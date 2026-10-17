"""Public execution policy contracts and resource detection at OS boundaries."""

from dataclasses import FrozenInstanceError, replace
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from MACS3.Utilities import Execution


@pytest.fixture
def cpu_host(monkeypatch):
    """Independent of the machine running the suite, including Linux quotas."""
    monkeypatch.setattr(Execution, 'os', SimpleNamespace(cpu_count=lambda: 32, environ=os.environ))
    monkeypatch.setattr(Execution, 'sys', SimpleNamespace(platform='darwin'))
    monkeypatch.delenv('MACS3_NUM_THREADS', raising=False)


def test_policy_precedence_clamping_and_snapshot(cpu_host, monkeypatch):
    assert Execution.resolve_execution_policy().threads == 8
    monkeypatch.setenv('MACS3_NUM_THREADS', '12')
    policy = Execution.resolve_execution_policy()
    assert policy.threads == 12 and policy.source == 'MACS3_NUM_THREADS'
    assert policy.process_workers == 11
    assert Execution.resolve_execution_policy(24).threads == 24
    assert Execution.resolve_execution_policy(64).threads == 32
    with pytest.raises(FrozenInstanceError):
        policy.threads = 1
    with pytest.raises(FrozenInstanceError):
        policy.availability.count = 1
    monkeypatch.setenv('MACS3_NUM_THREADS', 'bad')
    assert Execution.resolve_execution_policy(2).threads == 2
    assert Execution.resolve_execution_policy(execution_policy=policy) is policy
    monkeypatch.setattr(Execution.os, 'cpu_count', lambda: 1)
    assert Execution.resolve_execution_policy(execution_policy=policy) is policy
    assert policy.threads == 12
    assert Execution.resolve_execution_policy(8).threads == 1
    assert Execution.resolve_execution_policy(8).process_backend == 'serial'
    assert Execution.resolve_execution_policy(8).process_workers == 0


@pytest.mark.parametrize('value', [0, -1, True, False, 1.5, '2', 1 << 31])
def test_invalid_explicit_budget(cpu_host, value):
    with pytest.raises(ValueError, match='positive integer'):
        Execution.resolve_execution_policy(value)


@pytest.mark.parametrize('value', ['', '0', '-1', '1.5', 'yes', ' 2', '+2', '٢',
                                   str(1 << 31), pytest.param('9' * 5000, id='too-many-digits')])
def test_invalid_environment_budget(cpu_host, monkeypatch, value):
    monkeypatch.setenv('MACS3_NUM_THREADS', value)
    with pytest.raises(ValueError):
        Execution.resolve_execution_policy()


def test_policy_conflicts_and_backend(cpu_host):
    policy = Execution.resolve_execution_policy(4, 'serial')
    assert policy.threads == 4 and policy.process_workers == 0
    assert Execution.resolve_execution_policy(4, 'serial', execution_policy=policy) is policy
    for overrides in ({'threads': 2}, {'process_backend': 'auto'}):
        with pytest.raises(ValueError, match='conflicts'):
            Execution.resolve_execution_policy(execution_policy=policy, **overrides)
    for backend in (None, '', 'fork', 'spawn'):
        with pytest.raises(ValueError, match='process_backend'):
            Execution.resolve_execution_policy(process_backend=backend)
    with pytest.raises(TypeError):
        Execution.resolve_execution_policy(execution_policy=object())
    with pytest.raises(ValueError, match='clamped'):
        replace(policy, threads=2)
    with pytest.raises(ValueError, match='immutable'):
        Execution.CPUAvailability(2, [('host', 2)], ())
    assert Execution.resolve_execution_policy(default_threads=1).threads == 1


@pytest.mark.parametrize('host,process,affinity,expected', [
    (32, 8, set(range(4)), 4), (2, 8, set(range(4)), 2),
    (None, None, set(), 1), (None, 3, set(), 3),
    (None, None, {0, 1}, 2), (0, -2, set(), 1),
])
def test_resource_signals(cpu_host, monkeypatch, host, process, affinity, expected):
    monkeypatch.setattr(Execution.os, 'cpu_count', lambda: host)
    monkeypatch.setattr(Execution.os, 'process_cpu_count', lambda: process, raising=False)
    monkeypatch.setattr(Execution.os, 'sched_getaffinity', lambda pid: affinity, raising=False)
    assert Execution.available_cpus() == expected


@pytest.mark.parametrize('failure', [OSError, PermissionError, NotImplementedError])
def test_resource_query_failure(cpu_host, monkeypatch, failure):
    def unavailable(*args):
        raise failure('unavailable')
    monkeypatch.setattr(Execution.os, 'sched_getaffinity', unavailable, raising=False)
    monkeypatch.setattr(Execution.os, 'process_cpu_count', unavailable, raising=False)
    assert Execution.available_cpus() == 32
    monkeypatch.setattr(Execution.os, 'cpu_count', unavailable)
    assert Execution.available_cpus() == 1
    assert Execution.detect_cpu_resources().notes
    monkeypatch.delattr(Execution.os, 'cpu_count')
    assert Execution.available_cpus() == 1


@pytest.mark.parametrize('kind', ['v1', 'v2'])
@pytest.mark.parametrize('quota,expected', [
    ('200000 100000', 2), ('250000 100000', 2), ('50000 100000', 1),
    ('unlimited', 32), ('bad', 32), ('0 100000', 32), ('200000 0', 32),
    ('missing', 32),
])
def test_quota_detection(cpu_host, monkeypatch, tmp_path, kind, quota, expected):
    mount = tmp_path / 'cpu mount'
    child = mount / 'job'
    child.mkdir(parents=True)
    files = {
        '/proc/self/cgroup': ('0::/batch/job\n' if kind == 'v2'
                              else '7:cpu,cpuacct:/batch/job\n'),
        '/proc/self/mountinfo': (
            '5 1 0:1 /batch ' + str(mount).replace(' ', r'\040') +
            (' rw - cgroup2 cgroup rw\n' if kind == 'v2'
             else ' rw - cgroup cgroup rw,cpu,cpuacct\n')),
    }
    if quota != 'missing':
        if kind == 'v2':
            (child / 'cpu.max').write_text('max 100000' if quota == 'unlimited' else quota)
        else:
            values = ('-1 100000' if quota == 'unlimited' else quota).split()
            (child / 'cpu.cfs_quota_us').write_text(values[0])
            (child / 'cpu.cfs_period_us').write_text(values[1] if len(values) == 2 else 'bad')
    real_read = Path.read_text
    def read(path, *args, **kwargs):
        if str(path) in files:
            return files[str(path)]
        if str(path).startswith('/proc/'):
            raise PermissionError()
        return real_read(path, *args, **kwargs)
    monkeypatch.setattr(Path, 'read_text', read)
    monkeypatch.setattr(Execution, 'sys', SimpleNamespace(platform='linux'))
    assert Execution.available_cpus() == expected
    # A visible parent is more restrictive than a child/unlimited/invalid leaf.
    if kind == 'v2':
        (mount / 'cpu.max').write_text('100000 100000')
    else:
        (mount / 'cpu.cfs_quota_us').write_text('100000')
        (mount / 'cpu.cfs_period_us').write_text('100000')
    assert Execution.available_cpus() == 1
    # Membership in a cgroup namespace can name '/' below a non-root mount root.
    files['/proc/self/cgroup'] = '0::/\n' if kind == 'v2' else '7:cpu:/\n'
    assert Execution.available_cpus() == 1


def test_quota_metadata_unavailable(cpu_host, monkeypatch):
    monkeypatch.setattr(Execution, 'sys', SimpleNamespace(platform='linux'))
    def denied(*args, **kwargs):
        raise PermissionError()
    monkeypatch.setattr(Path, 'read_text', denied)
    assert Execution.available_cpus() == 32
    assert 'cgroup metadata unavailable' in Execution.detect_cpu_resources().notes


def test_resolution_does_not_mutate_environment(cpu_host, monkeypatch):
    monkeypatch.setenv('OMP_NUM_THREADS', '7')
    before = dict(os.environ)
    policies = [Execution.resolve_execution_policy(x) for x in (1, 2, 8, 16)]
    assert [p.threads for p in policies] == [1, 2, 8, 16]
    assert dict(os.environ) == before


def test_command_api_preserves_supplied_policy(cpu_host, macs3_argparser, monkeypatch):
    from MACS3.Utilities.OptValidator import opt_validate_callpeak
    policy = Execution.resolve_execution_policy(4, 'serial')
    options = macs3_argparser.parse_args(['callpeak', '-t', 't.bed'])
    options.execution_policy = policy
    monkeypatch.setenv('MACS3_NUM_THREADS', 'bad unused environment value')
    assert opt_validate_callpeak(options).execution_policy is policy


def test_clamped_callvar_options_can_be_revalidated(cpu_host, macs3_argparser):
    from MACS3.Utilities.OptValidator import opt_validate_callvar
    options = macs3_argparser.parse_args(
        ['callvar', '-t', 't.bam', '-b', 'p.bed', '-o', 'o.vcf', '--threads', '64'])
    opt_validate_callvar(options)
    policy = options.execution_policy
    assert options.np == 64 and policy.threads == 32
    assert opt_validate_callvar(options).execution_policy is policy


@pytest.mark.parametrize('backend,budget', [('auto', 1), ('serial', 4)])
def test_variant_serial_policy_starts_no_pool(cpu_host, monkeypatch, backend, budget):
    from MACS3.Commands import callvar_cmd
    policy = Execution.resolve_execution_policy(budget, backend)
    def forbidden(*args, **kwargs):
        raise AssertionError('serial policy tried to create a process context')
    monkeypatch.setattr(callvar_cmd.mp, 'get_context', forbidden)
    assert callvar_cmd.map_variant_ranges(sum, [(1, 2), (3, 4)], policy) == [3, 7]


def test_variant_host_fallback_starts_no_pool(cpu_host, monkeypatch):
    from MACS3.Commands import callvar_cmd
    policy = Execution.resolve_execution_policy(4)
    def forbidden(*args, **kwargs):
        raise AssertionError('unsafe host tried to create a process context')
    monkeypatch.setattr(callvar_cmd.mp, 'get_context', forbidden)
    assert callvar_cmd.map_variant_ranges(sum, [(1, 2), (3, 4)], policy) == [3, 7]


def test_variant_fork_worker_budget(monkeypatch):
    import multiprocessing
    import sys
    if sys.platform != 'linux' or 'fork' not in multiprocessing.get_all_start_methods():
        pytest.skip('real fork worker validation requires Linux')
    from MACS3.Commands import callvar_cmd
    policy = Execution.resolve_execution_policy(4)
    contexts = []
    real_context = multiprocessing.get_context
    class Context:
        def Pool(self, workers):
            contexts.append(workers)
            return real_context('fork').Pool(workers)
    monkeypatch.setattr(callvar_cmd.mp, 'get_context', lambda name: Context())
    expected = callvar_cmd.map_variant_ranges(sum, [(1, 2), (3, 4)], policy, timeout=10)
    assert expected == [3, 7]
    if policy.process_workers:
        assert contexts == [min(policy.process_workers, 2)]
    assert not multiprocessing.active_children()
