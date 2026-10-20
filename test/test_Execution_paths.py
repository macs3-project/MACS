"""Compiled public API and CLI regressions for invocation CPU budgets."""

import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from types import SimpleNamespace

import pytest

from conftest import macs3_env
from test_Parser_bam_truncation import HEADER, record, write_bam
from test_Parser_frag_eof import write_frag


@pytest.fixture
def linux_thread_counter(tmp_path, monkeypatch):
    """Select Linux observation against controlled procfs task records."""
    import execution_probe
    monkeypatch.setattr(execution_probe, 'sys', SimpleNamespace(platform='linux'))
    monkeypatch.setattr(execution_probe, 'Path', lambda path: tmp_path)
    return execution_probe.native_thread_count


def test_linux_thread_count_excludes_tasks_already_exiting(tmp_path, linux_thread_counter):
    # pthread_join can finish while the kernel still lists a PF_EXITING task.
    # Its state need not yet be Z; names can contain spaces and parentheses.
    for tid, state, flags in [(1, 'R', 4194304), (2, 'S', 4194304),
                              (3, 'R', 4194308), (4, 'Z', 4194308)]:
        task = tmp_path / str(tid)
        task.mkdir()
        (task / 'stat').write_text(f'{tid} (inflate ) work) {state} 1 1 1 0 -1 {flags} 0\n')
    # Linux comm is raw bytes, not necessarily valid UTF-8.
    (tmp_path / '2' / 'stat').write_bytes(b'2 (worker\xff ) A) S 1 1 1 0 -1 4194304 0\n')
    assert linux_thread_counter() == 2
    # An additional live worker must still fail a two-thread cleanup bound.
    (tmp_path / '3' / 'stat').write_text('3 (inflate worker) S 1 1 1 0 -1 4194304 0\n')
    assert linux_thread_counter() == 3


def test_linux_thread_count_handles_a_task_removed_during_observation(tmp_path, linux_thread_counter):
    (tmp_path / 'gone').mkdir()  # stat disappeared after directory enumeration
    live = tmp_path / 'live'
    live.mkdir()
    (live / 'stat').write_text('1 (main) R 1 1 1 0 -1 4194304 0\n')
    assert linux_thread_counter() == 1


@pytest.mark.parametrize('error', [ProcessLookupError, PermissionError])
def test_linux_thread_count_distinguishes_vanished_and_unreadable_tasks(
        tmp_path, linux_thread_counter, monkeypatch, error):
    task = tmp_path / '1'
    task.mkdir()

    def read_bytes(path):
        raise error('task stat read failed')

    monkeypatch.setattr(Path, 'read_bytes', read_bytes)
    if error is ProcessLookupError:
        assert linux_thread_counter() == 0
    else:
        with pytest.raises(PermissionError):
            linux_thread_counter()


def test_linux_thread_count_does_not_hide_invalid_procfs_data(tmp_path, linux_thread_counter):
    task = tmp_path / '1'
    task.mkdir()
    (task / 'stat').write_text('1 (main) R 1 1 1 0 -1 invalid 0\n')
    with pytest.raises(ValueError):
        linux_thread_counter()


def probe(tmp_path, scenario, budget, *args):
    script = tmp_path / 'execution_probe.py'
    script.write_text(Path(__file__).with_name('execution_probe.py').read_text())
    process = subprocess.Popen(
        [sys.executable, str(script), scenario, str(budget), *map(str, args)],
        cwd=tmp_path, env=macs3_env(tmp_path), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=45)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        pytest.fail('execution subprocess timed out\n' + stdout + stderr)
    assert process.returncode == 0, stdout + stderr
    return json.loads(stdout)


@pytest.mark.parametrize('encoding', ['plain', 'gzip', 'bgzf'])
@pytest.mark.parametrize('kind', ['frag', 'bam', 'bampe'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_parser_output_across_budgets(tmp_path, encoding, kind, operation):
    path = tmp_path / 'input'
    if kind == 'frag':
        data = b''.join(b'chr%d\t%d\t%d\tBC%d\t2\n' %
                        (i % 4, i * 3, i * 3 + 40, i % 3) for i in range(12000))
        write_frag(path, data, encoding)
        expected = 24000
    else:
        write_bam(path, HEADER + record() * 4000,
                  'raw' if encoding == 'plain' else encoding)
        expected = 4000
    results = [probe(tmp_path, 'parser', budget, kind, operation, path)
               for budget in (1, 2, 4)]
    assert results[0] == results[1] == results[2]
    assert results[0][0] == expected


@pytest.mark.skipif(sys.platform not in ('linux', 'darwin'),
                    reason='native thread observation uses Linux /proc or Mach')
def test_bgzf_native_threads_obey_budget(tmp_path):
    # ~5 MiB decompressed; enough to observe worker lifetime across windows.
    path = tmp_path / 'fragments.bgzf'
    data = b''.join(b'chr%d\t%d\t%d\tBC%d\t2\n' %
                    (i % 4, i * 3, i * 3 + 40, i % 3) for i in range(180000))
    write_frag(path, data, 'bgzf')
    for budget in (1, 2, 4):
        result = probe(tmp_path, 'native', budget, path)
        assert result['extra_native_threads'] <= result['effective'] - 1
        assert result['stream_open_samples'] > 0


@pytest.mark.parametrize('scenario', ['track', 'model', 'cleanup'])
def test_public_phase_output_across_budgets(tmp_path, scenario):
    results = [probe(tmp_path, scenario, budget) for budget in (1, 2, 4)]
    assert results[0] == results[1] == results[2]


@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('control', ['yes', 'no'])
@pytest.mark.parametrize('symbol', ['p', 'q'])
def test_caller_output_across_budgets(tmp_path, mode, control, symbol):
    results = [probe(tmp_path, 'caller', budget, mode, control, symbol)
               for budget in (1, 2, 4)]
    assert results[0] == results[1] == results[2]


@pytest.mark.parametrize('command', ['callpeak', 'predictd', 'pileup', 'filterdup',
                                    'randsample', 'refinepeak', 'hmmratac', 'callvar'])
def test_cli_budget_help(run_macs3, command):
    result = run_macs3([command, '--help'])
    assert result.returncode == 0
    assert '--threads' in result.stdout and 'MACS3_NUM_THREADS' in result.stdout


@pytest.mark.parametrize('value', ['0', '-1', '2147483648', '1.5'])
def test_cli_rejects_invalid_budget(run_macs3, tmp_path, value):
    result = run_macs3(['callpeak', '-t', 'missing', '--threads', value])
    assert result.returncode != 0
    assert 'positive integer' in result.stderr or 'invalid int value' in result.stderr


@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
def test_cli_peak_artifacts_across_budgets(tmp_path, run_macs3, mode):
    data = b''.join(b'chr%d\t%d\t%d\n' % (c, p, p + 150)
                    for c in range(3) for p in [1000 + c * 100] * 30 + [9000])
    path = tmp_path / 'reads.bedpe'
    path.write_bytes(data)
    outputs = []
    for budget in (1, 2, 4):
        outdir = tmp_path / str(budget)
        args = ['callpeak', '-t', path, '-f', 'BEDPE', '-g', '100000',
                '-n', 'peaks', '--outdir', outdir, '--threads', budget,
                '--keep-dup', 'all', '-B', '--cutoff-analysis']
        if mode != 'narrow':
            args.append('--broad' if mode == 'broad' else '--call-summits')
        result = run_macs3(args, timeout=45)
        assert result.returncode == 0, result.stderr
        assert 'Execution: requested=%d' % budget in result.stderr
        stable = {p.name: p.read_bytes() for p in outdir.iterdir() if p.suffix != '.xls'}
        xls = outdir / 'peaks_peaks.xls'
        stable['xls-data'] = b'\n'.join(line for line in xls.read_bytes().splitlines()
                                      if line and not line.startswith(b'#'))
        assert stable
        outputs.append(stable)
    assert outputs[0] == outputs[1] == outputs[2]


def test_validation_under_optimized_python(tmp_path):
    script = tmp_path / 'optimized.py'
    script.write_text('''
from MACS3.Utilities.Execution import resolve_execution_policy
for value in (0, -1, True, 1.5, 1 << 31):
    try:
        resolve_execution_policy(value)
    except ValueError:
        continue
    raise SystemExit('invalid budget accepted under python -O')
''')
    result = subprocess.run([sys.executable, '-O', str(script)], cwd=tmp_path,
                            env=macs3_env(tmp_path), capture_output=True, text=True,
                            timeout=15)
    assert result.returncode == 0, result.stdout + result.stderr


def test_cleanup_closes_unclaimed_files_on_startup_failure(tmp_path, monkeypatch):
    import threading
    from MACS3.Signal.CallPeakUnit import close_unlinked_files_later
    from MACS3.Utilities.Execution import ExecutionPolicy, CPUAvailability
    policy = ExecutionPolicy(4, 4, CPUAvailability(4, (), ()))
    fds = []
    for i in range(8):
        path = tmp_path / str(i)
        fds.append(os.open(path, os.O_WRONLY | os.O_CREAT, 0o600))
        path.unlink()
    release = threading.Event()
    started = []
    original = threading.Thread

    def factory(*args, **kwargs):
        if started:
            release.set()
            raise RuntimeError('injected cleanup startup failure')
        target = kwargs['target']
        def held():
            if not release.wait(5):
                raise RuntimeError('cleanup worker was not released')
            target()
        thread = original(target=held)
        started.append(thread)
        return thread

    monkeypatch.setattr(threading, 'Thread', factory)
    try:
        with pytest.raises(RuntimeError, match='injected cleanup startup failure'):
            close_unlinked_files_later(fds, [1 << 20] * len(fds), policy)
        assert started and all(not thread.is_alive() for thread in started)
        for fd in fds:
            with pytest.raises(OSError):
                os.fstat(fd)
    finally:
        release.set()
        for thread in started:
            thread.join(5)
        for fd in fds:
            try:
                os.close(fd)
            except OSError:
                pass
