"""Public caller/CLI regressions for bounded work and empty chromosomes."""

import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest

from conftest import macs3_env, no_core_dumps
from MACS3.Utilities.Performance import DEFAULTS


def probe(tmp_path, **config):
    script = tmp_path / 'chromosome_batching_probe.py'
    script.write_text(Path(__file__).with_name(script.name).read_text())
    proc = subprocess.Popen(
        [sys.executable, str(script), json.dumps(config)], cwd=tmp_path,
        env=macs3_env(tmp_path), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, start_new_session=True, preexec_fn=no_core_dumps)
    try:
        stdout, stderr = proc.communicate(timeout=45)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        stdout, stderr = proc.communicate()
        pytest.fail('chromosome caller timed out\n' + stdout + stderr)
    assert proc.returncode == 0, stdout + stderr
    return json.loads(stdout)


@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
@pytest.mark.parametrize('phase', ['histogram', 'calling'])
@pytest.mark.parametrize('workload', [
    'tiny', 'many_tiny', 'below', 'at', 'above', 'single', 'skewed', 'mixed', 'rounding'])
def test_public_work_gate_and_submission_bound(tmp_path, kind, phase, workload):
    probe(tmp_path, action='schedule', kind=kind, phase=phase, workload=workload)


@pytest.mark.parametrize('phase', ['histogram', 'calling'])
@pytest.mark.parametrize('budget,backend', [(1, 'auto'), (2, 'auto'), (4, 'serial')])
def test_scheduler_obeys_execution_policy(tmp_path, phase, budget, backend):
    probe(tmp_path, action='schedule', kind='fw', phase=phase, workload='skewed',
          budget=budget, backend=backend)


@pytest.mark.parametrize('phase', ['histogram', 'calling'])
@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
def test_distinct_control_records_count_toward_work(tmp_path, phase, kind):
    result = probe(tmp_path, action='schedule', kind=kind, phase=phase,
                   workload='below', control=True)
    assert result['tasks'] > 0


@pytest.mark.parametrize('phase', ['histogram', 'calling'])
def test_frag_multiplicities_do_not_count_as_stored_work(tmp_path, phase):
    result = probe(tmp_path, action='schedule', kind='frag', phase=phase,
                   workload='tiny', multiplicity=10000)
    assert result['locations'] == 90 and result['tasks'] == 0


@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('workload', ['empty', 'disjoint', 'tiny', 'mixed', 'normal', 'skewed'])
def test_caller_results_and_order_across_budgets(tmp_path, kind, mode, workload):
    probe(tmp_path, action='outputs', kind=kind, mode=mode, workload=workload)


@pytest.mark.parametrize('kind', ['fw', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
def test_retained_empty_chromosomes(tmp_path, kind, mode):
    # PETrackI.sample_percent divides by zero when the whole track is empty;
    # its genuine-empty and retained mixed cases are covered separately.
    probe(tmp_path, action='outputs', kind=kind, mode=mode, workload='retained')


@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('workload', ['empty', 'mixed'])
def test_empty_bedgraph_and_cutoff_artifacts(tmp_path, kind, mode, workload):
    probe(tmp_path, action='outputs', kind=kind, mode=mode, workload=workload,
          artifacts=True)


@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
def test_one_sided_empty_chromosomes(tmp_path, kind, mode):
    probe(tmp_path, action='one_sided', kind=kind, mode=mode)


@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('symbol', ['p', 'q'])
def test_controlled_outputs_across_budgets(tmp_path, kind, mode, symbol):
    probe(tmp_path, action='outputs', kind=kind, mode=mode, workload='normal',
          control=True, symbol=symbol)


@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('large', [False, True])
def test_cli_outputs_without_bedgraph_or_cutoff(tmp_path, run_macs3, mode, large):
    count = DEFAULTS.min_parallel_peak_locations // 3 + 1 if large else 30
    path = tmp_path / 'reads.bedpe'
    path.write_bytes(b''.join(
        b'chr%d\t%d\t%d\n' % (chrom, start, start + 200)
        for chrom in range(3) for start in [1000 + chrom * 400] * count + [9000]))
    expected = None
    for budget in (1, 4):
        out = tmp_path / str(budget)
        args = ['callpeak', '-t', path, '-f', 'BEDPE', '-g', '100000000',
                '-n', 'peaks', '--outdir', out, '--keep-dup', 'all',
                '--threads', str(budget)]
        if mode != 'narrow':
            args.append('--broad' if mode == 'broad' else '--call-summits')
        result = run_macs3(args, timeout=45)
        assert result.returncode == 0, result.stderr
        files = {p.name: p.read_bytes() for p in out.iterdir() if p.suffix != '.xls'}
        files['xls-data'] = b'\n'.join(line for line in (out / 'peaks_peaks.xls').read_bytes().splitlines()
                                      if line and not line.startswith(b'#'))
        assert files and any(data for data in files.values())
        if expected is None:
            expected = files
        else:
            assert files == expected
        if not large and budget == 4 and sys.platform == 'linux':
            # A one-CPU allocation selects the serial policy before work gating.
            assert ('threads=1 (including caller)' in result.stderr or
                    'serial (workload:' in result.stderr)
