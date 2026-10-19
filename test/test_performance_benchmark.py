"""Public benchmark CLI gates, using the same compiled build as pytest.

The end-to-end case takes ~10 s and <200 MiB on macOS. Fixtures are tiny;
the full 64 MiB log-cache probe deliberately exercises the production cap.
"""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from conftest import macs3_env


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'benchmark_performance.py'


def invoke(tmp_path, *args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                          cwd=tmp_path, env=macs3_env(tmp_path),
                          capture_output=True, text=True, timeout=90)


def test_generated_inputs_are_reproducible(tmp_path):
    manifests = []
    for name in ('one', 'two'):
        directory = tmp_path / name
        result = invoke(tmp_path, 'prepare', '--out', directory,
                        '--small', 300, '--large', 600, '--contigs', 9)
        assert result.returncode == 0, result.stdout + result.stderr
        manifests.append(json.loads((directory / 'manifest.json').read_text()))
    assert len(manifests[0]['cases']) == 9
    assert [case['sha256'] for case in manifests[0]['cases']] == [
        case['sha256'] for case in manifests[1]['cases']]


@pytest.mark.slow
def test_benchmark_cli_outputs_cache_and_failure_gates(tmp_path):
    pytest.importorskip('psutil', reason='benchmark runner requires optional psutil')
    data = tmp_path / 'data'
    result = invoke(tmp_path, 'prepare', '--out', data,
                    '--small', 900, '--large', 3000, '--contigs', 12)
    assert result.returncode == 0, result.stdout + result.stderr
    output = tmp_path / 'results'
    result = invoke(tmp_path, 'run', '--manifest', data / 'manifest.json',
                    '--out', output, '--label', 'test', '--revision', 'pytest-selected-build',
                    '--cli', SCRIPT.parents[1] / 'bin' / 'macs3',
                    '--budgets', '1,2', '--repeats', 1,
                    '--cases', 'small-bam,small-bampe,small-frag')
    assert result.returncode == 0, result.stdout + result.stderr
    metadata = json.loads((output / 'metadata.json').read_text())
    import MACS3
    assert Path(metadata['package']).resolve() == Path(MACS3.__file__).resolve()
    rows = json.loads((output / 'runs.json').read_text())
    assert len(rows) == 9
    assert all(row['self_peak_rss_bytes'] > 0 and row['samples'] > 0 for row in rows)
    assert all(row['execution_log'] and row['temp_remaining_bytes'] == 0 for row in rows)
    phases = json.loads((output / 'cache.json').read_text())['phases']
    assert phases[0]['live_table_bytes'] == 2 * 1024 * 8
    assert phases[0]['live_table_arrays'] == 2
    assert phases[-1]['live_table_bytes'] == 64 * 1024**2
    assert phases[-1]['live_table_arrays'] == 2
    comparison = tmp_path / 'comparison.md'
    result = invoke(tmp_path, 'compare', output, '--out', comparison)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'gates: PASS' in comparison.read_text()

    # A single numerical output change must fail, even with matching row counts.
    original = json.loads((output / 'runs.json').read_text())
    rows[0]['artifacts']['bench_peaks.narrowPeak']['sha256'] = '0' * 64
    (output / 'runs.json').write_text(json.dumps(rows))
    result = invoke(tmp_path, 'compare', output, '--out', comparison)
    assert result.returncode == 1
    assert 'scientific output bytes differ' in comparison.read_text()

    (output / 'runs.json').write_text(json.dumps(original[:-1]))
    result = invoke(tmp_path, 'compare', output, '--out', comparison)
    assert result.returncode == 1
    assert 'missing, extra or duplicate runs' in comparison.read_text()

    # An edited input cannot silently be compared against the recorded manifest.
    manifest = json.loads((data / 'manifest.json').read_text())
    Path(manifest['cases'][0]['treat']).write_bytes(b'changed')
    result = invoke(tmp_path, 'run', '--manifest', data / 'manifest.json',
                    '--out', tmp_path / 'invalid', '--label', 'invalid', '--revision', 'test')
    assert result.returncode != 0
    assert 'Input checksum mismatch' in result.stderr


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux deleted descriptor accounting')
def test_sample_counts_unlinked_temp_file_once(tmp_path):
    pytest.importorskip('psutil', reason='benchmark runner requires optional psutil')
    script = tmp_path / 'sample.py'
    script.write_text('''
import os
from pathlib import Path
import sys
sys.path.insert(0, sys.argv[1])
from benchmark_performance import process_sample
directory = Path(sys.argv[2])
path = directory / 'data'
path.write_bytes(b'x' * 8192)
first = os.open(path, os.O_RDONLY)
second = os.dup(first)
assert process_sample(os.getpid(), directory)[1] == 8192
path.unlink()
assert process_sample(os.getpid(), directory)[1] == 8192
os.close(first)
os.close(second)
assert process_sample(os.getpid(), directory)[1] == 0
''')
    directory = tmp_path / 'private-tmp'
    directory.mkdir()
    result = subprocess.run([sys.executable, str(script), str(SCRIPT.parent), str(directory)],
                            cwd=tmp_path, env=macs3_env(tmp_path),
                            capture_output=True, text=True, timeout=15)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.slow
def test_benchmark_timeout_is_a_failed_run(tmp_path):
    pytest.importorskip('psutil', reason='benchmark runner requires optional psutil')
    data = tmp_path / 'data'
    result = invoke(tmp_path, 'prepare', '--out', data, '--small', 900, '--large', 900)
    assert result.returncode == 0, result.stderr
    cli = tmp_path / 'hang.py'
    cli.write_text('''
import sys, time
if '--help' in sys.argv:
    print('--threads --process-backend')
else:
    time.sleep(30)
''')
    output = tmp_path / 'timeout'
    result = invoke(tmp_path, 'run', '--manifest', data / 'manifest.json',
                    '--out', output, '--label', 'timeout', '--revision', 'test',
                    '--cli', cli, '--budgets', '1', '--repeats', 1,
                    '--cases', 'small-frag', '--timeout', 3)
    assert result.returncode == 1, result.stdout + result.stderr
    row, = json.loads((output / 'runs.json').read_text())
    assert 'cli failed' in row['error']
    metrics, = output.glob('small-frag-*/cli.metrics.json')
    assert json.loads(metrics.read_text())['timed_out'] is True
    result = invoke(tmp_path, 'compare', output, '--out', tmp_path / 'failed.md')
    assert result.returncode == 1
