"""Public caller regressions run in subprocesses to bound native hangs."""
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest


SCRIPT = r'''
import sys
import threading
import logging
from concurrent.futures import ThreadPoolExecutor
import MACS3.Signal.CallPeakUnit as unit
from MACS3.Signal.FixWidthTrack import FWTrack

same_names, broad, failure, parallel = map(int, sys.argv[1:])
# Keep logging quiet while exercising overlapping callers.
logging.disable(logging.CRITICAL)


def caller(job):
    track = FWTrack(fw=50)
    prefix = 'chr' if same_names else 'job%d_chr' % job
    for index in range(2):
        chrom = (prefix + str(index)).encode()
        for pos in [1000 + job * 2000] * (12 + job * 7) + [9000]:
            track.add_loc(chrom, pos, 0)
    track.finalize()
    return unit.CallerFromAlignments(track, None, d=200, ctrl_d_s=[],
                                    ctrl_scaling_factor_s=[], lambda_bg=0.5)


def run(calc, job):
    if broad and job == 1:
        peaks = calc.call_broadpeaks(['q'], [3.0], [1.0], min_length=100)
    else:
        peaks = calc.call_peaks(['q'], [3.0], min_length=100)
    keys = ['start', 'end', 'score', 'pileup', 'pscore', 'qscore', 'fc']
    if broad and job == 1:
        keys += ['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
        chroms = sorted(peaks.peaks)
    else:
        keys += ['summit']
        chroms = sorted(peaks.get_chr_names())
    return [(chrom, [tuple(peak[key] for key in keys) for peak in
                    (peaks.peaks[chrom] if broad and job == 1
                     else peaks.get_data_from_chrom(chrom))])
            for chrom in chroms]



unit.n_available_workers = lambda: 1
baseline = []
for job in range(2):
    calc = caller(job)
    try:
        rows = run(calc, job)
        assert rows, 'fixture must call peaks'
        baseline.append(rows)
    finally:
        calc.destroy()
assert baseline[0] != baseline[1]

real_executor = unit.ProcessPoolExecutor
barrier = threading.Barrier(2, timeout=15)


class FailingExecutor(real_executor):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.submissions = 0

    def submit(self, *args, **kwargs):
        self.submissions += 1
        if self.submissions == 2:
            raise RuntimeError('injected submission failure')
        return super().submit(*args, **kwargs)


if failure:
    unit.n_available_workers = lambda: 2
    unit.ProcessPoolExecutor = FailingExecutor
    calc = caller(0)
    try:
        try:
            run(calc, 0)
        except RuntimeError as exc:
            assert str(exc) == 'injected submission failure'
        else:
            raise AssertionError('failure was not injected')
        unit.ProcessPoolExecutor = real_executor
        assert run(calc, 0) == baseline[0]
    finally:
        calc.destroy()
    later = caller(1)
    try:
        assert run(later, 1) == baseline[1]
    finally:
        later.destroy()
else:
    unit.ProcessPoolExecutor = real_executor
    local = threading.local()

    def workers():
        if not getattr(local, 'entered', False):
            local.entered = True
            barrier.wait()
        return 2 if parallel else 1

    unit.n_available_workers = workers
    for repeat in range(3):
        callers = [caller(job) for job in range(2)]
        try:
            with ThreadPoolExecutor(max_workers=2) as executor:
                results = list(executor.map(lambda job: run(callers[job], job), range(2)))
                assert results == baseline, (results, baseline)
                # Reuse each caller's q-table/cache as well as testing fresh instances.
                results = list(executor.map(lambda job: run(callers[job], job), range(2)))
                assert results == baseline, (results, baseline)
        finally:
            for calc in callers:
                calc.destroy()
'''


@pytest.mark.skipif(os.name != 'posix', reason='bounded process-group cleanup requires POSIX')
@pytest.mark.parametrize('same_names', [False, True])
@pytest.mark.parametrize('broad', [False, True])
@pytest.mark.parametrize('scenario', ['host_fallback', 'serial', 'failure'])
def test_independent_callers(tmp_path, same_names, broad, scenario):
    import multiprocessing
    import MACS3
    if scenario == 'failure' and (sys.platform != 'linux' or
                                 'fork' not in multiprocessing.get_all_start_methods()):
        pytest.skip('submission failure requires the supported Linux fork backend')
    env = os.environ.copy()
    env['TMPDIR'] = str(tmp_path)
    # Match the package pytest imported, including non-editable wheel installs.
    # Running -c from the checkout would put Cython source ahead of that build.
    package_root = str(Path(MACS3.__file__).resolve().parent.parent)
    env['PYTHONPATH'] = package_root + (
        os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    process = subprocess.Popen(
        [sys.executable, '-c', SCRIPT, str(int(same_names)), str(int(broad)),
         str(int(scenario == 'failure')), str(int(scenario == 'host_fallback'))],
        cwd=tmp_path, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=60)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        pytest.fail('caller subprocess timed out\n' + stdout + stderr)
    assert process.returncode == 0, stdout + stderr
