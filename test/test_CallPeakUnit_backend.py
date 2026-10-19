"""Exercise compiled public calls with bounded host-context subprocesses."""
import multiprocessing
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest


SCRIPT = r'''
import numpy as np
import multiprocessing
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

import MACS3.Signal.CallPeakUnit as unit
import MACS3.Signal.CallerExecution as execution
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Utilities.Execution import ExecutionPolicy, CPUAvailability
from MACS3.Utilities.Performance import DEFAULTS

def execution_policy(budget):
    return ExecutionPolicy(budget, budget, CPUAvailability(64, (("test host", 64),), ()))


def run(scenario, mode):
    pools = []
    real_executor = execution.ProcessPoolExecutor

    def executor(*args, **kwargs):
        assert kwargs['mp_context'].get_start_method() == 'fork'
        assert kwargs['initializer'] is unit.initialize_peak_worker
        assert len(kwargs['initargs']) == 1
        pools.append(kwargs['initargs'][0])
        return real_executor(*args, **kwargs)

    execution.ProcessPoolExecutor = executor

    def rows(peaks):
        keys = ['start', 'end', 'length', 'score', 'pileup',
                'pscore', 'qscore', 'fc', 'name']
        keys += (['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
                 if mode == 'broad' else ['summit'])
        return [(chrom, [tuple(row[key] for key in keys)
                         for row in peaks.peaks[chrom]])
                for chrom in sorted(peaks.peaks)]

    def call(backend):
        track = FWTrack(fw=50)
        for index in range(3):
            for pos in [1000 + index * 200] * (20 + index * 7) + [9000]:
                track.add_loc(('chr%d' % index).encode(), pos, 0)
        # Enough real work to reach pool/safety boundaries under the minimum.
        for index in range(3):
            track.add_loc_arrays(('chr%d' % index).encode(),
                np.arange(DEFAULTS.min_parallel_peak_locations // 3 + 1,
                          dtype='i4') * 300 + 10000, np.array([], dtype='i4'))
        track.finalize()
        caller = unit.CallerFromAlignments(
            track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[],
            lambda_bg=0.5, process_backend=backend,
            execution_policy=ExecutionPolicy(3, 3, CPUAvailability(64, (), ()),
                                             requested_process_backend=backend))
        try:
            if mode == 'broad':
                peaks = caller.call_broadpeaks(['q'], [3.0], [1.0], min_length=100)
            else:
                peaks = caller.call_peaks(['q'], [3.0], min_length=100,
                                          call_summits=mode == 'summits')
            result = rows(peaks)
            assert result and all(rows_for_chrom for _, rows_for_chrom in result)
            # Exercise reuse of the same caller and its precomputed q-table.
            if mode != 'broad':
                again = caller.call_peaks(['q'], [3.0], min_length=100,
                                          call_summits=mode == 'summits')
                assert result == rows(again)
            return result
        finally:
            caller.destroy()

    baseline = call('serial')
    assert not pools, 'explicit serial created a process pool'
    stop = threading.Event()
    background = None
    try:
        # Simulate Linux eligibility for individual rejection conditions.
        # The linux_fork case below alone validates real process execution.
        if scenario != 'linux_fork':
            sys.platform = 'linux'
        if scenario == 'lock_change':
            real_eligible = execution.peak_processes_allowed
            checks = 0
            def change_host_before_recheck(backend):
                nonlocal checks, background
                checks += 1
                if checks == 2:
                    background = threading.Thread(target=stop.wait)
                    background.start()
                return real_eligible(backend)
            execution.peak_processes_allowed = change_host_before_recheck
        if scenario == 'macos':
            sys.platform = 'darwin'
        elif scenario == 'unsupported':
            sys.platform = 'win32'
        elif scenario == 'no_fork':
            execution.multiprocessing.get_all_start_methods = lambda: ['spawn']
        elif scenario == 'background':
            background = threading.Thread(target=stop.wait)
            background.start()

        if scenario == 'non_main':
            with ThreadPoolExecutor(max_workers=1) as host:
                actual = host.submit(call, 'auto').result(timeout=30)
        else:
            actual = call('serial' if scenario == 'serial' else 'auto')
        assert actual == baseline, (actual, baseline)
        if scenario == 'linux_fork':
            assert len(pools) == (1 if mode == 'broad' else 3), len(pools)
        else:
            assert not pools, (scenario, len(pools))
    finally:
        stop.set()
        if background:
            background.join()


if __name__ == '__main__':
    scenario, mode = sys.argv[1:]
    if scenario == 'daemon':
        child = multiprocessing.get_context('spawn').Process(
            target=run, args=(scenario, mode), daemon=True)
        child.start()
        child.join(40)
        if child.is_alive():
            child.kill()
            child.join()
            raise AssertionError('daemon child timed out')
        assert child.exitcode == 0, child.exitcode
    else:
        run(scenario, mode)
'''


@pytest.mark.skipif(os.name != 'posix', reason='process-group cleanup requires POSIX')
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('scenario', [
    'serial', 'macos', 'unsupported', 'background', 'non_main', 'daemon',
    'no_fork', 'lock_change', 'linux_fork'])
def test_peak_process_backend(tmp_path, scenario, mode):
    import MACS3
    if scenario == 'linux_fork' and (
            sys.platform != 'linux' or
            'fork' not in multiprocessing.get_all_start_methods()):
        pytest.skip('real Linux fork validation requires Linux')
    env = os.environ.copy()
    env['TMPDIR'] = str(tmp_path)
    package_root = str(Path(MACS3.__file__).resolve().parent.parent)
    env['PYTHONPATH'] = package_root + (
        os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    script = tmp_path / 'backend_check.py'
    script.write_text(SCRIPT)
    process = subprocess.Popen(
        [sys.executable, str(script), scenario, mode], cwd=tmp_path, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=60)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        pytest.fail('backend subprocess timed out\n' + stdout + stderr)
    assert process.returncode == 0, stdout + stderr


@pytest.mark.parametrize('backend', ['', 'fork', 'spawn', 'forkserver', None])
def test_invalid_process_backend(backend):
    from MACS3.Signal.CallPeakUnit import CallerFromAlignments
    with pytest.raises((ValueError, TypeError), match='process_backend'):
        CallerFromAlignments(None, None, process_backend=backend)
