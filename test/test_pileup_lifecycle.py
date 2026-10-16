"""Bounded public peak-call probes for ownership of temporary pileup files."""
import multiprocessing
import os
import signal
import subprocess
import sys

import pytest

from conftest import macs3_env


SCRIPT = r'''
from concurrent.futures import Future
import multiprocessing
import os
from pathlib import Path
import sys
import tempfile
import signal
import threading
import time

import numpy as np
import MACS3.Signal.CallPeakUnit as unit
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Utilities.Execution import ExecutionPolicy, CPUAvailability

def execution_policy(budget):
    return ExecutionPolicy(budget, budget, CPUAvailability(64, (("test host", 64),), ()))

scenario, mode, backend = sys.argv[1:]
root = Path(tempfile.gettempdir())
sentinel = root / 'unrelated-user-file'
sentinel.write_bytes(b'keep me')
before = set(root.iterdir())
pools = []
original_executor = unit.ProcessPoolExecutor

def recording_executor(*args, **kwargs):
    pools.append(kwargs)
    return original_executor(*args, **kwargs)

unit.ProcessPoolExecutor = recording_executor

def make(background=0.5, selected_backend=None):
    selected_backend = backend if selected_backend is None else selected_backend
    track = FWTrack(fw=50)
    for chrom in [b'chr1', b'chr2']:
        for pos in [1000] * 20 + [9000]:
            track.add_loc(chrom, pos, 0)
    track.finalize()
    return unit.CallerFromAlignments(
        track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[],
        lambda_bg=background, process_backend=selected_backend,
        execution_policy=ExecutionPolicy(3, 3, CPUAvailability(64, (), ()),
                                         requested_process_backend=selected_backend),
        cutoff_analysis_filename=str(root.parent / 'cutoff.tsv'),
        save_bedGraph=scenario == 'bedgraph_retry',
        bedGraph_treat_filename=str(root.parent / 'treat.bdg'),
        bedGraph_control_filename=str(root.parent / 'ctrl.bdg'))

def call(caller, symbol='q'):
    if mode == 'broad':
        return caller.call_broadpeaks([symbol], [3.0], [1.0], min_length=100,
                                      cutoff_analysis=scenario == 'cutoff_retry')
    return caller.call_peaks([symbol], [3.0], min_length=100,
                             cutoff_analysis=scenario == 'cutoff_retry')

def rows(peaks):
    keys = ['start', 'end', 'length', 'score', 'pileup',
            'pscore', 'qscore', 'fc', 'name']
    keys += (['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
             if mode == 'broad' else ['summit'])
    return [(chrom, [tuple(row[key] for key in keys) for row in peaks.peaks[chrom]])
            for chrom in sorted(peaks.peaks)]

def clean():
    assert set(root.iterdir()) == before, list(root.rglob('*'))
    assert sentinel.read_bytes() == b'keep me'

baseline_caller = make()
expected = rows(call(baseline_caller))
assert expected and all(p for _, p in expected)
assert bool(pools) == (backend == 'auto'), 'unexpected production process selection'
output_names = (['cutoff.tsv'] if scenario == 'cutoff_retry' else
                ['treat.bdg', 'ctrl.bdg'] if scenario == 'bedgraph_retry' else [])
expected_outputs = {name: (root.parent / name).read_bytes() for name in output_names}
baseline_caller.destroy()
clean()
fault_backend = ('auto' if scenario.startswith('pool_') or scenario in (
    'interrupt', 'shutdown_interrupt', 'shutdown_signal', 'shutdown_custom_signal')
    else backend)
caller = make(0.0 if scenario == 'zero_lambda' else 0.5, fault_backend)
real_write = unit.write_pileup_file
real_executor = unit.ProcessPoolExecutor
real_eligible = unit.peak_processes_allowed
real_completed = unit.as_completed
real_stat_worker = unit.pscore_stat_worker
closed_fds = []
shutdowns = []
signal_deliveries = []
saved_handler = signal.getsignal(signal.SIGINT)
if scenario == 'shutdown_custom_signal':
    def host_handler(signum, frame):
        assert shutdowns and getattr(shutdowns[-1], 'finished', False)
        signal_deliveries.append(signum)
    signal.signal(signal.SIGINT, host_handler)
expected_handler = signal.getsignal(signal.SIGINT)
created_files = []
error = KeyboardInterrupt('injected interruption') if scenario in (
    'interrupt', 'write_interrupt', 'shutdown_interrupt') else RuntimeError('injected task failure')

if scenario.startswith('pool_') or scenario in (
        'interrupt', 'shutdown_interrupt', 'shutdown_signal', 'shutdown_custom_signal'):
    # Deterministic process-boundary failures. No private caller methods or
    # attributes are inspected. Real compiled precomputation and retry run.
    if scenario == 'pool_peak_partial':
        assert rows(call(caller)) == expected
    unit.peak_processes_allowed = lambda backend: True
    unit.as_completed = lambda futures: iter(futures)

    class FaultPool:
        def __init__(self, **kwargs):
            self.count = 0
            self.directory = next(p for p in root.iterdir() if p not in before)
            if scenario == 'pool_construct':
                self.new_file()
                raise error

        def new_file(self):
            path = self.directory / ('unreturned-%d' % len(created_files))
            path.write_bytes(b'worker pileup')
            created_files.append(path)
            return os.fsencode(path)

        def submit(self, fn, task):
            self.count += 1
            path = self.new_file()
            if scenario == 'pool_submit' and self.count == 2:
                raise error
            future = Future()
            if self.count == 1 and scenario in ('pool_partial', 'pool_peak_partial'):
                chrom = b'chr1'
                if scenario == 'pool_peak_partial':
                    future.set_result([(chrom, path, [])])
                else:
                    keys = np.array([(20 << 32) | int(
                        np.array([0.5], dtype='f4').view('u4')[0])], dtype='u8')
                    future.set_result(([(chrom, path)], keys,
                                       np.array([0], dtype='i8'), np.array([200], dtype='i8')))
            elif scenario == 'shutdown_interrupt':
                future.set_result(([], np.array([], dtype='u8'),
                                   np.array([], dtype='i8'), np.array([], dtype='i8')))
            else:
                future.set_exception(error)
            return future

        def shutdown(self, wait, cancel_futures):
            assert wait and cancel_futures
            shutdowns.append(self)
            if scenario == 'shutdown_interrupt' and len(shutdowns) == 1:
                raise error
            if scenario in ('shutdown_signal', 'shutdown_custom_signal') and len(shutdowns) == 1:
                os.kill(os.getpid(), signal.SIGINT)
            # Models a running worker completing during join. Cleanup must
            # happen after this write, including submission/result failures.
            self.new_file()
            self.finished = True

    unit.ProcessPoolExecutor = FaultPool

if scenario == 'real_shutdown_signal':
    def delayed_stat_worker(batch):
        if batch[0] == 0:
            deadline = time.monotonic() + 5
            while not (root.parent / 'slow-worker-started').exists():
                assert time.monotonic() < deadline, 'second worker did not start'
                time.sleep(0.01)
            raise RuntimeError('worker failure before shutdown')
        result = real_stat_worker(batch)
        (root.parent / 'slow-worker-started').write_text('started')
        time.sleep(0.5)
        (root.parent / 'slow-worker-finished').write_text('finished')
        return result

    class InterruptedShutdown(original_executor):
        def shutdown(self, **kwargs):
            timer = threading.Timer(0.1, lambda: os.kill(os.getpid(), signal.SIGINT))
            timer.start()
            try:
                return super().shutdown(**kwargs)
            finally:
                timer.join()

    unit.pscore_stat_worker = delayed_stat_worker
    unit.ProcessPoolExecutor = InterruptedShutdown

if scenario in ('write_error', 'write_interrupt'):
    def fail_write(fd, arrays):
        with open(fd, 'wb') as stream:
            stream.write(b'partial header')
        closed_fds.append(fd)
        raise error if scenario == 'write_interrupt' else OSError('disk full')
    unit.write_pileup_file = fail_write

try:
    if scenario in ('success', 'large_success', 'write_error'):
        assert rows(call(caller)) == expected
        if scenario == 'write_error':
            assert closed_fds
            assert not [p for p in root.rglob('*') if p.is_file() and p != sentinel]
        else:
            files = [p for p in root.rglob('*') if p.is_file() and p != sentinel]
            assert len(files) == 2
        if scenario == 'large_success':
            for path in root.rglob('*'):
                if path.is_file() and path != sentinel:
                    with path.open('ab') as stream:
                        stream.write(b'0' * (2 << 20))
        caller.destroy()
        clean()
    else:
        try:
            symbol = ('invalid' if scenario in
                      ('invalid_retry', 'cutoff_retry', 'bedgraph_retry') else 'q')
            call(caller, symbol)
        except BaseException as exc:
            if scenario == 'zero_lambda':
                assert isinstance(exc, AssertionError), repr(exc)
            elif scenario in ('shutdown_signal', 'real_shutdown_signal'):
                assert isinstance(exc, KeyboardInterrupt), repr(exc)
            elif scenario in ('invalid_retry', 'cutoff_retry', 'bedgraph_retry'):
                assert isinstance(exc, (IndexError, ValueError)), repr(exc)
            else:
                assert exc is error, repr(exc)
        else:
            raise AssertionError('expected injected failure')
        clean()  # failure cleanup is automatic, before destroy is called
        assert signal.getsignal(signal.SIGINT) is expected_handler
        if scenario == 'real_shutdown_signal':
            assert (root.parent / 'slow-worker-finished').read_text() == 'finished'
            assert not multiprocessing.active_children(), 'workers remain after API exit'
        if scenario in ('shutdown_signal', 'shutdown_custom_signal'):
            assert len(shutdowns) == 1
        if scenario == 'shutdown_custom_signal':
            assert signal_deliveries == [signal.SIGINT]
        if created_files:
            assert shutdowns or scenario == 'pool_construct'
            assert not any(p.exists() for p in created_files)
    caller.destroy()
    caller.destroy()
    clean()
    unit.write_pileup_file = real_write
    unit.ProcessPoolExecutor = real_executor
    unit.peak_processes_allowed = real_eligible
    unit.as_completed = real_completed
    unit.pscore_stat_worker = real_stat_worker
    if scenario != 'zero_lambda':
        assert rows(call(caller)) == expected
        for name, content in expected_outputs.items():
            assert (root.parent / name).read_bytes() == content
        caller.destroy()
        clean()
    for fd in closed_fds:
        try:
            os.fstat(fd)
        except OSError:
            pass
        else:
            raise AssertionError('write descriptor remains open')
finally:
    unit.write_pileup_file = real_write
    unit.ProcessPoolExecutor = real_executor
    signal.signal(signal.SIGINT, saved_handler)
    unit.pscore_stat_worker = real_stat_worker
    caller.destroy()
'''


@pytest.mark.skipif(os.name != 'posix', reason='bounded process-group probes require POSIX')
@pytest.mark.parametrize('mode', ['narrow', 'broad'])
@pytest.mark.parametrize('scenario', [
    'success', 'large_success', 'zero_lambda', 'invalid_retry',
    'cutoff_retry', 'bedgraph_retry',
    'write_error', 'write_interrupt',
    'pool_unreturned', 'pool_partial', 'pool_submit', 'pool_construct', 'interrupt',
    'shutdown_interrupt', 'shutdown_signal', 'shutdown_custom_signal',
])
def test_pileup_lifecycle(tmp_path, scenario, mode):
    run_probe(tmp_path, scenario, mode, 'serial')


def test_peak_worker_partial_failure(tmp_path):
    run_probe(tmp_path, 'pool_peak_partial', 'narrow', 'serial')


@pytest.mark.skipif(sys.platform != 'linux' or
                    'fork' not in multiprocessing.get_all_start_methods(),
                    reason='production process backend requires Linux fork')
@pytest.mark.parametrize('scenario', [
    'success', 'zero_lambda', 'invalid_retry', 'real_shutdown_signal'])
@pytest.mark.parametrize('mode', ['narrow', 'broad'])
def test_pileup_real_processes(tmp_path, scenario, mode):
    run_probe(tmp_path, scenario, mode, 'auto')


def run_probe(tmp_path, scenario, mode, backend):
    tempdir = tmp_path / 'temp'
    tempdir.mkdir()
    script = tmp_path / 'pileup_probe.py'
    script.write_text(SCRIPT)
    process = subprocess.Popen(
        [sys.executable, str(script), scenario, mode, backend],
        cwd=tmp_path, env=macs3_env(tempdir), start_new_session=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = process.communicate(timeout=30)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        pytest.fail('pileup probe timed out\n' + stdout + stderr)
    assert process.returncode == 0, stdout + stderr
