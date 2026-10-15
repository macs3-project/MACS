"""Bounded subprocess checks of worker cleanup through public workflows."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = r'''
import sys
import threading
import time
from unittest.mock import patch

import numpy as np
from concurrent.futures import Future
import MACS3.Signal.CallPeakUnit as unit
from MACS3.Utilities import Threads
from MACS3.Signal.FixWidthTrack import FWTrack
import MACS3.Signal.PeakModel as model_module
from test_PeakModel import model_opt, site_track

workflow, failure_kind, fail_at = sys.argv[1], sys.argv[2], int(sys.argv[3])
original_thread = threading.Thread
original_start = threading.Thread.start
workers = []
finished = []
observed_arrays = []
entered = threading.Event()
release = threading.Event()
failure = RuntimeError('injected lifecycle failure')
created = 0
starts = 0

# The timer is outside the patched worker factory. It bounds the wait and
# lets a correct API join its workers without needing the caller to return.
timer = threading.Timer(0.15, release.set)
if failure_kind not in ('start', 'construct'):
    release.set()

def inject():
    assert entered.wait(5)
    original_start(timer)
    raise failure

def thread_factory(*args, **kwargs):
    global created
    created += 1
    if failure_kind == 'construct' and created == fail_at:
        inject()
    target = kwargs['target']
    target_args = kwargs.get('args', ())
    observed_arrays.extend(arg for arg in target_args if isinstance(arg, np.ndarray))
    def held_work():
        entered.set()
        assert release.wait(5)
        try:
            target(*target_args)
        finally:
            finished.append(threading.current_thread())
    thread = original_thread(target=held_work)
    workers.append(thread)
    return thread

def start(thread):
    global starts
    starts += 1
    if failure_kind == 'start' and starts == fail_at:
        inject()
    original_start(thread)

def run():
    if workflow == 'map':
        def task(chrom):
            if failure_kind == 'task':
                raise failure
            return chrom * 2
        return Threads.map_chromosomes(task, [1, 2, 3, 4], [1, 4, 2, 3])
    if workflow == 'scoring':
        # Supply a process-boundary histogram large enough to require four
        # scoring threads, without forking a host that has a timer thread.
        # The real compiled cache merge, scoring and peak call still run.
        keys = (np.uint64(10) << np.uint64(32)) | np.linspace(
            0.5, 2.0, 32769, dtype='f4').view('u4').astype('u8')
        bg_bits = np.array([0.5], dtype='f4').view('u4')[0]
        keys = np.concatenate([keys, np.array([
            (n << 32) | int(bg_bits) for n in (0, 1, 20)], dtype='u8')])
        firsts = np.arange(len(keys), dtype='i8')
        lengths = np.ones(len(keys), dtype='i8')
        class HistogramPool:
            def __init__(self, **kwargs):
                pass
            def submit(self, fn, batch):
                future = Future()
                future.set_result(([(sorted(track.get_chr_names())[i], None) for i in batch],
                                   keys, firsts, lengths))
                return future
            def shutdown(self, **kwargs):
                pass
        track = FWTrack(fw=36)
        for chrom in (b'a', b'b', b'c', b'd'):
            for pos in [1000] * 20 + [9000]:
                track.add_loc(chrom, pos, 0)
        track.finalize()
        caller = unit.CallerFromAlignments(
            track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[],
            lambda_bg=0.5, process_backend='auto')
        eligibility = iter([True, True, False])
        try:
            with patch.object(unit, 'peak_processes_allowed', lambda mode: next(eligibility)), \
                 patch.object(unit, 'n_available_workers', lambda: 4), \
                 patch.object(unit, 'ProcessPoolExecutor', HistogramPool):
                caller.call_peaks(['q'], [3.0], min_length=100)
                from MACS3.Signal.Prob import poisson_cdf
                pscores = observed_arrays[3]
                for i in (0, 8191, 16384, len(keys) - 1):
                    key = int(keys[i])
                    lam = np.array([key & 0xffffffff], dtype='u4').view('f4')[0]
                    np.testing.assert_allclose(pscores[i],
                        -poisson_cdf(key >> 32, float(lam), False, True),
                        rtol=1e-6)
        finally:
            caller.destroy()
        return
    if workflow == 'track':
        track = FWTrack(fw=36)
        for chrom in (b'a', b'b', b'c', b'd'):
            for pos in range(65536, 0, -1):
                track.add_loc(chrom, pos, 0)
        track.finalize()
        assert track.total == 4 * 65536
        for chrom in track.get_chr_names():
            np.testing.assert_array_equal(track.get_locations_by_chr(chrom)[0],
                                          np.arange(1, 65537))
        return
    opt, _ = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 120,
                          chroms=(b'a', b'b', b'c', b'd'))
    model = model_module.PeakModel(opt, track)
    if failure_kind in ('initializer', 'task'):
        allocation = 'zeros_like' if failure_kind == 'initializer' else 'empty'
        real_zeros = getattr(np, allocation)
        def zeros(*args, **kwargs):
            if threading.current_thread() is not threading.main_thread():
                raise failure
            return real_zeros(*args, **kwargs)
        with patch.object(np, allocation, zeros):
            model.build()
    else:
        model.build()
    return model.d

expected_d = None
if workflow == 'model':
    opt, _ = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 120,
                          chroms=(b'a', b'b', b'c', b'd'))
    baseline = model_module.PeakModel(opt, track)
    baseline.build()
    expected_d = baseline.d

try:
    with patch.object(Threads, 'n_threads', lambda: 4), \
         patch.object(model_module, '_model_threads', lambda: 4), \
         patch.object(threading, 'Thread', thread_factory), \
         patch.object(original_thread, 'start', start):
        try:
            result = run()
        except RuntimeError as error:
            assert failure_kind != 'normal'
            assert error is failure
            assert workers
            assert all(not t.is_alive() for t in workers), 'worker outlived API failure'
            assert len(finished) == sum(t.ident is not None for t in workers)
            snapshot = list(finished)
            arrays = [array.tobytes() for array in observed_arrays]
            time.sleep(0.02)
            assert arrays == [array.tobytes() for array in observed_arrays]
            assert finished == snapshot, 'worker still changing state after failure'
        else:
            assert failure_kind == 'normal', 'injected failure did not reach caller'
            assert all(not t.is_alive() for t in workers)
            if workflow == 'map':
                assert result == [2, 4, 6, 8]
            if workflow == 'model':
                assert result == expected_d
finally:
    release.set()
    if timer.ident is not None:
        timer.join(5)
    for worker in workers:
        if worker.ident is not None:
            worker.join(5)
'''


@pytest.mark.parametrize('workflow', ['map', 'track', 'model', 'scoring'])
@pytest.mark.parametrize('failure_kind,fail_at', [
    ('normal', 0), ('start', 2), ('start', 3), ('construct', 2),
])
def test_public_worker_lifecycle(tmp_path, workflow, failure_kind, fail_at):
    run_probe(tmp_path, workflow, failure_kind, fail_at)


@pytest.mark.parametrize('workflow,failure_kind', [
    ('map', 'task'), ('model', 'task'), ('model', 'initializer'),
])
def test_public_worker_exceptions(tmp_path, workflow, failure_kind):
    run_probe(tmp_path, workflow, failure_kind, 0)


def run_probe(tmp_path, workflow, failure_kind, fail_at):
    import MACS3

    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    # Follow pytest's selected installation, which may be outside the checkout.
    # Only expose the test directory for the shared model fixture helpers.
    package_root = Path(MACS3.__file__).resolve().parent.parent
    env['PYTHONPATH'] = str(package_root) + os.pathsep + str(root / 'test')
    result = subprocess.run(
        [sys.executable, '-c', SCRIPT, workflow, failure_kind, str(fail_at)],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
