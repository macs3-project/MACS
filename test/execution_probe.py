"""Bounded subprocess driver exercising compiled public execution paths."""

import ctypes
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import threading
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from MACS3.IO.Parser import BAMParser, BAMPEParser, FragParser
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Signal.PeakModel import PeakModel
from MACS3.Signal.CallPeakUnit import CallerFromAlignments, close_unlinked_files_later
from MACS3.Utilities.Execution import resolve_execution_policy


def snapshot(track):
    track.finalize()
    values = {}
    for chrom in sorted(track.get_chr_names()):
        locations = track.get_locations_by_chr(chrom)
        arrays = locations if isinstance(locations, (list, tuple)) else [locations]
        values[chrom.decode()] = [hashlib.sha256(a.tobytes()).hexdigest() for a in arrays]
        if isinstance(track, PETrackII):
            values[chrom.decode()].append(hashlib.sha256(track.barcodes[chrom].tobytes()).hexdigest())
    return [track.total, values]


def native_thread_count():
    """Count live worker threads, excluding joined OS teardown remnants."""
    if sys.platform == 'linux':
        count = 0
        # Linux sets PF_EXITING before clearing the TID that wakes
        # pthread_join; /proc may retain that task until release_task.
        # Such a task cannot run the worker callback again. Do not count
        # its stale directory as a leaked worker or wait for it to vanish.
        # PF_EXITING is 0x4 in include/linux/sched.h; stat field 9 is flags.
        for task in Path('/proc/self/task').iterdir():
            try:
                stat = (task / 'stat').read_bytes()
            except (FileNotFoundError, ProcessLookupError):
                continue  # Task was removed after directory enumeration.
            # comm (field 2) can contain spaces and closing parentheses.
            flags = int(stat.rsplit(b')', 1)[1].split()[6])
            if not flags & 0x4:
                count += 1
        return count
    if sys.platform != 'darwin':
        raise RuntimeError('native thread observation requires Linux or macOS')
    lib = ctypes.CDLL(None)
    lib.mach_task_self.restype = ctypes.c_uint
    task = lib.mach_task_self()
    array = ctypes.POINTER(ctypes.c_uint)()
    count = ctypes.c_uint()
    lib.task_threads.argtypes = [ctypes.c_uint, ctypes.POINTER(ctypes.POINTER(ctypes.c_uint)),
                                ctypes.POINTER(ctypes.c_uint)]
    # task_threads may still list Mach threads after pthread_join returns.
    # The pthread list has already dropped these joined threads. Native BGZF
    # workers use pthread_create, so live workers (even waiting ones) still map.
    lib.pthread_from_mach_thread_np.argtypes = [ctypes.c_uint]
    lib.pthread_from_mach_thread_np.restype = ctypes.c_void_p
    assert lib.task_threads(task, ctypes.byref(array), ctypes.byref(count)) == 0
    try:
        return sum(bool(lib.pthread_from_mach_thread_np(array[i]))
                   for i in range(count.value))
    finally:
        for i in range(count.value):
            lib.mach_port_deallocate(task, array[i])
        lib.vm_deallocate.argtypes = [ctypes.c_uint, ctypes.c_size_t, ctypes.c_size_t]
        lib.vm_deallocate(task, ctypes.cast(array, ctypes.c_void_p).value,
                          count.value * ctypes.sizeof(ctypes.c_uint))


def run(scenario, budget, args):
    policy = resolve_execution_policy(budget)
    # A later environment change must not affect objects made with this policy.
    os.environ['MACS3_NUM_THREADS'] = 'bad unused environment setting'
    started = []
    real_start = threading.Thread.start

    def start(thread):
        started.append(thread)
        return real_start(thread)

    with patch.object(threading.Thread, 'start', start):
        if scenario == 'parser':
            kind, operation, path = args
            parser = {'frag': FragParser, 'bam': BAMParser, 'bampe': BAMPEParser}[kind](
                path, buffer_size=4096, execution_policy=policy)
            assert parser.execution_policy is policy
            pe = kind != 'bam'
            if operation == 'append':
                # Appending does not overwrite the receiving track's policy.
                receiving = resolve_execution_policy(1)
                track_type = PETrackII if kind == 'frag' else PETrackI if pe else FWTrack
                track = track_type(execution_policy=receiving)
                result = parser.append_petrack(track) if pe else parser.append_fwtrack(track)
                assert result is track and track.execution_policy is receiving
            else:
                track = parser.build_petrack() if pe else parser.build_fwtrack()
                assert track.execution_policy is policy
            result = snapshot(track)
            if kind == 'frag':
                subset = track.subset({b'BC1'})
                assert subset.execution_policy is track.execution_policy
                result += [dict((k.decode(), v) for k, v in track.barcode_dict.items()),
                           snapshot(subset)]
            parser.close()
        elif scenario == 'native':
            ready, stop = threading.Event(), threading.Event()
            counts = []
            def observe():
                counts.append(native_thread_count())
                ready.set()
                while not stop.wait(0.0005):
                    counts.append(native_thread_count())
            observer = threading.Thread(target=observe)
            real_start(observer)
            assert ready.wait(5)
            before = counts[0]
            stream_samples = []
            real_open = io.open

            def observed_open(filename, *open_args, **kwargs):
                # The public file-open boundary is reached after native pool
                # construction and before consumption. Sample synchronously
                # here so a fast parse cannot hide every worker between polls.
                # No private stream/pool symbol is patched or invoked.
                if kwargs.get('buffering') == 0 and os.fspath(filename) == args[0]:
                    count = native_thread_count()
                    counts.append(count)
                    stream_samples.append(count)
                return real_open(filename, *open_args, **kwargs)

            try:
                with patch.object(io, 'open', observed_open):
                    parser = FragParser(args[0], execution_policy=policy)
                    track = parser.build_petrack()
                    parser.close()
                assert parser.n == 180000
            finally:
                stop.set()
                observer.join(5)
            delta = max(counts) - before
            assert stream_samples, 'native stream file-open boundary was not observed'
            assert delta <= policy.threads - 1, (counts, policy)
            if policy.threads == 1:
                assert delta == 0, counts
            else:
                assert delta > 0, 'the fixture did not exercise native workers'
            result = {'effective': policy.threads, 'extra_native_threads': delta,
                      'stream_open_samples': len(stream_samples)}
        elif scenario == 'track':
            result = []
            for track_type in (FWTrack, PETrackI, PETrackII):
                track = track_type(execution_policy=policy)
                for chrom in (b'a', b'b', b'c'):
                    for pos in range(65537, 0, -1):
                        if track_type is FWTrack:
                            track.add_loc(chrom, pos // 2, 0)
                        elif track_type is PETrackI:
                            track.add_loc(chrom, pos // 2, pos // 2 + 40)
                        else:
                            track.add_loc(chrom, pos // 2, pos // 2 + 40, b'A', 2)
                result.append(snapshot(track))
                if track_type is not PETrackII:
                    track.filter_dup(1)
                    result.append(snapshot(track))
                else:
                    assert track.subset({b'A'}).execution_policy is policy
                if track_type is PETrackI:
                    copied = track.sample_percent_copy(0.5, seed=17)
                    assert copied.execution_policy is policy
                    result.append(snapshot(copied))
                if track_type is PETrackII:
                    copied = track.sample_percent_copy(0.5, seed=17)
                    assert copied.execution_policy is policy
                    result.append(snapshot(copied))
        elif scenario == 'model':
            track = FWTrack(fw=36, execution_policy=policy)
            for i in range(120):
                chrom = (b'a', b'b', b'c')[i % 3]
                center = 2000 + (i // 3) * 5000
                for _ in range(5):
                    track.add_loc(chrom, center - 50, 0)
                    track.add_loc(chrom, center + 50, 1)
            track.finalize()
            noop = lambda *a, **kw: None
            opt = SimpleNamespace(gsize=1e6, lmfold=5, umfold=200, bw=100,
                                  d_min=20, info=noop, debug=noop, warn=noop)
            model = PeakModel(opt, track, execution_policy=policy)
            assert model.execution_policy is policy
            model.build()
            result = [model.d, model.alternative_d,
                      [hashlib.sha256(a.tobytes()).hexdigest()
                       for a in (model.plus_line, model.minus_line, model.ycorr)]]
        elif scenario == 'caller':
            mode, control, symbol = args
            # Caller settings are independent of the supplied track's budget.
            track = FWTrack(fw=36, threads=1)
            ctrl = FWTrack(fw=36, threads=1) if control == 'yes' else None
            for i in range(3):
                chrom = ('chr%d' % i).encode()
                for pos in [1000 + i * 200] * (20 + i * 7) + [9000]:
                    track.add_loc(chrom, pos, 0)
                if ctrl is not None:
                    for pos in (1000, 1200, 9000):
                        ctrl.add_loc(chrom, pos, 0)
            track.finalize()
            if ctrl is not None:
                ctrl.finalize()
            caller = CallerFromAlignments(
                track, ctrl, d=200, ctrl_d_s=[200] if ctrl is not None else [],
                ctrl_scaling_factor_s=[1.0] if ctrl is not None else [],
                lambda_bg=0.5, execution_policy=policy)
            assert caller.execution_policy is policy
            try:
                if mode == 'broad':
                    peaks = caller.call_broadpeaks([symbol], [3.0], [1.0], min_length=100)
                else:
                    peaks = caller.call_peaks([symbol], [3.0], min_length=100,
                                              call_summits=mode == 'summits')
                assert peaks.peaks and all(peaks.peaks.values())
                keys = ['start', 'end', 'length', 'score', 'pileup',
                        'pscore', 'qscore', 'fc', 'name']
                keys += (['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
                         if mode == 'broad' else ['summit'])
                def rows(value):
                    return [(chrom.decode(), [[row[key] for key in keys] for row in data])
                            for chrom, data in sorted(value.peaks.items())]
                result = rows(peaks)
                if mode != 'broad':
                    again = caller.call_peaks([symbol], [3.0], min_length=100,
                                              call_summits=mode == 'summits')
                    assert result == rows(again)
            finally:
                caller.destroy()
        elif scenario == 'cleanup':
            fds = []
            for i in range(8):
                path = 'unlinked-%d' % i
                fds.append(os.open(path, os.O_WRONLY | os.O_CREAT, 0o600))
                os.unlink(path)
            close_unlinked_files_later(fds, [1 << 20] * len(fds), policy)
            for fd in fds:
                try:
                    os.fstat(fd)
                except OSError:
                    pass
                else:
                    raise AssertionError('descriptor was not closed')
            result = 'closed'
        else:
            raise AssertionError(scenario)
    if policy.threads == 1:
        assert not started, [t.name for t in started]
    assert all(not t.is_alive() for t in started), 'worker outlived API return'
    return result


if __name__ == '__main__':
    print(json.dumps(run(sys.argv[1], int(sys.argv[2]), sys.argv[3:]), sort_keys=True,
                     default=lambda value: value.decode() if isinstance(value, bytes) else value.item()))
