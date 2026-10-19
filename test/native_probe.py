"""Isolated drivers for public native compatibility/failure regressions."""

import faulthandler
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import zlib

import numpy as np

faulthandler.enable()


class FragmentReceiver:
    """Duck-typed public append receiver selecting the scalar parsing path."""

    def __init__(self, track):
        self.track = track

    def __getattr__(self, name):
        return getattr(self.track, name)


def streams(manifest, cleanup=None):
    from MACS3.IO.Parser import BAMParser, BAMPEParser, FragParser
    from MACS3.Signal.FixWidthTrack import FWTrack
    from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
    if cleanup and sys.platform in ('darwin', 'linux'):
        from execution_probe import native_thread_count
        before = native_thread_count()
    results = []
    for case in json.loads(Path(manifest).read_text()):
        path, kind, operation, budget, cap = case
        parser = None
        track = None
        try:
            parser = {'bam': BAMParser, 'bampe': BAMPEParser, 'frag': FragParser}[kind](
                path, threads=budget, buffer_size=2)
            if operation in ('append', 'fallback'):
                track = {'bam': FWTrack, 'bampe': PETrackI, 'frag': PETrackII}[kind](
                    threads=budget, buffer_size=2)
                if operation == 'fallback':
                    track = FragmentReceiver(track)
            kwargs = {'max_count': cap} if kind == 'frag' else {}
            if kind == 'bam':
                track = (parser.build_fwtrack() if track is None else
                         parser.append_fwtrack(track))
            else:
                track = (parser.build_petrack(**kwargs) if track is None else
                         parser.append_petrack(track, **kwargs))
            track.finalize()
            arrays = []
            small_records = {}
            for chrom in sorted(track.get_chr_names()):
                locs = track.get_locations_by_chr(chrom)
                for array in locs if kind == 'bam' else [locs]:
                    arrays.append(hashlib.sha256(array.tobytes()).hexdigest())
                if kind == 'frag':
                    arrays.append(hashlib.sha256(track.barcodes[chrom].tobytes()).hexdigest())
                    if len(locs) <= 16:
                        small_records[chrom.hex()] = locs.tolist()
            barcodes = ([(hashlib.sha256(bar).hexdigest(), code)
                         for bar, code in track.barcode_dict.items()]
                        if kind == 'frag' else [])
            results.append(['ok', track.total, arrays, barcodes, small_records])
        except (gzip.BadGzipFile, EOFError, OverflowError, ValueError, zlib.error) as exc:
            results.append([type(exc).__name__, str(exc)])
        finally:
            if parser is not None:
                parser.close()
            if track is not None:
                track.destroy()
            if cleanup and sys.platform in ('darwin', 'linux'):
                assert native_thread_count() <= before, 'stream retained native workers'
    print(json.dumps(results))


def allocator(expected):
    # NumPy's diagnostic reports each array's retained NEP 49 handler. This
    # does not read MACS3 internals or call its allocator functions directly.
    try:
        from numpy._core.multiarray import get_handler_name
    except ImportError:
        from numpy.core.multiarray import get_handler_name
    from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
    original = get_handler_name()
    observations = []
    for cls in (PETrackI, PETrackII):
        track = cls(buffer_size=17, threads=1)
        handler = expected if cls is PETrackI else original
        chrom = b'chr1'
        if cls is PETrackII:
            track.barcode_id(b'BC')
        # Incremental bulk appends cross the 2 MiB mapping threshold, then
        # reserve growth and finalize shrink it. Under 150 MiB live per track.
        offset = 0
        for size in (16, 270000, 1200000):
            starts = np.arange(offset, offset + size, dtype='i4')
            if cls is PETrackI:
                track.add_loc_arrays(chrom, starts, starts + 40)
            else:
                track.add_loc_arrays(chrom, starts, starts + 40,
                                     np.ones(size, dtype='u2'), np.zeros(size, dtype='i4'))
            offset += size
            array = track.get_locations_by_chr(chrom)
            assert get_handler_name(array) == handler, get_handler_name(array)
            np.testing.assert_array_equal(array['l'][:offset], np.arange(offset))
            observations.append([cls.__name__, offset, get_handler_name(array)])
            del array
            assert get_handler_name() == original, 'track allocation leaked handler'
        track.finalize()
        assert track.total == offset
        np.testing.assert_array_equal(track.get_locations_by_chr(chrom)['l'], np.arange(offset))
        copy = track.sample_percent_copy(0.25, seed=770)
        assert copy.total == int(offset * 0.25)
        copy.destroy()
        if cls is PETrackII:
            subset = track.subset({b'BC'})
            assert subset.total == offset
            subset.destroy()
        track.destroy()
        track.destroy()
        assert not track.get_chr_names()
    assert get_handler_name() == original
    print(json.dumps(observations))


def process_start_failure(eligibility, phase, fail_at):
    import multiprocessing
    from multiprocessing.process import BaseProcess
    from unittest.mock import patch
    import MACS3.Signal.CallPeakUnit as unit
    from MACS3.Signal.FixWidthTrack import FWTrack
    from MACS3.Utilities.Execution import CPUAvailability, ExecutionPolicy
    from MACS3.Utilities.Performance import DEFAULTS
    def rows(peaks):
        keys = ('start', 'end', 'length', 'summit', 'score', 'pileup',
                'pscore', 'qscore', 'fc', 'name')
        return [(chrom, [tuple(row[key] for key in keys) for row in values])
                for chrom, values in peaks.peaks.items()]
    if eligibility == 'controlled':
        # Exercise real fork/start/join under a controlled host-safety decision.
        # Production macOS fallback has separate backend tests.
        unit.peak_processes_allowed = lambda backend: backend == 'auto'
    track = FWTrack(fw=36, threads=1)
    for chrom in (b'chr1', b'chr2', b'chr3'):
        track.add_loc_arrays(chrom, np.full(DEFAULTS.min_parallel_peak_locations // 3 + 1,
                                          1000, dtype='i4'), np.array([], dtype='i4'))
    track.finalize()
    reference = unit.CallerFromAlignments(
        track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[], lambda_bg=0.5,
        execution_policy=ExecutionPolicy(1, 1, CPUAvailability(4, (), ())))
    expected = rows(reference.call_peaks(['q'], [2.0], min_length=50))
    reference.destroy()
    calc = unit.CallerFromAlignments(
        track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[], lambda_bg=0.5,
        execution_policy=ExecutionPolicy(4, 4, CPUAvailability(4, (), ())))
    before = set(Path(os.environ['TMPDIR']).iterdir())
    if phase == 'calling':
        with patch.object(unit, 'peak_processes_allowed', lambda backend: False):
            assert calc.call_peaks(['q'], [2.0], min_length=50).total > 0
    real_start = BaseProcess.start
    started = []
    def start(process):
        if len(started) == int(fail_at) - 1:
            raise OSError('injected process startup failure')
        real_start(process)
        started.append(process)
    try:
        with patch.object(BaseProcess, 'start', start):
            try:
                calc.call_peaks(['q'], [2.0], min_length=50)
            except OSError as exc:
                assert str(exc) == 'injected process startup failure'
            else:
                raise AssertionError('process failure path was not reached')
        assert len(started) == int(fail_at) - 1
        assert all(not process.is_alive() for process in started)
        assert not multiprocessing.active_children()
        assert set(Path(os.environ['TMPDIR']).iterdir()) == before
        with patch.object(unit, 'peak_processes_allowed', lambda backend: False):
            actual = rows(calc.call_peaks(['q'], [2.0], min_length=50))
            assert actual == expected, (actual, expected)
    finally:
        calc.destroy()
        track.destroy()
    assert set(Path(os.environ['TMPDIR']).iterdir()) == before
    print(json.dumps({'started': len(started), 'joined': True}))


if __name__ == '__main__':
    {'streams': streams, 'allocator': allocator,
     'process_start_failure': process_start_failure}[sys.argv[1]](*sys.argv[2:])
