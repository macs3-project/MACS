"""Bounded subprocess driver for compiled public chromosome scheduling calls."""

from concurrent.futures import Future
import faulthandler
import json
import logging
from pathlib import Path
import sys

import numpy as np
import MACS3.Signal.CallPeakUnit as unit
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.IO.PeakIO import BroadPeakIO
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Utilities.Execution import CPUAvailability, ExecutionPolicy
from MACS3.Utilities.Performance import DEFAULTS

faulthandler.enable()
logging.disable(logging.CRITICAL)
assert Path(unit.__file__).suffix in ('.so', '.pyd'), unit.__file__


def policy(budget, backend='auto'):
    return ExecutionPolicy(budget, budget, CPUAvailability(64, (), ()),
                           requested_process_backend=backend)


def make_track(kind, counts, multiplicity=1, prefix='chr'):
    """Zero counts retain empty chromosomes via public sampling."""
    track = {'fw': FWTrack, 'pe': PETrackI, 'frag': PETrackII}[kind](
        buffer_size=128, threads=1)
    sampled = 0 in counts
    if kind == 'frag':
        track.barcode_id(b'A')
    for i, count in enumerate(counts):
        size = max(1, count * 10) if sampled else count
        chrom = ('%s%05d' % (prefix, i)).encode()
        # Distinct hills and a distant singleton make output/order observable.
        starts = np.full(size, 1000 + (i % 3) * 400, dtype='i4')
        if size > 10:
            starts[-1] = 9000
        if kind == 'fw':
            track.add_loc_arrays(chrom, starts, np.array([], dtype='i4'))
        elif kind == 'pe':
            track.add_loc_arrays(chrom, starts, starts + 200)
        else:
            track.add_loc_arrays(chrom, starts, starts + 200,
                                 np.full(size, multiplicity, dtype='u2'),
                                 np.zeros(size, dtype='i4'))
    if counts:
        track.finalize()
    if sampled:
        track.sample_percent(0.1, seed=11)
    return track


def without_empty(track, kind):
    out = {'fw': FWTrack, 'pe': PETrackI, 'frag': PETrackII}[kind](threads=1)
    if kind == 'frag':
        out.barcode_id(b'A')
    for chrom in sorted(track.get_chr_names()):
        if unit.n_locations(track, chrom) == 0:
            continue
        locs = track.get_locations_by_chr(chrom)
        if kind == 'fw':
            out.add_loc_arrays(chrom, *locs)
        elif kind == 'pe':
            out.add_loc_arrays(chrom, locs['l'], locs['r'])
        else:
            out.add_loc_arrays(chrom, locs['l'], locs['r'], locs['c'],
                               np.zeros(len(locs), dtype='i4'))
    if out.get_chr_names():
        out.finalize()
    return out


def call(caller, mode, symbol='q', cutoff=False):
    if mode == 'broad':
        return caller.call_broadpeaks([symbol], [2.0], [0.5], min_length=50,
                                      cutoff_analysis=cutoff)
    return caller.call_peaks([symbol], [2.0], min_length=50,
                             call_summits=mode == 'summits', cutoff_analysis=cutoff)


def rows(peaks):
    # Retain returned chromosome/peak order and every public peak value.
    keys = ['start', 'end', 'length', 'score', 'pileup', 'pscore', 'qscore', 'fc', 'name']
    keys += (['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
             if isinstance(peaks, BroadPeakIO) else ['summit'])
    return [(chrom, [tuple(row[key] for key in keys) for row in values])
            for chrom, values in peaks.peaks.items()]


def caller(track, ctrl=None, budget=4, backend='auto', artifacts=False):
    return unit.CallerFromAlignments(
        track, ctrl, d=200, ctrl_d_s=[200, 1000, 10000] if ctrl is not None else [],
        ctrl_scaling_factor_s=[0.1, 0.02, 0.002] if ctrl is not None else [],
        lambda_bg=0.5, execution_policy=policy(budget, backend),
        save_bedGraph=artifacts, bedGraph_treat_filename='treat.bdg',
        bedGraph_control_filename='ctrl.bdg', cutoff_analysis_filename='cutoff.tsv')


def schedule(config):
    kind, workload, phase = (config[k] for k in ('kind', 'workload', 'phase'))
    minimum = DEFAULTS.min_parallel_peak_locations
    counts = {
        'tiny': [20, 30, 40], 'many_tiny': [1] * 447,
        'below': [minimum // 2, minimum - minimum // 2 - 1],
        'at': [minimum // 2, minimum - minimum // 2],
        'rounding': [minimum // 224] * 224 + [minimum % 224],
        'above': [minimum // 2, minimum - minimum // 2 + 1],
        'skewed': [minimum, 1000] + [1] * 1000,
        'mixed': [minimum // 2, minimum // 2] + [0] * 2000,
        'single': [minimum * 2],
    }[workload]
    track = make_track(kind, counts, config.get('multiplicity', 1))
    ctrl = make_track(kind, counts) if config.get('control') else None
    names = sorted(track.get_chr_names())
    weights = [unit.n_locations(track, c) +
               (unit.n_locations(ctrl, c) if ctrl is not None else 0) for c in names]
    positive = {i for i, w in enumerate(weights) if w}
    total = sum(weights)
    expected = (total >= minimum and len(positive) > 1 and
                config.get('budget', 8) > 1 and config.get('backend', 'auto') == 'auto')
    pools = []
    class Intercepted(Exception):
        pass
    class Recorder:
        def __init__(self, **kwargs):
            self.workers = kwargs['max_workers']
            self.batches = []
            pools.append(self)
        def submit(self, function, task):
            assert function is (unit.pscore_stat_worker if phase == 'histogram'
                                else unit.call_peaks_worker)
            indices = task if phase == 'histogram' else task[0]
            self.batches.append(list(indices))
            future = Future()
            future.set_exception(Intercepted())
            return future
        def shutdown(self, **kwargs):
            assert kwargs == {'wait': True, 'cancel_futures': True}
    checks = 0
    def eligible(backend):
        nonlocal checks
        checks += 1
        return backend == 'auto' and (phase == 'histogram' or checks > 1)
    unit.peak_processes_allowed = eligible
    unit.ProcessPoolExecutor = Recorder
    calc = caller(track, ctrl, config.get('budget', 8), config.get('backend', 'auto'))
    try:
        try:
            call(calc, 'narrow')
        except Intercepted:
            assert expected
        else:
            assert not expected
    finally:
        calc.destroy()
    assert len(pools) == int(expected), (pools, expected)
    if pools:
        pool = pools[0]
        assert 1 <= pool.workers <= config.get('budget', 8) - 1
        limit = DEFAULTS.chromosome_tasks_per_worker * min(
            config.get('budget', 8) - 1, len(positive))
        assert 1 <= len(pool.batches) <= limit
        flat = [i for batch in pool.batches for i in batch]
        assert len(flat) == len(set(flat)) and set(flat) == positive
        target = (total + limit - 1) // limit
        assert all(batch and sum(weights[i] for i in batch) > 0 for batch in pool.batches)
        assert all(len(batch) == 1 for batch in pool.batches
                   if any(weights[i] >= target for i in batch))
        assert all(sum(weights[i] for i in batch) >= target for batch in pool.batches[:-1])
        if workload in ('skewed', 'rounding'):
            assert max(map(len, pool.batches)) > 1
        return {'tasks': len(pool.batches), 'locations': total, 'workers': pool.workers}
    return {'tasks': 0, 'locations': total}


def outputs(config):
    kind, mode, workload = (config[k] for k in ('kind', 'mode', 'workload'))
    n = DEFAULTS.min_parallel_peak_locations
    counts = {'empty': [], 'disjoint': [20, 30], 'retained': [0] * 17,
              'mixed': [40, 50, 0, 0, 60], 'tiny': [20, 30, 40],
              'normal': [n // 2, n // 3, n // 4],
              'skewed': [n, 50] + [1] * 100}[workload]
    track = make_track(kind, counts)
    ctrl = make_track(kind, [max(1, c // 10) if c else 0 for c in counts]) if config.get('control') else None
    if workload == 'disjoint':
        ctrl = make_track(kind, counts, prefix='other')
    # No mocking of host eligibility: Linux uses real eligible workers.
    real_executor = unit.ProcessPoolExecutor
    pools = []
    completions = []
    real_completed = unit.as_completed
    def recording_executor(**kwargs):
        pools.append(kwargs['max_workers'])
        return real_executor(**kwargs)
    def reverse_completed(futures):
        # A deterministic completion order different from task submission.
        done = list(real_completed(futures))
        completions.append(len(done))
        return iter(reversed(list(futures)))
    unit.ProcessPoolExecutor = recording_executor
    unit.as_completed = reverse_completed
    expected = None
    artifacts = config.get('artifacts', False)
    for budget, backend in ((1, 'serial'), (4, 'auto')):
        calc = caller(track, ctrl, budget, backend, artifacts)
        try:
            result = call(calc, mode, config.get('symbol', 'q'), artifacts)
            snapshot = rows(result)
            files = {p: Path(p).read_bytes() for p in ('treat.bdg', 'ctrl.bdg', 'cutoff.tsv')} if artifacts else {}
            if expected is None:
                expected = (snapshot, files)
            else:
                assert (snapshot, files) == expected
            again = call(calc, mode, config.get('symbol', 'q'), artifacts)
            assert rows(again) == snapshot
            calc.destroy()
            assert rows(call(calc, mode, config.get('symbol', 'q'), artifacts)) == snapshot
        finally:
            calc.destroy()
    if workload in ('mixed', 'retained'):
        clean = without_empty(track, kind)
        clean_ctrl = without_empty(ctrl, kind) if ctrl is not None else None
        calc = caller(clean, clean_ctrl, 1, 'serial')
        try:
            assert rows(call(calc, mode, config.get('symbol', 'q'))) == expected[0]
        finally:
            calc.destroy()
    if workload in ('empty', 'retained', 'disjoint'):
        assert expected[0] == [] and not pools
        if artifacts:
            assert files['treat.bdg'] == files['ctrl.bdg'] == b''
            assert files['cutoff.tsv'] == b'pscore\tqscore\tnpeaks\tlpeaks\tavelpeak\n'
    elif workload in ('tiny', 'normal', 'mixed'):
        assert expected[0], 'fixture must call peaks'
    if unit.peak_processes_allowed('auto') and workload in ('normal', 'skewed') and not artifacts:
        assert pools and completions, 'fixture must exercise real Linux workers'
    return {'chromosomes_with_peaks': len(expected[0]), 'pools': len(pools)}


def one_sided(config):
    kind, mode = config['kind'], config['mode']
    track = make_track(kind, [50, 0])
    ctrl = make_track(kind, [0, 50])
    expected = None
    for budget in (1, 4):
        calc = caller(track, ctrl, budget, artifacts=True)
        try:
            peaks = call(calc, mode)
            result = rows(peaks)
            files = [Path(p).read_bytes() for p in ('treat.bdg', 'ctrl.bdg')]
            if expected is None:
                expected = result, files
            else:
                assert (result, files) == expected
            assert set(peaks.peaks) == {b'chr00000'}
            # The treatment-empty chromosome still contributes its control span.
            treatment_rows = [line.split(b'\t') for line in files[0].splitlines()]
            control_rows = [line.split(b'\t') for line in files[1].splitlines()]
            assert any(row[0] == b'chr00001' for row in treatment_rows)
            assert all(float(row[3]) == 0 for row in treatment_rows if row[0] == b'chr00001')
            assert all(float(row[3]) >= 0.5 for row in control_rows)
            assert any(float(row[3]) == 0.5 for row in control_rows if row[0] == b'chr00000')
        finally:
            calc.destroy()
    return {'one_sided': 'ok'}


if __name__ == '__main__':
    config = json.loads(sys.argv[1])
    print(json.dumps({'schedule': schedule, 'outputs': outputs, 'one_sided': one_sided}[config['action']](config)))
