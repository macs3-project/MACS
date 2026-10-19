"""Exercise independent caller storage through compiled public APIs."""

import tempfile

import pytest

from MACS3.Signal.CallPeakUnit import CallerFromAlignments
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Utilities.Execution import CPUAvailability, ExecutionPolicy


def make_caller(kind, policy, offset):
    track = {'fw': FWTrack, 'pe': PETrackI, 'frag': PETrackII}[kind](threads=1)
    for chrom in (b'chr1', b'chr2'):
        for _ in range(25):
            if kind == 'fw':
                track.add_loc(chrom, 1000 + offset, 0)
            elif kind == 'pe':
                track.add_loc(chrom, 1000 + offset, 1200 + offset)
            else:
                track.add_loc(chrom, 1000 + offset, 1200 + offset, b'cell', 2)
    track.finalize()
    return CallerFromAlignments(
        track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[],
        lambda_bg=0.5, execution_policy=policy)


def call(caller, mode):
    if mode == 'broad':
        peaks = caller.call_broadpeaks(['q'], [3.0], [1.0], min_length=100)
    else:
        peaks = caller.call_peaks(['q'], [3.0], min_length=100,
                                  call_summits=mode == 'summits')
    keys = ['start', 'end', 'length', 'score', 'pileup', 'pscore', 'qscore',
            'fc', 'name']
    keys += (['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
             if mode == 'broad' else ['summit'])
    return [(chrom, [tuple(row[key] for key in keys) for row in rows])
            for chrom, rows in peaks.peaks.items()]


@pytest.mark.parametrize('kind', ['fw', 'pe', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
@pytest.mark.parametrize('budget', [1, 4])
def test_destroy_and_reuse_leave_other_live_caller_intact(tmp_path, monkeypatch,
                                                         kind, mode, budget):
    # Small fixtures keep each case below 0.1 s; storage is observed on disk,
    # without accessing caller internals or substituting scoring functions.
    monkeypatch.setattr(tempfile, 'tempdir', str(tmp_path))
    sentinel = tmp_path / 'unrelated.txt'
    sentinel.write_bytes(b'keep')
    policy = ExecutionPolicy(budget, budget, CPUAvailability(4, (), ()))
    first = make_caller(kind, policy, 0)
    second = make_caller(kind, policy, 5000)
    try:
        expected_first = call(first, mode)
        first_dir, = [p for p in tmp_path.iterdir() if p.is_dir()]
        expected_second = call(second, mode)
        second_dir, = [p for p in tmp_path.iterdir() if p.is_dir() and p != first_dir]
        assert expected_first and expected_second and expected_first != expected_second
        second_files = {p.name: p.read_bytes() for p in second_dir.iterdir()}
        assert len(second_files) == 2

        first.destroy()
        first.destroy()
        assert not first_dir.exists()
        assert {p.name: p.read_bytes() for p in second_dir.iterdir()} == second_files
        assert call(second, mode) == expected_second
        assert call(first, mode) == expected_first
        fresh_dir, = [p for p in tmp_path.iterdir() if p.is_dir() and p != second_dir]
        assert fresh_dir != first_dir

        # A missing cached chromosome file must be rebuilt by the same caller.
        cached_file = next(fresh_dir.iterdir())
        cached_file.unlink()
        assert call(first, mode) == expected_first
        assert len(list(fresh_dir.iterdir())) == 2
        second.destroy()
        assert not second_dir.exists()
        assert call(first, mode) == expected_first
    finally:
        first.destroy()
        second.destroy()
    assert list(tmp_path.iterdir()) == [sentinel]
    assert sentinel.read_bytes() == b'keep'
