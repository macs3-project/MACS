"""Public track contracts sensitive to NumPy sorting and native allocation."""

import numpy as np
import pytest

from MACS3.Signal.PairedEndTrack import PETrackII
from MACS3.Utilities.Execution import CPUAvailability, ExecutionPolicy


@pytest.mark.parametrize('size', [1, 2, 15, 16, 17, 63, 64, 65, 1001,
                                      65535, 65536, 65537])
@pytest.mark.parametrize('budget', [1, 4])
@pytest.mark.parametrize('operation', ['finalize', 'sort'])
def test_tied_fragments_follow_installed_numpy(size, budget, operation):
    assert_tied_order(size, budget, operation, 'random')


@pytest.mark.parametrize('pattern', ['sorted', 'reverse', 'all-tied'])
@pytest.mark.parametrize('size', [16, 17, 257, 4097])
@pytest.mark.parametrize('operation', ['finalize', 'sort'])
def test_adversarial_tied_order(size, pattern, operation):
    assert_tied_order(size, 1, operation, pattern)


def assert_tied_order(size, budget, operation, pattern):
    # Two chromosomes reach the thread threshold. Normally <1 s per case.
    policy = ExecutionPolicy(budget, budget, CPUAvailability(4, (), ()))
    track = PETrackII(buffer_size=size, execution_policy=policy)
    rng = np.random.default_rng(770)
    expected = {}
    for chrom in (b'chr1', b'chr2'):
        for index in range(size):
            if pattern == 'random':
                left = int(rng.integers(0, 8))
                right = left + 40
                count = (1, 2, 65535)[int(rng.integers(0, 3))]
            elif pattern == 'all-tied':
                left, right, count = 10, 50, 2
            else:
                rank = index if pattern == 'sorted' else size - index - 1
                left = rank // 8
                right = left + 40 + (rank // 4) % 2
                count = 1 + (rank // 2) % 2
            # Equal coordinates, different counts, and exact l/r/count ties.
            barcode = str(index).encode() + (b'\xff', b'\x00', b'B' * 37)[index % 3]
            track.add_loc(chrom, left, right, barcode, count)
        locs = track.get_locations_by_chr(chrom).copy()
        bars = track.barcodes[chrom].copy()
        order = np.argsort(locs, order=['l', 'r'])
        expected[chrom] = (locs[order], bars[order])
    getattr(track, operation)()
    for chrom, (locs, bars) in expected.items():
        np.testing.assert_array_equal(track.get_locations_by_chr(chrom), locs)
        np.testing.assert_array_equal(track.barcodes[chrom], bars)
    assert track.is_sorted
    # A subset must retain barcode association even among indistinguishable rows.
    selected = next(iter(track.barcode_dict))
    code = track.barcode_dict[selected]
    subset = track.subset({selected})
    for chrom, (locs, bars) in expected.items():
        np.testing.assert_array_equal(subset.get_locations_by_chr(chrom), locs[bars == code])
        np.testing.assert_array_equal(subset.barcodes[chrom], bars[bars == code])


@pytest.mark.parametrize('populated', [False, True])
@pytest.mark.parametrize('finalized', [False, True])
def test_frag_destroy_releases_all_chromosomes(populated, finalized):
    track = PETrackII(buffer_size=2)
    if populated:
        for chrom in (b'chr1', b'chr2'):
            for i in range(9):
                track.add_loc(chrom, i, i + 40, b'BC', 1)
        if finalized:
            track.finalize()
    for _ in range(2):
        track.destroy()
        assert track.locations == track.barcodes == track.barcode_dict == {}
        assert track.get_chr_names() == set()
