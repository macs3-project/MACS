#!/usr/bin/env python
"""Module Description: Test functions in PileupV2 (weighted p-v pileups).

Every pileup is compared with a small numpy reference: the weighted
coverage of half-open intervals [start, end) accumulated on the
elementary segments between distinct coordinates, written as change
points (pos, value), where segment i covers [pos[i-1], pos[i]) with
pos[-1] read as 0, and the last segment ends at the largest end.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import numpy as np
import pytest

from MACS3.Signal.PileupV2 import (mapping_function_always_1,
                                   pileup_from_LR_hmmratac,
                                   pileup_from_LR,
                                   pileup_from_LRC,
                                   pileup_from_PN)

INT32_MAX = 2**31 - 1
INT32_MIN = -2**31
LR_DTYPE = [('l', 'i4'), ('r', 'i4')]
LRC_DTYPE = [('l', 'i4'), ('r', 'i4'), ('c', 'u2')]
PV_DTYPE = np.dtype([('p', 'i4'), ('v', 'f4')])


# ------------------------------------
# reference implementation and helpers
# ------------------------------------

def merge_runs(pos, vals):
    """Drop breakpoints whose value equals the value of the next segment."""
    pos = np.asarray(pos, dtype=np.int64)
    vals = np.asarray(vals, dtype=np.float32)
    if pos.size == 0:
        return pos, vals
    keep = np.ones(pos.size, dtype=bool)
    keep[:-1] = vals[:-1] != vals[1:]
    return pos[keep], vals[keep]


def ref_pileup(starts, ends, weights=None):
    """Reference weighted coverage of [start, end) intervals as change points.

    Coverage is a cumulative sum of +w at starts and -w at ends over the
    sorted distinct coordinates (a compressed coverage array, so
    coordinates near INT32_MAX need no large array).
    """
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if weights is None:
        weights = np.ones(starts.size)
    weights = np.asarray(weights, dtype=np.float64)
    keep = ends > starts
    starts, ends, weights = starts[keep], ends[keep], weights[keep]
    if starts.size == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.float32)
    assert starts.min() >= 0, "reference covers non-negative coordinates"
    coords = np.unique(np.concatenate(([0], starts, ends)))
    delta = np.zeros(coords.size)
    np.add.at(delta, np.searchsorted(coords, starts), weights)
    np.add.at(delta, np.searchsorted(coords, ends), -weights)
    cov = np.cumsum(delta)[:-1]
    return merge_runs(coords[1:], cov.astype(np.float32))


def assert_pv_equal(result, expected):
    """PileupV2 returns merged change points, so compare them raw."""
    assert result.dtype == PV_DTYPE
    np.testing.assert_array_equal(result['p'], expected[0])
    np.testing.assert_array_equal(result['v'], expected[1])


def make_lr(pairs):
    arr = np.zeros(len(pairs), dtype=LR_DTYPE)
    for i, (l, r) in enumerate(pairs):
        arr[i] = (l, r)
    return arr


def make_lrc(triples):
    arr = np.zeros(len(triples), dtype=LRC_DTYPE)
    for i, (l, r, c) in enumerate(triples):
        arr[i] = (l, r, c)
    return arr


def random_lr(seed, n, span=2000, maxlen=300):
    rng = np.random.default_rng(seed)
    arr = np.zeros(n, dtype=LR_DTYPE)
    arr['l'] = rng.integers(0, span, n)
    arr['r'] = arr['l'] + rng.integers(1, maxlen, n)
    return arr


def random_lrc(seed, n, maxcount=5, span=2000, maxlen=300):
    rng = np.random.default_rng(seed + 1000)
    lr = random_lr(seed, n, span, maxlen)
    arr = np.zeros(n, dtype=LRC_DTYPE)
    arr['l'] = lr['l']
    arr['r'] = lr['r']
    arr['c'] = rng.integers(1, maxcount + 1, n)
    return arr


def dyadic_weight(L, R):
    # multiples of 1/4 keep every partial sum exact in float32
    return 0.25 * ((R - L) % 7 + 1)


# ------------------------------------
# mapping_function_always_1
# ------------------------------------

@pytest.mark.parametrize("L, R", [(0, 0), (10, 5), (-5, 100),
                                  (INT32_MAX, INT32_MIN)])
def test_mapping_function_always_1_returns_one(L, R):
    w = mapping_function_always_1(L, R)
    assert isinstance(w, float)
    assert w == 1.0


@pytest.mark.parametrize("L, R", [(2**31, 0), (0, -2**31 - 1)])
def test_mapping_function_always_1_rejects_out_of_int32(L, R):
    with pytest.raises(OverflowError):
        mapping_function_always_1(L, R)


# ------------------------------------
# pileup_from_LR
# ------------------------------------

def test_pileup_from_LR_single_fragment():
    # [5, 10) -> 0 on [0, 5), 1 on [5, 10)
    assert_pv_equal(pileup_from_LR(make_lr([(5, 10)])),
                    ([5, 10], [0.0, 1.0]))


def test_pileup_from_LR_fragment_at_zero_has_no_leading_zero_segment():
    assert_pv_equal(pileup_from_LR(make_lr([(0, 10)])), ([10], [1.0]))


def test_pileup_from_LR_touching_fragments_merge():
    # [0, 10) and [10, 20) give one segment of value 1 over [0, 20)
    assert_pv_equal(pileup_from_LR(make_lr([(0, 10), (10, 20)])),
                    ([20], [1.0]))


def test_pileup_from_LR_gap_is_a_zero_segment():
    # [0, 10) value 1, [10, 30) value 0, [30, 40) value 1
    assert_pv_equal(pileup_from_LR(make_lr([(0, 10), (30, 40)])),
                    ([10, 30, 40], [1.0, 0.0, 1.0]))


def test_pileup_from_LR_duplicates_add_up():
    assert_pv_equal(pileup_from_LR(make_lr([(5, 10)] * 3)),
                    ([5, 10], [0.0, 3.0]))


def test_pileup_from_LR_empty():
    res = pileup_from_LR(np.zeros(0, dtype=LR_DTYPE))
    assert res.dtype == PV_DTYPE
    assert res.shape == (0,)


def test_pileup_from_LR_int32_extreme_positions():
    lr = make_lr([(INT32_MAX - 10, INT32_MAX)])
    assert_pv_equal(pileup_from_LR(lr),
                    ([INT32_MAX - 10, INT32_MAX], [0.0, 1.0]))


@pytest.mark.parametrize("seed, n", [(0, 1), (1, 7), (2, 50), (3, 300),
                                     (4, 1000)])
def test_pileup_from_LR_matches_reference(seed, n):
    lr = random_lr(seed, n)
    assert_pv_equal(pileup_from_LR(lr), ref_pileup(lr['l'], lr['r']))


def test_pileup_from_LR_input_order_does_not_matter():
    lr = random_lr(11, 200)
    shuffled = lr.copy()
    np.random.default_rng(5).shuffle(shuffled)
    a = pileup_from_LR(lr)
    b = pileup_from_LR(shuffled)
    np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_pileup_from_LR_mapping_function_weights(seed):
    lr = random_lr(seed, 150)
    w = [dyadic_weight(int(l), int(r)) for l, r in lr]
    assert_pv_equal(pileup_from_LR(lr, dyadic_weight),
                    ref_pileup(lr['l'], lr['r'], w))


def test_pileup_from_LR_mapping_function_gets_left_and_right():
    seen = []

    def record(L, R):
        seen.append((L, R))
        return 1.0
    pileup_from_LR(make_lr([(3, 9), (20, 25)]), record)
    # called once for the start and once for the end of each fragment
    assert sorted(set(seen)) == [(3, 9), (20, 25)]


def test_pileup_from_LR_mapping_function_error_propagates():
    def broken(L, R):
        raise ValueError("no weight")
    with pytest.raises(ValueError, match="no weight"):
        pileup_from_LR(make_lr([(3, 9)]), broken)


# ------------------------------------
# pileup_from_LRC
# ------------------------------------

@pytest.mark.parametrize("count", [1, 2, 7, 65535])
def test_pileup_from_LRC_count_is_the_weight(count):
    assert_pv_equal(pileup_from_LRC(make_lrc([(5, 10, count)])),
                    ([5, 10], [0.0, float(count)]))


def test_pileup_from_LRC_zero_count_adds_no_coverage():
    # [0, 10) with count 0: value 0 everywhere up to 10
    assert_pv_equal(pileup_from_LRC(make_lrc([(5, 10, 0)])), ([10], [0.0]))


def test_pileup_from_LRC_large_counts_stack_exactly():
    # 3 x 65535 = 196605 is exact in float32
    lrc = make_lrc([(0, 10, 65535), (2, 8, 65535), (4, 6, 65535)])
    assert_pv_equal(pileup_from_LRC(lrc),
                    ([2, 4, 6, 8, 10],
                     [65535.0, 131070.0, 196605.0, 131070.0, 65535.0]))


def test_pileup_from_LRC_empty():
    res = pileup_from_LRC(np.zeros(0, dtype=LRC_DTYPE))
    assert res.dtype == PV_DTYPE
    assert res.shape == (0,)


@pytest.mark.parametrize("seed, n, maxcount", [(0, 5, 3), (1, 100, 5),
                                               (2, 500, 2), (3, 50, 65535)])
def test_pileup_from_LRC_matches_reference(seed, n, maxcount):
    lrc = random_lrc(seed, n, maxcount)
    assert_pv_equal(pileup_from_LRC(lrc),
                    ref_pileup(lrc['l'], lrc['r'], lrc['c']))


def test_pileup_from_LRC_mapping_function_multiplies_count():
    lrc = random_lrc(9, 120)
    w = [int(c) * dyadic_weight(int(l), int(r)) for l, r, c in lrc]
    assert_pv_equal(pileup_from_LRC(lrc, dyadic_weight),
                    ref_pileup(lrc['l'], lrc['r'], w))


def test_pileup_from_LRC_int32_extreme_positions():
    lrc = make_lrc([(INT32_MAX - 50, INT32_MAX, 2)])
    assert_pv_equal(pileup_from_LRC(lrc),
                    ([INT32_MAX - 50, INT32_MAX], [0.0, 2.0]))


# ------------------------------------
# pileup_from_LR_hmmratac
# ------------------------------------

def test_pileup_from_LR_hmmratac_hand_example():
    # lengths 10 -> 0.5, 20 -> 2.0; [0, 10) 0.5 and [5, 25) 2.0:
    # [0, 5) 0.5, [5, 10) 2.5, [10, 25) 2.0
    lr = make_lr([(0, 10), (5, 25)])
    assert_pv_equal(pileup_from_LR_hmmratac(lr, {10: 0.5, 20: 2.0}),
                    ([5, 10, 25], [0.5, 2.5, 2.0]))


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_pileup_from_LR_hmmratac_matches_reference(seed):
    lr = random_lr(seed, 200)
    lengths = sorted(set((lr['r'] - lr['l']).tolist()))
    mapping = {L: 0.25 * (L % 5) for L in lengths}   # includes weight 0
    w = [mapping[int(r - l)] for l, r in lr]
    assert_pv_equal(pileup_from_LR_hmmratac(lr, mapping),
                    ref_pileup(lr['l'], lr['r'], w))


def test_pileup_from_LR_hmmratac_missing_length_raises_keyerror():
    with pytest.raises(KeyError):
        pileup_from_LR_hmmratac(make_lr([(0, 10)]), {20: 1.0})


def test_pileup_from_LR_hmmratac_requires_dict():
    with pytest.raises(TypeError):
        pileup_from_LR_hmmratac(make_lr([(0, 10)]), [1.0] * 20)


def test_pileup_from_LR_hmmratac_empty():
    res = pileup_from_LR_hmmratac(np.zeros(0, dtype=LR_DTYPE), {})
    assert res.dtype == PV_DTYPE
    assert res.shape == (0,)


# ------------------------------------
# pileup_from_PN
# ------------------------------------

def test_pileup_from_PN_hand_example():
    # plus 10 -> [10, 30), plus 50 -> [50, 70),
    # minus 30 -> [10, 30), minus 100 -> [80, 100)
    P = np.array([10, 50], dtype="i4")
    N = np.array([30, 100], dtype="i4")
    assert_pv_equal(pileup_from_PN(P, N, 20),
                    ([10, 30, 50, 70, 80, 100],
                     [0.0, 2.0, 0.0, 1.0, 0.0, 1.0]))


@pytest.mark.parametrize("seed, n, extsize", [(0, 1, 5), (1, 20, 50),
                                              (2, 200, 147), (3, 500, 1)])
def test_pileup_from_PN_matches_reference(seed, n, extsize):
    rng = np.random.default_rng(seed)
    P = np.sort(rng.integers(0, 3000, n)).astype("i4")
    N = np.sort(rng.integers(extsize, 3000, n)).astype("i4")
    starts = np.concatenate((P, N - extsize))
    ends = np.concatenate((P + extsize, N))
    assert_pv_equal(pileup_from_PN(P, N, extsize), ref_pileup(starts, ends))


def test_pileup_from_PN_empty():
    e = np.zeros(0, dtype="i4")
    res = pileup_from_PN(e, e, 10)
    assert res.dtype == PV_DTYPE
    assert res.shape == (0,)


def test_pileup_from_PN_unequal_strand_counts():
    """Regression test: pileup_from_PN asserted len(P) == len(N), so FWTrack
    strands with unequal read counts raised AssertionError.

    Fixed upstream in e0bd1e2 (#733), which rewrote pileup_from_PN to
    concatenate the plus and minus starts/ends instead of pairing them in
    make_PV_from_PN.
    """
    P = np.array([10, 20, 40], dtype="i4")
    N = np.array([60], dtype="i4")
    # plus: [10, 30), [20, 40), [40, 60); minus 60 -> [40, 60)
    expected = ref_pileup([10, 20, 40, 40], [30, 40, 60, 60])
    assert_pv_equal(pileup_from_PN(P, N, 20), expected)


