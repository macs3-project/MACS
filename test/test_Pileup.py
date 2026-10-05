#!/usr/bin/env python
"""Module Description: Test functions for pileup functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest
import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
from math import log2
from MACS3.Signal.Pileup import (se_all_in_one_pileup,
                                 quick_pileup,
                                 naive_quick_pileup,
                                 over_two_pv_array,
                                 naive_call_peaks as naive_call_peaks_v1,
                                 pileup_and_write_se as pileup_and_write_se_v1)
from MACS3.Signal.PileupV2 import (pileup_from_LR,
                                   pileup_from_LRC,
                                   pileup_from_PN,
                                   pileup_from_LR_as_list,
                                   pileup_from_PN_shifted,
                                   over_two_pv_array as over_two_pv_array_v2,
                                   naive_quick_pileup as naive_quick_pileup_v2,
                                   naive_call_peaks as naive_call_peaks_v2,
                                   pileup_and_write_pe as pileup_and_write_pe_v2,
                                   pileup_and_write_se as pileup_and_write_se_v2)
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.IO.BedGraphIO import bedGraphIO

import pytest
from MACS3.Signal.Pileup import (pileup_and_write_se,
                                 pileup_and_write_pe,
                                 naive_call_peaks)
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI

# ------------------------------------
# Main function
# ------------------------------------


class Test_SE_Pileup(unittest.TestCase):
    """Unittest for pileup functions in Pileup.pyx for single-end
    datasets.

    Function to test: se_all_in_one_pileup

    """
    def setUp(self):
        self.plus_pos = np.array((0, 1, 3), dtype="int32")
        self.minus_pos = np.array((8, 9, 10), dtype="int32")
        self.rlength = 100      # right end of coordinates
        # expected result from pileup_bdg_se: (start, end, value)
        # the actual fragment length is 1+five_shift+three_shift
        self.param_1 = {"five_shift": 0,
                        "three_shift": 5,
                        "scale_factor": 0.5,
                        "baseline": 0}
        self.expect_pileup_1 = [(0, 1, 0.5),
                                (1, 3, 1.0),
                                (3, 4, 2.0),
                                (4, 6, 2.5),
                                (6, 8, 2.0),
                                (8, 9, 1.0),
                                (9, 10, 0.5)]
        # expected result from pileup_w_multiple_d_bdg_se: (start, end, value)
        self.param_2 = {"five_shift": 0,
                        "three_shift": 10,
                        "scale_factor": 2,
                        "baseline": 8}
        self.expect_pileup_2 = [(0,  1, 8.0),
                                (1,  3, 10.0),
                                (3,  8, 12.0),
                                (8,  9, 10.0),
                                (9,  10, 8.0),
                                (10, 11, 8.0),
                                (11, 13, 8.0)]

    def test_pileup_1(self):
        pileup = se_all_in_one_pileup(self.plus_pos, self.minus_pos,
                                      self.param_1["five_shift"],
                                      self.param_1["three_shift"],
                                      self.rlength,
                                      self.param_1["scale_factor"],
                                      self.param_1["baseline"])
        result = []
        (p,v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pileup_1)

    def test_pileup_2(self):
        pileup = se_all_in_one_pileup(self.plus_pos, self.minus_pos,
                                      self.param_2["five_shift"],
                                      self.param_2["three_shift"],
                                      self.rlength,
                                      self.param_2["scale_factor"],
                                      self.param_2["baseline"])
        result = []
        (p, v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pileup_2)


class Test_Quick_Pileup(unittest.TestCase):
    """Unittest for pileup functions in Pileup.pyx for quick-pileup.

    Function to test: quick_pileup

    """
    def setUp(self):
        self.start_pos = np.array((0, 1, 3, 3, 4, 5), dtype="int32")
        self.end_pos = np.array((5, 6, 8, 8, 9, 10), dtype="int32")
        # expected result from pileup_bdg_se: (start, end, value)
        self.param_1 = {"scale_factor": 0.5,
                        "baseline": 0}
        self.expect_pileup_1 = [(0, 1, 0.5),
                                (1, 3, 1.0),
                                (3, 4, 2.0),
                                (4, 6, 2.5),
                                (6, 8, 2.0),
                                (8, 9, 1.0),
                                (9, 10, 0.5)]

    def test_pileup_1(self):
        pileup = quick_pileup(self.start_pos, self.end_pos,
                              self.param_1["scale_factor"],
                              self.param_1["baseline"])
        result = []
        (p, v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pileup_1)


class Test_Naive_Pileup(unittest.TestCase):
    """Unittest for pileup functions in Pileup.pyx for naive-quick-pileup.

    Function to test: naive_quick_pileup

    """
    def setUp(self):
        self.pos = np.array((2, 3, 5, 5, 6, 7), dtype="int32")
        # expected result from pileup_bdg_se: (start, end, value)
        self.param_1 = {"extension": 2}
        self.expect_pileup_1 = [(0, 1, 1.0),
                                (1, 3, 2.0),
                                (3, 7, 4.0),
                                (7, 8, 2.0),
                                (8, 9, 1.0)]

    def test_pileup_1(self):
        pileup = naive_quick_pileup(self.pos,
                                    self.param_1["extension"])
        result = []
        (p, v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pileup_1)


class Test_Over_Two_PV_Array(unittest.TestCase):
    """Unittest for over_two_pv_array function

    Function to test: over_two_pv_array

    """
    def setUp(self):
        self.pv1 = [np.array((2, 5, 7, 8, 9, 12), dtype="int32"),
                    np.array((1, 2, 3, 4, 3, 2), dtype="float32")]
        self.pv2 = [np.array((1, 4, 6, 8, 10, 11), dtype="int32"),
                    np.array((5, 3, 2, 1, 0, 3), dtype="float32")]
        # expected result from pileup_bdg_se: (start, end, value)
        self.expect_pv_max = [(0, 1, 5.0),
                              (1, 2, 3.0),
                              (2, 4, 3.0),
                              (4, 5, 2.0),
                              (5, 6, 3.0),
                              (6, 7, 3.0),
                              (7, 8, 4.0),
                              (8, 9, 3.0),
                              (9, 10, 2.0),
                              (10, 11, 3.0)]
        self.expect_pv_min = [(0,  1, 1.0),
                              (1, 2, 1.0),
                              (2, 4, 2.0),
                              (4, 5, 2.0),
                              (5, 6, 2.0),
                              (6, 7, 1.0),
                              (7, 8, 1.0),
                              (8, 9, 0.0),
                              (9, 10, 0.0),
                              (10, 11, 2.0)]
        self.expect_pv_mean = [(0, 1, 3.0),
                               (1, 2, 2.0),
                               (2, 4, 2.5),
                               (4, 5, 2.0),
                               (5, 6, 2.5),
                               (6, 7, 2.0),
                               (7, 8, 2.5),
                               (8, 9, 1.5),
                               (9, 10, 1.0),
                               (10, 11, 2.5)]

    def test_max(self):
        pileup = over_two_pv_array(self.pv1, self.pv2, func="max")
        result = []
        (p, v) = pileup
        # print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pv_max)

    def test_min(self):
        pileup = over_two_pv_array(self.pv1, self.pv2, func="min")
        result = []
        (p, v) = pileup
        # print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pv_min)

    def test_mean(self):
        pileup = over_two_pv_array(self.pv1, self.pv2, func="mean")
        result = []
        (p, v) = pileup
        # print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append((pre, pos, value))
            pre = pos
        # check result
        self.assertEqual(result, self.expect_pv_mean)


class Test_PileupV2_PE(unittest.TestCase):
    """Unittest for pileup functions in PileupV2.pyx.

    Function to test: pileup_from_LR

    """
    def setUp(self):
        self.LR_array1 = np.array([(1, 5), (2, 6),
                                   (4, 8), (4, 8),
                                   (5, 9), (6, 10),
                                   (12, 14), (13, 17),
                                   (14, 18), (17, 19)],
                                  dtype=[('l', 'int32'), ('r', 'int32')])
        # expected result from pileup_from_LR: (end, value)
        self.expect_pileup_1 = np.array([(1, 0.0),
                                         (2, 1.0),
                                         (4, 2.0),
                                         (8, 4.0),
                                         (9, 2.0),
                                         (10, 1.0),
                                         (12, 0.0),
                                         (13, 1.0),
                                         (18, 2.0),
                                         (19, 1.0)],
                                        dtype=[('p', 'uint32'), ('v', 'float32')])
        # with log2(length) as weight
        self.expect_pileup_2 = np.array([(1, 0.0),
                                         (2, 2.0),
                                         (4, 4.0),
                                         (8, 8.0),
                                         (9, 4.0),
                                         (10, 2.0),
                                         (12, 0.0),
                                         (13, 1.0),
                                         (14, 3.0),
                                         (17, 4.0),
                                         (18, 3.0),
                                         (19, 1.0)],
                                        dtype=[('p', 'uint32'), ('v', 'float32')])

    def test_pileup_1(self):
        pileup = pileup_from_LR(self.LR_array1)
        np.testing.assert_equal(pileup, self.expect_pileup_1)

    def test_pileup_2(self):
        pileup = pileup_from_LR(self.LR_array1, lambda x, y: log2(y-x))
        np.testing.assert_equal(pileup, self.expect_pileup_2)


class Test_PileupV2_LRC(unittest.TestCase):
    """Unittest for count-weighted PileupV2 pileup."""

    def test_count_weighted_unsorted_right_endpoints(self):
        lrc = np.array([(5, 9, 2),
                        (1, 4, 1),
                        (2, 6, 3),
                        (2, 6, 2)],
                       dtype=[('l', 'int32'), ('r', 'int32'), ('c', 'uint16')])
        expect = np.array([(1, 0.0),
                           (2, 1.0),
                           (4, 6.0),
                           (5, 5.0),
                           (6, 7.0),
                           (9, 2.0)],
                          dtype=[('p', 'uint32'), ('v', 'float32')])

        pileup = pileup_from_LRC(lrc)
        np.testing.assert_equal(pileup, expect)


class Test_PileupV2_SE(unittest.TestCase):
    """Unittest for pileup functions in PileupV2.pyx.

    Function to test: pileup_from_PN

    """
    def setUp(self):
        self.P = np.array((0, 1, 3, 3, 4, 5), dtype="int32")   #plus strand pos
        self.N = np.array((5, 6, 8, 8, 9, 10), dtype="int32")  #minus strand pos
        # expected result from pileup_bdg_se: (start, end, value)
        self.extsize = 2
        self.expect_pileup_1 = np.array([(1, 1.0),
                                         (2, 2.0),
                                         (3, 1.0),
                                         (4, 3.0),
                                         (5, 5.0),
                                         (8, 3.0),
                                         (9, 2.0),
                                         (10, 1.0)],
                                        dtype=[('p', 'uint32'), ('v', 'float32')])

    def test_pileup_1(self):
        pileup = pileup_from_PN(self.P, self.N, self.extsize)
        np.testing.assert_equal(pileup, self.expect_pileup_1)


class Test_PileupV2_Compatibility_Helpers(unittest.TestCase):
    """V2 wrappers that preserve legacy [p, v] helper behavior."""

    def test_lr_as_list_matches_quick_pileup(self):
        lr = np.array([(3, 7), (1, 4), (6, 8), (6, 10)],
                      dtype=[("l", "i4"), ("r", "i4")])
        expected = quick_pileup(np.sort(lr["l"]), np.sort(lr["r"]), 1.25, 0.5)
        observed = pileup_from_LR_as_list(lr, 1.25, 0.5)
        np.testing.assert_array_equal(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])

    def test_shifted_pn_matches_se_all_in_one(self):
        plus = np.array([0, 3, 3, 8], dtype="i4")
        minus = np.array([5, 6, 9, 10], dtype="i4")
        expected = se_all_in_one_pileup(plus, minus, 2, 6, 20, 0.75, 0.25)
        observed = pileup_from_PN_shifted(plus, minus, 2, 6, 20, 0.75, 0.25)
        np.testing.assert_array_equal(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])

    def test_over_two_pv_array_matches_v1(self):
        pv1 = [np.array([1, 3, 5, 8], dtype="i4"),
               np.array([0, 2, 1, 3], dtype="f4")]
        pv2 = [np.array([2, 4, 5, 9], dtype="i4"),
               np.array([1, 1, 4, 0], dtype="f4")]
        for func in ("max", "min", "mean"):
            expected = over_two_pv_array(pv1, pv2, func=func)
            observed = over_two_pv_array_v2(pv1, pv2, func=func)
            np.testing.assert_array_equal(observed[0], expected[0])
            np.testing.assert_allclose(observed[1], expected[1])

    def test_naive_helpers_match_v1(self):
        tags = np.array([1, 3, 5, 8, 13, 21], dtype="i4")
        expected = naive_quick_pileup(tags, 3)
        observed = naive_quick_pileup_v2(tags, 3)
        np.testing.assert_array_equal(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])
        self.assertEqual(naive_call_peaks_v2(observed, 1, min_length=1),
                         naive_call_peaks_v1(expected, 1, min_length=1))


class Test_PileupV2_Writers(unittest.TestCase):
    """Unittest for direct V2 bedGraph writer entrypoints."""

    def test_se_writer_matches_v1(self):
        track = FWTrack()
        for pos in (0, 1, 3, 3, 4, 5):
            track.add_loc(b"chr1", pos, 0)
        for pos in (5, 6, 8, 8, 9, 10):
            track.add_loc(b"chr1", pos, 1)
        track.finalize()
        track.set_rlengths({b"chr1": 100})

        with TemporaryDirectory() as tmpdir:
            v1_path = Path(tmpdir) / "v1.bdg"
            v2_path = Path(tmpdir) / "v2.bdg"
            pileup_and_write_se_v1(track, str(v1_path).encode(), 2, 1,
                                   directional=True, halfextension=False)
            pileup_and_write_se_v2(track, str(v2_path).encode(), 2, 1,
                                   directional=True, halfextension=False)

            self.assertEqual(v2_path.read_text(), v1_path.read_text())

    def test_pe_writer_matches_track_bedgraph(self):
        track = PETrackI()
        for start, end in ((1, 5), (2, 6), (4, 8), (4, 8), (5, 9)):
            track.add_loc(b"chr1", start, end)
        track.finalize()
        track.set_rlengths({b"chr1": 100})

        with TemporaryDirectory() as tmpdir:
            expected_path = Path(tmpdir) / "expected.bdg"
            v2_path = Path(tmpdir) / "v2.bdg"
            bdg = track.pileup_bdg()
            bedGraphIO(str(expected_path), data=bdg).write_bedGraph(trackline=False)
            pileup_and_write_pe_v2(track, str(v2_path).encode())

            self.assertEqual(v2_path.read_text(), expected_path.read_text())


# ------------------------------------
# Reference pileup used by the tests below
# ------------------------------------
#
# The reference is the coverage of half-open intervals [start, end)
# accumulated on the elementary segments between the sorted distinct
# coordinates (a compressed coverage array, so INT32_MAX coordinates
# need no large array). Values are max(coverage * scale, baseline) in
# float32, the type MACS3 stores. Results are compared as change points
# (pos, value): segment i covers [pos[i-1], pos[i]) with pos[-1] read
# as 0, consecutive equal values merged.

INT32_MAX = 2**31 - 1


def merge_runs(pos, vals):
    """Drop breakpoints whose value equals the value of the next segment."""
    pos = np.asarray(pos, dtype=np.int64)
    vals = np.asarray(vals, dtype=np.float32)
    if pos.size == 0:
        return pos, vals
    keep = np.ones(pos.size, dtype=bool)
    keep[:-1] = vals[:-1] != vals[1:]
    return pos[keep], vals[keep]


def ref_pileup(starts, ends, scale=1.0, baseline=0.0):
    """Reference coverage of [start, end) intervals as change points."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    keep = ends > starts
    starts, ends = starts[keep], ends[keep]
    if starts.size == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.float32)
    assert starts.min() >= 0, "reference covers non-negative coordinates"
    coords = np.unique(np.concatenate(([0], starts, ends)))
    delta = np.zeros(coords.size, dtype=np.int64)
    np.add.at(delta, np.searchsorted(coords, starts), 1)
    np.add.at(delta, np.searchsorted(coords, ends), -1)
    cov = np.cumsum(delta)[:-1]
    vals = np.maximum(cov.astype(np.float32) * np.float32(scale),
                      np.float32(baseline))
    return merge_runs(coords[1:], vals)


def se_intervals(plus, minus, five_shift, three_shift, rlength):
    """Fragments from 5' ends: a plus tag at p covers [p - five_shift,
    p + three_shift), a minus tag at m covers [m - three_shift,
    m + five_shift); both ends are clipped to [0, rlength]."""
    plus = np.asarray(plus, dtype=np.int64)
    minus = np.asarray(minus, dtype=np.int64)
    starts = np.concatenate((plus - five_shift, minus - three_shift))
    ends = np.concatenate((plus + three_shift, minus + five_shift))
    return np.clip(starts, 0, rlength), np.clip(ends, 0, rlength)


def assert_cp_equal(result, expected):
    """Check a [p, v] result: int32/float32, strictly increasing
    positions, and the same change points as the reference."""
    p, v = result
    assert p.dtype == np.int32 and v.dtype == np.float32
    assert len(p) == len(v)
    assert np.all(np.diff(p.astype(np.int64)) > 0)
    gp, gv = merge_runs(p, v)
    np.testing.assert_array_equal(gp, expected[0])
    np.testing.assert_array_equal(gv, expected[1])


def bdg_text(chrom, cp):
    """bedGraph lines for change points, starting at 0, '%.5f' values."""
    out = []
    pre = 0
    for p, v in zip(cp[0].tolist(), cp[1].tolist()):
        out.append("%s\t%d\t%d\t%.5f\n" % (chrom, pre, p, v))
        pre = p
    return "".join(out)


def i4(x):
    return np.array(x, dtype="i4")


def f4(x):
    return np.array(x, dtype="f4")


# ------------------------------------
# se_all_in_one_pileup
# ------------------------------------

def test_se_all_in_one_pileup_hand_example():
    # plus 10 -> [10, 35), minus 30 -> [5, 30):
    # [0, 5) 0, [5, 10) 1, [10, 30) 2, [30, 35) 1
    p, v = se_all_in_one_pileup(i4([10]), i4([30]), 0, 25, INT32_MAX,
                                1.0, 0.0)
    np.testing.assert_array_equal(p, [5, 10, 30, 35])
    np.testing.assert_array_equal(v, [0.0, 1.0, 2.0, 1.0])


def test_se_all_in_one_pileup_extension_past_zero_is_clipped():
    # minus 20 -> [-10, 20) clipped to [0, 20); plus 0 -> [0, 30)
    p, v = se_all_in_one_pileup(i4([0]), i4([20]), 0, 30, INT32_MAX,
                                1.0, 0.0)
    np.testing.assert_array_equal(p, [20, 30])
    np.testing.assert_array_equal(v, [2.0, 1.0])


def test_se_all_in_one_pileup_clipped_at_rlength():
    # plus 90 -> [90, 120) clipped to [90, 100)
    p, v = se_all_in_one_pileup(i4([90]), i4([]), 0, 30, 100, 1.0, 0.0)
    np.testing.assert_array_equal(p, [90, 100])
    np.testing.assert_array_equal(v, [0.0, 1.0])


def test_se_all_in_one_pileup_duplicates_stack():
    p, v = se_all_in_one_pileup(i4([10, 10, 10]), i4([]), 0, 5, INT32_MAX,
                                1.0, 0.0)
    np.testing.assert_array_equal(p, [10, 15])
    np.testing.assert_array_equal(v, [0.0, 3.0])


def test_se_all_in_one_pileup_baseline_above_all_values():
    # coverage 0 then 1 both lift to the baseline 4
    res = se_all_in_one_pileup(i4([10]), i4([]), 0, 5, INT32_MAX, 1.0, 4.0)
    assert_cp_equal(res, ([15], [4.0]))


@pytest.mark.parametrize("plus, minus", [([INT32_MAX - 100], []),
                                         ([], [INT32_MAX])])
def test_se_all_in_one_pileup_int32_extreme(plus, minus):
    # plus at MAX-100 -> [MAX-100, MAX); minus at MAX -> [MAX-100, MAX)
    p, v = se_all_in_one_pileup(i4(plus), i4(minus), 0, 100, INT32_MAX,
                                1.0, 0.0)
    np.testing.assert_array_equal(p, [INT32_MAX - 100, INT32_MAX])
    np.testing.assert_array_equal(v, [0.0, 1.0])


def test_se_all_in_one_pileup_empty():
    p, v = se_all_in_one_pileup(i4([]), i4([]), 0, 100, INT32_MAX, 1.0, 0.0)
    assert p.dtype == np.int32 and v.dtype == np.float32
    assert p.shape == (0,) and v.shape == (0,)


SE_CASES = [
    # seed, n_plus, n_minus, five_shift, three_shift, rlength, scale, baseline
    (0, 5, 5, 0, 100, INT32_MAX, 1.0, 0.0),
    (1, 50, 40, 0, 147, INT32_MAX, 1.0, 0.0),
    (2, 100, 120, 73, 74, INT32_MAX, 0.5, 0.0),
    (3, 80, 0, -25, 75, INT32_MAX, 1.0, 0.0),
    (4, 0, 80, 0, 200, INT32_MAX, 2.0, 0.0),
    (5, 200, 200, 0, 300, 1500, 0.1, 0.0),
    (6, 60, 60, 50, 50, INT32_MAX, 1.0, 1.5),
    (7, 300, 300, 0, 100, 1800, 0.333, 0.25),
    (8, 30, 30, 1000, 1000, 2000, 1.0, 0.0),
]


@pytest.mark.parametrize("seed, n_plus, n_minus, five, three, rlength, "
                         "scale, baseline", SE_CASES)
def test_se_all_in_one_pileup_matches_reference(seed, n_plus, n_minus, five,
                                                three, rlength, scale,
                                                baseline):
    rng = np.random.default_rng(seed)
    plus = np.sort(rng.integers(0, 2000, n_plus)).astype("i4")
    minus = np.sort(rng.integers(0, 2000, n_minus)).astype("i4")
    res = se_all_in_one_pileup(plus, minus, five, three, rlength, scale,
                               baseline)
    starts, ends = se_intervals(plus, minus, five, three, rlength)
    assert_cp_equal(res, ref_pileup(starts, ends, scale, baseline))


def test_se_all_in_one_pileup_does_not_modify_input():
    plus = i4([5, 50, 90])
    minus = i4([10, 60])
    se_all_in_one_pileup(plus, minus, 20, 30, 70, 1.0, 0.0)
    np.testing.assert_array_equal(plus, [5, 50, 90])
    np.testing.assert_array_equal(minus, [10, 60])


# ------------------------------------
# quick_pileup
# ------------------------------------

def test_quick_pileup_touching_fragments_merge():
    # [0, 10) and [10, 20): one segment of value 1 over [0, 20)
    p, v = quick_pileup(i4([0, 10]), i4([10, 20]), 1.0, 0.0)
    np.testing.assert_array_equal(p, [20])
    np.testing.assert_array_equal(v, [1.0])


def test_quick_pileup_empty():
    p, v = quick_pileup(i4([]), i4([]), 1.0, 0.0)
    assert p.dtype == np.int32 and v.dtype == np.float32
    assert p.shape == (0,) and v.shape == (0,)


def test_quick_pileup_int32_extreme():
    p, v = quick_pileup(i4([INT32_MAX - 5]), i4([INT32_MAX]), 1.0, 0.0)
    np.testing.assert_array_equal(p, [INT32_MAX - 5, INT32_MAX])
    np.testing.assert_array_equal(v, [0.0, 1.0])


@pytest.mark.parametrize("seed, n, scale, baseline", [
    (0, 1, 1.0, 0.0), (1, 10, 1.0, 0.0), (2, 200, 0.5, 0.0),
    (3, 500, 1.0, 2.0), (4, 300, 0.07, 0.03), (5, 50, 3.0, 0.0)])
def test_quick_pileup_matches_reference(seed, n, scale, baseline):
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, 3000, n)
    ends = starts + rng.integers(1, 400, n)
    res = quick_pileup(np.sort(starts).astype("i4"),
                       np.sort(ends).astype("i4"), scale, baseline)
    assert_cp_equal(res, ref_pileup(starts, ends, scale, baseline))


def test_quick_pileup_duplicate_fragments():
    # sorted starts [5, 5, 5, 8], sorted ends [9, 20, 20, 20]:
    # [5, 8) 3, [8, 9) 4, [9, 20) 3
    res = quick_pileup(i4([5, 5, 5, 8]), i4([9, 20, 20, 20]), 1.0, 0.0)
    assert_cp_equal(res, ([5, 8, 9, 20], [0.0, 3.0, 4.0, 3.0]))


# ------------------------------------
# naive_quick_pileup
# ------------------------------------

def test_naive_quick_pileup_clips_at_zero():
    # 0 -> [-5, 5) clipped [0, 5); 3 -> [0, 8)
    res = naive_quick_pileup(i4([0, 3]), 5)
    assert_cp_equal(res, ([5, 8], [2.0, 1.0]))


@pytest.mark.parametrize("seed, n, extension", [(0, 1, 10), (1, 30, 1),
                                                (2, 200, 50), (3, 400, 200)])
def test_naive_quick_pileup_matches_reference(seed, n, extension):
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.integers(0, 3000, n)).astype("i4")
    res = naive_quick_pileup(pos, extension)
    starts = np.maximum(pos.astype(np.int64) - extension, 0)
    ends = pos.astype(np.int64) + extension
    assert_cp_equal(res, ref_pileup(starts, ends))


def test_naive_quick_pileup_empty_raises():
    with pytest.raises(Exception, match="length is 0"):
        naive_quick_pileup(i4([]), 5)


# ------------------------------------
# over_two_pv_array
# ------------------------------------

def ref_over_two(pv1, pv2, func):
    """Reference merge on the union of breakpoints up to the shorter end
    (the documented 'overlap regions')."""
    p1, v1 = np.asarray(pv1[0], np.int64), np.asarray(pv1[1], np.float32)
    p2, v2 = np.asarray(pv2[0], np.int64), np.asarray(pv2[1], np.float32)
    if p1.size == 0 or p2.size == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.float32)
    end = min(p1[-1], p2[-1])
    bps = np.union1d(p1, p2)
    bps = bps[bps <= end]
    a = v1[np.searchsorted(p1, bps)]
    b = v2[np.searchsorted(p2, bps)]
    if func == "max":
        out = np.maximum(a, b)
    elif func == "min":
        out = np.minimum(a, b)
    else:
        out = ((a.astype(np.float64) + b.astype(np.float64)) / 2)
    return bps, out.astype(np.float32)


def random_pv(rng, n):
    p = np.cumsum(rng.integers(1, 10, n)).astype("i4")
    v = (rng.random(n) * 10).astype("f4")
    return [p, v]


@pytest.mark.parametrize("func", ["max", "min", "mean"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_over_two_pv_array_matches_reference(func, seed):
    rng = np.random.default_rng(seed)
    pv1 = random_pv(rng, int(rng.integers(1, 60)))
    pv2 = random_pv(rng, int(rng.integers(1, 60)))
    p, v = over_two_pv_array(pv1, pv2, func=func)
    ep, ev = ref_over_two(pv1, pv2, func)
    assert p.dtype == np.int32 and v.dtype == np.float32
    np.testing.assert_array_equal(p, ep)
    np.testing.assert_array_equal(v, ev)


@pytest.mark.parametrize("func", ["max", "min", "mean"])
def test_over_two_pv_array_stops_at_shorter_array(func):
    # pv2 extends to 20 but the result covers only the overlap [0, 10)
    pv1 = [i4([5, 10]), f4([1.0, 2.0])]
    pv2 = [i4([8, 20]), f4([4.0, 0.0])]
    # segments [0, 5): (1, 4), [5, 8): (2, 4), [8, 10): (2, 0)
    expect = {"max": [4.0, 4.0, 2.0], "min": [1.0, 2.0, 0.0],
              "mean": [2.5, 3.0, 1.0]}
    p, v = over_two_pv_array(pv1, pv2, func=func)
    np.testing.assert_array_equal(p, [5, 8, 10])
    np.testing.assert_array_equal(v, expect[func])


@pytest.mark.parametrize("func", ["max", "min", "mean"])
def test_over_two_pv_array_identical_inputs(func):
    pv = [i4([3, 9, 15]), f4([0.5, 7.25, 1.0])]
    p, v = over_two_pv_array(pv, [pv[0].copy(), pv[1].copy()], func=func)
    np.testing.assert_array_equal(p, pv[0])
    np.testing.assert_array_equal(v, pv[1])


def test_over_two_pv_array_empty_side():
    p, v = over_two_pv_array([i4([]), f4([])], [i4([5]), f4([1.0])])
    assert p.shape == (0,) and v.shape == (0,)


def test_over_two_pv_array_default_is_max():
    pv1 = [i4([5, 10]), f4([1.0, 2.0])]
    pv2 = [i4([5, 10]), f4([3.0, 0.0])]
    p, v = over_two_pv_array(pv1, pv2)
    np.testing.assert_array_equal(v, [3.0, 2.0])


def test_over_two_pv_array_invalid_function():
    pv = [i4([5]), f4([1.0])]
    with pytest.raises(Exception, match="Invalid function"):
        over_two_pv_array(pv, pv, func="sum")


def test_over_two_pv_array_requires_lists():
    pv = (i4([5]), f4([1.0]))
    with pytest.raises(TypeError):
        over_two_pv_array(pv, [i4([5]), f4([1.0])])


# ------------------------------------
# naive_call_peaks
# ------------------------------------
# Each case is a bedGraph-like [p, v] array (segment i is
# [p[i-1], p[i]) with value v[i]); the expected (summit, height) pairs
# are derived by hand in the comment of each case.

NAIVE_PEAK_CASES = [
    # one region [100, 400) above 2, length 300 >= 200: summit 250
    ("single", [100, 400, 500], [0, 5, 0], {}, [(250, 5.0)]),
    # region [100, 250) is 150 long < 200: dropped
    ("too_short", [100, 250, 500], [0, 5, 0], {}, []),
    # length exactly min_length 200 is kept: summit 200
    ("min_length_boundary", [100, 300, 500], [0, 5, 0], {}, [(200, 5.0)]),
    # gap [200, 230) of 30 <= 50 merges; highest segment [230, 400): 315
    ("merge_gap", [100, 200, 230, 400, 600], [0, 5, 0, 7, 0], {},
     [(315, 7.0)]),
    # gap [400, 500) of 100 > 50 splits into two peaks
    ("split_gap", [100, 400, 500, 800, 900], [0, 5, 0, 3, 0], {},
     [(250, 5.0), (650, 3.0)]),
    # gap of exactly max_gap 50 merges; summit stays in [100, 400)
    ("gap_boundary", [100, 400, 450, 800, 900], [0, 5, 0, 3, 0], {},
     [(250, 5.0)]),
    # height 5 is not < max_v 5: dropped
    ("max_v_drop", [100, 400, 500], [0, 5, 0], {"max_v": 5.0}, []),
    ("max_v_keep", [100, 400, 500], [0, 5, 0], {"max_v": 5.5}, [(250, 5.0)]),
    # two tied tops at midpoints 150, 250: index int(3/2)-1 = 0 -> 150
    ("tie2", [100, 200, 300, 600], [0, 5, 5, 0], {}, [(150, 5.0)]),
    # three tied tops 150, 250, 350: index int(4/2)-1 = 1 -> 250
    ("tie3", [100, 200, 300, 400, 600], [0, 5, 5, 5, 0], {}, [(250, 5.0)]),
    # four tied tops: index int(5/2)-1 = 1 -> 250
    ("tie4", [100, 200, 300, 400, 500, 700], [0, 5, 5, 5, 5, 0], {},
     [(250, 5.0)]),
    # first segment [0, 300) above threshold: summit 150
    ("from_zero", [300, 400], [5, 0], {}, [(150, 5.0)]),
    # value equal to min_v is not above it
    ("equal_min_v", [100, 400], [2, 2], {}, []),
    # region running to the end of the array is closed: summit 250
    ("open_end", [100, 400], [0, 5], {}, [(250, 5.0)]),
    # short first region is discarded, later one kept: summit 650
    ("short_then_long", [100, 150, 500, 800, 900], [0, 5, 0, 4, 0], {},
     [(650, 4.0)]),
    # summit midpoint is floored: int((101 + 400) / 2) = 250
    ("floor_midpoint", [101, 400, 500], [0, 5, 0], {}, [(250, 5.0)]),
    # custom min_length and max_gap
    ("custom_params", [10, 20, 25, 40, 100], [0, 3, 0, 4, 0],
     {"max_gap": 5, "min_length": 30}, [(32, 4.0)]),
]


@pytest.mark.parametrize("ps, vs, kwargs, expected",
                         [c[1:] for c in NAIVE_PEAK_CASES],
                         ids=[c[0] for c in NAIVE_PEAK_CASES])
def test_naive_call_peaks(ps, vs, kwargs, expected):
    peaks = naive_call_peaks([i4(ps), f4(vs)], 2.0, **kwargs)
    assert peaks == expected


def test_naive_call_peaks_empty():
    assert naive_call_peaks([i4([]), f4([])], 2.0) == []


# ------------------------------------
# pileup_and_write_se / pileup_and_write_pe
# ------------------------------------

def se_shifts(d, directional, halfextension):
    """(five_shift, three_shift) as defined in pileup_and_write_se."""
    if directional:
        return (d // -4, d * 3 // 4) if halfextension else (0, d)
    return (d // 4, d // 4) if halfextension else (d // 2, d - d // 2)


def make_fwtrack():
    rng = np.random.default_rng(42)
    fw = FWTrack(fw=50)
    reads = [(b"chr1", int(x), 0) for x in rng.integers(0, 3000, 40)]
    reads += [(b"chr1", int(x), 1) for x in rng.integers(0, 3000, 30)]
    reads += [(b"chr2", int(x), 0) for x in rng.integers(0, 800, 10)]
    reads += [(b"chr2", int(x), 1) for x in rng.integers(0, 800, 15)]
    # reads whose fragments reach past position 0
    reads += [(b"chr1", 0, 0), (b"chr1", 3, 1), (b"chr2", 10, 1)]
    for c, p, s in reads:
        fw.add_loc(c, p, s)
    fw.finalize()
    return fw


def expected_se_bdg(track, d, scale, baseline, directional, halfextension):
    five, three = se_shifts(d, directional, halfextension)
    text = ""
    rlengths = track.get_rlengths()
    for chrom in list(rlengths.keys()):
        plus, minus = track.get_locations_by_chr(chrom)
        starts, ends = se_intervals(plus, minus, five, three,
                                    rlengths[chrom])
        text += bdg_text(chrom.decode(),
                         ref_pileup(starts, ends, scale, baseline))
    return text


@pytest.mark.parametrize("d", [200, 201])
@pytest.mark.parametrize("directional, halfextension",
                         [(True, True), (True, False),
                          (False, True), (False, False)])
def test_pileup_and_write_se_matches_reference(tmp_path, d, directional,
                                               halfextension):
    fw = make_fwtrack()
    fw.set_rlengths({b"chr1": 2500})    # chr2 gets INT32_MAX
    out = tmp_path / "se.bdg"
    pileup_and_write_se(fw, str(out).encode(), d, 1.0,
                        directional=directional, halfextension=halfextension)
    assert out.read_text() == expected_se_bdg(fw, d, 1.0, 0.0, directional,
                                              halfextension)


def test_pileup_and_write_se_scale_and_baseline(tmp_path):
    fw = make_fwtrack()
    out = tmp_path / "se.bdg"
    pileup_and_write_se(fw, str(out).encode(), 150, 0.3, baseline_value=0.5,
                        directional=True, halfextension=False)
    assert out.read_text() == expected_se_bdg(fw, 150, 0.3, 0.5, True, False)


def test_pileup_and_write_se_default_is_directional_half_extension(tmp_path):
    fw = make_fwtrack()
    out = tmp_path / "se.bdg"
    pileup_and_write_se(fw, str(out).encode(), 200, 1.0)
    assert out.read_text() == expected_se_bdg(fw, 200, 1.0, 0.0, True, True)


def test_pileup_and_write_se_overwrites_existing_file(tmp_path):
    fw = make_fwtrack()
    out = tmp_path / "se.bdg"
    out.write_text("stale content\n")
    pileup_and_write_se(fw, str(out).encode(), 100, 1.0,
                        directional=True, halfextension=False)
    assert out.read_text() == expected_se_bdg(fw, 100, 1.0, 0.0, True, False)


def test_pileup_and_write_se_single_read_by_hand(tmp_path):
    # plus 100, d 50, directional, no half extension: [100, 150)
    fw = FWTrack(fw=50)
    fw.add_loc(b"chrA", 100, 0)
    fw.finalize()
    out = tmp_path / "one.bdg"
    pileup_and_write_se(fw, str(out).encode(), 50, 1.0,
                        directional=True, halfextension=False)
    assert out.read_text() == ("chrA\t0\t100\t0.00000\n"
                               "chrA\t100\t150\t1.00000\n")


def test_pileup_and_write_se_filename_must_be_bytes(tmp_path):
    fw = make_fwtrack()
    with pytest.raises(TypeError):
        pileup_and_write_se(fw, str(tmp_path / "x.bdg"), 100, 1.0)


def make_petrack():
    rng = np.random.default_rng(7)
    pe = PETrackI()
    for chrom, n, span in ((b"chr2", 60, 4000), (b"chr1", 25, 1500)):
        starts = rng.integers(0, span, n)
        lens = rng.integers(20, 400, n)
        for s, ln in zip(starts.tolist(), lens.tolist()):
            pe.add_loc(chrom, s, s + ln)
    pe.add_loc(b"chr1", 0, 30)          # fragment starting at 0
    pe.add_loc(b"chr1", 0, 30)          # and its duplicate
    pe.finalize()
    return pe


def expected_pe_bdg(track, scale, baseline):
    text = ""
    for chrom in list(track.get_rlengths().keys()):
        locs = track.get_locations_by_chr(chrom)
        text += bdg_text(chrom.decode(),
                         ref_pileup(locs['l'], locs['r'], scale, baseline))
    return text


@pytest.mark.parametrize("scale, baseline", [(1.0, 0.0), (0.25, 0.0),
                                             (2.0, 1.0), (0.1, 0.05)])
def test_pileup_and_write_pe_matches_reference(tmp_path, scale, baseline):
    pe = make_petrack()
    out = tmp_path / "pe.bdg"
    pileup_and_write_pe(pe, str(out).encode(), scale_factor=scale,
                        baseline_value=baseline)
    assert out.read_text() == expected_pe_bdg(pe, scale, baseline)


def test_pileup_and_write_pe_single_fragment_by_hand(tmp_path):
    pe = PETrackI()
    pe.add_loc(b"chrA", 10, 60)
    pe.finalize()
    out = tmp_path / "one.bdg"
    pileup_and_write_pe(pe, str(out).encode())
    assert out.read_text() == ("chrA\t0\t10\t0.00000\n"
                               "chrA\t10\t60\t1.00000\n")


def test_pileup_and_write_pe_overwrites_existing_file(tmp_path):
    pe = make_petrack()
    out = tmp_path / "pe.bdg"
    out.write_text("stale content\n")
    pileup_and_write_pe(pe, str(out).encode())
    assert out.read_text() == expected_pe_bdg(pe, 1.0, 0.0)


def test_pileup_and_write_pe_filename_must_be_bytes(tmp_path):
    with pytest.raises(TypeError):
        pileup_and_write_pe(make_petrack(), str(tmp_path / "x.bdg"))
