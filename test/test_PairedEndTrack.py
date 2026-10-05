#!/usr/bin/env python
# Time-stamp: <2025-09-29 14:57:48 Tao Liu>

import unittest
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Signal.PileupV2 import (pileup_from_LRC,
                                   pileup_from_LRC_centers_as_list,
                                   over_two_pv_array)
from MACS3.Signal.Region import Regions
import numpy as np
import anndata  # noqa: F401
import pandas  # noqa: F401
import scipy  # noqa: F401

import io
import itertools
import logging
from collections import Counter

import pytest



class Test_PETrackI(unittest.TestCase):
    def setUp(self):
        self.input_regions = [(b"chrY", 0, 100),  # 17 tags in chrY
                              (b"chrY", 70, 270),  # will exclude
                              (b"chrY", 70, 100),  # will exclude
                              (b"chrY", 80, 160),  # will exclude
                              (b"chrY", 80, 160),  # will exclude
                              (b"chrY", 80, 180),  # will exclude
                              (b"chrY", 80, 180),  # will exclude
                              (b"chrY", 85, 185),  # will exclude
                              (b"chrY", 85, 285),  # will exclude
                              (b"chrY", 85, 285),  # will exclude
                              (b"chrY", 85, 285),  # will exclude
                              (b"chrY", 85, 385),  # will exclude
                              (b"chrY", 90, 190),  # will exclude
                              (b"chrY", 90, 190),  # will exclude
                              (b"chrY", 90, 191),  # will exclude
                              (b"chrY", 150, 190),  # will exclude
                              (b"chrY", 150, 250),  # will exclude
                              (b"chrX", 0, 100),  # 9 tags in chrX
                              (b"chrX", 70, 270),  # will exclude
                              (b"chrX", 70, 100),  # will exclude
                              (b"chrX", 80, 160),
                              (b"chrX", 80, 180),
                              (b"chrX", 85, 185),
                              (b"chrX", 85, 285),
                              (b"chrX", 90, 190),
                              (b"chrX", 90, 191)
                              ]
        self.t = sum([x[2]-x[1] for x in self.input_regions])

        self.exclude_regions = [(b"chrY", 85, 200), (b"chrX", 50, 75)]
        self.x_regions = Regions()
        for r in self.exclude_regions:
            self.x_regions.add_loc(r[0], r[1], r[2])

        self.pe = PETrackI()
        for (c, l, r) in self.input_regions:
            self.pe.add_loc(c, l, r)
        self.pe.finalize()

    def test_add_loc(self):
        self.assertEqual(self.pe.total, 26)
        self.assertEqual(self.pe.length, self.t)

    def test_filter_dup(self):
        self.pe.filter_dup(3)
        self.assertEqual(self.pe.total, 26)

        self.pe.filter_dup(2)
        self.assertEqual(self.pe.total, 25)

        self.pe.filter_dup(1)
        self.assertEqual(self.pe.total, 21)

    def test_sample_num(self):
        self.pe.sample_num(10)  # percentage is 38.5%
        self.assertEqual(self.pe.total, 9)

    def test_sample_percent(self):
        self.pe.sample_percent(0.5)  # 12=int(0.5*17)+int(0.5*9)
        self.assertEqual(self.pe.total, 12)

    def test_pileupbdg(self):
        self.pe.pileup_bdg()

    def test_exclude(self):
        self.pe.exclude(self.x_regions)
        self.assertEqual(self.pe.total, 7)
        self.assertEqual(self.pe.length, 781)
        self.assertAlmostEqual(self.pe.average_template_length, 112, 0)


class Test_PETrackII(unittest.TestCase):
    def setUp(self):
        self.input_regions = [(b"chrY", 0, 100, b"0w#AAACGAAAGACTCGGA", 2),  # will be excluded
                              (b"chrY", 70, 170, b"0w#AAACGAAAGACTCGGA", 1),  # will be excluded
                              (b"chrY", 80, 190, b"0w#AAACGAAAGACTCGGA", 1),
                              (b"chrY", 85, 180, b"0w#AAACGAAAGACTCGGA", 3),
                              (b"chrY", 100, 190, b"0w#AAACGAAAGACTCGGA", 1),
                              (b"chrY", 0, 100, b"0w#AAACGAACAAGTAACA", 1),  # will be excluded
                              (b"chrY", 70, 170, b"0w#AAACGAACAAGTAACA", 2),  # will be excluded
                              (b"chrY", 80, 190, b"0w#AAACGAACAAGTAACA", 1),
                              (b"chrY", 85, 180, b"0w#AAACGAACAAGTAACA", 1),
                              (b"chrY", 100, 190, b"0w#AAACGAACAAGTAACA", 3),
                              (b"chrY", 10, 110, b"0w#AAACGAACAAGTAAGA", 1),  # will be excluded
                              (b"chrY", 50, 160, b"0w#AAACGAACAAGTAAGA", 2),  # will be excluded
                              (b"chrY", 100, 170, b"0w#AAACGAACAAGTAAGA", 3)
                              ]
        self.exclude_regions = [(b"chrY", 10, 75)]
        self.x_regions = Regions()
        for r in self.exclude_regions:
            self.x_regions.add_loc(r[0], r[1], r[2])

        self.pileup_p = np.array([10, 50, 70, 80, 85,
                                  100, 110, 160, 170, 180,
                                  190], dtype="i4")
        self.pileup_v = np.array([3.0, 4.0, 6.0, 9.0, 11.0,
                                  15.0, 19.0, 18.0, 16.0, 10.0,
                                  6.0], dtype="f4")
        self.peak_str = "chrom:chrY	start:80	end:180	name:peak_1	score:19	summit:105\n"
        self.subset_barcodes = {b'0w#AAACGAACAAGTAACA', b"0w#AAACGAACAAGTAAGA"}
        self.subset_pileup_p = np.array([10, 50, 70, 80, 85,
                                         100, 110, 160, 170, 180,
                                         190], dtype="i4")
        self.subset_pileup_v = np.array([1.0, 2.0, 4.0, 6.0, 7.0,
                                         8.0, 13.0, 12.0, 10.0, 5.0,
                                         4.0], dtype="f4")
        self.subset_peak_str = "chrom:chrY	start:100	end:170	name:peak_1	score:13	summit:105\n"
        self.t = sum([(x[2]-x[1]) * x[4] for x in self.input_regions])

        self.pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            self.pe.add_loc(c, l, r, b, C)
        self.pe.finalize()

    def test_add_frag(self):
        self.assertEqual(self.pe.total, 22)
        self.assertEqual(self.pe.length, self.t)

        pe_subset = self.pe.subset(self.subset_barcodes)
        self.assertEqual(pe_subset.total, 14)
        self.assertEqual(pe_subset.length, 1305)

    def test_pileup(self):
        bdg = self.pe.pileup_bdg()
        d = bdg.get_data_by_chr(b'chrY')
        np.testing.assert_array_equal(d[0], self.pileup_p)
        np.testing.assert_array_equal(d[1], self.pileup_v)

        pe_subset = self.pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg()
        d = bdg.get_data_by_chr(b'chrY')
        np.testing.assert_array_equal(d[0], self.subset_pileup_p)
        np.testing.assert_array_equal(d[1], self.subset_pileup_v)

    def test_pileup2(self):
        bdg = self.pe.pileup_bdg2()
        d = bdg.get_data_by_chr(b'chrY')
        np.testing.assert_array_equal(d['p'], self.pileup_p)
        np.testing.assert_array_equal(d['v'], self.pileup_v)

        pe_subset = self.pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg2()
        d = bdg.get_data_by_chr(b'chrY')
        np.testing.assert_array_equal(d['p'], self.subset_pileup_p)
        np.testing.assert_array_equal(d['v'], self.subset_pileup_v)

    def test_callpeak(self):
        bdg = self.pe.pileup_bdg()
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.peak_str)

        pe_subset = self.pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg()
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.subset_peak_str)

    def test_callpeak2(self):
        bdg = self.pe.pileup_bdg2()
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.peak_str)

        pe_subset = self.pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg2()
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.subset_peak_str)

    def test_control_pileup_matches_legacy_temp_arrays(self):
        locs = np.array([(10, 90, 2),
                         (10, 90, 1),
                         (25, 70, 3),
                         (40, 100, 2),
                         (60, 80, 4)],
                        dtype=[("l", "i4"), ("r", "i4"), ("c", "u2")])
        ds = [20, 45]
        scales = [1.5, 0.25]
        baseline = 0.5
        expected = None

        for d, scale in zip(ds, scales):
            half_d = d // 2
            tmp_arr_l = locs.copy()
            tmp_arr_l["l"] = tmp_arr_l["l"] - half_d
            tmp_arr_l["r"] = tmp_arr_l["l"] + d

            tmp_arr_r = locs.copy()
            tmp_arr_r["l"] = tmp_arr_r["r"] - half_d
            tmp_arr_r["r"] = tmp_arr_r["l"] + d

            tmp_arr = np.concatenate([tmp_arr_l, tmp_arr_r])
            pv = pileup_from_LRC(tmp_arr)
            v = pv["v"] * scale
            v[v < baseline] = baseline
            current = [pv["p"], v]

            if expected:
                expected = over_two_pv_array(expected, current, func="max")
            else:
                expected = current

        pe = PETrackII()
        for i, (left, right, count) in enumerate(locs):
            pe.add_loc(b"chr1", int(left), int(right), b"BC%d" % i, int(count))
        pe.finalize()
        observed = pe.pileup_a_chromosome_c(b"chr1", ds, scales, baseline)

        np.testing.assert_array_equal(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])

    def test_lrc_center_helper_handles_unsorted_right_endpoints(self):
        locs = np.array([(1, 100, 2),
                         (2, 20, 3),
                         (2, 10, 1),
                         (5, 30, 4)],
                        dtype=[("l", "i4"), ("r", "i4"), ("c", "u2")])
        d = 10
        scale = 0.75
        baseline = 0.25
        half_d = d // 2
        tmp_arr_l = locs.copy()
        tmp_arr_l["l"] = tmp_arr_l["l"] - half_d
        tmp_arr_l["r"] = tmp_arr_l["l"] + d
        tmp_arr_r = locs.copy()
        tmp_arr_r["l"] = tmp_arr_r["r"] - half_d
        tmp_arr_r["r"] = tmp_arr_r["l"] + d
        pv = pileup_from_LRC(np.concatenate([tmp_arr_l, tmp_arr_r]))
        expected_v = pv["v"] * scale
        expected_v[expected_v < baseline] = baseline

        observed = pileup_from_LRC_centers_as_list(locs, d, scale, baseline)
        np.testing.assert_array_equal(observed[0], pv["p"])
        np.testing.assert_allclose(observed[1], expected_v)

    def test_exclude(self):
        self.pe.exclude(self.x_regions)
        self.assertEqual(self.pe.total, 13)
        self.assertEqual(self.pe.length, 1170)
        self.assertAlmostEqual(self.pe.average_template_length, 90, 0)

    def test_exclude_merges_overlapping_regions(self):
        pe = PETrackII()
        pe.add_loc(b"chr1", 0, 10, b"A", 1)
        pe.add_loc(b"chr1", 12, 18, b"B", 2)
        pe.add_loc(b"chr1", 26, 30, b"C", 3)
        pe.finalize()

        regions = Regions()
        regions.add_loc(b"chr1", 5, 15)
        regions.add_loc(b"chr1", 14, 25)

        pe.exclude(regions)

        self.assertEqual(regions.regions[b"chr1"], [(5, 15), (14, 25)])
        self.assertEqual(pe.total, 3)
        self.assertEqual(pe.length, 12)
        remaining = pe.get_locations_by_chr(b"chr1")
        self.assertEqual(remaining.shape[0], 1)
        self.assertEqual(int(remaining[0]['l']), 26)
        self.assertEqual(int(remaining[0]['r']), 30)
        self.assertEqual(int(remaining[0]['c']), 3)

    def test_return_anndata(self):
        petrack = PETrackII()
        petrack.add_loc(b"chr1", 0, 100, barcode=b"A", count=2)   # peak_1
        petrack.add_loc(b"chr1", 70, 270, barcode=b"A", count=1)   # peak_2
        petrack.add_loc(b"chr1", 0, 100, barcode=b"B", count=3)   # peak_2
        petrack.add_loc(b"chr1", 175, 325, barcode=b"C", count=4)   # peak_2
        petrack.finalize()

        regions = Regions()
        regions.add_loc(b"chr1", 0, 100)    # peak_1
        regions.add_loc(b"chr1", 200, 300)  # peak_2
        regions.add_loc(b"chr1", 500, 600)  # peak_3


        adata = petrack.return_anndata(regions)
        self.assertEqual(adata.shape, (3, 3))
        self.assertEqual(list(adata.obs.index), ["A", "B", "C"])
        self.assertEqual(list(adata.var.index), ["peak_1", "peak_2", "peak_3"])

        X = adata.X.toarray()
        expected = np.array([[3, 1, 0],
                              [3, 0, 0],
                              [0, 4, 0]], dtype=np.int32)
        np.testing.assert_array_equal(X, expected)

        # Explicit peak-inside-fragment check with simpler labels
        petrack = PETrackII()
        petrack.add_loc(b"chr1", 0, 100, barcode=b"A", count=1)   # contains peak_1
        petrack.add_loc(b"chr1", 175, 325, barcode=b"B", count=4)  # contains peak_2
        petrack.finalize()

        regions = Regions()
        regions.add_loc(b"chr1", 10, 90)    # peak_1
        regions.add_loc(b"chr1", 200, 300)  # peak_2

        adata = petrack.return_anndata(regions)
        self.assertEqual(adata.shape, (2, 2))
        self.assertEqual(list(adata.obs.index), ["A", "B"])
        self.assertEqual(list(adata.var.index), ["peak_1", "peak_2"])

        X = adata.X.toarray()
        expected = np.array([[1, 0],
                              [0, 4]], dtype=np.int32)
        np.testing.assert_array_equal(X, expected)

    def test_return_anndata_merges_overlapping_regions(self):
        petrack = PETrackII()
        petrack.add_loc(b"chr1", 0, 150, barcode=b"A", count=2)
        petrack.add_loc(b"chr1", 190, 260, barcode=b"B", count=3)
        petrack.finalize()

        regions = Regions()
        regions.add_loc(b"chr1", 10, 90)
        regions.add_loc(b"chr1", 80, 120)
        regions.add_loc(b"chr1", 200, 240)

        adata = petrack.return_anndata(regions)

        self.assertEqual(adata.shape, (2, 2))
        self.assertEqual(list(adata.obs.index), ["A", "B"])
        self.assertEqual(list(adata.var.index), ["peak_1", "peak_2"])
        self.assertEqual(regions.regions[b"chr1"], [(10, 90), (80, 120), (200, 240)])

        X = adata.X.toarray()
        expected = np.array([[2, 0],
                              [0, 3]], dtype=np.int32)
        np.testing.assert_array_equal(X, expected)

    def test_return_anndata_merges_adjacent_regions(self):
        petrack = PETrackII()
        petrack.add_loc(b"chr1", 20, 90, barcode=b"A", count=2)
        petrack.add_loc(b"chr1", 100, 170, barcode=b"B", count=5)
        petrack.finalize()

        regions = Regions()
        regions.add_loc(b"chr1", 10, 50)
        regions.add_loc(b"chr1", 50, 80)
        regions.add_loc(b"chr1", 120, 160)

        adata = petrack.return_anndata(regions)
        self.assertEqual(adata.shape, (2, 2))
        self.assertEqual(list(adata.obs.index), ["A", "B"])
        self.assertEqual(list(adata.var.index), ["peak_1", "peak_2"])
        self.assertEqual(regions.regions[b"chr1"], [(10, 50), (50, 80), (120, 160)])

        X = adata.X.toarray()
        expected = np.array([[2, 0],
                              [0, 5]], dtype=np.int32)
        np.testing.assert_array_equal(X, expected)

class TestPETrackIISampling(unittest.TestCase):

    def setUp(self):
        # Create a PETrackII with two chromosomes, three fragments each, and known counts
        # Suppose: locs dtype = [('l', 'i4'), ('r', 'i4'), ('c', 'u2')]
        self.petrack = PETrackII()
        self.petrack.locations = {}
        self.petrack.barcodes = {}
        self.petrack.size = {}
        self.petrack.buf_size = {}
        chroms = [b'chr1', b'chr2']
        for chrom in chroms:
            locs = np.array([(0, 10, 5), (10, 20, 3), (20, 30, 2)],
                            dtype=[('l', 'i4'), ('r', 'i4'), ('c', 'u2')])
            bars = np.array([1, 2, 3], dtype='i4')
            self.petrack.locations[chrom] = locs.copy()
            self.petrack.barcodes[chrom] = bars.copy()
            self.petrack.size[chrom] = 3
            self.petrack.buf_size[chrom] = 3

    def test_sample_percent(self):
        # In-place, 50% downsampling
        total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        self.petrack.sample_percent(0.5, seed=42)
        new_total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        self.assertAlmostEqual(new_total, round(total * 0.5), delta=2)  # allow rounding error
        # Should not have counts greater than the originals for any fragment
        for k in self.petrack.get_chr_names():
            orig = np.array([5, 3, 2])
            new = np.zeros_like(orig)
            for i, loc in enumerate(self.petrack.locations[k]):
                new[i] = loc['c']
            self.assertTrue(np.all(new <= orig))

    def test_sample_percent_copy(self):
        # Copy, 30% downsampling
        total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        petrack2 = self.petrack.sample_percent_copy(0.3, seed=123)
        new_total = sum(petrack2.locations[k]['c'].sum() for k in petrack2.get_chr_names())
        self.assertAlmostEqual(new_total, round(total * 0.3), delta=2)
        # Originals should remain unchanged
        orig_total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        self.assertEqual(orig_total, total)

    def test_sample_num(self):
        # In-place, absolute downsampling
        total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        target = 4
        self.petrack.sample_num(target, seed=1)
        new_total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        self.assertAlmostEqual(new_total, target, delta=1)
        # Should not have more than original in any fragment
        for k in self.petrack.get_chr_names():
            orig = np.array([5, 3, 2])
            new = np.zeros_like(orig)
            for i, loc in enumerate(self.petrack.locations[k]):
                new[i] = loc['c']
            self.assertTrue(np.all(new <= orig))

    def test_sample_num_copy(self):
        # Copy, absolute downsampling
        total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        target = 7
        petrack2 = self.petrack.sample_num_copy(target, seed=99)
        new_total = sum(petrack2.locations[k]['c'].sum() for k in petrack2.get_chr_names())
        self.assertAlmostEqual(new_total, target, delta=1)
        # Originals should remain unchanged
        orig_total = sum(self.petrack.locations[k]['c'].sum() for k in self.petrack.get_chr_names())
        self.assertEqual(orig_total, total)


# ------------------------------------
# Reference implementation and helpers for the tests below
# ------------------------------------
#
# The reference pileup is the (count-weighted) coverage of half-open
# fragments [start, end) accumulated on the elementary segments between
# the sorted distinct coordinates, with values
# max(coverage * scale, baseline) in float32. Results are compared as
# change points (pos, value): segment i covers [pos[i-1], pos[i]) with
# pos[-1] read as 0, equal neighbours merged.

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


def ref_pileup(starts, ends, weights=None, scale=1.0, baseline=0.0):
    """Reference weighted coverage of [start, end) intervals."""
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
    cov = np.cumsum(delta)[:-1].astype(np.float32)
    vals = np.maximum(cov * np.float32(scale), np.float32(baseline))
    return merge_runs(coords[1:], vals)


def ref_max_truncated(cps):
    """Max over several pileups on the union of their breakpoints, up to
    the shortest pileup's end (over_two_pv_array merges the overlap)."""
    if any(cp[0].size == 0 for cp in cps):
        return np.zeros(0, np.int64), np.zeros(0, np.float32)
    end = min(cp[0][-1] for cp in cps)
    bps = np.unique(np.concatenate([cp[0] for cp in cps]))
    bps = bps[bps <= end]
    vals = [cp[1][np.searchsorted(cp[0], bps)] for cp in cps]
    return merge_runs(bps, np.maximum.reduce(vals))


def end_windows(locs, lo, width, rlength=None):
    """Windows [x - lo, x - lo + width) around both ends x of each
    fragment, clipped to [0, rlength] when rlength is given."""
    x = np.concatenate((locs['l'], locs['r'])).astype(np.int64)
    starts, ends = x - lo, x - lo + width
    if rlength is not None:
        starts, ends = np.clip(starts, 0, rlength), np.clip(ends, 0, rlength)
    return starts, ends


def assert_cp_equal(p, v, expected):
    """Check pileup arrays: int32/float32, strictly increasing positions,
    and the same change points as the reference."""
    p = np.asarray(p)
    v = np.asarray(v)
    assert p.dtype == np.int32 and v.dtype == np.float32
    assert len(p) == len(v)
    assert np.all(np.diff(p.astype(np.int64)) > 0)
    gp, gv = merge_runs(p, v)
    np.testing.assert_array_equal(gp, expected[0])
    np.testing.assert_array_equal(gv, expected[1])


def assert_pv_equal(pv, expected):
    """PileupV2 structured output is already merged: compare it raw."""
    assert pv.dtype == np.dtype([('p', 'i4'), ('v', 'f4')])
    np.testing.assert_array_equal(pv['p'], expected[0])
    np.testing.assert_array_equal(pv['v'], expected[1])


def f32_mean(length, total):
    # float(length) / total evaluated in float32, as in finalize()
    return float(np.float32(np.float32(length) / np.float32(total)))


def keep_at_most(values, maxnum):
    """Sorted values with every distinct value kept at most maxnum times."""
    return [v for v, n in sorted(Counter(values).items())
            for _ in range(min(n, maxnum))]


def overlaps(l, r, regions):
    return any(s < r and l < e for s, e in regions)


def make_regions(regions):
    rg = Regions()
    for chrom, items in regions.items():
        for s, e in items:
            rg.add_loc(chrom, s, e)
    return rg


# ------------------------------------
# PETrackI
# ------------------------------------

PE_FRAGS = [(b"chrA", 100, 250), (b"chrA", 30, 80), (b"chrA", 30, 80),
            (b"chrA", 30, 90), (b"chrA", 0, 40), (b"chrA", 500, 650),
            (b"chrA", 30, 80), (b"chrA", 120, 121),
            (b"chrB", 1000, 1300), (b"chrB", 1000, 1300), (b"chrB", 50, 400),
            (b"chrB", 2000, 2010),
            (b"chrC", 7, 77)]


def build_pe(frags, buffer_size=100000):
    pe = PETrackI(buffer_size=buffer_size)
    for c, l, r in frags:
        pe.add_loc(c, l, r)
    pe.finalize()
    return pe


def expected_pe(frags):
    """{chrom: sorted [(l, r), ...]}."""
    out = {}
    for c, l, r in frags:
        out.setdefault(c, []).append((l, r))
    return {c: sorted(v) for c, v in out.items()}


def pe_records(pe):
    return {c: [tuple(x) for x in pe.get_locations_by_chr(c).tolist()]
            for c in pe.get_chr_names()}


def random_pe(seed, n=80, span=4000, chroms=(b"chr1",)):
    rng = np.random.default_rng(seed)
    pe = PETrackI()
    for chrom in chroms:
        starts = rng.integers(0, span, n)
        lens = rng.integers(20, 400, n)
        for s, ln in zip(starts.tolist(), lens.tolist()):
            pe.add_loc(chrom, s, s + ln)
    pe.add_loc(chroms[0], 0, 25)        # a fragment starting at 0
    pe.finalize()
    return pe


def test_petrackI_init_defaults():
    pe = PETrackI()
    assert (pe.total, pe.length, pe.annotation, pe.buffer_size) == \
        (0, 0, "", 100000)
    assert pe.average_template_length == 0.0
    assert pe.is_sorted is False
    assert (pe.locations, pe.size, pe.rlengths) == ({}, {}, {})


def test_petrackI_init_arguments():
    pe = PETrackI(anno="ctrl", buffer_size=7)
    assert (pe.annotation, pe.buffer_size) == ("ctrl", 7)


@pytest.mark.parametrize("buffer_size", [1, 2, 3, 100000])
def test_petrackI_add_loc_finalize(buffer_size):
    pe = build_pe(PE_FRAGS, buffer_size)
    exp = expected_pe(PE_FRAGS)
    assert pe_records(pe) == exp
    assert pe.size == {c: len(v) for c, v in exp.items()}
    length = sum(r - l for _, l, r in PE_FRAGS)
    assert pe.total == len(PE_FRAGS)
    assert pe.length == length
    assert pe.average_template_length == f32_mean(length, len(PE_FRAGS))
    assert pe.is_sorted is True


def test_petrackI_add_loc_updates_length_before_finalize():
    pe = PETrackI()
    pe.add_loc(b"chrA", 10, 40)
    pe.add_loc(b"chrA", 50, 90)
    assert pe.length == 70
    assert pe.size == {b"chrA": 2}
    assert pe.total == 0                # counted by finalize()


def test_petrackI_add_loc_int32_extreme():
    pe = build_pe([(b"chrA", INT32_MAX - 10, INT32_MAX)])
    assert pe_records(pe) == {b"chrA": [(INT32_MAX - 10, INT32_MAX)]}
    assert pe.length == 10


@pytest.mark.parametrize("start, end", [(2**31, 2**31 + 5), (0, 2**31)])
def test_petrackI_add_loc_outside_int32_raises(start, end):
    with pytest.raises(OverflowError):
        PETrackI().add_loc(b"chrA", start, end)


def test_petrackI_add_loc_chromosome_must_be_bytes():
    with pytest.raises(TypeError, match="expected bytes"):
        PETrackI().add_loc("chrA", 10, 20)


def test_petrackI_destroy():
    pe = build_pe(PE_FRAGS)
    pe.destroy()
    assert pe.get_chr_names() == set()
    with pytest.raises(Exception, match="No such chromosome name"):
        pe.get_locations_by_chr(b"chrA")


def test_petrackI_rlengths():
    pe = build_pe(PE_FRAGS)
    assert pe.get_rlengths() == {b"chrA": INT32_MAX, b"chrB": INT32_MAX,
                                 b"chrC": INT32_MAX}
    pe2 = build_pe(PE_FRAGS)
    assert pe2.set_rlengths({b"chrA": 700, b"chrZ": 9}) is True
    assert pe2.get_rlengths() == {b"chrA": 700, b"chrB": INT32_MAX,
                                  b"chrC": INT32_MAX}


def test_petrackI_get_locations_by_chr_missing_raises():
    with pytest.raises(Exception,
                       match=r"No such chromosome name \(b'chrZ'\)"):
        build_pe(PE_FRAGS).get_locations_by_chr(b"chrZ")


def test_petrackI_get_chr_names():
    assert build_pe(PE_FRAGS).get_chr_names() == {b"chrA", b"chrB", b"chrC"}
    assert PETrackI().get_chr_names() == set()


def test_petrackI_sort():
    pe = build_pe(PE_FRAGS)
    for c in pe.get_chr_names():
        np.random.default_rng(1).shuffle(pe.get_locations_by_chr(c))
    pe.is_sorted = False
    pe.sort()
    assert pe.is_sorted is True
    assert pe_records(pe) == expected_pe(PE_FRAGS)


def test_petrackI_count_fraglengths():
    pe = build_pe(PE_FRAGS)
    assert pe.count_fraglengths() == dict(Counter(r - l for _, l, r
                                                  in PE_FRAGS))


def test_petrackI_fraglengths():
    sizes = build_pe(PE_FRAGS).fraglengths()
    assert sizes.dtype == np.int32
    assert sorted(sizes.tolist()) == sorted(r - l for _, l, r in PE_FRAGS)


def test_petrackI_fraglengths_single_chromosome():
    sizes = build_pe([(b"chrA", 5, 10), (b"chrA", 0, 100)]).fraglengths()
    assert sizes.tolist() == [100, 5]   # in (l, r) order


def check_excluded(pe, frags, regions):
    kept = [(c, l, r) for c, l, r in frags
            if not overlaps(l, r, regions.get(c, []))]
    assert pe_records(pe) == expected_pe(kept)
    length = sum(r - l for _, l, r in kept)
    assert pe.total == len(kept)
    assert pe.length == length
    assert pe.average_template_length == f32_mean(length, len(kept))


def test_petrackI_exclude():
    frags = [(b"chrA", l, r) for l, r in ((0, 50), (10, 60), (40, 90),
                                          (100, 150), (120, 170),
                                          (300, 400), (500, 600))]
    frags += [(b"chrB", 5, 50), (b"chrB", 60, 90)]
    regions = {b"chrA": [(55, 110), (130, 140)]}
    pe = build_pe(frags)
    pe.exclude(make_regions(regions))
    check_excluded(pe, frags, regions)


def test_petrackI_exclude_removes_emptied_chromosome():
    frags = [(b"chrA", 0, 50), (b"chrA", 100, 300), (b"chrC", 10, 20)]
    regions = {b"chrC": [(15, 16)]}
    pe = build_pe(frags)
    pe.exclude(make_regions(regions))
    assert pe.get_chr_names() == {b"chrA"}
    check_excluded(pe, frags, regions)


def test_petrackI_exclude_requires_regions():
    with pytest.raises(AssertionError):
        build_pe(PE_FRAGS).exclude([(b"chrA", 0, 10)])


@pytest.mark.parametrize("maxnum", [1, 2, 3])
def test_petrackI_filter_dup(maxnum):
    pe = build_pe(PE_FRAGS)
    assert pe.filter_dup(maxnum) is None
    exp = {c: keep_at_most(v, maxnum) for c, v in
           expected_pe(PE_FRAGS).items()}
    assert pe_records(pe) == exp
    total = sum(len(v) for v in exp.values())
    length = sum(r - l for v in exp.values() for l, r in v)
    assert pe.total == total
    assert pe.length == length
    # self.length / self.total in double, stored as float32
    assert pe.average_template_length == float(np.float32(length / total))


@pytest.mark.parametrize("maxnum", [-1, -5])
def test_petrackI_filter_dup_negative_keeps_all(maxnum):
    pe = build_pe(PE_FRAGS)
    pe.filter_dup(maxnum)
    assert pe_records(pe) == expected_pe(PE_FRAGS)
    assert pe.total == len(PE_FRAGS)


def expected_pe_sample_count(n, percent):
    return int(round(n * float(np.float32(percent)), 5))


def check_pe_sampled(pe, frags, percent):
    exp = expected_pe(frags)
    got = pe_records(pe)
    total = length = 0
    for c, orig in exp.items():
        assert len(got[c]) == expected_pe_sample_count(len(orig), percent)
        assert got[c] == sorted(got[c])
        assert not Counter(got[c]) - Counter(orig)
        total += len(got[c])
        length += sum(r - l for l, r in got[c])
    assert pe.total == total
    assert pe.length == length
    if total:
        assert pe.average_template_length == f32_mean(length, total)


def sample_pe_frags():
    rng = np.random.default_rng(17)
    frags = []
    for chrom, n in ((b"chrA", 40), (b"chrB", 23), (b"chrC", 6)):
        for s in rng.integers(0, 2000, n).tolist():
            frags.append((chrom, s, s + int(rng.integers(30, 300))))
    return frags


@pytest.mark.parametrize("percent", [0.3, 0.5, 0.77, 1.0])
def test_petrackI_sample_percent(percent):
    frags = sample_pe_frags()
    pe = build_pe(frags)
    pe.sample_percent(percent, seed=3)
    check_pe_sampled(pe, frags, percent)


def test_petrackI_sample_percent_without_seed():
    frags = sample_pe_frags()
    pe = build_pe(frags)
    pe.sample_percent(0.5)
    check_pe_sampled(pe, frags, 0.5)


@pytest.mark.parametrize("seed", [0, 42])
def test_petrackI_sample_percent_same_seed_same_sample(seed):
    frags = sample_pe_frags()
    a, b = build_pe(frags), build_pe(frags)
    a.sample_percent(0.4, seed=seed)
    b.sample_percent(0.4, seed=seed)
    assert pe_records(a) == pe_records(b)


def test_petrackI_sample_percent_logs_seed(caplog):
    caplog.set_level(logging.INFO)
    build_pe(sample_pe_frags()).sample_percent(0.5, seed=7)
    assert any(r.getMessage().endswith("#   A random seed 7 has been used")
               for r in caplog.records)


@pytest.mark.parametrize("percent", [0.3, 0.5, 1.0])
def test_petrackI_sample_percent_copy(percent):
    frags = sample_pe_frags()
    pe = build_pe(frags)
    pe.set_rlengths({b"chrA": 5000})
    before = pe_records(pe)
    sub = pe.sample_percent_copy(percent, seed=5)
    assert isinstance(sub, PETrackI)
    assert pe_records(pe) == before     # source untouched
    assert pe.total == len(frags)
    check_pe_sampled(sub, frags, percent)
    assert sub.get_rlengths() == pe.get_rlengths()


def test_petrackI_sample_percent_copy_same_seed_same_sample():
    pe = build_pe(sample_pe_frags())
    a = pe.sample_percent_copy(0.4, seed=9)
    b = pe.sample_percent_copy(0.4, seed=9)
    assert pe_records(a) == pe_records(b)


def test_petrackI_sample_percent_copy_logs_seed(caplog):
    caplog.set_level(logging.INFO)
    build_pe(sample_pe_frags()).sample_percent_copy(0.5, seed=8)
    assert any(r.getMessage().endswith(
        "# A random seed 8 has been used in the sampling function")
        for r in caplog.records)


@pytest.mark.parametrize("samplesize", [10, 23, 69])
def test_petrackI_sample_num(samplesize):
    frags = sample_pe_frags()
    pe = build_pe(frags)
    # percent = float32(samplesize) / total, computed in float32
    percent = np.float32(np.float32(samplesize) / np.float32(len(frags)))
    pe.sample_num(samplesize, seed=2)
    check_pe_sampled(pe, frags, float(percent))


@pytest.mark.parametrize("samplesize", [10, 23, 69])
def test_petrackI_sample_num_copy(samplesize):
    frags = sample_pe_frags()
    pe = build_pe(frags)
    percent = np.float32(np.float32(samplesize) / np.float32(len(frags)))
    sub = pe.sample_num_copy(samplesize, seed=2)
    assert pe.total == len(frags)
    check_pe_sampled(sub, frags, float(percent))


def test_petrackI_print_to_bed_single_chromosome():
    pe = build_pe([(b"chrA", 50, 90), (b"chrA", 10, 40), (b"chrA", 10, 30)])
    buf = io.StringIO()
    pe.print_to_bed(buf)
    assert buf.getvalue() == "chrA\t10\t30\nchrA\t10\t40\nchrA\t50\t90\n"


def test_petrackI_print_to_bed_many_chromosomes():
    pe = build_pe(PE_FRAGS)
    buf = io.StringIO()
    pe.print_to_bed(buf)
    lines = buf.getvalue().splitlines()
    names = [x.split("\t")[0] for x in lines]
    blocks = [k for k, _ in itertools.groupby(names)]
    assert len(blocks) == len(set(blocks)) == 3
    for c, v in expected_pe(PE_FRAGS).items():
        name = c.decode()
        assert [x for x in lines if x.split("\t")[0] == name] == \
            ["%s\t%d\t%d" % (name, l, r) for l, r in v]


def test_petrackI_print_to_bed_defaults_to_stdout(capsys):
    build_pe([(b"chrA", 1, 9)]).print_to_bed()
    assert capsys.readouterr().out == "chrA\t1\t9\n"


def test_petrackI_print_to_bed_requires_file_object():
    with pytest.raises(AssertionError):
        build_pe([(b"chrA", 1, 9)]).print_to_bed(object())


def test_petrackI_pileup_a_chromosome_hand_example():
    # [0, 40) and [30, 80): [0, 30) 1, [30, 40) 2, [40, 80) 1
    pe = build_pe([(b"chrA", 30, 80), (b"chrA", 0, 40)])
    p, v = pe.pileup_a_chromosome(b"chrA")
    np.testing.assert_array_equal(p, [30, 40, 80])
    np.testing.assert_array_equal(v, [1.0, 2.0, 1.0])


@pytest.mark.parametrize("scale, baseline", [(1.0, 0.0), (0.5, 0.0),
                                             (2.0, 3.0), (0.1, 0.05)])
def test_petrackI_pileup_a_chromosome_matches_reference(scale, baseline):
    pe = random_pe(4)
    p, v = pe.pileup_a_chromosome(b"chr1", scale_factor=scale,
                                  baseline_value=baseline)
    locs = pe.get_locations_by_chr(b"chr1")
    assert_cp_equal(p, v, ref_pileup(locs['l'], locs['r'], scale=scale,
                                     baseline=baseline))


def test_petrackI_pileup_a_chromosome_int32_extreme():
    pe = build_pe([(b"chrA", INT32_MAX - 10, INT32_MAX)])
    p, v = pe.pileup_a_chromosome(b"chrA")
    assert_cp_equal(p, v, ([INT32_MAX - 10, INT32_MAX], [0.0, 1.0]))


def test_petrackI_pileup_a_chromosome_missing_chromosome():
    with pytest.raises(KeyError):
        build_pe(PE_FRAGS).pileup_a_chromosome(b"chrZ")


PE_C_CASES = [
    # ds, scale factors, baseline, rlength
    ([200], [1.0], 0.0, None),
    ([201, 1000, 10000], [1.0, 0.2, 0.02], 0.5, None),
    ([300, 2000], [0.5, 0.25], 0.0, 3500),
]


@pytest.mark.parametrize("ds, sfs, baseline, rlength", PE_C_CASES)
def test_petrackI_pileup_a_chromosome_c_matches_reference(ds, sfs, baseline,
                                                          rlength):
    # each fragment end x contributes [x - d//2, x + d//2), as defined by
    # five_shift = three_shift = d//2, clipped to [0, rlength]
    pe = random_pe(6)
    if rlength is not None:
        pe.set_rlengths({b"chr1": rlength})
    else:
        rlength = INT32_MAX
    p, v = pe.pileup_a_chromosome_c(b"chr1", ds, sfs, baseline_value=baseline)
    locs = pe.get_locations_by_chr(b"chr1")
    cps = []
    for d, sf in zip(ds, sfs):
        starts, ends = end_windows(locs, d // 2, 2 * (d // 2), rlength)
        cps.append(ref_pileup(starts, ends, scale=sf, baseline=baseline))
    assert_cp_equal(p, v, ref_max_truncated(cps))


def test_petrackI_pileup_a_chromosome_c_length_mismatch():
    with pytest.raises(AssertionError, match="same length"):
        random_pe(6).pileup_a_chromosome_c(b"chr1", [100, 200], [1.0])


@pytest.mark.parametrize("scale, baseline", [(1.0, 0.0), (0.5, 0.25)])
def test_petrackI_pileup_bdg(scale, baseline):
    pe = random_pe(9, chroms=(b"chr1", b"chr2", b"chr3"))
    bdg = pe.pileup_bdg(scale_factor=scale, baseline_value=baseline)
    assert bdg.get_chr_names() == {b"chr1", b"chr2", b"chr3"}
    assert bdg.baseline_value == baseline
    for c in pe.get_chr_names():
        p, v = bdg.get_data_by_chr(c)
        locs = pe.get_locations_by_chr(c)
        assert_cp_equal(p, v, ref_pileup(locs['l'], locs['r'], scale=scale,
                                         baseline=baseline))


def hmmr_mapping(pe):
    lengths = sorted(pe.count_fraglengths())
    return [{L: 1.0 for L in lengths},
            {L: 0.25 * (L % 4) for L in lengths},
            {L: 0.5 if L < 150 else 2.0 for L in lengths}]


def test_petrackI_pileup_bdg_hmmr():
    pe = random_pe(12, chroms=(b"chr1", b"chr2"))
    mapping = hmmr_mapping(pe)
    res = pe.pileup_bdg_hmmr(mapping)
    assert len(res) == len(mapping)
    for m, by_chrom in zip(mapping, res):
        assert set(by_chrom) == {b"chr1", b"chr2"}
        for c, pv in by_chrom.items():
            locs = pe.get_locations_by_chr(c)
            w = [m[int(r - l)] for l, r in locs.tolist()]
            assert_pv_equal(pv, ref_pileup(locs['l'], locs['r'], w))


def test_petrackI_pileup_bdg_hmmr_ignores_baseline():
    pe = random_pe(12)
    mapping = hmmr_mapping(pe)[:1]
    a = pe.pileup_bdg_hmmr(mapping)
    b = pe.pileup_bdg_hmmr(mapping, baseline_value=5.0)
    np.testing.assert_array_equal(a[0][b"chr1"], b[0][b"chr1"])


def test_petrackI_pileup_bdg_hmmr_missing_length_raises():
    pe = build_pe([(b"chrA", 0, 100)])
    with pytest.raises(KeyError):
        pe.pileup_bdg_hmmr([{50: 1.0}])


# ------------------------------------
# PETrackII
# ------------------------------------

FRAGS2 = [(b"chrA", 0, 100, b"AAA", 2), (b"chrA", 70, 170, b"AAA", 1),
          (b"chrA", 80, 190, b"CCC", 1), (b"chrA", 85, 180, b"AAA", 3),
          (b"chrA", 100, 190, b"GGG", 1), (b"chrA", 0, 100, b"CCC", 1),
          (b"chrA", 70, 170, b"CCC", 2), (b"chrA", 500, 600, b"TTT", 4),
          (b"chrB", 10, 110, b"GGG", 1), (b"chrB", 50, 160, b"AAA", 2),
          (b"chrB", 100, 170, b"CCC", 3)]
# count totals: chrA 15, chrB 6; barcode TTT only on chrA


def build_pe2(frags, buffer_size=100000):
    pe = PETrackII(buffer_size=buffer_size)
    for c, l, r, b, n in frags:
        pe.add_loc(c, l, r, b, n)
    pe.finalize()
    return pe


def pe2_records(pe):
    """{chrom: sorted [(l, r, count, barcode), ...]} from a PETrackII."""
    names = {v: k for k, v in pe.barcode_dict.items()}
    out = {}
    for c in pe.get_chr_names():
        locs = pe.get_locations_by_chr(c)
        bars = pe.barcodes[c]
        out[c] = sorted((int(l), int(r), int(n), names[int(b)])
                        for (l, r, n), b in zip(locs.tolist(), bars.tolist()))
    return out


def expected_pe2(frags):
    out = {}
    for c, l, r, b, n in frags:
        out.setdefault(c, []).append((l, r, n, b))
    return {c: sorted(v) for c, v in out.items()}


def random_pe2(seed, n=80, span=4000, lo=0, chroms=(b"chr1",)):
    rng = np.random.default_rng(seed)
    pe = PETrackII()
    for chrom in chroms:
        starts = rng.integers(lo, span, n).tolist()
        lens = rng.integers(20, 400, n).tolist()
        counts = rng.integers(1, 6, n).tolist()
        bcs = rng.integers(0, 5, n).tolist()
        for s, ln, k, b in zip(starts, lens, counts, bcs):
            pe.add_loc(chrom, s, s + ln, b"BC%d" % b, k)
    pe.finalize()
    return pe


def test_petrackII_init_defaults():
    pe = PETrackII()
    assert (pe.total, pe.length, pe.annotation, pe.buffer_size) == \
        (0, 0, "", 100000)
    assert (pe.locations, pe.barcodes, pe.barcode_dict) == ({}, {}, {})
    assert pe.is_sorted is False


@pytest.mark.parametrize("buffer_size", [1, 2, 100000])
def test_petrackII_add_loc_finalize(buffer_size):
    pe = build_pe2(FRAGS2, buffer_size)
    assert pe2_records(pe) == expected_pe2(FRAGS2)
    # barcodes are numbered in order of first appearance
    assert pe.barcode_dict == {b"AAA": 0, b"CCC": 1, b"GGG": 2, b"TTT": 3}
    for c in pe.get_chr_names():
        locs = pe.get_locations_by_chr(c)
        assert locs.dtype == np.dtype([('l', 'i4'), ('r', 'i4'),
                                       ('c', 'u2')])
        assert locs[['l', 'r']].tolist() == sorted(locs[['l', 'r']].tolist())
    total = sum(n for *_, n in FRAGS2)
    length = sum((r - l) * n for _, l, r, _, n in FRAGS2)
    assert pe.total == total == 21
    assert pe.length == length
    assert pe.average_template_length == f32_mean(length, total)
    assert pe.is_sorted is True


@pytest.mark.parametrize("count", [1, 2, 1000, 65535])
def test_petrackII_add_loc_count_is_stored(count):
    pe = build_pe2([(b"chrA", 10, 30, b"B", count)])
    assert pe.get_locations_by_chr(b"chrA")['c'].tolist() == [count]
    assert pe.total == count
    assert pe.length == 20 * count


@pytest.mark.parametrize("count", [-1, 65536])
def test_petrackII_add_loc_count_outside_uint16_raises(count):
    with pytest.raises(OverflowError):
        PETrackII().add_loc(b"chrA", 10, 30, b"B", count)


def test_petrackII_add_loc_70000_is_not_silently_wrapped():
    # a Python caller cannot store a wrapped count: the conversion to the
    # uint16 count column is checked
    pe = PETrackII()
    with pytest.raises(OverflowError):
        pe.add_loc(b"chr1", 10, 100, b"BC1", 70000)


def test_petrackII_add_loc_int32_extreme():
    pe = build_pe2([(b"chrA", INT32_MAX - 50, INT32_MAX, b"B", 2)])
    assert pe2_records(pe) == {b"chrA": [(INT32_MAX - 50, INT32_MAX, 2,
                                          b"B")]}


def test_petrackII_add_loc_chromosome_and_barcode_must_be_bytes():
    with pytest.raises(TypeError):
        PETrackII().add_loc("chrA", 10, 30, b"B", 1)
    with pytest.raises(TypeError):
        PETrackII().add_loc(b"chrA", 10, 30, "B", 1)


def test_petrackII_rlengths():
    pe = build_pe2(FRAGS2)
    assert pe.get_rlengths() == {b"chrA": INT32_MAX, b"chrB": INT32_MAX}
    pe2 = build_pe2(FRAGS2)
    assert pe2.set_rlengths({b"chrB": 900, b"chrQ": 1}) is True
    assert pe2.get_rlengths() == {b"chrA": INT32_MAX, b"chrB": 900}


def test_petrackII_finalize_without_fragments_raises():
    with pytest.raises(AssertionError, match="no fragments in PETrackII"):
        PETrackII().finalize()


def test_petrackII_finalize_with_only_zero_counts_raises():
    pe = PETrackII()
    pe.add_loc(b"chrA", 10, 30, b"B", 0)
    with pytest.raises(AssertionError, match="no fragments in PETrackII"):
        pe.finalize()


def test_petrackII_get_locations_by_chr_missing_raises():
    with pytest.raises(Exception,
                       match=r"No such chromosome name \(b'chrZ'\)"):
        build_pe2(FRAGS2).get_locations_by_chr(b"chrZ")


def test_petrackII_get_chr_names():
    assert build_pe2(FRAGS2).get_chr_names() == {b"chrA", b"chrB"}


def test_petrackII_sort_keeps_barcodes_aligned():
    pe = build_pe2(FRAGS2)
    for c in pe.get_chr_names():
        order = np.random.default_rng(2).permutation(pe.size[c])
        pe.locations[c] = pe.locations[c][order]
        pe.barcodes[c] = pe.barcodes[c][order]
    pe.is_sorted = False
    pe.sort()
    assert pe.is_sorted is True
    assert pe2_records(pe) == expected_pe2(FRAGS2)
    for c in pe.get_chr_names():
        lr = pe.get_locations_by_chr(c)[['l', 'r']].tolist()
        assert lr == sorted(lr)


def test_petrackII_count_fraglengths():
    weighted = Counter()
    for _, l, r, _, n in FRAGS2:
        weighted[r - l] += n
    assert build_pe2(FRAGS2).count_fraglengths() == dict(weighted)


def test_petrackII_fraglengths():
    sizes = build_pe2(FRAGS2).fraglengths()
    assert sizes.dtype == np.int32
    assert sorted(sizes.tolist()) == sorted(
        r - l for _, l, r, _, n in FRAGS2 for _ in range(n))


def test_petrackII_fraglengths_empty_track():
    sizes = PETrackII().fraglengths()
    assert sizes.dtype == np.int32
    assert sizes.shape == (0,)


@pytest.mark.parametrize("selected", [{b"AAA"}, {b"CCC", b"GGG"},
                                      {b"AAA", b"CCC", b"GGG", b"TTT"},
                                      {b"TTT", b"NOPE"}])
def test_petrackII_subset(selected):
    pe = build_pe2(FRAGS2)
    sub = pe.subset(selected)
    kept = [f for f in FRAGS2 if f[3] in selected]
    assert pe2_records(sub) == expected_pe2(kept)
    assert sub.total == sum(f[4] for f in kept)
    assert sub.length == sum((f[2] - f[1]) * f[4] for f in kept)
    assert sub.barcode_dict == {b: pe.barcode_dict[b] for b in selected
                                if b in pe.barcode_dict}
    assert pe2_records(pe) == expected_pe2(FRAGS2)     # source untouched


def test_petrackII_subset_without_matching_barcodes_raises():
    with pytest.raises(AssertionError, match="no fragments in PETrackII"):
        build_pe2(FRAGS2).subset({b"NOPE"})


def test_petrackII_pileup_a_chromosome_hand_example():
    # [0, 100) x2 and [50, 150) x1: [0, 50) 2, [50, 100) 3, [100, 150) 1
    pe = build_pe2([(b"chrA", 0, 100, b"A", 2), (b"chrA", 50, 150, b"B", 1)])
    p, v = pe.pileup_a_chromosome(b"chrA")
    np.testing.assert_array_equal(p, [50, 100, 150])
    np.testing.assert_array_equal(v, [2.0, 3.0, 1.0])


@pytest.mark.parametrize("scale, baseline", [(1.0, 0.0), (0.5, 0.0),
                                             (2.0, 3.0), (0.1, 0.05)])
def test_petrackII_pileup_a_chromosome_matches_reference(scale, baseline):
    pe = random_pe2(3)
    p, v = pe.pileup_a_chromosome(b"chr1", scale_factor=scale,
                                  baseline_value=baseline)
    locs = pe.get_locations_by_chr(b"chr1")
    assert_cp_equal(p, v, ref_pileup(locs['l'], locs['r'], locs['c'],
                                     scale, baseline))


def test_petrackII_pileup_a_chromosome_large_count():
    pe = build_pe2([(b"chrA", 10, 20, b"A", 65535),
                    (b"chrA", 15, 25, b"B", 65535)])
    p, v = pe.pileup_a_chromosome(b"chrA")
    assert_cp_equal(p, v, ([10, 15, 20, 25],
                           [0.0, 65535.0, 131070.0, 65535.0]))


def test_petrackII_pileup_a_chromosome_zero_fragments():
    pe = build_pe2(FRAGS2 + [(b"chrC", 10, 60, b"AAA", 1)])
    pe.sample_percent(0.4, seed=1)      # chrC: round(1 * 0.4) -> 0
    p, v = pe.pileup_a_chromosome(b"chrC")
    assert p.shape == (0,) and v.shape == (0,)


def test_petrackII_pileup_a_chromosome_c_single_d():
    # each fragment end x contributes count x [x - d//2, x - d//2 + d);
    # all fragments start at >= 500 so no window reaches position 0
    pe = random_pe2(5, lo=500)
    p, v = pe.pileup_a_chromosome_c(b"chr1", [200], [0.5],
                                    baseline_value=0.25)
    locs = pe.get_locations_by_chr(b"chr1")
    starts, ends = end_windows(locs, 100, 200)
    w = np.concatenate((locs['c'], locs['c']))
    assert_cp_equal(p, v, ref_pileup(starts, ends, w, 0.5, 0.25))


def test_petrackII_pileup_a_chromosome_c_length_mismatch():
    with pytest.raises(AssertionError, match="same length"):
        random_pe2(5).pileup_a_chromosome_c(b"chr1", [100, 200], [1.0])


def test_petrackII_pileup_a_chromosome_c_several_ds():
    """Regression test: with several ds, pileup_a_chromosome_c passed the
    strided 'p' field of a structured pileup to over_two_pv_array, whose
    raw int pointer read alternating position and value bits.

    Fixed upstream in 5456b02 (#737), which builds each pileup with
    pileup_from_LRC_centers_as_list (contiguous arrays).
    """
    pe = random_pe2(7, lo=600)
    ds, sfs = [200, 1000], [1.0, 0.5]
    p, v = pe.pileup_a_chromosome_c(b"chr1", ds, sfs)
    locs = pe.get_locations_by_chr(b"chr1")
    w = np.concatenate((locs['c'], locs['c']))
    cps = []
    for d, sf in zip(ds, sfs):
        starts, ends = end_windows(locs, d // 2, d)
        cps.append(ref_pileup(starts, ends, w, sf, 0.0))
    assert_cp_equal(p, v, ref_max_truncated(cps))


@pytest.mark.parametrize("scale, baseline", [(1.0, 0.0), (0.5, 0.25)])
def test_petrackII_pileup_bdg(scale, baseline):
    pe = random_pe2(8, chroms=(b"chr1", b"chr2"))
    bdg = pe.pileup_bdg(scale_factor=scale, baseline_value=baseline)
    assert bdg.get_chr_names() == {b"chr1", b"chr2"}
    for c in pe.get_chr_names():
        p, v = bdg.get_data_by_chr(c)
        locs = pe.get_locations_by_chr(c)
        assert_cp_equal(p, v, ref_pileup(locs['l'], locs['r'], locs['c'],
                                         scale, baseline))


def test_petrackII_pileup_bdg2():
    pe = random_pe2(8, chroms=(b"chr1", b"chr2"))
    bdg = pe.pileup_bdg2()
    assert bdg.get_chr_names() == {b"chr1", b"chr2"}
    for c in pe.get_chr_names():
        locs = pe.get_locations_by_chr(c)
        assert_pv_equal(bdg.get_data_by_chr(c),
                        ref_pileup(locs['l'], locs['r'], locs['c']))


def check_pe2_excluded(pe, frags, regions):
    kept = [f for f in frags if not overlaps(f[1], f[2],
                                             regions.get(f[0], []))]
    assert pe2_records(pe) == expected_pe2(kept)
    total = sum(f[4] for f in kept)
    length = sum((f[2] - f[1]) * f[4] for f in kept)
    assert pe.total == total
    assert pe.length == length
    assert pe.average_template_length == f32_mean(length, total)


@pytest.mark.parametrize("regions", [
    {b"chrA": [(10, 75)]},
    {b"chrA": [(60, 75), (550, 560)], b"chrB": [(105, 106)]},
    {b"chrB": [(0, 10)]},
    {b"chrZ": [(0, 1000)]},
])
def test_petrackII_exclude(regions):
    pe = build_pe2(FRAGS2)
    pe.exclude(make_regions(regions))
    check_pe2_excluded(pe, FRAGS2, regions)


def check_pe2_sampled(pe, frags, percent):
    """Per chromosome the kept counts sum to round(n * percent) (float32
    product), every kept fragment exists in the input with at least its
    new count, and the arrays are sorted."""
    orig = {(c, l, r, b): n for c, l, r, b, n in frags}
    totals = Counter()
    for c, l, r, b, n in frags:
        totals[c] += n
    got = pe2_records(pe)
    total = length = 0
    for c, n in totals.items():
        want = round(float(np.float32(n) * np.float32(percent)))
        recs = got.get(c, [])
        assert sum(k for _, _, k, _ in recs) == want
        for l, r, k, b in recs:
            assert 0 < k <= orig[(c, l, r, b)]
        lr = pe.get_locations_by_chr(c)[['l', 'r']].tolist() if recs else []
        assert lr == sorted(lr)
        total += want
        length += sum((r - l) * k for l, r, k, _ in recs)
    assert pe.total == total
    assert pe.length == length


@pytest.mark.parametrize("percent", [0.0, 0.25, 0.35, 0.6, 1.0])
def test_petrackII_sample_percent(percent):
    pe = build_pe2(FRAGS2)
    pe.sample_percent(percent, seed=4)
    check_pe2_sampled(pe, FRAGS2, percent)


def test_petrackII_sample_percent_without_seed():
    pe = build_pe2(FRAGS2)
    pe.sample_percent(0.6)
    check_pe2_sampled(pe, FRAGS2, 0.6)


@pytest.mark.parametrize("seed", [0, 31])
def test_petrackII_sample_percent_same_seed_same_sample(seed):
    a, b = build_pe2(FRAGS2), build_pe2(FRAGS2)
    a.sample_percent(0.5, seed=seed)
    b.sample_percent(0.5, seed=seed)
    assert pe2_records(a) == pe2_records(b)


@pytest.mark.parametrize("percent", [-0.1, 1.5])
def test_petrackII_sample_percent_out_of_range(percent):
    with pytest.raises(AssertionError, match=r"percent must be in \[0, 1\]"):
        build_pe2(FRAGS2).sample_percent(percent)


def test_petrackII_sample_percent_logs_seed(caplog):
    caplog.set_level(logging.INFO)
    build_pe2(FRAGS2).sample_percent(0.5, seed=6)
    assert any(r.getMessage().endswith("#   A random seed 6 has been used")
               for r in caplog.records)


@pytest.mark.parametrize("percent", [0.25, 0.6, 1.0])
def test_petrackII_sample_percent_copy(percent):
    pe = build_pe2(FRAGS2)
    pe.set_rlengths({b"chrA": 800})
    sub = pe.sample_percent_copy(percent, seed=4)
    assert isinstance(sub, PETrackII)
    assert pe2_records(pe) == expected_pe2(FRAGS2)     # source untouched
    assert pe.total == 21
    check_pe2_sampled(sub, FRAGS2, percent)
    assert sub.get_rlengths() == {b"chrA": 800, b"chrB": INT32_MAX}


def test_petrackII_sample_percent_copy_same_seed_same_sample():
    pe = build_pe2(FRAGS2)
    a = pe.sample_percent_copy(0.5, seed=12)
    b = pe.sample_percent_copy(0.5, seed=12)
    assert pe2_records(a) == pe2_records(b)


def test_petrackII_sample_percent_copy_out_of_range():
    with pytest.raises(AssertionError, match=r"percent must be in \[0, 1\]"):
        build_pe2(FRAGS2).sample_percent_copy(2.0)


@pytest.mark.parametrize("samplesize, percent", [(7, 7 / 21), (10, 10 / 21),
                                                 (21, 1.0), (100, 1.0)])
def test_petrackII_sample_num(samplesize, percent):
    # percent = min(samplesize / total count, 1)
    pe = build_pe2(FRAGS2)
    pe.sample_num(samplesize, seed=3)
    check_pe2_sampled(pe, FRAGS2, percent)


@pytest.mark.parametrize("samplesize, percent", [(7, 7 / 21), (10, 10 / 21),
                                                 (100, 1.0)])
def test_petrackII_sample_num_copy(samplesize, percent):
    pe = build_pe2(FRAGS2)
    sub = pe.sample_num_copy(samplesize, seed=3)
    assert pe.total == 21
    check_pe2_sampled(sub, FRAGS2, percent)


def test_petrackII_pileup_bdg_hmmr():
    pe = random_pe2(10, chroms=(b"chr1", b"chr2"))
    lengths = sorted(pe.count_fraglengths())
    mapping = [{L: 1.0 for L in lengths},
               {L: 0.25 * (L % 4) for L in lengths}]
    res = pe.pileup_bdg_hmmr(mapping)
    assert len(res) == 2
    for m, by_chrom in zip(mapping, res):
        assert set(by_chrom) == {b"chr1", b"chr2"}
        for c, pv in by_chrom.items():
            locs = pe.get_locations_by_chr(c)
            w = [n * m[int(r - l)] for l, r, n in locs.tolist()]
            assert_pv_equal(pv, ref_pileup(locs['l'], locs['r'], w))


def test_petrackII_pileup_bdg_hmmr_empty_chromosome():
    pe = build_pe2(FRAGS2 + [(b"chrC", 10, 60, b"AAA", 1)])
    pe.sample_percent(0.4, seed=1)      # chrC is left without fragments
    lengths = sorted(Counter(f[2] - f[1] for f in FRAGS2))
    res = pe.pileup_bdg_hmmr([{L: 1.0 for L in lengths}])
    assert res[0][b"chrC"].shape == (0,)


if __name__ == '__main__':
    unittest.main()
