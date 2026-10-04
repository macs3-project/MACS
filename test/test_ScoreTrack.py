#!/usr/bin/env python
# Time-stamp: <2024-10-18 15:30:21 Tao Liu>

import io
import unittest
from numpy.testing import assert_array_equal  # assert_equal, assert_almost_equal

import numpy as np
from MACS3.Signal.ScoreTrack import ScoreTrackII, TwoConditionScores
from MACS3.Signal.BedGraph import bedGraphTrackI

import math

import pytest
from scipy.stats import poisson

from MACS3.IO.PeakIO import (PeakIO,
                             BroadPeakIO)


def make_two_summit_profile(length=260, right_apex_distance=15):
    """Return a peak profile with a summit close to the right boundary."""
    signal = np.full(length, 10.0)
    half_width = 70
    for center in (60, length - right_apex_distance):
        positions = np.arange(max(center - half_width + 1, 0),
                              min(center + half_width, length))
        signal[positions] = np.maximum(
            signal[positions],
            10 + 30 * (1 - np.abs(positions - center) / half_width),
        )
    return np.rint(signal).astype(int)


def make_score_track_from_profile(signal, start):
    """Encode a coverage profile as a finalized fold-enrichment track."""
    score_track = ScoreTrackII(1, 1)
    score_track.add_chromosome(b"chrSynthetic", len(signal) + 1)
    score_track.add(b"chrSynthetic", start, 0, 1)
    for offset, value in enumerate(signal, 1):
        score_track.add(b"chrSynthetic", start + offset, value, 1)
    score_track.finalize()
    score_track.change_score_method(ord("F"))
    return score_track


class Test_TwoConditionScores(unittest.TestCase):
    def setUp(self):
        self.t1bdg = bedGraphTrackI()
        self.t2bdg = bedGraphTrackI()
        self.c1bdg = bedGraphTrackI()
        self.c2bdg = bedGraphTrackI()
        self.test_regions1 = [(b"chrY", 0, 70, 0.00, 0.01),
                              (b"chrY", 70, 80, 7.00, 0.5),
                              (b"chrY", 80, 150, 0.00, 0.02)]
        self.test_regions2 = [(b"chrY", 0, 75, 20.0, 4.00),
                              (b"chrY", 75, 90, 35.0, 6.00),
                              (b"chrY", 90, 150, 10.0, 15.00)]
        for a in self.test_regions1:
            self.t1bdg.add_loc(a[0], a[1], a[2], a[3])
            self.c1bdg.add_loc(a[0], a[1], a[2], a[4])

        for a in self.test_regions2:
            self.t2bdg.add_loc(a[0], a[1], a[2], a[3])
            self.c2bdg.add_loc(a[0], a[1], a[2], a[4])

        self.twoconditionscore = TwoConditionScores(self.t1bdg,
                                                    self.c1bdg,
                                                    self.t2bdg,
                                                    self.c2bdg,
                                                    1.0,
                                                    1.0)
        self.twoconditionscore.build()
        self.twoconditionscore.finalize()
        (self.cat1, self.cat2, self.cat3) = self.twoconditionscore.call_peaks(min_length=10, max_gap=10, cutoff=3)

    def test_call_peaks_score_is_mean_logLR(self):
        # Regression test for issue #715: the score of a bdgdiff region
        # must be the length-weighted mean of the log10 likelihood
        # ratios of its intervals, not the mean of their integer parts.
        # Two adjacent cond1 intervals of equal length with
        # log10LR(t1 vs t2) = 0.3271 and 6.1784 (pseudocount 0.01)
        # must give (0.3271 + 6.1784) / 2 = 3.2527, whereas truncating
        # each value to an integer gives (0 + 6) / 2 = 3.0.
        regions = [(0, 1000, 20.0, 2.0, 18.0, 2.0),
                   (1000, 2000, 20.0, 2.0, 15.0, 2.0),
                   (2000, 3000, 40.0, 2.0, 15.0, 2.0),
                   (3000, 4000, 2.0, 2.0, 2.0, 2.0)]
        bdgs = [bedGraphTrackI() for i in range(4)]
        for (s, e, t1, c1, t2, c2) in regions:
            for (bdg, v) in zip(bdgs, (t1, c1, t2, c2)):
                bdg.add_loc(b"chr1", s, e, v)
        tcs = TwoConditionScores(bdgs[0], bdgs[1], bdgs[2], bdgs[3], 1.0, 1.0)
        tcs.build()
        tcs.finalize()
        (cat1, cat2, _cat3) = tcs.call_peaks(min_length=200, max_gap=100, cutoff=0.2)
        self.assertEqual(cat1.total, 1)
        self.assertEqual(cat2.total, 0)
        peaks = cat1.get_data_from_chrom(b"chr1")
        self.assertAlmostEqual(peaks[0]["score"], 3.2527, places=3)


class Test_ScoreTrackII(unittest.TestCase):

    def setUp(self):
        # for initiate scoretrack
        self.test_regions1 = [(b"chrY", 10, 100, 10),
                              (b"chrY", 60, 10, 10),
                              (b"chrY", 110, 15, 20),
                              (b"chrY", 160, 5, 20),
                              (b"chrY", 210, 20, 5)]
        self.treat_edm = 10
        self.ctrl_edm = 5
        # for different scoring method
        self.p_result = [60.49, 0.38, 0.08, 0.0, 6.41]  # -log10(p-value), pseudo count 1 added
        self.q_result = [58.17, 0.0, 0.0, 0.0, 5.13]  # -log10(q-value) from BH, pseudo count 1 added
        self.l_result = [58.17, 0.0, -0.28, -3.25, 4.91]  # log10 likelihood ratio, pseudo count 1 added
        self.f_result = [0.96, 0.00, -0.12, -0.54, 0.54]  # note, pseudo count 1 would be introduced.
        self.d_result = [90.00, 0, -5.00, -15.00, 15.00]
        self.m_result = [10.00, 1.00, 1.50, 0.50, 2.00]
        # for norm
        self.norm_T = np.array([[ 10, 100,  20,   0],
                                [ 60,  10,  20,   0],
                                [110,  15,  40,   0],
                                [160,   5,  40,   0],
                                [210,  20,  10,   0]]).transpose()
        self.norm_C = np.array([[ 10,  50,  10,   0],
                                [ 60,   5,  10,   0],
                                [110,   7.5,  20,   0],
                                [160,   2.5,  20,   0],
                                [210,  10,   5,   0]]).transpose()
        self.norm_M = np.array([[ 10,  10,   2,   0],
                                [ 60,   1,   2,   0],
                                [110,   1.5,   4,   0],
                                [160,   0.5,   4,   0],
                                [210,   2,   1,   0]]).transpose()
        self.norm_N = np.array([[ 10, 100,  10,   0],  # note precision lost
                                [ 60,  10,  10,   0],
                                [110,  15,  20,   0],
                                [160,   5,  20,   0],
                                [210,  20,   5,   0]]).transpose()

        # for write_bedGraph
        self.bdg1 = """chrY	0	10	100.00000
chrY	10	60	10.00000
chrY	60	110	15.00000
chrY	110	160	5.00000
chrY	160	210	20.00000
"""
        self.bdg2 = """chrY	0	60	10.00000
chrY	60	160	20.00000
chrY	160	210	5.00000
"""
        self.bdg3 = """chrY	0	10	60.48912
chrY	10	60	0.37599
chrY	60	110	0.07723
chrY	110	160	0.00006
chrY	160	210	6.40804
"""
        # for peak calls
        self.peak1 = """chrY	0	60	MACS_peak_1	60.4891
chrY	160	210	MACS_peak_2	6.40804
"""
        self.summit1 = """chrY	5	6	MACS_peak_1	60.4891
chrY	185	186	MACS_peak_2	6.40804
"""
        self.xls1    ="""chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name
chrY	1	60	60	6	100	63.2725	9.18182	-1	MACS_peak_1
chrY	161	210	50	186	20	7.09102	3.5	-1	MACS_peak_2
"""

    def assertListAlmostEqual(self, a, b, places=2):
        return all([self.assertAlmostEqual(x, y, places=places) for (x, y) in zip(a, b)])

    def test_compute_scores(self):
        s1 = ScoreTrackII(self.treat_edm, self.ctrl_edm)
        s1.add_chromosome(b"chrY", 5)
        for a in self.test_regions1:
            s1.add(a[0], a[1], a[2], a[3])

        s1.set_pseudocount(1.0)

        s1.change_score_method(ord('p'))
        r = s1.get_data_by_chr(b"chrY")
        self.assertListAlmostEqual([round(x, 2) for x in r[3]], self.p_result)

        s1.change_score_method(ord('q'))
        r = s1.get_data_by_chr(b"chrY")
        self.assertListAlmostEqual([round(x, 2) for x in list(r[3])], self.q_result)

        s1.change_score_method(ord('l'))
        r = s1.get_data_by_chr(b"chrY")
        self.assertListAlmostEqual([round(x, 2) for x in list(r[3])], self.l_result)

        s1.change_score_method(ord('f'))
        r = s1.get_data_by_chr(b"chrY")
        self.assertListAlmostEqual([round(x, 2) for x in list(r[3])], self.f_result)

        s1.change_score_method(ord('d'))
        r = s1.get_data_by_chr(b"chrY")
        self.assertListAlmostEqual([round(x, 2) for x in list(r[3])], self.d_result)

        s1.change_score_method(ord('m'))
        r = s1.get_data_by_chr(b"chrY")
        self.assertListAlmostEqual([round(x, 2) for x in list(r[3])], self.m_result)

    def test_normalize(self):
        s1 = ScoreTrackII(self.treat_edm, self.ctrl_edm)
        s1.add_chromosome(b"chrY", 5)
        for a in self.test_regions1:
            s1.add(a[0], a[1], a[2], a[3])

        s1.change_normalization_method(ord('T'))
        r = s1.get_data_by_chr(b"chrY")
        assert_array_equal(r, self.norm_T)

        s1.change_normalization_method(ord('C'))
        r = s1.get_data_by_chr(b"chrY")
        assert_array_equal(r, self.norm_C)

        s1.change_normalization_method(ord('M'))
        r = s1.get_data_by_chr(b"chrY")
        assert_array_equal(r, self.norm_M)

        s1.change_normalization_method(ord('N'))
        r = s1.get_data_by_chr(b"chrY")
        assert_array_equal(r, self.norm_N)

    def test_writebedgraph(self):
        s1 = ScoreTrackII(self.treat_edm, self.ctrl_edm)
        s1.add_chromosome(b"chrY", 5)
        for a in self.test_regions1:
            s1.add(a[0], a[1], a[2], a[3])

        s1.change_score_method(ord('p'))

        strio = io.StringIO()
        s1.write_bedGraph(strio, "NAME", "DESC", 1)
        self.assertEqual(strio.getvalue(), self.bdg1)
        strio = io.StringIO()
        s1.write_bedGraph(strio, "NAME", "DESC", 2)
        self.assertEqual(strio.getvalue(), self.bdg2)
        strio = io.StringIO()
        s1.write_bedGraph(strio, "NAME", "DESC", 3)
        self.assertEqual(strio.getvalue(), self.bdg3)

    def test_callpeak(self):
        s1 = ScoreTrackII(self.treat_edm, self.ctrl_edm)
        s1.add_chromosome(b"chrY", 5)
        for a in self.test_regions1:
            s1.add(a[0], a[1], a[2], a[3])

        s1.change_score_method(ord('p'))
        p = s1.call_peaks(cutoff=0.10, min_length=10, max_gap=10)
        strio = io.StringIO()
        p.write_to_bed(strio, trackline=False)
        self.assertEqual(strio.getvalue(), self.peak1)

        strio = io.StringIO()
        p.write_to_summit_bed(strio, trackline=False)
        self.assertEqual(strio.getvalue(), self.summit1)

        strio = io.StringIO()
        p.write_to_xls(strio)
        self.assertEqual(strio.getvalue(), self.xls1)

    def test_call_summits_rejects_maximum_in_below_cutoff_gap(self):
        signal = np.rint(np.interp(np.arange(123),
                                   np.linspace(0, 122, 8),
                                   [27, 1, 27, 8,
                                    14, 10, 27, 12])).astype(int)
        score_track = ScoreTrackII(1, 1)
        score_track.add_chromosome(b"chrSynthetic", len(signal) + 1)
        score_track.add(b"chrSynthetic", 1000, 0, 1)
        for offset, value in enumerate(signal, 1):
            score_track.add(b"chrSynthetic", 1000 + offset, value, 1)
        score_track.finalize()
        score_track.change_score_method(ord("F"))

        peaks = score_track.call_peaks(cutoff=9, min_length=75,
                                       max_gap=150, call_summits=True)

        peak = peaks.peaks[b"chrSynthetic"][0]
        actual_pileup = signal[peak["summit"] - 1000]
        self.assertEqual(peak["summit"], 1035)
        self.assertEqual(peak["pileup"], actual_pileup)
        self.assertEqual(actual_pileup, 27)

    def test_call_summits_keeps_right_edge_candidate(self):
        """ScoreTrackII uses the same corrected coordinates as callpeak."""
        signal = make_two_summit_profile()

        for start in (0, 1, 5, 9, 10, 1000):
            score_track = make_score_track_from_profile(signal, start)
            peaks = score_track.call_peaks(cutoff=2, min_length=50,
                                           max_gap=100,
                                           call_summits=True)
            rows = peaks.peaks[b"chrSynthetic"]

            relative_summits = [row["summit"] - start for row in rows]
            self.assertEqual(len(relative_summits), 2)
            self.assertLessEqual(abs(relative_summits[0] - 60), 1)
            self.assertLessEqual(abs(relative_summits[1] - 236), 1)
            self.assertTrue(all(row["start"] == start for row in rows))
            self.assertTrue(all(row["end"] == start + len(signal)
                                for row in rows))
            self.assertTrue(all(row["start"] <= row["summit"] < row["end"]
                                for row in rows))


# ------------------------------------
# Helpers for the tests below
# ------------------------------------
#
# A ScoreTrackII holds, per chromosome, four arrays: region end
# positions, treatment pileup, control pileup and score; region i spans
# [pos[i-1], pos[i]) with pos[-1] taken as 0.

LN10 = math.log(10)


def make_st(rows, treat_depth=1.0, ctrl_depth=1.0, pseudocount=1.0,
            chrom=b"chrY"):
    """Finalized ScoreTrackII from (end, treat, ctrl) rows."""
    s = ScoreTrackII(treat_depth, ctrl_depth, pseudocount)
    s.add_chromosome(chrom, len(rows))
    for (end, t, c) in rows:
        s.add(chrom, end, t, c)
    s.finalize()
    return s


def make_st_multi(chrom_rows, treat_depth=1.0, ctrl_depth=1.0,
                  pseudocount=1.0):
    s = ScoreTrackII(treat_depth, ctrl_depth, pseudocount)
    for chrom, rows in chrom_rows.items():
        s.add_chromosome(chrom, len(rows))
        for (end, t, c) in rows:
            s.add(chrom, end, t, c)
    s.finalize()
    return s


def column(s, i, chrom=b"chrY"):
    return s.get_data_by_chr(chrom)[i].tolist()


def f32(x):
    return float(np.float32(x))


def ref_pscore(observed, lam):
    """-log10 P(X > observed) for X ~ Poisson(lam)."""
    return -poisson.logsf(observed, lam) / LN10


def poisson_llr(x, y):
    """Poisson log-likelihood ratio x*ln(x/y) + y - x in log10 units."""
    return (x * math.log(x / y) + y - x) / LN10


def ref_logLR_asym(x, y):
    # signed by the direction of enrichment of x over y
    if x > y:
        return poisson_llr(x, y)
    if x < y:
        return -poisson_llr(x, y)
    return 0.0


def ref_logLR_sym(x, y):
    # LLR of the larger value over the smaller one, negative when y > x
    if x > y:
        return poisson_llr(x, y)
    if y > x:
        return -poisson_llr(y, x)
    return 0.0


def bh_qscores(pscore_lengths):
    """Benjamini-Hochberg over base pairs as documented in make_pq_table:
    {pscore: qscore} for a {pscore: total length} dict. A p-score v
    whose block starts at rank k (1 + bp of all higher p-scores) out of
    N bp gets v + log10(k / N), made monotone and clipped at 0."""
    n = sum(pscore_lengths.values())
    k = 1
    pre_q = math.inf
    out = {}
    for v in sorted(pscore_lengths, reverse=True):
        q = max(min(pre_q, v + math.log10(k / n)), 0.0)
        out[v] = q
        pre_q = q
        k += pscore_lengths[v]
    return out


def region_lengths(pos):
    return [e - s for (s, e) in zip([0] + list(pos[:-1]), pos)]


def peak_fields(peakio):
    """(chrom, start, end, summit, score, pileup, pscore, fc, qscore)."""
    rows = []
    for chrom in sorted(peakio.get_chr_names()):
        for p in peakio.get_data_from_chrom(chrom):
            rows.append((chrom, p["start"], p["end"], p["summit"],
                         p["score"], p["pileup"], p["pscore"], p["fc"],
                         p["qscore"]))
    return rows


# Rows for the score-method tests: (end, treat, ctrl). Treat 2.7 checks
# the integer truncation of the observed count in the p-score.
SCORE_ROWS = [(10, 0.0, 1.0), (20, 3.0, 1.0), (30, 2.7, 5.0),
              (40, 10.0, 2.0), (50, 5.0, 5.0), (60, 1.0, 4.0)]


# ------------------------------------
# ScoreTrackII: construction, add_chromosome, add, finalize,
# get_data_by_chr, get_chr_names
# ------------------------------------

def test_ScoreTrackII_new_track_is_empty():
    s = ScoreTrackII(1.0, 1.0)
    assert s.get_chr_names() == set()
    assert s.get_data_by_chr(b"chrY") is None


def test_ScoreTrackII_add_chromosome_allocates_zero_arrays():
    s = ScoreTrackII(1.0, 1.0)
    assert s.add_chromosome(b"chr1", 3) is None
    (pos, t, c, v) = s.get_data_by_chr(b"chr1")
    assert [a.dtype for a in (pos, t, c, v)] == \
        [np.int32, np.float32, np.float32, np.float32]
    assert [a.tolist() for a in (pos, t, c, v)] == [[0, 0, 0]] + \
        [[0.0, 0.0, 0.0]] * 3


def test_ScoreTrackII_add_chromosome_twice_keeps_data():
    s = ScoreTrackII(1.0, 1.0)
    s.add_chromosome(b"chr1", 2)
    s.add(b"chr1", 10, 1.0, 2.0)
    s.add_chromosome(b"chr1", 5)
    s.add(b"chr1", 20, 3.0, 4.0)
    s.finalize()
    assert column(s, 0, b"chr1") == [10, 20]
    assert column(s, 1, b"chr1") == [1.0, 3.0]
    assert column(s, 2, b"chr1") == [2.0, 4.0]


def test_ScoreTrackII_add_stores_float32():
    s = make_st([(10, 0.1, 1e39), (2 ** 31 - 1, -2.5, 3.4e38)])
    assert column(s, 0) == [10, 2 ** 31 - 1]
    assert column(s, 1) == [f32(0.1), -2.5]
    assert column(s, 2) == [math.inf, f32(3.4e38)]


def test_ScoreTrackII_add_unknown_chromosome_raises():
    s = ScoreTrackII(1.0, 1.0)
    s.add_chromosome(b"chr1", 1)
    with pytest.raises(KeyError, match="chr2"):
        s.add(b"chr2", 10, 1.0, 1.0)


def test_ScoreTrackII_add_beyond_capacity_raises():
    s = ScoreTrackII(1.0, 1.0)
    s.add_chromosome(b"chr1", 1)
    s.add(b"chr1", 10, 1.0, 1.0)
    with pytest.raises(IndexError, match="index 1 is out of bounds"):
        s.add(b"chr1", 20, 1.0, 1.0)


def test_ScoreTrackII_add_end_beyond_int32_raises():
    s = ScoreTrackII(1.0, 1.0)
    s.add_chromosome(b"chr1", 1)
    with pytest.raises(OverflowError, match="value too large to convert to int"):
        s.add(b"chr1", 2 ** 31, 1.0, 1.0)


def test_ScoreTrackII_finalize_trims_to_added_rows():
    s = ScoreTrackII(1.0, 1.0)
    s.add_chromosome(b"chr1", 5)
    s.add_chromosome(b"chr2", 4)
    s.add(b"chr1", 10, 1.0, 2.0)
    s.add(b"chr1", 20, 3.0, 4.0)
    assert s.finalize() is None
    assert [len(a) for a in s.get_data_by_chr(b"chr1")] == [2, 2, 2, 2]
    assert [len(a) for a in s.get_data_by_chr(b"chr2")] == [0, 0, 0, 0]
    assert s.get_chr_names() == {b"chr1", b"chr2"}


# ------------------------------------
# ScoreTrackII.set_pseudocount
# ------------------------------------

def test_ScoreTrackII_set_pseudocount():
    # linear fold enrichment (t + pc) / (c + pc) for t=3, c=1
    s = make_st([(10, 3.0, 1.0)], pseudocount=1.0)
    s.change_score_method(ord("F"))
    assert column(s, 3) == [2.0]
    s.set_pseudocount(0.0)
    s.change_score_method(ord("F"))
    assert column(s, 3) == [3.0]
    s.set_pseudocount(3.0)
    s.change_score_method(ord("F"))
    assert column(s, 3) == [1.5]


# ------------------------------------
# ScoreTrackII.change_normalization_method (covers normalize)
# ------------------------------------
#
# treat_depth 4, ctrl_depth 2; raw treat [8, 16], raw ctrl [8, 4].
# Each state's expected arrays, from the documented meaning:
#   N raw; T control scaled to treatment depth (ctrl * 4/2);
#   C treatment scaled to control depth (treat * 2/4);
#   M both per million reads (treat / 4, ctrl / 2).
NORM_STATE = {"N": ([8.0, 16.0], [8.0, 4.0]),
              "T": ([8.0, 16.0], [16.0, 8.0]),
              "C": ([4.0, 8.0], [8.0, 4.0]),
              "M": ([2.0, 4.0], [4.0, 2.0])}


# Switching back to 'N' from 'T' or 'C' is left out: in this version it
# multiplies by the depth instead of undoing the scaling.
@pytest.mark.parametrize("first, second",
                         [(a, b) for a in "NTCM" for b in "NTCM"
                          if (a, b) not in (("T", "N"), ("C", "N"))])
def test_ScoreTrackII_change_normalization_method(first, second):
    s = make_st([(10, 8.0, 8.0), (20, 16.0, 4.0)], treat_depth=4.0,
                ctrl_depth=2.0)
    s.change_normalization_method(ord(first))
    assert (column(s, 1), column(s, 2)) == NORM_STATE[first]
    s.change_normalization_method(ord(second))
    assert (column(s, 1), column(s, 2)) == NORM_STATE[second]


def test_ScoreTrackII_change_normalization_method_unknown_code():
    s = make_st([(10, 8.0, 8.0)], treat_depth=4.0, ctrl_depth=2.0)
    try:
        s.change_normalization_method(ord("X"))
    except NotImplementedError:
        pass
    assert (column(s, 1), column(s, 2)) == ([8.0], [8.0])


def test_ScoreTrackII_change_normalization_method_str_code_raises():
    s = make_st([(10, 8.0, 8.0)])
    with pytest.raises(TypeError, match="an integer is required"):
        s.change_normalization_method("T")


# ------------------------------------
# ScoreTrackII.change_score_method (covers compute_pvalue,
# compute_qvalue, compute_likelihood, compute_sym_likelihood,
# compute_logFE, compute_foldenrichment, compute_subtraction,
# compute_SPMR, compute_max and the module functions get_pscore,
# logLR_asym, logLR_sym, get_logFE, get_subtraction)
# ------------------------------------

@pytest.mark.parametrize("pseudocount", [0.0, 1.0, 2.5])
def test_ScoreTrackII_score_p(pseudocount):
    # -log10 Poisson upper tail with observed int(treat + pc) and
    # lambda ctrl + pc
    s = make_st(SCORE_ROWS, pseudocount=pseudocount)
    s.change_score_method(ord("p"))
    expected = [ref_pscore(int(t + pseudocount), c + pseudocount)
                for (_, t, c) in SCORE_ROWS]
    assert column(s, 3) == pytest.approx(expected, abs=1e-4)


@pytest.mark.parametrize("from_p", [False, True])
def test_ScoreTrackII_score_q(from_p):
    s = make_st(SCORE_ROWS)
    if from_p:
        s.change_score_method(ord("p"))
    s.change_score_method(ord("q"))
    pscores = [ref_pscore(int(t + 1.0), c + 1.0) for (_, t, c) in SCORE_ROWS]
    stat = {}
    for (v, ln) in zip(pscores, region_lengths([r[0] for r in SCORE_ROWS])):
        stat[v] = stat.get(v, 0) + ln
    q = bh_qscores(stat)
    # the lower p-scores reach q = 0, so the table ends with a cut-off
    assert sorted(q.values())[0] == 0.0
    assert column(s, 3) == pytest.approx([q[v] for v in pscores], abs=1e-4)


@pytest.mark.parametrize("pseudocount", [0.5, 1.0])
def test_ScoreTrackII_score_l(pseudocount):
    s = make_st(SCORE_ROWS, pseudocount=pseudocount)
    s.change_score_method(ord("l"))
    expected = [ref_logLR_asym(f32(t + pseudocount), f32(c + pseudocount))
                for (_, t, c) in SCORE_ROWS]
    assert column(s, 3) == pytest.approx(expected, rel=1e-5, abs=1e-6)


def test_ScoreTrackII_score_l_equal_values_is_zero():
    s = make_st([(10, 4.0, 4.0)])
    s.change_score_method(ord("l"))
    assert column(s, 3) == [0.0]


@pytest.mark.parametrize("pseudocount", [0.5, 1.0])
def test_ScoreTrackII_score_s(pseudocount):
    s = make_st(SCORE_ROWS, pseudocount=pseudocount)
    s.change_score_method(ord("s"))
    expected = [ref_logLR_sym(f32(t + pseudocount), f32(c + pseudocount))
                for (_, t, c) in SCORE_ROWS]
    assert column(s, 3) == pytest.approx(expected, rel=1e-5, abs=1e-6)


def test_ScoreTrackII_score_s_is_antisymmetric():
    # swapping treatment and control flips the sign
    a = make_st([(10, 9.0, 2.0), (20, 1.0, 6.0)])
    b = make_st([(10, 2.0, 9.0), (20, 6.0, 1.0)])
    a.change_score_method(ord("s"))
    b.change_score_method(ord("s"))
    assert column(a, 3) == pytest.approx([-x for x in column(b, 3)],
                                         rel=1e-6)


@pytest.mark.parametrize("code, func", [
    ("f", lambda t, c, pc: math.log10((t + pc) / (c + pc))),
    ("F", lambda t, c, pc: (t + pc) / (c + pc)),
    # subtraction and maximum do not use the pseudocount
    ("d", lambda t, c, pc: t - c),
    ("M", lambda t, c, pc: max(t, c)),
])
@pytest.mark.parametrize("pseudocount", [0.5, 1.0])
def test_ScoreTrackII_score_simple(code, func, pseudocount):
    s = make_st(SCORE_ROWS, pseudocount=pseudocount)
    s.change_score_method(ord(code))
    expected = [func(t, c, pseudocount) for (_, t, c) in SCORE_ROWS]
    assert column(s, 3) == pytest.approx(expected, rel=1e-6, abs=1e-6)


@pytest.mark.parametrize("norm", ["N", "T", "C", "M"])
def test_ScoreTrackII_score_m_is_treatment_per_million(norm):
    # whatever the normalization, SPMR is raw treatment / treatment depth
    s = make_st([(10, 8.0, 8.0), (20, 16.0, 4.0)], treat_depth=4.0,
                ctrl_depth=2.0)
    s.change_normalization_method(ord(norm))
    s.change_score_method(ord("m"))
    assert column(s, 3) == [2.0, 4.0]


def test_ScoreTrackII_score_many_chromosomes():
    s = make_st_multi({b"chr2": [(10, 5.0, 1.0)],
                       b"chr1": [(10, 2.0, 1.0), (30, 1.0, 4.0)]})
    s.change_score_method(ord("d"))
    assert column(s, 3, b"chr1") == [1.0, -3.0]
    assert column(s, 3, b"chr2") == [4.0]


@pytest.mark.parametrize("code", [ord("x"), ord("P"), ord("N"), 0])
def test_ScoreTrackII_change_score_method_unknown_code(code):
    s = make_st(SCORE_ROWS)
    with pytest.raises(NotImplementedError):
        s.change_score_method(code)


def test_ScoreTrackII_change_score_method_str_code_raises():
    with pytest.raises(TypeError, match="an integer is required"):
        make_st(SCORE_ROWS).change_score_method("p")


@pytest.mark.parametrize("code", ["p", "q"])
def test_ScoreTrackII_score_p_zero_lambda_raises(code):
    # pseudocount 0 and control 0 give a Poisson lambda of 0
    s = make_st([(10, 3.0, 0.0)], pseudocount=0.0)
    with pytest.raises(AssertionError, match="Lambda must > 0, however we got 0"):
        s.change_score_method(ord(code))


def test_ScoreTrackII_score_f_zero_control_raises():
    s = make_st([(10, 3.0, 0.0)], pseudocount=0.0)
    with pytest.raises(ZeroDivisionError, match="float division"):
        s.change_score_method(ord("f"))


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_ScoreTrackII_score_F_zero_control():
    s = make_st([(10, 0.0, 0.0), (20, 3.0, 0.0)], pseudocount=0.0)
    s.change_score_method(ord("F"))
    v = column(s, 3)
    assert math.isnan(v[0]) and v[1] == math.inf


@pytest.mark.parametrize("code", ["l", "s"])
def test_ScoreTrackII_score_llr_zero_control(code):
    # x == y gives 0; x > y = 0 gives an infinite ratio
    s = make_st([(10, 0.0, 0.0), (20, 3.0, 0.0)], pseudocount=0.0)
    s.change_score_method(ord(code))
    assert column(s, 3) == [0.0, math.inf]


# ------------------------------------
# ScoreTrackII.make_pq_table
# ------------------------------------

def test_ScoreTrackII_make_pq_table():
    s = make_st(SCORE_ROWS)
    s.change_score_method(ord("p"))
    table = dict(s.make_pq_table().items())
    stat = {}
    for (v, ln) in zip(column(s, 3),
                       region_lengths([r[0] for r in SCORE_ROWS])):
        stat[v] = stat.get(v, 0) + ln
    ref = bh_qscores(stat)
    assert sorted(table) == sorted(ref)
    for v in ref:
        assert table[v] == pytest.approx(ref[v], abs=1e-5)


def test_ScoreTrackII_make_pq_table_existing_example():
    # p-scores 60.48912 (10 bp), 0.37599, 0.07723, 0.00006 (50 bp each)
    # and 6.40804 (50 bp), N = 210: q(60.49) = 60.49 - log10(210),
    # q(6.41) = 6.41 + log10(11/210), the rest reach 0
    s = make_st([(10, 100, 10), (60, 10, 10), (110, 15, 20), (160, 5, 20),
                 (210, 20, 5)], treat_depth=10, ctrl_depth=5)
    s.change_score_method(ord("p"))
    table = dict(s.make_pq_table().items())
    p = column(s, 3)
    assert table[p[0]] == pytest.approx(p[0] - math.log10(210), abs=1e-5)
    assert table[p[4]] == pytest.approx(p[4] + math.log10(11 / 210), abs=1e-5)
    assert [table[p[i]] for i in (1, 2, 3)] == [0.0, 0.0, 0.0]


@pytest.mark.parametrize("prepare", ["new", "q"])
def test_ScoreTrackII_make_pq_table_requires_p(prepare):
    s = make_st(SCORE_ROWS)
    if prepare == "q":
        s.change_score_method(ord("q"))
    with pytest.raises(AssertionError):
        s.make_pq_table()


# ------------------------------------
# ScoreTrackII.write_bedGraph / enable_trackline
# ------------------------------------

# treat [1, 1, 1.000002, 3, 3]: 1.000002 is within 1e-5 of 1 and merged
WB_ROWS = [(10, 1.0, 2.0), (20, 1.0, 2.0), (30, 1.000002, 2.0),
           (40, 3.0, 2.5), (50, 3.0, 2.5)]


@pytest.mark.parametrize("col, prepare, expected", [
    (1, None, "chrY\t0\t30\t1.00000\nchrY\t30\t50\t3.00000\n"),
    (2, None, "chrY\t0\t30\t2.00000\nchrY\t30\t50\t2.50000\n"),
    (3, "d", "chrY\t0\t30\t-1.00000\nchrY\t30\t50\t0.50000\n"),
])
def test_ScoreTrackII_write_bedGraph(col, prepare, expected):
    s = make_st(WB_ROWS)
    if prepare:
        s.change_score_method(ord(prepare))
    fhd = io.StringIO()
    assert s.write_bedGraph(fhd, "NAME", "DESC", col) is True
    assert fhd.getvalue() == expected


def test_ScoreTrackII_write_bedGraph_default_column_is_score():
    s = make_st(WB_ROWS)
    s.change_score_method(ord("M"))
    fhd = io.StringIO()
    s.write_bedGraph(fhd, "NAME", "DESC")
    assert fhd.getvalue() == "chrY\t0\t30\t2.00000\nchrY\t30\t50\t3.00000\n"


def test_ScoreTrackII_write_bedGraph_splits_beyond_5_digits():
    s = make_st([(10, 1.0, 1.0), (20, 1.0001, 1.0), (25, 0.1, 1.0)])
    fhd = io.StringIO()
    s.write_bedGraph(fhd, "NAME", "DESC", 1)
    assert fhd.getvalue() == ("chrY\t0\t10\t1.00000\nchrY\t10\t20\t1.00010\n"
                              "chrY\t20\t25\t0.10000\n")


def test_ScoreTrackII_write_bedGraph_many_chromosomes_sorted():
    s = make_st_multi({b"chr2": [(10, 5.0, 1.0)],
                       b"chr1": [(10, 2.0, 1.0), (30, 1.0, 4.0)],
                       b"chr3": []})
    fhd = io.StringIO()
    s.write_bedGraph(fhd, "NAME", "DESC", 1)
    # chr3 has no rows and is skipped
    assert fhd.getvalue() == ("chr1\t0\t10\t2.00000\nchr1\t10\t30\t1.00000\n"
                              "chr2\t0\t10\t5.00000\n")


def test_ScoreTrackII_write_bedGraph_empty_track():
    fhd = io.StringIO()
    assert ScoreTrackII(1.0, 1.0).write_bedGraph(fhd, "N", "D", 3) is True
    assert fhd.getvalue() == ""


@pytest.mark.parametrize("col", [0, 4, -1])
def test_ScoreTrackII_write_bedGraph_bad_column(col):
    s = make_st(WB_ROWS)
    with pytest.raises(AssertionError,
                       match="column should be between 1, 2 or 3."):
        s.write_bedGraph(io.StringIO(), "NAME", "DESC", col)


def test_ScoreTrackII_enable_trackline_returns_none():
    s = make_st(WB_ROWS)
    assert s.enable_trackline() is None


# ------------------------------------
# ScoreTrackII.cutoff_analysis
# ------------------------------------
#
# Same sweep as bedGraphTrackI.cutoff_analysis, on the score column.
# Rows give 'd' scores [0, 4, 0, 2, 0] over
# [0,100) [100,200) [200,300) [300,350) [350,1000).

CA_HEADER = "score\tnpeaks\tlpeaks\tavelpeak\n"
CA_ROWS = [(100, 1.0, 1.0), (200, 5.0, 1.0), (300, 1.0, 1.0),
           (350, 3.0, 1.0), (1000, 1.0, 1.0)]


def make_ca_track(rows=CA_ROWS):
    s = make_st(rows)
    s.change_score_method(ord("d"))
    return s


@pytest.mark.parametrize("kwargs, body", [
    (dict(max_gap=0, min_length=0, steps=4),
     "3.00\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.00\t2\t150\t75.00\n0.00\t2\t150\t75.00\n"),
    (dict(max_gap=100, min_length=0, steps=4),
     "3.00\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.00\t1\t250\t250.00\n0.00\t1\t250\t250.00\n"),
    (dict(max_gap=0, min_length=100, steps=4),
     "3.00\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.00\t1\t100\t100.00\n0.00\t1\t100\t100.00\n"),
    (dict(max_gap=0, min_length=200, steps=4), ""),
    (dict(max_gap=0, min_length=0, steps=4, min_score=1, max_score=3),
     "2.50\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.50\t2\t150\t75.00\n1.00\t2\t150\t75.00\n"),
    # header only: no steps, or an empty score range
    (dict(max_gap=0, min_length=0, steps=0), ""),
    (dict(max_gap=0, min_length=0, steps=-1), ""),
    (dict(max_gap=0, min_length=0, steps=4, min_score=10), ""),
])
def test_ScoreTrackII_cutoff_analysis(kwargs, body):
    assert make_ca_track().cutoff_analysis(**kwargs) == CA_HEADER + body


def test_ScoreTrackII_cutoff_analysis_defaults():
    # defaults max_gap 50, min_length 200, steps 100 over scores 0..4:
    # no peak is 200 bp long
    assert make_ca_track().cutoff_analysis() == CA_HEADER


def test_ScoreTrackII_cutoff_analysis_many_chromosomes():
    s = make_st_multi({b"chrY": CA_ROWS,
                       b"chr2": [(50, 1.0, 1.0), (150, 4.0, 1.0),
                                 (200, 1.0, 1.0)]})
    s.change_score_method(ord("d"))
    assert s.cutoff_analysis(max_gap=0, min_length=0, steps=4) == \
        CA_HEADER + ("3.00\t1\t100\t100.00\n2.00\t2\t200\t100.00\n"
                     "1.00\t3\t250\t83.33\n0.00\t3\t250\t83.33\n")


@pytest.mark.parametrize("rows", [[], [(100, 2.0, 1.0), (200, 2.0, 1.0)]])
def test_ScoreTrackII_cutoff_analysis_empty_or_constant(rows):
    s = make_st(rows) if rows else ScoreTrackII(1.0, 1.0)
    s.change_score_method(ord("d"))
    assert s.cutoff_analysis(max_gap=0, min_length=0, steps=4) == CA_HEADER


# ------------------------------------
# ScoreTrackII.call_peaks
# ------------------------------------
#
# Rows give 'd' scores [0, 10, 3, 10, 0] over 10 bp regions from 0 to
# 50 and treatment pileup [1, 11, 4, 11, 1] (control 1). The summit is
# the middle of the region with the highest treatment pileup; ties
# take index int((n + 1) / 2) - 1. pileup, pscore (without
# pseudocount) and fold change (with pseudocount 1) are those of the
# summit region.

CP_ROWS = [(10, 1.0, 1.0), (20, 11.0, 1.0), (30, 4.0, 1.0),
           (40, 11.0, 1.0), (50, 1.0, 1.0)]


def cp_peak(start, end, summit):
    return (b"chrY", start, end, summit, 10.0, 11.0,
            pytest.approx(ref_pscore(11, 1.0), abs=1e-4), 6.0, -1.0)


@pytest.mark.parametrize("cutoff, min_length, max_gap, expected", [
    # one peak [10,40); pileup 11 tied in [10,20) and [30,40): first
    (1, 0, 0, [cp_peak(10, 40, 15)]),
    # score >= 5 only in [10,20) and [30,40), 10 bp apart
    (5, 0, 0, [cp_peak(10, 20, 15), cp_peak(30, 40, 35)]),
    (5, 0, 9, [cp_peak(10, 20, 15), cp_peak(30, 40, 35)]),
    (5, 0, 10, [cp_peak(10, 40, 15)]),
    # min_length is inclusive
    (5, 10, 0, [cp_peak(10, 20, 15), cp_peak(30, 40, 35)]),
    (5, 11, 0, []),
    (1, 30, 0, [cp_peak(10, 40, 15)]),
    # cutoff is inclusive
    (10, 0, 0, [cp_peak(10, 20, 15), cp_peak(30, 40, 35)]),
    (10.5, 0, 0, []),
])
def test_ScoreTrackII_call_peaks(cutoff, min_length, max_gap, expected):
    s = make_st(CP_ROWS)
    s.change_score_method(ord("d"))
    peaks = s.call_peaks(cutoff=cutoff, min_length=min_length,
                         max_gap=max_gap)
    assert isinstance(peaks, PeakIO)
    assert peak_fields(peaks) == expected


def test_ScoreTrackII_call_peaks_summit_follows_treatment_pileup():
    # scores [0, 9, 4, 0]; pileup 12 in [20,30) beats 10 in [10,20), so
    # the summit is there and the peak score is the score there (4).
    # Surprising: __close_peak's docstring says "the region with the
    # highest score", but the code ranks regions by treatment pileup, as
    # CallPeakUnit does (its comments call the pileup the "general score
    # to find summit"), so the docstring's "score" is ambiguous. The
    # reported pscore has no pseudocount, matching CallPeakUnit and the
    # class comment that the pseudocount is for logLR, FE and logFE.
    s = make_st([(10, 1.0, 1.0), (20, 10.0, 1.0), (30, 12.0, 8.0),
                 (40, 1.0, 1.0)])
    s.change_score_method(ord("d"))
    assert peak_fields(s.call_peaks(cutoff=1, min_length=0, max_gap=0)) == \
        [(b"chrY", 10, 30, 25, 4.0, 12.0,
          pytest.approx(ref_pscore(12, 8.0), abs=1e-4),
          pytest.approx(13 / 9, rel=1e-6), -1.0)]


def test_ScoreTrackII_call_peaks_first_region_starts_at_zero():
    s = make_st([(10, 5.0, 1.0), (20, 1.0, 1.0)])
    s.change_score_method(ord("d"))
    assert [p[1:4] for p in peak_fields(
        s.call_peaks(cutoff=1, min_length=0, max_gap=0))] == [(0, 10, 5)]


def test_ScoreTrackII_call_peaks_qscore_reported_for_q_scores():
    # q-scores: [1000,1100) has the top p-score (k = 1, N = 2000), the
    # rest has q = 0; the peak's qscore is its score
    s = make_st([(1000, 1.0, 5.0), (1100, 30.0, 5.0), (2000, 1.0, 5.0)])
    s.change_score_method(ord("q"))
    q = ref_pscore(31, 6.0) - math.log10(2000)
    rows = peak_fields(s.call_peaks(cutoff=1, min_length=0, max_gap=0))
    assert rows == [(b"chrY", 1000, 1100, 1050, pytest.approx(q, abs=1e-4),
                     30.0, pytest.approx(ref_pscore(30, 5.0), abs=1e-4),
                     pytest.approx(31 / 6, rel=1e-6),
                     pytest.approx(q, abs=1e-4))]


def test_ScoreTrackII_call_peaks_many_chromosomes_and_no_peaks():
    s = make_st_multi({b"chr2": [(10, 1.0, 1.0), (20, 5.0, 1.0)],
                       b"chr1": [(30, 3.0, 1.0)],
                       b"chr3": [(30, 1.0, 1.0)]})
    s.change_score_method(ord("d"))
    assert [p[:4] for p in peak_fields(
        s.call_peaks(cutoff=1, min_length=0, max_gap=0))] == \
        [(b"chr1", 0, 30, 15), (b"chr2", 10, 20, 15)]
    assert s.call_peaks(cutoff=100, min_length=0, max_gap=0).total == 0
    assert ScoreTrackII(1.0, 1.0).call_peaks().total == 0


def hump_track(centers):
    """Treatment pileup made of triangles of height 40 and half-width
    100 bp sampled in 5 bp regions from 1000 to 1400; control 1. With
    'd' scores and cutoff 2, the region above cutoff is [1005, 1395) for
    two humps at 1100 and 1300 (the dip between them is 10 bp), and the
    highest pileup (39) of a hump at c is in [c-5, c+5)."""
    rows = [(1000, 0.0, 1.0)]
    for x in range(1000, 1400, 5):
        v = max(max(0.0, 40 - abs(x + 2.5 - c) * 0.4) for c in centers)
        rows.append((x + 5, float(round(v)), 1.0))
    rows.append((3000, 0.0, 1.0))
    s = make_st(rows)
    s.change_score_method(ord("d"))
    return s


def test_ScoreTrackII_call_peaks_two_humps_without_summits():
    # four regions share the top pileup 39: index 1 -> [1100,1105)
    s = hump_track([1100, 1300])
    rows = peak_fields(s.call_peaks(cutoff=2, min_length=10, max_gap=100))
    assert [r[:6] for r in rows] == [(b"chrY", 1005, 1395, 1102, 38.0, 39.0)]


@pytest.mark.parametrize("centers", [[1200], [1100, 1300]])
def test_ScoreTrackII_call_peaks_call_summits(centers):
    # one summit per hump, at the hump's top
    s = hump_track(centers)
    peaks = s.call_peaks(cutoff=2, min_length=10, max_gap=100,
                         call_summits=True)
    rows = peak_fields(peaks)
    assert len(rows) == len(centers)
    for row, c in zip(rows, centers):
        assert c - 5 <= row[3] < c + 5
        assert (row[4], row[5]) == (38.0, 39.0)


def test_ScoreTrackII_call_peaks_call_summits_boundaries():
    """Regression test: call_peaks(call_summits=True) widened peak
    boundaries by the 10 bp summit-search padding (call_summits=False and
    CallPeakUnit report the region above the cutoff).

    Fixed upstream in 1622eb8 (#749, issue #747), which reports
    peak_start/peak_end instead of the padded start/end.
    """
    s = hump_track([1100, 1300])
    peaks = s.call_peaks(cutoff=2, min_length=10, max_gap=100,
                         call_summits=True)
    assert [r[1:3] for r in peak_fields(peaks)] == [(1005, 1395),
                                                    (1005, 1395)]


# ------------------------------------
# ScoreTrackII.call_broadpeaks
# ------------------------------------

def test_ScoreTrackII_call_broadpeaks_without_strong_peaks():
    # scores [0, 2, 0]: no lvl1 peak, so nothing is reported
    s = make_st([(10, 1.0, 1.0), (20, 3.0, 1.0), (30, 1.0, 1.0)])
    s.change_score_method(ord("d"))
    bp = s.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=0,
                           lvl1_max_gap=0, lvl2_max_gap=5)
    assert isinstance(bp, BroadPeakIO)
    assert bp.peaks == {}


@pytest.mark.parametrize("kwargs, message", [
    (dict(lvl1_cutoff=1, lvl2_cutoff=1),
     "level 1 cutoff should be larger than level 2."),
    (dict(lvl1_cutoff=5, lvl2_cutoff=1, lvl1_max_gap=400, lvl2_max_gap=400),
     "level 2 maximum gap should be larger than level 1."),
])
def test_ScoreTrackII_call_broadpeaks_bad_arguments(kwargs, message):
    with pytest.raises(AssertionError, match=message):
        make_st(CP_ROWS).call_broadpeaks(**kwargs)


# ------------------------------------
# TwoConditionScores (build covers build_chromosome, get_common_chrs,
# add_chromosome, add; call_peaks covers the peak merging and
# mean_from_peakcontent)
# ------------------------------------

def make_bdg(rows, chrom=b"chr1"):
    t = bedGraphTrackI()
    for (s, e, v) in rows:
        t.add_loc(chrom, s, e, v)
    return t


# t1: [0,10)=10 [10,30)=2      c1: [0,20)=1 [20,30)=3
# t2: [0,30)=1                 c2: [0,5)=1  [5,30)=2
# union of breakpoints: 5, 10, 20, 30, giving the intervals
#   [0,5)=(10,1,1,1) [5,10)=(10,1,1,2) [10,20)=(2,1,1,2) [20,30)=(2,3,1,2)
TC_TRACKS = ([(0, 10, 10.0), (10, 30, 2.0)],
             [(0, 20, 1.0), (20, 30, 3.0)],
             [(0, 30, 1.0)],
             [(0, 5, 1.0), (5, 30, 2.0)])
TC_INTERVALS = [(10.0, 1.0, 1.0, 1.0), (10.0, 1.0, 1.0, 2.0),
                (2.0, 1.0, 1.0, 2.0), (2.0, 3.0, 1.0, 2.0)]


def make_tc(tracks=TC_TRACKS, f1=1.0, f2=1.0, pseudocount=None,
            chroms=(b"chr1",) * 4, build=True, finalize=True):
    bdgs = [make_bdg(rows, chrom) for rows, chrom in zip(tracks, chroms)]
    if pseudocount is None:
        tc = TwoConditionScores(bdgs[0], bdgs[1], bdgs[2], bdgs[3], f1, f2)
    else:
        tc = TwoConditionScores(bdgs[0], bdgs[1], bdgs[2], bdgs[3], f1, f2,
                                pseudocount)
    if build:
        tc.build()
    if finalize:
        tc.finalize()
    return tc


def tc_ref_scores(intervals, f1, f2, pc):
    """The three logLR columns for each interval, with the float32
    arithmetic of TwoConditionScores.add."""
    pc = f32(pc)
    out = []
    for (t1, c1, t2, c2) in intervals:
        x1 = f32(f32(t1 + pc) * f32(f1))
        y1 = f32(f32(c1 + pc) * f32(f1))
        x2 = f32(f32(t2 + pc) * f32(f2))
        y2 = f32(f32(c2 + pc) * f32(f2))
        out.append((ref_logLR_asym(x1, y1), ref_logLR_asym(x2, y2),
                    ref_logLR_sym(x1, x2)))
    return [list(col) for col in zip(*out)]


@pytest.mark.parametrize("f1, f2, pseudocount", [
    (1.0, 1.0, None),            # default pseudocount 0.01
    (0.5, 2.0, 1.0),
    (1.0, 0.25, 0.5),
])
def test_TwoConditionScores_build_scores(f1, f2, pseudocount):
    tc = make_tc(f1=f1, f2=f2, pseudocount=pseudocount)
    ref = tc_ref_scores(TC_INTERVALS, f1, f2,
                        0.01 if pseudocount is None else pseudocount)
    data = tc.get_data_by_chr(b"chr1")
    for i in (1, 2, 3):
        assert data[i].tolist() == pytest.approx(ref[i - 1], rel=1e-5,
                                                 abs=1e-6)


def test_TwoConditionScores_set_pseudocount():
    tc = make_tc(build=False, finalize=False)
    assert tc.set_pseudocount(1.0) is None
    tc.build()
    tc.finalize()
    ref = tc_ref_scores(TC_INTERVALS, 1.0, 1.0, 1.0)
    assert tc.get_data_by_chr(b"chr1")[3].tolist() == \
        pytest.approx(ref[2], rel=1e-5, abs=1e-6)


def test_TwoConditionScores_finalize_trims_capacity():
    # capacity is the sum of the region counts of the four tracks
    # (2 + 2 + 1 + 2 = 7); four intervals are stored
    tc = make_tc(finalize=False)
    assert [len(a) for a in tc.get_data_by_chr(b"chr1")] == [7, 7, 7, 7]
    assert tc.finalize() is None
    assert [len(a) for a in tc.get_data_by_chr(b"chr1")] == [4, 4, 4, 4]
    assert [a.dtype for a in tc.get_data_by_chr(b"chr1")] == \
        [np.int32, np.float32, np.float32, np.float32]


def test_TwoConditionScores_only_common_chromosomes():
    t1 = make_bdg([(0, 10, 5.0)], b"chr1")
    t1.add_loc(b"chr2", 0, 10, 5.0)
    t1.add_loc(b"chr3", 0, 10, 5.0)
    others = []
    for _ in range(3):
        b = make_bdg([(0, 10, 1.0)], b"chr1")
        b.add_loc(b"chr2", 0, 10, 1.0)
        others.append(b)
    # chr2 is missing from c2
    others[2] = make_bdg([(0, 10, 1.0)], b"chr1")
    tc = TwoConditionScores(t1, others[0], others[1], others[2])
    tc.build()
    tc.finalize()
    assert tc.get_chr_names() == {b"chr1"}
    assert tc.get_data_by_chr(b"chr2") is None


def test_TwoConditionScores_no_common_chromosome():
    tc = make_tc(chroms=(b"chr1", b"chr1", b"chr1", b"chr2"))
    assert tc.get_chr_names() == set()
    fhd = io.StringIO()
    assert tc.write_bedGraph(fhd, "N", "D", 3) is True
    assert tc.write_matrix(fhd, "N", "D") is True
    assert fhd.getvalue() == ""
    cats = tc.call_peaks()
    assert len(cats) == 3
    assert [c.total for c in cats] == [0, 0, 0]


# write_bedGraph / write_matrix / call_peaks are checked against the
# arrays the object stores (get_data_by_chr), so that they test the
# writers and the peak caller on their own, independent of the build
# bug above.

def ref_bdg_text(chrom, pos, val):
    """bedGraph lines: consecutive values closer than 1e-6 are merged."""
    lines = []
    pre = 0
    pre_v = val[0]
    for i in range(1, len(pos)):
        if abs(pre_v - val[i]) >= 1e-6:
            lines.append("%s\t%d\t%d\t%.5f\n" % (chrom, pre, pos[i - 1], pre_v))
            pre_v = val[i]
            pre = pos[i - 1]
    lines.append("%s\t%d\t%d\t%.5f\n" % (chrom, pre, pos[-1], pre_v))
    return "".join(lines)


@pytest.mark.parametrize("col", [1, 2, 3])
def test_TwoConditionScores_write_bedGraph(col):
    tc = make_tc()
    data = tc.get_data_by_chr(b"chr1")
    fhd = io.StringIO()
    assert tc.write_bedGraph(fhd, "NAME", "DESC", col) is True
    assert fhd.getvalue() == ref_bdg_text("chr1", data[0].tolist(),
                                          data[col].tolist())


def test_TwoConditionScores_write_bedGraph_merges_equal_values():
    # t1 constant: column 3 (t1 vs t2) is equal on [0,5) and [5,10)
    tc = make_tc(tracks=([(0, 30, 4.0)], [(0, 30, 1.0)], [(0, 30, 1.0)],
                         [(0, 5, 1.0), (5, 30, 2.0)]))
    data = tc.get_data_by_chr(b"chr1")
    assert len(set(data[3].tolist())) == 1
    fhd = io.StringIO()
    tc.write_bedGraph(fhd, "NAME", "DESC", 3)
    assert fhd.getvalue().count("\n") == 1


@pytest.mark.parametrize("col", [0, 4])
def test_TwoConditionScores_write_bedGraph_bad_column(col):
    with pytest.raises(AssertionError,
                       match="column should be between 1, 2 or 3."):
        make_tc().write_bedGraph(io.StringIO(), "NAME", "DESC", col)


def test_TwoConditionScores_write_matrix():
    tc = make_tc()
    (pos, v1, v2, v3) = [a.tolist() for a in tc.get_data_by_chr(b"chr1")]
    fhd = io.StringIO()
    assert tc.write_matrix(fhd, "NAME", "DESC") is True
    starts = [0] + pos[:-1]
    assert fhd.getvalue() == "".join(
        "chr1:%d_%d\t%.5f\t%.5f\t%.5f\n" % row
        for row in zip(starts, pos, v1, v2, v3))
    assert fhd.getvalue().count("\n") == 4


# Differential design over 100 bp regions (all intervals of the four
# tracks fall on multiples of 100):
#   [100,200) t1=100, [200,300) t1=50  -> condition 1 stronger (cat1)
#   [400,500) t2=100                   -> condition 2 stronger (cat2)
#   [600,700) t1=t2=100                -> both enriched (cat3)
TC_DIFF = ([(0, 100, 1.0), (100, 200, 100.0), (200, 300, 50.0),
            (300, 600, 1.0), (600, 700, 100.0), (700, 800, 1.0)],
           [(0, 800, 1.0)],
           [(0, 400, 1.0), (400, 500, 100.0), (500, 600, 1.0),
            (600, 700, 100.0), (700, 800, 1.0)],
           [(0, 800, 1.0)])


def ref_tc_peaks(pos, c1, c2, c3, cutoff, min_length, max_gap):
    """(start, end, mean score) per category from the stored arrays, as
    documented in call_peaks: entries are merged when the gap is at most
    max_gap, kept when at least min_length long, and scored by the
    length-weighted mean of t1-vs-t2 (negated for category 2, absolute
    for category 3)."""
    cats = ([], [], [])
    tests = (lambda i: c1[i] >= cutoff and c3[i] >= cutoff,
             lambda i: c2[i] >= cutoff and c3[i] <= -cutoff,
             lambda i: (c1[i] >= cutoff and c2[i] >= cutoff and
                        -cutoff <= c3[i] <= cutoff))
    scores = (lambda i: c3[i], lambda i: -c3[i], lambda i: abs(c3[i]))
    for k in range(3):
        content = []
        for i in range(len(pos)):
            if not tests[k](i):
                continue
            item = (pos[i - 1] if i > 0 else 0, pos[i], scores[k](i))
            if content and item[0] - content[-1][1] > max_gap:
                cats[k].append(content)
                content = []
            content.append(item)
        if content:
            cats[k].append(content)
    out = []
    for k in range(3):
        peaks = []
        for content in cats[k]:
            length = content[-1][1] - content[0][0]
            if length >= min_length:
                mean = sum(v * (e - s) for (s, e, v) in content) / length
                peaks.append((content[0][0], content[-1][1], mean))
        out.append(peaks)
    return out


@pytest.mark.parametrize("cutoff, min_length, max_gap", [
    (3, 1, 0),
    (3, 150, 0),
    (3, 1, 250),
    (100, 1, 0),
])
def test_TwoConditionScores_call_peaks(cutoff, min_length, max_gap):
    tc = make_tc(tracks=TC_DIFF)
    data = [a.tolist() for a in tc.get_data_by_chr(b"chr1")]
    ref = ref_tc_peaks(*data, cutoff, min_length, max_gap)
    cats = tc.call_peaks(cutoff=cutoff, min_length=min_length,
                         max_gap=max_gap)
    assert isinstance(cats, tuple) and len(cats) == 3
    for peaks, expected in zip(cats, ref):
        got = [(p["start"], p["end"], p["summit"], p["pileup"], p["pscore"],
                p["fc"], p["qscore"])
               for p in peaks.get_data_from_chrom(b"chr1")]
        assert got == [(s, e, -1, 0.0, 0.0, 0.0, 0.0)
                       for (s, e, _) in expected]


def test_TwoConditionScores_call_peaks_one_peak_per_category():
    # independent of positions: each category is one block of entries
    tc = make_tc(tracks=TC_DIFF)
    cats = tc.call_peaks(cutoff=3, min_length=1, max_gap=0)
    assert [c.total for c in cats] == [1, 1, 1]


def test_TwoConditionScores_call_peaks_score_is_mean_logLR():
    """Regression test: mean_from_peakcontent stored each score in a C
    long, truncating it to an integer before averaging.

    Fixed upstream in 9dbb3c7 (#739, issue #715), which declares the score
    and accumulator as cython.double.
    """
    tc = make_tc(tracks=TC_DIFF)
    data = [a.tolist() for a in tc.get_data_by_chr(b"chr1")]
    ref = ref_tc_peaks(*data, 3, 1, 0)
    cats = tc.call_peaks(cutoff=3, min_length=1, max_gap=0)
    for peaks, expected in zip(cats, ref):
        got = [p["score"] for p in peaks.get_data_from_chrom(b"chr1")]
        assert got == pytest.approx([m for (_, _, m) in expected], rel=1e-5)
