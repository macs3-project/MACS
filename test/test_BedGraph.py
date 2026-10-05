# Time-stamp: <2025-02-10 08:46:48 Tao Liu>

import pytest

from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.IO.PeakIO import PeakIO

import math
from array import array

import numpy as np
from scipy.stats import chi2

from MACS3.Signal.BedGraph import (bedGraphTrackII,)
from MACS3.Signal.ScoreTrack import (ScoreTrackII,)
from MACS3.IO.PeakIO import (BroadPeakIO,)

test_regions1 = [(b"chrY", 0, 10593155, 0.0),
                 (b"chrY", 10593155, 10597655, 0.0066254580149)]
overlie_max = [(b"chrY", 0, 75, 20.0),
               (b"chrY", 75, 85, 35.0),
               (b"chrY", 85, 90, 75.0),
               (b"chrY", 90, 150, 10.0),
               (b"chrY", 150, 155, 0.0)]
overlie_mean = [(b"chrY", 0, 70, 6.66667),
                (b"chrY", 70, 75, 9.0),
                (b"chrY", 75, 80, 14.0),
                (b"chrY", 80, 85, 11.66667),
                (b"chrY", 85, 90, 36.66667),
                (b"chrY", 90, 150, 3.33333),
                (b"chrY", 150, 155, 0.0)]
overlie_fisher = [(b"chrY", 0, 70, (0, 0, 20), 92.10340371976183, 1.1074313239555578e-17, 16.9557),
                  (b"chrY", 70, 75, (7, 0, 20), 124.33959502167846, 1.9957116587802055e-24, 23.6999),
                  (b"chrY", 75, 80, (7, 0, 35), 193.41714781149986, 4.773982707347631e-39, 38.3211),
                  (b"chrY", 80, 85, (0, 0, 35), 161.1809565095832, 3.329003070922764e-32, 31.4777),
                  (b"chrY", 85, 90, (0, 75, 35), 506.56872045869005, 3.233076792862357e-106, 105.4904),
                  (b"chrY", 90, 150, (0, 0, 10), 46.051701859880914, 2.8912075645386016e-08, 7.5389),
                  (b"chrY", 150, 155, (0, 0, 0), 0.0, 1.0, 0.0)]
some_peak = [(b"chrY", 50, 80),
             (b"chrY", 95, 152)]


@pytest.fixture
def define_regions():
    test_regions1 = [(b"chrY", 0, 70, 0.0),
                     (b"chrY", 70, 80, 7.0),
                     (b"chrY", 80, 150, 0.0)]
    test_regions2 = [(b"chrY", 0, 85, 0.0),
                     (b"chrY", 85, 90, 75.0),
                     (b"chrY", 90, 155, 0.0)]
    test_regions3 = [(b"chrY", 0, 75, 20.0),
                     (b"chrY", 75, 90, 35.0),
                     (b"chrY", 90, 150, 10.0)]
    bdg1 = bedGraphTrackI()
    bdg2 = bedGraphTrackI()
    bdg3 = bedGraphTrackI()
    for a in test_regions1:
        bdg1.add_loc(a[0], a[1], a[2], a[3])
    for a in test_regions2:
        bdg2.add_loc(a[0], a[1], a[2], a[3])
    for a in test_regions3:
        bdg3.add_loc(a[0], a[1], a[2], a[3])
    return (bdg1, bdg2, bdg3)


def test_add_loc1():
    bdg = bedGraphTrackI()
    for a in test_regions1:
        bdg.add_loc(a[0], a[1], a[2], a[3])


def test_refine_peak():
    bdg = bedGraphTrackI()
    for a in overlie_mean:
        bdg.add_loc(a[0], a[1], a[2], a[3])
    peak = PeakIO()
    for a in some_peak:
        peak.add(a[0], a[1], a[2])
    new_peak = bdg.refine_peaks(peak)
    out = str(new_peak)
    std = "chrom:chrY\tstart:50\tend:80\tname:peak_1\tscore:14\tsummit:77\nchrom:chrY\tstart:95\tend:152\tname:peak_2\tscore:3.33333\tsummit:122\n"
    assert out == std


def test_refine_peak_nested_input_peaks():
    # regression test: IndexError raised from __close_peak() when
    # refine_peaks() is given two input peaks where the second peak's
    # end has already been passed by the bedGraph cursor by the time
    # the cursor would otherwise overlap it (e.g. two closely-spaced
    # HMMRATAC "open" peaks against a coarse fold-change bedGraph track).
    bdg_regions = [(b"chrY", 0, 15, 1.0),
                   (b"chrY", 15, 60, 2.0),
                   (b"chrY", 60, 100, 3.0)]
    nested_peaks = [(b"chrY", 10, 50),
                    (b"chrY", 20, 30)]
    bdg = bedGraphTrackI()
    for a in bdg_regions:
        bdg.add_loc(a[0], a[1], a[2], a[3])
    peak = PeakIO()
    for a in nested_peaks:
        peak.add(a[0], a[1], a[2])

    # this used to raise IndexError: list index out of range
    new_peak = bdg.refine_peaks(peak)

    result = new_peak.get_data_from_chrom(b"chrY")
    assert len(result) == 1
    assert result[0]["start"] == 10
    assert result[0]["end"] == 50


# ------------------------------------
# Helpers for the tests below
# ------------------------------------
#
# A bedGraphTrackI stores, per chromosome, the end positions of
# consecutive regions and the value of each region; region i spans
# [p[i-1], p[i]) with p[-1] taken as 0. Expected values below are
# derived by hand from that representation.

def make_track(rows, baseline=0):
    """bedGraphTrackI built with add_loc from (chrom, start, end, value)."""
    t = bedGraphTrackI(baseline_value=baseline)
    for (chrom, s, e, v) in rows:
        t.add_loc(chrom, s, e, v)
    return t


def track_content(track):
    """{chrom: (positions, values)} as plain lists, for exact comparison."""
    out = {}
    for chrom in sorted(track.get_chr_names()):
        (p, v) = track.get_data_by_chr(chrom)
        out[chrom] = (list(p), list(v))
    return out


def f32(x):
    """Round a Python float to float32 precision (the storage type)."""
    return float(np.float32(x))


def peak_rows(peakio):
    """(chrom, start, end, summit, score) of every peak, chromosomes sorted."""
    rows = []
    for chrom in sorted(peakio.get_chr_names()):
        for p in peakio.get_data_from_chrom(chrom):
            rows.append((chrom, p["start"], p["end"], p["summit"], p["score"]))
    return rows


def broad_rows(bpeaks):
    """Every field that the gappedPeak writer uses, chromosomes sorted."""
    rows = []
    for chrom in sorted(bpeaks.peaks.keys()):
        for p in bpeaks.peaks[chrom]:
            rows.append((chrom, p["start"], p["end"], p["thickStart"],
                         p["thickEnd"], p["blockNum"], p["blockSizes"],
                         p["blockStarts"], p["score"]))
    return rows


def make_track2(chrom_rows):
    """bedGraphTrackII from {chrom: [(end, value), ...]} via add_chrom_data,
    then finalize (add_chrom_data rather than add_loc)."""
    t = bedGraphTrackII()
    for chrom, rows in chrom_rows.items():
        t.add_chrom_data(chrom, np.array(rows, dtype=[("p", "u4"),
                                                      ("v", "f4")]))
    t.finalize()
    return t


def fisher_ref(values):
    """-log10 of Fisher's combined p-value for -log10 p-values ``values``:
    statistic 2 * sum(-ln p) = 2 * ln(10) * sum(values), chi2 with 2k df."""
    stat = 2.0 * math.log(10) * sum(values)
    if stat <= 0:
        return 0.0
    return -math.log10(chi2.sf(stat, 2 * len(values)))


# Peak with ties used by both track classes: a 0 region [0,10), then
# 2n-1 regions of 10 bp alternating 9 and 2 starting at 10, then a 0
# region. The n regions with value 9 start at 10, 30, 50, 70, so their
# midpoints are 15, 35, 55, 75; the summit is the midpoint with index
# int((n + 1) / 2) - 1 among them.
TIE_SUMMIT = {1: 15, 2: 15, 3: 35, 4: 35}


def tie_values(n):
    return [9.0 if i % 2 == 0 else 2.0 for i in range(2 * n - 1)]


# Track used by several call_peaks/cutoff tests (one chromosome):
# [0,10)=0, [10,20)=5, [20,30)=0, [30,40)=7, [40,100)=0
PEAK_ROWS = [(b"chr1", 0, 10, 0.0), (b"chr1", 10, 20, 5.0),
             (b"chr1", 20, 30, 0.0), (b"chr1", 30, 40, 7.0),
             (b"chr1", 40, 100, 0.0)]


# ------------------------------------
# bedGraphTrackI.__init__ and public attributes
# ------------------------------------

def test_bedGraphTrackI_init_defaults():
    t = bedGraphTrackI()
    assert t.get_chr_names() == set()
    assert t.total() == 0
    # sentinels chosen so that the first add_loc updates them
    assert t.maxvalue == -10000000
    assert t.minvalue == 10000000
    assert t.baseline_value == 0


@pytest.mark.parametrize("baseline", [0.0, 2.5, -1.0, 0.1])
def test_bedGraphTrackI_init_baseline_is_float32(baseline):
    t = bedGraphTrackI(baseline_value=baseline)
    assert t.baseline_value == f32(baseline)


# ------------------------------------
# bedGraphTrackI.add_loc
# ------------------------------------

@pytest.mark.parametrize("rows, baseline, expected", [
    # identical neighbours are merged into one region
    ([(0, 10, 1.0), (10, 20, 1.0)], 0, ([20], [1.0])),
    ([(0, 10, 2.0), (10, 20, 2.0), (20, 30, 2.0)], 0, ([30], [2.0])),
    # different neighbours are kept apart
    ([(0, 10, 1.0), (10, 20, 2.0)], 0, ([10, 20], [1.0, 2.0])),
    ([(0, 10, 1.0), (10, 20, 2.0), (20, 30, 1.0)], 0,
     ([10, 20, 30], [1.0, 2.0, 1.0])),
    # a first region starting after 0 gets a baseline block [0, start)
    ([(5, 10, 1.0)], 0, ([5, 10], [0.0, 1.0])),
    ([(5, 10, 1.0)], 3, ([5, 10], [3.0, 1.0])),
    # a negative start is clamped to 0
    ([(-5, 10, 1.0)], 0, ([10], [1.0])),
    ([(-2 ** 31, 10, 1.0)], 0, ([10], [1.0])),
    # a single interval
    ([(0, 1, 4.0)], 0, ([1], [4.0])),
])
def test_bedGraphTrackI_add_loc(rows, baseline, expected):
    t = make_track([(b"chr1",) + r for r in rows], baseline=baseline)
    assert track_content(t) == {b"chr1": expected}


@pytest.mark.parametrize("end", [0, -1, -100])
def test_bedGraphTrackI_add_loc_ignores_nonpositive_end(end):
    t = bedGraphTrackI()
    t.add_loc(b"chr1", -10, end, 1.0)
    assert t.get_chr_names() == set()
    # max/min are untouched because the call returns early
    assert t.maxvalue == -10000000 and t.minvalue == 10000000


def test_bedGraphTrackI_add_loc_updates_max_min():
    t = make_track([(b"chr1", 0, 10, 3.0), (b"chr1", 10, 20, -2.0),
                    (b"chr2", 0, 5, 7.5)])
    assert t.maxvalue == 7.5
    assert t.minvalue == -2.0


def test_bedGraphTrackI_add_loc_many_chromosomes_are_independent():
    t = make_track([(b"chr2", 0, 10, 1.0), (b"chr1", 0, 5, 2.0),
                    (b"chr2", 10, 30, 3.0), (b"chrX", 7, 9, 4.0)])
    assert t.get_chr_names() == {b"chr1", b"chr2", b"chrX"}
    assert track_content(t) == {b"chr1": ([5], [2.0]),
                                b"chr2": ([10, 30], [1.0, 3.0]),
                                b"chrX": ([7, 9], [0.0, 4.0])}


@pytest.mark.parametrize("value, stored", [
    (0.1, f32(0.1)),                    # rounded to float32
    (3.4e38, f32(3.4e38)),              # near float32 maximum
    (1e39, math.inf),                   # beyond float32 range
    (-1e39, -math.inf),
    (1e-46, 0.0),                       # below the smallest subnormal
])
def test_bedGraphTrackI_add_loc_float32_values(value, stored):
    t = bedGraphTrackI()
    t.add_loc(b"chr1", 0, 10, value)
    assert list(t.get_data_by_chr(b"chr1")[1]) == [stored]


def test_bedGraphTrackI_add_loc_int32_max_end():
    t = bedGraphTrackI()
    t.add_loc(b"chr1", 0, 2 ** 31 - 1, 1.0)
    assert track_content(t) == {b"chr1": ([2 ** 31 - 1], [1.0])}


def test_bedGraphTrackI_add_loc_end_beyond_int32_raises():
    t = bedGraphTrackI()
    with pytest.raises(OverflowError, match="value too large to convert to int"):
        t.add_loc(b"chr1", 0, 2 ** 31, 1.0)


def test_bedGraphTrackI_add_loc_str_chromosome_raises():
    t = bedGraphTrackI()
    with pytest.raises(TypeError, match="expected bytes, got str"):
        t.add_loc("chr1", 0, 10, 1.0)


# ------------------------------------
# bedGraphTrackI.add_loc_wo_merge
# ------------------------------------

@pytest.mark.parametrize("rows, baseline, expected", [
    # identical neighbours are NOT merged
    ([(0, 10, 1.0), (10, 20, 1.0)], 0, ([10, 20], [1.0, 1.0])),
    # values below the baseline are raised to the baseline
    ([(0, 10, 1.0), (10, 20, 5.0)], 2, ([10, 20], [2.0, 5.0])),
    # first region after 0 gets a baseline block
    ([(5, 10, 1.0)], 0, ([5, 10], [0.0, 1.0])),
    ([(-5, 10, 1.0)], 0, ([10], [1.0])),
    ([(0, 10, -3.0)], 0, ([10], [0.0])),
])
def test_bedGraphTrackI_add_loc_wo_merge(rows, baseline, expected):
    t = bedGraphTrackI(baseline_value=baseline)
    for (s, e, v) in rows:
        t.add_loc_wo_merge(b"chr1", s, e, v)
    assert track_content(t) == {b"chr1": expected}


def test_bedGraphTrackI_add_loc_wo_merge_ignores_nonpositive_end():
    t = bedGraphTrackI()
    t.add_loc_wo_merge(b"chr1", -5, 0, 1.0)
    assert t.get_chr_names() == set()


def test_bedGraphTrackI_add_loc_wo_merge_max_min_use_clamped_value():
    t = bedGraphTrackI(baseline_value=1)
    t.add_loc_wo_merge(b"chr1", 0, 10, -5.0)   # stored as 1.0
    t.add_loc_wo_merge(b"chr1", 10, 20, 4.0)
    assert t.maxvalue == 4.0
    assert t.minvalue == 1.0


# ------------------------------------
# bedGraphTrackI.add_chrom_data / add_chrom_data_PV
# ------------------------------------

def test_bedGraphTrackI_add_chrom_data_stores_arrays():
    p = array("i", [10, 20, 30])
    v = array("f", [1.0, -2.0, 6.0])
    t = bedGraphTrackI()
    t.add_chrom_data(b"chr1", p, v)
    data = t.get_data_by_chr(b"chr1")
    # the arrays are stored as given, not copied
    assert data[0] is p and data[1] is v
    assert t.maxvalue == 6.0 and t.minvalue == -2.0
    assert t.total() == 3


def test_bedGraphTrackI_add_chrom_data_replaces_chromosome():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0)])
    t.add_chrom_data(b"chr1", array("i", [5]), array("f", [9.0]))
    assert track_content(t) == {b"chr1": ([5], [9.0])}
    assert t.maxvalue == 9.0


def test_bedGraphTrackI_add_chrom_data_empty_raises():
    t = bedGraphTrackI()
    with pytest.raises(ValueError, match="empty"):
        t.add_chrom_data(b"chr1", array("i"), array("f"))


def test_bedGraphTrackI_add_chrom_data_PV():
    pv = np.array([(10, 1.5), (25, -0.5), (40, 3.0)],
                  dtype=[("p", "i4"), ("v", "f4")])
    t = bedGraphTrackI()
    t.add_chrom_data_PV(b"chr1", pv)
    (p, v) = t.get_data_by_chr(b"chr1")
    assert (p.typecode, v.typecode) == ("i", "f")
    assert (list(p), list(v)) == ([10, 25, 40], [1.5, -0.5, 3.0])
    assert t.maxvalue == 3.0 and t.minvalue == -0.5


def test_bedGraphTrackI_add_chrom_data_PV_missing_field_raises():
    pv = np.array([(10, 1.5)], dtype=[("x", "i4"), ("v", "f4")])
    with pytest.raises(ValueError, match="no field of name p"):
        bedGraphTrackI().add_chrom_data_PV(b"chr1", pv)


# ------------------------------------
# bedGraphTrackI.destroy / get_data_by_chr / get_chr_names / total
# ------------------------------------

def test_bedGraphTrackI_destroy():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr2", 0, 10, 2.0)])
    assert t.destroy() is True
    assert t.get_chr_names() == set()
    assert t.total() == 0
    # the track is usable again afterwards
    t.add_loc(b"chr3", 0, 5, 1.0)
    assert track_content(t) == {b"chr3": ([5], [1.0])}


def test_bedGraphTrackI_destroy_empty():
    assert bedGraphTrackI().destroy() is True


def test_bedGraphTrackI_get_data_by_chr_missing_is_empty_list():
    t = make_track([(b"chr1", 0, 10, 1.0)])
    assert t.get_data_by_chr(b"chr2") == []


def test_bedGraphTrackI_get_chr_names_returns_new_set():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr2", 0, 10, 1.0)])
    names = t.get_chr_names()
    names.add(b"chrZ")
    assert t.get_chr_names() == {b"chr1", b"chr2"}


@pytest.mark.parametrize("rows, expected", [
    ([], 0),
    ([(b"chr1", 0, 10, 1.0)], 1),
    # merged neighbours count once
    ([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 1.0)], 1),
    # baseline block counts as a region
    ([(b"chr1", 5, 10, 1.0)], 2),
    ([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0),
      (b"chr2", 0, 10, 1.0)], 3),
])
def test_bedGraphTrackI_total(rows, expected):
    assert make_track(rows).total() == expected


# ------------------------------------
# bedGraphTrackI.filter_score / reset_baseline
# ------------------------------------

def test_bedGraphTrackI_filter_score_sets_low_regions_to_baseline():
    # baseline 1: regions below cutoff 1 become 1, and consecutive
    # baseline regions are merged: [0,10)=5 [10,30)=1 [30,40)=6
    t = make_track([(b"chr1", 0, 10, 5.0), (b"chr1", 10, 20, 0.5),
                    (b"chr1", 20, 30, 0.2), (b"chr1", 30, 40, 6.0)],
                   baseline=1)
    assert t.filter_score(cutoff=1) is True
    assert track_content(t) == {b"chr1": ([10, 30, 40], [5.0, 1.0, 6.0])}


def test_bedGraphTrackI_filter_score_keeps_values_equal_to_cutoff():
    # value == cutoff is kept ("value < cutoff" is filtered)
    t = make_track([(b"chr1", 0, 10, 3.0), (b"chr1", 10, 20, 2.0),
                    (b"chr1", 20, 30, 3.0)], baseline=1)
    t.filter_score(cutoff=2)
    assert track_content(t) == {b"chr1": ([10, 20, 30], [3.0, 2.0, 3.0])}


def test_bedGraphTrackI_filter_score_all_below():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 1.5)],
                   baseline=2)
    t.filter_score(cutoff=2)
    assert track_content(t) == {b"chr1": ([20], [2.0])}


def test_bedGraphTrackI_filter_score_default_keeps_nonnegative_values():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 0.0),
                    (b"chr2", 0, 5, 2.0)], baseline=0)
    t.filter_score()
    assert track_content(t) == {b"chr1": ([10, 20], [1.0, 0.0]),
                                b"chr2": ([5], [2.0])}


def test_bedGraphTrackI_reset_baseline():
    # values < 2 become 2; then equal neighbours are merged:
    # [0,10)=1->2, [10,20)=5, [20,30)=2, [30,40)=1.5->2
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 5.0),
                    (b"chr1", 20, 30, 2.0), (b"chr1", 30, 40, 1.5)])
    assert t.reset_baseline(2.0) is None
    assert t.baseline_value == 2.0
    assert track_content(t) == {b"chr1": ([10, 20, 40], [2.0, 5.0, 2.0])}


def test_bedGraphTrackI_reset_baseline_many_chromosomes():
    t = make_track([(b"chr1", 0, 10, 4.0), (b"chr1", 10, 20, 3.0),
                    (b"chr2", 0, 10, 0.5), (b"chr2", 10, 30, 0.7),
                    (b"chr2", 30, 40, 9.0)])
    t.reset_baseline(3.5)
    assert track_content(t) == {b"chr1": ([10, 20], [4.0, 3.5]),
                                b"chr2": ([30, 40], [3.5, 9.0])}


# ------------------------------------
# bedGraphTrackI.summary
# ------------------------------------

def test_bedGraphTrackI_summary_single_region():
    # sum = 3 * 10, length 10, mean 3, all values equal -> std 0
    assert make_track([(b"chr1", 0, 10, 3.0)]).summary() == \
        (30.0, 10, 3.0, 3.0, 3.0, 0.0)


def test_bedGraphTrackI_summary_sum_length_max_min_mean():
    # chr1: [0,10)=1 [10,20)=3; chr2: [0,5)=2
    # sum = 10 + 30 + 10 = 50; length 25; mean 2
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 3.0),
                    (b"chr2", 0, 5, 2.0)])
    assert t.summary()[:5] == (50.0, 25, 3.0, 1.0, 2.0)


def test_bedGraphTrackI_summary_empty_raises():
    with pytest.raises(ZeroDivisionError, match="float division"):
        bedGraphTrackI().summary()


# ------------------------------------
# bedGraphTrackI.call_peaks
# ------------------------------------

@pytest.mark.parametrize("cutoff, min_length, max_gap, expected", [
    # gap [20,30) = 10 <= 10: one peak [10,40), summit in the 7 region
    (1, 0, 10, [(b"chr1", 10, 40, 35, 7.0)]),
    # gap 10 > 9: two peaks
    (1, 0, 9, [(b"chr1", 10, 20, 15, 5.0), (b"chr1", 30, 40, 35, 7.0)]),
    # min_length is inclusive (10 >= 10)
    (1, 10, 9, [(b"chr1", 10, 20, 15, 5.0), (b"chr1", 30, 40, 35, 7.0)]),
    (1, 11, 9, []),
    (1, 30, 10, [(b"chr1", 10, 40, 35, 7.0)]),
    (1, 31, 10, []),
    # cutoff is inclusive (5 >= 5)
    (5, 0, 0, [(b"chr1", 10, 20, 15, 5.0), (b"chr1", 30, 40, 35, 7.0)]),
    (5.5, 0, 100, [(b"chr1", 30, 40, 35, 7.0)]),
    (7.5, 0, 100, []),
    # a cutoff of 0 takes the whole chromosome
    (0, 0, 0, [(b"chr1", 0, 100, 35, 7.0)]),
])
def test_bedGraphTrackI_call_peaks(cutoff, min_length, max_gap, expected):
    t = make_track(PEAK_ROWS)
    peaks = t.call_peaks(cutoff=cutoff, min_length=min_length,
                         max_gap=max_gap)
    assert peak_rows(peaks) == expected


def test_bedGraphTrackI_call_peaks_other_fields_are_zero():
    peaks = make_track(PEAK_ROWS).call_peaks(cutoff=1, min_length=0,
                                             max_gap=0)
    for p in peaks.get_data_from_chrom(b"chr1"):
        assert (p["pileup"], p["pscore"], p["fc"], p["qscore"],
                p["length"]) == (0.0, 0.0, 0.0, 0.0, 10)


def test_bedGraphTrackI_call_peaks_default_arguments():
    # defaults: cutoff 1, min_length 200, max_gap 50
    t = make_track([(b"chr1", 0, 100, 0.0), (b"chr1", 100, 300, 1.0),
                    (b"chr1", 300, 350, 0.0), (b"chr1", 350, 360, 2.0),
                    (b"chr1", 360, 1000, 0.0)])
    # gap 50 <= 50: [100,360), summit (350 + 360) / 2
    assert peak_rows(t.call_peaks()) == [(b"chr1", 100, 360, 355, 2.0)]


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_bedGraphTrackI_call_peaks_summit_ties(n):
    rows = [(b"chr1", 0, 10, 0.0)]
    pos = 10
    for v in tie_values(n):
        rows.append((b"chr1", pos, pos + 10, v))
        pos += 10
    rows.append((b"chr1", pos, pos + 10, 0.0))
    peaks = make_track(rows).call_peaks(cutoff=1, min_length=0, max_gap=0)
    assert peak_rows(peaks) == [(b"chr1", 10, pos, TIE_SUMMIT[n], 9.0)]


def test_bedGraphTrackI_call_peaks_summit_truncates_half_position():
    # summit of [10,15) is int((10 + 15) / 2) = 12
    t = make_track([(b"chr1", 0, 10, 0.0), (b"chr1", 10, 15, 5.0),
                    (b"chr1", 15, 30, 0.0)])
    assert peak_rows(t.call_peaks(cutoff=1, min_length=0, max_gap=0)) == \
        [(b"chr1", 10, 15, 12, 5.0)]


@pytest.mark.parametrize("rows, expected", [
    # first region above the cutoff: the peak starts at 0
    ([(0, 10, 5.0), (10, 20, 0.0)], [(b"chr1", 0, 10, 5, 5.0)]),
    # last region above the cutoff: the peak ends at the chromosome end
    ([(0, 10, 0.0), (10, 20, 5.0)], [(b"chr1", 10, 20, 15, 5.0)]),
    # whole chromosome above
    ([(0, 20, 5.0)], [(b"chr1", 0, 20, 10, 5.0)]),
    # a region starting after 0 has a baseline block before it
    ([(30, 40, 5.0)], [(b"chr1", 30, 40, 35, 5.0)]),
])
def test_bedGraphTrackI_call_peaks_chromosome_ends(rows, expected):
    t = make_track([(b"chr1",) + r for r in rows])
    assert peak_rows(t.call_peaks(cutoff=1, min_length=0, max_gap=0)) == \
        expected


def test_bedGraphTrackI_call_peaks_many_chromosomes():
    t = make_track([(b"chr2", 0, 10, 0.0), (b"chr2", 10, 20, 3.0),
                    (b"chr1", 0, 50, 4.0), (b"chr3", 0, 10, 0.0)])
    peaks = t.call_peaks(cutoff=1, min_length=0, max_gap=0)
    assert peak_rows(peaks) == [(b"chr1", 0, 50, 25, 4.0),
                                (b"chr2", 10, 20, 15, 3.0)]
    assert peaks.get_chr_names() == {b"chr1", b"chr2"}


def test_bedGraphTrackI_call_peaks_empty_track():
    peaks = bedGraphTrackI().call_peaks(cutoff=1, min_length=0, max_gap=0)
    assert isinstance(peaks, PeakIO)
    assert peaks.total == 0


# ------------------------------------
# bedGraphTrackI.call_broadpeaks
# ------------------------------------

def test_bedGraphTrackI_call_broadpeaks_two_cores_in_one_link():
    # lvl1 (>=5, gap 50): [100,200) and [300,400), apart by 100 > 50
    # lvl2 (>=1, gap 400): [100,400) + [500,600) merged (gap 100) ->
    # [100,600). Cores start at the lvl2 start (no left 1-bp block) but
    # end before the lvl2 end (right 1-bp block at 599).
    t = make_track([(b"chr1", 0, 100, 0.0), (b"chr1", 100, 200, 10.0),
                    (b"chr1", 200, 300, 3.0), (b"chr1", 300, 400, 10.0),
                    (b"chr1", 400, 500, 0.0), (b"chr1", 500, 600, 3.0),
                    (b"chr1", 600, 700, 0.0)])
    bp = t.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=50,
                           lvl1_max_gap=50, lvl2_max_gap=400)
    assert isinstance(bp, BroadPeakIO)
    assert broad_rows(bp) == [(b"chr1", 100, 600, b"100", b"600", 3,
                               b"100,100,1", b"0,200,499", 10.0)]
    p = bp.peaks[b"chr1"][0]
    assert (p["pileup"], p["pscore"], p["fc"], p["qscore"]) == (0, 0, 0, 0)


def test_bedGraphTrackI_call_broadpeaks_core_inside_link():
    # lvl2 [100,400), lvl1 [200,300): 1-bp blocks on both sides
    t = make_track([(b"chr1", 0, 100, 0.0), (b"chr1", 100, 200, 3.0),
                    (b"chr1", 200, 300, 10.0), (b"chr1", 300, 400, 3.0),
                    (b"chr1", 400, 500, 0.0)])
    bp = t.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=50,
                           lvl1_max_gap=50, lvl2_max_gap=400)
    assert broad_rows(bp) == [(b"chr1", 100, 400, b"100", b"400", 3,
                               b"1,100,1", b"0,100,299", 10.0)]


def test_bedGraphTrackI_call_broadpeaks_link_without_core():
    # lvl2 peaks [100,200) and [1000,1100) (gap 800 > 400); only the
    # first contains a lvl1 peak (and equals it). The second is written
    # with two 1-bp blocks at its ends.
    t = make_track([(b"chr1", 0, 100, 0.0), (b"chr1", 100, 200, 10.0),
                    (b"chr1", 200, 1000, 0.0), (b"chr1", 1000, 1100, 3.0),
                    (b"chr1", 1100, 1200, 0.0)])
    bp = t.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=50,
                           lvl1_max_gap=50, lvl2_max_gap=400)
    assert broad_rows(bp) == [
        (b"chr1", 100, 200, b"100", b"200", 1, b"100", b"0", 10.0),
        (b"chr1", 1000, 1100, b"1000", b"1100", 2, b"1,1", b"0,99", 3.0)]


def test_bedGraphTrackI_call_broadpeaks_chromosome_without_core_is_skipped():
    # chromosomes are taken from the lvl1 peaks (as in CallPeakUnit), so
    # chr2, which has only a lvl2 peak, gives no broad peak
    t = make_track([(b"chr1", 0, 100, 10.0), (b"chr2", 0, 100, 3.0)])
    bp = t.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=50,
                           lvl1_max_gap=50, lvl2_max_gap=400)
    assert broad_rows(bp) == [(b"chr1", 0, 100, b"0", b"100", 1, b"100",
                               b"0", 10.0)]


def test_bedGraphTrackI_call_broadpeaks_empty_track():
    assert bedGraphTrackI().call_broadpeaks().peaks == {}


@pytest.mark.parametrize("kwargs, message", [
    (dict(lvl1_cutoff=1, lvl2_cutoff=1),
     "level 1 cutoff should be larger than level 2."),
    (dict(lvl1_cutoff=1, lvl2_cutoff=2),
     "level 1 cutoff should be larger than level 2."),
    (dict(lvl1_cutoff=5, lvl2_cutoff=1, lvl1_max_gap=50, lvl2_max_gap=50),
     "level 2 maximum gap should be larger than level 1."),
])
def test_bedGraphTrackI_call_broadpeaks_bad_arguments(kwargs, message):
    t = make_track([(b"chr1", 0, 10, 1.0)])
    with pytest.raises(AssertionError, match=message):
        t.call_broadpeaks(**kwargs)


# ------------------------------------
# bedGraphTrackI.refine_peaks
# ------------------------------------

@pytest.mark.parametrize("rows, peaks, expected", [
    # peak inside one region: summit at the peak middle
    ([(0, 100, 5.0), (100, 200, 1.0)], [(10, 20)],
     [(b"chr1", 10, 20, 15, 5.0)]),
    # peak over three regions: parts are clipped to the peak and the
    # summit is the middle of the highest part
    ([(0, 10, 1.0), (10, 20, 4.0), (20, 30, 2.0)], [(5, 25)],
     [(b"chr1", 5, 25, 15, 4.0)]),
    # a peak running past the end of the data is cut at the data end
    ([(0, 50, 2.0)], [(40, 80)], [(b"chr1", 40, 50, 45, 2.0)]),
    # a highest part clipped by the peak start: summit of [12,20)
    ([(0, 10, 0.0), (10, 20, 5.0), (20, 30, 0.0)], [(12, 25)],
     [(b"chr1", 12, 25, 16, 5.0)]),
])
def test_bedGraphTrackI_refine_peaks(rows, peaks, expected):
    t = make_track([(b"chr1",) + r for r in rows])
    pk = PeakIO()
    for (s, e) in peaks:
        pk.add(b"chr1", s, e)
    assert peak_rows(t.refine_peaks(pk)) == expected


def test_bedGraphTrackI_refine_peaks_sorts_input():
    t = make_track([(b"chr1", 0, 40, 1.0), (b"chr1", 40, 100, 4.0)])
    pk = PeakIO()
    pk.add(b"chr1", 50, 60)
    pk.add(b"chr1", 10, 20)
    new = t.refine_peaks(pk)
    assert peak_rows(new) == [(b"chr1", 10, 20, 15, 1.0),
                              (b"chr1", 50, 60, 55, 4.0)]
    # the input PeakIO is sorted in place
    assert pk.CO_sorted is True
    assert [p["start"] for p in pk.get_data_from_chrom(b"chr1")] == [10, 50]


def test_bedGraphTrackI_refine_peaks_only_common_chromosomes():
    t = make_track([(b"chr1", 0, 100, 2.0), (b"chr2", 0, 100, 3.0)])
    pk = PeakIO()
    pk.add(b"chr2", 10, 30)
    pk.add(b"chr3", 10, 20)
    assert peak_rows(t.refine_peaks(pk)) == [(b"chr2", 10, 30, 20, 3.0)]


def test_bedGraphTrackI_refine_peaks_empty_peaks():
    t = make_track([(b"chr1", 0, 100, 2.0)])
    assert t.refine_peaks(PeakIO()).total == 0


def test_bedGraphTrackI_refine_peaks_not_peakio_raises():
    t = make_track([(b"chr1", 0, 100, 2.0)])
    # a list has .sort(), so the isinstance assertion is what fails
    with pytest.raises(AssertionError):
        t.refine_peaks([1, 2])
    with pytest.raises(AttributeError, match="has no attribute 'sort'"):
        t.refine_peaks(None)


# ------------------------------------
# bedGraphTrackI.set_single_value
# ------------------------------------

def test_bedGraphTrackI_set_single_value():
    t = make_track([(b"chr1", 5, 10, 1.0), (b"chr1", 10, 30, 2.0),
                    (b"chr2", 0, 7, 3.0)])
    new = t.set_single_value(4.0)
    assert isinstance(new, bedGraphTrackI)
    # one region [0, last end) per chromosome
    assert track_content(new) == {b"chr1": ([30], [4.0]),
                                  b"chr2": ([7], [4.0])}
    assert (new.maxvalue, new.minvalue) == (4.0, 4.0)
    # the original is unchanged
    assert track_content(t) == {b"chr1": ([5, 10, 30], [0.0, 1.0, 2.0]),
                                b"chr2": ([7], [3.0])}


def test_bedGraphTrackI_set_single_value_empty():
    assert bedGraphTrackI().set_single_value(1.0).get_chr_names() == set()


# ------------------------------------
# bedGraphTrackI.overlie (covers mean_func, fisher_func, subtract_func,
# divide_func, product_func)
# ------------------------------------

# A: [0,10)=1 [10,30)=3 ; B: [0,20)=2 [20,25)=4
# shared intervals (stop at the shorter track, 25):
#   [0,10)=(1,2)  [10,20)=(3,2)  [20,25)=(3,4)
OV_A = [(b"chr1", 0, 10, 1.0), (b"chr1", 10, 30, 3.0)]
OV_B = [(b"chr1", 0, 20, 2.0), (b"chr1", 20, 25, 4.0)]
# C: [0,5)=1 [5,40)=2 ; with A and B:
#   [0,5)=(1,2,1) [5,10)=(1,2,2) [10,20)=(3,2,2) [20,25)=(3,4,2)
OV_C = [(b"chr1", 0, 5, 1.0), (b"chr1", 5, 40, 2.0)]


@pytest.mark.parametrize("func, expected", [
    ("max", ([10, 20, 25], [2.0, 3.0, 4.0])),
    ("sum", ([10, 20, 25], [3.0, 5.0, 7.0])),
    ("product", ([10, 20, 25], [2.0, 6.0, 12.0])),
    ("mean", ([10, 20, 25], [1.5, 2.5, 3.5])),
    # subtract is (other track) - (self)
    ("subtract", ([10, 20, 25], [1.0, -1.0, 1.0])),
    ("fisher", ([10, 20, 25], [f32(fisher_ref([1, 2])),
                               f32(fisher_ref([3, 2])),
                               f32(fisher_ref([3, 4]))])),
])
def test_bedGraphTrackI_overlie_two_tracks(func, expected):
    a = make_track(OV_A)
    b = make_track(OV_B)
    got = track_content(a.overlie([b], func=func))
    assert list(got) == [b"chr1"]
    assert got[b"chr1"][0] == expected[0]
    assert got[b"chr1"][1] == pytest.approx(expected[1], rel=1e-6)


@pytest.mark.parametrize("func, expected", [
    # the first two intervals both give 2 and are merged
    ("max", ([10, 20, 25], [2.0, 3.0, 4.0])),
    ("sum", ([5, 10, 20, 25], [4.0, 5.0, 7.0, 9.0])),
    ("product", ([5, 10, 20, 25], [2.0, 4.0, 12.0, 24.0])),
    ("mean", ([5, 10, 20, 25], [f32(4 / 3), f32(5 / 3), f32(7 / 3), 3.0])),
    ("fisher", ([5, 10, 20, 25], [f32(fisher_ref([1, 2, 1])),
                                  f32(fisher_ref([1, 2, 2])),
                                  f32(fisher_ref([3, 2, 2])),
                                  f32(fisher_ref([3, 4, 2]))])),
])
def test_bedGraphTrackI_overlie_three_tracks(func, expected):
    a = make_track(OV_A)
    got = track_content(a.overlie([make_track(OV_B), make_track(OV_C)],
                                  func=func))
    assert got[b"chr1"][0] == expected[0]
    assert got[b"chr1"][1] == pytest.approx(expected[1], rel=1e-6)


def test_bedGraphTrackI_overlie_docstring_example():
    # the example in the overlie docstring, func "max"
    a = make_track([(b"chr1", 0, 100, 0.0), (b"chr1", 100, 200, 3.0),
                    (b"chr1", 200, 300, 4.0)])
    b = make_track([(b"chr1", 0, 150, 1.0), (b"chr1", 150, 250, 2.0),
                    (b"chr1", 250, 300, 4.0)])
    assert track_content(a.overlie([b])) == \
        {b"chr1": ([100, 200, 300], [1.0, 3.0, 4.0])}


def test_bedGraphTrackI_overlie_shared_breakpoints_and_merging():
    # both tracks change at 10; sums are 3 and 3, merged into one region
    a = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0)])
    b = make_track([(b"chr1", 0, 10, 2.0), (b"chr1", 10, 20, 1.0)])
    assert track_content(a.overlie([b], func="sum")) == \
        {b"chr1": ([20], [3.0])}


def test_bedGraphTrackI_overlie_fisher_of_zeros_is_zero():
    a = make_track([(b"chr1", 0, 10, 0.0)])
    b = make_track([(b"chr1", 0, 10, 0.0)])
    assert track_content(a.overlie([b], func="fisher")) == \
        {b"chr1": ([10], [0.0])}


def test_bedGraphTrackI_overlie_only_common_chromosomes():
    a = make_track([(b"chr1", 0, 10, 1.0), (b"chr2", 0, 10, 1.0)])
    b = make_track([(b"chr2", 0, 20, 5.0), (b"chr3", 0, 10, 1.0)])
    assert track_content(a.overlie([b], func="sum")) == \
        {b"chr2": ([10], [6.0])}


def test_bedGraphTrackI_overlie_no_common_chromosome():
    a = make_track([(b"chr1", 0, 10, 1.0)])
    b = make_track([(b"chr2", 0, 10, 1.0)])
    assert a.overlie([b]).get_chr_names() == set()


def test_bedGraphTrackI_overlie_leaves_inputs_unchanged():
    a = make_track(OV_A)
    b = make_track(OV_B)
    a.overlie([b], func="sum")
    assert track_content(a) == {b"chr1": ([10, 30], [1.0, 3.0])}
    assert track_content(b) == {b"chr1": ([20, 25], [2.0, 4.0])}


@pytest.mark.parametrize("tracks, func, exc, message", [
    ([], "max", AssertionError, "Specify at least one more bdg objects."),
    ([1], "max", AssertionError, "bdgTrack1 is not a bedGraphTrackI object"),
    (["T", 1], "max", AssertionError,
     "bdgTrack2 is not a bedGraphTrackI object"),
    (["T", "T"], "subtract", Exception,
     "Only one more bdg object is allowed, but provided 2"),
    (["T", "T"], "divide", Exception,
     "Only one more bdg object is allowed, but provided 2"),
    (["T"], "median", Exception, "Invalid function"),
])
def test_bedGraphTrackI_overlie_errors(tracks, func, exc, message):
    a = make_track(OV_A)
    tracks = [make_track(OV_B) if x == "T" else x for x in tracks]
    with pytest.raises(exc, match=message):
        a.overlie(tracks, func=func)


def test_bedGraphTrackI_overlie_requires_a_list():
    a = make_track(OV_A)
    with pytest.raises(TypeError, match="has no len"):
        a.overlie(make_track(OV_B))


# ------------------------------------
# bedGraphTrackI.apply_func
# ------------------------------------

def test_bedGraphTrackI_apply_func():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0),
                    (b"chr2", 0, 5, -3.0)])
    assert t.apply_func(lambda x: x * 2 + 1) is True
    assert track_content(t) == {b"chr1": ([10, 20], [3.0, 5.0]),
                                b"chr2": ([5], [-5.0])}
    assert (t.maxvalue, t.minvalue) == (5.0, -5.0)


def test_bedGraphTrackI_apply_func_does_not_merge():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0)])
    t.apply_func(lambda x: 7.0)
    assert track_content(t) == {b"chr1": ([10, 20], [7.0, 7.0])}


def test_bedGraphTrackI_apply_func_result_stored_as_float32():
    t = make_track([(b"chr1", 0, 10, 1.0)])
    t.apply_func(lambda x: x / 3)
    assert list(t.get_data_by_chr(b"chr1")[1]) == [f32(1 / 3)]


def test_bedGraphTrackI_apply_func_error_propagates():
    t = make_track([(b"chr1", 0, 10, 0.0)])
    with pytest.raises(ZeroDivisionError):
        t.apply_func(lambda x: 1 / x)


# ------------------------------------
# bedGraphTrackI.p2q
# ------------------------------------
#
# Reference: Benjamini-Hochberg on base pairs, as in
# ScoreTrackII.make_pq_table: p-scores are sorted in decreasing order;
# a p-score v whose block starts at rank k (1 + bp of all higher
# p-scores) out of N bp gets q = v + log10(k / N), made monotone and
# clipped at 0.

def bh_qscores(pscore_lengths):
    """{pscore: qscore} for a {pscore: total length} dict."""
    n = sum(pscore_lengths.values())
    k = 1
    pre_q = math.inf
    out = {}
    for v in sorted(pscore_lengths, reverse=True):
        q = min(pre_q, v + math.log10(k / n))
        q = max(q, 0.0)
        out[v] = q
        pre_q = q
        k += pscore_lengths[v]
    return out


def test_bedGraphTrackI_p2q_single_value():
    # N = 1000, k = 1: q = 5 - 3 = 2
    t = make_track([(b"chr1", 0, 1000, 5.0)])
    assert bh_qscores({5.0: 1000}) == {5.0: 2.0}
    assert t.p2q() is None
    assert track_content(t) == {b"chr1": ([1000], [2.0])}


def test_bedGraphTrackI_p2q_low_scores_clip_to_zero_and_merge():
    # q(10) = 10 - 3 = 7; q(0.5) and q(0.2) are below 0 -> 0, and the
    # two zero regions are merged
    t = make_track([(b"chr1", 0, 100, 10.0), (b"chr1", 100, 600, 0.5),
                    (b"chr1", 600, 1000, 0.2)])
    ref = bh_qscores({10.0: 100, 0.5: 500, 0.2: 400})
    assert ref == {10.0: 7.0, 0.5: 0.0, 0.2: 0.0}
    t.p2q()
    assert track_content(t) == {b"chr1": ([100, 1000], [7.0, 0.0])}


def test_bedGraphTrackI_p2q_counts_all_chromosomes():
    # N = 100 + 900 = 1000 over two chromosomes
    t = make_track([(b"chr1", 0, 100, 10.0), (b"chr2", 0, 900, 0.0)])
    t.p2q()
    assert track_content(t) == {b"chr1": ([100], [7.0]),
                                b"chr2": ([900], [0.0])}


def test_bedGraphTrackI_p2q_empty_track():
    t = bedGraphTrackI()
    t.p2q()
    assert t.get_chr_names() == set()


# ------------------------------------
# bedGraphTrackI.extract_value
# ------------------------------------

# self: [0,100)=5 [100,200)=1 ; regions: [0,50)=1 [50,150)=0 [150,200)=2
# regions with value > 0 are [0,50) and [150,200); self is 5 and 1 there
EV_SELF = [(b"chr1", 0, 100, 5.0), (b"chr1", 100, 200, 1.0)]
EV_REGIONS = [(b"chr1", 0, 50, 1.0), (b"chr1", 50, 150, 0.0),
              (b"chr1", 150, 200, 2.0)]


def test_bedGraphTrackI_extract_value_values_and_lengths():
    ret = make_track(EV_SELF).extract_value(make_track(EV_REGIONS))
    assert len(ret) == 3
    assert list(ret[1]) == [5.0, 1.0]
    assert list(ret[2]) == [50, 50]
    assert len(ret[0]) == 2


def test_bedGraphTrackI_extract_value_region_names():
    """Pins the current output.

    Region names are built with str(chrom) on the bytes chromosome name,
    so they read "b'chr1'.0.50". The name format is not documented and
    extract_value has no caller in MACS3, so whether a decoded name was
    intended cannot be told.
    """
    ret = make_track(EV_SELF).extract_value(make_track(EV_REGIONS))
    assert ret[0] == ["b'chr1'.0.50", "b'chr1'.150.200"]


def test_bedGraphTrackI_extract_value_no_common_chromosome():
    ret = make_track(EV_SELF).extract_value(
        make_track([(b"chr2", 0, 10, 1.0)]))
    assert (ret[0], list(ret[1]), list(ret[2])) == ([], [], [])


def test_bedGraphTrackI_extract_value_not_a_track_raises():
    with pytest.raises(AssertionError, match="not a bedGraphTrackI object"):
        make_track(EV_SELF).extract_value(None)


# ------------------------------------
# bedGraphTrackI.extract_value_hmmr
# ------------------------------------

def test_bedGraphTrackI_extract_value_hmmr():
    # signal chr1: [0,10)=1 [10,20)=2 [20,30)=3 ; chrd: [0,30)=9
    # bins chr1 (add_loc_wo_merge): [0,5)=0 [5,15)=2 [15,25)=2
    #      chrd: [0,10)=1
    # Each bin end with a nonzero bin value is reported with the signal
    # value of the signal region containing it; chromosomes come out in
    # reverse sorted order.
    sig = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0),
                      (b"chr1", 20, 30, 3.0), (b"chrd", 0, 30, 9.0)])
    bins = bedGraphTrackI()
    for (s, e, v) in [(0, 5, 0.0), (5, 15, 2.0), (15, 25, 2.0)]:
        bins.add_loc_wo_merge(b"chr1", s, e, v)
    bins.add_loc_wo_merge(b"chrd", 0, 10, 1.0)
    (pos, val, nbin) = sig.extract_value_hmmr(bins)
    assert pos == [(b"chrd", 10), (b"chr1", 15), (b"chr1", 25)]
    assert list(val) == [9.0, 2.0, 3.0]
    assert list(nbin) == [1, 2, 2]
    assert (val.typecode, nbin.typecode) == ("f", "i")


def test_bedGraphTrackI_extract_value_hmmr_negative_bin_value_is_kept():
    sig = make_track([(b"chr1", 0, 20, 4.0)])
    bins = bedGraphTrackI(baseline_value=-5)
    bins.add_loc_wo_merge(b"chr1", 0, 10, -1.0)
    bins.add_loc_wo_merge(b"chr1", 10, 15, 0.0)
    (pos, val, nbin) = sig.extract_value_hmmr(bins)
    assert pos == [(b"chr1", 10)]
    assert (list(val), list(nbin)) == ([4.0], [-1])


def test_bedGraphTrackI_extract_value_hmmr_not_a_track_raises():
    with pytest.raises(AssertionError, match="not a bedGraphTrackI object"):
        make_track(EV_SELF).extract_value_hmmr([])


# ------------------------------------
# bedGraphTrackI.make_ScoreTrackII_for_macs
# ------------------------------------

def test_bedGraphTrackI_make_ScoreTrackII_for_macs():
    # treat chr1: [0,10)=1 [10,30)=2 ; ctrl chr1: [0,20)=5 [20,25)=6
    # union of ends up to the shorter track: 10, 20, 25
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 30, 2.0),
                    (b"only_t", 0, 10, 1.0)])
    c = make_track([(b"chr1", 0, 20, 5.0), (b"chr1", 20, 25, 6.0),
                    (b"only_c", 0, 10, 1.0)])
    s = t.make_ScoreTrackII_for_macs(c, depth1=2.0, depth2=4.0)
    assert isinstance(s, ScoreTrackII)
    assert s.get_chr_names() == {b"chr1"}
    (pos, treat, ctrl, score) = s.get_data_by_chr(b"chr1")
    assert pos.tolist() == [10, 20, 25]
    assert treat.tolist() == [1.0, 2.0, 2.0]
    assert ctrl.tolist() == [5.0, 5.0, 6.0]
    assert score.tolist() == [0.0, 0.0, 0.0]
    assert (pos.dtype, treat.dtype) == (np.int32, np.float32)
    # the depths are used for SPMR: treat / depth1
    s.change_score_method(ord("m"))
    assert s.get_data_by_chr(b"chr1")[3].tolist() == [0.5, 1.0, 1.0]


def test_bedGraphTrackI_make_ScoreTrackII_for_macs_shared_breakpoints():
    t = make_track([(b"chr1", 0, 10, 1.0), (b"chr1", 10, 20, 2.0)])
    c = make_track([(b"chr1", 0, 10, 3.0), (b"chr1", 10, 20, 4.0)])
    (pos, treat, ctrl, _) = t.make_ScoreTrackII_for_macs(c).get_data_by_chr(
        b"chr1")
    assert (pos.tolist(), treat.tolist(), ctrl.tolist()) == \
        ([10, 20], [1.0, 2.0], [3.0, 4.0])


def test_bedGraphTrackI_make_ScoreTrackII_for_macs_not_a_track_raises():
    with pytest.raises(AssertionError,
                       match="bdgTrack2 is not a bedGraphTrackI object"):
        make_track(OV_A).make_ScoreTrackII_for_macs(None)


# ------------------------------------
# bedGraphTrackI.cutoff_analysis
# ------------------------------------
#
# Cutoffs are np.arange(max(min_score, minvalue), min(maxvalue,
# max_score), step) with step = range / steps; a region is above a
# cutoff when its value is strictly greater. Rows are written from the
# highest cutoff down, and cutoffs with no peak are left out.

CA_HEADER = "score\tnpeaks\tlpeaks\tavelpeak\n"
# [0,100)=0 [100,200)=4 [200,300)=0 [300,350)=2 [350,1000)=0
CA_ROWS = [(b"chr1", 0, 100, 0.0), (b"chr1", 100, 200, 4.0),
           (b"chr1", 200, 300, 0.0), (b"chr1", 300, 350, 2.0),
           (b"chr1", 350, 1000, 0.0)]


@pytest.mark.parametrize("kwargs, body", [
    # cutoffs 0,1,2,3; at 0 and 1 both regions, at 2 and 3 only [100,200)
    (dict(max_gap=0, min_length=0, steps=4),
     "3.00\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.00\t2\t150\t75.00\n0.00\t2\t150\t75.00\n"),
    # gap 300 - 200 = 100 <= max_gap: one peak [100,350)
    (dict(max_gap=100, min_length=0, steps=4),
     "3.00\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.00\t1\t250\t250.00\n0.00\t1\t250\t250.00\n"),
    # min_length 100 drops the 50 bp peak
    (dict(max_gap=0, min_length=100, steps=4),
     "3.00\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.00\t1\t100\t100.00\n0.00\t1\t100\t100.00\n"),
    # nothing is long enough: header only
    (dict(max_gap=0, min_length=200, steps=4), ""),
    # sweep from 1 to 3 in steps of 0.5
    (dict(max_gap=0, min_length=0, steps=4, min_score=1, max_score=3),
     "2.50\t1\t100\t100.00\n2.00\t1\t100\t100.00\n"
     "1.50\t2\t150\t75.00\n1.00\t2\t150\t75.00\n"),
])
def test_bedGraphTrackI_cutoff_analysis(kwargs, body):
    assert make_track(CA_ROWS).cutoff_analysis(**kwargs) == CA_HEADER + body


def test_bedGraphTrackI_cutoff_analysis_many_chromosomes():
    # chr2 adds a 100 bp peak of value 3 (above cutoffs 0, 1, 2)
    t = make_track(CA_ROWS + [(b"chr2", 0, 50, 0.0), (b"chr2", 50, 150, 3.0),
                              (b"chr2", 150, 200, 0.0)])
    assert t.cutoff_analysis(max_gap=0, min_length=0, steps=4) == \
        CA_HEADER + ("3.00\t1\t100\t100.00\n2.00\t2\t200\t100.00\n"
                     "1.00\t3\t250\t83.33\n0.00\t3\t250\t83.33\n")


# ------------------------------------
# bedGraphTrackII.__init__
# ------------------------------------

def test_bedGraphTrackII_init_defaults():
    t = bedGraphTrackII()
    assert t.get_chr_names() == set()
    assert t.total() == 0
    assert (t.maxvalue, t.minvalue, t.baseline_value) == \
        (-10000000, 10000000, 0)
    assert bedGraphTrackII(baseline_value=1.5).baseline_value == 1.5


# ------------------------------------
# bedGraphTrackII.add_loc / add_loc_wo_merge
# ------------------------------------

@pytest.mark.parametrize("method", ["add_loc", "add_loc_wo_merge"])
def test_bedGraphTrackII_add_loc_ignores_nonpositive_end(method):
    t = bedGraphTrackII()
    getattr(t, method)(b"chr1", -10, 0, 1.0)
    assert t.get_chr_names() == set()


# ------------------------------------
# bedGraphTrackII.add_chrom_data / finalize / get_data_by_chr /
# get_chr_names / destroy / total
# ------------------------------------

def test_bedGraphTrackII_add_chrom_data_and_finalize():
    pv = np.array([(30, 3.0), (10, -1.0), (20, 2.0)],
                  dtype=[("p", "u4"), ("v", "f4")])
    t = bedGraphTrackII()
    t.add_chrom_data(b"chr1", pv)
    # stored as given until finalize
    assert t.get_data_by_chr(b"chr1") is pv
    assert t.total() == 3
    assert t.finalize() is None
    # finalize sorts by position and sets max/min
    assert t.get_data_by_chr(b"chr1").tolist() == \
        [(10, -1.0), (20, 2.0), (30, 3.0)]
    assert (t.maxvalue, t.minvalue) == (3.0, -1.0)


def test_bedGraphTrackII_many_chromosomes():
    t = make_track2({b"chr1": [(10, 1.0)], b"chr2": [(5, 2.0), (9, 0.0)]})
    assert t.get_chr_names() == {b"chr1", b"chr2"}
    assert t.total() == 3
    assert (t.maxvalue, t.minvalue) == (2.0, 0.0)


def test_bedGraphTrackII_get_data_by_chr_missing_is_none():
    assert make_track2({b"chr1": [(10, 1.0)]}).get_data_by_chr(b"chr2") is None


def test_bedGraphTrackII_finalize_empty_chromosome_raises():
    t = bedGraphTrackII()
    t.add_chrom_data(b"chr1", np.zeros(0, dtype=[("p", "u4"), ("v", "f4")]))
    with pytest.raises(ValueError, match="zero-size array"):
        t.finalize()


def test_bedGraphTrackII_destroy():
    t = make_track2({b"chr1": [(10, 1.0)], b"chr2": [(5, 2.0)]})
    assert t.destroy() is True
    assert t.get_chr_names() == set()
    assert t.total() == 0


# ------------------------------------
# bedGraphTrackII.filter_score
# ------------------------------------

def test_bedGraphTrackII_filter_score_keeps_rows_above_cutoff():
    t = make_track2({b"chr1": [(10, 0.0), (20, 5.0), (30, 1.0), (40, 6.0)],
                     b"chr2": [(10, 2.0), (20, 1.5)]})
    assert t.filter_score(cutoff=1.0) is True
    # strictly greater rows are kept; the others are removed
    assert [v for _, v in t.get_data_by_chr(b"chr1").tolist()] == [5.0, 6.0]
    assert t.get_data_by_chr(b"chr2").tolist() == [(10, 2.0), (20, 1.5)]
    assert t.maxvalue == 6.0


# ------------------------------------
# bedGraphTrackII.summary
# ------------------------------------

def test_bedGraphTrackII_summary_single_region():
    assert make_track2({b"chr1": [(10, 3.0)]}).summary() == \
        (30.0, 10, 3.0, 3.0, 3.0, 0.0)


# ------------------------------------
# bedGraphTrackII.call_peaks
# ------------------------------------

PEAK_ROWS2 = {b"chr1": [(10, 0.0), (20, 5.0), (30, 0.0), (40, 7.0),
                        (100, 0.0)]}


@pytest.mark.parametrize("cutoff, min_length, max_gap, expected", [
    (1, 0, 10, [(b"chr1", 10, 40, 35, 7.0)]),
    (1, 0, 9, [(b"chr1", 10, 20, 15, 5.0), (b"chr1", 30, 40, 35, 7.0)]),
    (1, 10, 9, [(b"chr1", 10, 20, 15, 5.0), (b"chr1", 30, 40, 35, 7.0)]),
    (1, 11, 9, []),
    (5, 0, 0, [(b"chr1", 10, 20, 15, 5.0), (b"chr1", 30, 40, 35, 7.0)]),
    (5.5, 0, 100, [(b"chr1", 30, 40, 35, 7.0)]),
    (7.5, 0, 100, []),
])
def test_bedGraphTrackII_call_peaks(cutoff, min_length, max_gap, expected):
    peaks = make_track2(PEAK_ROWS2).call_peaks(cutoff=cutoff,
                                               min_length=min_length,
                                               max_gap=max_gap)
    assert peak_rows(peaks) == expected


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_bedGraphTrackII_call_peaks_summit_ties(n):
    rows = [(10, 0.0)]
    pos = 10
    for v in tie_values(n):
        pos += 10
        rows.append((pos, v))
    rows.append((pos + 10, 0.0))
    peaks = make_track2({b"chr1": rows}).call_peaks(cutoff=1, min_length=0,
                                                    max_gap=0)
    assert peak_rows(peaks) == [(b"chr1", 10, pos, TIE_SUMMIT[n], 9.0)]


def test_bedGraphTrackII_call_peaks_many_chromosomes_and_empty():
    t = make_track2({b"chr2": [(10, 0.0), (20, 3.0)],
                     b"chr1": [(10, 0.0), (50, 4.0), (60, 0.0)],
                     b"chr3": [(10, 0.0)]})
    assert peak_rows(t.call_peaks(cutoff=1, min_length=0, max_gap=0)) == \
        [(b"chr1", 10, 50, 30, 4.0), (b"chr2", 10, 20, 15, 3.0)]
    assert bedGraphTrackII().call_peaks().total == 0


# ------------------------------------
# bedGraphTrackII.call_broadpeaks
# ------------------------------------

def test_bedGraphTrackII_call_broadpeaks():
    # same signal as test_bedGraphTrackI_call_broadpeaks_two_cores_in_one_link
    t = make_track2({b"chr1": [(100, 0.0), (200, 10.0), (300, 3.0),
                               (400, 10.0), (500, 0.0), (600, 3.0),
                               (700, 0.0)]})
    bp = t.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=50,
                           lvl1_max_gap=50, lvl2_max_gap=400)
    assert broad_rows(bp) == [(b"chr1", 100, 600, b"100", b"600", 3,
                               b"100,100,1", b"0,200,499", 10.0)]


def test_bedGraphTrackII_call_broadpeaks_link_without_core():
    t = make_track2({b"chr1": [(100, 0.0), (200, 3.0), (300, 10.0),
                               (400, 3.0), (500, 0.0), (1000, 0.0),
                               (1100, 3.0), (1200, 0.0)]})
    bp = t.call_broadpeaks(lvl1_cutoff=5, lvl2_cutoff=1, min_length=50,
                           lvl1_max_gap=50, lvl2_max_gap=400)
    assert broad_rows(bp) == [
        (b"chr1", 100, 400, b"100", b"400", 3, b"1,100,1", b"0,100,299",
         10.0),
        (b"chr1", 1000, 1100, b"1000", b"1100", 2, b"1,1", b"0,99", 3.0)]


@pytest.mark.parametrize("kwargs, message", [
    (dict(lvl1_cutoff=1, lvl2_cutoff=1),
     "level 1 cutoff should be larger than level 2."),
    (dict(lvl1_cutoff=5, lvl2_cutoff=1, lvl1_max_gap=60, lvl2_max_gap=50),
     "level 2 maximum gap should be larger than level 1."),
])
def test_bedGraphTrackII_call_broadpeaks_bad_arguments(kwargs, message):
    with pytest.raises(AssertionError, match=message):
        make_track2(PEAK_ROWS2).call_broadpeaks(**kwargs)


# ------------------------------------
# bedGraphTrackII.refine_peaks
# ------------------------------------

