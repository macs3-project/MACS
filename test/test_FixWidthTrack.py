#!/usr/bin/env python
# Time-stamp: <2025-09-29 14:24:50 Tao Liu>

import unittest

from MACS3.Signal.FixWidthTrack import FWTrack

import io
import itertools
from collections import Counter

import numpy as np
import pytest

from MACS3.IO.PeakIO import PeakIO


class Test_FWTrack(unittest.TestCase):

    def setUp(self):

        self.input_regions = [(b"chrY", 0, 0),
                              (b"chrY", 90, 0),
                              (b"chrY", 150, 0),
                              (b"chrY", 70, 0),
                              (b"chrY", 80, 0),
                              (b"chrY", 85, 0),
                              (b"chrY", 85, 0),
                              (b"chrY", 85, 0),
                              (b"chrY", 85, 0),
                              (b"chrY", 90, 1),
                              (b"chrY", 150, 1),
                              (b"chrY", 70, 1),
                              (b"chrY", 80, 1),
                              (b"chrY", 80, 1),
                              (b"chrY", 80, 1),
                              (b"chrY", 85, 1),
                              (b"chrY", 90, 1),
                              ]
        self.fw = 50

    def test_add_loc(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for (c, p, s) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual(fw.total, 17)
        self.assertEqual(fw.length, 17*self.fw)

    def test_filter_dup(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for (c, p, s) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual(fw.total, 17)
        self.assertEqual(fw.length, 17*self.fw)

        # filter out more than 3 tags
        fw.filter_dup(3)
        # one chrY:85:0 should be removed
        self.assertEqual(fw.total, 16)

        # filter out more than 2 tags
        fw.filter_dup(2)
        # then, one chrY:85:0 and one chrY:80:- should be removed
        self.assertEqual(fw.total, 14)

        # filter out more than 1 tag
        fw.filter_dup(1)
        # then, one chrY:85:0 and one chrY:80:1, one chrY:90:1 should be removed
        self.assertEqual(fw.total, 11)

    def test_sample_num(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for (c, p, s) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual(fw.total, 17)
        self.assertEqual(fw.length, 17*self.fw)

        fw.sample_num(10)
        self.assertEqual(fw.total, 9)

    def test_sample_percent(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for (c, p, s) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual(fw.total, 17)
        self.assertEqual(fw.length, 17*self.fw)

        fw.sample_percent(0.5)
        self.assertEqual(fw.total, 8)


# ------------------------------------
# Reference implementation and helpers for the tests below
# ------------------------------------
#
# The reference pileup is the coverage of half-open fragments
# [start, end) accumulated on the elementary segments between the
# sorted distinct coordinates (a compressed coverage array), with
# values max(coverage * scale, baseline) in float32. Results are
# compared as change points (pos, value): segment i covers
# [pos[i-1], pos[i]) with pos[-1] read as 0, equal neighbours merged.

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


def fw_intervals(plus, minus, d, directional=True, end_shift=0,
                 rlength=INT32_MAX):
    """Fragments of 5' cuts. A plus cut at p moves to p + end_shift and
    extends d bp to the right; a minus cut at m moves to m - end_shift
    and extends d bp to the left. Without direction the fragment is
    centred on the moved cut, d//2 bp on its 5' side and d - d//2 bp on
    its 3' side. Both ends are clipped to [0, rlength]."""
    plus = np.asarray(plus, dtype=np.int64) + end_shift
    minus = np.asarray(minus, dtype=np.int64) - end_shift
    if directional:
        starts = np.concatenate((plus, minus - d))
        ends = np.concatenate((plus + d, minus))
    else:
        starts = np.concatenate((plus - d // 2, minus - (d - d // 2)))
        ends = np.concatenate((plus + (d - d // 2), minus + d // 2))
    return np.clip(starts, 0, rlength), np.clip(ends, 0, rlength)


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


def build_track(reads, fw=50, buffer_size=100000):
    t = FWTrack(fw=fw, buffer_size=buffer_size)
    for c, p, s in reads:
        t.add_loc(c, p, s)
    t.finalize()
    return t


def expected_strands(reads):
    """{chrom: (sorted plus positions, sorted minus positions)}."""
    out = {}
    for c, p, s in reads:
        out.setdefault(c, ([], []))[s].append(p)
    return {c: (sorted(a), sorted(b)) for c, (a, b) in out.items()}


def keep_at_most(values, maxnum):
    """Sorted values with every distinct value kept at most maxnum times."""
    return [v for v, n in sorted(Counter(values).items())
            for _ in range(min(n, maxnum))]


def random_track(seed, n_plus=60, n_minus=50, span=3000, chrom=b"chr1",
                 fw=50):
    rng = np.random.default_rng(seed)
    t = FWTrack(fw=fw)
    for p in rng.integers(0, span, n_plus).tolist():
        t.add_loc(chrom, p, 0)
    for p in rng.integers(0, span, n_minus).tolist():
        t.add_loc(chrom, p, 1)
    for p, s in ((0, 0), (2, 1), (15, 1)):     # fragments reaching past 0
        t.add_loc(chrom, p, s)
    t.finalize()
    return t


READS = [(b"chrY", 0, 0), (b"chrY", 90, 0), (b"chrY", 150, 0),
         (b"chrY", 70, 0), (b"chrY", 85, 0), (b"chrY", 85, 0),
         (b"chrY", 90, 1), (b"chrY", 150, 1), (b"chrY", 70, 1),
         (b"chrX", 500, 1), (b"chrX", 20, 1), (b"chrX", 300, 0),
         (b"chr1", 7, 0), (b"chr1", 3, 0), (b"chr1", 7, 0),
         (b"chr1", 1000, 1), (b"chr1", 999, 1), (b"chr1", 7, 1)]

FILTER_READS = ([(b"chrA", p, 0) for p in (20, 5, 10, 5, 20, 5)] +
                [(b"chrA", p, 1) for p in (40, 30, 40, 30, 40, 40)] +
                [(b"chrB", p, 0) for p in (7, 7)] +
                [(b"chrC", p, 0) for p in (3, 1, 3, 2)] +
                [(b"chrC", p, 1) for p in (9, 9, 9)])


# ------------------------------------
# FWTrack.__init__
# ------------------------------------

def test_fwtrack_init_defaults():
    t = FWTrack()
    assert (t.fw, t.total, t.length, t.annotation, t.buffer_size) == \
        (0, 0, 0, "", 100000)
    assert t.get_chr_names() == set()


def test_fwtrack_init_arguments():
    t = FWTrack(fw=36, anno="treat", buffer_size=10)
    assert (t.fw, t.annotation, t.buffer_size) == (36, "treat", 10)


# ------------------------------------
# FWTrack.add_loc and FWTrack.finalize
# ------------------------------------

@pytest.mark.parametrize("buffer_size", [1, 2, 3, 100000])
def test_add_loc_finalize_sorts_each_strand(buffer_size):
    t = build_track(READS, buffer_size=buffer_size)
    exp = expected_strands(READS)
    assert t.get_chr_names() == set(exp)
    for c, (plus, minus) in exp.items():
        gp, gm = t.get_locations_by_chr(c)
        assert gp.dtype == np.int32 and gm.dtype == np.int32
        assert gp.tolist() == plus
        assert gm.tolist() == minus
    assert t.total == len(READS)
    assert t.length == 50 * len(READS)


def test_add_loc_one_strand_only():
    t = build_track([(b"chrA", 30, 1), (b"chrA", 10, 1)])
    plus, minus = t.get_locations_by_chr(b"chrA")
    assert plus.tolist() == []
    assert minus.tolist() == [10, 30]
    assert t.total == 2


def test_add_loc_int32_extremes():
    t = build_track([(b"chrA", INT32_MAX, 0), (b"chrA", -2**31, 1)])
    plus, minus = t.get_locations_by_chr(b"chrA")
    assert plus.tolist() == [INT32_MAX]
    assert minus.tolist() == [-2**31]


@pytest.mark.parametrize("pos", [2**31, -2**31 - 1])
def test_add_loc_position_outside_int32_raises(pos):
    t = FWTrack(fw=50)
    with pytest.raises(OverflowError):
        t.add_loc(b"chrA", pos, 0)


def test_add_loc_strand_other_than_0_or_1_raises():
    t = FWTrack(fw=50)
    with pytest.raises(IndexError):
        t.add_loc(b"chrA", 10, 2)       # new chromosome
    t2 = FWTrack(fw=50)
    t2.add_loc(b"chrA", 10, 0)
    with pytest.raises(IndexError):
        t2.add_loc(b"chrA", 20, 2)      # existing chromosome


def test_add_loc_chromosome_must_be_bytes():
    t = FWTrack(fw=50)
    with pytest.raises(TypeError, match="expected bytes"):
        t.add_loc("chrA", 10, 0)


def test_finalize_empty_track():
    t = FWTrack(fw=50)
    t.finalize()
    assert t.total == 0
    assert t.length == 0


# ------------------------------------
# FWTrack.destroy
# ------------------------------------

def test_destroy_releases_all_chromosomes():
    t = build_track(READS)
    t.destroy()
    assert t.get_chr_names() == set()
    with pytest.raises(Exception, match="No such chromosome name"):
        t.get_locations_by_chr(b"chrY")


# ------------------------------------
# FWTrack.set_rlengths and FWTrack.get_rlengths
# ------------------------------------

def test_get_rlengths_defaults_to_int32_max():
    t = build_track(READS)
    assert t.get_rlengths() == {b"chrY": INT32_MAX, b"chrX": INT32_MAX,
                                b"chr1": INT32_MAX}


def test_set_rlengths_fills_missing_and_ignores_extra():
    t = build_track(READS)
    assert t.set_rlengths({b"chrY": 1000, b"chrZ": 5}) is True
    assert t.get_rlengths() == {b"chrY": 1000, b"chrX": INT32_MAX,
                                b"chr1": INT32_MAX}


def test_get_rlengths_empty_track():
    assert FWTrack().get_rlengths() == {}


# ------------------------------------
# FWTrack.get_locations_by_chr and FWTrack.get_chr_names
# ------------------------------------

def test_get_locations_by_chr_missing_raises():
    t = build_track(READS)
    with pytest.raises(Exception,
                       match=r"No such chromosome name \(b'chrZ'\)"):
        t.get_locations_by_chr(b"chrZ")


def test_get_chr_names_returns_set_of_bytes():
    names = build_track(READS).get_chr_names()
    assert isinstance(names, set)
    assert names == {b"chrY", b"chrX", b"chr1"}


# ------------------------------------
# FWTrack.sort
# ------------------------------------

def test_sort_restores_ascending_order():
    t = build_track(READS)
    for c in t.get_chr_names():
        for arr in t.get_locations_by_chr(c):
            arr[:] = arr[::-1].copy()
    t.sort()
    for c, (plus, minus) in expected_strands(READS).items():
        gp, gm = t.get_locations_by_chr(c)
        assert gp.tolist() == plus
        assert gm.tolist() == minus


# ------------------------------------
# FWTrack.filter_dup
# ------------------------------------

@pytest.mark.parametrize("maxnum", [1, 2, 3, 10])
def test_filter_dup_keeps_at_most_maxnum(maxnum):
    t = build_track(FILTER_READS, fw=20)
    ret = t.filter_dup(maxnum)
    total = 0
    for c, (plus, minus) in expected_strands(FILTER_READS).items():
        gp, gm = t.get_locations_by_chr(c)
        assert gp.tolist() == keep_at_most(plus, maxnum)
        assert gm.tolist() == keep_at_most(minus, maxnum)
        total += len(keep_at_most(plus, maxnum)) + \
            len(keep_at_most(minus, maxnum))
    assert ret == total
    assert t.total == total
    assert t.length == 20 * total


def test_filter_dup_negative_keeps_all():
    # maxnum < 0 is the 'keep all duplicates' setting
    t = build_track(FILTER_READS)
    assert t.filter_dup(-1) == len(FILTER_READS)
    assert t.total == len(FILTER_READS)
    for c, (plus, minus) in expected_strands(FILTER_READS).items():
        gp, gm = t.get_locations_by_chr(c)
        assert gp.tolist() == plus
        assert gm.tolist() == minus


def test_filter_dup_default_keeps_all():
    t = build_track(FILTER_READS)
    assert t.filter_dup() == len(FILTER_READS)


# ------------------------------------
# FWTrack.sample_percent and FWTrack.sample_num
# ------------------------------------

def sample_track():
    rng = np.random.default_rng(3)
    reads = [(b"chr1", int(p), 0) for p in rng.integers(0, 500, 40)]
    reads += [(b"chr1", int(p), 1) for p in rng.integers(0, 500, 33)]
    reads += [(b"chr2", int(p), 0) for p in rng.integers(0, 50, 7)]
    return reads, build_track(reads, fw=30)


def expected_sample_count(n, percent):
    # num = int(round(n * percent, 5)) with percent passed as float32
    return int(round(n * float(np.float32(percent)), 5))


def check_sampled(t, reads, percent):
    total = 0
    for c, strands in expected_strands(reads).items():
        got = t.get_locations_by_chr(c)
        for orig, arr in zip(strands, got):
            values = arr.tolist()
            assert len(values) == expected_sample_count(len(orig), percent)
            assert values == sorted(values)
            assert not Counter(values) - Counter(orig)   # a sub-multiset
            total += len(values)
    assert t.total == total
    assert t.length == t.fw * total


@pytest.mark.parametrize("percent", [0.0, 0.3, 0.5, 0.77, 1.0])
def test_sample_percent_counts_and_subset(percent):
    reads, t = sample_track()
    t.sample_percent(percent, seed=11)
    check_sampled(t, reads, percent)


def test_sample_percent_without_seed_counts():
    reads, t = sample_track()
    t.sample_percent(0.5)
    check_sampled(t, reads, 0.5)


@pytest.mark.parametrize("seed", [0, 7, 12345])
def test_sample_percent_same_seed_same_sample(seed):
    _, t1 = sample_track()
    _, t2 = sample_track()
    t1.sample_percent(0.4, seed=seed)
    t2.sample_percent(0.4, seed=seed)
    for c in t1.get_chr_names():
        for a, b in zip(t1.get_locations_by_chr(c),
                        t2.get_locations_by_chr(c)):
            assert a.tolist() == b.tolist()


def test_sample_percent_one_keeps_everything():
    reads, t = sample_track()
    t.sample_percent(1.0, seed=1)
    for c, (plus, minus) in expected_strands(reads).items():
        gp, gm = t.get_locations_by_chr(c)
        assert gp.tolist() == plus
        assert gm.tolist() == minus


@pytest.mark.parametrize("samplesize", [1, 10, 40, 80])
def test_sample_num_counts(samplesize):
    reads, t = sample_track()
    total = t.total
    # percent = float32(samplesize) / total, computed in float32
    percent = np.float32(np.float32(samplesize) / np.float32(total))
    t.sample_num(samplesize, seed=5)
    check_sampled(t, reads, float(percent))


def test_sample_num_same_seed_same_sample():
    _, t1 = sample_track()
    _, t2 = sample_track()
    t1.sample_num(30, seed=99)
    t2.sample_num(30, seed=99)
    for c in t1.get_chr_names():
        for a, b in zip(t1.get_locations_by_chr(c),
                        t2.get_locations_by_chr(c)):
            assert a.tolist() == b.tolist()


# ------------------------------------
# FWTrack.print_to_bed
# ------------------------------------

def test_print_to_bed_single_chromosome():
    # plus cut p -> [p, p + fw), minus cut m -> [m - fw, m)
    t = build_track([(b"chrA", 30, 0), (b"chrA", 10, 0), (b"chrA", 100, 1),
                     (b"chrA", 60, 1)], fw=25)
    buf = io.StringIO()
    t.print_to_bed(buf)
    assert buf.getvalue() == ("chrA\t10\t35\t.\t.\t+\n"
                              "chrA\t30\t55\t.\t.\t+\n"
                              "chrA\t35\t60\t.\t.\t-\n"
                              "chrA\t75\t100\t.\t.\t-\n")


def test_print_to_bed_many_chromosomes():
    t = build_track(READS, fw=5)
    buf = io.StringIO()
    t.print_to_bed(buf)
    lines = buf.getvalue().splitlines()
    names = [x.split("\t")[0] for x in lines]
    # each chromosome is written as one block
    blocks = [k for k, _ in itertools.groupby(names)]
    assert len(blocks) == len(set(blocks)) == 3
    for c, (plus, minus) in expected_strands(READS).items():
        name = c.decode()
        exp = (["%s\t%d\t%d\t.\t.\t+" % (name, p, p + 5) for p in plus] +
               ["%s\t%d\t%d\t.\t.\t-" % (name, m - 5, m) for m in minus])
        assert [x for x in lines if x.split("\t")[0] == name] == exp


def test_print_to_bed_defaults_to_stdout(capsys):
    t = build_track([(b"chrA", 10, 0)], fw=5)
    t.print_to_bed()
    assert capsys.readouterr().out == "chrA\t10\t15\t.\t.\t+\n"


def test_print_to_bed_requires_positive_fw():
    t = build_track([(b"chrA", 10, 0)], fw=0)
    with pytest.raises(AssertionError, match="should be set larger than 0"):
        t.print_to_bed(io.StringIO())


def test_print_to_bed_requires_file_object():
    t = build_track([(b"chrA", 10, 0)], fw=5)
    with pytest.raises(AssertionError):
        t.print_to_bed(object())


# ------------------------------------
# FWTrack.extract_region_tags
# ------------------------------------

@pytest.mark.parametrize("start, end", [(0, 10), (5, 30), (31, 39),
                                        (70, 90), (151, 200), (-50, 5000)])
def test_extract_region_tags_inclusive_window(start, end):
    reads = [(b"chrA", p, 0) for p in (0, 5, 10, 10, 30, 70, 90, 150)]
    reads += [(b"chrA", p, 1) for p in (4, 5, 30, 40, 90, 91)]
    t = build_track(reads)
    plus, minus = t.extract_region_tags(b"chrA", start, end)
    exp = expected_strands(reads)[b"chrA"]
    assert list(plus) == [p for p in exp[0] if start <= p <= end]
    assert list(minus) == [m for m in exp[1] if start <= m <= end]


def test_extract_region_tags_missing_chromosome():
    t = build_track(READS)
    with pytest.raises(AssertionError, match="can't be found"):
        t.extract_region_tags(b"chrZ", 0, 100)


# ------------------------------------
# FWTrack.compute_region_tags_from_peaks
# ------------------------------------

def collect_tags(chrom, plus, minus, startpos, endpos, name=None,
                 window_size=None, cutoff=None):
    return (chrom, plus.tolist(), minus.tolist(), startpos, endpos, name,
            window_size, cutoff)


def brute_force_tags(track, peaks, window_size, cutoff):
    """Expected callback arguments for each peak, by direct filtering."""
    out = []
    for chrom in sorted(peaks):
        plus, minus = track.get_locations_by_chr(chrom)
        for start, end, name in sorted(peaks[chrom]):
            s, e = start - window_size, end + window_size
            out.append((chrom, [p for p in plus.tolist() if s <= p <= e],
                        [m for m in minus.tolist() if s <= m <= e], s, e,
                        name, window_size, cutoff))
    return out


def make_peakio(peaks):
    pio = PeakIO()
    for chrom, items in peaks.items():
        for start, end, name in items:
            pio.add(chrom, start, end, summit=(start + end) // 2, name=name)
    return pio


SEPARATED_PLUS = [10, 60, 120, 180, 240, 260, 900, 960, 1050, 1149, 1150,
                  1151, 1300, 5100, 5150, 5151]
SEPARATED_MINUS = [55, 130, 199, 251, 940, 1000, 1101, 1200, 4949, 4950,
                   5000]


def separated_track():
    reads = [(b"chr1", p, 0) for p in SEPARATED_PLUS]
    reads += [(b"chr1", m, 1) for m in SEPARATED_MINUS]
    reads += [(b"chr2", 100, 0), (b"chr2", 200, 0), (b"chr2", 150, 1),
              (b"chr2", 300, 1)]
    return build_track(reads)


def test_compute_region_tags_from_peaks_separated_peaks():
    t = separated_track()
    peaks = {b"chr1": [(5000, 5100, b"p3"), (100, 200, b"p1"),
                       (1000, 1100, b"p2")],
             b"chr2": [(120, 180, b"q1")]}
    got = t.compute_region_tags_from_peaks(make_peakio(peaks), collect_tags,
                                           window_size=50, cutoff=2.5)
    assert got == brute_force_tags(t, peaks, 50, 2.5)


def test_compute_region_tags_from_peaks_default_arguments():
    t = separated_track()
    peaks = {b"chr2": [(120, 180, b"q1")]}
    got = t.compute_region_tags_from_peaks(make_peakio(peaks), collect_tags)
    # default window_size 100 and cutoff 5.0: window [20, 280]
    assert got == [(b"chr2", [100, 200], [150], 20, 280, b"q1", 100, 5.0)]


def test_compute_region_tags_from_peaks_missing_chromosome():
    t = separated_track()
    pio = make_peakio({b"chrZ": [(10, 20, b"z")]})
    with pytest.raises(AssertionError, match="can't be found"):
        t.compute_region_tags_from_peaks(pio, collect_tags)


# ------------------------------------
# FWTrack.pileup_a_chromosome
# ------------------------------------

def test_pileup_a_chromosome_hand_example():
    # plus 10 -> [10, 35), minus 30 -> [5, 30)
    t = build_track([(b"chrA", 10, 0), (b"chrA", 30, 1)])
    p, v = t.pileup_a_chromosome(b"chrA", 25)
    np.testing.assert_array_equal(p, [5, 10, 30, 35])
    np.testing.assert_array_equal(v, [0.0, 1.0, 2.0, 1.0])


def test_pileup_a_chromosome_not_directional_hand_example():
    # d 50 centred: plus 100 -> [75, 125), minus 300 -> [275, 325)
    t = build_track([(b"chrA", 100, 0), (b"chrA", 300, 1)])
    res = t.pileup_a_chromosome(b"chrA", 50, directional=False)
    assert_cp_equal(res, ([75, 125, 275, 325], [0.0, 1.0, 0.0, 1.0]))


def test_pileup_a_chromosome_end_shift_hand_example():
    # end_shift 10 moves cuts 3'-ward: plus 100 -> [110, 160),
    # minus 300 -> [240, 290)
    t = build_track([(b"chrA", 100, 0), (b"chrA", 300, 1)])
    res = t.pileup_a_chromosome(b"chrA", 50, end_shift=10)
    assert_cp_equal(res, ([110, 160, 240, 290], [0.0, 1.0, 0.0, 1.0]))


def test_pileup_a_chromosome_extension_past_zero_is_clipped():
    # minus 20 with d 30 -> [-10, 20) clipped to [0, 20); plus 0 -> [0, 30)
    t = build_track([(b"chrA", 0, 0), (b"chrA", 20, 1)])
    res = t.pileup_a_chromosome(b"chrA", 30)
    assert_cp_equal(res, ([20, 30], [2.0, 1.0]))


PILEUP_CASES = [
    # seed, d, directional, end_shift, scale, baseline, rlength
    (0, 200, True, 0, 1.0, 0.0, None),
    (1, 147, True, 0, 0.5, 0.0, None),
    (2, 200, False, 0, 1.0, 0.0, None),
    (3, 201, False, 0, 1.0, 0.0, None),
    (4, 100, True, 20, 1.0, 0.0, None),
    (5, 100, True, -20, 1.0, 0.0, None),
    (6, 151, False, 13, 0.2, 0.1, None),
    (7, 300, True, 0, 1.0, 2.0, 2000),
    (8, 1000, False, 0, 0.01, 0.0, 1500),
]


@pytest.mark.parametrize("seed, d, directional, end_shift, scale, "
                         "baseline, rlength", PILEUP_CASES)
def test_pileup_a_chromosome_matches_reference(seed, d, directional,
                                               end_shift, scale, baseline,
                                               rlength):
    t = random_track(seed)
    if rlength is not None:
        t.set_rlengths({b"chr1": rlength})
    else:
        rlength = INT32_MAX
    res = t.pileup_a_chromosome(b"chr1", d, scale_factor=scale,
                                baseline_value=baseline,
                                directional=directional, end_shift=end_shift)
    plus, minus = t.get_locations_by_chr(b"chr1")
    starts, ends = fw_intervals(plus, minus, d, directional, end_shift,
                                rlength)
    assert_cp_equal(res, ref_pileup(starts, ends, scale, baseline))


def test_pileup_a_chromosome_duplicates_stack():
    t = build_track([(b"chrA", 10, 0)] * 4)
    res = t.pileup_a_chromosome(b"chrA", 5)
    assert_cp_equal(res, ([10, 15], [0.0, 4.0]))


def test_pileup_a_chromosome_int32_extreme():
    t = build_track([(b"chrA", INT32_MAX - 100, 0)])
    res = t.pileup_a_chromosome(b"chrA", 100)
    assert_cp_equal(res, ([INT32_MAX - 100, INT32_MAX], [0.0, 1.0]))


def test_pileup_a_chromosome_independent_of_other_chromosomes():
    one = random_track(4)
    many = random_track(4)
    for p in (5, 50, 500):
        many.add_loc(b"chr2", p, 0)
    many.finalize()
    a = one.pileup_a_chromosome(b"chr1", 120)
    b = many.pileup_a_chromosome(b"chr1", 120)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])


def test_pileup_a_chromosome_missing_chromosome():
    t = build_track(READS)
    with pytest.raises(KeyError):
        t.pileup_a_chromosome(b"chrZ", 100)


# ------------------------------------
# FWTrack.pileup_a_chromosome_c
# ------------------------------------

PILEUP_C_CASES = [
    # ds, scale factors, baseline, directional, end_shift, rlength
    ([200, 1000, 4000], [1.0, 0.2, 0.05], 0.3, False, 0, None),
    ([150, 300], [1.0, 0.5], 0.0, True, 0, None),
    ([50, 51, 52], [1.0, 1.0, 1.0], 0.0, False, 5, None),
    ([100, 2000], [1.0, 0.1], 0.0, False, 0, 2500),
    ([100], [2.0], 0.0, True, 0, None),
]


@pytest.mark.parametrize("ds, sfs, baseline, directional, end_shift, "
                         "rlength", PILEUP_C_CASES)
def test_pileup_a_chromosome_c_matches_reference(ds, sfs, baseline,
                                                 directional, end_shift,
                                                 rlength):
    t = random_track(21)
    if rlength is not None:
        t.set_rlengths({b"chr1": rlength})
    else:
        rlength = INT32_MAX
    res = t.pileup_a_chromosome_c(b"chr1", ds, sfs, baseline_value=baseline,
                                  directional=directional,
                                  end_shift=end_shift)
    plus, minus = t.get_locations_by_chr(b"chr1")
    cps = []
    for d, sf in zip(ds, sfs):
        starts, ends = fw_intervals(plus, minus, d, directional, end_shift,
                                    rlength)
        cps.append(ref_pileup(starts, ends, sf, baseline))
    assert_cp_equal(res, ref_max_truncated(cps))


def test_pileup_a_chromosome_c_single_d_equals_pileup_a_chromosome():
    t = random_track(8)
    a = t.pileup_a_chromosome_c(b"chr1", [180], [0.5], baseline_value=0.2,
                                directional=False)
    b = t.pileup_a_chromosome(b"chr1", 180, scale_factor=0.5,
                              baseline_value=0.2, directional=False)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])


def test_pileup_a_chromosome_c_length_mismatch():
    t = random_track(8)
    with pytest.raises(AssertionError, match="same length"):
        t.pileup_a_chromosome_c(b"chr1", [100, 200], [1.0])
