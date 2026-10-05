#!/usr/bin/env python
# Time-stamp: <2025-09-29 13:48:03 Tao Liu>

import io
import os
import pickle
import random
import string
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np
import pytest

import MACS3
from MACS3.IO.PeakIO import PeakIO
from MACS3.IO.PeakIO import (PeakContent,
                             RegionIO,
                             BroadPeakContent,
                             BroadPeakIO)


class Test_PeakIO(unittest.TestCase):
    def setUp(self):
        self.test_peaks1 = [(b"chrY", 0, 100),
                            (b"chrY", 300, 500),
                            (b"chrY", 700, 900),
                            (b"chrY", 1000, 1200),
                            (b"chrY", 1250, 1450),
                            (b"chrX", 1000, 2000),
                            (b"chrX", 3000, 4000),
                            (b"chr1", 100, 10000),  # chr1 only has one region but overlapping with peaks2
                            (b"chr2", 1000, 2000),  # only peaks1 has chr2
                            (b"chr4", 500, 800),    # chr4 only one region, and not overlapping with peaks2
                            ]
        self.test_peaks2 = [(b"chrY", 100, 200),
                            (b"chrY", 300, 400),
                            (b"chrY", 600, 800),
                            (b"chrY", 1100, 1300),
                            (b"chrY", 1700, 1800),
                            (b"chrX", 1100, 1200),
                            (b"chrX", 1300, 1400),
                            (b"chr1", 2000, 3000),
                            (b"chr3", 1000, 5000),  # only peaks2 has chr3
                            (b"chr4", 1000, 2000),
                            ]
        self.result_exclude2from1 = [(b"chrY", 0, 100),
                                     (b"chrX", 3000, 4000),
                                     (b"chr2", 1000, 2000),
                                     (b"chr4", 500, 800),
                                     ]
        self.exclude2from1 = PeakIO()
        for a in self.result_exclude2from1:
            self.exclude2from1.add(a[0], a[1], a[2])

    def test_exclude(self):
        r1 = PeakIO()
        for a in self.test_peaks1:
            r1.add(a[0], a[1], a[2])
        r2 = PeakIO()
        for a in self.test_peaks2:
            r2.add(a[0], a[1], a[2])
        r1.exclude(r2)
        result = str(r1)
        expected = str(self.exclude2from1)
        # print( "result:\n", result )
        # print( "expected:\n", expected )
        self.assertEqual(result, expected)


# ------------------------------------
# Helpers
# ------------------------------------

def _f32(x):
    """``x`` after a round trip through a C float (the peak fields are float32)."""
    return float(np.float32(x))


def _three_peaks(name=b"MACS3"):
    """Two chromosomes, values exactly representable in float32."""
    p = PeakIO(name=name)
    p.add(b"chr1", 100, 200, summit=150, peak_score=2.5, pileup=7.0,
          pscore=3.25, fold_change=4.5, qscore=2.5)
    p.add(b"chr1", 300, 400, summit=310, peak_score=12.75, pileup=20.0,
          pscore=13.5, fold_change=8.0, qscore=12.75)
    p.add(b"chr2", 0, 50, summit=0, peak_score=0.5, pileup=1.0,
          pscore=0.75, fold_change=1.25, qscore=0.5)
    return p


def _subpeaks(n, end=1000):
    """One peak region (0, end) with ``n`` summits, as call-summits makes."""
    p = PeakIO()
    for i in range(n):
        p.add(b"chr1", 0, end, summit=i, peak_score=1.0)
    return p


def _text(write, *args, **kwargs):
    fh = io.StringIO()
    write(fh, *args, **kwargs)
    return fh.getvalue()


def _run_python(code, cwd):
    """Run ``code`` in a fresh interpreter that imports this MACS3 build."""
    env = os.environ.copy()
    pkg_parent = str(Path(MACS3.__file__).resolve().parent.parent)
    env["PYTHONPATH"] = pkg_parent + (os.pathsep + env["PYTHONPATH"]
                                      if env.get("PYTHONPATH") else "")
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[var] = "1"
    return subprocess.run([sys.executable, "-c", code], cwd=str(cwd), env=env,
                          capture_output=True, text=True, timeout=60)


def _bases(intervals):
    covered = set()
    for s, e in intervals:
        covered.update(range(s, e))
    return covered


XLS_HEADER = ("chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\t"
              "fold_enrichment\t-log10(qvalue)\tname\n")
BROAD_XLS_HEADER = ("chr\tstart\tend\tlength\tpileup\t-log10(pvalue)\t"
                    "fold_enrichment\t-log10(qvalue)\tname\n")


# ------------------------------------
# PeakContent
# ------------------------------------

def test_PeakContent_fields():
    pc = PeakContent(b"chr1", 100, 250, 180, 2.3, 7.5, 3.1, 4.2, 1.9,
                     name=b"p1")
    assert pc["chrom"] == b"chr1"
    assert (pc["start"], pc["end"], pc["length"], pc["summit"]) == \
        (100, 250, 150, 180)
    # the floating fields are stored as C floats
    assert pc["score"] == _f32(2.3)
    assert pc["pileup"] == _f32(7.5)
    assert pc["pscore"] == _f32(3.1)
    assert pc["fc"] == _f32(4.2)
    assert pc["qscore"] == _f32(1.9)
    assert pc["name"] == b"p1"


def test_PeakContent_default_name_and_unknown_key():
    pc = PeakContent(b"chr1", 0, 10, 5, 0, 0, 0, 0, 0)
    assert pc["name"] == b""
    assert pc["no_such_field"] is None


@pytest.mark.parametrize("key,value,expected", [
    ("chrom", b"chrX", b"chrX"),
    ("start", 5, 5),
    ("end", 99, 99),
    ("length", 42, 42),
    ("summit", 7, 7),
    ("score", 0.1, _f32(0.1)),
    ("pileup", 2.25, 2.25),
    ("pscore", 1e-3, _f32(1e-3)),
    ("fc", 3.0, 3.0),
    ("qscore", 0.7, _f32(0.7)),
    ("name", b"x", b"x"),
])
def test_PeakContent_setitem(key, value, expected):
    pc = PeakContent(b"chr1", 0, 10, 5, 1.0, 1.0, 1.0, 1.0, 1.0)
    pc[key] = value
    assert pc[key] == expected


def test_PeakContent_setitem_unknown_key_is_ignored():
    pc = PeakContent(b"chr1", 0, 10, 5, 1.0, 1.0, 1.0, 1.0, 1.0)
    before = pc.__getstate__()
    pc["no_such_field"] = 3
    assert pc.__getstate__() == before


def test_PeakContent_str():
    pc = PeakContent(b"chr1", 10, 20, 15, 2.5, 0, 0, 0, 0)
    # "%s" of a bytes object gives its repr, "%f" gives 6 decimals
    assert str(pc) == "chrom:b'chr1';start:10;end:20;score:2.500000"


def test_PeakContent_getstate_order():
    pc = PeakContent(b"chr1", 10, 30, 15, 2.5, 3.5, 4.5, 5.5, 6.5, name=b"n")
    assert pc.__getstate__() == (b"chr1", 10, 30, 20, 15, 2.5, 3.5, 4.5,
                                 5.5, 6.5, b"n")


def test_PeakContent_pickle_round_trip():
    pc = PeakContent(b"chr1", 10, 30, 15, 2.3, 3.5, 4.5, 5.5, 6.5, name=b"n")
    back = pickle.loads(pickle.dumps(pc))
    assert back.__getstate__() == pc.__getstate__()


def test_PeakContent_int32_extremes():
    pc = PeakContent(b"chr1", 0, 2147483647, 2147483646, 0, 0, 0, 0, 0)
    assert pc["length"] == 2147483647
    with pytest.raises(OverflowError):
        PeakContent(b"chr1", 0, 2147483648, 0, 0, 0, 0, 0, 0)


def test_PeakContent_rejects_str_chrom():
    with pytest.raises(TypeError):
        PeakContent("chr1", 0, 10, 5, 0, 0, 0, 0, 0)


# ------------------------------------
# PeakIO construction, add, add_PeakContent, get_data_from_chrom,
# get_chr_names, sort
# ------------------------------------

def test_PeakIO_defaults():
    p = PeakIO()
    assert p.peaks == {}
    assert p.total == 0
    assert p.CO_sorted is False
    assert p.name == b"MACS3"
    assert PeakIO(name=b"exp").name == b"exp"


def test_add_defaults_and_counts():
    p = PeakIO()
    p.add(b"chr1", 10, 20)
    p.add(b"chr2", 5, 9, summit=7, peak_score=1.5, name=b"x")
    p.add(b"chr1", 0, 4)
    assert p.total == 3
    assert p.get_chr_names() == {b"chr1", b"chr2"}
    first = p.get_data_from_chrom(b"chr1")[0]
    assert first.__getstate__() == (b"chr1", 10, 20, 10, 0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, b"")
    assert [x["start"] for x in p.get_data_from_chrom(b"chr1")] == [10, 0]
    assert p.get_data_from_chrom(b"chr2")[0]["name"] == b"x"


def test_add_resets_sorted_flag():
    p = PeakIO()
    p.add(b"chr1", 10, 20)
    p.sort()
    assert p.CO_sorted is True
    p.add(b"chr1", 0, 5)
    assert p.CO_sorted is False


def test_add_rejects_str_chromosome():
    with pytest.raises(TypeError):
        PeakIO().add("chr1", 0, 10)


def test_add_PeakContent_stores_the_object():
    p = PeakIO()
    pc = PeakContent(b"chr1", 100, 200, 150, 10.0, 5.0, 3.0, 2.0, 1.0)
    p.add_PeakContent(b"chr1", pc)
    assert p.get_data_from_chrom(b"chr1")[0] is pc
    assert p.total == 1
    assert p.CO_sorted is False


def test_add_PeakContent_rejects_other_types():
    with pytest.raises(TypeError):
        PeakIO().add_PeakContent(b"chr1", (b"chr1", 0, 10))


def test_get_data_from_chrom_creates_missing_chromosome():
    p = PeakIO()
    data = p.get_data_from_chrom(b"chr9")
    assert data == []
    assert p.get_chr_names() == {b"chr9"}
    assert p.total == 0
    data.append("sentinel")
    assert p.peaks[b"chr9"] == ["sentinel"]


def test_get_chr_names_empty():
    assert PeakIO().get_chr_names() == set()


def test_sort_by_start_stable_for_ties():
    p = PeakIO()
    p.add(b"chr1", 300, 400)
    p.add(b"chr1", 100, 900, summit=1)
    p.add(b"chr1", 100, 200, summit=2)
    p.add(b"chr2", 50, 60)
    p.add(b"chr2", 10, 20)
    p.sort()
    assert [(x["start"], x["summit"]) for x in p.peaks[b"chr1"]] == \
        [(100, 1), (100, 2), (300, 0)]
    assert [x["start"] for x in p.peaks[b"chr2"]] == [10, 50]
    assert p.CO_sorted is True
    assert p.total == 5


def test_sort_again_after_add():
    p = PeakIO()
    p.add(b"chr1", 300, 400)
    p.sort()
    p.add(b"chr1", 100, 200)
    p.sort()
    assert [x["start"] for x in p.peaks[b"chr1"]] == [100, 300]


def test_str_exact_with_subpeaks():
    p = PeakIO()
    p.add(b"chr1", 10, 20, summit=15, peak_score=1.5)
    p.add(b"chr1", 10, 20, summit=17, peak_score=2.5)
    p.add(b"chr1", 30, 40, summit=35, peak_score=3.0)
    p.add(b"chr2", 5, 9, summit=6, peak_score=0.25)
    assert str(p) == (
        "chrom:chr1\tstart:10\tend:20\tname:peak_1a\tscore:1.5\tsummit:15\n"
        "chrom:chr1\tstart:10\tend:20\tname:peak_1b\tscore:2.5\tsummit:17\n"
        "chrom:chr1\tstart:30\tend:40\tname:peak_2\tscore:3\tsummit:35\n"
        "chrom:chr2\tstart:5\tend:9\tname:peak_3\tscore:0.25\tsummit:6\n")


# ------------------------------------
# PeakIO.randomly_pick
# ------------------------------------

def _pick_pool():
    p = PeakIO()
    for i in range(6):
        p.add(b"chr2", 10 * i, 10 * i + 5)
    for i in range(5):
        p.add(b"chr1", 100 * i, 100 * i + 5)
    return p


def _reference_pick(p, n, seed):
    # all peaks, chromosomes in sorted order, then random.shuffle with seed
    pool = [(c, x["start"]) for c in sorted(p.peaks) for x in p.peaks[c]]
    state = random.getstate()
    try:
        random.seed(seed)
        random.shuffle(pool)
    finally:
        random.setstate(state)
    return pool[:n]


@pytest.mark.parametrize("n,seed", [(4, 12345), (4, 1), (1, 7), (11, 3),
                                    (50, 3)])
def test_randomly_pick_matches_seeded_shuffle(n, seed):
    p = _pick_pool()
    got = p.randomly_pick(n, seed=seed)
    expected = _reference_pick(p, n, seed)
    assert got.total == len(expected)
    for chrom in (b"chr1", b"chr2"):
        assert [x["start"] for x in got.peaks.get(chrom, [])] == \
            [s for c, s in expected if c == chrom]


def test_randomly_pick_reproducible_and_default_seed():
    p = _pick_pool()
    a = p.randomly_pick(5)
    b = p.randomly_pick(5, seed=12345)
    assert str(a) == str(b)
    assert p.total == 11


def test_randomly_pick_shares_peak_objects():
    p = _pick_pool()
    got = p.randomly_pick(11, seed=5)
    originals = {id(x) for c in p.peaks for x in p.peaks[c]}
    assert {id(x) for c in got.peaks for x in got.peaks[c]} == originals


@pytest.mark.parametrize("n", [0, -1])
def test_randomly_pick_requires_positive_n(n):
    with pytest.raises(AssertionError):
        _pick_pool().randomly_pick(n)


# ------------------------------------
# PeakIO.filter_pscore / filter_qscore / filter_fc / filter_score
# ------------------------------------

FILTER_FIELDS = [
    ("filter_pscore", "pscore", "pscore"),
    ("filter_qscore", "qscore", "qscore"),
    ("filter_fc", "fold_change", "fc"),
    ("filter_score", "peak_score", "score"),
]


def _filter_peaks(kwarg, values, chroms=(b"chr1", b"chr2")):
    p = PeakIO()
    for chrom in chroms:
        for i, v in enumerate(values):
            p.add(chrom, 100 * i, 100 * i + 50, **{kwarg: v})
    return p


@pytest.mark.parametrize("method,kwarg,field", FILTER_FIELDS)
def test_filter_lower_bound_inclusive(method, kwarg, field):
    p = _filter_peaks(kwarg, [4.0, 1.0, 2.5])
    getattr(p, method)(2.5)
    assert p.total == 4
    for chrom in (b"chr1", b"chr2"):
        assert [x[field] for x in p.peaks[chrom]] == [4.0, 2.5]


@pytest.mark.parametrize("method,kwarg,field", FILTER_FIELDS)
def test_filter_removes_everything(method, kwarg, field):
    p = _filter_peaks(kwarg, [1.0, 2.0])
    getattr(p, method)(100.0)
    assert p.total == 0
    assert _text(p.write_to_bed, trackline=False) == ""


@pytest.mark.parametrize("method", ["filter_pscore", "filter_qscore",
                                    "filter_fc", "filter_score"])
def test_filter_empty_peakio(method):
    p = PeakIO()
    getattr(p, method)(1.0)
    assert p.total == 0
    assert p.peaks == {}


@pytest.mark.parametrize("method,kwarg,kept", [
    ("filter_pscore", "pscore", False),
    ("filter_qscore", "qscore", False),
    ("filter_fc", "fold_change", True),
    ("filter_score", "peak_score", True),
])
def test_filter_cut_precision(method, kwarg, kept):
    # The stored value is float32(2.3) = 2.2999999523... < 2.3. The pscore
    # and qscore cuts are C doubles (2.3 itself), so the peak is below the
    # cut; the fc and score cuts are C floats (float32(2.3)), so it is equal.
    assert _f32(2.3) < 2.3
    p = PeakIO()
    p.add(b"chr1", 0, 10, **{kwarg: 2.3})
    getattr(p, method)(2.3)
    assert p.total == (1 if kept else 0)


@pytest.mark.parametrize("method,kwarg,field", [
    ("filter_fc", "fold_change", "fc"),
    ("filter_score", "peak_score", "score"),
])
@pytest.mark.parametrize("low,up,expected", [
    pytest.param(2.0, 4.0, [2.0, 3.0], id="upper_exclusive"),
    pytest.param(2.0, 0, [2.0, 3.0, 4.0], id="zero_upper_ignored"),
    pytest.param(2.0, 2.0, [2.0, 3.0, 4.0], id="upper_equal_lower_ignored"),
    pytest.param(3.0, 2.0, [3.0, 4.0], id="upper_below_lower_ignored"),
    pytest.param(5.0, 0, [], id="none_left"),
])
def test_filter_range(method, kwarg, field, low, up, expected):
    p = _filter_peaks(kwarg, [1.0, 2.0, 3.0, 4.0], chroms=(b"chr1",))
    getattr(p, method)(low, up)
    assert [x[field] for x in p.peaks[b"chr1"]] == expected
    assert p.total == len(expected)


# ------------------------------------
# PeakIO.write_to_bed (and the C-only subpeak_letters)
# ------------------------------------

def test_write_to_bed_defaults_exact(tmp_path):
    path = tmp_path / "p.bed"
    with open(path, "w") as fh:
        _three_peaks().write_to_bed(fh)
    assert path.read_text() == (
        'track name="MACS (peaks)" description="MACS" visibility=1\n'
        "chr1\t100\t200\tMACS_peak_1\t2.5\n"
        "chr1\t300\t400\tMACS_peak_2\t12.75\n"
        "chr2\t0\t50\tMACS_peak_3\t0.5\n")


def test_write_to_bed_no_trackline_and_name():
    text = _text(_three_peaks().write_to_bed, name=b"exp", trackline=False)
    assert text == ("chr1\t100\t200\texp_peak_1\t2.5\n"
                    "chr1\t300\t400\texp_peak_2\t12.75\n"
                    "chr2\t0\t50\texp_peak_3\t0.5\n")


@pytest.mark.parametrize("column,values", [
    ("score", ["2.5", "12.75", "0.5"]),
    ("pileup", ["7", "20", "1"]),
    ("pscore", ["3.25", "13.5", "0.75"]),
    ("fc", ["4.5", "8", "1.25"]),
    ("qscore", ["2.5", "12.75", "0.5"]),
])
def test_write_to_bed_score_column(column, values):
    text = _text(_three_peaks().write_to_bed, score_column=column,
                 trackline=False)
    assert [line.split("\t")[4] for line in text.splitlines()] == values


def test_write_to_bed_trackline_escapes_quotes():
    text = _text(_three_peaks().write_to_bed, name=b'my "x"',
                 description=b"about %s")
    assert text.splitlines()[0] == \
        'track name="my \\"x\\" (peaks)" description="about my \\"x\\"" visibility=1'


def test_write_to_bed_prefix_and_description_without_placeholder():
    text = _text(_three_peaks().write_to_bed, name_prefix=b"peak",
                 description=b"plain")
    lines = text.splitlines()
    assert lines[0] == 'track name="MACS (peaks)" description="plain" visibility=1'
    assert [line.split("\t")[3] for line in lines[1:]] == \
        ["peak1", "peak2", "peak3"]


def test_write_to_bed_undecodable_name_falls_back_to_default_trackline():
    text = _text(_three_peaks().write_to_bed, name_prefix=b"p",
                 name=b"\xff", trackline=True)
    assert text == ("track name=MACS description=Unknown\n"
                    "chr1\t100\t200\tp1\t2.5\n"
                    "chr1\t300\t400\tp2\t12.75\n"
                    "chr2\t0\t50\tp3\t0.5\n")


def test_write_to_bed_chromosomes_in_bytes_order():
    p = PeakIO()
    p.add(b"chrX", 1, 2, peak_score=1)
    p.add(b"chr2", 1, 2, peak_score=1)
    p.add(b"chr10", 1, 2, peak_score=1)
    text = _text(p.write_to_bed, trackline=False)
    assert [line.split("\t")[0] for line in text.splitlines()] == \
        ["chr10", "chr2", "chrX"]
    assert [line.split("\t")[3] for line in text.splitlines()] == \
        ["MACS_peak_1", "MACS_peak_2", "MACS_peak_3"]


@pytest.mark.parametrize("score,text", [
    (1234567.0, "1.23457e+06"),
    (2.5e-05, "2.5e-05"),
    (0.0001, "0.0001"),
    (-3.5, "-3.5"),
    (0.0, "0"),
])
def test_write_to_bed_score_format(score, text):
    p = PeakIO()
    p.add(b"chr1", 0, 10, peak_score=score)
    assert _text(p.write_to_bed, trackline=False) == \
        "chr1\t0\t10\tMACS_peak_1\t%s\n" % text


def test_write_to_bed_empty():
    assert _text(PeakIO().write_to_bed, trackline=False) == ""
    assert _text(PeakIO().write_to_bed) == \
        'track name="MACS (peaks)" description="MACS" visibility=1\n'


def test_subpeak_letters_two_summits_then_single():
    p = _subpeaks(2, end=100)
    p.add(b"chr1", 200, 300, summit=250, peak_score=1.0)
    text = _text(p.write_to_bed, trackline=False)
    assert [line.split("\t")[3] for line in text.splitlines()] == \
        ["MACS_peak_1a", "MACS_peak_1b", "MACS_peak_2"]


def test_subpeak_letters_first_26():
    text = _text(_subpeaks(26).write_to_bed, trackline=False)
    assert [line.split("\t")[3] for line in text.splitlines()] == \
        ["MACS_peak_1" + c for c in string.ascii_lowercase]


def test_subpeak_grouping_uses_end_only():
    """Pins the current output.

    Consecutive peaks are grouped as summits of one peak when their ends
    are equal, even if their starts differ. Whether distinct starts should
    split the group is not specified anywhere, so no independent value
    exists.
    """
    p = PeakIO()
    p.add(b"chr1", 0, 100, peak_score=1)
    p.add(b"chr1", 50, 100, peak_score=1)
    text = _text(p.write_to_bed, trackline=False)
    assert text == ("chr1\t0\t100\tMACS_peak_1a\t1\n"
                    "chr1\t50\t100\tMACS_peak_1b\t1\n")


# ------------------------------------
# PeakIO.write_to_summit_bed
# ------------------------------------

def test_write_to_summit_bed_defaults_exact():
    assert _text(_three_peaks().write_to_summit_bed) == (
        "chr1\t150\t151\tMACS_peak_1\t2.5\n"
        "chr1\t310\t311\tMACS_peak_2\t12.75\n"
        "chr2\t0\t1\tMACS_peak_3\t0.5\n")


def test_write_to_summit_bed_trackline_and_column():
    text = _text(_three_peaks().write_to_summit_bed, name=b"exp",
                 score_column="qscore", trackline=True)
    assert text == (
        'track name="exp (summits)" description="exp" visibility=1\n'
        "chr1\t150\t151\texp_peak_1\t2.5\n"
        "chr1\t310\t311\texp_peak_2\t12.75\n"
        "chr2\t0\t1\texp_peak_3\t0.5\n")


def test_write_to_summit_bed_subpeaks():
    text = _text(_subpeaks(3).write_to_summit_bed)
    assert text == ("chr1\t0\t1\tMACS_peak_1a\t1\n"
                    "chr1\t1\t2\tMACS_peak_1b\t1\n"
                    "chr1\t2\t3\tMACS_peak_1c\t1\n")


# ------------------------------------
# PeakIO.tobed / to_summits_bed (write to sys.stdout bound at import)
# ------------------------------------

_STDOUT_SCRIPT = """
from MACS3.IO.PeakIO import PeakIO
p = PeakIO(%s)
p.add(b"chr1", 10, 20, summit=12, peak_score=2.5)
p.add(b"chr1", 10, 20, summit=18, peak_score=1.5)
p.add(b"chr2", 0, 5, summit=3, peak_score=4.0)
r = p.%s()
assert r is None
"""


@pytest.mark.parametrize("ctor,method,expected", [
    ("", "tobed", "chr1\t10\t20\tMACS3_peak_1a\t2.5\n"
                  "chr1\t10\t20\tMACS3_peak_1b\t1.5\n"
                  "chr2\t0\t5\tMACS3_peak_2\t4\n"),
    ("name=b'exp'", "tobed", "chr1\t10\t20\texp_peak_1a\t2.5\n"
                             "chr1\t10\t20\texp_peak_1b\t1.5\n"
                             "chr2\t0\t5\texp_peak_2\t4\n"),
    ("", "to_summits_bed", "chr1\t12\t13\tMACS3_peak_1a\t2.5\n"
                           "chr1\t18\t19\tMACS3_peak_1b\t1.5\n"
                           "chr2\t3\t4\tMACS3_peak_2\t4\n"),
])
def test_stdout_writers(tmp_path, ctor, method, expected):
    res = _run_python(_STDOUT_SCRIPT % (ctor, method), tmp_path)
    assert res.returncode == 0, res.stderr
    assert res.stdout == expected


# ------------------------------------
# PeakIO.write_to_narrowPeak
# ------------------------------------

def test_write_to_narrowPeak_defaults_exact(tmp_path):
    path = tmp_path / "p.narrowPeak"
    with open(path, "w") as fh:
        _three_peaks().write_to_narrowPeak(fh)
    # score column: int(10 * score); int(127.5) = 127; summit offset from start
    assert path.read_text() == (
        "chr1\t100\t200\tMACS_peak_1\t25\t.\t4.5\t3.25\t2.5\t50\n"
        "chr1\t300\t400\tMACS_peak_2\t127\t.\t8\t13.5\t12.75\t10\n"
        "chr2\t0\t50\tMACS_peak_3\t5\t.\t1.25\t0.75\t0.5\t0\n")


def test_write_to_narrowPeak_trackline_name_and_column():
    text = _text(_three_peaks().write_to_narrowPeak, name=b"exp",
                 score_column="pscore", trackline=True)
    assert text == (
        'track type=narrowPeak name="exp" description="exp" nextItemButton=on\n'
        "chr1\t100\t200\texp_peak_1\t32\t.\t4.5\t3.25\t2.5\t50\n"
        "chr1\t300\t400\texp_peak_2\t135\t.\t8\t13.5\t12.75\t10\n"
        "chr2\t0\t50\texp_peak_3\t7\t.\t1.25\t0.75\t0.5\t0\n")


@pytest.mark.parametrize("score,column5", [
    pytest.param(150.5, 1505, id="not_capped_at_1000"),
    pytest.param(2.3, 22, id="float32_truncated"),  # 10*float32(2.3)=22.99999952
    pytest.param(2.5, 25, id="exact"),
    pytest.param(-1.25, -12, id="negative_truncates_toward_zero"),
    pytest.param(0.0, 0, id="zero"),
])
def test_write_to_narrowPeak_score_scaling(score, column5):
    assert int(10 * _f32(score)) == column5
    p = PeakIO()
    p.add(b"chr1", 0, 10, summit=5, peak_score=score)
    line = _text(p.write_to_narrowPeak)
    assert line.split("\t")[4] == str(column5)


def test_write_to_narrowPeak_summit_minus_one():
    p = PeakIO()
    p.add(b"chr1", 100, 200, summit=-1, peak_score=1.0)
    assert _text(p.write_to_narrowPeak) == \
        "chr1\t100\t200\tMACS_peak_1\t10\t.\t0\t0\t0\t-1\n"


def test_write_to_narrowPeak_subpeaks():
    p = PeakIO()
    p.add(b"chr1", 100, 200, summit=120, peak_score=1.0)
    p.add(b"chr1", 100, 200, summit=170, peak_score=2.0)
    assert _text(p.write_to_narrowPeak, name_prefix=b"%s_p", name=b"x") == (
        "chr1\t100\t200\tx_p1a\t10\t.\t0\t0\t0\t20\n"
        "chr1\t100\t200\tx_p1b\t20\t.\t0\t0\t0\t70\n")


def test_write_to_narrowPeak_empty():
    assert _text(PeakIO().write_to_narrowPeak) == ""


# ------------------------------------
# PeakIO.write_to_xls
# ------------------------------------

def test_write_to_xls_exact(tmp_path):
    path = tmp_path / "p.xls"
    with open(path, "w") as fh:
        _three_peaks().write_to_xls(fh)
    # 1-based start and summit, BED end
    assert path.read_text() == XLS_HEADER + (
        "chr1\t101\t200\t100\t151\t7\t3.25\t4.5\t2.5\tMACS_peak_1\n"
        "chr1\t301\t400\t100\t311\t20\t13.5\t8\t12.75\tMACS_peak_2\n"
        "chr2\t1\t50\t50\t1\t1\t0.75\t1.25\t0.5\tMACS_peak_3\n")


@pytest.mark.parametrize("pileup,text", [
    pytest.param(12.345, "12.35", id="float32_above_half"),  # 12.34500027
    pytest.param(2.675, "2.67", id="float32_below_half"),    # 2.67499995
    pytest.param(7.0, "7", id="integer"),
    pytest.param(0.004, "0", id="rounds_to_zero"),
    pytest.param(1234567.0, "1.23457e+06", id="large"),
])
def test_write_to_xls_pileup_rounding(pileup, text):
    p = PeakIO()
    p.add(b"chr1", 0, 10, summit=4, pileup=pileup)
    lines = _text(p.write_to_xls).splitlines()
    assert lines[1].split("\t")[5] == text


def test_write_to_xls_subpeaks_and_prefix():
    p = _subpeaks(2, end=100)
    assert _text(p.write_to_xls, name_prefix=b"%s_", name=b"s") == \
        XLS_HEADER + ("chr1\t1\t100\t100\t1\t0\t0\t0\t0\ts_1a\n"
                      "chr1\t1\t100\t100\t2\t0\t0\t0\t0\ts_1b\n")


def test_write_to_xls_empty():
    assert _text(PeakIO().write_to_xls) == XLS_HEADER


# ------------------------------------
# PeakIO.exclude
# ------------------------------------

EXCLUDE_CASES = [
    pytest.param([(0, 100), (200, 300)], [(50, 60)], id="inner_hit"),
    pytest.param([(0, 100)], [(100, 200)], id="touching_right_kept"),
    pytest.param([(100, 200)], [(0, 100)], id="touching_left_kept"),
    pytest.param([(0, 100), (150, 160), (300, 400)], [(90, 155)],
                 id="one_hits_two"),
    pytest.param([(0, 10), (20, 30)], [(0, 100)], id="all_removed"),
    pytest.param([(300, 400), (0, 100), (150, 200)], [(350, 360)],
                 id="unsorted"),
    pytest.param([(0, 10), (500, 600), (700, 800)], [(5, 6)],
                 id="tail_after_other_exhausted"),
    pytest.param([(0, 10), (40, 50), (90, 100)], [(5, 45), (20, 30)],
                 id="nested_in_other"),
]


@pytest.mark.parametrize("a,b", EXCLUDE_CASES)
def test_exclude_matches_reference(a, b):
    covered = _bases(b)
    expected = sorted(iv for iv in a if not (_bases([iv]) & covered))
    p1, p2 = PeakIO(), PeakIO()
    for s, e in a:
        p1.add(b"chr1", s, e)
    for s, e in b:
        p2.add(b"chr1", s, e)
    p1.exclude(p2)
    assert [(x["start"], x["end"]) for x in p1.peaks[b"chr1"]] == expected
    assert p1.total == len(expected)


def test_exclude_keeps_peak_objects_and_other_chromosomes():
    p1 = PeakIO()
    p1.add(b"chr1", 0, 10, summit=5, peak_score=3.5, name=b"keep")
    p1.add(b"chr1", 20, 30, summit=25)
    p1.add(b"chr2", 0, 10)
    kept = p1.peaks[b"chr1"][0]
    p2 = PeakIO()
    p2.add(b"chr1", 25, 26)
    p2.add(b"chr3", 0, 100)
    p1.exclude(p2)
    assert p1.peaks[b"chr1"] == [kept]
    assert kept["name"] == b"keep" and kept["score"] == 3.5
    assert [(x["start"], x["end"]) for x in p1.peaks[b"chr2"]] == [(0, 10)]
    assert p1.total == 2
    assert p1.get_chr_names() == {b"chr1", b"chr2"}


def test_exclude_with_empty_other_keeps_all():
    p1 = _three_peaks()
    before = _text(p1.write_to_bed)
    p1.exclude(PeakIO())
    assert _text(p1.write_to_bed) == before
    assert p1.total == 3


def test_exclude_rejects_non_peakio():
    with pytest.raises(AssertionError):
        _three_peaks().exclude(RegionIO())


# ------------------------------------
# PeakIO.read_from_xls
# ------------------------------------


# ------------------------------------
# RegionIO
# ------------------------------------

def _regionio(spec):
    r = RegionIO()
    for chrom, intervals in spec.items():
        for s, e in intervals:
            r.add_loc(chrom, s, e)
    return r


def _runs(bases):
    out = []
    for b in sorted(bases):
        if out and out[-1][1] == b:
            out[-1][1] = b + 1
        else:
            out.append([b, b + 1])
    return [tuple(x) for x in out]


def _bed(spec):
    return "".join("%s\t%d\t%d\n" % (c.decode(), s, e)
                   for c in sorted(spec) for s, e in spec[c])


def test_RegionIO_add_loc_write_to_bed_insertion_order():
    r = _regionio({b"chr2": [(5, 6)], b"chr1": [(30, 40), (0, 10)]})
    assert _text(r.write_to_bed) == "chr1\t30\t40\nchr1\t0\t10\nchr2\t5\t6\n"


def test_RegionIO_sort():
    r = _regionio({b"chr2": [(5, 6), (1, 2)], b"chr1": [(30, 40), (0, 50),
                                                        (0, 10)]})
    r.sort()
    assert _text(r.write_to_bed) == \
        "chr1\t0\t10\nchr1\t0\t50\nchr1\t30\t40\nchr2\t1\t2\nchr2\t5\t6\n"


@pytest.mark.parametrize("spec,expected", [
    ({}, set()),
    ({b"chr1": [(0, 1)]}, {b"chr1"}),
    ({b"chr2": [(0, 1)], b"chr1": [(0, 1)]}, {b"chr1", b"chr2"}),
], ids=["empty", "one", "many"])
def test_RegionIO_get_chr_names(spec, expected):
    assert _regionio(spec).get_chr_names() == expected


@pytest.mark.parametrize("intervals", [
    pytest.param([(0, 10), (20, 30)], id="disjoint"),
    pytest.param([(0, 10), (5, 20)], id="overlapping"),
    pytest.param([(0, 10), (10, 20)], id="touching"),
    pytest.param([(50, 60), (0, 10), (5, 20)], id="unsorted"),
    pytest.param([(0, 5), (4, 9), (8, 12)], id="chain"),
    pytest.param([(3, 7)], id="single"),
])
def test_RegionIO_merge_overlap(intervals):
    r = _regionio({b"chr1": intervals, b"chr2": [(7, 8)]})
    r.merge_overlap()
    assert _text(r.write_to_bed) == \
        _bed({b"chr1": _runs(_bases(intervals)), b"chr2": [(7, 8)]})


def test_RegionIO_merge_overlap_empty():
    r = RegionIO()
    r.merge_overlap()
    assert _text(r.write_to_bed) == ""
    assert r.get_chr_names() == set()


def test_RegionIO_write_to_bed_file(tmp_path):
    r = _regionio({b"chr1": [(0, 2147483647)]})
    path = tmp_path / "r.bed"
    with open(path, "w") as fh:
        r.write_to_bed(fh)
    assert path.read_text() == "chr1\t0\t2147483647\n"


def test_RegionIO_add_loc_errors():
    with pytest.raises(OverflowError):
        RegionIO().add_loc(b"chr1", 0, 2147483648)
    with pytest.raises(TypeError):
        RegionIO().add_loc("chr1", 0, 10)


# ------------------------------------
# BroadPeakContent
# ------------------------------------

def test_BroadPeakContent_fields():
    bp = BroadPeakContent(100, 500, 2.3, b"150", b"300", 3, b"1,50,1",
                          b"0,200,399", 7.5, 3.1, 4.2, 1.9)
    assert (bp["start"], bp["end"], bp["length"]) == (100, 500, 400)
    assert bp["score"] == _f32(2.3)
    assert (bp["thickStart"], bp["thickEnd"]) == (b"150", b"300")
    assert bp["blockNum"] == 3
    assert (bp["blockSizes"], bp["blockStarts"]) == (b"1,50,1", b"0,200,399")
    assert bp["pileup"] == _f32(7.5)
    assert bp["pscore"] == _f32(3.1)
    assert bp["fc"] == _f32(4.2)
    assert bp["qscore"] == _f32(1.9)
    assert bp["name"] == b"MACS3"
    assert bp["chrom"] is None
    assert str(bp) == "start:100;end:500;score:2.300000"


# ------------------------------------
# BroadPeakIO
# ------------------------------------

def _broad_peaks():
    bp = BroadPeakIO()
    bp.add(b"chr1", 100, 500, score=2.5, thickStart=b"150", thickEnd=b"300",
           blockNum=3, blockSizes=b"1,50,1", blockStarts=b"0,200,399",
           pileup=7.0, pscore=3.25, fold_change=4.5, qscore=2.5)
    bp.add(b"chr1", 1000, 1200, score=1.5, pileup=2.0, pscore=1.75,
           fold_change=2.25, qscore=1.5)
    bp.add(b"chr2", 0, 300, score=150.5, thickStart=b"10", thickEnd=b"20",
           blockNum=2, blockSizes=b"1,1", blockStarts=b"0,299", pileup=30.0,
           pscore=160.25, fold_change=9.5, qscore=150.5)
    return bp


def test_BroadPeakIO_add_defaults_and_total():
    bp = BroadPeakIO()
    assert bp.peaks == {} and bp.total() == 0
    bp.add(b"chr1", 10, 20)
    x = bp.peaks[b"chr1"][0]
    assert (x["start"], x["end"], x["length"], x["score"]) == (10, 20, 10, 0)
    assert (x["thickStart"], x["thickEnd"], x["blockNum"]) == (b".", b".", 0)
    assert (x["blockSizes"], x["blockStarts"]) == (b".", b".")
    assert (x["pileup"], x["pscore"], x["fc"], x["qscore"]) == (0, 0, 0, 0)
    assert x["name"] == b"NA"
    assert _broad_peaks().total() == 3


def test_BroadPeakIO_add_rejects_str_fields():
    with pytest.raises(TypeError):
        BroadPeakIO().add(b"chr1", 0, 10, thickStart="1")


@pytest.mark.parametrize("method,kwarg,field", [
    ("filter_pscore", "pscore", "pscore"),
    ("filter_qscore", "qscore", "qscore"),
    ("filter_fc", "fold_change", "fc"),
])
def test_BroadPeakIO_filter_lower_bound_inclusive(method, kwarg, field):
    bp = BroadPeakIO()
    for i, v in enumerate([4.0, 1.0, 2.5]):
        bp.add(b"chr1", 100 * i, 100 * i + 50, **{kwarg: v})
    getattr(bp, method)(2.5)
    assert [x[field] for x in bp.peaks[b"chr1"]] == [4.0, 2.5]
    assert bp.total() == 2


@pytest.mark.parametrize("method,kwarg,kept", [
    ("filter_pscore", "pscore", True),
    ("filter_qscore", "qscore", True),
    ("filter_fc", "fold_change", False),
])
def test_BroadPeakIO_filter_cut_precision(method, kwarg, kept):
    # The stored value is float32(2.3) < 2.3. The pscore and qscore cuts are
    # declared cython.float, so the cut is float32(2.3) and the peak is kept;
    # the fc cut is annotated ``float``, which Cython maps to a C double, so
    # the cut stays 2.3 and the peak is below it.
    bp = BroadPeakIO()
    bp.add(b"chr1", 0, 10, **{kwarg: 2.3})
    getattr(bp, method)(2.3)
    assert bp.total() == (1 if kept else 0)


@pytest.mark.parametrize("low,up,expected", [
    pytest.param(2.0, 4.0, [2.0, 3.0], id="upper_exclusive"),
    pytest.param(2.0, -1, [2.0, 3.0, 4.0], id="negative_upper_ignored"),
    pytest.param(2.0, 0, [], id="zero_upper_applies"),
    pytest.param(3.0, 2.0, [], id="upper_below_lower"),
])
def test_BroadPeakIO_filter_fc_range(low, up, expected):
    bp = BroadPeakIO()
    for i, v in enumerate([1.0, 2.0, 3.0, 4.0]):
        bp.add(b"chr1", 100 * i, 100 * i + 50, fold_change=v)
    bp.filter_fc(low, up)
    assert [x["fc"] for x in bp.peaks[b"chr1"]] == expected
    assert bp.total() == len(expected)


def test_BroadPeakIO_filter_fc_default_upper():
    bp = BroadPeakIO()
    for i, v in enumerate([1.0, 2.0, 3.0]):
        bp.add(b"chr1", 100 * i, 100 * i + 50, fold_change=v)
    bp.filter_fc(2.0)
    assert [x["fc"] for x in bp.peaks[b"chr1"]] == [2.0, 3.0]


def test_BroadPeakIO_write_to_gappedPeak_exact(tmp_path):
    path = tmp_path / "b.gappedPeak"
    with open(path, "w") as fh:
        _broad_peaks().write_to_gappedPeak(fh)
    # the peak without blocks (thickStart ".") is skipped but still numbered
    assert path.read_text() == (
        'track name="peak" description="peak" type=gappedPeak nextItemButton=on\n'
        "chr1\t100\t500\tpeak_1\t25\t.\t0\t0\t0\t3\t1,50,1\t0,200,399\t4.5\t3.25\t2.5\n"
        "chr2\t0\t300\tpeak_3\t1505\t.\t0\t0\t0\t2\t1,1\t0,299\t9.5\t160.25\t150.5\n")


def test_BroadPeakIO_write_to_gappedPeak_options():
    text = _text(_broad_peaks().write_to_gappedPeak, name_prefix=b"%s_b",
                 name=b"exp", description=b"d %s", score_column="pscore",
                 trackline=True)
    assert text == (
        'track name="exp" description="d exp" type=gappedPeak nextItemButton=on\n'
        "chr1\t100\t500\texp_b1\t32\t.\t0\t0\t0\t3\t1,50,1\t0,200,399\t4.5\t3.25\t2.5\n"
        "chr2\t0\t300\texp_b3\t1602\t.\t0\t0\t0\t2\t1,1\t0,299\t9.5\t160.25\t150.5\n")
    assert _text(_broad_peaks().write_to_gappedPeak,
                 trackline=False).startswith("chr1\t100\t500\tpeak_1\t")


def test_BroadPeakIO_write_to_Bed12_exact():
    assert _text(_broad_peaks().write_to_Bed12) == (
        'track name="peak" description="peak" type=bed nextItemButton=on\n'
        "chr1\t100\t500\tpeak_1\t25\t.\t150\t300\t0\t3\t1,50,1\t0,200,399\n"
        "chr1\t1000\t1200\tpeak_2\t15\t.\n"
        "chr2\t0\t300\tpeak_3\t1505\t.\t10\t20\t0\t2\t1,1\t0,299\n")


def test_BroadPeakIO_write_to_Bed12_no_trackline_prefix():
    assert _text(_broad_peaks().write_to_Bed12, name_prefix=b"%s_",
                 name=b"e", trackline=False) == (
        "chr1\t100\t500\te_1\t25\t.\t150\t300\t0\t3\t1,50,1\t0,200,399\n"
        "chr1\t1000\t1200\te_2\t15\t.\n"
        "chr2\t0\t300\te_3\t1505\t.\t10\t20\t0\t2\t1,1\t0,299\n")


def test_BroadPeakIO_write_to_broadPeak_exact():
    assert _text(_broad_peaks().write_to_broadPeak) == (
        'track type=broadPeak name="peak" description="peak" nextItemButton=on\n'
        "chr1\t100\t500\tpeak_1\t25\t.\t4.5\t3.25\t2.5\n"
        "chr1\t1000\t1200\tpeak_2\t15\t.\t2.25\t1.75\t1.5\n"
        "chr2\t0\t300\tpeak_3\t1505\t.\t9.5\t160.25\t150.5\n")


def test_BroadPeakIO_write_to_broadPeak_options():
    assert _text(_broad_peaks().write_to_broadPeak, name_prefix=b"%s_broad_",
                 name=b"exp", score_column="qscore", trackline=False) == (
        "chr1\t100\t500\texp_broad_1\t25\t.\t4.5\t3.25\t2.5\n"
        "chr1\t1000\t1200\texp_broad_2\t15\t.\t2.25\t1.75\t1.5\n"
        "chr2\t0\t300\texp_broad_3\t1505\t.\t9.5\t160.25\t150.5\n")


def test_BroadPeakIO_write_to_xls_exact(tmp_path):
    path = tmp_path / "b.xls"
    with open(path, "w") as fh:
        _broad_peaks().write_to_xls(fh)
    assert path.read_text() == BROAD_XLS_HEADER + (
        "chr1\t101\t500\t400\t7\t3.25\t4.5\t2.5\tMACS_peak_1\n"
        "chr1\t1001\t1200\t200\t2\t1.75\t2.25\t1.5\tMACS_peak_2\n"
        "chr2\t1\t300\t300\t30\t160.25\t9.5\t150.5\tMACS_peak_3\n")


def test_BroadPeakIO_writers_empty():
    bp = BroadPeakIO()
    assert _text(bp.write_to_gappedPeak, trackline=False) == ""
    assert _text(bp.write_to_Bed12, trackline=False) == ""
    assert _text(bp.write_to_broadPeak, trackline=False) == ""
    assert _text(bp.write_to_xls) == BROAD_XLS_HEADER
