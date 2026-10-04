import sys
import tempfile
import types
import pytest
import numpy as np

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.CallPeakUnit import CallerFromAlignments
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.IO.PeakIO import PeakIO


# ---------------------------------------------------------------------------
# Provide tiny Cython stubs so the Python sources import cleanly.
# ---------------------------------------------------------------------------
cython_stub = sys.modules.get("cython")
if cython_stub is None:
    cython_stub = types.ModuleType("cython")
    sys.modules["cython"] = cython_stub


def _identity_decorator(*args, **kwargs):
    def decorate(func):
        return func
    if args and callable(args[0]) and len(args) == 1 and not kwargs:
        return args[0]
    return decorate


for name in ("cfunc", "ccall", "cclass", "locals", "inline", "returns",
             "boundscheck", "wraparound"):
    setattr(cython_stub, name, getattr(cython_stub, name, _identity_decorator))

for attr, default in [
    ("declare", lambda *a, **k: None),
    ("short", int),
    ("float", float),
    ("double", float),
    ("int", int),
    ("long", int),
    ("ulong", int),
    ("bint", bool),
]:
    setattr(cython_stub, attr, getattr(cython_stub, attr, default))

cimports_mod = sys.modules.get("cython.cimports")
if cimports_mod is None:
    cimports_mod = types.ModuleType("cython.cimports")
    sys.modules["cython.cimports"] = cimports_mod
cython_stub.cimports = cimports_mod

if "cython.cimports.cpython" not in sys.modules:
    cpython_mod = types.ModuleType("cython.cimports.cpython")
    cpython_mod.bool = bool
    sys.modules["cython.cimports.cpython"] = cpython_mod
    cimports_mod.cpython = cpython_mod

if "cython.cimports.numpy" not in sys.modules:
    numpy_mod = types.ModuleType("cython.cimports.numpy")
    numpy_mod.ndarray = lambda *args, **kwargs: None
    sys.modules["cython.cimports.numpy"] = numpy_mod
    cimports_mod.numpy = numpy_mod

def make_fwtrack(layout, fw=50):
    track = FWTrack(fw=fw)
    for chrom, positions in layout.items():
        for pos in positions:
            track.add_loc(chrom, pos, strand=0)
    track.finalize()
    return track


def make_tracks():
    treat = make_fwtrack({b"chr1": [10, 30]})
    ctrl = make_fwtrack({b"chr1": [15], b"chr2": [5]})
    return treat, ctrl


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


def make_petrack_from_profile(signal, start):
    """Encode an integer coverage profile as a paired-end track."""
    track = PETrackI(buffer_size=200000)
    for level in range(1, int(signal.max()) + 1):
        edges = np.diff(np.r_[False, signal >= level, False].astype(int))
        for left, right in zip(np.flatnonzero(edges == 1),
                               np.flatnonzero(edges == -1)):
            track.add_loc(b"chrSynthetic", start + int(left),
                          start + int(right))
    track.finalize()
    return track


def test_constructor_rejects_unknown_track_type():
    with pytest.raises(Exception):
        CallerFromAlignments(object(), object())


def test_destroy_accepts_empty_state(tmp_path):
    calc = CallerFromAlignments(*make_tracks())
    if hasattr(calc, "pileup_data_files"):
        dummy = tmp_path / "pile.tmp"
        dummy.write_text("tmp")
        calc.pileup_data_files = {b"chr1": str(dummy)}
    calc.destroy()  # should not raise even if files missing or attribute hidden


def test_set_pseudocount_no_exception():
    calc = CallerFromAlignments(*make_tracks())
    calc.set_pseudocount(5.0)
    if hasattr(calc, "pseudocount"):
        assert calc.pseudocount == pytest.approx(5.0)


def test_enable_trackline_idempotent():
    calc = CallerFromAlignments(*make_tracks())
    calc.enable_trackline()
    calc.enable_trackline()
    if hasattr(calc, "trackline"):
        assert calc.trackline is True


def test_call_peaks_with_no_chromosomes_returns_peakio():
    treat, ctrl = make_tracks()
    calc = CallerFromAlignments(treat, ctrl)
    if hasattr(calc, "chromosomes"):
        calc.chromosomes = []
    if not hasattr(calc, "pqtable"):
        pytest.skip("pqtable not exposed in current build")
    calc.pqtable[0.0] = 0.0
    peaks = calc.call_peaks(['p'], [1.0], min_length=10, max_gap=5, call_summits=False, cutoff_analysis=False)
    assert isinstance(peaks, PeakIO)


def test_call_summits_rejects_maximum_in_below_cutoff_gap(tmp_path,
                                                           monkeypatch):
    """A smoothed summit must map to an above-cutoff signal chunk."""
    signal = np.rint(np.interp(np.arange(123),
                               np.linspace(0, 122, 8),
                               [27, 1, 27, 8, 14, 10, 27, 12])).astype(int)
    track = PETrackI(buffer_size=10000)

    for level in range(1, int(signal.max()) + 1):
        edges = np.diff(np.r_[False, signal >= level, False].astype(int))
        starts = np.flatnonzero(edges == 1)
        ends = np.flatnonzero(edges == -1)
        for start, end in zip(starts, ends):
            track.add_loc(b"chrSynthetic", 1000 + int(start),
                          1000 + int(end))
    track.finalize()

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    caller = CallerFromAlignments(track, None, ctrl_d_s=[],
                                  ctrl_scaling_factor_s=[], lambda_bg=1.0)
    try:
        peaks = caller.call_peaks(["f"], [9.0], min_length=75,
                                  max_gap=100, call_summits=True)
    finally:
        caller.destroy()

    peak = peaks.peaks[b"chrSynthetic"][0]
    actual_pileup = signal[peak["summit"] - 1000]
    assert peak["summit"] == 1035
    assert peak["pileup"] == actual_pileup == 27
    assert peak["fc"] == pytest.approx((actual_pileup + 1) / 2)


@pytest.mark.parametrize("start", [0, 1, 5, 9, 10, 1000])
def test_call_summits_keeps_right_edge_candidate(tmp_path, monkeypatch,
                                                  start):
    """Regression test for the padding-coordinate mismatch in issue #747."""
    signal = make_two_summit_profile()
    track = make_petrack_from_profile(signal, start)

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    caller = CallerFromAlignments(track, None, ctrl_d_s=[],
                                  ctrl_scaling_factor_s=[], lambda_bg=1.0)
    try:
        peaks = caller.call_peaks(["f"], [2.0], min_length=50,
                                  max_gap=100, call_summits=True)
    finally:
        caller.destroy()

    rows = peaks.peaks[b"chrSynthetic"]
    relative_summits = [row["summit"] - start for row in rows]
    assert len(relative_summits) == 2
    assert abs(relative_summits[0] - 60) <= 1
    assert abs(relative_summits[1] - 236) <= 1
    assert all(row["start"] == start for row in rows)
    assert all(row["end"] == start + len(signal) for row in rows)
    assert all(row["start"] <= row["summit"] < row["end"] for row in rows)


# ===========================================================================
# Comprehensive tests for CallerFromAlignments (added below the original
# tests above, which are kept unchanged).
#
# Expected values are derived independently of CallPeakUnit:
#
# * pileups: each read or fragment becomes a half-open interval following
#   the extension rule documented at each ``*_intervals`` helper, and the
#   coverage is counted directly;
# * a chromosome is cut at every position where the treatment coverage or
#   the coverage of any control window changes, from 0 to the end of the
#   treatment pileup (MACS3 keeps these cuts even where the scaled value
#   does not change, which decides the summit of a flat-topped peak);
# * treatment value = coverage * treat_scaling_factor; control lambda =
#   max(lambda_bg, max over windows of coverage * window scaling factor),
#   both in float32;
# * p-score = -log10(scipy.stats.poisson.sf(int(t), lambda)) =
#   -log10 P(X > t);
# * q-scores follow the Benjamini-Hochberg variant MACS3 implements
#   (``ref_qtable``);
# * peaks follow the cutoff / max_gap / min_length rules (``RefCaller``).
#
# C-only (cfunc) helpers and the Python-visible callers covering them:
#   get_pscore, __cal_pscore, __cal_pvalue_qvalue_table, __cal_qscore
#       -> call_peaks / call_broadpeaks with 'p' and 'q' cutoffs
#   __cal_FE, __cal_subtraction -> call_peaks with 'f' and 's' cutoffs
#   apply_multiple_cutoffs -> test_call_peaks_multiple_cutoffs_*
#   pileup_treat_ctrl_a_chromosome, __chrom_pair_treat_ctrl,
#   __write_bedGraph_for_a_chromosome, clean_up_ndarray
#       -> bedGraph tests and repeated call_peaks
#   __chrom_call_peak_using_certain_criteria, __close_peak_wo_subpeaks,
#   __close_peak_with_subpeaks -> call_peaks(call_summits=False/True)
#   __pre_computes -> call_peaks/call_broadpeaks(cutoff_analysis=True)
#   __chrom_call_broadpeak_using_certain_criteria,
#   __close_peak_for_broad_region, __add_broadpeak,
#   mean_from_value_length, getitem_then_subtract -> call_broadpeaks
# C-only helpers that nothing in the module calls (unreachable, so not
# tested): get_logLR_asym, chi2_k1_cdf, log10_chi2_k1_cdf, chi2_k2_cdf,
# log10_chi2_k2_cdf, chi2_k4_cdf, log10_chi2_k4_CDF,
# get_from_multiple_scores, get_logFE, get_subtraction, left_sum,
# right_sum, left_forward, right_forward, median_from_value_length,
# find_optimal_cutoff (its only call is commented out), __cal_logLR,
# __cal_logFE, __cal_score.
# ===========================================================================

import functools
import logging
import re
import tempfile

import numpy as np
from scipy.stats import poisson

from MACS3.Signal.PairedEndTrack import (PETrackI,
                                         PETrackII)
from MACS3.IO.PeakIO import BroadPeakIO

F32 = np.float32
CUTOFF_LIST = [round(x, 5)
               for x in sorted(list(np.arange(0.3, 10.0, 0.3)), reverse=True)]


# ------------------------------------
# reference model
# ------------------------------------

def _clip0(iv):
    """MACS3 clips negative pileup coordinates of FWTrack/PETrackI to 0."""
    s, e, w = iv
    return (max(s, 0), max(e, 0), w)


def se_treat_intervals(plus, minus, d, end_shift=0):
    """SE treatment fragment: a + tag at p covers [p + s, p + s + d), a - tag
    at p covers [p - s - d, p - s), where s is end_shift."""
    ivs = [(p + end_shift, p + end_shift + d, 1) for p in plus]
    ivs += [(p - end_shift - d, p - end_shift, 1) for p in minus]
    return [_clip0(x) for x in ivs]


def se_ctrl_intervals(plus, minus, w):
    """SE control window of size w around each tag: a + tag at p covers
    [p - w//2, p + w - w//2), a - tag covers [p - (w - w//2), p + w//2)."""
    ivs = [(p - w // 2, p + w - w // 2, 1) for p in plus]
    ivs += [(p - (w - w // 2), p + w // 2, 1) for p in minus]
    return [_clip0(x) for x in ivs]


def pe_treat_intervals(frags):
    """A fragment (l, r[, count]) covers [l, r) with weight count."""
    return [(f[0], f[1], f[2] if len(f) > 2 else 1) for f in frags]


def pe1_ctrl_intervals(frags, w):
    """PETrackI control: each fragment end x covers [x - w//2, x + w//2)."""
    ivs = []
    for f in frags:
        for x in (f[0], f[1]):
            ivs.append((x - w // 2, x + w // 2, 1))
    return [_clip0(x) for x in ivs]


def pe2_ctrl_intervals(frags, w):
    """PETrackII control: each fragment end x covers [x - w//2, x - w//2 + w)
    with the fragment count as weight."""
    ivs = []
    for l, r, c in frags:
        for x in (l, r):
            ivs.append((x - w // 2, x - w // 2 + w, c))
    return [_clip0(x) for x in ivs]


def _breakpoints(ivs):
    """Positions where the weighted coverage of ``ivs`` changes."""
    net = {}
    for s, e, w in ivs:
        net[s] = net.get(s, 0) + w
        net[e] = net.get(e, 0) - w
    return {p for p, v in net.items() if v != 0}


def _coverage(ivs, pos):
    """Weighted coverage of ``ivs`` at each position of the array ``pos``."""
    cov = np.zeros(pos.shape[0], dtype=np.int64)
    for s, e, w in ivs:
        cov += w * ((pos >= s) & (pos < e))
    return cov


@functools.lru_cache(maxsize=None)
def ref_pscore(t, lam):
    """-log10 P(X > int(t)) for X ~ Poisson(lam), rounded to 5 decimals
    (as MACS3 rounds it, so that lambdas differing in the last float32 bit
    share one value) and stored as float32."""
    return float(F32(round(-np.log10(poisson.sf(int(t), float(lam))), 5)))


def ref_qtable(stat):
    """MACS3's q-score table from {pscore: total length}.

    P-scores are visited in decreasing order; k is 1 plus the total length
    of all strictly larger p-scores and N the total length.  The q-score is
    p + log10(k) - log10(N), never larger than the q-score above it, and 0
    from the first value <= 0 on.
    """
    n = sum(stat.values())
    table = {}
    k = 1
    pre_q = float("inf")
    zero = False
    for v in sorted(stat, reverse=True):
        if zero:
            table[v] = 0.0
            continue
        q = min(v + np.log10(k) - np.log10(n), pre_q)
        if q <= 0:
            table[v] = 0.0
            zero = True
            continue
        table[v] = q
        pre_q = q
        k += stat[v]
    return table


def _rle_lines(chrom, starts, ends, values):
    """bedGraph lines merging consecutive equal values ("%.5f")."""
    lines = []
    cur_s, cur_v = starts[0], values[0]
    for i in range(1, len(values)):
        if values[i] != cur_v:
            lines.append("%s\t%d\t%d\t%.5f" % (chrom.decode(), cur_s,
                                               ends[i - 1], cur_v))
            cur_s, cur_v = starts[i], values[i]
    lines.append("%s\t%d\t%d\t%.5f" % (chrom.decode(), cur_s, ends[-1],
                                       cur_v))
    return lines


class RefCaller:
    """Independent model of CallerFromAlignments.

    ``treat`` maps chromosome -> treatment intervals, ``windows`` maps
    chromosome -> [(control intervals, scaling factor), ...].  Only the
    chromosomes given are modelled (the common ones).
    """

    def __init__(self, treat, windows, treat_scale=1.0, lambda_bg=0.0,
                 pseudocount=1.0, no_lambda=False):
        self.pc = pseudocount
        self.segs = {}
        for chrom in sorted(treat):
            tiv = treat[chrom]
            bps = _breakpoints(tiv)
            end = max(bps)         # where the treatment coverage returns to 0
            wins = [] if no_lambda else windows[chrom]
            for iv, _ in wins:
                bps |= _breakpoints(iv)
            ends = sorted(p for p in bps | {end} if 0 < p <= end)
            starts = [0] + ends[:-1]
            pos = np.array(starts)
            t = np.maximum(_coverage(tiv, pos).astype(F32) * F32(treat_scale),
                           F32(0))
            if no_lambda:
                c = np.full(len(ends), F32(lambda_bg), dtype=F32)
            else:
                c = None
                for iv, scale in wins:
                    cw = np.maximum(_coverage(iv, pos).astype(F32) * F32(scale),
                                    F32(lambda_bg))
                    c = cw if c is None else np.maximum(c, cw)
            self.segs[chrom] = list(zip(starts, ends, t.tolist(), c.tolist()))
        self.stat = {}
        for chrom, segs in self.segs.items():
            for s, e, t, c in segs:
                p = ref_pscore(t, c)
                self.stat[p] = self.stat.get(p, 0) + (e - s)
        self.qtable = ref_qtable(self.stat)

    # scores ------------------------------------------------------------
    def score(self, sym, t, c):
        if sym == "p":
            return ref_pscore(t, c)
        if sym == "q":
            return self.qtable[ref_pscore(t, c)]
        if sym == "f":
            return float(F32((t + self.pc) / (c + self.pc)))
        if sym == "s":
            return float(F32(t) - F32(c))
        raise ValueError(sym)

    def regions(self, chrom, syms, cutoffs, max_gap, min_length):
        """Groups of segment indices above any cutoff, merged across gaps
        <= max_gap and kept when at least min_length long."""
        segs = self.segs[chrom]
        groups = []
        for i, (s, e, t, c) in enumerate(segs):
            if not any(self.score(y, t, c) > k for y, k in zip(syms, cutoffs)):
                continue
            if groups and s - segs[groups[-1][-1]][1] <= max_gap:
                groups[-1].append(i)
            else:
                groups.append([i])
        return [g for g in groups
                if segs[g[-1]][1] - segs[g[0]][0] >= min_length]

    def peaks(self, syms, cutoffs, min_length=200, max_gap=50):
        """Narrow peaks without sub-peaks: the summit is the midpoint of the
        middle one of the segments holding the highest treatment value."""
        out = []
        for chrom, segs in self.segs.items():
            for g in self.regions(chrom, syms, cutoffs, max_gap, min_length):
                tmax = max(segs[i][2] for i in g)
                tied = [i for i in g if segs[i][2] == tmax]
                si = tied[(len(tied) + 1) // 2 - 1]
                s, e, t, c = segs[si]
                if any(self.score(y, t, c) < k for y, k in zip(syms, cutoffs)):
                    continue
                p = ref_pscore(t, c)
                out.append(dict(chrom=chrom, start=segs[g[0]][0],
                                end=segs[g[-1]][1], summit=(s + e) // 2,
                                pileup=t, pscore=p, qscore=self.qtable[p],
                                fc=(t + self.pc) / (c + self.pc)))
        return out

    def broad_regions(self, chrom, syms, cutoffs, max_gap, min_length):
        """Regions with length-weighted mean scores of their segments."""
        segs = self.segs[chrom]
        out = []
        for g in self.regions(chrom, syms, cutoffs, max_gap, min_length):
            ln = np.array([segs[i][1] - segs[i][0] for i in g], dtype=float)

            def wmean(vals):
                return float(np.sum(np.array(vals) * ln) / ln.sum())
            ts = [segs[i][2] for i in g]
            cs = [segs[i][3] for i in g]
            ps = [ref_pscore(t, c) for t, c in zip(ts, cs)]
            out.append(dict(start=segs[g[0]][0], end=segs[g[-1]][1],
                            pileup=wmean(ts), pscore=wmean(ps),
                            qscore=wmean([self.qtable[p] for p in ps]),
                            fc=wmean([(t + self.pc) / (c + self.pc)
                                      for t, c in zip(ts, cs)])))
        return out

    def broadpeaks(self, syms, lvl1, lvl2, min_length=200, gap1=50,
                   gap2=400):
        """Broad regions (level 2) with their level-1 blocks; a region
        without a level-1 block gets two 1-bp blocks at its ends, a region
        whose blocks do not reach its ends gets 1-bp end blocks."""
        out = []
        for chrom in sorted(self.segs):
            strong = self.broad_regions(chrom, syms, lvl1, gap1, min_length)
            if not strong:
                continue
            for weak in self.broad_regions(chrom, syms, lvl2, gap2,
                                           min_length):
                inside = [x for x in strong if weak["start"] <= x["start"]
                          and x["end"] <= weak["end"]]
                s, e = weak["start"], weak["end"]
                if not inside:
                    blocks = [(s, s + 1), (e - 1, e)]
                else:
                    blocks = [(x["start"], x["end"]) for x in inside]
                    if blocks[0][0] != s:
                        blocks.insert(0, (s, s + 1))
                    if blocks[-1][1] != e:
                        blocks.append((e - 1, e))
                # thick part always spans the whole region (end blocks are
                # added when the level-1 blocks do not reach the ends)
                out.append(dict(
                    chrom=chrom, start=s, end=e, score=weak["qscore"],
                    pileup=weak["pileup"], pscore=weak["pscore"],
                    fc=weak["fc"], qscore=weak["qscore"],
                    thickStart=b"%d" % s, thickEnd=b"%d" % e,
                    blockNum=len(blocks),
                    blockSizes=b",".join(b"%d" % (y - x) for x, y in blocks),
                    blockStarts=b",".join(b"%d" % (x - s) for x, y in blocks)))
        return out

    def bedgraph(self, which, denominator=1.0):
        lines = []
        for chrom, segs in self.segs.items():
            col = 2 if which == "treat" else 3
            vals = [float(F32(x[col]) / F32(denominator)) for x in segs]
            lines += _rle_lines(chrom, [x[0] for x in segs],
                                [x[1] for x in segs], vals)
        return lines

    def cutoff_rows(self, max_gap, min_length):
        """Expected rows of the cutoff-analysis file."""
        stat = dict(self.stat)
        for k in CUTOFF_LIST:
            stat.setdefault(k, 0)
        qtab = ref_qtable(stat)
        rows = []
        for k in CUTOFF_LIST:
            n = 0
            total = 0
            for chrom, segs in self.segs.items():
                for g in self.regions(chrom, ["p"], [k], max_gap, min_length):
                    n += 1
                    total += segs[g[-1]][1] - segs[g[0]][0]
            if n > 0:
                rows.append((k, qtab[k], n, total, total / n))
        return rows


# ------------------------------------
# track builders and comparison helpers
# ------------------------------------

def build_fw(plus, minus, fw=50):
    """FWTrack from {chrom: [positions]} for each strand."""
    track = FWTrack(fw=fw)
    for chrom in sorted(set(plus) | set(minus)):
        for p in plus.get(chrom, []):
            track.add_loc(chrom, p, 0)
        for p in minus.get(chrom, []):
            track.add_loc(chrom, p, 1)
    track.finalize()
    return track


def build_pe1(frags):
    track = PETrackI()
    for chrom in sorted(frags):
        for l, r in frags[chrom]:
            track.add_loc(chrom, l, r)
    track.finalize()
    return track


def build_pe2(frags):
    track = PETrackII()
    for chrom in sorted(frags):
        for i, (l, r, c) in enumerate(frags[chrom]):
            track.add_loc(chrom, l, r, b"BC%d" % (i % 3), c)
    track.finalize()
    return track


def peak_rows(peakio):
    rows = []
    for chrom in sorted(peakio.get_chr_names()):
        for p in peakio.get_data_from_chrom(chrom):
            rows.append(dict(chrom=chrom, start=p["start"], end=p["end"],
                             summit=p["summit"], pileup=p["pileup"],
                             pscore=p["pscore"], qscore=p["qscore"],
                             fc=p["fc"], score=p["score"]))
    return rows


def assert_peaks_match(got, expected):
    """Coordinates exactly; pileup and fold change to float32 precision;
    p/q-scores to 1e-4 (MACS3 rounds -log10 p to 5 decimals and stops its
    Poisson series at a 1e-5 relative change)."""
    assert [(p["chrom"], p["start"], p["end"], p["summit"]) for p in got] == \
        [(p["chrom"], p["start"], p["end"], p["summit"]) for p in expected]
    for g, e in zip(got, expected):
        assert g["pileup"] == pytest.approx(e["pileup"], rel=1e-6)
        assert g["fc"] == pytest.approx(e["fc"], rel=1e-6)
        assert g["pscore"] == pytest.approx(e["pscore"], abs=1e-4)
        assert g["qscore"] == pytest.approx(e["qscore"], abs=1e-4)
        assert g["score"] == pytest.approx(e["qscore"], abs=1e-4)


def read_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


@pytest.fixture
def caller_tmp(tmp_path, monkeypatch):
    """Send CallerFromAlignments' per-chromosome pickles (mkstemp) into a
    directory of this test instead of the system temporary directory."""
    d = tmp_path / "pileup_tmp"
    d.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(d))
    return d


# ------------------------------------
# single-end scenario: chr1 and chr2 in both samples, chr3 only in the
# treatment and chr4 only in the control (both ignored).
# ------------------------------------

SE_D = 100
SE_CTRL_D = [100, 400, 1000]
SE_CTRL_SCALE = [1.2, 0.3, 0.12]
SE_LAMBDA_BG = 0.3
SE_TREAT_PLUS = {b"chr1": [200, 1000, 1005, 1010, 1020, 1030, 1040, 1050,
                           1060, 1070, 1400, 1405, 1410, 1420, 2000, 2010,
                           2030, 2050, 2080, 2100, 2150, 2200, 2250, 2600,
                           3990, 20000],
                 b"chr2": [100, 400, 410, 420, 430, 440, 450],
                 b"chr3": [50]}
SE_TREAT_MINUS = {b"chr1": [1150, 1160, 1170, 1180, 1190, 1200, 1210, 1520,
                            1530, 1535, 2300, 2320, 2350, 2400, 2450, 3000,
                            4500],
                  b"chr2": [520, 530, 540, 550, 900],
                  b"chr3": [300]}
SE_CTRL_PLUS = {b"chr1": [300, 900, 2500, 4100, 4700, 20500],
                b"chr2": [150, 700, 1000],
                b"chr4": [500]}
SE_CTRL_MINUS = {b"chr1": [1250, 3300, 4800],
                 b"chr2": [600, 1100]}
COMMON = (b"chr1", b"chr2")


def se_tracks():
    return (build_fw(SE_TREAT_PLUS, SE_TREAT_MINUS),
            build_fw(SE_CTRL_PLUS, SE_CTRL_MINUS))


def _clip_to(ivs, length):
    """Coordinates beyond the chromosome length are set to the length."""
    return [(min(s, length), min(e, length), w) for s, e, w in ivs]


def se_ref(d=SE_D, ctrl_d_s=SE_CTRL_D, scales=SE_CTRL_SCALE,
           lambda_bg=SE_LAMBDA_BG, treat_scale=1.0, end_shift=0,
           pseudocount=1.0, no_lambda=False, with_control=True,
           rlengths=None):
    treat = {c: se_treat_intervals(SE_TREAT_PLUS[c], SE_TREAT_MINUS[c], d,
                                   end_shift) for c in COMMON}
    if rlengths:
        treat = {c: _clip_to(iv, rlengths[c]) for c, iv in treat.items()}
    cp, cm = ((SE_CTRL_PLUS, SE_CTRL_MINUS) if with_control
              else (SE_TREAT_PLUS, SE_TREAT_MINUS))
    chroms = COMMON if with_control else (b"chr1", b"chr2", b"chr3")
    if not with_control:
        treat[b"chr3"] = se_treat_intervals(SE_TREAT_PLUS[b"chr3"],
                                            SE_TREAT_MINUS[b"chr3"], d,
                                            end_shift)
    windows = {c: [(se_ctrl_intervals(cp[c], cm[c], w), s)
                   for w, s in zip(ctrl_d_s, scales)] for c in chroms}
    if rlengths:
        windows = {c: [(_clip_to(iv, rlengths[c]), s) for iv, s in ws]
                   for c, ws in windows.items()}
    return RefCaller(treat, windows, treat_scale=treat_scale,
                     lambda_bg=lambda_bg, pseudocount=pseudocount,
                     no_lambda=no_lambda)


def se_caller(tmp_path, save_bedGraph=False, **kw):
    treat, ctrl = se_tracks()
    args = dict(d=SE_D, ctrl_d_s=list(SE_CTRL_D),
                ctrl_scaling_factor_s=list(SE_CTRL_SCALE),
                lambda_bg=SE_LAMBDA_BG, save_bedGraph=save_bedGraph,
                bedGraph_filename_prefix="NAME",
                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                bedGraph_control_filename=str(tmp_path / "c.bdg"),
                cutoff_analysis_filename=str(tmp_path / "cutoff.txt"))
    args.update(kw)
    return CallerFromAlignments(treat, ctrl, **args)


def test_reference_scenario_is_valid():
    """The scenario keeps every control window past the end of the
    treatment pileup (so no truncation is involved) and has enough
    background for the q-score table to reach 0."""
    for c in COMMON:
        t_end = max(e for _, e, _ in se_treat_intervals(
            SE_TREAT_PLUS[c], SE_TREAT_MINUS[c], SE_D))
        for w in SE_CTRL_D:
            w_end = max(e for _, e, _ in se_ctrl_intervals(
                SE_CTRL_PLUS[c], SE_CTRL_MINUS[c], w))
            assert w_end >= t_end
    ref = se_ref()
    assert min(ref.qtable.values()) == 0.0


# ------------------------------------
# CallerFromAlignments: constructor
# ------------------------------------

def test_constructor_rejects_unknown_track_message():
    with pytest.raises(Exception, match="Should be FWTrack or PETrackI/II object!"):
        CallerFromAlignments(None, None)


@pytest.mark.parametrize("bad", [1, "chr1", [1, 2], {b"chr1": []}])
def test_constructor_rejects_other_types(bad):
    with pytest.raises(Exception, match="Should be FWTrack or PETrackI/II object!"):
        CallerFromAlignments(bad, None)


def test_defaults_with_zero_lambda_raise(caller_tmp):
    """With the default lambda_bg of 0 a region without control reads has
    lambda 0 and the Poisson score refuses it."""
    treat = build_fw({b"chr1": [100]}, {})
    ctrl = build_fw({b"chr1": [50000]}, {})
    calc = CallerFromAlignments(treat, ctrl)
    with pytest.raises(AssertionError, match="Lambda must > 0"):
        calc.call_peaks(["p"], [5.0])


def test_stderr_on_is_accepted_and_ignored(tmp_path, caller_tmp):
    a = se_caller(tmp_path, stderr_on=True).call_peaks(["p"], [3.0], min_length=100)
    b = se_caller(tmp_path, stderr_on=False).call_peaks(["p"], [3.0], min_length=100)
    assert peak_rows(a) == peak_rows(b)


def test_only_common_chromosomes_are_called(tmp_path, caller_tmp):
    """chr3 (treatment only) and chr4 (control only) are skipped, even with
    a cutoff of 0 that every region passes."""
    peaks = se_caller(tmp_path).call_peaks(["s"], [-1e9], min_length=1)
    assert sorted(peaks.get_chr_names()) == [b"chr1", b"chr2"]


def test_no_control_uses_treatment_as_control(tmp_path, caller_tmp):
    treat = build_fw(SE_TREAT_PLUS, SE_TREAT_MINUS)
    calc = CallerFromAlignments(treat, None, d=SE_D, ctrl_d_s=[1000],
                                ctrl_scaling_factor_s=[0.1], lambda_bg=0.25,
                                save_bedGraph=True,
                                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                                bedGraph_control_filename=str(tmp_path / "c.bdg"))
    peaks = calc.call_peaks(["p"], [3.0], min_length=100)
    ref = se_ref(ctrl_d_s=[1000], scales=[0.1], lambda_bg=0.25,
                 with_control=False)
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    assert_peaks_match(peak_rows(peaks), ref.peaks(["p"], [3.0], 100, 50))


# ------------------------------------
# CallerFromAlignments: bedGraph output (pileup and control lambda)
# ------------------------------------

def test_bedgraph_treat_pileup(tmp_path, caller_tmp):
    se_caller(tmp_path, save_bedGraph=True).call_peaks(["p"], [3.0])
    assert read_lines(tmp_path / "t.bdg") == se_ref().bedgraph("treat")


def test_bedgraph_control_lambda_is_max_of_windows(tmp_path, caller_tmp):
    se_caller(tmp_path, save_bedGraph=True).call_peaks(["p"], [3.0])
    assert read_lines(tmp_path / "c.bdg") == se_ref().bedgraph("ctrl")


@pytest.mark.parametrize("ctrl_d_s,scales,lambda_bg", [
    ([100], [1.2], 0.3),
    ([400], [0.3], 0.05),
    ([100, 1000], [1.2, 0.12], 0.3),
    ([1000, 100, 400], [0.12, 1.2, 0.3], 0.3),
    ([100, 400, 1000], [1.2, 0.3, 0.12], 2.0),
    ([101, 399, 1001], [1.0, 0.5, 0.25], 0.1),
])
def test_bedgraph_control_lambda_windows(tmp_path, caller_tmp, ctrl_d_s,
                                         scales, lambda_bg):
    """Each window extends a control tag by w//2 to the left and w - w//2
    to the right (odd sizes included); lambda_bg floors the result."""
    se_caller(tmp_path, save_bedGraph=True, ctrl_d_s=ctrl_d_s,
              ctrl_scaling_factor_s=scales,
              lambda_bg=lambda_bg).call_peaks(["p"], [3.0])
    ref = se_ref(ctrl_d_s=ctrl_d_s, scales=scales, lambda_bg=lambda_bg)
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")


@pytest.mark.parametrize("d", [50, 100, 147, 250])
def test_bedgraph_treat_extension_d(tmp_path, caller_tmp, d):
    se_caller(tmp_path, save_bedGraph=True, d=d).call_peaks(["p"], [3.0])
    assert read_lines(tmp_path / "t.bdg") == se_ref(d=d).bedgraph("treat")


@pytest.mark.parametrize("end_shift", [-30, 0, 20, 150])
def test_bedgraph_treat_end_shift(tmp_path, caller_tmp, end_shift):
    """end_shift moves the 5' end towards 3' before extending by d."""
    se_caller(tmp_path, save_bedGraph=True,
              end_shift=end_shift).call_peaks(["p"], [3.0])
    assert read_lines(tmp_path / "t.bdg") == \
        se_ref(end_shift=end_shift).bedgraph("treat")


def test_bedgraph_clipped_at_chromosome_length(tmp_path, caller_tmp):
    """With chromosome lengths set on both tracks, coordinates past the end
    are moved to the end: on chr1 (4300 bp) the - read at 4500 vanishes,
    so the treatment pileup ends at 4090, and the wider control windows
    stop at 4300; on chr2 (1200 bp) only control windows are cut."""
    treat, ctrl = se_tracks()
    rl = {b"chr1": 4300, b"chr2": 1200, b"chr3": 1000, b"chr4": 1000}
    treat.set_rlengths(rl)
    ctrl.set_rlengths(rl)
    calc = CallerFromAlignments(treat, ctrl, d=SE_D, ctrl_d_s=list(SE_CTRL_D),
                                ctrl_scaling_factor_s=list(SE_CTRL_SCALE),
                                lambda_bg=SE_LAMBDA_BG, save_bedGraph=True,
                                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                                bedGraph_control_filename=str(tmp_path / "c.bdg"))
    peaks = calc.call_peaks(["p"], [3.0], min_length=100)
    ref = se_ref(rlengths=rl)
    t = read_lines(tmp_path / "t.bdg")
    assert t == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    assert [x for x in t if x.startswith("chr1\t")][-1].split("\t")[2] == "4090"
    assert_peaks_match(peak_rows(peaks), ref.peaks(["p"], [3.0], 100, 50))


def test_bedgraph_treat_scaling_factor(tmp_path, caller_tmp):
    se_caller(tmp_path, save_bedGraph=True,
              treat_scaling_factor=0.4).call_peaks(["p"], [3.0])
    ref = se_ref(treat_scale=0.4)
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")


@pytest.mark.parametrize("treat_scale", [1.0, 0.5])
def test_bedgraph_spmr(tmp_path, caller_tmp, treat_scale):
    """With save_SPMR both files are divided by the depth in millions:
    the treatment's when treat_scaling_factor is 1, else the control's."""
    se_caller(tmp_path, save_bedGraph=True, save_SPMR=True,
              treat_scaling_factor=treat_scale).call_peaks(["p"], [3.0])
    n_treat = len(sum(SE_TREAT_PLUS.values(), [])) + len(sum(SE_TREAT_MINUS.values(), []))
    n_ctrl = len(sum(SE_CTRL_PLUS.values(), [])) + len(sum(SE_CTRL_MINUS.values(), []))
    denom = (n_treat if treat_scale == 1.0 else n_ctrl) / 1e6
    ref = se_ref(treat_scale=treat_scale)
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat", denom)
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl", denom)


def test_bedgraph_no_lambda_is_constant_lambda_bg(tmp_path, caller_tmp):
    """Empty ctrl_d_s (or empty scaling factors) disables the local lambda:
    the control is lambda_bg over the whole treatment pileup."""
    for kw in (dict(ctrl_d_s=[]), dict(ctrl_scaling_factor_s=[])):
        se_caller(tmp_path, save_bedGraph=True, lambda_bg=0.7,
                  **kw).call_peaks(["p"], [3.0])
        ref = se_ref(lambda_bg=0.7, no_lambda=True)
        assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
        assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")


def test_bedgraph_written_once(tmp_path, caller_tmp):
    """The bedGraph files are written by the first call only."""
    calc = se_caller(tmp_path, save_bedGraph=True)
    calc.call_peaks(["p"], [3.0])
    (tmp_path / "t.bdg").write_text("sentinel\n")
    calc.call_peaks(["p"], [3.0])
    assert read_lines(tmp_path / "t.bdg") == ["sentinel"]


def test_bedgraph_not_written_without_save(tmp_path, caller_tmp):
    se_caller(tmp_path, save_bedGraph=False).call_peaks(["p"], [3.0])
    assert not (tmp_path / "t.bdg").exists()
    assert not (tmp_path / "c.bdg").exists()


def test_bedgraph_written_by_call_broadpeaks(tmp_path, caller_tmp):
    se_caller(tmp_path, save_bedGraph=True).call_broadpeaks(["p"], [5.0], [2.0])
    ref = se_ref()
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")


def test_enable_trackline_writes_track_lines(tmp_path, caller_tmp):
    """The first line of each file is a UCSC track line."""
    calc = se_caller(tmp_path, save_bedGraph=True)
    calc.enable_trackline()
    calc.call_peaks(["p"], [3.0])
    t = read_lines(tmp_path / "t.bdg")
    c = read_lines(tmp_path / "c.bdg")
    assert t[0].startswith('track type=bedGraph name="treatment pileup"')
    assert c[0].startswith('track type=bedGraph name="control lambda"')
    ref = se_ref()
    assert t[1:] == ref.bedgraph("treat")
    assert c[1:] == ref.bedgraph("ctrl")


def test_trackline_off_by_default(tmp_path, caller_tmp):
    se_caller(tmp_path, save_bedGraph=True).call_peaks(["p"], [3.0])
    assert not read_lines(tmp_path / "t.bdg")[0].startswith("track")


# ------------------------------------
# CallerFromAlignments.call_peaks: narrow peaks
# ------------------------------------

@pytest.mark.parametrize("cutoff", [2.0, 3.0, 5.0, 8.0, 12.0])
def test_call_peaks_pscore(tmp_path, caller_tmp, cutoff):
    peaks = se_caller(tmp_path).call_peaks(["p"], [cutoff], min_length=100,
                                           max_gap=50)
    expected = se_ref().peaks(["p"], [cutoff], 100, 50)
    assert_peaks_match(peak_rows(peaks), expected)


@pytest.mark.parametrize("min_length,max_gap", [
    (1, 0), (50, 10), (100, 50), (200, 50), (300, 100), (400, 400),
    (1000, 2000)])
def test_call_peaks_min_length_max_gap(tmp_path, caller_tmp, min_length,
                                       max_gap):
    peaks = se_caller(tmp_path).call_peaks(["p"], [2.0],
                                           min_length=min_length,
                                           max_gap=max_gap)
    expected = se_ref().peaks(["p"], [2.0], min_length, max_gap)
    assert_peaks_match(peak_rows(peaks), expected)


def test_call_peaks_defaults_min_length_200_max_gap_50(tmp_path, caller_tmp):
    peaks = se_caller(tmp_path).call_peaks(["p"], [2.0])
    assert_peaks_match(peak_rows(peaks), se_ref().peaks(["p"], [2.0], 200, 50))


@pytest.mark.parametrize("cutoff", [1.0, 2.0, 5.0, 10.0])
def test_call_peaks_qscore(tmp_path, caller_tmp, cutoff):
    peaks = se_caller(tmp_path).call_peaks(["q"], [cutoff], min_length=100)
    expected = se_ref().peaks(["q"], [cutoff], 100, 50)
    assert expected or cutoff == 10.0
    assert_peaks_match(peak_rows(peaks), expected)


@pytest.mark.parametrize("cutoff", [2.0, 5.0, 10.0])
def test_call_peaks_fold_change(tmp_path, caller_tmp, cutoff):
    peaks = se_caller(tmp_path).call_peaks(["f"], [cutoff], min_length=100)
    assert_peaks_match(peak_rows(peaks), se_ref().peaks(["f"], [cutoff], 100, 50))


@pytest.mark.parametrize("cutoff", [0.5, 3.0, 8.0])
def test_call_peaks_subtraction(tmp_path, caller_tmp, cutoff):
    peaks = se_caller(tmp_path).call_peaks(["s"], [cutoff], min_length=100)
    assert_peaks_match(peak_rows(peaks), se_ref().peaks(["s"], [cutoff], 100, 50))


@pytest.mark.parametrize("syms,cutoffs", [
    (["p", "f"], [3.0, 4.0]),
    (["p", "q"], [5.0, 2.0]),
    (["f", "s"], [8.0, 2.0]),
    (["s", "p"], [9.0, 3.0]),
])
def test_call_peaks_multiple_cutoffs(tmp_path, caller_tmp, syms, cutoffs):
    """Regions are where any score exceeds its cutoff; a peak is kept only
    when its summit reaches every cutoff."""
    peaks = se_caller(tmp_path).call_peaks(syms, cutoffs, min_length=100)
    assert_peaks_match(peak_rows(peaks), se_ref().peaks(syms, cutoffs, 100, 50))


def test_call_peaks_high_cutoff_gives_no_peaks(tmp_path, caller_tmp):
    peaks = se_caller(tmp_path).call_peaks(["p"], [1000.0])
    assert isinstance(peaks, PeakIO)
    assert peaks.total == 0
    assert peak_rows(peaks) == []


def test_call_peaks_mismatched_cutoffs(tmp_path, caller_tmp):
    with pytest.raises(AssertionError,
                       match="number of functions and cutoffs should be the same!"):
        se_caller(tmp_path).call_peaks(["p", "q"], [3.0])


def test_call_peaks_unknown_symbol(tmp_path, caller_tmp):
    """An unknown symbol yields no score array, so the cutoffs have nothing
    to be applied to."""
    with pytest.raises(IndexError):
        se_caller(tmp_path).call_peaks(["x"], [3.0])


def test_call_peaks_repeatable(tmp_path, caller_tmp):
    """A second call reuses the q-score table and the cached pileups."""
    calc = se_caller(tmp_path)
    a = peak_rows(calc.call_peaks(["p"], [3.0], min_length=100))
    b = peak_rows(calc.call_peaks(["p"], [3.0], min_length=100))
    c = peak_rows(calc.call_peaks(["q"], [2.0], min_length=100))
    assert a == b
    assert_peaks_match(c, se_ref().peaks(["q"], [2.0], 100, 50))


LARGE_BASE = 2_100_000_000


def _large_coordinate_case(tmp_path):
    """Reads near the top of the int32 range (2.1e9) on a second
    chromosome."""
    base = LARGE_BASE
    tp = {b"chr1": [100, 5000], b"chrB": [base + i for i in range(0, 60, 10)]}
    tm = {b"chrB": [base + 300]}
    cp = {b"chr1": [6000], b"chrB": [base + 2000]}
    calc = CallerFromAlignments(build_fw(tp, tm), build_fw(cp, {}), d=150,
                                ctrl_d_s=[150, 1500],
                                ctrl_scaling_factor_s=[1.0, 0.1],
                                lambda_bg=0.2, save_bedGraph=True,
                                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                                bedGraph_control_filename=str(tmp_path / "c.bdg"))
    peaks = calc.call_peaks(["p"], [3.0], min_length=50)
    ref = RefCaller({c: se_treat_intervals(tp[c], tm.get(c, []), 150)
                     for c in tp},
                    {c: [(se_ctrl_intervals(cp[c], [], 150), 1.0),
                         (se_ctrl_intervals(cp[c], [], 1500), 0.1)]
                     for c in tp}, lambda_bg=0.2)
    return peak_rows(peaks), ref


def test_call_peaks_large_coordinates(tmp_path, caller_tmp):
    """Pileups, lambda and peak boundaries keep coordinates near 2.1e9."""
    rows, ref = _large_coordinate_case(tmp_path)
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    expected = ref.peaks(["p"], [3.0], 50, 50)
    assert [p["chrom"] for p in expected] == [b"chrB"]
    assert expected[0]["start"] > LARGE_BASE
    assert [(r["chrom"], r["start"], r["end"]) for r in rows] == \
        [(p["chrom"], p["start"], p["end"]) for p in expected]
    assert rows[0]["pileup"] == expected[0]["pileup"]
    assert rows[0]["pscore"] == pytest.approx(expected[0]["pscore"], abs=1e-4)


@pytest.mark.parametrize("pc", [0.0, 0.5, 1.0, 5.0])
def test_set_pseudocount_changes_fold_change(tmp_path, caller_tmp, pc):
    calc = se_caller(tmp_path)
    calc.set_pseudocount(pc)
    peaks = calc.call_peaks(["p"], [3.0], min_length=100)
    assert_peaks_match(peak_rows(peaks),
                       se_ref(pseudocount=pc).peaks(["p"], [3.0], 100, 50))


@pytest.mark.parametrize("pc", [0.5, 2.0])
def test_pseudocount_argument(tmp_path, caller_tmp, pc):
    peaks = se_caller(tmp_path, pseudocount=pc).call_peaks(["f"], [3.0],
                                                           min_length=100)
    assert_peaks_match(peak_rows(peaks),
                       se_ref(pseudocount=pc).peaks(["f"], [3.0], 100, 50))


def test_call_peaks_no_lambda(tmp_path, caller_tmp):
    peaks = se_caller(tmp_path, ctrl_d_s=[], lambda_bg=0.5).call_peaks(
        ["p"], [3.0], min_length=100)
    assert_peaks_match(peak_rows(peaks),
                       se_ref(lambda_bg=0.5, no_lambda=True).peaks(["p"], [3.0], 100, 50))


def test_call_peaks_end_shift_and_scaling(tmp_path, caller_tmp):
    peaks = se_caller(tmp_path, end_shift=20, treat_scaling_factor=0.6).call_peaks(
        ["p"], [2.0], min_length=100)
    assert_peaks_match(peak_rows(peaks),
                       se_ref(end_shift=20, treat_scale=0.6).peaks(["p"], [2.0], 100, 50))


def test_call_peaks_first_segment_peak_starts_at_zero(tmp_path, caller_tmp):
    """Reads starting at position 0 make the first segment significant;
    the peak then starts at 0."""
    treat = build_fw({b"chr1": [0, 0, 0, 0, 0, 0, 3000]}, {})
    ctrl = build_fw({b"chr1": [5000]}, {})
    calc = CallerFromAlignments(treat, ctrl, d=100, ctrl_d_s=[100, 1000],
                                ctrl_scaling_factor_s=[1.0, 0.1],
                                lambda_bg=0.2)
    rows = peak_rows(calc.call_peaks(["p"], [5.0], min_length=50))
    assert [(r["start"], r["end"], r["summit"]) for r in rows] == [(0, 100, 50)]


def test_call_peaks_flat_top_summit_is_middle_segment(tmp_path, caller_tmp):
    """Three segments share the highest pileup (cut by control windows);
    the summit is the midpoint of the middle one."""
    # treatment: 4 reads at 1000 -> pileup 4 on [1000, 1100); control tags
    # at 1020 and 1040 with a 2 bp window cut that range at 1019/1021 and
    # 1039/1041 without changing lambda (lambda_bg dominates).
    treat = build_fw({b"chr1": [1000] * 4 + [3000]}, {})
    ctrl = build_fw({b"chr1": [1020, 1040, 3500]}, {})
    calc = CallerFromAlignments(treat, ctrl, d=100, ctrl_d_s=[2],
                                ctrl_scaling_factor_s=[0.01], lambda_bg=0.2)
    rows = peak_rows(calc.call_peaks(["p"], [3.0], min_length=50))
    # segments with pileup 4: [1000,1019) [1019,1021) [1021,1039)
    # [1039,1041) [1041,1100); the middle one is [1021,1039) -> 1030
    assert [(r["start"], r["end"], r["summit"]) for r in rows] == \
        [(1000, 1100, 1030)]
    assert rows[0]["pileup"] == 4.0
    assert rows[0]["pscore"] == pytest.approx(
        ref_pscore(4, float(F32(0.2))), abs=1e-4)


# ------------------------------------
# CallerFromAlignments.call_peaks(call_summits=True)
# ------------------------------------

def _two_hill_tracks():
    """One region with two separated hills of different heights, built from
    staggered + tags (pileup rises one read at a time)."""
    plus = list(range(1000, 1200, 20)) + list(range(1500, 1660, 20))
    minus = list(range(1300, 1500, 20)) + list(range(1760, 1920, 20))
    treat = build_fw({b"chr1": plus + [5000]}, {b"chr1": minus})
    ctrl = build_fw({b"chr1": [6000]}, {})
    return treat, ctrl, plus, minus


def test_call_summits_two_hills(tmp_path, caller_tmp):
    treat, ctrl, plus, minus = _two_hill_tracks()
    calc = CallerFromAlignments(treat, ctrl, d=200, ctrl_d_s=[200, 2000],
                                ctrl_scaling_factor_s=[1.0, 0.1],
                                lambda_bg=0.5)
    rows = peak_rows(calc.call_peaks(["p"], [3.0], min_length=200,
                                     max_gap=200, call_summits=True))
    ref = RefCaller({b"chr1": se_treat_intervals(plus + [5000], minus, 200)},
                    {b"chr1": [(se_ctrl_intervals([6000], [], 200), 1.0),
                               (se_ctrl_intervals([6000], [], 2000), 0.1)]},
                    lambda_bg=0.5)
    (region,) = ref.peaks(["p"], [3.0], 200, 200)
    # every sub-peak shares the boundaries of the region
    assert len(rows) == 2
    assert {(r["start"], r["end"]) for r in rows} == \
        {(region["start"], region["end"])}
    # each summit sits on a local maximum of the treatment pileup, one per
    # hill, and reports that segment's values
    segs = ref.segs[b"chr1"]
    for r in rows:
        (seg,) = [s for s in segs if s[0] <= r["summit"] < s[1]]
        assert r["pileup"] == seg[2]
        assert r["pscore"] == pytest.approx(ref_pscore(seg[2], seg[3]), abs=1e-4)
        assert r["fc"] == pytest.approx((seg[2] + 1) / (seg[3] + 1), rel=1e-6)
    assert rows[0]["summit"] < 1450 < rows[1]["summit"]
    hill1 = max(s[2] for s in segs if s[1] <= 1450)
    hill2 = max(s[2] for s in segs if s[0] >= 1450)
    assert [rows[0]["pileup"], rows[1]["pileup"]] == [hill1, hill2]


def test_call_summits_flat_peak(tmp_path, caller_tmp):
    """A flat-topped peak (6 reads at 1000, d 150) is reported once with
    the same boundaries, values and summit (1075, the middle) as without
    call_summits.

    Pins the current output for the sub-peak summit. Upstream 1622eb8
    (#749, issue #747) removed the double 10-bp padding shift of the
    summit-search signal; the plateau is now symmetric in the padded
    window, maxima() returns two candidates (start + 49 and start + 99),
    enforce_peakyness rejects both as too flat, and the caller falls back
    to the summit without sub-peaks (the middle). Before that commit the
    single candidate start + 49 (1049) was reported.
    """
    kw = dict(d=150, ctrl_d_s=[150, 1500], ctrl_scaling_factor_s=[1.0, 0.1],
              lambda_bg=0.3)
    a = peak_rows(CallerFromAlignments(
        build_fw({b"chr1": [1000] * 6 + [3000]}, {}),
        build_fw({b"chr1": [4000]}, {}), **kw).call_peaks(
        ["p"], [3.0], min_length=100, call_summits=True))
    b = peak_rows(CallerFromAlignments(
        build_fw({b"chr1": [1000] * 6 + [3000]}, {}),
        build_fw({b"chr1": [4000]}, {}), **kw).call_peaks(
        ["p"], [3.0], min_length=100, call_summits=False))
    assert [(r["start"], r["end"], r["summit"]) for r in b] == [(1000, 1150, 1075)]
    assert [(r["start"], r["end"], r["summit"]) for r in a] == [(1000, 1150, 1075)]
    for k in ("pileup", "pscore", "qscore", "fc"):
        assert a[0][k] == b[0][k]
    assert b[0]["pileup"] == 6.0


def _uneven_hills():
    """Hill A (max pileup 12) and hill B (max 22) 140 bp apart, d = 200,
    lambda 0.5 everywhere (lambda_bg)."""
    plus = list(range(1000, 1120, 20)) + list(range(1500, 1620, 10))
    minus = list(range(1260, 1380, 20)) + list(range(1800, 1920, 10))
    treat = build_fw({b"chr1": plus}, {b"chr1": minus})
    ctrl = build_fw({b"chr1": [6000]}, {})
    calc = CallerFromAlignments(treat, ctrl, d=200, ctrl_d_s=[200, 2000],
                                ctrl_scaling_factor_s=[1.0, 0.1],
                                lambda_bg=0.5)
    ref = RefCaller({b"chr1": se_treat_intervals(plus, minus, 200)},
                    {b"chr1": [(se_ctrl_intervals([6000], [], 200), 1.0),
                               (se_ctrl_intervals([6000], [], 2000), 0.1)]},
                    lambda_bg=0.5)
    return calc, ref


def test_call_summits_uneven_hills_single_cutoff(caller_tmp):
    calc, ref = _uneven_hills()
    rows = peak_rows(calc.call_peaks(["p"], [3.0], min_length=100,
                                     max_gap=300, call_summits=True))
    (region,) = ref.peaks(["p"], [3.0], 100, 300)
    assert len(rows) == 2
    assert {(r["start"], r["end"]) for r in rows} == \
        {(region["start"], region["end"])}
    segs = ref.segs[b"chr1"]
    hills = [max(s[2] for s in segs if s[1] <= 1450),
             max(s[2] for s in segs if s[0] >= 1450)]
    assert hills == [12.0, 22.0]
    assert [r["pileup"] for r in rows] == hills
    for r, top in zip(rows, hills):
        (seg,) = [s for s in segs if s[0] <= r["summit"] < s[1]]
        assert seg[2] == top
        assert r["pscore"] == pytest.approx(ref_pscore(seg[2], seg[3]), abs=1e-4)
        assert r["qscore"] == pytest.approx(
            ref.qtable[ref_pscore(seg[2], seg[3])], abs=1e-4)


def test_call_summits_positions_pinned(caller_tmp):
    """Exact sub-peak summits of the two hill layouts.

    Pins the current output.  The summits come from the zero crossings of
    a Savitzky-Golay derivative over the padded pileup followed by the
    peakyness filter; the tests above check independently that each lies
    on the top of its hill, but the exact base within the top is not
    practical to derive by hand.
    """
    calc, _ = _uneven_hills()
    rows = peak_rows(calc.call_peaks(["p"], [3.0], min_length=100,
                                     max_gap=300, call_summits=True))
    assert [(r["start"], r["end"], r["summit"]) for r in rows] == \
        [(1060, 1880, 1179), (1060, 1880, 1704)]


# ------------------------------------
# q-scores
# ------------------------------------

def _uniform_tracks():
    """Ten + reads tiling [0, 1000) with d = 100: pileup 1 everywhere, so
    the chromosome is one segment and the p-score table has one value."""
    return build_fw({b"chr1": list(range(0, 1000, 100))}, {})


def test_qscore_uniform_coverage_pscore(caller_tmp):
    calc = CallerFromAlignments(_uniform_tracks(), None, d=100, ctrl_d_s=[],
                                ctrl_scaling_factor_s=[], lambda_bg=0.01)
    rows = peak_rows(calc.call_peaks(["p"], [2.0], min_length=100))
    assert [(r["start"], r["end"]) for r in rows] == [(0, 1000)]
    assert rows[0]["pscore"] == pytest.approx(ref_pscore(1, F32(0.01)), abs=1e-4)


def test_qscore_of_every_significant_segment(tmp_path, caller_tmp):
    """With min_length 1 and max_gap -1 (adjacent segments are not merged)
    every segment above the cutoff is a peak, so each reported q-score is
    checked against the table."""
    ref = se_ref()
    vals = sorted(ref.qtable.items(), reverse=True)
    assert all(a[1] >= b[1] for a, b in zip(vals, vals[1:]))
    rows = peak_rows(se_caller(tmp_path).call_peaks(["p"], [2.0], min_length=1,
                                                    max_gap=-1))
    expected = ref.peaks(["p"], [2.0], 1, -1)
    assert len(expected) > 5
    assert_peaks_match(rows, expected)


# ------------------------------------
# CallerFromAlignments.destroy and temporary pileup files
# ------------------------------------

def test_destroy_removes_pileup_files(tmp_path, caller_tmp):
    calc = se_caller(tmp_path)
    calc.call_peaks(["p"], [3.0])
    assert len(list(caller_tmp.iterdir())) == 2      # one per chromosome
    calc.destroy()
    assert list(caller_tmp.iterdir()) == []
    calc.destroy()                                    # nothing left: no error


def test_call_peaks_after_destroy_recomputes(tmp_path, caller_tmp):
    calc = se_caller(tmp_path)
    a = peak_rows(calc.call_peaks(["p"], [3.0], min_length=100))
    calc.destroy()
    b = peak_rows(calc.call_peaks(["p"], [3.0], min_length=100))
    assert a == b
    assert len(list(caller_tmp.iterdir())) == 2


def test_destroy_before_calling(tmp_path, caller_tmp):
    se_caller(tmp_path).destroy()
    assert list(caller_tmp.iterdir()) == []


def test_enable_trackline_returns_none(tmp_path):
    calc = se_caller(tmp_path)
    assert calc.enable_trackline() is None
    assert calc.set_pseudocount(2.0) is None


# ------------------------------------
# log messages
# ------------------------------------

def _messages(caplog):
    return [re.sub(r"^\[\d+ MB\] ", "", r.getMessage()) for r in caplog.records
            if r.name == "MACS3.Signal.CallPeakUnit"]


@pytest.mark.parametrize("spmr,treat_scale,line", [
    (False, 1.0, "#3   Pileup will be based on sequencing depth in treatment."),
    (False, 0.5, "#3   Pileup will be based on sequencing depth in control."),
    (True, 1.0, "#3   --SPMR is requested, so pileup will be normalized by "
                "sequencing depth in million reads."),
])
def test_call_peaks_log_messages(tmp_path, caller_tmp, caplog, spmr,
                                 treat_scale, line):
    caplog.set_level(logging.INFO, logger="MACS3.Signal.CallPeakUnit")
    se_caller(tmp_path, save_bedGraph=True, save_SPMR=spmr,
              treat_scaling_factor=treat_scale).call_peaks(["p"], [3.0])
    assert _messages(caplog) == [
        "#3 Pre-compute pvalue-qvalue table...",
        "#3 In the peak calling step, the following will be performed simultaneously:",
        "#3   Write bedGraph files for treatment pileup (after scaling if necessary)... NAME_treat_pileup.bdg",
        "#3   Write bedGraph files for control lambda (after scaling if necessary)... NAME_control_lambda.bdg",
        line,
        "#3 Call peaks for each chromosome...",
    ]


def test_call_peaks_log_messages_cached_table(tmp_path, caller_tmp, caplog):
    calc = se_caller(tmp_path)
    calc.call_peaks(["p"], [3.0])
    caplog.set_level(logging.INFO, logger="MACS3.Signal.CallPeakUnit")
    calc.call_peaks(["p"], [3.0])
    assert _messages(caplog) == ["#3 Call peaks for each chromosome..."]


def test_call_broadpeaks_log_messages(tmp_path, caller_tmp, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.Signal.CallPeakUnit")
    se_caller(tmp_path, save_bedGraph=True).call_broadpeaks(["p"], [5.0], [2.0])
    assert _messages(caplog) == [
        "#3 Pre-compute pvalue-qvalue table...",
        "#3 In the peak calling step, the following will be performed simultaneously:",
        "#3   Write bedGraph files for treatment pileup (after scaling if necessary)... NAME_treat_pileup.bdg",
        "#3   Write bedGraph files for control lambda (after scaling if necessary)... NAME_control_lambda.bdg",
        "#3 Call peaks for each chromosome...",
    ]


# ------------------------------------
# cutoff analysis (call_peaks / call_broadpeaks with cutoff_analysis=True)
# ------------------------------------

def _cutoff_lines(rows):
    return ["%.2f\t%d\t%d\t%.2f" % (k, n, ln, ave) for k, q, n, ln, ave in rows]


def _read_cutoff(path):
    lines = read_lines(path)
    assert lines[0] == "pscore\tqscore\tnpeaks\tlpeaks\tavelpeak"
    rows = [x.split("\t") for x in lines[1:]]
    return ["\t".join([r[0]] + r[2:]) for r in rows], [float(r[1]) for r in rows]


@pytest.mark.parametrize("max_gap,min_length", [(50, 100), (0, 1), (200, 300)])
def test_cutoff_analysis_file(tmp_path, caller_tmp, max_gap, min_length):
    """Rows for every -log10 p cutoff from 9.9 down to 0.3 (step 0.3) that
    calls at least one peak.  This scenario has no significant first
    segment."""
    calc = se_caller(tmp_path)
    calc.call_peaks(["p"], [3.0], min_length=min_length, max_gap=max_gap,
                    cutoff_analysis=True)
    rows = se_ref().cutoff_rows(max_gap, min_length)
    rows = [r for r in rows if r[0] > 0.6]  # first segment is significant below
    got, q = _read_cutoff(tmp_path / "cutoff.txt")
    got = [g for g in got if float(g.split("\t")[0]) > 0.6]
    assert got == _cutoff_lines(rows)
    assert q[:len(rows)] == pytest.approx([r[1] for r in rows], abs=0.006)


def test_cutoff_analysis_same_peaks(tmp_path, caller_tmp):
    a = peak_rows(se_caller(tmp_path).call_peaks(["p"], [3.0], min_length=100,
                                                 cutoff_analysis=True))
    b = peak_rows(se_caller(tmp_path).call_peaks(["p"], [3.0], min_length=100))
    assert a == b


def test_cutoff_analysis_broad_uses_lvl2_gap(tmp_path, caller_tmp):
    se_caller(tmp_path).call_broadpeaks(["p"], [5.0], [2.0], min_length=100,
                                        lvl2_max_gap=300, cutoff_analysis=True)
    rows = [r for r in se_ref().cutoff_rows(300, 100) if r[0] > 0.6]
    got, _ = _read_cutoff(tmp_path / "cutoff.txt")
    got = [g for g in got if float(g.split("\t")[0]) > 0.6]
    assert got == _cutoff_lines(rows)


def test_cutoff_analysis_log_message(tmp_path, caller_tmp, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.Signal.CallPeakUnit")
    se_caller(tmp_path).call_peaks(["p"], [3.0], cutoff_analysis=True)
    msgs = _messages(caplog)
    assert msgs[:2] == ["#3 Pre-compute pvalue-qvalue table...",
                        "#3 Cutoff vs peaks called will be analyzed!"]


def test_cutoff_analysis_unwritable_file(tmp_path, caller_tmp):
    calc = se_caller(tmp_path, cutoff_analysis_filename=str(tmp_path / "no" / "x.txt"))
    with pytest.raises(FileNotFoundError):
        calc.call_peaks(["p"], [3.0], cutoff_analysis=True)


# ------------------------------------
# CallerFromAlignments.call_broadpeaks
# ------------------------------------

def broad_rows(bpeaks):
    rows = []
    for chrom in sorted(bpeaks.peaks):
        for p in bpeaks.peaks[chrom]:
            rows.append(dict(chrom=chrom, start=p["start"], end=p["end"],
                             score=p["score"], pileup=p["pileup"],
                             pscore=p["pscore"], fc=p["fc"],
                             qscore=p["qscore"], thickStart=p["thickStart"],
                             thickEnd=p["thickEnd"], blockNum=p["blockNum"],
                             blockSizes=p["blockSizes"],
                             blockStarts=p["blockStarts"]))
    return rows


def assert_broad_match(got, expected):
    keys = ("chrom", "start", "end", "thickStart", "thickEnd", "blockNum",
            "blockSizes", "blockStarts")
    assert [tuple(g[k] for k in keys) for g in got] == \
        [tuple(e[k] for k in keys) for e in expected]
    for g, e in zip(got, expected):
        assert g["pileup"] == pytest.approx(e["pileup"], rel=1e-5)
        assert g["fc"] == pytest.approx(e["fc"], rel=1e-5)
        assert g["pscore"] == pytest.approx(e["pscore"], abs=1e-4)
        assert g["qscore"] == pytest.approx(e["qscore"], abs=1e-4)
        assert g["score"] == pytest.approx(e["qscore"], abs=1e-4)


@pytest.mark.parametrize("lvl1,lvl2,min_length,gap1,gap2", [
    (5.0, 2.0, 100, 50, 400),
    (8.0, 2.0, 100, 50, 400),
    (5.0, 1.0, 100, 30, 2000),
    (3.0, 2.0, 50, 0, 100),
    (12.0, 3.0, 50, 50, 200),
])
def test_call_broadpeaks_pscore(tmp_path, caller_tmp, lvl1, lvl2, min_length,
                                gap1, gap2):
    bp = se_caller(tmp_path).call_broadpeaks(["p"], [lvl1], [lvl2],
                                             min_length=min_length,
                                             lvl1_max_gap=gap1,
                                             lvl2_max_gap=gap2)
    assert isinstance(bp, BroadPeakIO)
    expected = se_ref().broadpeaks(["p"], [lvl1], [lvl2], min_length, gap1, gap2)
    assert expected
    assert_broad_match(broad_rows(bp), expected)


def test_call_broadpeaks_qscore(tmp_path, caller_tmp):
    bp = se_caller(tmp_path).call_broadpeaks(["q"], [5.0], [1.5], min_length=100)
    assert_broad_match(broad_rows(bp),
                       se_ref().broadpeaks(["q"], [5.0], [1.5], 100, 50, 400))


def test_call_broadpeaks_defaults(tmp_path, caller_tmp):
    """Defaults: min_length 200, lvl1_max_gap 50, lvl2_max_gap 400."""
    bp = se_caller(tmp_path).call_broadpeaks(["p"], [5.0], [2.0])
    assert_broad_match(broad_rows(bp),
                       se_ref().broadpeaks(["p"], [5.0], [2.0], 200, 50, 400))


def test_call_broadpeaks_nothing_above_level1(tmp_path, caller_tmp):
    bp = se_caller(tmp_path).call_broadpeaks(["p"], [1000.0], [2.0])
    assert bp.peaks == {}


@pytest.mark.parametrize("lvl1,lvl2", [([5.0], [2.0, 1.0]), ([5.0, 1.0], [2.0])])
def test_call_broadpeaks_mismatched_cutoffs(tmp_path, caller_tmp, lvl1, lvl2):
    with pytest.raises(AssertionError,
                       match="number of functions and cutoffs should be the same!"):
        se_caller(tmp_path).call_broadpeaks(["p"], lvl1, lvl2)


def test_call_broadpeaks_weak_region_without_strong_block(tmp_path, caller_tmp):
    """On a chromosome with a level-1 region, a level-2 region holding no
    level-1 region is reported with two 1-bp blocks at its ends."""
    # cluster A (strong) near 1000 and cluster B (weak) near 3000
    plus = [1000] * 8 + [3000, 3000, 3000]
    treat = build_fw({b"chr1": plus}, {})
    ctrl = build_fw({b"chr1": [4000]}, {})
    calc = CallerFromAlignments(treat, ctrl, d=200, ctrl_d_s=[200, 2000],
                                ctrl_scaling_factor_s=[1.0, 0.1],
                                lambda_bg=0.3)
    rows = broad_rows(calc.call_broadpeaks(["p"], [8.0], [2.0], min_length=100))
    assert [(r["start"], r["end"], r["blockNum"], r["blockSizes"],
             r["blockStarts"], r["thickStart"], r["thickEnd"]) for r in rows] == [
        (1000, 1200, 1, b"200", b"0", b"1000", b"1200"),
        (3000, 3200, 2, b"1,1", b"0,199", b"3000", b"3200")]


# ------------------------------------
# paired-end tracks (PETrackI, BEDPE/BAMPE)
# ------------------------------------

PE_TREAT = {b"chr1": [(100, 300), (950, 1150), (960, 1160), (980, 1200),
                      (1000, 1210), (1010, 1190), (1020, 1230), (2500, 2700),
                      (4000, 4180)],
            b"chr2": [(200, 380), (400, 610), (410, 600), (430, 640),
                      (900, 1100)]}
PE_CTRL = {b"chr1": [(300, 500), (1200, 1400), (2600, 2800), (4300, 4500),
                     (4700, 4900)],
           b"chr2": [(100, 300), (700, 900), (1200, 1400)]}
PE_D = [200, 1000]
PE_SCALE = [0.5, 0.1]


def pe_ref(lambda_bg=0.2, with_control=True):
    treat = {c: pe_treat_intervals(PE_TREAT[c]) for c in PE_TREAT}
    src = PE_CTRL if with_control else PE_TREAT
    windows = {c: [(pe1_ctrl_intervals(src[c], w), s)
                   for w, s in zip(PE_D, PE_SCALE)] for c in PE_TREAT}
    return RefCaller(treat, windows, lambda_bg=lambda_bg)


def test_pe_bedgraph_and_peaks(tmp_path, caller_tmp):
    """Treatment pileup is fragment coverage; the control piles up both
    ends of each control fragment in each window."""
    calc = CallerFromAlignments(build_pe1(PE_TREAT), build_pe1(PE_CTRL),
                                ctrl_d_s=PE_D, ctrl_scaling_factor_s=PE_SCALE,
                                lambda_bg=0.2, save_bedGraph=True,
                                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                                bedGraph_control_filename=str(tmp_path / "c.bdg"))
    peaks = calc.call_peaks(["p"], [3.0], min_length=100)
    ref = pe_ref()
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    assert_peaks_match(peak_rows(peaks), ref.peaks(["p"], [3.0], 100, 50))


def test_pe_without_control(tmp_path, caller_tmp):
    calc = CallerFromAlignments(build_pe1(PE_TREAT), None,
                                ctrl_d_s=PE_D, ctrl_scaling_factor_s=PE_SCALE,
                                lambda_bg=0.2, save_bedGraph=True,
                                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                                bedGraph_control_filename=str(tmp_path / "c.bdg"))
    peaks = calc.call_peaks(["p"], [3.0], min_length=100)
    ref = pe_ref(with_control=False)
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    assert_peaks_match(peak_rows(peaks), ref.peaks(["p"], [3.0], 100, 50))


def test_pe_broadpeaks(tmp_path, caller_tmp):
    calc = CallerFromAlignments(build_pe1(PE_TREAT), build_pe1(PE_CTRL),
                                ctrl_d_s=PE_D, ctrl_scaling_factor_s=PE_SCALE,
                                lambda_bg=0.2)
    bp = calc.call_broadpeaks(["p"], [5.0], [2.0], min_length=100)
    assert_broad_match(broad_rows(bp),
                       pe_ref().broadpeaks(["p"], [5.0], [2.0], 100, 50, 400))


# ------------------------------------
# fragment files (PETrackII, FRAG)
# ------------------------------------

FRAG_TREAT = {b"chr1": [(5100, 5300, 1), (5950, 6150, 2), (5960, 6160, 1),
                        (5980, 6200, 3), (6000, 6210, 1), (7500, 7700, 1),
                        (9000, 9180, 2)]}
FRAG_CTRL = {b"chr1": [(5200, 5400, 1), (6200, 6400, 2), (7600, 7800, 1),
                       (9300, 9500, 1), (9700, 9900, 1)]}


def frag_ref(with_control, ds, scales, lambda_bg=0.2):
    treat = {b"chr1": pe_treat_intervals(FRAG_TREAT[b"chr1"])}
    src = FRAG_CTRL if with_control else FRAG_TREAT
    windows = {b"chr1": [(pe2_ctrl_intervals(src[b"chr1"], w), s)
                         for w, s in zip(ds, scales)]}
    return RefCaller(treat, windows, lambda_bg=lambda_bg)


def _frag_caller(tmp_path, with_control, ds, scales):
    ctrl = build_pe2(FRAG_CTRL) if with_control else None
    return CallerFromAlignments(build_pe2(FRAG_TREAT), ctrl, ctrl_d_s=ds,
                                ctrl_scaling_factor_s=scales, lambda_bg=0.2,
                                save_bedGraph=True,
                                bedGraph_treat_filename=str(tmp_path / "t.bdg"),
                                bedGraph_control_filename=str(tmp_path / "c.bdg"))


def test_frag_single_window_without_control(tmp_path, caller_tmp):
    """Fragment counts weight the pileup; with one control window the
    lambda is that window over both fragment ends."""
    peaks = _frag_caller(tmp_path, False, [2000], [0.1]).call_peaks(
        ["p"], [3.0], min_length=100)
    ref = frag_ref(False, [2000], [0.1])
    assert read_lines(tmp_path / "t.bdg") == ref.bedgraph("treat")
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")
    assert_peaks_match(peak_rows(peaks), ref.peaks(["p"], [3.0], 100, 50))


def test_frag_single_window_with_control(tmp_path, caller_tmp):
    _frag_caller(tmp_path, True, [400], [0.5]).call_peaks(["p"], [3.0])
    ref = frag_ref(True, [400], [0.5])
    assert read_lines(tmp_path / "c.bdg") == ref.bedgraph("ctrl")


def test_frag_control_lambda_max_over_windows(tmp_path, caller_tmp):
    """Regression test: with FRAG input the max over control windows read
    the strided positions through a raw pointer, mixing position and value
    bits.

    Fixed upstream in 5456b02 (#737): PETrackII.pileup_a_chromosome_c now
    builds each window's pileup as contiguous arrays
    (pileup_from_LRC_centers_as_list) before over_two_pv_array.
    """
    ds, scales = [200, 1000, 4000], [0.5, 0.1, 0.025]
    _frag_caller(tmp_path, True, ds, scales).call_peaks(["p"], [3.0])
    assert read_lines(tmp_path / "c.bdg") == frag_ref(True, ds, scales).bedgraph("ctrl")


