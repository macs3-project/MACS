#!/usr/bin/env python

"""Module Description: Test the bdgcmp subcommand
(MACS3/Commands/bdgcmp_cmd.py) and its option checks in
MACS3/Utilities/OptValidator.py (opt_validate_bdgcmp).

Every scoring method is compared with a reference formula computed in
the test (scipy.stats.poisson for ppois, the MACS3 Benjamini-Hochberg
procedure for qpois, log10 ratios and likelihood ratios written out)
on a tiny treatment/control pair. Output content is checked in-process
(the real argparse parser plus ``run``); exit codes, logs and argparse
errors are checked with ``macs3`` in a subprocess.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import logging
import os
from math import log, log10
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import poisson

from MACS3.Commands.bdgcmp_cmd import (run as bdgcmp_run)

# ------------------------------------
# helpers
# ------------------------------------

LOG10_E = 0.43429448190325176
METHODS = ["ppois", "qpois", "subtract", "logFE", "FE", "logLR", "slogLR",
           "max"]

# Treatment and control, continuous from 0. chr3 exists only in the
# treatment and chrX only in the control, so neither is written. The
# control chr2 is longer than the treatment chr2, so chr2 stops at 150.
T_ROWS = [("chr1", 0, 100, 1), ("chr1", 100, 200, 5),
          ("chr1", 200, 300, 12), ("chr1", 300, 400, 3),
          ("chr2", 0, 50, 2), ("chr2", 50, 150, 8),
          ("chr3", 0, 100, 4)]
C_ROWS = [("chr1", 0, 150, 2), ("chr1", 150, 400, 4),
          ("chr2", 0, 200, 3),
          ("chrX", 0, 100, 1)]
# Paired intervals over the common chromosomes (union of breakpoints):
# chr1 [0,100) (1,2) [100,150) (5,2) [150,200) (5,4) [200,300) (12,4)
#      [300,400) (3,4); chr2 [0,50) (2,3) [50,150) (8,3)


@pytest.fixture(autouse=True)
def _restore_optvalidator_level():
    """opt_validate_bdgcmp sets the level of the OptValidator logger from
    --verbose; restore it so in-process runs leave no trace."""
    lg = logging.getLogger("MACS3.Utilities.OptValidator")
    level = lg.level
    yield
    lg.setLevel(level)


def call(argparser, argv):
    """Parse ``macs3 bdgcmp <argv>`` with the real parser and run it
    in-process. Returns the options namespace."""
    options = argparser.parse_args(["bdgcmp"] + [str(a) for a in argv])
    bdgcmp_run(options)
    return options


def lines(path):
    """All lines of a text file without line endings."""
    return Path(path).read_text().splitlines()


def rows_of(path):
    """(chrom, start, end, value) from a bedGraph without a track line."""
    out = []
    for x in lines(path):
        c, s, e, v = x.split("\t")
        out.append((c, int(s), int(e), float(v)))
    return out


def paired_intervals(trows, crows):
    """(chrom, start, end, t, c) over the chromosomes in both tracks, on
    the union of breakpoints, stopping at the shorter track's end."""
    out = []
    common = sorted({r[0] for r in trows} & {r[0] for r in crows})
    for chrom in common:
        t = [(e, v) for c, s, e, v in trows if c == chrom]
        k = [(e, v) for c, s, e, v in crows if c == chrom]
        end = min(t[-1][0], k[-1][0])
        pre = 0
        for b in sorted({e for e, _ in t + k if e <= end}):
            tv = next(v for e, v in t if e >= b)
            cv = next(v for e, v in k if e >= b)
            out.append((chrom, pre, b, tv, cv))
            pre = b
    return out


def logLR_asym(x, y):
    """log10 likelihood ratio of x over y (negative for depletion)."""
    if x > y:
        return (x * (log(x) - log(y)) + y - x) * LOG10_E
    if x < y:
        return (x * (-log(x) + log(y)) - y + x) * LOG10_E
    return 0.0


def logLR_sym(x, y):
    """Symmetric log10 likelihood ratio between two ChIP signals."""
    if x > y:
        return (x * (log(x) - log(y)) + y - x) * LOG10_E
    if y > x:
        return (y * (log(x) - log(y)) + y - x) * LOG10_E
    return 0.0


def pscore(obs, lam):
    """-log10 P(X > obs) for X ~ Poisson(lam), from scipy."""
    return -log10(poisson.sf(obs, lam))


def macs_qscores(scored):
    """MACS3's q-value transform of -log10 p-scores.

    ``scored`` is a list of (length, pscore). Following MACS3's
    definition: N = total length; going through distinct p-scores from
    the largest, q = p + log10(k) - log10(N) with k = 1 + total length
    of larger p-scores, q is capped by the previous q (monotone), and
    once q <= 0 it and all smaller p-scores get 0.
    """
    stat = {}
    for ln, p in scored:
        key = round(p, 5)
        stat[key] = stat.get(key, 0) + ln
    n = sum(stat.values())
    k = 1
    pre_q = float("inf")
    table = {key: 0.0 for key in stat}
    for key in sorted(stat, reverse=True):
        q = min(key + log10(k) - log10(n), pre_q)
        if q <= 0:
            break
        table[key] = q
        pre_q = q
        k += stat[key]
    return [table[round(p, 5)] for _, p in scored]


def reference_scores(method, trows, crows, sfactor=1.0, pseudocount=0.0):
    """Expected (chrom, start, end, value) rows of ``bdgcmp -m method``.

    Values are scaled by -S first (only when it differs from 1), in
    float32 as MACS3 stores them; the pseudocount is added to both
    tracks for ppois, logFE, FE, logLR and slogLR (not for subtract or
    max); ppois uses the integer part of the treatment value. Adjacent
    intervals whose scores differ by <= 1e-5 are merged into one line
    carrying the first value, as ScoreTrackII.write_bedGraph does.
    """
    iv = paired_intervals(trows, crows)
    f32 = np.float32
    scale = abs(sfactor - 1) > 1e-6
    pc = float(f32(pseudocount))
    vals = []
    for chrom, s, e, t, c in iv:
        t = f32(t) * f32(sfactor) if scale else f32(t)
        c = f32(c) * f32(sfactor) if scale else f32(c)
        tp = float(f32(t + f32(pc)))
        cp = float(f32(c + f32(pc)))
        if method in ("ppois", "qpois"):
            v = pscore(int(tp), cp)
        elif method == "subtract":
            v = float(f32(t - c))
        elif method == "logFE":
            v = log10(tp / cp)
        elif method == "FE":
            v = float(f32(t + f32(pc)) / f32(c + f32(pc)))
        elif method == "logLR":
            v = logLR_asym(tp, cp)
        elif method == "slogLR":
            v = logLR_sym(tp, cp)
        elif method == "max":
            v = float(max(t, c))
        vals.append(v)
    if method == "qpois":
        vals = macs_qscores([(e - s, v) for (_, s, e, _, _), v in
                             zip(iv, vals)])
    out = []
    for (chrom, s, e, _, _), v in zip(iv, vals):
        v = float(f32(v))
        if out and out[-1][0] == chrom and abs(out[-1][3] - v) <= 1e-5:
            out[-1] = (chrom, out[-1][1], e, out[-1][3])
        else:
            out.append((chrom, s, e, v))
    return out


def assert_rows_match(got, expected, tol):
    """Same intervals; values within ``tol``."""
    assert [r[:3] for r in got] == [r[:3] for r in expected]
    for g, x in zip(got, expected):
        assert g[3] == pytest.approx(x[3], abs=tol), (g, x)


TOL = {"ppois": 1e-4, "qpois": 1e-4, "subtract": 6e-6, "logFE": 2e-5,
       "FE": 6e-6, "logLR": 2e-5, "slogLR": 2e-5, "max": 6e-6}


@pytest.fixture
def tc_files(write_bedgraph):
    """The treatment and control bedGraphs above."""
    return (write_bedgraph(T_ROWS, name="t.bdg"),
            write_bedgraph(C_ROWS, name="c.bdg"))


# ------------------------------------
# every method against its reference
# ------------------------------------

@pytest.mark.parametrize("sfactor, pseudocount", [
    (1.0, 0.0), (1.0, 1.0), (0.5, 0.0), (2.0, 0.5)])
@pytest.mark.parametrize("method", METHODS)
def test_method_matches_reference(macs3_argparser, tc_files, tmp_path,
                                  method, sfactor, pseudocount):
    """Each -m with -S and -p against the reference formula."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", method, "-S", sfactor,
                           "-p", pseudocount, "--outdir", tmp_path,
                           "--o-prefix", "X"])
    got = rows_of(tmp_path / ("X_%s.bdg" % method))
    exp = reference_scores(method, T_ROWS, C_ROWS, sfactor, pseudocount)
    assert_rows_match(got, exp, TOL[method])


def test_subtract_exact(macs3_argparser, tc_files, tmp_path):
    """subtract, by hand: chr1 1-2, 5-2, 5-4, 12-4, 3-4; chr2 2-3, 8-3.
    The output has no track line and stops where the treatment chr2
    ends (150); chr3 and chrX are not common and are not written."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "subtract", "--outdir",
                           tmp_path, "--o-prefix", "X"])
    assert lines(tmp_path / "X_subtract.bdg") == [
        "chr1\t0\t100\t-1.00000", "chr1\t100\t150\t3.00000",
        "chr1\t150\t200\t1.00000", "chr1\t200\t300\t8.00000",
        "chr1\t300\t400\t-1.00000",
        "chr2\t0\t50\t-1.00000", "chr2\t50\t150\t5.00000"]


def test_max_merges_equal_neighbours(macs3_argparser, tc_files, tmp_path):
    """max, by hand: chr1 2, 5, 5, 12, 4 -> the two 5s ([100,150) and
    [150,200)) become one line; chr2 3, 8."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "max", "--outdir",
                           tmp_path, "--o-prefix", "X"])
    assert lines(tmp_path / "X_max.bdg") == [
        "chr1\t0\t100\t2.00000", "chr1\t100\t200\t5.00000",
        "chr1\t200\t300\t12.00000", "chr1\t300\t400\t4.00000",
        "chr2\t0\t50\t3.00000", "chr2\t50\t150\t8.00000"]


def test_fe_with_pseudocount_exact(macs3_argparser, tc_files, tmp_path):
    """FE -p 1, by hand: (t+1)/(c+1) = 2/3, 6/3, 6/5, 13/5, 4/5; chr2
    3/4, 9/4."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "FE", "-p", "1",
                           "--outdir", tmp_path, "--o-prefix", "X"])
    assert lines(tmp_path / "X_FE.bdg") == [
        "chr1\t0\t100\t0.66667", "chr1\t100\t150\t2.00000",
        "chr1\t150\t200\t1.20000", "chr1\t200\t300\t2.60000",
        "chr1\t300\t400\t0.80000",
        "chr2\t0\t50\t0.75000", "chr2\t50\t150\t2.25000"]


def test_qpois_by_hand(macs3_argparser, tc_files, tmp_path):
    """qpois with N = 550 bp: p-scores in decreasing order are 3.5627
    (12 vs 4, 100 bp), 2.4199 (8 vs 3, 100), 1.7808 (5 vs 2, 50),
    0.6678 (5 vs 4, 50), 0.2468 (3 vs 4, 100), ...
    q = p + log10(k) - log10(550) with k = 1, 101, 201, 251, 301 gives
    0.8223, 1.6842 -> capped 0.8223, 1.3436 -> capped 0.8223, 0.3272,
    and -0.015 -> 0 (and 0 for everything smaller)."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "qpois", "--outdir",
                           tmp_path, "--o-prefix", "X"])
    got = rows_of(tmp_path / "X_qpois.bdg")
    q1 = pscore(12, 4) - log10(550)
    q4 = pscore(5, 4) + log10(251) - log10(550)
    exp = [("chr1", 0, 100, 0.0), ("chr1", 100, 150, q1),
           ("chr1", 150, 200, q4), ("chr1", 200, 300, q1),
           ("chr1", 300, 400, 0.0), ("chr2", 0, 50, 0.0),
           ("chr2", 50, 150, q1)]
    assert_rows_match(got, exp, 1e-4)
    assert q1 == pytest.approx(0.8223, abs=1e-3)
    assert q4 == pytest.approx(0.3272, abs=1e-3)


def test_ppois_truncates_treatment(macs3_argparser, write_bedgraph,
                                   tmp_path):
    """ppois uses int(treatment + pseudocount): 2.9 scores like 2 and
    3.0 like 3 against lambda 1."""
    t = write_bedgraph([("chr1", 0, 100, 2.9), ("chr1", 100, 200, 2),
                        ("chr1", 200, 300, 3)], name="t.bdg")
    c = write_bedgraph([("chr1", 0, 300, 1)], name="c.bdg")
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "ppois", "--outdir",
                           tmp_path, "--o-prefix", "X"])
    got = rows_of(tmp_path / "X_ppois.bdg")
    # [0,100) and [100,200) have the same score and are one line
    assert_rows_match(got, [("chr1", 0, 200, pscore(2, 1)),
                            ("chr1", 200, 300, pscore(3, 1))], 1e-4)


def test_pseudocount_handles_zeros(macs3_argparser, write_bedgraph,
                                   tmp_path):
    """With -p 1, zero treatment and zero control give finite scores:
    logLR(1, 3) < 0, logLR(4, 1) > 0, ppois(int(0+1), 0+1)."""
    t = write_bedgraph([("chr1", 0, 100, 0), ("chr1", 100, 200, 3)],
                       name="t.bdg")
    c = write_bedgraph([("chr1", 0, 100, 2), ("chr1", 100, 200, 0)],
                       name="c.bdg")
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "logLR", "ppois", "-p",
                           "1", "--outdir", tmp_path, "--o-prefix", "Z"])
    assert_rows_match(rows_of(tmp_path / "Z_logLR.bdg"),
                      [("chr1", 0, 100, logLR_asym(1, 3)),
                       ("chr1", 100, 200, logLR_asym(4, 1))], 2e-5)
    assert_rows_match(rows_of(tmp_path / "Z_ppois.bdg"),
                      [("chr1", 0, 100, pscore(1, 3)),
                       ("chr1", 100, 200, pscore(4, 1))], 1e-4)


@pytest.mark.parametrize("method", ["subtract", "max"])
def test_pseudocount_not_used(macs3_argparser, tc_files, tmp_path, method):
    """subtract and max ignore -p."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", method, "--outdir",
                           tmp_path, "-o", "a.bdg"])
    call(macs3_argparser, ["-t", t, "-c", c, "-m", method, "-p", "5",
                           "--outdir", tmp_path, "-o", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")


@pytest.mark.parametrize("sfactor, expected", [
    # (t - c) * S, by hand
    ("0.5", ["-0.50000", "1.50000", "0.50000", "4.00000", "-0.50000",
             "-0.50000", "2.50000"]),
    ("2", ["-2.00000", "6.00000", "2.00000", "16.00000", "-2.00000",
           "-2.00000", "10.00000"]),
    # within 1e-6 of 1: no scaling
    ("1.0000001", ["-1.00000", "3.00000", "1.00000", "8.00000", "-1.00000",
                   "-1.00000", "5.00000"]),
])
def test_scaling_factor_subtract(macs3_argparser, tc_files, tmp_path,
                                 sfactor, expected):
    """-S multiplies both tracks before scoring."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "subtract", "-S",
                           sfactor, "--outdir", tmp_path, "-o", "s.bdg"])
    assert [x.split("\t")[3] for x in lines(tmp_path / "s.bdg")] == expected


def test_fe_scaling_cancels_without_pseudocount(macs3_argparser, tc_files,
                                                tmp_path):
    """FE = tS/cS: -S changes nothing when -p is 0."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "FE", "--outdir",
                           tmp_path, "-o", "a.bdg"])
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "FE", "-S", "4",
                           "--outdir", tmp_path, "-o", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")


def test_identical_tracks(macs3_argparser, tc_files, tmp_path):
    """Treatment == control: subtract 0, FE 1, logFE 0, logLR 0, slogLR
    0, one merged line per chromosome."""
    t, _ = tc_files
    call(macs3_argparser, ["-t", t, "-c", t, "-m", "subtract", "FE",
                           "logFE", "logLR", "slogLR", "--outdir",
                           tmp_path, "--o-prefix", "I"])
    for method, val in [("subtract", "0.00000"), ("FE", "1.00000"),
                        ("logFE", "0.00000"), ("logLR", "0.00000"),
                        ("slogLR", "0.00000")]:
        assert lines(tmp_path / ("I_%s.bdg" % method)) == [
            "chr1\t0\t400\t" + val, "chr2\t0\t150\t" + val,
            "chr3\t0\t100\t" + val], method


def test_long_option_names(macs3_argparser, tc_files, tmp_path):
    """--tfile/--cfile/--scaling-factor/--pseudocount/--method/--ofile
    are the same options as -t/-c/-S/-p/-m/-o."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-S", "2", "-p", "1", "-m",
                           "FE", "logLR", "--outdir", tmp_path, "-o", "a1",
                           "a2"])
    call(macs3_argparser, ["--tfile", t, "--cfile", c, "--scaling-factor",
                           "2", "--pseudocount", "1", "--method", "FE",
                           "logLR", "--outdir", tmp_path, "--ofile", "b1",
                           "b2"])
    assert lines(tmp_path / "a1") == lines(tmp_path / "b1")
    assert lines(tmp_path / "a2") == lines(tmp_path / "b2")
    assert_rows_match(rows_of(tmp_path / "a1"),
                      reference_scores("FE", T_ROWS, C_ROWS, 2.0, 1.0), 6e-6)


def test_empty_treatment(macs3_argparser, tc_files, tmp_path):
    """An empty treatment bedGraph has no chromosome in common with the
    control: an empty output file."""
    _, c = tc_files
    t = tmp_path / "empty.bdg"
    t.write_text("")
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "ppois", "--outdir",
                           tmp_path, "-o", "e.bdg"])
    assert (tmp_path / "e.bdg").read_text() == ""


def test_no_common_chromosome(macs3_argparser, write_bedgraph, tmp_path):
    """No chromosome in common: an empty output file."""
    t = write_bedgraph([("chrA", 0, 100, 3)], name="t.bdg")
    c = write_bedgraph([("chrB", 0, 100, 1)], name="c.bdg")
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "FE", "--outdir",
                           tmp_path, "-o", "e.bdg"])
    assert (tmp_path / "e.bdg").read_text() == ""


# ------------------------------------
# several methods, file naming
# ------------------------------------

def test_all_methods_with_prefix(macs3_argparser, tc_files, tmp_path):
    """--o-prefix X writes X_<method>.bdg for each -m, each equal to the
    reference (qpois right after ppois reuses its p-scores)."""
    t, c = tc_files
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-t", t, "-c", c, "-m"] + METHODS
         + ["--outdir", outdir, "--o-prefix", "X"])
    assert sorted(os.listdir(outdir)) == sorted("X_%s.bdg" % m
                                                for m in METHODS)
    for m in METHODS:
        assert_rows_match(rows_of(outdir / ("X_%s.bdg" % m)),
                          reference_scores(m, T_ROWS, C_ROWS), TOL[m])


def test_ofile_pairs_with_methods_in_order(macs3_argparser, tc_files,
                                           tmp_path):
    """-o names are used in the order of -m."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "max", "subtract", "FE",
                           "--outdir", tmp_path, "-o", "one.bdg", "two.bdg",
                           "three.bdg"])
    assert_rows_match(rows_of(tmp_path / "one.bdg"),
                      reference_scores("max", T_ROWS, C_ROWS), 6e-6)
    assert_rows_match(rows_of(tmp_path / "two.bdg"),
                      reference_scores("subtract", T_ROWS, C_ROWS), 6e-6)
    assert_rows_match(rows_of(tmp_path / "three.bdg"),
                      reference_scores("FE", T_ROWS, C_ROWS), 6e-6)


def test_repeated_method_with_prefix_written_once(macs3_argparser,
                                                  tc_files, tmp_path):
    """-m FE FE --o-prefix X writes X_FE.bdg once."""
    t, c = tc_files
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-t", t, "-c", c, "-m", "FE", "FE", "--outdir",
                           outdir, "--o-prefix", "X"])
    assert os.listdir(outdir) == ["X_FE.bdg"]


@pytest.mark.parametrize("order", [["qpois"], ["ppois", "qpois"],
                                   ["qpois", "ppois"], ["FE", "qpois"]])
def test_qpois_independent_of_order(macs3_argparser, tc_files, tmp_path,
                                    order):
    """qpois equals the reference whether or not ppois was computed
    first in the same run."""
    t, c = tc_files
    call(macs3_argparser, ["-t", t, "-c", c, "-m"] + order
         + ["--outdir", tmp_path, "--o-prefix", "X"])
    assert_rows_match(rows_of(tmp_path / "X_qpois.bdg"),
                      reference_scores("qpois", T_ROWS, C_ROWS), 1e-4)


# ------------------------------------
# option validation and the command line
# ------------------------------------

def test_cli_ofile_count_must_match_methods(run_macs3, parse_log, tc_files,
                                            tmp_path):
    """Two -m and one -o: opt_validate_bdgcmp logs an error, exits 1 and
    nothing is written."""
    t, c = tc_files
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "ppois", "FE",
                     "--outdir", tmp_path, "-o", "a.bdg"], timeout=60)
    assert res.returncode == 1
    assert parse_log(res.stderr) == [
        ("ERROR", "The number and the order of arguments for --ofile must "
                  "be the same as for -m.")]
    assert not (tmp_path / "a.bdg").exists()


def test_cli_more_ofiles_than_methods(run_macs3, parse_log, tc_files,
                                      tmp_path):
    """One -m and two -o: the same error."""
    t, c = tc_files
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "FE", "--outdir",
                     tmp_path, "-o", "a.bdg", "b.bdg"], timeout=60)
    assert res.returncode == 1
    assert parse_log(res.stderr)[-1] == (
        "ERROR", "The number and the order of arguments for --ofile must "
                 "be the same as for -m.")


def test_invalid_method_in_process(macs3_argparser, tc_files, tmp_path,
                                   caplog):
    """A method outside the list (only reachable without argparse's
    choices check) is rejected by opt_validate_bdgcmp: exit 1."""
    t, c = tc_files
    options = macs3_argparser.parse_args(["bdgcmp", "-t", t, "-c", c, "-m",
                                          "FE", "--outdir", str(tmp_path),
                                          "--o-prefix", "X"])
    options.method = ["FE", "bogus"]
    with caplog.at_level(logging.INFO):
        with pytest.raises(SystemExit) as exc:
            bdgcmp_run(options)
    assert exc.value.code == 1
    assert caplog.records[-1].getMessage().endswith("Invalid method: bogus")
    assert not os.path.exists(tmp_path / "X_FE.bdg")


def test_cli_log_messages(run_macs3, parse_log, tc_files, tmp_path):
    """A two-method run with -S logs exactly these INFO messages."""
    t, c = tc_files
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "FE", "subtract",
                     "-S", "2", "--outdir", tmp_path, "--o-prefix", "L"],
                    timeout=60)
    assert res.returncode == 0, res.stderr
    fe = os.path.join(str(tmp_path), "L_FE.bdg")
    sub = os.path.join(str(tmp_path), "L_subtract.bdg")
    assert parse_log(res.stderr) == [
        ("INFO", "Read and build treatment bedGraph..."),
        ("INFO", "Read and build control bedGraph..."),
        ("INFO", "Build ScoreTrackII..."),
        ("INFO", "Values in your input bedGraph files will be multiplied "
                 "by 2.000000 ..."),
        ("INFO", "Calculate scores comparing treatment and control by "
                 "'FE'..."),
        ("INFO", "Write bedGraph of scores..."),
        ("INFO", "Finished 'FE'! Please check '%s'!" % fe),
        ("INFO", "Calculate scores comparing treatment and control by "
                 "'subtract'..."),
        ("INFO", "Write bedGraph of scores..."),
        ("INFO", "Finished 'subtract'! Please check '%s'!" % sub)]
    assert res.stdout == ""


@pytest.mark.parametrize("verbose, n_info", [("0", 0), ("1", 0), ("2", 6),
                                             ("3", 6)])
def test_cli_verbose(run_macs3, parse_log, tc_files, tmp_path, verbose,
                     n_info):
    """--verbose 0/1 hide INFO messages; 2 and 3 show them (bdgcmp has no
    debug messages)."""
    t, c = tc_files
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "FE", "--verbose",
                     verbose, "--outdir", tmp_path, "--o-prefix", "V"],
                    timeout=60)
    assert res.returncode == 0
    log = parse_log(res.stderr)
    assert [lv for lv, _ in log] == ["INFO"] * n_info
    assert (tmp_path / "V_FE.bdg").exists()


def test_cli_verbose_zero_still_shows_errors(run_macs3, parse_log,
                                             tc_files, tmp_path):
    """At --verbose 0 the -o/-m count error is still logged."""
    t, c = tc_files
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "FE", "max",
                     "--verbose", "0", "--outdir", tmp_path, "-o", "a"],
                    timeout=60)
    assert res.returncode == 1
    assert [lv for lv, _ in parse_log(res.stderr)] == ["ERROR"]


def test_cli_zero_control_without_pseudocount(run_macs3, write_bedgraph,
                                              tmp_path):
    """ppois with a control value of 0 and no pseudocount stops with
    Poisson's 'Lambda must > 0' assertion (exit 1)."""
    t = write_bedgraph([("chr1", 0, 100, 3)], name="t.bdg")
    c = write_bedgraph([("chr1", 0, 100, 0)], name="c.bdg")
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "ppois", "--outdir",
                     tmp_path, "--o-prefix", "Z"], timeout=60)
    assert res.returncode == 1
    assert "Lambda must > 0" in res.stderr


def test_cli_outdir_is_created(run_macs3, tc_files, tmp_path):
    """A missing --outdir is created and holds the outputs."""
    t, c = tc_files
    outdir = tmp_path / "x" / "y"
    res = run_macs3(["bdgcmp", "-t", t, "-c", c, "-m", "max", "FE",
                     "--outdir", outdir, "--o-prefix", "O"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert sorted(os.listdir(outdir)) == ["O_FE.bdg", "O_max.bdg"]


@pytest.mark.parametrize("which", ["-t", "-c"])
def test_cli_missing_input_file(run_macs3, tc_files, tmp_path, which):
    """A nonexistent -t or -c fails with FileNotFoundError (exit 1)."""
    t, c = tc_files
    missing = str(tmp_path / "nope.bdg")
    args = {"-t": t, "-c": c}
    args[which] = missing
    res = run_macs3(["bdgcmp", "-t", args["-t"], "-c", args["-c"], "-m",
                     "FE", "--outdir", tmp_path, "--o-prefix", "M"],
                    timeout=60)
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr
    assert missing in res.stderr
    assert not (tmp_path / "M_FE.bdg").exists()


@pytest.mark.parametrize("argv, message", [
    (["-c", "c.bdg", "-m", "FE", "-o", "a"],
     "the following arguments are required: -t/--tfile"),
    (["-t", "t.bdg", "-m", "FE", "-o", "a"],
     "the following arguments are required: -c/--cfile"),
    (["-t", "t.bdg", "-c", "c.bdg", "-m", "FE"],
     "one of the arguments --o-prefix -o/--ofile is required"),
    (["-t", "t.bdg", "-c", "c.bdg", "-m", "FE", "--o-prefix", "X", "-o",
      "a"], "argument -o/--ofile: not allowed with argument --o-prefix"),
    (["-t", "t.bdg", "-c", "c.bdg", "-m", "bogus", "-o", "a"],
     "argument -m/--method: invalid choice: 'bogus'"),
    (["-t", "t.bdg", "-c", "c.bdg", "-m", "fe", "-o", "a"],
     "argument -m/--method: invalid choice: 'fe'"),
    (["-t", "t.bdg", "-c", "c.bdg", "-S", "big", "-o", "a"],
     "argument -S/--scaling-factor: invalid float value: 'big'"),
    (["-t", "t.bdg", "-c", "c.bdg", "-p", "one", "-o", "a"],
     "argument -p/--pseudocount: invalid float value: 'one'"),
    (["-t", "t.bdg", "-c", "c.bdg", "-m", "-o", "a"],
     "argument -m/--method: expected at least one argument"),
])
def test_cli_argparse_errors(run_macs3, argv, message):
    """Bad or missing arguments: exit 2, usage and the argparse error."""
    res = run_macs3(["bdgcmp"] + argv, timeout=60)
    assert res.returncode == 2
    assert res.stderr.startswith("usage: macs3 bdgcmp")
    assert res.stderr.rstrip().splitlines()[-1].startswith(
        "macs3 bdgcmp: error: " + message)


# ------------------------------------
# realistic run on the CTCF example
# ------------------------------------

def test_ctcf_pileups_match_standard(macs3_argparser, test_dir, tmp_path):
    """bdgcmp -m ppois FE -p 1 on the CTCF chr22 ChIP and control
    pileups (macs3 pileup --extsize 200) reproduces both upstream
    bedGraphs.

    Pins the current output. The reference files are MACS3's own
    standard results; scores over 190k intervals cannot be derived by
    hand line by line (ppois and FE are checked against scipy and the
    formula above).
    """
    pdir = test_dir / "standard_results_pileup"
    sdir = test_dir / "standard_results_bdgcmp"
    call(macs3_argparser, ["-t", pdir / "run_pileup_ChIP.bed.bdg",
                           "-c", pdir / "run_pileup_CTRL.bed.bdg",
                           "-m", "ppois", "FE", "-p", "1",
                           "--outdir", tmp_path, "--o-prefix", "run_bdgcmp"])
    for name in ("run_bdgcmp_ppois.bdg", "run_bdgcmp_FE.bdg"):
        assert lines(tmp_path / name) == lines(sdir / name), name
