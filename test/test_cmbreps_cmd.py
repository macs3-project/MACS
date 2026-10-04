#!/usr/bin/env python

"""Module Description: Test the cmbreps subcommand
(MACS3/Commands/cmbreps_cmd.py) and its option checks in
MACS3/Utilities/OptValidator.py (opt_validate_cmbreps): combining
replicate score tracks by Fisher's method, max or mean.

Output content is checked in-process (the real argparse parser plus
``run``) on tiny bedGraphs; Fisher's method is compared with
scipy.stats.chi2 on the combined -log10 p-values. Exit codes, logs and
argparse errors are checked with ``macs3`` in a subprocess.

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
from scipy.stats import chi2

from MACS3.Commands.cmbreps_cmd import (run as cmbreps_run)

# ------------------------------------
# helpers
# ------------------------------------

# Three replicates. chr3 exists only in R3 and is never written. The
# chr2 ends differ (200, 250, 200), so chr2 stops at 200.
R1 = [("chr1", 0, 100, 2), ("chr1", 100, 300, 4), ("chr2", 0, 200, 1)]
R2 = [("chr1", 0, 150, 3), ("chr1", 150, 300, 1), ("chr2", 0, 100, 5),
      ("chr2", 100, 250, 2)]
R3 = [("chr1", 0, 300, 6), ("chr2", 0, 200, 0), ("chr3", 0, 100, 1)]
# Paired intervals for R1, R2:
# chr1 [0,100) (2,3) [100,150) (4,3) [150,300) (4,1)
# chr2 [0,100) (1,5) [100,200) (1,2)


@pytest.fixture(autouse=True)
def _restore_optvalidator_level():
    """opt_validate_cmbreps sets the level of the OptValidator logger
    from --verbose; restore it so in-process runs leave no trace."""
    lg = logging.getLogger("MACS3.Utilities.OptValidator")
    level = lg.level
    yield
    lg.setLevel(level)


@pytest.fixture
def reps(write_bedgraph):
    """Paths of R1, R2, R3."""
    return [write_bedgraph(R1, name="r1.bdg"),
            write_bedgraph(R2, name="r2.bdg"),
            write_bedgraph(R3, name="r3.bdg")]


def call(argparser, argv):
    """Parse ``macs3 cmbreps <argv>`` with the real parser and run it
    in-process. Returns the options namespace."""
    options = argparser.parse_args(["cmbreps"] + [str(a) for a in argv])
    cmbreps_run(options)
    return options


def lines(path):
    """All lines of a text file without line endings."""
    return Path(path).read_text().splitlines()


def track(method):
    m = method.upper()
    return ('track type=bedGraph name="%s_combined_scores" '
            'description="Scores calculated by %s" visibility=2 '
            'alwaysZero=on' % (m, m))


def fisher(scores):
    """Fisher's combined -log10 p-value of -log10 p-values ``scores``:
    X = -2 sum(ln p) = 2 ln(10) sum(scores) ~ chi2 with 2k df."""
    x = 2 * log(10) * sum(scores)
    if x <= 0:
        return 0.0
    return -log10(chi2.sf(x, 2 * len(scores)))


def combine(tracks, func):
    """Expected (chrom, start, end, value) of combining ``tracks``:
    chromosomes in all tracks, the union of breakpoints up to the
    shortest track's end, ``func`` of the float32 values; neighbours
    with equal float32 results are one line."""
    common = sorted(set.intersection(*[{r[0] for r in t} for t in tracks]))
    out = []
    for chrom in common:
        per = [[(e, v) for c, s, e, v in t if c == chrom] for t in tracks]
        end = min(p[-1][0] for p in per)
        pre = 0
        for b in sorted({e for p in per for e, _ in p if e <= end}):
            vals = [float(np.float32(next(v for e, v in p if e >= b)))
                    for p in per]
            v = float(np.float32(func(vals)))
            if out and out[-1][0] == chrom and out[-1][3] == v:
                out[-1] = (chrom, out[-1][1], b, v)
            else:
                out.append((chrom, pre, b, v))
            pre = b
    return out


def rows_of(path):
    out = []
    for x in lines(path)[1:]:
        c, s, e, v = x.split("\t")
        out.append((c, int(s), int(e), float(v)))
    return out


# ------------------------------------
# max and mean, by hand
# ------------------------------------

@pytest.mark.parametrize("method, which, expected", [
    # max of R1, R2: chr1 3, 4, 4 (the two 4s become one line); chr2 5, 2
    ("max", [0, 1], ["chr1\t0\t100\t3.00000", "chr1\t100\t300\t4.00000",
                     "chr2\t0\t100\t5.00000", "chr2\t100\t200\t2.00000"]),
    # mean of R1, R2: chr1 2.5, 3.5, 2.5; chr2 3, 1.5
    ("mean", [0, 1], ["chr1\t0\t100\t2.50000", "chr1\t100\t150\t3.50000",
                      "chr1\t150\t300\t2.50000",
                      "chr2\t0\t100\t3.00000", "chr2\t100\t200\t1.50000"]),
    # max of R1, R2, R3: chr1 is 6 throughout; chr2 5, 2
    ("max", [0, 1, 2], ["chr1\t0\t300\t6.00000", "chr2\t0\t100\t5.00000",
                        "chr2\t100\t200\t2.00000"]),
    # mean of three: chr1 11/3, 13/3, 11/3; chr2 6/3, 3/3
    ("mean", [0, 1, 2], ["chr1\t0\t100\t3.66667", "chr1\t100\t150\t4.33333",
                         "chr1\t150\t300\t3.66667",
                         "chr2\t0\t100\t2.00000", "chr2\t100\t200\t1.00000"]),
])
def test_max_and_mean_exact(macs3_argparser, reps, tmp_path, method, which,
                            expected):
    """max and mean of 2 and 3 replicates, values taken as they are."""
    call(macs3_argparser, ["-i"] + [reps[k] for k in which]
         + ["-m", method, "--outdir", tmp_path, "-o", "cmb.bdg"])
    assert lines(tmp_path / "cmb.bdg") == [track(method)] + expected


# ------------------------------------
# Fisher's method against scipy
# ------------------------------------

@pytest.mark.parametrize("which", [[0, 1], [0, 1, 2], [2, 0, 1],
                                   [0, 1, 2, 0]])
def test_fisher_matches_chi2(macs3_argparser, reps, tmp_path, which):
    """Fisher's method for 2, 3 and 4 replicates (in any order) equals
    -log10 chi2.sf(2 ln10 sum, 2k)."""
    tracks = [[R1, R2, R3][k] for k in which]
    call(macs3_argparser, ["-i"] + [reps[k] for k in which]
         + ["-m", "fisher", "--outdir", tmp_path, "-o", "f.bdg"])
    got = rows_of(tmp_path / "f.bdg")
    exp = combine(tracks, fisher)
    assert [r[:3] for r in got] == [r[:3] for r in exp]
    for g, x in zip(got, exp):
        assert g[3] == pytest.approx(x[3], abs=2e-5), (g, x)


def test_fisher_two_by_hand(macs3_argparser, reps, tmp_path):
    """With 4 df the tail is exp(-x/2)(1 + x/2), so the combined score
    of sum S is S - log10(1 + S ln10): S = 5 -> 3.90264, S = 7 ->
    5.76654, S = 6 -> 4.82928, S = 3 -> 2.10195."""
    call(macs3_argparser, ["-i", reps[0], reps[1], "--outdir", tmp_path,
                           "-o", "f.bdg"])

    def f4(s):
        return s - log10(1 + s * log(10))

    got = rows_of(tmp_path / "f.bdg")
    exp = [("chr1", 0, 100, f4(5)), ("chr1", 100, 150, f4(7)),
           ("chr1", 150, 300, f4(5)), ("chr2", 0, 100, f4(6)),
           ("chr2", 100, 200, f4(3))]
    assert [r[:3] for r in got] == [r[:3] for r in exp]
    for g, x in zip(got, exp):
        assert g[3] == pytest.approx(x[3], abs=2e-5)
    assert lines(tmp_path / "f.bdg")[0] == track("fisher")


def test_fisher_of_zero_scores_is_zero(macs3_argparser, write_bedgraph,
                                       tmp_path):
    """-log10 p = 0 in every replicate (p = 1) combines to 0."""
    a = write_bedgraph([("chr1", 0, 100, 0)], name="a.bdg")
    b = write_bedgraph([("chr1", 0, 100, 0)], name="b.bdg")
    call(macs3_argparser, ["-i", a, b, "-m", "fisher", "--outdir", tmp_path,
                           "-o", "f.bdg"])
    assert lines(tmp_path / "f.bdg")[1:] == ["chr1\t0\t100\t0.00000"]


def test_default_method_is_fisher(macs3_argparser, reps, tmp_path):
    """Without -m the method is fisher."""
    call(macs3_argparser, ["-i", reps[0], reps[1], "--outdir", tmp_path,
                           "-o", "a.bdg"])
    call(macs3_argparser, ["-i", reps[0], reps[1], "-m", "fisher",
                           "--outdir", tmp_path, "-o", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")
    assert lines(tmp_path / "a.bdg")[0] == track("fisher")


@pytest.mark.parametrize("method, func", [
    ("max", max), ("mean", lambda v: sum(v) / len(v)), ("fisher", fisher)])
def test_same_replicate_twice(macs3_argparser, reps, tmp_path, method,
                              func):
    """-i R1 R1: max and mean give R1 back; fisher gives the 4-df score
    of 2v."""
    call(macs3_argparser, ["-i", reps[0], reps[0], "-m", method,
                           "--outdir", tmp_path, "-o", "s.bdg"])
    got = rows_of(tmp_path / "s.bdg")
    exp = combine([R1, R1], func)
    assert [r[:3] for r in got] == [r[:3] for r in exp]
    for g, x in zip(got, exp):
        assert g[3] == pytest.approx(x[3], abs=2e-5)


def test_long_option_names(macs3_argparser, reps, tmp_path):
    """--method/--ofile are the same options as -m/-o."""
    call(macs3_argparser, ["-i", reps[0], reps[1], "-m", "mean", "--outdir",
                           tmp_path, "-o", "a.bdg"])
    call(macs3_argparser, ["-i", reps[0], reps[1], "--method", "mean",
                           "--outdir", tmp_path, "--ofile", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")


def test_empty_replicate(macs3_argparser, reps, tmp_path):
    """An empty replicate shares no chromosome: track line only."""
    empty = tmp_path / "empty.bdg"
    empty.write_text("")
    call(macs3_argparser, ["-i", reps[0], empty, "-m", "max", "--outdir",
                           tmp_path, "-o", "e.bdg"])
    assert lines(tmp_path / "e.bdg") == [track("max")]


def test_no_common_chromosome(macs3_argparser, write_bedgraph, tmp_path):
    """Replicates without a common chromosome give the track line only."""
    a = write_bedgraph([("chrA", 0, 100, 1)], name="a.bdg")
    b = write_bedgraph([("chrB", 0, 100, 1)], name="b.bdg")
    call(macs3_argparser, ["-i", a, b, "-m", "max", "--outdir", tmp_path,
                           "-o", "e.bdg"])
    assert lines(tmp_path / "e.bdg") == [track("max")]


def test_output_in_outdir(macs3_argparser, reps, tmp_path):
    """-o is written inside --outdir and nothing else is written there."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", reps[0], reps[1], "-m", "max", "--outdir",
                           outdir, "-o", "m.bdg"])
    assert os.listdir(outdir) == ["m.bdg"]


# ------------------------------------
# option validation and the command line
# ------------------------------------

def test_cli_needs_two_replicates(run_macs3, parse_log, reps, tmp_path):
    """One -i file: opt_validate_cmbreps logs an error and exits 1."""
    res = run_macs3(["cmbreps", "-i", reps[0], "-m", "max", "--outdir",
                     tmp_path, "-o", "c.bdg"], timeout=60)
    assert res.returncode == 1
    assert parse_log(res.stderr) == [
        ("ERROR", "Combining replicates needs at least two replicates!")]
    assert not (tmp_path / "c.bdg").exists()


def test_invalid_method_in_process(macs3_argparser, reps, tmp_path, caplog):
    """A method outside the list (only reachable without argparse's
    choices check) is rejected by opt_validate_cmbreps: exit 1."""
    options = macs3_argparser.parse_args(["cmbreps", "-i", reps[0], reps[1],
                                          "--outdir", str(tmp_path), "-o",
                                          "c.bdg"])
    options.method = "median"
    with caplog.at_level(logging.INFO):
        with pytest.raises(SystemExit) as exc:
            cmbreps_run(options)
    assert exc.value.code == 1
    assert caplog.records[-1].getMessage().endswith("Invalid method: median")
    assert not (tmp_path / "c.bdg").exists()


@pytest.mark.parametrize("n", [2, 3])
def test_cli_log_messages(run_macs3, parse_log, reps, tmp_path, n):
    """A successful run logs one 'Read file' line per replicate and these
    INFO messages."""
    res = run_macs3(["cmbreps", "-i"] + reps[:n]
                    + ["-m", "mean", "--outdir", tmp_path, "-o", "c.bdg"],
                    timeout=60)
    assert res.returncode == 0, res.stderr
    out = os.path.join(str(tmp_path), "c.bdg")
    assert parse_log(res.stderr) == (
        [("INFO", "Read and build bedGraph for each replicate...")]
        + [("INFO", "Read file #%d" % (k + 1)) for k in range(n)]
        + [("INFO", "combining tracks 1-%d with method 'mean'" % n),
           ("INFO", "Write bedGraph of combined scores..."),
           ("INFO", "Finished 'mean'! Please check '%s'!" % out)])


@pytest.mark.parametrize("verbose, n_info", [("0", 0), ("1", 0), ("2", 6),
                                             ("3", 6)])
def test_cli_verbose(run_macs3, parse_log, reps, tmp_path, verbose, n_info):
    """--verbose 0/1 hide INFO messages; 2 and 3 show them."""
    res = run_macs3(["cmbreps", "-i", reps[0], reps[1], "-m", "max",
                     "--verbose", verbose, "--outdir", tmp_path, "-o",
                     "c.bdg"], timeout=60)
    assert res.returncode == 0
    assert [lv for lv, _ in parse_log(res.stderr)] == ["INFO"] * n_info
    assert (tmp_path / "c.bdg").exists()


def test_cli_outdir_is_created(run_macs3, reps, tmp_path):
    """A missing --outdir is created."""
    outdir = tmp_path / "p" / "q"
    res = run_macs3(["cmbreps", "-i", reps[0], reps[1], "-m", "max",
                     "--outdir", outdir, "-o", "c.bdg"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert os.listdir(outdir) == ["c.bdg"]


def test_cli_missing_replicate_file(run_macs3, reps, tmp_path):
    """A nonexistent replicate fails with FileNotFoundError (exit 1)."""
    missing = tmp_path / "nope.bdg"
    res = run_macs3(["cmbreps", "-i", reps[0], missing, "-m", "max",
                     "--outdir", tmp_path, "-o", "c.bdg"], timeout=60)
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr
    assert str(missing) in res.stderr
    assert not (tmp_path / "c.bdg").exists()


@pytest.mark.parametrize("argv, message", [
    (["-m", "max", "-o", "a"], "the following arguments are required: -i"),
    (["-i", "-m", "max", "-o", "a"], "argument -i: expected at least one "
                                     "argument"),
    (["-i", "a.bdg", "b.bdg", "-m", "max"],
     "the following arguments are required: -o/--ofile"),
    (["-i", "a.bdg", "b.bdg", "-m", "median", "-o", "a"],
     "argument -m/--method: invalid choice: 'median'"),
    (["-i", "a.bdg", "b.bdg", "-m", "max", "-o", "a", "--verbose", "x"],
     "argument --verbose: invalid int value: 'x'"),
])
def test_cli_argparse_errors(run_macs3, argv, message):
    """Bad or missing arguments: exit 2, usage and the argparse error."""
    res = run_macs3(["cmbreps"] + argv, timeout=60)
    assert res.returncode == 2
    assert res.stderr.startswith("usage: macs3 cmbreps")
    assert res.stderr.rstrip().splitlines()[-1].startswith(
        "macs3 cmbreps: error: " + message)


# ------------------------------------
# realistic run on the CTCF example
# ------------------------------------

@pytest.mark.parametrize("method", ["max", "mean", "fisher"])
def test_ctcf_tracks_match_standard(macs3_argparser, test_dir, tmp_path,
                                    method):
    """cmbreps on the CTCF chr22 treatment pileup, control lambda and
    bdgcmp ppois track (as in MACS3's command-line test) reproduces the
    upstream bedGraph.

    Pins the current output. The reference files are MACS3's own
    standard results over ~400k intervals (max, mean and fisher are
    checked against hand values and scipy's chi2 above).
    """
    narrow = test_dir / "standard_results_callpeak_narrow"
    files = [narrow / "run_callpeak_narrow0_treat_pileup.bdg",
             narrow / "run_callpeak_narrow0_control_lambda.bdg",
             test_dir / "standard_results_bdgcmp" / "run_bdgcmp_ppois.bdg"]
    name = "run_cmbreps_%s.bdg" % method
    call(macs3_argparser, ["-i"] + files + ["-m", method, "--outdir",
                                            tmp_path, "-o", name])
    assert (lines(tmp_path / name)
            == lines(test_dir / "standard_results_cmbreps" / name))
