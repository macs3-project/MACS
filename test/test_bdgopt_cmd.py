#!/usr/bin/env python

"""Module Description: Test the bdgopt subcommand
(MACS3/Commands/bdgopt_cmd.py) and its option checks in
MACS3/Utilities/OptValidator.py (opt_validate_bdgopt): operations on
the score column of a bedGraph (multiply, add, max, min, p2q).

Output content is checked in-process (the real argparse parser plus
``run``) on tiny bedGraphs with values worked out by hand; exit codes,
logs and argparse errors are checked with ``macs3`` in a subprocess.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import logging
import os
from pathlib import Path

import pytest

from MACS3.Commands.bdgopt_cmd import (run as bdgopt_run)

# ------------------------------------
# helpers
# ------------------------------------

# Input O. chr2 starts at 50, so [0,50) is filled with the baseline 0
# when the file is read.
O_ROWS = [("chr1", 0, 100, 1.5), ("chr1", 100, 200, -2),
          ("chr1", 200, 300, 4),
          ("chr2", 50, 150, 3),
          ("chr3", 0, 100, 0)]
# The intervals as read (including chr2's leading baseline block).
O_IV = [("chr1", 0, 100), ("chr1", 100, 200), ("chr1", 200, 300),
        ("chr2", 0, 50), ("chr2", 50, 150), ("chr3", 0, 100)]
O_VALS = [1.5, -2.0, 4.0, 0.0, 3.0, 0.0]


@pytest.fixture(autouse=True)
def _restore_optvalidator_level():
    """opt_validate_bdgopt sets the level of the OptValidator logger from
    --verbose; restore it so in-process runs leave no trace."""
    lg = logging.getLogger("MACS3.Utilities.OptValidator")
    level = lg.level
    yield
    lg.setLevel(level)


def call(argparser, argv):
    """Parse ``macs3 bdgopt <argv>`` with the real parser and run it
    in-process. Returns the options namespace."""
    options = argparser.parse_args(["bdgopt"] + [str(a) for a in argv])
    bdgopt_run(options)
    return options


def lines(path):
    """All lines of a text file without line endings."""
    return Path(path).read_text().splitlines()


def track(method):
    m = method.upper()
    return ('track type=bedGraph name="%s_modified_scores" '
            'description="Scores calculated by %s" visibility=2 '
            'alwaysZero=on' % (m, m))


def bdg_lines(iv, vals):
    return ["%s\t%d\t%d\t%.5f" % (c, s, e, v) for (c, s, e), v in
            zip(iv, vals)]


# ------------------------------------
# multiply, add, max, min
# ------------------------------------

@pytest.mark.parametrize("method, param, func", [
    ("multiply", "2", lambda x: x * 2),
    ("multiply", "0.5", lambda x: x * 0.5),
    ("multiply", "-1", lambda x: x * -1.0),
    ("multiply", "0", lambda x: 0.0 * x),
    ("add", "-1.5", lambda x: x - 1.5),
    ("add", "10", lambda x: x + 10),
    ("max", "1", lambda x: max(x, 1)),
    ("max", "-3", lambda x: max(x, -3)),
    ("min", "1", lambda x: min(x, 1)),
    ("min", "0", lambda x: min(x, 0)),
])
def test_score_operation(macs3_argparser, write_bedgraph, tmp_path, method,
                         param, func):
    """Each value v becomes v*p, v+p, max(v, p) or min(v, p); the
    intervals (including chr2's leading baseline block) are unchanged and
    a track line named after the method comes first."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    call(macs3_argparser, ["-i", bdg, "-m", method, "-p", param, "--outdir",
                           tmp_path, "-o", "out.bdg"])
    assert lines(tmp_path / "out.bdg") == (
        [track(method)] + bdg_lines(O_IV, [func(v) for v in O_VALS]))


def test_equal_neighbours_not_merged(macs3_argparser, write_bedgraph,
                                     tmp_path):
    """max -p 5 makes chr1 all 5: the three intervals stay three lines."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    call(macs3_argparser, ["-i", bdg, "-m", "max", "-p", "5", "--outdir",
                           tmp_path, "-o", "out.bdg"])
    assert lines(tmp_path / "out.bdg")[1:4] == [
        "chr1\t0\t100\t5.00000", "chr1\t100\t200\t5.00000",
        "chr1\t200\t300\t5.00000"]


def test_extra_params_beyond_first_ignored(macs3_argparser, write_bedgraph,
                                           tmp_path):
    """-p takes any number of values; only the first is used."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    call(macs3_argparser, ["-i", bdg, "-m", "multiply", "-p", "2",
                           "--outdir", tmp_path, "-o", "a.bdg"])
    call(macs3_argparser, ["-i", bdg, "-m", "multiply", "-p", "2", "3", "4",
                           "--outdir", tmp_path, "-o", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")


def test_float32_storage(macs3_argparser, write_bedgraph, tmp_path):
    """Values are stored as float32: 0.1 * 3 is written from the float32
    product, 0.30000001 -> 0.30000; 1e-6 * 3 rounds to 0.00000."""
    bdg = write_bedgraph([("chr1", 0, 10, 0.1), ("chr1", 10, 20, 1e-6)],
                         name="f.bdg")
    call(macs3_argparser, ["-i", bdg, "-m", "multiply", "-p", "3",
                           "--outdir", tmp_path, "-o", "out.bdg"])
    assert lines(tmp_path / "out.bdg")[1:] == [
        "chr1\t0\t10\t0.30000", "chr1\t10\t20\t0.00000"]


def test_float32_overflow_to_inf(macs3_argparser, write_bedgraph,
                                 tmp_path):
    """Scores are float32: 3e38 * 10 overflows to inf, -3e38 * 10 to
    -inf; values near the int32 coordinate limit are kept."""
    bdg = write_bedgraph([("chr1", 0, 100, 3e38), ("chr1", 100, 200, -3e38),
                          ("chr1", 200, 2147483000, 1)], name="big.bdg")
    call(macs3_argparser, ["-i", bdg, "-m", "multiply", "-p", "10",
                           "--outdir", tmp_path, "-o", "out.bdg"])
    assert lines(tmp_path / "out.bdg")[1:] == [
        "chr1\t0\t100\tinf", "chr1\t100\t200\t-inf",
        "chr1\t200\t2147483000\t10.00000"]


def test_long_option_names(macs3_argparser, write_bedgraph, tmp_path):
    """--ifile/--method/--extra-param/--ofile are the same options as
    -i/-m/-p/-o."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    call(macs3_argparser, ["-i", bdg, "-m", "add", "-p", "2.5", "--outdir",
                           tmp_path, "-o", "a.bdg"])
    call(macs3_argparser, ["--ifile", bdg, "--method", "add",
                           "--extra-param", "2.5", "--outdir", tmp_path,
                           "--ofile", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")
    assert lines(tmp_path / "a.bdg")[1:] == bdg_lines(
        O_IV, [v + 2.5 for v in O_VALS])


def test_output_in_outdir_with_sorted_chromosomes(macs3_argparser,
                                                  write_bedgraph, tmp_path):
    """-o is written inside --outdir; chromosomes come out sorted."""
    bdg = write_bedgraph([("chrZ", 0, 10, 1), ("chr10", 0, 10, 2),
                          ("chr2", 0, 10, 3)], name="z.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "-m", "add", "-p", "1", "--outdir",
                           outdir, "-o", "r.bdg"])
    assert os.listdir(outdir) == ["r.bdg"]
    assert lines(outdir / "r.bdg")[1:] == [
        "chr10\t0\t10\t3.00000", "chr2\t0\t10\t4.00000",
        "chrZ\t0\t10\t2.00000"]


def test_empty_input(macs3_argparser, tmp_path):
    """An empty bedGraph gives the track line only."""
    bdg = tmp_path / "empty.bdg"
    bdg.write_text("")
    call(macs3_argparser, ["-i", bdg, "-m", "multiply", "-p", "2",
                           "--outdir", tmp_path, "-o", "out.bdg"])
    assert lines(tmp_path / "out.bdg") == [track("multiply")]


# ------------------------------------
# p2q
# ------------------------------------

def test_p2q_single_enriched_value(macs3_argparser, write_bedgraph,
                                   tmp_path):
    """p2q (the default -m) with N = 1000 bp: the top p-score 10 has rank
    1, q = 10 + log10(1) - log10(1000) = 7; 0.5 and 0 give q <= 0 -> 0,
    and the two zero blocks on chr1 merge. -p is not needed."""
    bdg = write_bedgraph([("chr1", 0, 100, 10), ("chr1", 100, 200, 0.5),
                          ("chr1", 200, 500, 0), ("chr2", 0, 500, 0)],
                         name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "q.bdg"])
    assert lines(tmp_path / "q.bdg") == [
        track("p2q"), "chr1\t0\t100\t7.00000", "chr1\t100\t500\t0.00000",
        "chr2\t0\t500\t0.00000"]


def test_p2q_explicit_method_ignores_param(macs3_argparser, write_bedgraph,
                                           tmp_path):
    """-m p2q -p 3 gives the same output as the default p2q."""
    bdg = write_bedgraph([("chr1", 0, 100, 10), ("chr1", 100, 1000, 0)],
                         name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "a.bdg"])
    call(macs3_argparser, ["-i", bdg, "-m", "p2q", "-p", "3", "--outdir",
                           tmp_path, "-o", "b.bdg"])
    assert lines(tmp_path / "a.bdg") == lines(tmp_path / "b.bdg")


# ------------------------------------
# option validation and the command line
# ------------------------------------

@pytest.mark.parametrize("method, param", [
    ("multiply", []), ("add", []), ("multiply", ["-p"]), ("add", ["-p"])])
def test_cli_missing_extra_param(run_macs3, parse_log, write_bedgraph,
                                 tmp_path, method, param):
    """multiply and add without a -p value: opt_validate_bdgopt logs an
    error and exits 1 before writing anything."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    res = run_macs3(["bdgopt", "-i", bdg, "-m", method] + param
                    + ["--outdir", tmp_path, "-o", "out.bdg"], timeout=60)
    assert res.returncode == 1
    assert parse_log(res.stderr) == [
        ("ERROR", "Need EXTRAPARAM for method multiply or add!")]
    assert not (tmp_path / "out.bdg").exists()


def test_invalid_method_in_process(macs3_argparser, write_bedgraph,
                                   tmp_path, caplog):
    """A method outside the list (only reachable without argparse's
    choices check) is rejected by opt_validate_bdgopt: exit 1."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    options = macs3_argparser.parse_args(["bdgopt", "-i", bdg, "-m", "add",
                                          "-p", "1", "--outdir",
                                          str(tmp_path), "-o", "out.bdg"])
    options.method = "divide"
    with caplog.at_level(logging.INFO):
        with pytest.raises(SystemExit) as exc:
            bdgopt_run(options)
    assert exc.value.code == 1
    assert caplog.records[-1].getMessage().endswith("Invalid method: divide")
    assert not (tmp_path / "out.bdg").exists()


def test_cli_log_messages(run_macs3, parse_log, write_bedgraph, tmp_path):
    """A successful run exits 0 and logs exactly these INFO messages."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    res = run_macs3(["bdgopt", "-i", bdg, "-m", "add", "-p", "1",
                     "--outdir", tmp_path, "-o", "out.bdg"], timeout=60)
    assert res.returncode == 0, res.stderr
    out = os.path.join(str(tmp_path), "out.bdg")
    assert parse_log(res.stderr) == [
        ("INFO", "Read and build bedGraph..."),
        ("INFO", "Modify bedGraph..."),
        ("INFO", "Write bedGraph of modified scores..."),
        ("INFO", "Finished 'add'! Please check '%s'!" % out)]


@pytest.mark.parametrize("verbose, n_info", [("0", 0), ("1", 0), ("2", 4),
                                             ("3", 4)])
def test_cli_verbose(run_macs3, parse_log, write_bedgraph, tmp_path,
                     verbose, n_info):
    """--verbose 0/1 hide INFO messages; 2 and 3 show them."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    res = run_macs3(["bdgopt", "-i", bdg, "-m", "min", "-p", "1",
                     "--verbose", verbose, "--outdir", tmp_path, "-o",
                     "out.bdg"], timeout=60)
    assert res.returncode == 0
    assert [lv for lv, _ in parse_log(res.stderr)] == ["INFO"] * n_info
    assert (tmp_path / "out.bdg").exists()


def test_cli_outdir_is_created(run_macs3, write_bedgraph, tmp_path):
    """A missing --outdir is created."""
    bdg = write_bedgraph(O_ROWS, name="o.bdg")
    outdir = tmp_path / "n" / "m"
    res = run_macs3(["bdgopt", "-i", bdg, "-m", "add", "-p", "1",
                     "--outdir", outdir, "-o", "out.bdg"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert os.listdir(outdir) == ["out.bdg"]


def test_cli_missing_input_file(run_macs3, tmp_path):
    """A nonexistent -i fails with FileNotFoundError (exit 1)."""
    missing = tmp_path / "nope.bdg"
    res = run_macs3(["bdgopt", "-i", missing, "-m", "add", "-p", "1",
                     "--outdir", tmp_path, "-o", "out.bdg"], timeout=60)
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr
    assert str(missing) in res.stderr
    assert not (tmp_path / "out.bdg").exists()


@pytest.mark.parametrize("argv, message", [
    (["-m", "add", "-p", "1", "-o", "a"],
     "the following arguments are required: -i/--ifile"),
    (["-i", "x.bdg", "-m", "add", "-p", "1"],
     "the following arguments are required: -o/--ofile"),
    (["-i", "x.bdg", "-m", "divide", "-o", "a"],
     "argument -m/--method: invalid choice: 'divide'"),
    (["-i", "x.bdg", "-m", "MULTIPLY", "-p", "2", "-o", "a"],
     "argument -m/--method: invalid choice: 'MULTIPLY'"),
    (["-i", "x.bdg", "-m", "add", "-p", "one", "-o", "a"],
     "argument -p/--extra-param: invalid float value: 'one'"),
    (["-i", "x.bdg", "-m", "add", "-p", "1", "-o", "a", "--verbose", "v"],
     "argument --verbose: invalid int value: 'v'"),
])
def test_cli_argparse_errors(run_macs3, argv, message):
    """Bad or missing arguments: exit 2, usage and the argparse error."""
    res = run_macs3(["bdgopt"] + argv, timeout=60)
    assert res.returncode == 2
    assert res.stderr.startswith("usage: macs3 bdgopt")
    assert res.stderr.rstrip().splitlines()[-1].startswith(
        "macs3 bdgopt: error: " + message)


# ------------------------------------
# realistic run on the CTCF example
# ------------------------------------

@pytest.mark.parametrize("method, param, std", [
    ("min", "10", "run_bdgopt_min.bdg"), ("max", "2", "run_bdgopt_max.bdg")])
def test_ctcf_treat_pileup_matches_standard(macs3_argparser, test_dir,
                                            tmp_path, method, param, std):
    """bdgopt -m min -p 10 / -m max -p 2 on the CTCF chr22 callpeak
    treatment pileup reproduces the upstream bedGraphs.

    Pins the current output. The reference files are MACS3's own
    standard results (min/max themselves are checked by hand above).
    """
    pileup = (test_dir / "standard_results_callpeak_narrow"
              / "run_callpeak_narrow0_treat_pileup.bdg")
    call(macs3_argparser, ["-i", pileup, "-m", method, "-p", param,
                           "--outdir", tmp_path, "-o", std])
    assert (lines(tmp_path / std)
            == lines(test_dir / "standard_results_bdgopt" / std))
