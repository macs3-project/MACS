#!/usr/bin/env python

"""Module Description: Test the bdgbroadcall subcommand
(MACS3/Commands/bdgbroadcall_cmd.py): two-level (nested) broad peak
calling on a bedGraph score track, written as gappedPeak lines.

Output content is checked in-process (the real argparse parser plus
``run``) on a tiny bedGraph whose expected broad regions are derived
by hand; the command line surface is checked with ``macs3`` in a
subprocess.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from MACS3.Commands.bdgbroadcall_cmd import (run as bdgbroadcall_run)

# ------------------------------------
# helpers
# ------------------------------------

# Input B (continuous from 0). Defaults: -c 2 (level 1), -C 1 (level
# 2), -l 200, -g 30 (level 1 max gap), -G 800 (level 2 max gap).
# chr1 structure 1: cores [100,400)=3 and [500,800)=4 joined by a
#                   linking block [400,500)=1.5.
# chr1 structure 2: core [2100,2400)=2.5 inside links [2000,2100) and
#                   [2400,2500) at 1.5 (1200 bp after structure 1).
# chr1 structure 3: [4000,4250)=1.5, a level-2 region without a core.
# chr3: core [0,300)=5 at the chromosome start.
B_ROWS = [("chr1", 0, 100, 0), ("chr1", 100, 400, 3),
          ("chr1", 400, 500, 1.5), ("chr1", 500, 800, 4),
          ("chr1", 800, 2000, 0), ("chr1", 2000, 2100, 1.5),
          ("chr1", 2100, 2400, 2.5), ("chr1", 2400, 2500, 1.5),
          ("chr1", 2500, 4000, 0), ("chr1", 4000, 4250, 1.5),
          ("chr1", 4250, 5000, 0),
          ("chr3", 0, 300, 5), ("chr3", 300, 1000, 0)]

TRACK = ('track name="peak" description="peak" type=gappedPeak '
         'nextItemButton=on')


def call(argparser, argv):
    """Parse ``macs3 bdgbroadcall <argv>`` with the real parser and run
    it in-process. Returns the options namespace."""
    options = argparser.parse_args(["bdgbroadcall"] + [str(a) for a in argv])
    bdgbroadcall_run(options)
    return options


def lines(path):
    """All lines of a text file without line endings."""
    return Path(path).read_text().splitlines()


def gp_line(chrom, start, end, name, score10, nblocks, sizes, starts):
    """A gappedPeak line as bdgbroadcall writes it: score int(10*max of
    the level-2 region), strand '.', thickStart/thickEnd/itemRgb 0,
    blocks, then fold change/-log10p/-log10q 0."""
    return ("%s\t%d\t%d\t%s\t%d\t.\t0\t0\t0\t%d\t%s\t%s\t0\t0\t0"
            % (chrom, start, end, name, score10, nblocks, sizes, starts))


def named(prefix, rows):
    """gappedPeak lines from (chrom, start, end, score10, n, sizes,
    starts) rows, numbered <prefix>_broadRegion1.. across chromosomes."""
    return [gp_line(r[0], r[1], r[2], "%s_broadRegion%d" % (prefix, i + 1),
                    *r[3:]) for i, r in enumerate(rows)]


# Expected rows for the defaults, derived by hand:
# structure 1: level-2 region [100,800) (all >= 1, no gap), max 4 ->
#   score 40; level-1 cores [100,400) and [500,800) are 100 bp apart
#   (> 30) so they are two blocks: sizes 300,300 at offsets 0,400.
#   Both region ends are covered by a core, so no 1 bp end blocks.
# structure 2: level-2 [2000,2500), max 2.5 -> 25; one core
#   [2100,2400) plus 1 bp blocks at both ends: sizes 1,300,1, offsets
#   0,100,499.
# structure 3: level-2 [4000,4250), max 1.5 -> 15, no core: two 1 bp
#   blocks at offsets 0 and 249.
# chr3: level-2 = level-1 = [0,300), max 5 -> 50, one block.
S1 = ("chr1", 100, 800, 40, 2, "300,300", "0,400")
S2 = ("chr1", 2000, 2500, 25, 3, "1,300,1", "0,100,499")
S3 = ("chr1", 4000, 4250, 15, 2, "1,1", "0,249")
S4 = ("chr3", 0, 300, 50, 1, "300", "0")


# ------------------------------------
# broad peak calling: output content
# ------------------------------------

def test_default_output_exact(macs3_argparser, write_bedgraph, tmp_path):
    """Defaults on input B (derivation above the expected rows)."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "--o-prefix", "B"])
    out = tmp_path / "B_c2.0_C1.00_l200_g30_G800_broad.bed12"
    assert lines(out) == [TRACK] + named("B", [S1, S2, S3, S4])


@pytest.mark.parametrize("argv, fname", [
    ([], "B_c2.0_C1.00_l200_g30_G800_broad.bed12"),
    (["-c", "3.5", "-C", "1.6", "-l", "260", "-g", "100", "-G", "1200"],
     "B_c3.5_C1.60_l260_g100_G1200_broad.bed12"),
    (["-c", "2.25", "-C", "0.125"], "B_c2.2_C0.12_l200_g30_G800_broad.bed12"),
])
def test_prefix_file_name(macs3_argparser, write_bedgraph, tmp_path, argv,
                          fname):
    """--o-prefix names the file PREFIX_c<%.1f>_C<%.2f>_l<minlen>
    _g<lvl1 gap>_G<lvl2 gap>_broad.bed12 inside --outdir, and nothing
    else is written."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", outdir, "--o-prefix", "B"]
         + argv)
    assert sorted(os.listdir(outdir)) == [fname]


def test_ofile_names_file_and_regions(macs3_argparser, write_bedgraph,
                                      tmp_path):
    """-o gives the exact file name and the region name prefix
    <ofile>_broadRegion<n>; the track line does not change."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", outdir, "-o", "br.gp"])
    assert sorted(os.listdir(outdir)) == ["br.gp"]
    assert lines(outdir / "br.gp") == [TRACK] + named("br.gp",
                                                      [S1, S2, S3, S4])


@pytest.mark.parametrize("outopt", [["--o-prefix", "B"], ["-o", "B.gp"]])
def test_no_trackline(macs3_argparser, write_bedgraph, tmp_path, outopt):
    """--no-trackline drops only the track line."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    d1 = tmp_path / "with"
    d2 = tmp_path / "without"
    d1.mkdir()
    d2.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", d1] + outopt)
    call(macs3_argparser, ["-i", bdg, "--outdir", d2, "--no-trackline"]
         + outopt)
    (f1,) = os.listdir(d1)
    (f2,) = os.listdir(d2)
    assert f1 == f2
    assert lines(d1 / f1)[0] == TRACK
    assert lines(d2 / f2) == lines(d1 / f1)[1:]


@pytest.mark.parametrize("argv, expected", [
    # -c 3.5: only [500,800)=4 and chr3 are cores. Structure 1 keeps one
    # core, so a 1 bp block is added on the left only (sizes 1,300 at
    # 0,400); structure 2 has no core any more (1 bp blocks at 0, 499);
    # it is still reported because chr1 has a core elsewhere.
    (["-c", "3.5"],
     [("chr1", 100, 800, 40, 2, "1,300", "0,400"),
      ("chr1", 2000, 2500, 25, 2, "1,1", "0,499"), S3, S4]),
    # -c 2.5 at equality: [2100,2400)=2.5 is still a core
    (["-c", "2.5"],
     [("chr1", 100, 800, 40, 2, "300,300", "0,400"), S2, S3, S4]),
    # -c 4: [500,800)=4 at equality is the only chr1 core
    (["-c", "4"],
     [("chr1", 100, 800, 40, 2, "1,300", "0,400"),
      ("chr1", 2000, 2500, 25, 2, "1,1", "0,499"), S3, S4]),
])
def test_cutoff_peak(macs3_argparser, write_bedgraph, tmp_path, argv,
                     expected):
    """-c sets the level-1 (core) cutoff; equality counts as enriched."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp",
                           "--no-trackline"] + argv)
    assert lines(tmp_path / "o.gp") == named("o.gp", expected)


@pytest.mark.parametrize("argv, expected", [
    # -C 1.5 at equality keeps the 1.5 linking blocks: same as default
    (["-C", "1.5"], [S1, S2, S3, S4]),
    # -C 1.6 drops the 1.5 blocks: structure 1 is still one level-2
    # region (cores 100 bp apart <= -G 800); structure 2 shrinks to its
    # core [2100,2400) (one 300 bp block); structure 3 disappears.
    (["-C", "1.6"],
     [S1, ("chr1", 2100, 2400, 25, 1, "300", "0"), S4]),
])
def test_cutoff_link(macs3_argparser, write_bedgraph, tmp_path, argv,
                     expected):
    """-C sets the level-2 (linking) cutoff."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp",
                           "--no-trackline"] + argv)
    assert lines(tmp_path / "o.gp") == named("o.gp", expected)


@pytest.mark.parametrize("minlen, expected", [
    ("250", [S1, S2, S3, S4]),     # structure 3 is exactly 250 bp
    ("251", [S1, S2, S4]),
    ("300", [S1, S2, S4]),         # every core is exactly 300 bp
])
def test_min_length(macs3_argparser, write_bedgraph, tmp_path, minlen,
                    expected):
    """-l applies to level-1 and level-2 regions; equality is kept."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp",
                           "--no-trackline", "-l", minlen])
    assert lines(tmp_path / "o.gp") == named("o.gp", expected)


@pytest.mark.parametrize("gap, expected", [
    ("99", [S1, S2, S3, S4]),
    # -g 100: the two cores of structure 1 merge into one 700 bp core
    ("100", [("chr1", 100, 800, 40, 1, "700", "0"), S2, S3, S4]),
])
def test_lvl1_max_gap(macs3_argparser, write_bedgraph, tmp_path, gap,
                      expected):
    """-g merges level-1 cores whose gap is <= the value."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp",
                           "--no-trackline", "-g", gap])
    assert lines(tmp_path / "o.gp") == named("o.gp", expected)


@pytest.mark.parametrize("gap, expected", [
    ("1199", [S1, S2, S3, S4]),
    # -G 1200: structures 1 and 2 (1200 bp apart) form one level-2 region
    # [100,2500) with three cores; it starts on a core, so only a right
    # 1 bp block is added: sizes 300,300,300,1 at 0,400,2000,2399.
    # Structure 3 is 1500 bp further and stays separate.
    ("1200", [("chr1", 100, 2500, 40, 4, "300,300,300,1",
               "0,400,2000,2399"), S3, S4]),
])
def test_lvl2_max_gap(macs3_argparser, write_bedgraph, tmp_path, gap,
                      expected):
    """-G merges level-2 regions whose gap is <= the value."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp",
                           "--no-trackline", "-G", gap])
    assert lines(tmp_path / "o.gp") == named("o.gp", expected)


def test_long_option_names(macs3_argparser, write_bedgraph, tmp_path):
    """--ifile/--cutoff-peak/--cutoff-link/--min-length/--lvl1-max-gap/
    --lvl2-max-gap/--ofile are the same options as the short forms."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    call(macs3_argparser, ["-i", bdg, "-c", "3.5", "-C", "1.6", "-l", "250",
                           "-g", "100", "-G", "1200", "--outdir", tmp_path,
                           "-o", "a.gp", "--no-trackline"])
    call(macs3_argparser, ["--ifile", bdg, "--cutoff-peak", "3.5",
                           "--cutoff-link", "1.6", "--min-length", "250",
                           "--lvl1-max-gap", "100", "--lvl2-max-gap", "1200",
                           "--outdir", tmp_path, "--ofile", "b.gp",
                           "--no-trackline"])
    a = lines(tmp_path / "a.gp")
    assert a == [x.replace("b.gp", "a.gp") for x in lines(tmp_path / "b.gp")]
    # -c 3.5 -C 1.6 -g 100 -G 1200: one level-2 region [100,800) with the
    # single core [500,800) (the 2.5 block is no longer a core and the
    # 1.6 link cutoff drops the 1.5 blocks), plus chr3.
    assert a == named("a.gp", [("chr1", 100, 800, 40, 2, "1,300", "0,400"),
                               ("chr1", 2100, 2400, 25, 2, "1,1", "0,299"),
                               S4])


def test_empty_input(macs3_argparser, tmp_path):
    """An empty bedGraph gives the track line only."""
    bdg = tmp_path / "empty.bdg"
    bdg.write_text("")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp"])
    assert lines(tmp_path / "o.gp") == [TRACK]


def test_no_enrichment(macs3_argparser, write_bedgraph, tmp_path):
    """A track below both cutoffs gives no regions."""
    bdg = write_bedgraph([("chr1", 0, 1000, 0.5)], name="z.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.gp",
                           "--no-trackline"])
    assert lines(tmp_path / "o.gp") == []


# ------------------------------------
# option checks inside call_broadpeaks
# ------------------------------------

@pytest.mark.parametrize("argv, message", [
    (["-c", "1", "-C", "2"], "level 1 cutoff should be larger than level 2."),
    (["-c", "2", "-C", "2"], "level 1 cutoff should be larger than level 2."),
    (["-g", "800", "-G", "800"],
     "level 2 maximum gap should be larger than level 1."),
    (["-g", "900", "-G", "800"],
     "level 2 maximum gap should be larger than level 1."),
])
def test_cli_cutoff_and_gap_order_checks(run_macs3, write_bedgraph,
                                         tmp_path, argv, message):
    """-c must exceed -C and -G must exceed -g: an AssertionError ends the
    run with exit 1 and no output file is left."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    res = run_macs3(["bdgbroadcall", "-i", bdg, "--outdir", tmp_path,
                     "-o", "o.gp"] + argv, timeout=60)
    assert res.returncode == 1
    assert res.stderr.rstrip().splitlines()[-1] == "AssertionError: " + message
    assert not (tmp_path / "o.gp").exists()


# ------------------------------------
# command line: logs, exit codes
# ------------------------------------

def test_cli_log_messages(run_macs3, parse_log, write_bedgraph, tmp_path):
    """A successful run exits 0 and logs exactly these INFO messages."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    res = run_macs3(["bdgbroadcall", "-i", bdg, "--outdir", tmp_path,
                     "--o-prefix", "B"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert parse_log(res.stderr) == [
        ("INFO", "Read and build bedGraph..."),
        ("INFO", "Call peaks from bedGraph..."),
        ("INFO", "Write peaks..."),
        ("INFO", "Done")]
    assert lines(tmp_path / "B_c2.0_C1.00_l200_g30_G800_broad.bed12") == (
        [TRACK] + named("B", [S1, S2, S3, S4]))


def test_cli_outdir_is_created(run_macs3, write_bedgraph, tmp_path):
    """A missing --outdir is created."""
    bdg = write_bedgraph(B_ROWS, name="b.bdg")
    outdir = tmp_path / "a" / "b"
    res = run_macs3(["bdgbroadcall", "-i", bdg, "--outdir", outdir,
                     "-o", "o.gp"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert os.listdir(outdir) == ["o.gp"]


def test_cli_missing_input_file(run_macs3, tmp_path):
    """A nonexistent -i fails with FileNotFoundError (exit 1)."""
    missing = tmp_path / "nope.bdg"
    res = run_macs3(["bdgbroadcall", "-i", missing, "--outdir", tmp_path,
                     "-o", "o.gp"], timeout=60)
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr
    assert str(missing) in res.stderr
    assert not (tmp_path / "o.gp").exists()


@pytest.mark.parametrize("argv, message", [
    (["-o", "a"], "the following arguments are required: -i/--ifile"),
    (["-i", "x.bdg"],
     "one of the arguments -o/--ofile --o-prefix is required"),
    (["-i", "x.bdg", "--o-prefix", "P", "-o", "a"],
     "argument -o/--ofile: not allowed with argument --o-prefix"),
    (["-i", "x.bdg", "-o", "a", "-c", "high"],
     "argument -c/--cutoff-peak: invalid float value: 'high'"),
    (["-i", "x.bdg", "-o", "a", "-C", "low"],
     "argument -C/--cutoff-link: invalid float value: 'low'"),
    (["-i", "x.bdg", "-o", "a", "-l", "2e2"],
     "argument -l/--min-length: invalid int value: '2e2'"),
    (["-i", "x.bdg", "-o", "a", "-g", "3.0"],
     "argument -g/--lvl1-max-gap: invalid int value: '3.0'"),
    (["-i", "x.bdg", "-o", "a", "-G", "far"],
     "argument -G/--lvl2-max-gap: invalid int value: 'far'"),
])
def test_cli_argparse_errors(run_macs3, argv, message):
    """Bad or missing arguments: exit 2, usage and the argparse error."""
    res = run_macs3(["bdgbroadcall"] + argv, timeout=60)
    assert res.returncode == 2
    assert res.stderr.startswith("usage: macs3 bdgbroadcall")
    assert res.stderr.rstrip().splitlines()[-1] == (
        "macs3 bdgbroadcall: error: " + message)


# ------------------------------------
# realistic run on the CTCF example
# ------------------------------------

def test_ctcf_fe_track_matches_standard(macs3_argparser, test_dir,
                                        tmp_path):
    """bdgbroadcall -c 2 -C 1.5 on the CTCF chr22 FE track reproduces the
    upstream gappedPeak file in every column except the score; the
    score column is int(10 * the largest FE inside each region).

    Pins the current output. The reference file is MACS3's own
    standard result; broad regions over 192k intervals cannot be derived
    by hand. The score column of the stored standard file is not
    correct: it was written when BroadPeakIO.add declared ``score`` as
    a C long, so the maximum FE was truncated to an integer before the
    x10 (every stored score is a multiple of 10; 55 of 6709 regions
    differ, e.g. 30 where the maximum FE is 3.5). The current code
    keeps the float and writes int(10 * max FE), as bdgpeakcall does
    for narrowPeak, so the score is checked against the FE track
    directly instead.
    """
    fe = test_dir / "standard_results_bdgcmp" / "run_bdgcmp_FE.bdg"
    std = (test_dir / "standard_results_bdgbroadcall"
           / "run_bdgbroadcall_w_prefix_c2.0_C1.50_l200_g30_G800_broad.bed12")
    call(macs3_argparser, ["-i", fe, "-c", "2", "-C", "1.5", "--outdir",
                           tmp_path, "--o-prefix",
                           "run_bdgbroadcall_w_prefix"])
    got = lines(tmp_path / std.name)
    exp = lines(std)
    assert got[0] == exp[0] == TRACK
    got_rows = [x.split("\t") for x in got[1:]]
    exp_rows = [x.split("\t") for x in exp[1:]]
    assert ([r[:4] + r[5:] for r in got_rows]
            == [r[:4] + r[5:] for r in exp_rows])
    # score = int(10 * max FE over the intervals inside the region)
    fe_rows = [x.split("\t") for x in lines(fe)]
    starts = np.array([int(r[1]) for r in fe_rows])
    ends = np.array([int(r[2]) for r in fe_rows])
    vals = np.array([float(r[3]) for r in fe_rows], dtype="float32")
    for r in got_rows:
        s, e = int(r[1]), int(r[2])
        i = np.searchsorted(starts, s)
        j = np.searchsorted(ends, e, side="right")
        assert int(r[4]) == int(10 * float(vals[i:j].max())), r
