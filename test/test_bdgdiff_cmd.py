#!/usr/bin/env python

"""Module Description: Test the bdgdiff subcommand
(MACS3/Commands/bdgdiff_cmd.py): differential regions from treatment
and control bedGraphs of two conditions, scored by log10 likelihood
ratios.

The inputs are tiny tracks on a 100 bp grid; a reference
implementation in the test computes the three log10 likelihood ratio
tracks (t1 vs c1, t2 vs c2, t1 vs t2, pseudocount 0.01, depth scaling
from --d1/--d2), the three categories, the merging (-g), the length
filter (-l) and the length-weighted mean score. Output content is
checked in-process (the real argparse parser plus ``run``); exit
codes, logs and argparse errors are checked with ``macs3`` in a
subprocess.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import os
from math import log
from pathlib import Path

import numpy as np
import pytest

from MACS3.Commands.bdgdiff_cmd import (run as bdgdiff_run)

# ------------------------------------
# helpers
# ------------------------------------

LOG10_E = 0.43429448190325176
PSEUDOCOUNT = 0.01      # TwoConditionScores' default, used by bdgdiff
BIN = 100
CATS = ("cond1", "cond2", "common")


def alternating(n):
    """Control values 1.0, 1.25, 1.0, ... so every bin edge is a
    breakpoint of the combined tracks."""
    return [1.0 if i % 2 == 0 else 1.25 for i in range(n)]


def blocks(n, spans, high, low=1.0):
    v = [low] * n
    for a, b in spans:
        for i in range(a, b):
            v[i] = high
    return v


# Per-bin values. chr1: t1 high in bins 5-10 (condition 1 only), both
# high in bins 15-20 (common), t2 high in bins 24-27 (condition 2
# only). chr2: t1 moderately high in bins 3-5 and 7-9 with one
# background bin between them.
BINS = {
    "chr1": {"t1": blocks(30, [(5, 11), (15, 21)], 40.0),
             "t2": blocks(30, [(15, 21), (24, 28)], 40.0),
             "c1": alternating(30), "c2": alternating(30)},
    "chr2": {"t1": blocks(15, [(3, 6), (7, 10)], 20.0),
             "t2": [1.0] * 15,
             "c1": alternating(15), "c2": alternating(15)},
}


def track_rows(key, bins=BINS):
    rows = []
    for chrom in sorted(bins):
        for i, v in enumerate(bins[chrom][key]):
            rows.append((chrom, i * BIN, (i + 1) * BIN, v))
    return rows


@pytest.fixture
def four(write_bedgraph):
    """Paths of the t1, c1, t2, c2 bedGraphs built from BINS."""
    return {k: write_bedgraph(track_rows(k), name="%s.bdg" % k)
            for k in ("t1", "c1", "t2", "c2")}


def argv4(f):
    return ["--t1", f["t1"], "--c1", f["c1"], "--t2", f["t2"], "--c2",
            f["c2"]]


def call(argparser, argv):
    """Parse ``macs3 bdgdiff <argv>`` with the real parser and run it
    in-process. Returns the options namespace."""
    options = argparser.parse_args(["bdgdiff"] + [str(a) for a in argv])
    bdgdiff_run(options)
    return options


def lines(path):
    """All lines of a text file without line endings."""
    return Path(path).read_text().splitlines()


TRACKS = {
    "cond1": 'track name="condition 1 (peaks)" description="unique '
             'regions in condition 1" visibility=1',
    "cond2": 'track name="condition 2 (peaks)" description="unique '
             'regions in condition 2" visibility=1',
    "common": 'track name="common (peaks)" description="common regions '
              'in both conditions" visibility=1',
}


def logLR_asym(x, y):
    if x > y:
        return (x * (log(x) - log(y)) + y - x) * LOG10_E
    if x < y:
        return (x * (-log(x) + log(y)) - y + x) * LOG10_E
    return 0.0


def logLR_sym(x, y):
    if x > y:
        return (x * (log(x) - log(y)) + y - x) * LOG10_E
    if y > x:
        return (y * (log(x) - log(y)) + y - x) * LOG10_E
    return 0.0


def depth_factors(d1, d2):
    """--d1/--d2: the deeper condition is scaled down to the other."""
    if d1 > d2:
        return d2 / d1, 1.0
    if d1 < d2:
        return 1.0, d1 / d2
    return 1.0, 1.0


def reference(cutoff=3.0, minlen=200, maxgap=100, d1=1.0, d2=1.0,
              bins=BINS):
    """Expected regions per category: {cat: [(chrom, start, end,
    score)]}.

    Per bin: a = logLR(t1 vs c1), b = logLR(t2 vs c2), s = symmetric
    logLR(t1 vs t2), each value plus 0.01 and times its condition's
    depth factor (float32 as in MACS3). cond1: a >= C and s >= C;
    cond2: b >= C and s <= -C; common: a >= C, b >= C and |s| <= C.
    Bins of a category are merged when the gap is <= maxgap; regions
    shorter than minlen are dropped; the score is the length-weighted
    mean of s, -s or |s|.
    """
    f32 = np.float32
    f1, f2 = depth_factors(d1, d2)
    out = {k: [] for k in CATS}
    for chrom in sorted(bins):
        b_ = bins[chrom]
        n = len(b_["t1"])

        def adj(v, f):
            return float(f32(f32(f32(v) + f32(PSEUDOCOUNT)) * f32(f)))

        a = [float(f32(logLR_asym(adj(b_["t1"][i], f1), adj(b_["c1"][i], f1))))
             for i in range(n)]
        b = [float(f32(logLR_asym(adj(b_["t2"][i], f2), adj(b_["c2"][i], f2))))
             for i in range(n)]
        s = [float(f32(logLR_sym(adj(b_["t1"][i], f1), adj(b_["t2"][i], f2))))
             for i in range(n)]
        sel = {"cond1": [a[i] >= cutoff and s[i] >= cutoff
                         for i in range(n)],
               "cond2": [b[i] >= cutoff and s[i] <= -cutoff
                         for i in range(n)],
               "common": [a[i] >= cutoff and b[i] >= cutoff
                          and -cutoff <= s[i] <= cutoff for i in range(n)]}
        score = {"cond1": s, "cond2": [-x for x in s],
                 "common": [abs(x) for x in s]}
        for cat in CATS:
            cur = []
            for i in range(n):
                if not sel[cat][i]:
                    continue
                piece = (i * BIN, (i + 1) * BIN, score[cat][i])
                if cur and piece[0] - cur[-1][1] > maxgap:
                    if cur[-1][1] - cur[0][0] >= minlen:
                        out[cat].append(close(chrom, cur))
                    cur = []
                cur.append(piece)
            if cur and cur[-1][1] - cur[0][0] >= minlen:
                out[cat].append(close(chrom, cur))
    return out


def close(chrom, pieces):
    tot = sum(e - s for s, e, _ in pieces)
    mean = sum(v * (e - s) for s, e, v in pieces) / tot
    return (chrom, pieces[0][0], pieces[-1][1], mean)


def regions(path):
    """(chrom, start, end, name, score) rows after the track line."""
    return [tuple(x.split("\t")) for x in lines(path)[1:]]


def lengths(rows):
    """(chrom, length) of each region: the part of the output that does
    not depend on where regions are placed."""
    return [(r[0], int(r[2]) - int(r[1])) for r in rows]


def prefix_files(outdir, prefix, cutoff):
    return {cat: Path(outdir) / ("%s_c%.1f_%s.bed" % (prefix, cutoff, cat))
            for cat in CATS}


# ------------------------------------
# output files and their layout
# ------------------------------------

def test_prefix_files_tracklines_and_names(macs3_argparser, four, tmp_path):
    """--o-prefix D writes D_c3.0_cond1.bed, D_c3.0_cond2.bed and
    D_c3.0_common.bed, each with its track line and regions named
    D_<cat>_<n>. Defaults give: cond1 chr1 600 bp and chr2 700 bp (two
    300 bp blocks 100 bp apart, merged by -g 100), cond2 chr1 400 bp,
    common chr1 600 bp (its score, |t1 vs t2|, is 0)."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, argv4(four) + ["--outdir", outdir, "--o-prefix",
                                         "D"])
    files = prefix_files(outdir, "D", 3.0)
    assert sorted(os.listdir(outdir)) == sorted(f.name for f in
                                                files.values())
    exp_len = {"cond1": [("chr1", 600), ("chr2", 700)],
               "cond2": [("chr1", 400)], "common": [("chr1", 600)]}
    for cat in CATS:
        assert lines(files[cat])[0] == TRACKS[cat]
        rows = regions(files[cat])
        assert [r[3] for r in rows] == ["D_%s_%d" % (cat, i + 1)
                                        for i in range(len(rows))]
        assert lengths(rows) == exp_len[cat]
    assert [r[4] for r in regions(files["common"])] == ["0"]


def test_ofile_three_names(macs3_argparser, four, tmp_path):
    """-o A B C writes the three files in the order cond1, cond2, common
    and uses each file name (as given) as the region name prefix."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, argv4(four) + ["--outdir", tmp_path, "--o-prefix",
                                         "D"])
    call(macs3_argparser, argv4(four) + ["--outdir", outdir, "-o", "u1.bed",
                                         "u2.bed", "both.bed"])
    assert sorted(os.listdir(outdir)) == ["both.bed", "u1.bed", "u2.bed"]
    ref = prefix_files(tmp_path, "D", 3.0)
    for cat, name in zip(CATS, ("u1.bed", "u2.bed", "both.bed")):
        got = lines(outdir / name)
        assert got[0] == TRACKS[cat]
        rows = regions(outdir / name)
        assert [r[3] for r in rows] == ["%s%d" % (name, i + 1)
                                        for i in range(len(rows))]
        # same regions and scores as the --o-prefix run
        assert ([r[:3] + r[4:] for r in rows]
                == [r[:3] + r[4:] for r in regions(ref[cat])])


@pytest.mark.parametrize("cutoff, fname", [
    ("3", "D_c3.0_cond1.bed"), ("10", "D_c10.0_cond1.bed"),
    ("2.25", "D_c2.2_cond1.bed")])
def test_cutoff_in_file_names(macs3_argparser, four, tmp_path, cutoff,
                              fname):
    """The -C value is written with one decimal in the file names."""
    call(macs3_argparser, argv4(four) + ["-C", cutoff, "--outdir", tmp_path,
                                         "--o-prefix", "D"])
    assert (tmp_path / fname).exists()
    assert (tmp_path / fname.replace("cond1", "cond2")).exists()
    assert (tmp_path / fname.replace("cond1", "common")).exists()


# ------------------------------------
# -C, -l, -g, --d1/--d2 against the reference
# ------------------------------------

def assert_lengths_match(outdir, prefix, cutoff, ref):
    files = prefix_files(outdir, prefix, cutoff)
    for cat in CATS:
        exp = [(c, e - s) for c, s, e, _ in ref[cat]]
        assert lengths(regions(files[cat])) == exp, cat


@pytest.mark.parametrize("cutoff", [3.0, 10.0, 20.0, 46.0, 47.0])
def test_cutoff(macs3_argparser, four, tmp_path, cutoff):
    """-C applies to all three likelihood ratios. chr2's condition-1
    blocks (logLR about 16-18) go above 20; at 46 only bins with control
    1.0 pass (every other bin, merged across the 100 bp gaps); at 47
    nothing passes."""
    call(macs3_argparser, argv4(four) + ["-C", cutoff, "--outdir", tmp_path,
                                         "--o-prefix", "D"])
    assert_lengths_match(tmp_path, "D", cutoff, reference(cutoff=cutoff))


@pytest.mark.parametrize("minlen, n_regions", [
    ("300", 4), ("400", 4), ("401", 3), ("700", 1), ("701", 0)])
def test_min_len(macs3_argparser, four, tmp_path, minlen, n_regions):
    """-l drops regions shorter than it (equality is kept)."""
    call(macs3_argparser, argv4(four) + ["-l", minlen, "--outdir", tmp_path,
                                         "--o-prefix", "D"])
    ref = reference(minlen=int(minlen))
    assert sum(len(v) for v in ref.values()) == n_regions
    assert_lengths_match(tmp_path, "D", 3.0, ref)


@pytest.mark.parametrize("maxgap, chr2_lengths", [
    ("100", [700]), ("99", [300, 300]), ("0", [300, 300]),
    ("199", [700])])
def test_max_gap(macs3_argparser, four, tmp_path, maxgap, chr2_lengths):
    """-g merges regions of a category whose gap is <= it: the two chr2
    condition-1 blocks are 100 bp apart."""
    call(macs3_argparser, argv4(four) + ["-g", maxgap, "--outdir", tmp_path,
                                         "--o-prefix", "D"])
    ref = reference(maxgap=int(maxgap))
    assert [e - s for c, s, e, _ in ref["cond1"] if c == "chr2"] == \
        chr2_lengths
    assert_lengths_match(tmp_path, "D", 3.0, ref)


@pytest.mark.parametrize("d1, d2, counts", [
    ("1", "1", (2, 1, 1)),
    ("2", "2", (2, 1, 1)),
    # condition 1 halved: the common region (40 vs 40 -> 20 vs 40,
    # symmetric logLR -3.36) becomes condition-2 specific
    ("2", "1", (2, 2, 0)),
    # condition 2 halved: it becomes condition-1 specific (+3.36)
    ("1", "2", (3, 1, 0)),
    ("10", "20", (3, 1, 0)),
    # 2/3 scaling is not enough (logLR -1.25): it stays common
    ("1.5", "1", (2, 1, 1)),
])
def test_depths(macs3_argparser, four, tmp_path, d1, d2, counts):
    """--d1/--d2 scale the deeper condition down to the other one."""
    call(macs3_argparser, argv4(four) + ["--d1", d1, "--d2", d2, "--outdir",
                                         tmp_path, "--o-prefix", "D"])
    ref = reference(d1=float(d1), d2=float(d2))
    assert tuple(len(ref[c]) for c in CATS) == counts
    assert_lengths_match(tmp_path, "D", 3.0, ref)


def test_depth_long_option_names(macs3_argparser, four, tmp_path):
    """--depth1/--depth2 are the same options as --d1/--d2."""
    call(macs3_argparser, argv4(four) + ["--d1", "2", "--d2", "1",
                                         "--outdir", tmp_path,
                                         "--o-prefix", "A"])
    call(macs3_argparser, argv4(four) + ["--depth1", "2", "--depth2", "1",
                                         "--outdir", tmp_path,
                                         "--o-prefix", "B"])
    fa = prefix_files(tmp_path, "A", 3.0)
    fb = prefix_files(tmp_path, "B", 3.0)
    for cat in CATS:
        assert ([r[:3] + r[4:] for r in regions(fa[cat])]
                == [r[:3] + r[4:] for r in regions(fb[cat])])


def test_long_option_names(macs3_argparser, four, tmp_path):
    """--cutoff/--min-len/--max-gap/--ofile are the same options as
    -C/-l/-g/-o."""
    call(macs3_argparser, argv4(four) + ["-C", "10", "-l", "300", "-g", "50",
                                         "--outdir", tmp_path, "-o", "a1",
                                         "a2", "a3"])
    call(macs3_argparser, argv4(four) + ["--cutoff", "10", "--min-len",
                                         "300", "--max-gap", "50",
                                         "--outdir", tmp_path, "--ofile",
                                         "b1", "b2", "b3"])
    for a, b in (("a1", "b1"), ("a2", "b2"), ("a3", "b3")):
        ra = regions(tmp_path / a)
        rb = regions(tmp_path / b)
        assert [r[:3] + r[4:] for r in ra] == [r[:3] + r[4:] for r in rb]
    # -g 50 splits chr2's condition-1 region into two 300 bp regions
    assert lengths(regions(tmp_path / "a1")) == [
        ("chr1", 600), ("chr2", 300), ("chr2", 300)]


def test_no_common_chromosome(macs3_argparser, write_bedgraph, four,
                              tmp_path):
    """When one track has none of the others' chromosomes, all three
    files hold only their track line."""
    other = write_bedgraph([("chrZ", 0, 1000, 5)], name="z.bdg")
    f = dict(four)
    f["c2"] = other
    call(macs3_argparser, argv4(f) + ["--outdir", tmp_path, "--o-prefix",
                                      "N"])
    for cat, path in prefix_files(tmp_path, "N", 3.0).items():
        assert lines(path) == [TRACKS[cat]]


def test_chromosome_missing_from_one_track_is_skipped(macs3_argparser,
                                                      write_bedgraph,
                                                      tmp_path):
    """chr2 absent from c1: only chr1 regions are reported."""
    f = {k: write_bedgraph(track_rows(k), name="%s.bdg" % k)
         for k in ("t1", "t2", "c2")}
    f["c1"] = write_bedgraph([r for r in track_rows("c1") if r[0] == "chr1"],
                             name="c1.bdg")
    call(macs3_argparser, argv4(f) + ["--outdir", tmp_path, "--o-prefix",
                                      "S"])
    files = prefix_files(tmp_path, "S", 3.0)
    assert lengths(regions(files["cond1"])) == [("chr1", 600)]


def test_region_scores(macs3_argparser, four, tmp_path):
    """The score is the length-weighted mean logLR: about 46.99 for the
    chr1 condition-1 region (t1 40 vs t2 1 in every bin), not 46.

    Regression test: the score truncated each log10 likelihood ratio to an
    integer (cython.long accumulator). Fixed upstream in 9dbb3c7 (#739,
    issue #715)."""
    call(macs3_argparser, argv4(four) + ["--outdir", tmp_path, "--o-prefix",
                                         "D"])
    ref = reference()
    for cat, path in prefix_files(tmp_path, "D", 3.0).items():
        got = [float(r[4]) for r in regions(path)]
        assert got == pytest.approx([x[3] for x in ref[cat]], abs=1e-4), cat


# ------------------------------------
# command line: checks, logs, exit codes
# ------------------------------------

def test_cli_log_messages(run_macs3, parse_log, four, tmp_path):
    """A successful run exits 0 and logs exactly these INFO messages."""
    res = run_macs3(["bdgdiff"] + argv4(four) + ["--outdir", tmp_path,
                                                 "--o-prefix", "D"],
                    timeout=60)
    assert res.returncode == 0, res.stderr
    assert parse_log(res.stderr) == [
        ("INFO", "Read and build treatment 1 bedGraph..."),
        ("INFO", "Read and build control 1 bedGraph..."),
        ("INFO", "Read and build treatment 2 bedGraph..."),
        ("INFO", "Read and build control 2 bedGraph..."),
        ("INFO", "Write peaks..."),
        ("INFO", "Done")]


@pytest.mark.parametrize("gap, minlen", [("200", "200"), ("300", "200")])
def test_cli_maxgap_not_below_minlen_is_logged(run_macs3, parse_log, four,
                                               tmp_path, gap, minlen):
    """-g >= -l logs a CRITICAL message first."""
    res = run_macs3(["bdgdiff"] + argv4(four) + ["-g", gap, "-l", minlen,
                                                 "--outdir", tmp_path,
                                                 "--o-prefix", "D"],
                    timeout=60)
    assert parse_log(res.stderr)[0] == (
        "CRITICAL", "MAXGAP should be smaller than MINLEN! Your input is "
                    "MAXGAP = %s and MINLEN = %s" % (gap, minlen))


def test_cli_maxgap_not_below_minlen_still_runs(run_macs3, four, tmp_path):
    """-g >= -l is only advice: the run goes on, exits 0 and writes the
    three files."""
    # Surprising but documented as advice: the help says the maximum gap
    # "should be smaller" than the minimum length, and the CRITICAL
    # message says the same without stopping the run; the output is the
    # ordinary result for these -g/-l values.
    res = run_macs3(["bdgdiff"] + argv4(four) + ["-g", "200", "-l", "200",
                                                 "--outdir", tmp_path,
                                                 "--o-prefix", "D"],
                    timeout=60)
    assert res.returncode == 0, res.stderr
    for path in prefix_files(tmp_path, "D", 3.0).values():
        assert path.exists()


def test_cli_outdir_is_created(run_macs3, four, tmp_path):
    """A missing --outdir is created and holds the three files."""
    outdir = tmp_path / "d" / "e"
    res = run_macs3(["bdgdiff"] + argv4(four) + ["--outdir", outdir, "-o",
                                                 "a.bed", "b.bed", "c.bed"],
                    timeout=60)
    assert res.returncode == 0, res.stderr
    assert sorted(os.listdir(outdir)) == ["a.bed", "b.bed", "c.bed"]


@pytest.mark.parametrize("which", ["t1", "c1", "t2", "c2"])
def test_cli_missing_input_file(run_macs3, four, tmp_path, which):
    """Any nonexistent input fails with FileNotFoundError (exit 1) and no
    output file is written."""
    f = dict(four)
    f[which] = str(tmp_path / "nope.bdg")
    res = run_macs3(["bdgdiff"] + argv4(f) + ["--outdir", tmp_path,
                                              "--o-prefix", "M"],
                    timeout=60)
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr
    assert f[which] in res.stderr
    assert not prefix_files(tmp_path, "M", 3.0)["cond1"].exists()


@pytest.mark.parametrize("argv, message", [
    (["--c1", "c1", "--t2", "t2", "--c2", "c2", "--o-prefix", "D"],
     "the following arguments are required: --t1"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--o-prefix", "D"],
     "the following arguments are required: --c2"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2"],
     "one of the arguments --o-prefix -o/--ofile is required"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2", "-o", "a",
      "b"], "argument -o/--ofile: expected 3 arguments"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2", "--o-prefix",
      "D", "-o", "a", "b", "c"],
     "argument -o/--ofile: not allowed with argument --o-prefix"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2", "--o-prefix",
      "D", "-C", "x"], "argument -C/--cutoff: invalid float value: 'x'"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2", "--o-prefix",
      "D", "-l", "1.5"], "argument -l/--min-len: invalid int value: '1.5'"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2", "--o-prefix",
      "D", "-g", "g"], "argument -g/--max-gap: invalid int value: 'g'"),
    (["--t1", "t1", "--c1", "c1", "--t2", "t2", "--c2", "c2", "--o-prefix",
      "D", "--d1", "deep"],
     "argument --d1/--depth1: invalid float value: 'deep'"),
])
def test_cli_argparse_errors(run_macs3, argv, message):
    """Bad or missing arguments: exit 2, usage and the argparse error."""
    res = run_macs3(["bdgdiff"] + argv, timeout=60)
    assert res.returncode == 2
    assert res.stderr.startswith("usage: macs3 bdgdiff")
    assert res.stderr.rstrip().splitlines()[-1] == (
        "macs3 bdgdiff: error: " + message)


# ------------------------------------
# realistic run on the CTCF example
# ------------------------------------

CTCF_OUTPUTS = [
    (["--o-prefix", "run_bdgdiff_prefix"],
     ["run_bdgdiff_prefix_c3.0_cond1.bed", "run_bdgdiff_prefix_c3.0_cond2.bed",
      "run_bdgdiff_prefix_c3.0_common.bed"]),
    (["-o", "cond1.bed", "cond2.bed", "common.bed"],
     ["cond1.bed", "cond2.bed", "common.bed"]),
]

CTCF_TRACKS = {
    "t1": ("standard_results_callpeak_narrow",
           "run_callpeak_narrow0_treat_pileup.bdg"),
    "c1": ("standard_results_callpeak_narrow",
           "run_callpeak_narrow0_control_lambda.bdg"),
    "t2": ("standard_results_callpeak_narrow_revert",
           "run_callpeak_narrow_revert_treat_pileup.bdg"),
    "c2": ("standard_results_callpeak_narrow_revert",
           "run_callpeak_narrow_revert_control_lambda.bdg"),
}


def run_ctcf(argparser, test_dir, outdir, outopt):
    """bdgdiff on the CTCF chr22 callpeak tracks: ChIP vs control as
    condition 1, the reverted run (control vs ChIP) as condition 2."""
    argv = []
    for key in ("t1", "c1", "t2", "c2"):
        sub, name = CTCF_TRACKS[key]
        argv += ["--" + key, test_dir / sub / name]
    call(argparser, argv + ["--outdir", outdir] + outopt)


@pytest.mark.parametrize("outopt, names", CTCF_OUTPUTS)
def test_ctcf_matches_standard(macs3_argparser, test_dir, tmp_path, outopt,
                               names):
    """The three files have the upstream track lines, the region names
    follow <prefix><n>, and there is no condition-2 or common region
    (the independent reference below finds none either).

    The coordinates and scores of the stored standard files carry a
    one-interval shift (they were made with float64 score sums, so their
    scores are not truncated), so they are not compared here.
    """
    run_ctcf(macs3_argparser, test_dir, tmp_path, outopt)
    std_dir = test_dir / "standard_results_bdgdiff"
    for name, cat in zip(names, CATS):
        got = lines(tmp_path / name)
        assert got[0] == lines(std_dir / name)[0] == TRACKS[cat], name
    prefix = (outopt[1] + "_cond1_" if outopt[0] == "--o-prefix"
              else names[0])
    rows = regions(tmp_path / names[0])
    assert len(rows) > 600
    assert [r[3] for r in rows] == ["%s%d" % (prefix, i + 1)
                                    for i in range(len(rows))]
    assert len(lines(tmp_path / names[1])) == 1
    assert len(lines(tmp_path / names[2])) == 1


