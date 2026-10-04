#!/usr/bin/env python

"""Module Description: Test the bdgpeakcall subcommand
(MACS3/Commands/bdgpeakcall_cmd.py): naive peak calling on a bedGraph
score track and the --cutoff-analysis report.

Output content is checked in-process (the real argparse parser plus
``run``) on tiny bedGraphs whose expected peaks are derived by hand in
the comments; the command line surface (exit codes, argparse errors,
log messages, --outdir creation) is checked with ``macs3`` in a
subprocess.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from MACS3.Commands.bdgpeakcall_cmd import (run as bdgpeakcall_run)

# ------------------------------------
# helpers
# ------------------------------------

# Input P, continuous from 0 on every chromosome.
# chr1: [100,300) is exactly at the default cutoff 5; the gap
#       [300,330) is exactly the default max gap 30; [330,500) holds
#       the maximum 8.5.
# chr2: three blocks of the maximum 6 separated by 5s, so the region
#       above 5 is [0,330) without a gap.
# chr3: never reaches 5.
P_ROWS = [("chr1", 0, 100, 1), ("chr1", 100, 300, 5),
          ("chr1", 300, 330, 1), ("chr1", 330, 500, 8.5),
          ("chr1", 500, 1000, 1),
          ("chr2", 0, 100, 6), ("chr2", 100, 110, 5),
          ("chr2", 110, 210, 6), ("chr2", 210, 230, 5),
          ("chr2", 230, 330, 6), ("chr2", 330, 1000, 0),
          ("chr3", 0, 1000, 2)]


def call(argparser, argv):
    """Parse ``macs3 bdgpeakcall <argv>`` with the real parser and run it
    in-process. Returns the options namespace."""
    options = argparser.parse_args(["bdgpeakcall"] + [str(a) for a in argv])
    bdgpeakcall_run(options)
    return options


def lines(path):
    """All lines of a text file without line endings."""
    return Path(path).read_text().splitlines()


def np_line(chrom, start, end, name, score10, summit_offset):
    """A narrowPeak line as bdgpeakcall writes it: score is int(10*max),
    strand '.', fold change/-log10p/-log10q are 0, then the summit
    offset from the peak start."""
    return "%s\t%d\t%d\t%s\t%d\t.\t0\t0\t0\t%d" % (chrom, start, end, name,
                                                score10, summit_offset)


def np_track(name):
    return ('track type=narrowPeak name="%s" description="%s" '
            'nextItemButton=on' % (name, name))


def ref_peaks(rows, cutoff, minlen, maxgap, strict=False):
    """Reference peak caller on a continuous bedGraph.

    Regions with score >= cutoff (> cutoff when ``strict``) are merged
    when the gap between them is <= maxgap; merged regions shorter than
    minlen are dropped. Scores and the cutoff are compared as float32,
    the type MACS3 stores bedGraph values in. Returns a list of
    (chrom, start, end, [(s, e, v), ...]) in chromosome order.
    """
    cut = np.float32(cutoff)
    out = []
    for chrom in sorted({r[0] for r in rows}):
        blocks = [(s, e, np.float32(v)) for c, s, e, v in rows if c == chrom]
        above = [b for b in blocks if (b[2] > cut if strict else b[2] >= cut)]
        cur = []
        for b in above:
            if cur and b[0] - cur[-1][1] > maxgap:
                if cur[-1][1] - cur[0][0] >= minlen:
                    out.append((chrom, cur[0][0], cur[-1][1], cur))
                cur = []
            cur.append(b)
        if cur and cur[-1][1] - cur[0][0] >= minlen:
            out.append((chrom, cur[0][0], cur[-1][1], cur))
    return out


def ref_cutoff_table(rows, minlen, maxgap, max_score=100, steps=100):
    """Reference --cutoff-analysis report.

    Following the option help, the score range from the smallest score
    to min(largest score, max_score) is cut into ``steps`` intervals
    (cutoffs rounded to 3 decimals); for each cutoff, regions with score
    strictly above it (as in callpeak's cutoff analysis) are merged
    (gap <= maxgap) and kept when >= minlen. Rows run from the highest
    cutoff down and rows with no peak are left out. Inputs used with
    this reference have no enriched region starting at position 0 and
    no negative score.
    """
    vals = [float(np.float32(r[3])) for r in rows]
    minv = min(vals)
    maxv = min(max(vals), max_score)
    # the step is held as a float32
    s = float(np.float32((maxv - minv) / steps))
    cutoffs = [round(x, 3) for x in np.arange(minv, maxv, s)]
    out = ["score\tnpeaks\tlpeaks\tavelpeak"]
    for c in reversed(cutoffs):
        c32 = float(np.float32(c))
        peaks = ref_peaks(rows, c32, minlen, maxgap, strict=True)
        n = len(peaks)
        tot = sum(e - s0 for _, s0, e, _ in peaks)
        if n > 0:
            out.append("%.2f\t%d\t%d\t%.2f" % (c32, n, tot, tot / n))
    return out


# ------------------------------------
# peak calling: output file content
# ------------------------------------

def test_default_output_exact(macs3_argparser, write_bedgraph, tmp_path):
    """Defaults -c 5 -l 200 -g 30 on input P.

    chr1: [100,300) (=5, kept at equality) and [330,500) (8.5) are 30 bp
    apart (= max gap) -> one peak [100,500); summit is the midpoint of
    the 8.5 block, int((330+500)/2) = 415, offset 315; score
    int(10*8.5) = 85.
    chr2: [0,330) is all >= 5; the maximum 6 occurs in three blocks with
    midpoints 50, 160, 280 and the middle one (160) is the summit;
    score 60. chr3 never reaches 5. Peak numbers run on across
    chromosomes.
    """
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "--o-prefix", "P"])
    out = tmp_path / "P_c5.0_l200_g30_peaks.narrowPeak"
    assert lines(out) == [
        np_track("P"),
        np_line("chr1", 100, 500, "P_narrowPeak1", 85, 315),
        np_line("chr2", 0, 330, "P_narrowPeak2", 60, 160),
    ]


@pytest.mark.parametrize("argv, fname", [
    ([], "P_c5.0_l200_g30_peaks.narrowPeak"),
    (["-c", "7.5", "-l", "150", "-g", "10"],
     "P_c7.5_l150_g10_peaks.narrowPeak"),
    # the cutoff is printed with one decimal in the name
    (["-c", "5.04"], "P_c5.0_l200_g30_peaks.narrowPeak"),
    (["-c", "0", "-l", "0", "-g", "0"], "P_c0.0_l0_g0_peaks.narrowPeak"),
])
def test_prefix_file_name(macs3_argparser, write_bedgraph, tmp_path,
                          argv, fname):
    """--o-prefix names the file PREFIX_c<cutoff %.1f>_l<minlen>_g<maxgap>
    _peaks.narrowPeak inside --outdir, and nothing else is written."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", outdir, "--o-prefix", "P"]
         + argv)
    assert sorted(os.listdir(outdir)) == [fname]


def test_ofile_names_file_and_peaks(macs3_argparser, write_bedgraph,
                                    tmp_path):
    """-o gives the exact file name; the track name and the peak name
    prefix (<ofile>_narrowPeak<n>) are taken from it."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", outdir, "-o", "my.np"])
    assert sorted(os.listdir(outdir)) == ["my.np"]
    assert lines(outdir / "my.np") == [
        np_track("my.np"),
        np_line("chr1", 100, 500, "my.np_narrowPeak1", 85, 315),
        np_line("chr2", 0, 330, "my.np_narrowPeak2", 60, 160),
    ]


@pytest.mark.parametrize("outopt", [["--o-prefix", "P"], ["-o", "P.np"]])
def test_no_trackline(macs3_argparser, write_bedgraph, tmp_path, outopt):
    """--no-trackline drops only the first (track) line."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
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
    with_tl = lines(d1 / f1)
    assert with_tl[0].startswith("track type=narrowPeak")
    assert lines(d2 / f2) == with_tl[1:]


@pytest.mark.parametrize("cutoff, minlen, expected", [
    # 5 is kept at equality (see test_default_output_exact)
    ("5", "200", [("chr1", 100, 500, 85, 315), ("chr2", 0, 330, 60, 160)]),
    # just above 5: chr1 keeps only [330,500) (170 bp < 200); chr2 keeps
    # the three 6-blocks [0,100) [110,210) [230,330), gaps 10 and 20,
    # merged to [0,330) with the same middle summit 160
    ("5.01", "200", [("chr2", 0, 330, 60, 160)]),
    # 8.5 at equality with -l 100: [330,500), summit 415 (offset 85)
    ("8.5", "100", [("chr1", 330, 500, 85, 85)]),
    ("8.51", "100", []),
])
def test_cutoff_at_equality(macs3_argparser, write_bedgraph, tmp_path,
                            cutoff, minlen, expected):
    """A score equal to -c is enriched; anything below is not."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "-c", cutoff, "-l", minlen, "--no-trackline"])
    assert lines(tmp_path / "o.np") == [
        np_line(c, s, e, "o.np_narrowPeak%d" % (i + 1), sc, off)
        for i, (c, s, e, sc, off) in enumerate(expected)]


@pytest.mark.parametrize("cutoff, n_peaks", [("2.3", 1), ("2.31", 0)])
def test_cutoff_equality_in_float32(macs3_argparser, write_bedgraph,
                                    tmp_path, cutoff, n_peaks):
    """Scores and the cutoff are both float32, so a score of 2.3 equals
    -c 2.3 although 2.3 is not exact in binary."""
    bdg = write_bedgraph([("chr1", 0, 100, 0), ("chr1", 100, 400, 2.3),
                          ("chr1", 400, 600, 0)], name="f.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "-c", cutoff, "--no-trackline"])
    got = lines(tmp_path / "o.np")
    assert len(got) == n_peaks
    if n_peaks:
        # summit int((100+400)/2) = 250 -> offset 150; int(10*2.3f) = 22
        assert got == [np_line("chr1", 100, 400, "o.np_narrowPeak1", 22, 150)]


@pytest.mark.parametrize("minlen, expected_chroms", [
    ("330", ["chr1", "chr2"]),   # chr2 peak is exactly 330 bp
    ("331", ["chr1"]),
    ("400", ["chr1"]),           # chr1 peak is exactly 400 bp
    ("401", []),
])
def test_min_length(macs3_argparser, write_bedgraph, tmp_path, minlen,
                    expected_chroms):
    """-l keeps a peak whose length equals it and drops shorter ones."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "-l", minlen, "--no-trackline"])
    got = lines(tmp_path / "o.np")
    assert [x.split("\t")[0] for x in got] == expected_chroms


@pytest.mark.parametrize("maxgap, expected", [
    # the 30 bp gap on chr1 is merged at -g 30
    ("30", [("chr1", 100, 500, 85, 315), ("chr2", 0, 330, 60, 160)]),
    # -g 29: chr1 splits; [100,300) is 200 bp (kept, summit 200, score
    # 50); [330,500) is 170 bp (dropped). chr2 has no gap at all.
    ("29", [("chr1", 100, 300, 50, 100), ("chr2", 0, 330, 60, 160)]),
    ("0", [("chr1", 100, 300, 50, 100), ("chr2", 0, 330, 60, 160)]),
    # a large gap changes nothing more here
    ("1000", [("chr1", 100, 500, 85, 315), ("chr2", 0, 330, 60, 160)]),
])
def test_max_gap(macs3_argparser, write_bedgraph, tmp_path, maxgap,
                 expected):
    """-g merges enriched regions whose gap is <= the value."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "-g", maxgap, "--no-trackline"])
    assert lines(tmp_path / "o.np") == [
        np_line(c, s, e, "o.np_narrowPeak%d" % (i + 1), sc, off)
        for i, (c, s, e, sc, off) in enumerate(expected)]


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_summit_middle_of_equal_maxima(macs3_argparser, write_bedgraph,
                                       tmp_path, n_max):
    """With n blocks at the maximum, the summit is the middle block's
    midpoint, the left one of the two middles when n is even (index
    (n+1)//2 - 1)."""
    # blocks of 100 bp at 9, separated by 50 bp at 6, starting at 1000;
    # the k-th maximum block is [1000+150k, 1100+150k)
    rows = [("chr1", 0, 1000, 0)]
    pos = 1000
    mids = []
    for k in range(n_max):
        rows.append(("chr1", pos, pos + 100, 9))
        mids.append(pos + 50)
        pos += 100
        if k < n_max - 1:
            rows.append(("chr1", pos, pos + 50, 6))
            pos += 50
    rows.append(("chr1", pos, pos + 1000, 0))
    bdg = write_bedgraph(rows, name="m.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "-c", "5", "-l", "50", "--no-trackline"])
    summit = mids[(n_max + 1) // 2 - 1]
    assert lines(tmp_path / "o.np") == [
        np_line("chr1", 1000, pos, "o.np_narrowPeak1", 90, summit - 1000)]


def test_summit_midpoint_rounds_down(macs3_argparser, write_bedgraph,
                                     tmp_path):
    """The summit of block [101,202) is int(303/2) = 151; the peak is
    [51,252), so the offset is 100."""
    bdg = write_bedgraph([("chr1", 0, 51, 0), ("chr1", 51, 101, 6),
                          ("chr1", 101, 202, 7), ("chr1", 202, 252, 6),
                          ("chr1", 252, 400, 0)], name="r.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "-l", "100", "--no-trackline"])
    assert lines(tmp_path / "o.np") == [
        np_line("chr1", 51, 252, "o.np_narrowPeak1", 70, 100)]


@pytest.mark.parametrize("value, score10", [
    (8.5, 85), (6.25, 62), (12, 120), (1000.5, 10005)])
def test_score_column_is_int_of_ten_times_max(macs3_argparser,
                                              write_bedgraph, tmp_path,
                                              value, score10):
    """Column 5 is int(10 * the peak's maximum score) (truncated)."""
    bdg = write_bedgraph([("chr1", 0, 100, 0), ("chr1", 100, 400, value),
                          ("chr1", 400, 500, 0)], name="s.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "--no-trackline"])
    assert lines(tmp_path / "o.np") == [
        np_line("chr1", 100, 400, "o.np_narrowPeak1", score10, 150)]


def test_long_option_names(macs3_argparser, write_bedgraph, tmp_path):
    """--ifile/--cutoff/--min-length/--max-gap/--ofile are the same
    options as -i/-c/-l/-g/-o."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "-c", "5.01", "-l", "100", "-g", "10",
                           "--outdir", tmp_path, "-o", "a.np"])
    call(macs3_argparser, ["--ifile", bdg, "--cutoff", "5.01",
                           "--min-length", "100", "--max-gap", "10",
                           "--outdir", tmp_path, "--ofile", "a2.np"])
    a = lines(tmp_path / "a.np")
    assert a[1:] == [x.replace("a2.np", "a.np") for x in
                     lines(tmp_path / "a2.np")[1:]]
    # -g 10: chr2's 6-blocks are 10 and 20 bp apart -> [0,210) only
    # (the 100 bp block after the 20 bp gap stands alone and is kept)
    assert [tuple(x.split("\t")[:3]) for x in a[1:]] == [
        ("chr1", "330", "500"), ("chr2", "0", "210"), ("chr2", "230", "330")]


def test_call_summits_flag_has_no_effect(macs3_argparser, write_bedgraph,
                                         tmp_path):
    """--call-summits (help suppressed) is a reserved flag: the summit is
    always computed, so the output is identical with and without it."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "a.np"])
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "b.np",
                           "--call-summits"])
    a = lines(tmp_path / "a.np")
    b = lines(tmp_path / "b.np")
    assert a[1:] == [x.replace("b.np", "a.np") for x in b[1:]]
    assert len(a) == 3


def test_header_lines_are_skipped(macs3_argparser, write_bedgraph,
                                  tmp_path):
    """Lines starting with 'track', 'browser' or '#' are not data."""
    plain = write_bedgraph(P_ROWS, name="p.bdg")
    headed = write_bedgraph(["track type=bedGraph name=x",
                             "browser position chr1:1-100",
                             "# a comment"] + P_ROWS, name="h.bdg")
    call(macs3_argparser, ["-i", plain, "--outdir", tmp_path, "-o", "a.np"])
    call(macs3_argparser, ["-i", headed, "--outdir", tmp_path, "-o", "a2.np"])
    assert (lines(tmp_path / "a.np")[1:]
            == [x.replace("a2.np", "a.np") for x in
                lines(tmp_path / "a2.np")[1:]])


def test_empty_input_writes_trackline_only(macs3_argparser, tmp_path):
    """An empty bedGraph gives a narrowPeak file with only the track line,
    and nothing at all with --no-trackline."""
    bdg = tmp_path / "empty.bdg"
    bdg.write_text("")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "a.np"])
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "b.np",
                           "--no-trackline"])
    assert lines(tmp_path / "a.np") == [np_track("a.np")]
    assert (tmp_path / "b.np").read_text() == ""


def test_single_interval_chromosome(macs3_argparser, write_bedgraph,
                                    tmp_path):
    """A chromosome made of one interval above the cutoff is one peak
    covering it, with the summit at its midpoint."""
    bdg = write_bedgraph([("chrM", 0, 300, 7)], name="one.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "--no-trackline"])
    assert lines(tmp_path / "o.np") == [
        np_line("chrM", 0, 300, "o.np_narrowPeak1", 70, 150)]


def test_chromosomes_written_in_sorted_order(macs3_argparser,
                                             write_bedgraph, tmp_path):
    """Chromosomes are written in byte-sorted order, not input order."""
    rows = [("chrB", 0, 300, 7), ("chr10", 0, 300, 7), ("chr2", 0, 300, 7)]
    bdg = write_bedgraph(rows, name="o.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "o.np",
                           "--no-trackline"])
    assert [x.split("\t")[0] for x in lines(tmp_path / "o.np")] == [
        "chr10", "chr2", "chrB"]


def test_peaks_match_reference_caller(macs3_argparser, write_bedgraph,
                                      tmp_path):
    """Coordinates on a busier track match the reference caller for a
    sweep of settings (summits and scores are covered elsewhere)."""
    rng = np.random.default_rng(7)
    rows = []
    for chrom in ("chr1", "chr2", "chr3"):
        pos = 0
        for _ in range(60):
            w = int(rng.integers(5, 80))
            rows.append((chrom, pos, pos + w, int(rng.integers(0, 10))))
            pos += w
    bdg = write_bedgraph(rows, name="busy.bdg")
    for cutoff, minlen, maxgap in [(5, 50, 20), (3, 100, 0), (8, 10, 60)]:
        name = "o_%d_%d_%d.np" % (cutoff, minlen, maxgap)
        call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", name,
                               "-c", cutoff, "-l", minlen, "-g", maxgap,
                               "--no-trackline"])
        got = [tuple(x.split("\t")[:3]) for x in lines(tmp_path / name)]
        exp = [(c, str(s), str(e)) for c, s, e, _ in
               ref_peaks(rows, cutoff, minlen, maxgap)]
        assert got == exp


# ------------------------------------
# --cutoff-analysis
# ------------------------------------

# Input CA: every chromosome starts with a 0 block, and no score
# (0, 3, 7, 9) equals a cutoff of the grids used below.
CA_ROWS = [("chr1", 0, 100, 0), ("chr1", 100, 400, 3),
           ("chr1", 400, 420, 0), ("chr1", 420, 700, 7),
           ("chr1", 700, 1000, 0),
           ("chr2", 0, 200, 0), ("chr2", 200, 260, 9),
           ("chr2", 260, 600, 3), ("chr2", 600, 1000, 0),
           ("chr3", 0, 500, 0), ("chr3", 500, 560, 7),
           ("chr3", 560, 1000, 0)]


def test_cutoff_analysis_by_hand(macs3_argparser, write_bedgraph,
                                 tmp_path):
    """Derivation for -l 200 -g 30 --cutoff-analysis-steps 4 (max 100):
    range [0, 9], step 2.25, cutoffs 0, 2.25, 4.5, 6.75.
    6.75: chr1 [420,700) 280 bp; chr2 [200,260) 60 (dropped); chr3
          [500,560) 60 (dropped) -> 1 peak, 280 bp.
    4.5:  same as 6.75 -> 1, 280.
    2.25: chr1 [100,400)+[420,700) (gap 20) = [100,700) 600; chr2
          [200,600) 400; chr3 60 (dropped) -> 2, 1000, mean 500.
    0:    same as 2.25 (0 is not above 0) -> 2, 1000.
    """
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "ca.txt",
                           "--cutoff-analysis", "--cutoff-analysis-steps",
                           "4"])
    assert lines(tmp_path / "ca.txt") == [
        "score\tnpeaks\tlpeaks\tavelpeak",
        "6.75\t1\t280\t280.00",
        "4.50\t1\t280\t280.00",
        "2.25\t2\t1000\t500.00",
        "0.00\t2\t1000\t500.00",
    ]


def test_cutoff_analysis_default_name_and_only_output(macs3_argparser,
                                                      write_bedgraph,
                                                      tmp_path):
    """With --o-prefix the report is PREFIX_l<minlen>_g<maxgap>
    _cutoff_analysis.txt; no narrowPeak file is written; the default
    max (100) and steps (100) give the reference table."""
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", outdir, "--o-prefix", "Q",
                           "--cutoff-analysis"])
    assert sorted(os.listdir(outdir)) == ["Q_l200_g30_cutoff_analysis.txt"]
    assert (lines(outdir / "Q_l200_g30_cutoff_analysis.txt")
            == ref_cutoff_table(CA_ROWS, 200, 30))


def test_cutoff_analysis_ofile(macs3_argparser, write_bedgraph, tmp_path):
    """With -o the report goes to exactly that file name."""
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    outdir = tmp_path / "out"
    outdir.mkdir()
    call(macs3_argparser, ["-i", bdg, "--outdir", outdir, "-o", "rep.tsv",
                           "--cutoff-analysis"])
    assert sorted(os.listdir(outdir)) == ["rep.tsv"]
    assert lines(outdir / "rep.tsv") == ref_cutoff_table(CA_ROWS, 200, 30)


@pytest.mark.parametrize("steps, cmax", [
    (4, 100), (4, 8), (10, 5), (100, 3), (7, 100), (1, 100)])
def test_cutoff_analysis_steps_and_max(macs3_argparser, write_bedgraph,
                                       tmp_path, steps, cmax):
    """--cutoff-analysis-steps sets the number of cutoffs and
    --cutoff-analysis-max caps the top of the range."""
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "ca.txt",
                           "--cutoff-analysis",
                           "--cutoff-analysis-steps", steps,
                           "--cutoff-analysis-max", cmax])
    assert (lines(tmp_path / "ca.txt")
            == ref_cutoff_table(CA_ROWS, 200, 30, max_score=cmax,
                                steps=steps))


def test_cutoff_analysis_max_changes_grid(macs3_argparser, write_bedgraph,
                                          tmp_path):
    """--cutoff-analysis-max 8 with 4 steps: range [0, 8], cutoffs 0, 2,
    4, 6 (the 9 on chr2 is above the cap but still counted)."""
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "ca.txt",
                           "--cutoff-analysis", "--cutoff-analysis-steps",
                           "4", "--cutoff-analysis-max", "8"])
    assert [x.split("\t")[0] for x in lines(tmp_path / "ca.txt")] == [
        "score", "6.00", "4.00", "2.00", "0.00"]


@pytest.mark.parametrize("minlen, maxgap", [(50, 30), (50, 0), (300, 50)])
def test_cutoff_analysis_minlen_maxgap(macs3_argparser, write_bedgraph,
                                       tmp_path, minlen, maxgap):
    """-l and -g apply to every cutoff of the report and appear in the
    default file name."""
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "--o-prefix",
                           "Q", "--cutoff-analysis", "-l", minlen, "-g",
                           maxgap, "--cutoff-analysis-steps", "20"])
    f = tmp_path / ("Q_l%d_g%d_cutoff_analysis.txt" % (minlen, maxgap))
    assert lines(f) == ref_cutoff_table(CA_ROWS, minlen, maxgap, steps=20)


def test_cutoff_analysis_ignores_cutoff_option(macs3_argparser,
                                               write_bedgraph, tmp_path):
    """-c does not change the report."""
    bdg = write_bedgraph(CA_ROWS, name="ca.bdg")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "a.txt",
                           "--cutoff-analysis"])
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "b.txt",
                           "--cutoff-analysis", "-c", "8"])
    assert lines(tmp_path / "a.txt") == lines(tmp_path / "b.txt")


def test_cutoff_analysis_empty_input(macs3_argparser, tmp_path):
    """An empty bedGraph gives the header row only."""
    bdg = tmp_path / "empty.bdg"
    bdg.write_text("")
    call(macs3_argparser, ["-i", bdg, "--outdir", tmp_path, "-o", "ca.txt",
                           "--cutoff-analysis"])
    assert lines(tmp_path / "ca.txt") == ["score\tnpeaks\tlpeaks\tavelpeak"]


# ------------------------------------
# command line: logs, exit codes, files
# ------------------------------------

@pytest.mark.parametrize("extra, messages", [
    ([], ["Read and build bedGraph...", "Call peaks from bedGraph...",
          "Write peaks...", "Done"]),
    (["--cutoff-analysis"],
     ["Read and build bedGraph...",
      "Analyze cutoff vs number of peaks/total length of peaks/average "
      "length of peak", "Write report...", "Done"]),
])
def test_cli_log_messages(run_macs3, parse_log, write_bedgraph, tmp_path,
                          extra, messages):
    """A successful run exits 0 and logs exactly these INFO messages."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    res = run_macs3(["bdgpeakcall", "-i", bdg, "--outdir", tmp_path,
                     "--o-prefix", "P"] + extra, timeout=60)
    assert res.returncode == 0, res.stderr
    assert parse_log(res.stderr) == [("INFO", m) for m in messages]
    assert res.stdout == ""


def test_cli_outdir_is_created(run_macs3, write_bedgraph, tmp_path):
    """A missing --outdir (with missing parents) is created."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    outdir = tmp_path / "new" / "sub"
    res = run_macs3(["bdgpeakcall", "-i", bdg, "--outdir", outdir,
                     "--o-prefix", "P"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert sorted(os.listdir(outdir)) == ["P_c5.0_l200_g30_peaks.narrowPeak"]


def test_cli_outdir_cannot_be_created(run_macs3, write_bedgraph, tmp_path):
    """An --outdir below a regular file cannot be created: exit 1 with
    main()'s message."""
    bdg = write_bedgraph(P_ROWS, name="p.bdg")
    blocker = tmp_path / "afile"
    blocker.write_text("x")
    outdir = blocker / "sub"
    res = run_macs3(["bdgpeakcall", "-i", bdg, "--outdir", outdir,
                     "--o-prefix", "P"], timeout=60)
    assert res.returncode == 1
    assert res.stderr.strip() == ("Output directory (%s) could not be "
                                  "created. Terminate program." % outdir)


def test_cli_missing_input_file(run_macs3, tmp_path):
    """A nonexistent -i fails with FileNotFoundError (exit 1) and writes
    no output."""
    missing = tmp_path / "nope.bdg"
    res = run_macs3(["bdgpeakcall", "-i", missing, "--outdir", tmp_path,
                     "--o-prefix", "P"], timeout=60)
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr
    assert str(missing) in res.stderr
    assert not (tmp_path / "P_c5.0_l200_g30_peaks.narrowPeak").exists()


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root can read a mode-000 file")
def test_cli_unreadable_input_file(run_macs3, write_bedgraph, tmp_path):
    """An unreadable -i fails with PermissionError (exit 1)."""
    bdg = Path(write_bedgraph(P_ROWS, name="p.bdg"))
    bdg.chmod(0)
    try:
        res = run_macs3(["bdgpeakcall", "-i", bdg, "--outdir", tmp_path,
                         "--o-prefix", "P"], timeout=60)
    finally:
        bdg.chmod(0o644)
    assert res.returncode == 1
    assert "PermissionError" in res.stderr


@pytest.mark.parametrize("argv, message", [
    (["--o-prefix", "P"],
     "the following arguments are required: -i/--ifile"),
    (["-i", "x.bdg"],
     "one of the arguments -o/--ofile --o-prefix is required"),
    (["-i", "x.bdg", "-o", "a", "--o-prefix", "P"],
     "argument --o-prefix: not allowed with argument -o/--ofile"),
    (["-i", "x.bdg", "-o", "a", "-c", "abc"],
     "argument -c/--cutoff: invalid float value: 'abc'"),
    (["-i", "x.bdg", "-o", "a", "-l", "1.5"],
     "argument -l/--min-length: invalid int value: '1.5'"),
    (["-i", "x.bdg", "-o", "a", "-g", "x"],
     "argument -g/--max-gap: invalid int value: 'x'"),
    (["-i", "x.bdg", "-o", "a", "--cutoff-analysis-max", "2.5"],
     "argument --cutoff-analysis-max: invalid int value: '2.5'"),
    (["-i", "x.bdg", "-o", "a", "--cutoff-analysis-steps", "many"],
     "argument --cutoff-analysis-steps: invalid int value: 'many'"),
    (["-i", "x.bdg", "-o", "a", "--verbose", "loud"],
     "argument --verbose: invalid int value: 'loud'"),
])
def test_cli_argparse_errors(run_macs3, tmp_path, argv, message):
    """Bad or missing arguments: exit 2, usage and the argparse error."""
    res = run_macs3(["bdgpeakcall"] + argv, timeout=60)
    assert res.returncode == 2
    assert res.stderr.startswith("usage: macs3 bdgpeakcall")
    assert res.stderr.rstrip().splitlines()[-1] == (
        "macs3 bdgpeakcall: error: " + message)


def test_cli_unknown_option(run_macs3):
    """An unknown option is rejected by the top-level parser (exit 2)."""
    res = run_macs3(["bdgpeakcall", "-i", "x.bdg", "-o", "a", "--bogus"],
                    timeout=60)
    assert res.returncode == 2
    assert res.stderr.rstrip().splitlines()[-1] == (
        "macs3: error: unrecognized arguments: --bogus")


# ------------------------------------
# realistic run on the CTCF example
# ------------------------------------

def test_ctcf_fe_track_matches_standard(macs3_argparser, test_dir,
                                        tmp_path):
    """bdgpeakcall -c 2 on the CTCF chr22 FE track (bdgcmp output of the
    MACS3 command-line test) reproduces the upstream narrowPeak.

    Pins the current output. The reference file is MACS3's own
    standard result; a peak set over 192k intervals cannot be derived
    by hand. The FE track has no gaps, no negative value and no
    coordinate near 2^30, so none of the bugs marked above applies.
    """
    fe = test_dir / "standard_results_bdgcmp" / "run_bdgcmp_FE.bdg"
    std = (test_dir / "standard_results_bdgpeakcall"
           / "run_bdgpeakcall_w_prefix_c2.0_l200_g30_peaks.narrowPeak")
    call(macs3_argparser, ["-i", fe, "-c", "2", "--outdir", tmp_path,
                           "--o-prefix", "run_bdgpeakcall_w_prefix"])
    out = tmp_path / std.name
    assert lines(out) == lines(std)
