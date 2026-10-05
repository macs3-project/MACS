#!/usr/bin/env python
"""Module Description: Test functions of the refinepeak subcommand
(MACS3/Commands/refinepeak_cmd.py), called directly and through
``macs3 refinepeak``.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import logging
import re

import numpy as np
import pytest

from MACS3.Commands.refinepeak_cmd import (run,
                                           find_summit,
                                           load_tag_files_options)
from MACS3.IO.Parser import BEDParser

# ------------------------------------
# local helpers
# ------------------------------------

READLEN = 36
REFS = (("chr1", 100000), ("chr2", 100000))
SE_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BOWTIE", "BAM"]
AUTO_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BAM"]
MEM_PREFIX = re.compile(r"^\[\d+ MB\] ")


def reads_at(chrom, plus5, minus5, readlen=READLEN):
    """Reads whose 5' ends are ``plus5`` (plus strand) and ``minus5``
    (minus strand); a minus read's 5' end is its BED end."""
    return [(chrom, p, "+") for p in plus5] + \
        [(chrom, m - readlen, "-") for m in minus5]


# five plus-strand 5' ends at 1010 and five minus-strand 5' ends at 1160
REF_READS = reads_at("chr1", [1010] * 5, [1160] * 5)
PEAK1 = "chr1\t1000\t1200\tpeak1"

# three peaks on two chromosomes
MULTI_READS = (REF_READS +
               reads_at("chr1", [5010] * 3, [5100] * 3) +
               reads_at("chr2", [300] * 4, [420] * 4))
MULTI_PEAKS = ["chr2\t250\t450\tpeak3", "chr1\t5000\t5200\tpeak2", PEAK1]


def bed_lines(reads, readlen=READLEN):
    return ["%s\t%d\t%d\tr%d\t0\t%s" % (c, s, s + readlen, i, st)
            for i, (c, s, st) in enumerate(reads)]


def eland_lines(reads, readlen=READLEN):
    return [">r%d\t%s\tU0\t1\t0\t0\t%s.fa\t%d\t%s\t..\t26A" %
            (i, "A" * readlen, c, s + 1, "F" if st == "+" else "R")
            for i, (c, s, st) in enumerate(reads)]


def elandmulti_lines(reads, readlen=READLEN):
    return [">r%d\t%s\t1:0:0\t%s.fa:%d%s0" %
            (i, "A" * readlen, c, s + 1, "F" if st == "+" else "R")
            for i, (c, s, st) in enumerate(reads)]


def elandexport_lines(reads, readlen=READLEN):
    return ["\t".join(["HWUSI", "1", "1", "1", "1000", str(1000 + i), "0",
                       "1", "A" * readlen, "I" * readlen, c, "",
                       str(s + 1), "F" if st == "+" else "R", "", "", "",
                       "", "", "", "", "Y"])
            for i, (c, s, st) in enumerate(reads)]


def bowtie_lines(reads, readlen=READLEN):
    return ["r%d\t%s\t%s\t%d\t%s\t%s\t0\t" %
            (i, st, c, s, "A" * readlen, "I" * readlen)
            for i, (c, s, st) in enumerate(reads)]


TEXT_WRITERS = {"BED": bed_lines, "ELAND": eland_lines,
                "ELANDMULTI": elandmulti_lines,
                "ELANDEXPORT": elandexport_lines, "BOWTIE": bowtie_lines}


def write_se(fmt, reads, write_bed, make_alignments, refs=REFS,
             readlen=READLEN):
    if fmt in ("SAM", "BAM"):
        recs = [dict(name="r%d" % i, ref=c, pos=s,
                     flag=0 if st == "+" else 16, cigar="%dM" % readlen)
                for i, (c, s, st) in enumerate(reads)]
        return make_alignments(recs, refs=refs, fmt=fmt.lower(),
                               name="reads." + fmt.lower())
    return write_bed(TEXT_WRITERS[fmt](reads, readlen),
                     name="reads." + fmt.lower())


def read_text(path):
    with open(path) as fh:
        return fh.read()


def cmd_messages(caplog):
    return [(r.levelname, MEM_PREFIX.sub("", r.getMessage()))
            for r in caplog.records
            if r.name == "MACS3.Utilities.OptValidator"]


@pytest.fixture(autouse=True)
def _restore_optvalidator_level():
    """run() sets the level of the OptValidator logger from --verbose;
    restore it after each test."""
    lg = logging.getLogger("MACS3.Utilities.OptValidator")
    level = lg.level
    yield
    lg.setLevel(level)


@pytest.fixture
def run_macs3(run_macs3):
    """The conftest runner with single-threaded OpenMP/BLAS pools."""
    def _run(args, env=None, **kwargs):
        one = {v: "1" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                "MKL_NUM_THREADS")}
        one.update(env or {})
        return run_macs3(args, env=one, **kwargs)
    return _run


def run_refinepeak(macs3_argparser, args):
    options = macs3_argparser.parse_args(["refinepeak"] +
                                         [str(a) for a in args])
    run(options)
    return options


@pytest.fixture
def ref_bed(write_bed):
    return write_bed(bed_lines(REF_READS), name="reads.bed")


@pytest.fixture
def peak1(write_bed):
    return write_bed([PEAK1], name="peak1.bed")


# ------------------------------------
# argument parser of `macs3 refinepeak`
# ------------------------------------

@pytest.mark.parametrize("dest, value", [
    ("format", "AUTO"), ("cutoff", 5.0), ("windowsize", 200),
    ("buffer_size", 100000), ("verbose", 2), ("outdir", ""),
    ("ofile", None), ("oprefix", "p")])
def test_parser_defaults(macs3_argparser, dest, value):
    options = macs3_argparser.parse_args(["refinepeak", "-b", "x", "-i", "y",
                                          "--o-prefix", "p"])
    assert getattr(options, dest) == value


@pytest.mark.parametrize("args, message", [
    (["-i", "y", "-o", "z"], "the following arguments are required: -b"),
    (["-b", "x", "-o", "z"],
     "the following arguments are required: -i/--ifile"),
    (["-b", "x", "-i", "y"],
     "one of the arguments -o/--ofile --o-prefix is required"),
    (["-b", "x", "-i", "y", "-o", "z", "--o-prefix", "p"],
     "argument --o-prefix: not allowed with argument -o/--ofile"),
    (["-b", "x", "-i", "y", "-o", "z", "-f", "BAMPE"],
     "argument -f/--format: invalid choice: 'BAMPE'"),
    (["-b", "x", "-i", "y", "-o", "z", "-w", "1.5"],
     "argument -w/--window-size: invalid int value: '1.5'"),
    (["-b", "x", "-i", "y", "-o", "z", "-c", "high"],
     "argument -c/--cutoff: invalid float value: 'high'"),
])
def test_parser_errors(macs3_argparser, capsys, args, message):
    with pytest.raises(SystemExit) as exc:
        macs3_argparser.parse_args(["refinepeak"] + args)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_cli_output_option_required(run_macs3, ref_bed, peak1):
    proc = run_macs3(["refinepeak", "-b", peak1, "-i", ref_bed], timeout=60)
    assert proc.returncode == 2
    assert proc.stderr.startswith("usage: macs3 refinepeak")
    assert "one of the arguments -o/--ofile --o-prefix is required" in \
        proc.stderr


# ------------------------------------
# find_summit
# ------------------------------------
# The score at position j is
#   2 * sqrt(W_left * C_right) - W_right - C_left
# where W/C count plus/minus 5' ends in the window_size-wide window left
# (ending at j) or right (starting at j) of j.  The summit is the first
# position with the largest score; "_R" if the score > cutoff, else "_F".

def i4(values):
    return np.array(values, dtype="i4")


@pytest.mark.parametrize("plus, minus, start, end, w, expect_pos, expect", [
    # 5 plus at 1010, 5 minus at 1160: both windows hold 5 tags for j in
    # [1011, 1160] -> 2 * sqrt(25) = 10, first at 1011
    ([1010] * 5, [1160] * 5, 800, 1400, 200, 1011, 10.0),
    # 4 plus and 9 minus -> 2 * sqrt(36) = 12
    ([1010] * 4, [1160] * 9, 800, 1400, 200, 1011, 12.0),
    # with 50-bp windows the strands never pair: best score 0 at the start
    ([1010] * 5, [1160] * 5, 950, 1250, 50, 950, 0.0),
    # strands 40 bp apart and 50-bp windows: 10 for j in [1021, 1060]
    ([1020] * 5, [1060] * 5, 950, 1150, 50, 1021, 10.0),
    # no tags at all: score 0 everywhere
    ([], [], 500, 900, 200, 500, 0.0),
    # plus tags only: never positive, 0 where the right window is empty
    ([1010] * 3, [], 800, 1400, 200, 800, 0.0),
])
def test_find_summit(plus, minus, start, end, w, expect_pos, expect):
    chrom, pos, pos1, name, score = find_summit(
        b"chr1", i4(plus), i4(minus), start, end, name=b"pk",
        window_size=w, cutoff=5)
    assert (chrom, pos, pos1) == (b"chr1", expect_pos, expect_pos + 1)
    assert score == expect
    assert name == (b"pk_R" if expect > 5 else b"pk_F")


@pytest.mark.parametrize("cutoff, suffix", [(9.99, b"_R"), (10, b"_F"),
                                            (10.01, b"_F"), (-1, b"_R")])
def test_find_summit_cutoff_is_strict(cutoff, suffix):
    result = find_summit(b"chr1", i4([1010] * 5), i4([1160] * 5), 800, 1400,
                         name=b"p", window_size=200, cutoff=cutoff)
    assert result[3] == b"p" + suffix
    assert result[4] == 10.0


def test_find_summit_default_name_and_types():
    result = find_summit(b"chrX", i4([1010] * 5), i4([1160] * 5), 800, 1400,
                         window_size=200)
    assert result == (b"chrX", 1011, 1012, b"peak_R", 10.0)
    assert isinstance(result[4], float)


def test_find_summit_scores_unequal_clusters():
    """Two plus clusters: 2 tags at 1000 and 3 at 1050; 6 minus at 1150.
    For j in [1051, 1150] both left plus tags and all minus tags count:
    2 * sqrt(5 * 6) = 10.954..."""
    result = find_summit(b"chr1", i4([1000] * 2 + [1050] * 3),
                         i4([1150] * 6), 801, 1400, window_size=200)
    assert result[1] == 1051
    assert result[4] == pytest.approx(2 * (5 * 6) ** 0.5, rel=1e-12)


# ------------------------------------
# -w/--window-size, -c/--cutoff, -o/--o-prefix
# ------------------------------------

@pytest.mark.parametrize("w, line", [
    # scan [800, 1400]: summit 1011, score 10
    ("200", "chr1\t1011\t1012\tpeak1_R\t10.00"),
    # scan [950, 1250] with 50-bp windows: the strands never pair
    ("50", "chr1\t950\t951\tpeak1_F\t0.00"),
    # 150-bp windows reach from 1011 to 1160: same summit
    ("150", "chr1\t1011\t1012\tpeak1_R\t10.00"),
])
def test_window_size(macs3_argparser, ref_bed, peak1, tmp_path, w, line):
    run_refinepeak(macs3_argparser, ["-b", peak1, "-i", ref_bed, "-f", "BED",
                                     "-w", w, "--outdir", tmp_path,
                                     "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == [line]


@pytest.mark.parametrize("cutoff, line", [
    ("5", "chr1\t1011\t1012\tpeak1_R\t10.00"),
    ("9.99", "chr1\t1011\t1012\tpeak1_R\t10.00"),
    ("10", "chr1\t1011\t1012\tpeak1_F\t10.00"),
    ("50", "chr1\t1011\t1012\tpeak1_F\t10.00"),
])
def test_cutoff(macs3_argparser, ref_bed, peak1, tmp_path, cutoff, line):
    run_refinepeak(macs3_argparser, ["-b", peak1, "-i", ref_bed, "-f", "BED",
                                     "-c", cutoff, "--outdir", tmp_path,
                                     "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == [line]


def test_cli_ofile(run_macs3, parse_log, ref_bed, peak1, tmp_path):
    proc = run_macs3(["refinepeak", "-b", peak1, "-i", ref_bed, "-f", "BED",
                      "--outdir", tmp_path / "res", "-o", "refined.bed"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert [f.name for f in (tmp_path / "res").iterdir()] == ["refined.bed"]
    assert read_text(tmp_path / "res" / "refined.bed").splitlines() == \
        ["chr1\t1011\t1012\tpeak1_R\t10.00"]
    assert parse_log(proc.stderr) == [("INFO", "read tag files..."),
                                      ("INFO", "# read treatment tags..."),
                                      ("INFO", "Done!")]


def test_cli_o_prefix(run_macs3, ref_bed, peak1, tmp_path):
    proc = run_macs3(["refinepeak", "-b", peak1, "-i", ref_bed,
                      "--outdir", tmp_path / "res", "--o-prefix", "sample"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert [f.name for f in (tmp_path / "res").iterdir()] == \
        ["sample_refinepeak.bed"]
    out = tmp_path / "res" / "sample_refinepeak.bed"
    assert read_text(out).splitlines() == ["chr1\t1011\t1012\tpeak1_R\t10.00"]


def test_peaks_sorted_by_chromosome_then_start(macs3_argparser, write_bed,
                                               tmp_path):
    """Hand-derived summits: chr1 1011 (5 + 5 tags -> 10), chr1 5011
    (3 + 3 -> 6), chr2 301 (4 + 4 -> 8); peaks come out sorted."""
    reads = write_bed(bed_lines(MULTI_READS), name="multi.bed")
    peaks = write_bed(MULTI_PEAKS, name="peaks.bed")
    run_refinepeak(macs3_argparser, ["-b", peaks, "-i", reads, "-f", "BED",
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == [
        "chr1\t1011\t1012\tpeak1_R\t10.00",
        "chr1\t5011\t5012\tpeak2_R\t6.00",
        "chr2\t301\t302\tpeak3_R\t8.00"]


def test_narrowpeak_input_and_space_separated_peaks(macs3_argparser,
                                                    write_bed, ref_bed,
                                                    tmp_path):
    peaks = write_bed(["chr1 1000 1200 peak1 100 . 5.0 10.0 8.0 100"],
                      name="peaks.narrowPeak")
    run_refinepeak(macs3_argparser, ["-b", peaks, "-i", ref_bed,
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == \
        ["chr1\t1011\t1012\tpeak1_R\t10.00"]


def test_peak_without_tags_fails(macs3_argparser, write_bed, ref_bed,
                                 tmp_path):
    peaks = write_bed([PEAK1, "chr1\t50000\t50300\tempty"], name="p.bed")
    run_refinepeak(macs3_argparser, ["-b", peaks, "-i", ref_bed,
                                     "--outdir", tmp_path, "-o", "out.bed"])
    # no tag within [49800, 50500]: score 0 at the first scanned base
    assert read_text(tmp_path / "out.bed").splitlines() == [
        "chr1\t1011\t1012\tpeak1_R\t10.00",
        "chr1\t49800\t49801\tempty_F\t0.00"]


def test_single_peaks_of_the_overlap_case(macs3_argparser, write_bed,
                                          tmp_path):
    reads = write_bed(bed_lines(reads_at("chr1", [860] * 5 + [1301],
                                         [1010] * 5)), name="r.bed")
    for name, line in (("pa", "chr1\t1000\t1100\tpa"),
                       ("pb", "chr1\t1050\t1300\tpb")):
        peaks = write_bed([line], name=name + ".bed")
        run_refinepeak(macs3_argparser, ["-b", peaks, "-i", reads,
                                         "--outdir", tmp_path,
                                         "-o", name + ".out"])
        assert read_text(tmp_path / (name + ".out")).splitlines() == \
            ["chr1\t861\t862\t%s_R\t10.00" % name]


# ------------------------------------
# -f/--format
# ------------------------------------

@pytest.mark.parametrize("fmt", SE_FORMATS)
def test_formats_explicit(macs3_argparser, write_bed, make_alignments, peak1,
                          tmp_path, fmt):
    path = write_se(fmt, REF_READS, write_bed, make_alignments)
    run_refinepeak(macs3_argparser, ["-b", peak1, "-i", path, "-f", fmt,
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == \
        ["chr1\t1011\t1012\tpeak1_R\t10.00"]


@pytest.mark.parametrize("fmt", AUTO_FORMATS)
def test_formats_auto(macs3_argparser, write_bed, make_alignments, peak1,
                      tmp_path, fmt):
    path = write_se(fmt, REF_READS, write_bed, make_alignments)
    run_refinepeak(macs3_argparser, ["-b", peak1, "-i", path,
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == \
        ["chr1\t1011\t1012\tpeak1_R\t10.00"]


def test_cli_auto_log(run_macs3, parse_log, write_bed, peak1, tmp_path):
    path = write_bed(bed_lines(REF_READS), name="reads.bed.gz")
    proc = run_macs3(["refinepeak", "-b", peak1, "-i", path,
                      "--outdir", tmp_path, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [("INFO", "read tag files..."),
                                      ("INFO", "# read treatment tags..."),
                                      ("INFO", "Detected format is: BED"),
                                      ("INFO", "* Input file is gzipped."),
                                      ("INFO", "Done!")]


# ------------------------------------
# several input files, --buffer-size, --verbose
# ------------------------------------

def test_multiple_input_files(macs3_argparser, write_bed, tmp_path):
    a = write_bed(bed_lines(MULTI_READS[:9]), name="a.bed")
    b = write_bed(bed_lines(MULTI_READS[9:]), name="b.bed")
    peaks = write_bed(MULTI_PEAKS, name="peaks.bed")
    run_refinepeak(macs3_argparser, ["-b", peaks, "-i", a, b, "-f", "BED",
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == [
        "chr1\t1011\t1012\tpeak1_R\t10.00",
        "chr1\t5011\t5012\tpeak2_R\t6.00",
        "chr2\t301\t302\tpeak3_R\t8.00"]


@pytest.mark.parametrize("buffer_size", [1, 3])
def test_buffer_size_does_not_change_output(macs3_argparser, ref_bed, peak1,
                                            tmp_path, buffer_size):
    run_refinepeak(macs3_argparser, ["-b", peak1, "-i", ref_bed,
                                     "--buffer-size", buffer_size,
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert read_text(tmp_path / "out.bed").splitlines() == \
        ["chr1\t1011\t1012\tpeak1_R\t10.00"]


@pytest.mark.parametrize("verbose, n_messages", [(0, 0), (1, 0), (2, 3),
                                                 (3, 3)])
def test_verbose_levels(macs3_argparser, caplog, ref_bed, peak1, tmp_path,
                        verbose, n_messages):
    run_refinepeak(macs3_argparser, ["-b", peak1, "-i", ref_bed, "-f", "BED",
                                     "--verbose", verbose,
                                     "--outdir", tmp_path, "-o", "out.bed"])
    assert len(cmd_messages(caplog)) == n_messages


# ------------------------------------
# upstream test data (test/cmdlinetest)
# ------------------------------------

CTCF_OUTPUTS = [
    (["-o", "run_refinepeak_w_ofile.bed"], "run_refinepeak_w_ofile.bed",
     "run_refinepeak_w_ofile.bed"),
    (["--o-prefix", "run_refinepeak_w_prefix"],
     "run_refinepeak_w_prefix_refinepeak.bed",
     "run_refinepeak_w_prefix_refinepeak.bed"),
]


def run_ctcf(macs3_argparser, test_dir, tmp_path, outopt):
    """The cmdlinetest run: the CTCF chr22 narrowPeak refined with the
    CTCF chr22 ChIP reads and the defaults (-w 200, -c 5)."""
    peaks = test_dir / "standard_results_callpeak_narrow" / \
        "run_callpeak_narrow0_peaks.narrowPeak"
    run_refinepeak(macs3_argparser,
                   ["-b", peaks, "-i",
                    test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz",
                    "--outdir", tmp_path] + outopt)


@pytest.mark.parametrize("outopt, outname, std", CTCF_OUTPUTS)
def test_standard_results(macs3_argparser, test_dir, tmp_path, outopt,
                          outname, std):
    """One line per peak, in the order of upstream's
    test/standard_results_refinepeak, with the same chromosome and peak
    name. The summits, scores and _R/_F suffixes of the stored files
    depend on find_summit's window counting, which is not checked here,
    so they are not compared."""
    run_ctcf(macs3_argparser, test_dir, tmp_path, outopt)
    got = [x.split("\t") for x in read_text(tmp_path / outname).splitlines()]
    exp = [x.split("\t") for x in read_text(
        test_dir / "standard_results_refinepeak" / std).splitlines()]
    assert len(got) == len(exp) == 730
    assert ([(r[0], r[3].rsplit("_", 1)[0]) for r in got]
            == [(r[0], r[3].rsplit("_", 1)[0]) for r in exp])


# ------------------------------------
# load_tag_files_options
# ------------------------------------

class _Opts:
    def __init__(self, parser, ifile, buffer_size=100000):
        self.parser = parser
        self.ifile = ifile
        self.buffer_size = buffer_size
        self.messages = []
        self.info = self.messages.append


def test_load_tag_files_options(write_bed):
    a = write_bed(bed_lines(MULTI_READS[:10]), name="a.bed")
    b = write_bed(bed_lines(MULTI_READS[10:]), name="b.bed")
    opts = _Opts(BEDParser, [a, b])
    track = load_tag_files_options(opts)
    assert opts.messages == ["# read treatment tags..."]
    assert track.total == len(MULTI_READS)
    plus, minus = track.get_locations_by_chr(b"chr2")
    assert list(plus) == [300] * 4 and list(minus) == [420] * 4


# ------------------------------------
# opt_validate_refinepeak checks that argparse makes unreachable
# ------------------------------------

def test_validator_rejects_unknown_format(macs3_argparser, caplog, ref_bed,
                                          peak1, tmp_path):
    """-f has fixed choices, so only a direct run() call reaches this."""
    options = macs3_argparser.parse_args(["refinepeak", "-b", peak1,
                                          "-i", ref_bed,
                                          "--outdir", str(tmp_path),
                                          "-o", "out.bed"])
    options.format = "bampe"
    with pytest.raises(SystemExit) as exc:
        run(options)
    assert exc.value.code == 1
    assert cmd_messages(caplog) == [
        ("ERROR", "Format \"BAMPE\" cannot be recognized!")]
    assert not (tmp_path / "out.bed").exists()


@pytest.mark.parametrize("fmt", ["bed", "auto"])
def test_validator_uppercases_format(macs3_argparser, ref_bed, peak1,
                                     tmp_path, fmt):
    options = macs3_argparser.parse_args(["refinepeak", "-b", peak1,
                                          "-i", ref_bed,
                                          "--outdir", str(tmp_path),
                                          "-o", "out.bed"])
    options.format = fmt
    run(options)
    assert options.format == fmt.upper()
    assert read_text(tmp_path / "out.bed").splitlines() == \
        ["chr1\t1011\t1012\tpeak1_R\t10.00"]
