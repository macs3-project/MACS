#!/usr/bin/env python
"""Module Description: Test functions of the randsample subcommand
(MACS3/Commands/randsample_cmd.py), called directly and through
``macs3 randsample``.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import gzip
import logging
import re
from collections import Counter

import numpy as np
import pytest

from MACS3.Commands.randsample_cmd import (run,
                                           load_tag_files_options,
                                           load_frag_files_options)
from MACS3.IO.Parser import (BEDParser,
                             BEDPEParser)

# ------------------------------------
# local helpers
# ------------------------------------

READLEN = 36
REFS = (("chr1", 100000), ("chr2", 100000))
SE_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BOWTIE", "BAM"]
AUTO_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BAM"]
MEM_PREFIX = re.compile(r"^\[\d+ MB\] ")

# (chrom, 0-based leftmost start, strand), every read READLEN long.
# chr1: 10 plus, 6 minus; chr2: 4 plus, 3 minus; 23 reads in total.
SE_READS = ([("chr1", 100 * i, "+") for i in range(1, 11)] +
            [("chr1", 100 * i + 50, "-") for i in range(1, 7)] +
            [("chr2", 200 * i, "+") for i in range(1, 5)] +
            [("chr2", 200 * i + 70, "-") for i in range(1, 4)])
STRAND_COUNTS = [10, 6, 4, 3]       # chr1 +, chr1 -, chr2 +, chr2 -

# 200 reads at distinct positions on chr1, for seed comparisons
BIG_READS = ([("chr1", 50 * i, "+") for i in range(1, 101)] +
             [("chr1", 50 * i + 20, "-") for i in range(1, 101)])

# paired-end fragments: 10 on chr1, 5 on chr2
PE_FRAGS = ([("chr1", 100 * i, 100 * i + 150 + i) for i in range(1, 11)] +
            [("chr2", 300 * i, 300 * i + 200) for i in range(1, 6)])


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


def write_bampe(frags, make_alignments, make_pe_pair, refs=REFS):
    recs = []
    for i, (c, s, e) in enumerate(frags):
        recs += make_pe_pair("p%d" % i, c, s, e)
    return make_alignments(recs, refs=refs, name="frags.bam")


def tag_lines(reads, tsize=READLEN, readlen=READLEN):
    """The BED line randsample writes for each read."""
    out = []
    for c, s, st in reads:
        if st == "+":
            out.append("%s\t%d\t%d\t.\t.\t+" % (c, s, s + tsize))
        else:
            out.append("%s\t%d\t%d\t.\t.\t-" % (c, s + readlen - tsize,
                                                s + readlen))
    return out


def frag_lines(frags):
    return ["%s\t%d\t%d" % f for f in frags]


def kept(n, fraction):
    """Reads kept from n on one strand of one chromosome: the fraction is
    passed as a C float (single precision), multiplied by n in double
    precision, rounded to 5 decimals, then truncated."""
    return int(round(n * float(np.float32(fraction)), 5))


def read_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def is_sub_multiset(sub, full):
    return not (Counter(sub) - Counter(full))


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


def run_randsample(macs3_argparser, args):
    options = macs3_argparser.parse_args(["randsample"] +
                                         [str(a) for a in args])
    run(options)
    return options


@pytest.fixture
def se_bed(write_bed):
    return write_bed(bed_lines(SE_READS), name="reads.bed")


@pytest.fixture
def big_bed(write_bed):
    return write_bed(bed_lines(BIG_READS), name="big.bed")


@pytest.fixture
def pe_bedpe(write_bedpe):
    return write_bedpe(PE_FRAGS, name="frags.bedpe")


# ------------------------------------
# argument parser of `macs3 randsample`
# ------------------------------------

@pytest.mark.parametrize("dest, value", [
    ("percentage", 10.0), ("number", None), ("seed", -1),
    ("outputfile", None), ("outdir", ""), ("tsize", None),
    ("format", "AUTO"), ("buffer_size", 100000), ("verbose", 2)])
def test_parser_defaults(macs3_argparser, dest, value):
    options = macs3_argparser.parse_args(["randsample", "-i", "x", "-p", "10"])
    assert getattr(options, dest) == value


def test_parser_number_is_float(macs3_argparser):
    options = macs3_argparser.parse_args(["randsample", "-i", "x",
                                          "-n", "8e+6"])
    assert options.number == 8e6 and options.percentage is None


@pytest.mark.parametrize("args, message", [
    (["-p", "10"], "the following arguments are required: -i/--ifile"),
    (["-i", "x"],
     "one of the arguments -p/--percentage -n/--number is required"),
    (["-i", "x", "-p", "10", "-n", "5"],
     "argument -n/--number: not allowed with argument -p/--percentage"),
    (["-i", "x", "-p", "ten"],
     "argument -p/--percentage: invalid float value: 'ten'"),
    (["-i", "x", "-p", "10", "--seed", "1.5"],
     "argument --seed: invalid int value: '1.5'"),
    (["-i", "x", "-p", "10", "-f", "FRAG"],
     "argument -f/--format: invalid choice: 'FRAG'"),
])
def test_parser_errors(macs3_argparser, capsys, args, message):
    with pytest.raises(SystemExit) as exc:
        macs3_argparser.parse_args(["randsample"] + args)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


@pytest.mark.parametrize("args", [[], ["-p", "10", "-n", "5"]])
def test_cli_p_or_n_required(run_macs3, se_bed, args):
    proc = run_macs3(["randsample", "-i", se_bed] + args, timeout=60)
    assert proc.returncode == 2
    assert proc.stderr.startswith("usage: macs3 randsample")


# ------------------------------------
# opt_validate_randsample and run() checks through the CLI
# ------------------------------------

def test_cli_percentage_above_100(run_macs3, parse_log, se_bed):
    proc = run_macs3(["randsample", "-i", se_bed, "-p", "100.5"], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "Percentage can't be bigger than 100.0. Please check your "
                  "options and retry!")]
    assert proc.stdout == ""


def test_cli_negative_number(run_macs3, parse_log, se_bed):
    proc = run_macs3(["randsample", "-i", se_bed, "-n", "-5"], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "Number of tags can't be smaller than or equal to 0. "
                  "Please check your options and retry!")]


def test_cli_number_above_total(run_macs3, parse_log, se_bed, tmp_path):
    proc = run_macs3(["randsample", "-i", se_bed, "-n", "100",
                      "--outdir", tmp_path, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 1
    log = parse_log(proc.stderr)
    assert log[-2:] == [
        ("CRITICAL", " Number you want is bigger than total number of tags "
                     "in alignment file! Please specify a smaller number and "
                     "try again!"),
        ("CRITICAL", " 1.00e+02 > 2.30e+01")]
    assert ("INFO", " total tags in alignment file: 23") in log


# ------------------------------------
# -p/--percentage
# ------------------------------------

@pytest.mark.parametrize("percentage", [100, 50, 33.3, 10, 1, 0])
def test_percentage_number_of_reads(macs3_argparser, caplog, se_bed,
                                    tmp_path, percentage):
    run_randsample(macs3_argparser,
                   ["-i", se_bed, "-p", percentage, "--seed", "3",
                    "--outdir", tmp_path, "-o", "out.bed"])
    lines = read_lines(tmp_path / "out.bed")
    n = sum(kept(c, percentage / 100.0) for c in STRAND_COUNTS)
    assert len(lines) == n
    assert is_sub_multiset(lines, tag_lines(SE_READS))
    msgs = [m for _, m in cmd_messages(caplog)]
    assert " Percentage of tags you want to keep: %.2f%%" % percentage \
        in msgs
    assert " tags after random sampling in alignment file: %d" % n in msgs


def test_percentage_100_keeps_every_read(macs3_argparser, se_bed, tmp_path):
    run_randsample(macs3_argparser,
                   ["-i", se_bed, "-p", "100", "--outdir", tmp_path,
                    "-o", "out.bed"])
    lines = read_lines(tmp_path / "out.bed")
    assert sorted(lines) == sorted(tag_lines(SE_READS))


def test_per_strand_counts_at_50_percent(macs3_argparser, se_bed, tmp_path):
    # 10 -> 5, 6 -> 3, 4 -> 2, 3 -> 1 (1.5 truncated)
    run_randsample(macs3_argparser,
                   ["-i", se_bed, "-p", "50", "--seed", "11",
                    "--outdir", tmp_path, "-o", "out.bed"])
    counts = Counter((x.split("\t")[0], x.split("\t")[5])
                     for x in read_lines(tmp_path / "out.bed"))
    assert counts == {("chr1", "+"): 5, ("chr1", "-"): 3,
                      ("chr2", "+"): 2, ("chr2", "-"): 1}


def test_output_is_sorted_within_strand(macs3_argparser, big_bed, tmp_path):
    run_randsample(macs3_argparser,
                   ["-i", big_bed, "-p", "40", "--seed", "5",
                    "--outdir", tmp_path, "-o", "out.bed"])
    lines = read_lines(tmp_path / "out.bed")
    plus = [int(x.split("\t")[1]) for x in lines if x.endswith("+")]
    minus = [int(x.split("\t")[2]) for x in lines if x.endswith("-")]
    assert len(plus) == len(minus) == 40
    assert plus == sorted(plus) and minus == sorted(minus)
    # plus-strand reads come first
    assert lines == [x for x in lines if x.endswith("+")] + \
        [x for x in lines if x.endswith("-")]


# ------------------------------------
# -n/--number
# ------------------------------------

@pytest.mark.parametrize("number", [23, 11, 5, 1.5])
def test_number_of_reads(macs3_argparser, caplog, se_bed, tmp_path, number):
    run_randsample(macs3_argparser,
                   ["-i", se_bed, "-n", number, "--seed", "2",
                    "--outdir", tmp_path, "-o", "out.bed"])
    fraction = float(number) / 23 * 100 / 100.0
    n = sum(kept(c, fraction) for c in STRAND_COUNTS)
    lines = read_lines(tmp_path / "out.bed")
    assert len(lines) == n
    assert n <= number
    assert is_sub_multiset(lines, tag_lines(SE_READS))
    msgs = [m for _, m in cmd_messages(caplog)]
    assert " Number of tags you want to keep: %.2e" % number in msgs
    assert " Percentage of tags you want to keep: %.2f%%" % \
        (float(number) / 23 * 100) in msgs


def test_cli_log_of_number_run(run_macs3, parse_log, se_bed, tmp_path):
    """-n 11 of 23 reads: 47.83% per strand gives 4 + 2 + 1 + 1 = 8."""
    proc = run_macs3(["randsample", "-i", se_bed, "-f", "BED", "-n", "11",
                      "--seed", "1", "--outdir", tmp_path, "-o", "out.bed"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "read tag files..."),
        ("INFO", "# read treatment tags..."),
        ("INFO", "tag size is determined as 36 bps"),
        ("INFO", "tag size = 36"),
        ("INFO", " total tags in alignment file: 23"),
        ("INFO", " Number of tags you want to keep: 1.10e+01"),
        ("INFO", " Percentage of tags you want to keep: 47.83%"),
        ("INFO", " Random seed has been set as: 1"),
        ("INFO", " tags after random sampling in alignment file: 8"),
        ("INFO", "Write to BED file"),
        ("INFO", "finished! Check out.bed.")]
    assert len(read_lines(tmp_path / "out.bed")) == 8


# ------------------------------------
# --seed
# ------------------------------------

def test_same_seed_same_output(macs3_argparser, big_bed, tmp_path):
    for name in ("a.bed", "b.bed"):
        run_randsample(macs3_argparser,
                       ["-i", big_bed, "-p", "50", "--seed", "12345",
                        "--outdir", tmp_path, "-o", name])
    assert (tmp_path / "a.bed").read_text() == \
        (tmp_path / "b.bed").read_text()


def test_cli_same_seed_same_output_across_processes(run_macs3, big_bed):
    outs = [run_macs3(["randsample", "-i", big_bed, "-p", "50",
                       "--seed", "7"], timeout=60).stdout
            for _ in range(2)]
    assert outs[0] == outs[1]
    assert len(outs[0].splitlines()) == 100


@pytest.mark.parametrize("seed", ["0", "2147483647"])
def test_seed_range(macs3_argparser, big_bed, tmp_path, seed):
    for name in ("a.bed", "b.bed"):
        run_randsample(macs3_argparser,
                       ["-i", big_bed, "-p", "50", "--seed", seed,
                        "--outdir", tmp_path, "-o", name])
    a = read_lines(tmp_path / "a.bed")
    assert len(a) == 100 and a == read_lines(tmp_path / "b.bed")


def test_different_seeds_differ(macs3_argparser, big_bed, tmp_path):
    for seed in (1, 2):
        run_randsample(macs3_argparser,
                       ["-i", big_bed, "-p", "50", "--seed", seed,
                        "--outdir", tmp_path, "-o", "s%d.bed" % seed])
    a = read_lines(tmp_path / "s1.bed")
    b = read_lines(tmp_path / "s2.bed")
    assert len(a) == len(b) == 100
    assert a != b


@pytest.mark.parametrize("seed", ["-1", "-7"])
def test_no_seed_differs_between_runs(macs3_argparser, caplog, big_bed,
                                      tmp_path, seed):
    """A negative --seed (the default -1) leaves the generator unseeded."""
    for name in ("a.bed", "b.bed"):
        run_randsample(macs3_argparser,
                       ["-i", big_bed, "-p", "50", "--seed", seed,
                        "--outdir", tmp_path, "-o", name])
    a = read_lines(tmp_path / "a.bed")
    b = read_lines(tmp_path / "b.bed")
    assert len(a) == len(b) == 100
    assert a != b
    assert not any("Random seed" in m for _, m in cmd_messages(caplog))


def test_seed_matches_numpy_shuffle(macs3_argparser, write_bed, tmp_path):
    """With one chromosome and one strand the kept reads are the first
    k of the 5' positions after np.random.seed(seed); shuffle."""
    starts = [100 * i for i in range(1, 21)]
    path = write_bed(bed_lines([("chr1", s, "+") for s in starts]),
                     name="one.bed")
    run_randsample(macs3_argparser,
                   ["-i", path, "-p", "25", "--seed", "42",
                    "--outdir", tmp_path, "-o", "out.bed"])
    arr = np.array(starts, dtype="i4")
    np.random.seed(42)
    np.random.shuffle(arr)
    expect = sorted(arr[:5].tolist())
    assert read_lines(tmp_path / "out.bed") == \
        ["chr1\t%d\t%d\t.\t.\t+" % (s, s + 36) for s in expect]


# ------------------------------------
# -o/--ofile, --outdir, -s/--tsize
# ------------------------------------

def test_cli_stdout_by_default(run_macs3, se_bed, tmp_path):
    proc = run_macs3(["randsample", "-i", se_bed, "-p", "100"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert sorted(proc.stdout.splitlines()) == sorted(tag_lines(SE_READS))


def test_cli_ofile_in_new_outdir(run_macs3, se_bed, tmp_path):
    proc = run_macs3(["randsample", "-i", se_bed, "-p", "100",
                      "--outdir", tmp_path / "sub", "-o", "r.bed"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""
    assert sorted(read_lines(tmp_path / "sub" / "r.bed")) == \
        sorted(tag_lines(SE_READS))


@pytest.mark.parametrize("tsize", [20, 50])
def test_tsize_sets_output_width(macs3_argparser, caplog, se_bed, tmp_path,
                                 tsize):
    run_randsample(macs3_argparser,
                   ["-i", se_bed, "-p", "100", "-s", tsize,
                    "--outdir", tmp_path, "-o", "out.bed"])
    assert sorted(read_lines(tmp_path / "out.bed")) == \
        sorted(tag_lines(SE_READS, tsize=tsize))
    msgs = [m for _, m in cmd_messages(caplog)]
    assert "tag size is determined as %d bps" % tsize in msgs
    assert "tag size = %d" % tsize in msgs


# ------------------------------------
# -f/--format
# ------------------------------------

@pytest.mark.parametrize("fmt", SE_FORMATS)
def test_formats_explicit(macs3_argparser, write_bed, make_alignments,
                          tmp_path, fmt):
    path = write_se(fmt, SE_READS, write_bed, make_alignments)
    run_randsample(macs3_argparser,
                   ["-i", path, "-f", fmt, "-p", "100",
                    "--outdir", tmp_path, "-o", "out.bed"])
    assert sorted(read_lines(tmp_path / "out.bed")) == \
        sorted(tag_lines(SE_READS))


@pytest.mark.parametrize("fmt", AUTO_FORMATS)
def test_formats_auto(macs3_argparser, write_bed, make_alignments, tmp_path,
                      fmt):
    path = write_se(fmt, SE_READS, write_bed, make_alignments)
    run_randsample(macs3_argparser,
                   ["-i", path, "-p", "100", "--outdir", tmp_path,
                    "-o", "out.bed"])
    assert sorted(read_lines(tmp_path / "out.bed")) == \
        sorted(tag_lines(SE_READS))


@pytest.mark.parametrize("fmt", ["BEDPE", "BAMPE"])
@pytest.mark.parametrize("percentage", [100, 50])
def test_paired_end_formats_keep_pairs(macs3_argparser, write_bedpe,
                                       make_alignments, make_pe_pair,
                                       tmp_path, fmt, percentage):
    if fmt == "BEDPE":
        path = write_bedpe(PE_FRAGS)
    else:
        path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair)
    run_randsample(macs3_argparser,
                   ["-i", path, "-f", fmt, "-p", percentage, "--seed", "9",
                    "--outdir", tmp_path, "-o", "out.bedpe"])
    lines = read_lines(tmp_path / "out.bedpe")
    # whole fragments are sampled per chromosome: 10 -> 5 and 5 -> 2
    n = kept(10, percentage / 100.0) + kept(5, percentage / 100.0)
    assert len(lines) == n
    assert all(len(x.split("\t")) == 3 for x in lines)
    assert is_sub_multiset(lines, frag_lines(PE_FRAGS))


def test_cli_paired_end_log(run_macs3, parse_log, pe_bedpe, tmp_path):
    proc = run_macs3(["randsample", "-i", pe_bedpe, "-f", "BEDPE", "-p", "50",
                      "--seed", "1", "--outdir", tmp_path,
                      "-o", "out.bedpe"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "# read input file in Paired-end mode."),
        ("INFO", "# read treatment fragments..."),
        ("INFO", "# total fragments/pairs in alignment file: 15"),
        ("INFO", " Percentage of tags you want to keep: 50.00%"),
        ("INFO", " Random seed has been set as: 1"),
        ("INFO", "#   A random seed 1 has been used"),
        ("INFO", " tags after random sampling in alignment file: 7"),
        ("INFO", "Write to BED file"),
        ("INFO", "finished! Check out.bedpe.")]


def test_paired_end_seed_reproducible(macs3_argparser, pe_bedpe, tmp_path):
    for name in ("a.bedpe", "b.bedpe"):
        run_randsample(macs3_argparser,
                       ["-i", pe_bedpe, "-f", "BEDPE", "-p", "50",
                        "--seed", "4", "--outdir", tmp_path, "-o", name])
    assert (tmp_path / "a.bedpe").read_text() == \
        (tmp_path / "b.bedpe").read_text()


def test_paired_end_number(macs3_argparser, caplog, pe_bedpe, tmp_path):
    run_randsample(macs3_argparser,
                   ["-i", pe_bedpe, "-f", "BEDPE", "-n", "9", "--seed", "4",
                    "--outdir", tmp_path, "-o", "out.bedpe"])
    fraction = 9.0 / 15 * 100 / 100.0
    n = kept(10, fraction) + kept(5, fraction)
    assert len(read_lines(tmp_path / "out.bedpe")) == n
    assert " tags after random sampling in alignment file: %d" % n in \
        [m for _, m in cmd_messages(caplog)]


# ------------------------------------
# several input files, --buffer-size, --verbose
# ------------------------------------

def test_multiple_input_files(macs3_argparser, write_bed, write_bedpe,
                              tmp_path):
    a = write_bed(bed_lines(SE_READS[:12]), name="a.bed")
    b = write_bed(bed_lines(SE_READS[12:]), name="b.bed")
    run_randsample(macs3_argparser, ["-i", a, b, "-p", "100",
                                     "--outdir", tmp_path, "-o", "se.bed"])
    assert sorted(read_lines(tmp_path / "se.bed")) == \
        sorted(tag_lines(SE_READS))
    c = write_bedpe(PE_FRAGS[:4], name="c.bedpe")
    d = write_bedpe(PE_FRAGS[4:], name="d.bedpe")
    run_randsample(macs3_argparser, ["-i", c, d, "-f", "BEDPE", "-p", "100",
                                     "--outdir", tmp_path, "-o", "pe.bedpe"])
    assert sorted(read_lines(tmp_path / "pe.bedpe")) == \
        sorted(frag_lines(PE_FRAGS))


@pytest.mark.parametrize("buffer_size", [1, 5])
def test_buffer_size_does_not_change_output(macs3_argparser, big_bed,
                                            tmp_path, buffer_size):
    run_randsample(macs3_argparser, ["-i", big_bed, "-p", "30", "--seed", "8",
                                     "--outdir", tmp_path, "-o", "ref.bed"])
    run_randsample(macs3_argparser, ["-i", big_bed, "-p", "30", "--seed", "8",
                                     "--buffer-size", buffer_size,
                                     "--outdir", tmp_path, "-o", "buf.bed"])
    assert (tmp_path / "buf.bed").read_text() == \
        (tmp_path / "ref.bed").read_text()


@pytest.mark.parametrize("verbose, n_messages", [(0, 0), (1, 0), (2, 10),
                                                 (3, 10)])
def test_verbose_levels(macs3_argparser, caplog, se_bed, tmp_path, verbose,
                        n_messages):
    run_randsample(macs3_argparser,
                   ["-i", se_bed, "-p", "50", "--seed", "1",
                    "--verbose", verbose, "--outdir", tmp_path,
                    "-o", "out.bed"])
    assert len(cmd_messages(caplog)) == n_messages


# ------------------------------------
# upstream test data
# ------------------------------------

def strand_counts_of(path):
    counts = Counter()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            counts[(f[0], f[5])] += 1
    return counts


def test_cli_upstream_se_number(run_macs3, test_dir, tmp_path):
    """The cmdlinetest run (-n 10000 --seed 31415926): the number of
    reads follows from the per-strand counts of the input, and every
    output line is an input line (all reads are 101 bp)."""
    chip = test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz"
    proc = run_macs3(["randsample", "-i", chip, "-n", "10000",
                      "--seed", "31415926", "--outdir", tmp_path,
                      "-o", "run_randsample.bed"], timeout=120)
    assert proc.returncode == 0, proc.stderr
    counts = strand_counts_of(chip)
    total = sum(counts.values())
    fraction = 10000.0 / total * 100 / 100.0
    out = read_lines(tmp_path / "run_randsample.bed")
    assert len(out) == sum(kept(n, fraction) for n in counts.values())
    with gzip.open(chip, "rt") as fh:
        full = fh.read().splitlines()
    assert is_sub_multiset(out, full)


def test_cli_upstream_contigs50k(run_macs3, test_dir, tmp_path):
    """cmdlinetest 14.4: 50k contigs, --buffer-size 1000. The BED file
    has its strand in column 5, so every read is read as plus strand;
    the count follows from the per-contig counts."""
    contigs = test_dir / "contigs50k.bed.gz"
    proc = run_macs3(["randsample", "-i", contigs, "-n", "100000",
                      "--seed", "31415926", "--outdir", tmp_path,
                      "-o", "run_randsample.bed", "--buffer-size", "1000"],
                     timeout=300)
    assert proc.returncode == 0, proc.stderr
    with gzip.open(contigs, "rt") as fh:
        rows = [line.split("\t")[:3] for line in fh if line.strip()]
    per_contig = Counter(r[0] for r in rows)
    fraction = 100000.0 / len(rows) * 100 / 100.0
    out = [x.split("\t") for x in read_lines(tmp_path / "run_randsample.bed")]
    assert len(out) == sum(kept(n, fraction) for n in per_contig.values())
    assert all(x[3:] == [".", ".", "+"] for x in out)
    assert is_sub_multiset(["\t".join(x[:3]) for x in out],
                           ["\t".join(r) for r in rows])


def test_upstream_bampe_number(macs3_argparser, test_dir, tmp_path):
    bam = test_dir / "CTCF_PE_ChIP_chr22_50k.bam"
    run_randsample(macs3_argparser,
                   ["-f", "BAMPE", "-i", bam, "-n", "10000",
                    "--seed", "31415926",
                    "--outdir", tmp_path, "-o", "pe.bedpe"])
    with gzip.open(test_dir / "CTCF_PE_ChIP_chr22_50k.bedpe.gz", "rt") as fh:
        full = fh.read().splitlines()
    out = read_lines(tmp_path / "pe.bedpe")
    assert len(out) == kept(len(full), 10000.0 / len(full) * 100 / 100.0)
    assert is_sub_multiset(out, full)


# ------------------------------------
# load_tag_files_options / load_frag_files_options
# ------------------------------------

class _Opts:
    def __init__(self, parser, ifile, tsize=None, buffer_size=100000):
        self.parser = parser
        self.ifile = ifile
        self.tsize = tsize
        self.buffer_size = buffer_size
        self.messages = []
        self.info = self.messages.append


def test_load_tag_files_options(se_bed):
    opts = _Opts(BEDParser, [se_bed])
    track = load_tag_files_options(opts)
    assert track.total == 23
    assert opts.tsize == 36
    assert opts.messages == ["# read treatment tags...",
                             "tag size is determined as 36 bps"]
    assert sorted(track.get_chr_names()) == [b"chr1", b"chr2"]


def test_load_tag_files_options_given_tsize_and_two_files(write_bed):
    a = write_bed(bed_lines(SE_READS[:5]), name="a.bed")
    b = write_bed(bed_lines(SE_READS[5:]), name="b.bed")
    opts = _Opts(BEDParser, [a, b], tsize=40)
    track = load_tag_files_options(opts)
    assert track.total == 23 and opts.tsize == 40
    assert opts.messages[-1] == "tag size is determined as 40 bps"


def test_load_frag_files_options(write_bedpe):
    a = write_bedpe(PE_FRAGS[:7], name="a.bedpe")
    b = write_bedpe(PE_FRAGS[7:], name="b.bedpe")
    opts = _Opts(BEDPEParser, [a, b])
    track = load_frag_files_options(opts)
    assert track.total == 15
    assert opts.messages == ["# read treatment fragments..."]
    assert sorted(track.get_chr_names()) == [b"chr1", b"chr2"]


def test_run_in_process_sets_pe_mode(macs3_argparser, pe_bedpe, tmp_path):
    options = run_randsample(macs3_argparser,
                             ["-i", pe_bedpe, "-f", "BEDPE", "-p", "100",
                              "--outdir", tmp_path, "-o", "x.bedpe"])
    assert options.PE_MODE is True
    assert sorted(read_lines(tmp_path / "x.bedpe")) == \
        sorted(frag_lines(PE_FRAGS))


# ------------------------------------
# opt_validate_randsample checks that argparse makes unreachable
# ------------------------------------

def test_validator_rejects_unknown_format(macs3_argparser, caplog, se_bed,
                                          tmp_path):
    """-f has fixed choices, so only a direct run() call reaches this."""
    options = macs3_argparser.parse_args(["randsample", "-i", se_bed,
                                          "-p", "50",
                                          "--outdir", str(tmp_path)])
    options.format = "bogus"
    with pytest.raises(SystemExit) as exc:
        run(options)
    assert exc.value.code == 1
    assert cmd_messages(caplog) == [
        ("ERROR", "Format \"BOGUS\" cannot be recognized!")]


@pytest.mark.parametrize("fmt", ["bed", "auto"])
def test_validator_uppercases_format(macs3_argparser, se_bed, tmp_path, fmt):
    options = macs3_argparser.parse_args(["randsample", "-i", se_bed,
                                          "-p", "100",
                                          "--outdir", str(tmp_path),
                                          "-o", "out.bed"])
    options.format = fmt
    run(options)
    assert options.format == fmt.upper()
    assert sorted(read_lines(tmp_path / "out.bed")) == \
        sorted(tag_lines(SE_READS))


# ------------------------------------
# boundaries: one read, one fragment
# ------------------------------------

@pytest.mark.parametrize("percentage, n", [(100, 1), (99.9, 0), (50, 0)])
def test_single_read(macs3_argparser, write_bed, tmp_path, percentage, n):
    se = write_bed(["chr1\t500\t536\tr\t0\t-"], name="one.bed")
    run_randsample(macs3_argparser, ["-i", se, "-p", percentage,
                                     "--seed", "1", "--outdir", tmp_path,
                                     "-o", "se.bed"])
    assert read_lines(tmp_path / "se.bed") == \
        ["chr1\t500\t536\t.\t.\t-"][:n]


@pytest.mark.parametrize("percentage, n", [(100, 1)])
def test_single_fragment(macs3_argparser, write_bedpe, tmp_path, percentage,
                         n):
    pe = write_bedpe([("chr1", 500, 750)], name="one.bedpe")
    run_randsample(macs3_argparser, ["-i", pe, "-f", "BEDPE",
                                     "-p", percentage, "--seed", "1",
                                     "--outdir", tmp_path, "-o", "pe.bedpe"])
    assert read_lines(tmp_path / "pe.bedpe") == ["chr1\t500\t750"][:n]


