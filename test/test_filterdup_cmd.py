#!/usr/bin/env python
"""Module Description: Test functions of the filterdup subcommand
(MACS3/Commands/filterdup_cmd.py), called directly and through
``macs3 filterdup``.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import gzip
import logging
import os
import re

import pytest
from scipy.stats import binom

from MACS3.Commands.filterdup_cmd import (run,
                                          cal_max_dup_tags,
                                          load_tag_files_options,
                                          load_frag_files_options)
from MACS3.IO.Parser import (BEDParser,
                             BEDPEParser)
from MACS3.Utilities.Constants import EFFECTIVEGS

# ------------------------------------
# local helpers
# ------------------------------------

READLEN = 36
REFS = (("chr1", 100000),)
SE_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BOWTIE", "BAM"]

# small single-end data set on chr1: (chrom, 0-based leftmost start,
# strand); every read is READLEN bp long.  5' ends: plus = start,
# minus = start + READLEN.
PLUS = [100, 100, 100, 250, 400, 400, 1000, 5000]
MINUS = [300, 300, 520, 700, 700, 700, 2000]
SE_READS = sorted([("chr1", s, "+") for s in PLUS] +
                  [("chr1", s, "-") for s in MINUS], key=lambda r: r[1])

# paired-end fragments (chrom, left, right)
PE_FRAGS = [("chr1", 100, 300), ("chr1", 100, 300), ("chr1", 100, 300),
            ("chr1", 100, 350), ("chr1", 150, 400), ("chr1", 500, 650),
            ("chr1", 500, 650), ("chr1", 1000, 1301)]

MEM_PREFIX = re.compile(r"^\[\d+ MB\] ")


def bed_lines(reads, readlen=READLEN):
    return ["%s\t%d\t%d\tr%d\t0\t%s" % (c, s, s + readlen, i, st)
            for i, (c, s, st) in enumerate(reads)]


def eland_lines(reads, readlen=READLEN):
    # name, sequence, match code, #U0, #U1, #U2, chrom file,
    # 1-based leftmost position, strand F/R, ...
    return [">r%d\t%s\tU0\t1\t0\t0\t%s.fa\t%d\t%s\t..\t26A" %
            (i, "A" * readlen, c, s + 1, "F" if st == "+" else "R")
            for i, (c, s, st) in enumerate(reads)]


def elandmulti_lines(reads, readlen=READLEN):
    # name, sequence, hit counts, chrom.fa:<1-based pos><F|R><mismatches>
    return [">r%d\t%s\t1:0:0\t%s.fa:%d%s0" %
            (i, "A" * readlen, c, s + 1, "F" if st == "+" else "R")
            for i, (c, s, st) in enumerate(reads)]


def elandexport_lines(reads, readlen=READLEN):
    # 22 columns; 9th sequence, 11th chrom, 13th 1-based position,
    # 14th strand
    return ["\t".join(["HWUSI", "1", "1", "1", "1000", str(1000 + i), "0",
                       "1", "A" * readlen, "I" * readlen, c, "",
                       str(s + 1), "F" if st == "+" else "R", "", "", "",
                       "", "", "", "", "Y"])
            for i, (c, s, st) in enumerate(reads)]


def bowtie_lines(reads, readlen=READLEN):
    # name, strand, chrom, 0-based leftmost offset, sequence, qualities,
    # count, mismatches
    return ["r%d\t%s\t%s\t%d\t%s\t%s\t0\t" %
            (i, st, c, s, "A" * readlen, "I" * readlen)
            for i, (c, s, st) in enumerate(reads)]


TEXT_WRITERS = {"BED": bed_lines, "ELAND": eland_lines,
                "ELANDMULTI": elandmulti_lines,
                "ELANDEXPORT": elandexport_lines, "BOWTIE": bowtie_lines}


def write_se(fmt, reads, write_bed, make_alignments, refs=REFS,
             readlen=READLEN, name=None):
    """Write ``reads`` in format ``fmt`` and return the path."""
    if fmt in ("SAM", "BAM"):
        recs = [dict(name="r%d" % i, ref=c, pos=s,
                     flag=0 if st == "+" else 16, cigar="%dM" % readlen)
                for i, (c, s, st) in enumerate(reads)]
        return make_alignments(recs, refs=refs, fmt=fmt.lower(),
                               name=name or "reads." + fmt.lower())
    return write_bed(TEXT_WRITERS[fmt](reads, readlen),
                     name=name or "reads." + fmt.lower())


def write_bampe(frags, make_alignments, make_pe_pair, refs=REFS,
                name="frags.bam"):
    recs = []
    for i, (c, s, e) in enumerate(frags):
        recs += make_pe_pair("p%d" % i, c, s, e)
    return make_alignments(recs, refs=refs, name=name)


def cap(values, maxdup):
    """Keep at most ``maxdup`` copies of each value (input sorted)."""
    if maxdup is None:
        return list(values)
    out, seen = [], {}
    for v in values:
        seen[v] = seen.get(v, 0) + 1
        if seen[v] <= maxdup:
            out.append(v)
    return out


def expected_bed(reads, tsize=READLEN, maxdup=None, readlen=READLEN,
                 chroms=None):
    """BED lines written by filterdup: per chromosome, plus-strand tags
    sorted by 5' end, then minus-strand tags sorted by 5' end, each
    ``tsize`` long and anchored at the 5' end."""
    out = []
    for chrom in chroms or sorted(set(r[0] for r in reads)):
        plus = cap(sorted(s for c, s, st in reads
                          if c == chrom and st == "+"), maxdup)
        minus = cap(sorted(s + readlen for c, s, st in reads
                           if c == chrom and st == "-"), maxdup)
        out += ["%s\t%d\t%d\t.\t.\t+" % (chrom, p, p + tsize) for p in plus]
        out += ["%s\t%d\t%d\t.\t.\t-" % (chrom, q - tsize, q) for q in minus]
    return out


def expected_bedpe(frags, maxdup=None):
    out = []
    for chrom in sorted(set(f[0] for f in frags)):
        locs = cap(sorted((s, e) for c, s, e in frags if c == chrom),
                   maxdup)
        out += ["%s\t%d\t%d" % (chrom, s, e) for s, e in locs]
    return out


def binom_cutoff(n, gsize, p=1e-5):
    """Smallest x with P(X <= x) > 1 - p for X ~ Binomial(n, 1/gsize)."""
    x = 0
    while binom.cdf(x, n, 1.0 / gsize) <= 1 - p:
        x += 1
    return x


def read_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def cmd_messages(caplog):
    """(level, message) of the records logged by the command itself."""
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


def run_filterdup(macs3_argparser, args):
    """Parse ``macs3 filterdup <args>`` and call run() in-process."""
    options = macs3_argparser.parse_args(["filterdup"] +
                                         [str(a) for a in args])
    run(options)
    return options


@pytest.fixture
def se_bed(write_bed):
    return write_bed(bed_lines(SE_READS), name="reads.bed")


@pytest.fixture
def pe_bedpe(write_bedpe):
    return write_bedpe(PE_FRAGS, name="frags.bedpe")


# ------------------------------------
# argument parser of `macs3 filterdup`
# ------------------------------------

@pytest.mark.parametrize("dest, value", [
    ("format", "AUTO"), ("gsize", "hs"), ("tsize", None), ("pvalue", None),
    ("keepduplicates", "auto"), ("buffer_size", 100000), ("verbose", 2),
    ("outdir", ""), ("outputfile", "stdout"), ("dryrun", False)])
def test_parser_defaults(macs3_argparser, dest, value):
    options = macs3_argparser.parse_args(["filterdup", "-i", "x.bed"])
    assert getattr(options, dest) == value


def test_parser_ifile_takes_several_files(macs3_argparser):
    options = macs3_argparser.parse_args(["filterdup", "-i", "a", "b", "c"])
    assert options.ifile == ["a", "b", "c"]


@pytest.mark.parametrize("args, message", [
    ([], "the following arguments are required: -i/--ifile"),
    (["-i", "x", "-f", "FRAG"],
     "argument -f/--format: invalid choice: 'FRAG'"),
    (["-i", "x", "-f", "bed"], "argument -f/--format: invalid choice: 'bed'"),
    (["-i", "x", "-s", "abc"],
     "argument -s/--tsize: invalid int value: 'abc'"),
    (["-i", "x", "-p", "abc"],
     "argument -p/--pvalue: invalid float value: 'abc'"),
    (["-i", "x", "--buffer-size", "1.5"],
     "argument --buffer-size: invalid int value: '1.5'"),
    (["-i", "x", "--verbose", "high"],
     "argument --verbose: invalid int value: 'high'"),
])
def test_parser_errors(macs3_argparser, capsys, args, message):
    with pytest.raises(SystemExit) as exc:
        macs3_argparser.parse_args(["filterdup"] + args)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_cli_argparse_error_exit_code_and_usage(run_macs3):
    proc = run_macs3(["filterdup"], timeout=60)
    assert proc.returncode == 2
    assert proc.stderr.startswith("usage: macs3 filterdup")
    assert "error: the following arguments are required: -i/--ifile" \
        in proc.stderr


# ------------------------------------
# opt_validate_filterdup reached through the CLI
# ------------------------------------

def test_cli_invalid_gsize(run_macs3, parse_log, se_bed, tmp_path):
    proc = run_macs3(["filterdup", "-i", se_bed, "-g", "xyz",
                      "--outdir", tmp_path, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "Error when interpreting --gsize option: xyz"),
        ("ERROR", "Available shortcuts of effective genome sizes are "
                  "hs,mm,ce,dm")]
    assert not (tmp_path / "out.bed").exists()


def test_cli_invalid_keep_dup(run_macs3, parse_log, se_bed, tmp_path):
    proc = run_macs3(["filterdup", "-i", se_bed, "--keep-dup", "abc",
                      "--outdir", tmp_path], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "--keep-dup should be 'auto', 'all' or an integer!")]
    assert proc.stdout == ""


@pytest.mark.parametrize("keepdup", ["-1", "1.5", "AUTO", "All", ""])
def test_invalid_keep_dup_values(macs3_argparser, caplog, se_bed, tmp_path,
                                 keepdup):
    with pytest.raises(SystemExit) as exc:
        run_filterdup(macs3_argparser, ["-i", se_bed, "--keep-dup", keepdup,
                                        "--outdir", tmp_path])
    assert exc.value.code == 1
    assert cmd_messages(caplog) == [
        ("ERROR", "--keep-dup should be 'auto', 'all' or an integer!")]


@pytest.mark.parametrize("gsize, value", [
    ("hs", EFFECTIVEGS["hs"]), ("mm", EFFECTIVEGS["mm"]),
    ("ce", EFFECTIVEGS["ce"]), ("dm", EFFECTIVEGS["dm"]),
    ("1e6", 1e6), ("52000000", 52000000.0)])
def test_gsize_values(macs3_argparser, se_bed, tmp_path, gsize, value):
    options = run_filterdup(macs3_argparser,
                            ["-i", se_bed, "-g", gsize, "--dry-run",
                             "--outdir", tmp_path])
    assert options.gsize == value
    assert EFFECTIVEGS == {"hs": 2913022398, "mm": 2652783500,
                           "ce": 100286401, "dm": 142573017}


# ------------------------------------
# --keep-dup
# ------------------------------------

@pytest.mark.parametrize("keepdup, gsize, maxdup", [
    ("all", "hs", None),
    ("1", "hs", 1),
    ("2", "hs", 2),
    ("3", "hs", 3),
    ("auto", "1000", binom_cutoff(15, 1000)),       # = 2
    ("auto", "100000", binom_cutoff(15, 100000)),   # = 1
])
def test_keep_dup_values(macs3_argparser, se_bed, tmp_path, keepdup, gsize,
                         maxdup):
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "-f", "BED", "--keep-dup", keepdup,
                   "-g", gsize, "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == \
        expected_bed(SE_READS, maxdup=maxdup)


def test_keep_dup_int32_max_keeps_all(macs3_argparser, se_bed, pe_bedpe,
                                      tmp_path):
    for path, fmt, out in ((se_bed, "BED", "se.bed"),
                           (pe_bedpe, "BEDPE", "pe.bedpe")):
        run_filterdup(macs3_argparser,
                      ["-i", path, "-f", fmt, "--keep-dup", "2147483647",
                       "--outdir", tmp_path, "-o", out])
    assert read_lines(tmp_path / "se.bed") == expected_bed(SE_READS)
    assert read_lines(tmp_path / "pe.bedpe") == expected_bedpe(PE_FRAGS)


def test_binomial_cutoffs_of_the_small_data_set():
    # the two auto cases above use different cutoffs
    assert binom_cutoff(15, 1000) == 2
    assert binom_cutoff(15, 100000) == 1


def test_cli_keep_dup_auto_log(run_macs3, parse_log, se_bed, tmp_path):
    """Exact log of an auto run: 15 tags, binomial cutoff 2 for
    -g 1000, 13 tags kept (one extra copy at 100+ and at 700-), redundant
    rate 2/15."""
    maxdup = binom_cutoff(15, 1000)
    kept = len(expected_bed(SE_READS, maxdup=maxdup))
    proc = run_macs3(["filterdup", "-i", se_bed, "-f", "BED", "-g", "1000",
                      "--outdir", tmp_path, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "# read tag files..."),
        ("INFO", "# read treatment tags..."),
        ("INFO", "# tag size = 36"),
        ("INFO", "# total tags in alignment file: 15"),
        ("INFO", "calculate max duplicate tags in single position based "
                 "on binomal distribution..."),
        ("INFO", " max_dup_tags based on binomal = %d" % maxdup),
        ("INFO", "filter out redundant tags at the same location and the "
                 "same strand by allowing at most %d tag(s)" % maxdup),
        ("INFO", " tags after filtering in alignment file: %d" % kept),
        ("INFO", " Redundant rate of alignment file: %.2f" %
         ((15 - kept) / 15)),
        ("INFO", "Write to BED file"),
        ("INFO", "finished! Check out.bed.")]
    assert proc.stdout == ""


def test_keep_dup_integer_log(macs3_argparser, caplog, se_bed, tmp_path):
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "-f", "BED", "--keep-dup", "1",
                   "--outdir", tmp_path, "-o", "out.bed"])
    msgs = [m for _, m in cmd_messages(caplog)]
    assert "user defined the maximum tags..." in msgs
    assert ("filter out redundant tags at the same location and the same "
            "strand by allowing at most 1 tag(s)") in msgs
    assert " tags after filtering in alignment file: 9" in msgs
    assert " Redundant rate of alignment file: 0.40" in msgs


def test_keep_dup_all_skips_filtering_log(macs3_argparser, caplog, se_bed,
                                          tmp_path):
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "-f", "BED", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert [m for _, m in cmd_messages(caplog)] == [
        "# read tag files...", "# read treatment tags...",
        "# tag size = 36", "# total tags in alignment file: 15",
        "Write to BED file", "finished! Check out.bed."]


# ------------------------------------
# -o/--ofile, --outdir and --dry-run
# ------------------------------------

def test_cli_stdout_default(run_macs3, parse_log, se_bed, tmp_path):
    proc = run_macs3(["filterdup", "-i", se_bed, "--keep-dup", "all"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == expected_bed(SE_READS)
    assert proc.stdout.endswith("\n")
    log = parse_log(proc.stderr)
    assert log[-1] == ("INFO", "finished! Check stdout.")
    assert ("INFO", "Detected format is: BED") in log
    assert sorted(os.listdir(tmp_path)) == ["TMPDIR", "reads.bed"]


def test_cli_ofile_in_outdir(run_macs3, se_bed, tmp_path):
    outdir = tmp_path / "out"
    proc = run_macs3(["filterdup", "-i", se_bed, "--keep-dup", "all",
                      "--outdir", outdir, "-o", "result.bed"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""
    assert os.listdir(outdir) == ["result.bed"]
    assert read_lines(outdir / "result.bed") == expected_bed(SE_READS)


def test_absolute_ofile_ignores_outdir(macs3_argparser, se_bed, tmp_path):
    target = tmp_path / "abs.bed"
    (tmp_path / "other").mkdir()
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "--keep-dup", "all",
                   "--outdir", tmp_path / "other", "-o", target])
    assert read_lines(target) == expected_bed(SE_READS)
    assert os.listdir(tmp_path / "other") == []


def test_stdout_in_process(macs3_argparser, capsys, se_bed, tmp_path):
    run_filterdup(macs3_argparser, ["-i", se_bed, "--keep-dup", "1",
                                    "--outdir", tmp_path])
    assert capsys.readouterr().out.splitlines() == \
        expected_bed(SE_READS, maxdup=1)


def test_cli_dry_run(run_macs3, parse_log, se_bed, tmp_path):
    proc = run_macs3(["filterdup", "-i", se_bed, "-f", "BED", "-g", "1000",
                      "--dry-run", "--outdir", tmp_path, "-o", "dry.bed"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""
    log = [m for _, m in parse_log(proc.stderr)]
    assert log[-5:] == [" max_dup_tags based on binomal = 2",
                        "filter out redundant tags at the same location and "
                        "the same strand by allowing at most 2 tag(s)",
                        " tags after filtering in alignment file: 13",
                        " Redundant rate of alignment file: 0.13",
                        "Dry-run is finished!"]
    assert "Write to BED file" not in log
    # the output file is opened (and left empty) but nothing is written
    dry = tmp_path / "dry.bed"
    assert not dry.exists() or dry.read_text() == ""


def test_dry_run_keep_dup_all(macs3_argparser, caplog, se_bed, tmp_path,
                              capsys):
    run_filterdup(macs3_argparser, ["-i", se_bed, "-f", "BED", "-d",
                                    "--keep-dup", "all",
                                    "--outdir", tmp_path])
    assert [m for _, m in cmd_messages(caplog)] == [
        "# read tag files...", "# read treatment tags...",
        "# tag size = 36", "# total tags in alignment file: 15",
        "Dry-run is finished!"]
    assert capsys.readouterr().out == ""


# ------------------------------------
# -s/--tsize
# ------------------------------------

@pytest.mark.parametrize("tsize", [20, 36, 50])
def test_tsize_sets_output_width(macs3_argparser, caplog, se_bed, tmp_path,
                                 tsize):
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "-f", "BED", "--keep-dup", "all",
                   "-s", tsize, "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == \
        expected_bed(SE_READS, tsize=tsize)
    assert "# tag size = %d" % tsize in [m for _, m in cmd_messages(caplog)]


def test_tsize_detected_as_mean_of_first_ten_reads(macs3_argparser, caplog,
                                                   write_bed, tmp_path):
    # read lengths 30, 31, ..., 41 -> the first ten average 34.5 -> 34
    lines = ["chr1\t%d\t%d\tr\t0\t+" % (100 * i, 100 * i + 30 + i)
             for i in range(12)]
    path = write_bed(lines, name="var.bed")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BED", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert "# tag size = 34" in [m for _, m in cmd_messages(caplog)]
    assert read_lines(tmp_path / "out.bed") == \
        ["chr1\t%d\t%d\t.\t.\t+" % (100 * i, 100 * i + 34) for i in range(12)]


# ------------------------------------
# -f/--format: every single-end format gives the same tags
# ------------------------------------

@pytest.mark.parametrize("fmt", SE_FORMATS)
def test_formats_explicit(macs3_argparser, write_bed, make_alignments,
                          tmp_path, fmt):
    path = write_se(fmt, SE_READS, write_bed, make_alignments)
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", fmt, "--keep-dup", "1",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == \
        expected_bed(SE_READS, maxdup=1)


@pytest.mark.parametrize("fmt", ["BAM", "BED", "ELAND", "ELANDEXPORT"])
def test_formats_auto_detected(macs3_argparser, caplog, write_bed,
                               make_alignments, tmp_path, fmt):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    path = write_se(fmt, SE_READS, write_bed, make_alignments)
    run_filterdup(macs3_argparser,
                  ["-i", path, "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == expected_bed(SE_READS)
    parser_msgs = [MEM_PREFIX.sub("", r.getMessage()) for r in caplog.records
                   if r.name == "MACS3.IO.Parser"]
    assert "Detected format is: %s" % fmt in parser_msgs


@pytest.mark.parametrize("fmt", ["SAM", "BAM"])
def test_sam_bam_plus_strand_reads(macs3_argparser, write_bed,
                                   make_alignments, tmp_path, fmt):
    """Plus-strand reads only (unaffected by the SAM minus-strand bug),
    both with -f and with AUTO."""
    reads = [r for r in SE_READS if r[2] == "+"] + \
        [("chr1", 7000 + 10 * i, "+") for i in range(4)]
    path = write_se(fmt, reads, write_bed, make_alignments)
    for args in (["-f", fmt], []):
        run_filterdup(macs3_argparser,
                      ["-i", path, "--keep-dup", "all", "--outdir", tmp_path,
                       "-o", "out.bed"] + args)
        assert read_lines(tmp_path / "out.bed") == expected_bed(reads)


def test_cli_auto_gzipped_bed(run_macs3, parse_log, write_bed, tmp_path):
    path = write_bed(bed_lines(SE_READS), name="reads.bed.gz")
    proc = run_macs3(["filterdup", "-i", path, "--keep-dup", "all",
                      "--outdir", tmp_path, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    log = parse_log(proc.stderr)
    assert log[:5] == [("INFO", "# read tag files..."),
                       ("INFO", "# read treatment tags..."),
                       ("INFO", "Detected format is: BED"),
                       ("INFO", "* Input file is gzipped."),
                       ("INFO", "# tag size = 36")]
    assert read_lines(tmp_path / "out.bed") == expected_bed(SE_READS)


def test_bed_header_lines_are_skipped(macs3_argparser, write_bed, tmp_path):
    lines = ["track name=reads", "browser position chr1:1-1000",
             "# a comment"] + bed_lines(SE_READS)
    path = write_bed(lines, name="header.bed")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BED", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == expected_bed(SE_READS)


def test_bed_without_strand_column_is_plus(macs3_argparser, write_bed,
                                           tmp_path):
    path = write_bed(["chr1\t100\t136", "chr1\t200\t236\tr2"],
                     name="bed3.bed")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BED", "--keep-dup", "all", "-s", "36",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == [
        "chr1\t100\t136\t.\t.\t+", "chr1\t200\t236\t.\t.\t+"]


def test_sam_bam_flag_filters(macs3_argparser, make_alignments, tmp_path):
    """Unmapped (4), secondary (256), QC-fail (512) and supplementary
    (2048) records are dropped; of a proper pair only read 1 is kept."""
    keep = [dict(name="k%d" % i, ref="chr1", pos=1000 * i, flag=0)
            for i in range(1, 9)]
    drop = [dict(name="u", ref="chr1", pos=50, flag=4),
            dict(name="s", ref="chr1", pos=60, flag=256),
            dict(name="q", ref="chr1", pos=70, flag=512),
            dict(name="x", ref="chr1", pos=80, flag=2048)]
    pair = [dict(name="p", ref="chr1", pos=20000, flag=99, next_ref="chr1",
                 next_pos=20200, tlen=236),
            dict(name="p", ref="chr1", pos=20200, flag=147, next_ref="chr1",
                 next_pos=20000, tlen=-236)]
    path = make_alignments(keep + drop + pair, refs=REFS, name="flags.bam")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BAM", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == \
        ["chr1\t%d\t%d\t.\t.\t+" % (1000 * i, 1000 * i + 36)
         for i in range(1, 9)] + ["chr1\t20000\t20036\t.\t.\t+"]


@pytest.mark.parametrize("fmt", ["BAM"])
def test_minus_strand_cigar_sets_five_prime_end(macs3_argparser,
                                                make_alignments, tmp_path,
                                                fmt):
    """A minus-strand read's 5' end is its start plus the reference length
    of M/D/N/=/X operations (here 20M5D16M -> 41; 2S34M -> 34)."""
    recs = [dict(name="a%d" % i, ref="chr1", pos=100 * i, flag=0)
            for i in range(1, 11)]
    recs += [dict(name="d", ref="chr1", pos=5000, flag=16,
                  cigar="20M5D16M"),
             dict(name="s", ref="chr1", pos=6000, flag=16, cigar="2S34M")]
    path = make_alignments(recs, refs=REFS, name="cigar." + fmt.lower(),
                           fmt=fmt.lower())
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", fmt, "--keep-dup", "all", "-s", "36",
                   "--outdir", tmp_path, "-o", "out.bed"])
    lines = read_lines(tmp_path / "out.bed")
    assert lines[-2:] == ["chr1\t5005\t5041\t.\t.\t-",
                          "chr1\t5998\t6034\t.\t.\t-"]


# ------------------------------------
# paired-end input: BEDPE and BAMPE give BEDPE output
# ------------------------------------

@pytest.mark.parametrize("keepdup, gsize, maxdup", [
    ("all", "hs", None),
    ("1", "hs", 1),
    ("2", "hs", 2),
    ("auto", "2000", binom_cutoff(len(PE_FRAGS), 2000)),
])
@pytest.mark.parametrize("fmt", ["BEDPE", "BAMPE"])
def test_paired_end_formats(macs3_argparser, write_bedpe, make_alignments,
                            make_pe_pair, tmp_path, fmt, keepdup, gsize,
                            maxdup):
    if fmt == "BEDPE":
        path = write_bedpe(PE_FRAGS)
    else:
        path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair)
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", fmt, "--keep-dup", keepdup, "-g", gsize,
                   "--outdir", tmp_path, "-o", "out.bedpe"])
    assert read_lines(tmp_path / "out.bedpe") == \
        expected_bedpe(PE_FRAGS, maxdup)


def test_cli_paired_end_log(run_macs3, parse_log, make_alignments,
                            make_pe_pair, tmp_path):
    path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair)
    proc = run_macs3(["filterdup", "-i", path, "-f", "BAMPE",
                      "--keep-dup", "1", "--outdir", tmp_path,
                      "-o", "out.bedpe"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "# read input file in Paired-end mode."),
        ("INFO", "# read treatment fragments..."),
        ("INFO", "8 fragments have been read."),
        ("INFO", "# total fragments/pairs in alignment file: 8"),
        ("INFO", "user defined the maximum tags..."),
        ("INFO", "filter out redundant tags at the same location and the "
                 "same strand by allowing at most 1 tag(s)"),
        ("INFO", " tags after filtering in alignment file: 5"),
        ("INFO", " Redundant rate of alignment file: 0.38"),
        ("INFO", "Write to BED file"),
        ("INFO", "finished! Check out.bedpe.")]


def test_single_end_mode_of_paired_bam(macs3_argparser, make_alignments,
                                       make_pe_pair, tmp_path):
    """The same pairs read with -f BAM give BED6 tags from read 1 only."""
    path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair)
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BAM", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    reads = [("chr1", s, "+") for _, s, _ in PE_FRAGS]
    assert read_lines(tmp_path / "out.bed") == expected_bed(reads)


# ------------------------------------
# several input files and --buffer-size
# ------------------------------------

def test_multiple_input_files_se(macs3_argparser, write_bed, tmp_path):
    a = write_bed(bed_lines(SE_READS[:6]), name="a.bed")
    b = write_bed(bed_lines(SE_READS[6:]), name="b.bed")
    run_filterdup(macs3_argparser,
                  ["-i", a, b, "-f", "BED", "--keep-dup", "1",
                   "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == \
        expected_bed(SE_READS, maxdup=1)


def test_multiple_input_files_pe(macs3_argparser, write_bedpe, tmp_path):
    a = write_bedpe(PE_FRAGS[:3], name="a.bedpe")
    b = write_bedpe(PE_FRAGS[3:], name="b.bedpe")
    run_filterdup(macs3_argparser,
                  ["-i", a, b, "-f", "BEDPE", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bedpe"])
    assert read_lines(tmp_path / "out.bedpe") == expected_bedpe(PE_FRAGS)


@pytest.mark.parametrize("buffer_size", [1, 3, 100000])
def test_buffer_size_does_not_change_output(macs3_argparser, se_bed,
                                            pe_bedpe, tmp_path, buffer_size):
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "-f", "BED", "--keep-dup", "2",
                   "--buffer-size", buffer_size,
                   "--outdir", tmp_path, "-o", "se.bed"])
    run_filterdup(macs3_argparser,
                  ["-i", pe_bedpe, "-f", "BEDPE", "--keep-dup", "2",
                   "--buffer-size", buffer_size,
                   "--outdir", tmp_path, "-o", "pe.bedpe"])
    assert read_lines(tmp_path / "se.bed") == expected_bed(SE_READS, maxdup=2)
    assert read_lines(tmp_path / "pe.bedpe") == expected_bedpe(PE_FRAGS, 2)


# ------------------------------------
# --verbose
# ------------------------------------

@pytest.mark.parametrize("verbose, n_messages", [(0, 0), (1, 0), (2, 11),
                                                 (3, 11)])
def test_verbose_levels(macs3_argparser, caplog, se_bed, tmp_path, verbose,
                        n_messages):
    run_filterdup(macs3_argparser,
                  ["-i", se_bed, "-f", "BED", "-g", "1000",
                   "--verbose", verbose, "--outdir", tmp_path,
                   "-o", "out.bed"])
    msgs = cmd_messages(caplog)
    assert len(msgs) == n_messages
    assert all(level == "INFO" for level, _ in msgs)
    assert read_lines(tmp_path / "out.bed") == \
        expected_bed(SE_READS, maxdup=2)


# ------------------------------------
# chromosome order depends on the hash seed
# ------------------------------------


def test_chromosomes_are_contiguous_blocks(macs3_argparser, write_bed,
                                           tmp_path):
    reads = [("chr%d" % k, 100 * i, "+-"[i % 2])
             for k in range(1, 6) for i in range(1, 5)]
    path = write_bed(bed_lines(reads), name="multi.bed")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BED", "--keep-dup", "all",
                   "--outdir", tmp_path, "-o", "out.bed"])
    lines = read_lines(tmp_path / "out.bed")
    order = []
    for line in lines:
        chrom = line.split("\t")[0]
        if not order or order[-1] != chrom:
            order.append(chrom)
    assert sorted(order) == ["chr%d" % k for k in range(1, 6)]
    for chrom in order:
        assert [x for x in lines if x.split("\t")[0] == chrom] == \
            expected_bed(reads, chroms=[chrom])


# ------------------------------------
# upstream test data (test/cmdlinetest)
# ------------------------------------

def count_gz_lines(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for line in fh if line.strip())


def test_cli_standard_result_se(run_macs3, parse_log, test_dir, tmp_path):
    """Pins the current output. Compared with
    test/standard_results_filterdup/run_filterdup_result.bed, the
    upstream expected output (one chromosome, so the order is fixed).
    The cutoff is derived from scipy; the output itself has no
    independent derivation at this size."""
    chip = test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz"
    std = test_dir / "standard_results_filterdup" / "run_filterdup_result.bed"
    proc = run_macs3(["filterdup", "-g", "52000000", "-i", chip,
                      "--outdir", tmp_path, "-o", "run_filterdup_result.bed"],
                     timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "run_filterdup_result.bed").read_bytes() == \
        std.read_bytes()
    t0 = count_gz_lines(chip)
    t1 = len(read_lines(std))
    log = [m for _, m in parse_log(proc.stderr)]
    assert "# total tags in alignment file: %d" % t0 in log
    assert " max_dup_tags based on binomal = %d" % \
        binom_cutoff(t0, 52000000) in log
    assert " tags after filtering in alignment file: %d" % t1 in log
    assert " Redundant rate of alignment file: %.2f" % ((t0 - t1) / t0) in log


def test_upstream_se_matches_reference(macs3_argparser, test_dir, tmp_path):
    """Independent check of the run above: every read is 101 bp, so the
    output is each strand's 5' ends sorted and capped at the scipy
    binomial cutoff, written as 101-bp BED6 records."""
    chip = test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz"
    with gzip.open(chip, "rt") as fh:
        rows = [line.split("\t") for line in fh if line.strip()]
    reads = [(c, int(s) if st.strip() == "+" else int(e) - 101,
              st.strip()) for c, s, e, _, _, st in rows]
    maxdup = binom_cutoff(len(reads), 52000000)
    run_filterdup(macs3_argparser,
                  ["-g", "52000000", "-i", chip, "--outdir", tmp_path,
                   "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == \
        expected_bed(reads, tsize=101, maxdup=maxdup, readlen=101)


def test_upstream_bedpe_matches_reference(macs3_argparser, test_dir,
                                          tmp_path):
    bedpe = test_dir / "CTCF_PE_ChIP_chr22_50k.bedpe.gz"
    with gzip.open(bedpe, "rt") as fh:
        frags = [(c, int(s), int(e)) for c, s, e in
                 (line.split() for line in fh if line.strip())]
    maxdup = binom_cutoff(len(frags), 52000000)
    run_filterdup(macs3_argparser,
                  ["-g", "52000000", "-f", "BEDPE", "-i", bedpe,
                   "--outdir", tmp_path, "-o", "out.bedpe"])
    assert read_lines(tmp_path / "out.bedpe") == expected_bedpe(frags, maxdup)


def test_standard_result_pe(macs3_argparser, test_dir, tmp_path):
    """Pins the current output. Compared with
    test/standard_results_filterdup/run_filterdup_result_pe.bedpe from
    upstream (BAMPE, -g 52000000, one chromosome; 50k pairs are too
    many to derive by hand)."""
    std = test_dir / "standard_results_filterdup" / \
        "run_filterdup_result_pe.bedpe"
    run_filterdup(macs3_argparser,
                  ["-g", "52000000", "-f", "BAMPE",
                   "-i", test_dir / "CTCF_PE_ChIP_chr22_50k.bam",
                   "--outdir", tmp_path, "-o", "pe.bedpe"])
    assert (tmp_path / "pe.bedpe").read_bytes() == std.read_bytes()


def test_bedpe_and_bampe_of_upstream_data_agree(macs3_argparser, test_dir,
                                                tmp_path):
    """Pins the current output. The BEDPE copy of the upstream pairs gives
    the BAMPE standard result (50k pairs, too many to derive by hand)."""
    run_filterdup(macs3_argparser,
                  ["-g", "52000000", "-f", "BEDPE",
                   "-i", test_dir / "CTCF_PE_ChIP_chr22_50k.bedpe.gz",
                   "--outdir", tmp_path, "-o", "bedpe.bedpe"])
    std = test_dir / "standard_results_filterdup" / \
        "run_filterdup_result_pe.bedpe"
    assert read_lines(tmp_path / "bedpe.bedpe") == read_lines(std)


def test_contigs50k_buffer_size(run_macs3, test_dir, tmp_path):
    """Pins the current output. 50k contigs with --buffer-size 1000
    against upstream's standard_results_50kcontigs (compared as sorted
    lines because the contig order depends on the hash seed)."""
    std = test_dir / "standard_results_50kcontigs" / "run_filterdup_result.bed"
    proc = run_macs3(["filterdup", "-g", "10000000",
                      "-i", test_dir / "contigs50k.bed.gz",
                      "--outdir", tmp_path, "-o", "out.bed",
                      "--buffer-size", "1000"], timeout=300)
    assert proc.returncode == 0, proc.stderr
    assert sorted(read_lines(tmp_path / "out.bed")) == \
        sorted(read_lines(std))


# ------------------------------------
# cal_max_dup_tags
# ------------------------------------

@pytest.mark.parametrize("gsize, n, p", [
    (1000, 15, 1e-5), (100000, 15, 1e-5), (1000, 15, 0.01),
    (52000000, 50000, 1e-5), (1e6, 200000, 1e-5), (2913022398, 10, 1e-5),
    (5000, 1000, 1e-3)])
def test_cal_max_dup_tags_matches_scipy(gsize, n, p):
    assert cal_max_dup_tags(gsize, n, p) == binom_cutoff(n, gsize, p)


def test_cal_max_dup_tags_default_p():
    assert cal_max_dup_tags(1000, 15) == binom_cutoff(15, 1000, 1e-5) == 2


def test_cal_max_dup_tags_zero_tags():
    assert cal_max_dup_tags(1000, 0) == 0


def test_cal_max_dup_tags_invalid_p():
    with pytest.raises(Exception, match="CDF must >= 0 or <= 1"):
        cal_max_dup_tags(1000, 15, 2.0)


# ------------------------------------
# load_tag_files_options / load_frag_files_options
# ------------------------------------

class _Opts:
    """Minimal options object for the load_* functions."""

    def __init__(self, parser, ifile, tsize=None, buffer_size=100000):
        self.parser = parser
        self.ifile = ifile
        self.tsize = tsize
        self.buffer_size = buffer_size
        self.messages = []
        self.info = self.messages.append


def test_load_tag_files_options_detects_tsize(se_bed):
    opts = _Opts(BEDParser, [se_bed])
    track = load_tag_files_options(opts)
    assert opts.tsize == 36
    assert track.total == 15
    assert opts.messages == ["# read treatment tags..."]
    plus, minus = track.get_locations_by_chr(b"chr1")
    assert list(plus) == sorted(PLUS)
    assert list(minus) == sorted(s + READLEN for s in MINUS)


def test_load_tag_files_options_keeps_given_tsize(se_bed):
    opts = _Opts(BEDParser, [se_bed], tsize=50)
    load_tag_files_options(opts)
    assert opts.tsize == 50


def test_load_tag_files_options_several_files(write_bed):
    a = write_bed(bed_lines(SE_READS[:4]), name="a.bed")
    b = write_bed(bed_lines(SE_READS[4:]), name="b.bed")
    track = load_tag_files_options(_Opts(BEDParser, [a, b]))
    assert track.total == len(SE_READS)


def test_load_frag_files_options(write_bedpe):
    a = write_bedpe(PE_FRAGS[:5], name="a.bedpe")
    b = write_bedpe(PE_FRAGS[5:], name="b.bedpe")
    opts = _Opts(BEDPEParser, [a, b])
    track = load_frag_files_options(opts)
    assert opts.messages == ["# read treatment fragments..."]
    assert track.total == len(PE_FRAGS)
    lengths = [e - s for _, s, e in PE_FRAGS]
    assert track.average_template_length == \
        pytest.approx(sum(lengths) / len(lengths), rel=1e-6)
    locs = track.get_locations_by_chr(b"chr1")
    assert [tuple(x) for x in locs] == sorted((s, e) for _, s, e in PE_FRAGS)


def test_run_returns_none_and_sets_pe_mode(macs3_argparser, pe_bedpe,
                                           tmp_path):
    options = macs3_argparser.parse_args(
        ["filterdup", "-i", pe_bedpe, "-f", "BEDPE", "--keep-dup", "all",
         "--outdir", str(tmp_path), "-o", "x.bedpe"])
    assert run(options) is None
    assert options.PE_MODE is True
    assert options.format == "BEDPE"


# ------------------------------------
# inputs that break the parsers
# ------------------------------------


# ------------------------------------
# opt_validate_filterdup checks that argparse makes unreachable
# ------------------------------------

def test_validator_rejects_unknown_format(macs3_argparser, caplog, se_bed,
                                          tmp_path):
    """-f has fixed choices, so only a direct run() call reaches this."""
    options = macs3_argparser.parse_args(["filterdup", "-i", se_bed,
                                          "--outdir", str(tmp_path)])
    options.format = "bogus"
    with pytest.raises(SystemExit) as exc:
        run(options)
    assert exc.value.code == 1
    assert cmd_messages(caplog) == [
        ("ERROR", "Format \"BOGUS\" cannot be recognized!")]


@pytest.mark.parametrize("fmt", ["bed", "Bed", "auto"])
def test_validator_uppercases_format(macs3_argparser, se_bed, tmp_path, fmt):
    options = macs3_argparser.parse_args(["filterdup", "-i", se_bed,
                                          "--keep-dup", "all",
                                          "--outdir", str(tmp_path),
                                          "-o", "out.bed"])
    options.format = fmt
    run(options)
    assert options.format == fmt.upper()
    assert read_lines(tmp_path / "out.bed") == expected_bed(SE_READS)


# ------------------------------------
# boundaries: one read, one fragment
# ------------------------------------

@pytest.mark.parametrize("keepdup", ["all", "1", "auto"])
def test_single_read(macs3_argparser, write_bed, tmp_path, keepdup):
    path = write_bed(["chr1\t500\t536\tr\t0\t-"], name="one.bed")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BED", "-g", "1000",
                   "--keep-dup", keepdup, "--outdir", tmp_path,
                   "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == ["chr1\t500\t536\t.\t.\t-"]


@pytest.mark.parametrize("keepdup", ["all", "1", "auto"])
def test_single_fragment(macs3_argparser, write_bedpe, tmp_path, keepdup):
    path = write_bedpe([("chr1", 500, 750)], name="one.bedpe")
    run_filterdup(macs3_argparser,
                  ["-i", path, "-f", "BEDPE", "-g", "1000",
                   "--keep-dup", keepdup, "--outdir", tmp_path,
                   "-o", "out.bedpe"])
    assert read_lines(tmp_path / "out.bedpe") == ["chr1\t500\t750"]


def test_positions_near_int32_max(macs3_argparser, write_bed, write_bedpe,
                                  tmp_path):
    top = 2 ** 31 - 1
    path = write_bed(["chr1\t%d\t%d\tr\t0\t+" % (top - 36, top),
                      "chr1\t%d\t%d\tr\t0\t-" % (top - 36, top)],
                     name="top.bed")
    run_filterdup(macs3_argparser, ["-i", path, "-f", "BED",
                                    "--keep-dup", "all",
                                    "--outdir", tmp_path, "-o", "out.bed"])
    assert read_lines(tmp_path / "out.bed") == [
        "chr1\t%d\t%d\t.\t.\t+" % (top - 36, top),
        "chr1\t%d\t%d\t.\t.\t-" % (top - 36, top)]
    pe = write_bedpe([("chr1", top - 300, top)], name="top.bedpe")
    run_filterdup(macs3_argparser, ["-i", pe, "-f", "BEDPE",
                                    "--keep-dup", "all",
                                    "--outdir", tmp_path, "-o", "out.bedpe"])
    assert read_lines(tmp_path / "out.bedpe") == \
        ["chr1\t%d\t%d" % (top - 300, top)]
