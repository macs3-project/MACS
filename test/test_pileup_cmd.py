#!/usr/bin/env python
"""Module Description: Test functions of the pileup subcommand
(MACS3/Commands/pileup_cmd.py), called directly and through
``macs3 pileup``.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import gzip
import logging
import re

import numpy as np
import pytest

from MACS3.Commands.pileup_cmd import (run,
                                       load_tag_files_options,
                                       load_frag_files_options)
from MACS3.IO.Parser import (BEDParser,
                             BEDPEParser,
                             FragParser)

# ------------------------------------
# local helpers
# ------------------------------------

READLEN = 36
REFS = (("chr1", 100000),)
SE_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BOWTIE", "BAM"]
MEM_PREFIX = re.compile(r"^\[\d+ MB\] ")

# (chrom, 0-based leftmost start, strand); every read READLEN long
PLUS = [100, 100, 100, 250, 400, 400, 1000, 5000]
MINUS = [300, 300, 520, 700, 700, 700, 2000]
SE_READS = sorted([("chr1", s, "+") for s in PLUS] +
                  [("chr1", s, "-") for s in MINUS], key=lambda r: r[1])

PE_FRAGS = [("chr1", 100, 300), ("chr1", 100, 300), ("chr1", 150, 400),
            ("chr1", 300, 350), ("chr1", 500, 650), ("chr1", 900, 1000),
            ("chr2", 50, 250), ("chr2", 60, 200)]

# fragments file rows: (chrom, start, end, barcode, count)
FRAG_ROWS = [("chr1", 100, 300, "AAAC", 1), ("chr1", 150, 250, "CCCG", 3),
             ("chr1", 200, 500, "AAAC", 2), ("chr1", 600, 700, "GGGT", 4),
             ("chr2", 50, 120, "CCCG", 1), ("chr2", 80, 200, "GGGT", 2)]


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
             readlen=READLEN, name=None):
    if fmt in ("SAM", "BAM"):
        recs = [dict(name="r%d" % i, ref=c, pos=s,
                     flag=0 if st == "+" else 16, cigar="%dM" % readlen)
                for i, (c, s, st) in enumerate(reads)]
        return make_alignments(recs, refs=refs, fmt=fmt.lower(),
                               name=name or "reads." + fmt.lower())
    return write_bed(TEXT_WRITERS[fmt](reads, readlen),
                     name=name or "reads." + fmt.lower())


def write_bampe(frags, make_alignments, make_pe_pair, refs):
    recs = []
    for i, (c, s, e) in enumerate(frags):
        recs += make_pe_pair("p%d" % i, c, s, e)
    return make_alignments(recs, refs=refs, name="frags.bam")


def ref_bedgraph(chrom, intervals, rlen=None):
    """Reference pileup: weighted intervals (start, end, weight) are
    clipped to [0, rlen], added on a numpy coverage array, and written
    as runs of equal value from 0 to the last covered base."""
    clipped = [(max(s, 0), e if rlen is None else min(e, rlen), w)
               for s, e, w in intervals]
    end = max(e for _, e, _ in clipped)
    cov = np.zeros(end, dtype=float)
    for s, e, w in clipped:
        cov[s:e] += w
    rows, start = [], 0
    for i in range(1, end + 1):
        if i == end or cov[i] != cov[start]:
            rows.append("%s\t%d\t%d\t%.5f" % (chrom, start, i, cov[start]))
            start = i
    return rows


def sweep_bedgraph(chrom, intervals):
    """Reference for large inputs without a dense array: +w at every
    start, -w at every end, swept in position order; runs of equal
    coverage are merged and written from 0 to the last end."""
    events = {}
    for s, e, w in intervals:
        events[s] = events.get(s, 0) + w
        events[e] = events.get(e, 0) - w
    runs, cov, prev = [], 0, 0
    for pos in sorted(events):
        if pos > prev:
            if runs and runs[-1][2] == cov:
                runs[-1][1] = pos
            else:
                runs.append([prev, pos, cov])
        cov += events[pos]
        prev = pos
    return ["%s\t%d\t%d\t%.5f" % (chrom, s, e, v) for s, e, v in runs]


def read_gz_rows(path):
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def se_intervals(reads, chrom, extsize, both=False, readlen=READLEN):
    """Fragments made from single-end reads: plus reads extend
    [5', 5' + extsize), minus reads (5' = end) [5' - extsize, 5');
    with -B every read becomes [5' - extsize, 5' + extsize)."""
    out = []
    for c, s, st in reads:
        if c != chrom:
            continue
        five = s if st == "+" else s + readlen
        if both:
            out.append((five - extsize, five + extsize, 1))
        elif st == "+":
            out.append((five, five + extsize, 1))
        else:
            out.append((five - extsize, five, 1))
    return out


def read_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


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


def run_pileup(macs3_argparser, args):
    options = macs3_argparser.parse_args(["pileup"] + [str(a) for a in args])
    run(options)
    return options


@pytest.fixture
def se_bed(write_bed):
    return write_bed(bed_lines(SE_READS), name="reads.bed")


@pytest.fixture
def pe_bedpe(write_bedpe):
    return write_bedpe(PE_FRAGS, name="frags.bedpe")


@pytest.fixture
def frag_file(write_frag):
    return write_frag(FRAG_ROWS, name="frags.tsv")


# ------------------------------------
# argument parser of `macs3 pileup`
# ------------------------------------

@pytest.mark.parametrize("dest, value", [
    ("outdir", ""), ("format", "AUTO"), ("barcodefile", ""),
    ("maxcount", None), ("bothdirection", False), ("extsize", 200),
    ("buffer_size", 100000), ("verbose", 2)])
def test_parser_defaults(macs3_argparser, dest, value):
    options = macs3_argparser.parse_args(["pileup", "-i", "x", "-o", "y"])
    assert getattr(options, dest) == value


@pytest.mark.parametrize("args, message", [
    (["-o", "y"], "the following arguments are required: -i/--ifile"),
    (["-i", "x"], "the following arguments are required: -o/--ofile"),
    (["-i", "x", "-o", "y", "--extsize", "1.5"],
     "argument --extsize: invalid int value: '1.5'"),
    (["-i", "x", "-o", "y", "--max-count", "a"],
     "argument --max-count: invalid int value: 'a'"),
    (["-i", "x", "-o", "y", "-f", "BEDGRAPH"],
     "argument -f/--format: invalid choice: 'BEDGRAPH'"),
])
def test_parser_errors(macs3_argparser, capsys, args, message):
    with pytest.raises(SystemExit) as exc:
        macs3_argparser.parse_args(["pileup"] + args)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_cli_missing_ofile(run_macs3, se_bed):
    proc = run_macs3(["pileup", "-i", se_bed], timeout=60)
    assert proc.returncode == 2
    assert proc.stderr.startswith("usage: macs3 pileup")


# ------------------------------------
# opt_validate_pileup reached through the CLI
# ------------------------------------

@pytest.mark.parametrize("extsize", ["0", "-5"])
def test_cli_extsize_must_be_positive(run_macs3, parse_log, se_bed, tmp_path,
                                      extsize):
    proc = run_macs3(["pileup", "-i", se_bed, "-f", "BED", "--extsize",
                      extsize, "--outdir", tmp_path, "-o", "out.bdg"],
                     timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [("ERROR", "--extsize must > 0!")]
    assert not (tmp_path / "out.bdg").exists()


def test_cli_negative_max_count(run_macs3, parse_log, frag_file, tmp_path):
    proc = run_macs3(["pileup", "-i", frag_file, "-f", "FRAG",
                      "--max-count", "-1", "--outdir", tmp_path,
                      "-o", "out.bdg"], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "--max-count can't be a negative value")]


def test_negative_max_count_ignored_for_bed(macs3_argparser, se_bed,
                                            tmp_path):
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED",
                                 "--max-count", "-1", "--outdir", tmp_path,
                                 "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))


# ------------------------------------
# --extsize and -B/--both-direction on single-end reads
# ------------------------------------

@pytest.mark.parametrize("extsize", [1, 36, 50, 147, 200, 400, 1000])
def test_extsize(macs3_argparser, se_bed, tmp_path, extsize):
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED",
                                 "--extsize", extsize,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", extsize))


@pytest.mark.parametrize("extsize", [50, 200, 400])
def test_both_direction(macs3_argparser, se_bed, tmp_path, extsize):
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED", "-B",
                                 "--extsize", extsize,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", extsize,
                                          both=True))


def test_exact_bedgraph_text_small_case(macs3_argparser, write_bed,
                                        tmp_path):
    """Hand-derived: plus read at 100 and minus read ending at 180 with
    --extsize 50 cover [100,150) and [130,180)."""
    path = write_bed(["chr1\t100\t136\tr1\t0\t+",
                      "chr1\t144\t180\tr2\t0\t-"], name="two.bed")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED", "--extsize", "50",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert (tmp_path / "out.bdg").read_text() == (
        "chr1\t0\t100\t0.00000\n"
        "chr1\t100\t130\t1.00000\n"
        "chr1\t130\t150\t2.00000\n"
        "chr1\t150\t180\t1.00000\n")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED", "--extsize", "50",
                                 "-B", "--outdir", tmp_path,
                                 "-o", "both.bdg"])
    # -B: [50,150) and [130,230)
    assert (tmp_path / "both.bdg").read_text() == (
        "chr1\t0\t50\t0.00000\n"
        "chr1\t50\t130\t1.00000\n"
        "chr1\t130\t150\t2.00000\n"
        "chr1\t150\t230\t1.00000\n")


def test_extsize_reaching_int32_max(macs3_argparser, write_bed, tmp_path):
    """A text input has no chromosome length, so positions are clipped at
    2**31 - 1; a fragment ending exactly there is kept whole."""
    path = write_bed(["chr1\t100\t136\tr1\t0\t+"], name="one.bed")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED",
                                 "--extsize", str(2 ** 31 - 1 - 100),
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == [
        "chr1\t0\t100\t0.00000", "chr1\t100\t2147483647\t1.00000"]


def test_fragment_starting_at_zero_has_no_zero_row(macs3_argparser,
                                                   write_bed, tmp_path):
    path = write_bed(["chr1\t0\t36\tr1\t0\t+"], name="zero.bed")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED", "--extsize", "10",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == ["chr1\t0\t10\t1.00000"]


def test_minus_read_near_start_is_clipped_at_zero(macs3_argparser, write_bed,
                                                  tmp_path):
    path = write_bed(["chr1\t10\t46\tr1\t0\t-"], name="clip.bed")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED", "--extsize", "100",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == ["chr1\t0\t46\t1.00000"]


def test_bam_pileup_is_clipped_at_chromosome_end(macs3_argparser,
                                                 make_alignments, write_bed,
                                                 tmp_path):
    reads = [("chr1", 100 * i, "+") for i in range(1, 10)] + \
        [("chr1", 5900, "+")]
    refs = (("chr1", 6000),)
    bam = write_se("BAM", reads, write_bed, make_alignments, refs=refs)
    bed = write_se("BED", reads, write_bed, make_alignments, refs=refs)
    run_pileup(macs3_argparser, ["-i", bam, "-f", "BAM", "--extsize", "300",
                                 "--outdir", tmp_path, "-o", "bam.bdg"])
    run_pileup(macs3_argparser, ["-i", bed, "-f", "BED", "--extsize", "300",
                                 "--outdir", tmp_path, "-o", "bed.bdg"])
    ints = se_intervals(reads, "chr1", 300)
    assert read_lines(tmp_path / "bam.bdg") == \
        ref_bedgraph("chr1", ints, rlen=6000)
    # text formats carry no chromosome length, so nothing is clipped
    assert read_lines(tmp_path / "bed.bdg") == ref_bedgraph("chr1", ints)
    assert read_lines(tmp_path / "bed.bdg")[-1] == "chr1\t5900\t6200\t1.00000"
    assert read_lines(tmp_path / "bam.bdg")[-1] == "chr1\t5900\t6000\t1.00000"


def test_cli_log_single_end(run_macs3, parse_log, se_bed, tmp_path):
    """``macs3 pileup`` runs pileup_v2_cmd since upstream 88af6b3 (#734),
    which logs "extend each read downstream by N bps" instead of
    pileup_cmd's "extend each read towards downstream direction with N
    bps"."""
    proc = run_macs3(["pileup", "-i", se_bed, "-f", "BED", "--outdir",
                      tmp_path, "-o", "out.bdg"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "# read alignment files..."),
        ("INFO", "# read tags..."),
        ("INFO", "tag size is determined as 36 bps"),
        ("INFO", "# tag size = 36"),
        ("INFO", "# total tags in alignment file: 15"),
        ("INFO", "# Pileup alignment file, extend each read downstream by "
                 "200 bps"),
        ("INFO", "# Done! Check out.bdg")]
    assert proc.stdout == ""


def test_both_direction_log(macs3_argparser, caplog, se_bed, tmp_path):
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED", "-B",
                                 "--extsize", "77", "--outdir", tmp_path,
                                 "-o", "out.bdg"])
    assert ("INFO", "# Pileup alignment file, extend each read towards "
                    "up/downstream direction with 77 bps") in \
        cmd_messages(caplog)


# ------------------------------------
# -f/--format: every single-end format gives the same pileup
# ------------------------------------

@pytest.mark.parametrize("fmt", SE_FORMATS)
def test_formats_single_end(macs3_argparser, write_bed, make_alignments,
                            tmp_path, fmt):
    path = write_se(fmt, SE_READS, write_bed, make_alignments)
    run_pileup(macs3_argparser, ["-i", path, "-f", fmt, "--extsize", "150",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 150))


def test_sam_plus_strand_reads(macs3_argparser, write_bed, make_alignments,
                               tmp_path):
    """SAM input whose reads are all on the plus strand (the minus-strand
    SAM bug does not apply)."""
    reads = [r for r in SE_READS if r[2] == "+"]
    path = write_se("SAM", reads, write_bed, make_alignments)
    run_pileup(macs3_argparser, ["-i", path, "-f", "SAM", "--extsize", "150",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        ref_bedgraph("chr1", se_intervals(reads, "chr1", 150))


def test_gzipped_bed(macs3_argparser, write_bed, tmp_path):
    path = write_bed(bed_lines(SE_READS), name="reads.bed.gz")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))


def test_single_end_chromosome_order(macs3_argparser, write_bed,
                                     make_alignments, tmp_path):
    """Text input: chromosomes in order of first appearance in the file;
    BAM input: chromosomes sorted by name."""
    reads = [("chrB", 100, "+"), ("chrB", 400, "-"), ("chrA", 50, "+"),
             ("chrA", 900, "-"), ("chrB", 700, "+"), ("chrA", 300, "+"),
             ("chrB", 1000, "+"), ("chrA", 600, "-"), ("chrB", 1200, "-"),
             ("chrA", 1500, "+")]
    refs = (("chrB", 10000), ("chrA", 10000))
    bed = write_se("BED", reads, write_bed, make_alignments, refs=refs)
    bam = write_se("BAM", reads, write_bed, make_alignments, refs=refs)
    run_pileup(macs3_argparser, ["-i", bed, "-f", "BED", "--extsize", "100",
                                 "--outdir", tmp_path, "-o", "bed.bdg"])
    run_pileup(macs3_argparser, ["-i", bam, "-f", "BAM", "--extsize", "100",
                                 "--outdir", tmp_path, "-o", "bam.bdg"])
    block_a = ref_bedgraph("chrA", se_intervals(reads, "chrA", 100))
    block_b = ref_bedgraph("chrB", se_intervals(reads, "chrB", 100))
    assert read_lines(tmp_path / "bed.bdg") == block_b + block_a
    assert read_lines(tmp_path / "bam.bdg") == block_a + block_b


# ------------------------------------
# paired-end formats: BEDPE, BAMPE, FRAG
# ------------------------------------

def pe_expected(frags):
    out = []
    for chrom in sorted(set(f[0] for f in frags)):
        out += ref_bedgraph(chrom, [(s, e, 1) for c, s, e in frags
                                    if c == chrom])
    return out


@pytest.mark.parametrize("fmt", ["BEDPE", "BAMPE"])
def test_paired_end_formats(macs3_argparser, write_bedpe, make_alignments,
                            make_pe_pair, tmp_path, fmt):
    if fmt == "BEDPE":
        path = write_bedpe(PE_FRAGS)
    else:
        path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair,
                           refs=(("chr2", 100000), ("chr1", 100000)))
    run_pileup(macs3_argparser, ["-i", path, "-f", fmt,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == pe_expected(PE_FRAGS)


def test_paired_end_ignores_extsize_and_both_direction(macs3_argparser,
                                                       pe_bedpe, tmp_path):
    run_pileup(macs3_argparser, ["-i", pe_bedpe, "-f", "BEDPE", "-B",
                                 "--extsize", "17", "--outdir", tmp_path,
                                 "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == pe_expected(PE_FRAGS)


def test_paired_end_exact_text(macs3_argparser, write_bedpe, tmp_path):
    path = write_bedpe([("chr1", 10, 40), ("chr1", 20, 40),
                        ("chr1", 40, 60)], name="tiny.bedpe")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BEDPE",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert (tmp_path / "out.bdg").read_text() == (
        "chr1\t0\t10\t0.00000\n"
        "chr1\t10\t20\t1.00000\n"
        "chr1\t20\t40\t2.00000\n"
        "chr1\t40\t60\t1.00000\n")


def test_cli_log_paired_end(run_macs3, parse_log, pe_bedpe, tmp_path):
    """``macs3 pileup`` runs pileup_v2_cmd since upstream 88af6b3 (#734),
    which logs "# Pileup paired-end alignment file with PileupV2."."""
    proc = run_macs3(["pileup", "-i", pe_bedpe, "-f", "BEDPE",
                      "--outdir", tmp_path, "-o", "out.bdg"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "# read alignment files..."),
        ("INFO", "# read input file in Paired-end mode."),
        ("INFO", "# read fragments..."),
        ("INFO", "# total fragments/pairs in alignment file: 8"),
        ("INFO", "# Pileup paired-end alignment file with PileupV2."),
        ("INFO", "# Done! Check out.bdg")]


def frag_expected(rows, max_count=None, barcodes=None):
    out = []
    for chrom in sorted(set(r[0] for r in rows)):
        ints = [(s, e, min(n, max_count) if max_count else n)
                for c, s, e, b, n in rows
                if c == chrom and (barcodes is None or b in barcodes)]
        if ints:
            out += ref_bedgraph(chrom, ints)
    return out


def test_frag_counts_weight_the_pileup(macs3_argparser, frag_file, tmp_path):
    run_pileup(macs3_argparser, ["-i", frag_file, "-f", "FRAG",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == frag_expected(FRAG_ROWS)


@pytest.mark.parametrize("max_count", [0, 1, 2, 3, 100])
def test_frag_max_count(macs3_argparser, frag_file, tmp_path, max_count):
    run_pileup(macs3_argparser, ["-i", frag_file, "-f", "FRAG",
                                 "--max-count", max_count,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        frag_expected(FRAG_ROWS, max_count=max_count)


def test_frag_max_count_1_equals_bedpe(macs3_argparser, frag_file,
                                       write_bedpe, tmp_path):
    bedpe = write_bedpe([r[:3] for r in FRAG_ROWS], name="same.bedpe")
    run_pileup(macs3_argparser, ["-i", frag_file, "-f", "FRAG",
                                 "--max-count", "1",
                                 "--outdir", tmp_path, "-o", "frag.bdg"])
    run_pileup(macs3_argparser, ["-i", bedpe, "-f", "BEDPE",
                                 "--outdir", tmp_path, "-o", "bedpe.bdg"])
    assert (tmp_path / "frag.bdg").read_text() == \
        (tmp_path / "bedpe.bdg").read_text()


@pytest.mark.parametrize("barcodes", [["AAAC"], ["CCCG", "GGGT"],
                                      ["AAAC", "NOTTHERE"]])
def test_frag_barcodes(macs3_argparser, caplog, frag_file, tmp_path,
                       barcodes):
    bc = tmp_path / "barcodes.txt"
    bc.write_text("".join(b + "\n" for b in barcodes))
    run_pileup(macs3_argparser, ["-i", frag_file, "-f", "FRAG",
                                 "--barcodes", bc,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        frag_expected(FRAG_ROWS, barcodes=set(barcodes))
    n = sum(r[4] for r in FRAG_ROWS if r[3] in barcodes)
    msgs = [m for _, m in cmd_messages(caplog)]
    assert "# total fragments/pairs in alignment file: %d" % \
        sum(r[4] for r in FRAG_ROWS) in msgs
    assert "# extract fragments with given barcodes" in msgs
    assert "#   extracted %d fragments" % n in msgs


def test_frag_barcodes_and_max_count(macs3_argparser, frag_file, tmp_path):
    bc = tmp_path / "barcodes.txt"
    bc.write_text("CCCG\nGGGT\n")
    run_pileup(macs3_argparser, ["-i", frag_file, "-f", "FRAG",
                                 "--barcodes", bc, "--max-count", "2",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        frag_expected(FRAG_ROWS, max_count=2, barcodes={"CCCG", "GGGT"})


def test_barcodes_ignored_for_bedpe(macs3_argparser, pe_bedpe, tmp_path):
    bc = tmp_path / "barcodes.txt"
    bc.write_text("AAAC\n")
    run_pileup(macs3_argparser, ["-i", pe_bedpe, "-f", "BEDPE",
                                 "--barcodes", bc, "--max-count", "1",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == pe_expected(PE_FRAGS)


def test_frag_multiple_files_with_max_count(macs3_argparser, write_frag,
                                            tmp_path):
    a = write_frag(FRAG_ROWS[:3], name="a.tsv")
    b = write_frag(FRAG_ROWS[3:], name="b.tsv")
    run_pileup(macs3_argparser, ["-i", a, b, "-f", "FRAG", "--max-count", "2",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == \
        frag_expected(FRAG_ROWS, max_count=2)


# ------------------------------------
# output file handling, several inputs, --buffer-size, --verbose
# ------------------------------------

def test_existing_output_is_replaced(macs3_argparser, caplog, se_bed,
                                     tmp_path):
    out = tmp_path / "out.bdg"
    out.write_text("old content\n")
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(out) == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))
    assert ("INFO", "# Existing file %s will be replaced!" %
            str(tmp_path / "out.bdg")) in cmd_messages(caplog)


def test_cli_outdir_created(run_macs3, se_bed, tmp_path):
    proc = run_macs3(["pileup", "-i", se_bed, "-f", "BED",
                      "--outdir", tmp_path / "a" / "b", "-o", "p.bdg"],
                     timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert read_lines(tmp_path / "a" / "b" / "p.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))


def test_multiple_input_files(macs3_argparser, write_bed, write_bedpe,
                              tmp_path):
    a = write_bed(bed_lines(SE_READS[:7]), name="a.bed")
    b = write_bed(bed_lines(SE_READS[7:]), name="b.bed")
    run_pileup(macs3_argparser, ["-i", a, b, "-f", "BED",
                                 "--outdir", tmp_path, "-o", "se.bdg"])
    assert read_lines(tmp_path / "se.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))
    c = write_bedpe(PE_FRAGS[:3], name="c.bedpe")
    d = write_bedpe(PE_FRAGS[3:], name="d.bedpe")
    run_pileup(macs3_argparser, ["-i", c, d, "-f", "BEDPE",
                                 "--outdir", tmp_path, "-o", "pe.bdg"])
    assert read_lines(tmp_path / "pe.bdg") == pe_expected(PE_FRAGS)


@pytest.mark.parametrize("buffer_size", [1, 4])
def test_buffer_size_does_not_change_output(macs3_argparser, se_bed,
                                            frag_file, tmp_path, buffer_size):
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED",
                                 "--buffer-size", buffer_size,
                                 "--outdir", tmp_path, "-o", "se.bdg"])
    run_pileup(macs3_argparser, ["-i", frag_file, "-f", "FRAG",
                                 "--buffer-size", buffer_size,
                                 "--outdir", tmp_path, "-o", "frag.bdg"])
    assert read_lines(tmp_path / "se.bdg") == \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))
    assert read_lines(tmp_path / "frag.bdg") == frag_expected(FRAG_ROWS)


@pytest.mark.parametrize("verbose, n_messages", [(0, 0), (1, 0), (2, 7),
                                                 (3, 7)])
def test_verbose_levels(macs3_argparser, caplog, se_bed, tmp_path, verbose,
                        n_messages):
    run_pileup(macs3_argparser, ["-i", se_bed, "-f", "BED",
                                 "--verbose", verbose,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert len(cmd_messages(caplog)) == n_messages


# ------------------------------------
# upstream test data (test/cmdlinetest)
# ------------------------------------

@pytest.mark.parametrize("fmt, infile, std", [
    ("BED", "CTCF_SE_ChIP_chr22_50k.bed.gz", "run_pileup_ChIP.bed.bdg"),
    ("BED", "CTCF_SE_CTRL_chr22_50k.bed.gz", "run_pileup_CTRL.bed.bdg"),
    ("BAMPE", "CTCF_PE_ChIP_chr22_50k.bam", "run_pileup_ChIPPE.bampe.bdg"),
    ("BAMPE", "CTCF_PE_CTRL_chr22_50k.bam", "run_pileup_CTRLPE.bampe.bdg"),
    ("BEDPE", "CTCF_PE_ChIP_chr22_50k.bedpe.gz",
     "run_pileup_ChIPPE.bedpe.bdg"),
    ("BEDPE", "CTCF_PE_CTRL_chr22_50k.bedpe.gz",
     "run_pileup_CTRLPE.bedpe.bdg"),
])
def test_standard_results(macs3_argparser, test_dir, tmp_path, fmt, infile,
                          std):
    """Pins the current output. Byte comparison with upstream's
    test/standard_results_pileup files (50k reads; no independent
    derivation at this size beyond the small cases above)."""
    args = ["-f", fmt, "-i", test_dir / infile, "--outdir", tmp_path,
            "-o", "out.bdg"]
    if fmt == "BED":
        args += ["--extsize", "200"]
    run_pileup(macs3_argparser, args)
    assert (tmp_path / "out.bdg").read_bytes() == \
        (test_dir / "standard_results_pileup" / std).read_bytes()


def test_cli_standard_result_bed(run_macs3, test_dir, tmp_path):
    """Pins the current output. Upstream standard_results_pileup."""
    proc = run_macs3(["pileup", "-f", "BED",
                      "-i", test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz",
                      "--extsize", "200", "--outdir", tmp_path,
                      "-o", "run_pileup_ChIP.bed.bdg"], timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "run_pileup_ChIP.bed.bdg").read_bytes() == \
        (test_dir / "standard_results_pileup" /
         "run_pileup_ChIP.bed.bdg").read_bytes()


def test_sweep_reference_agrees_with_dense_reference():
    ints = [(100, 300, 1), (150, 250, 3), (250, 400, 2), (0, 50, 1),
            (600, 700, 4), (400, 450, 1)]
    assert sweep_bedgraph("c", ints) == ref_bedgraph("c", ints)


@pytest.mark.parametrize("infile", ["CTCF_SE_ChIP_chr22_50k.bed.gz",
                                    "CTCF_SE_CTRL_chr22_50k.bed.gz"])
def test_upstream_bed_matches_reference(macs3_argparser, test_dir, tmp_path,
                                        infile):
    """The upstream single-end runs (--extsize 200) against the sweep
    reference built from the BED records: an independent check of the
    standard_results_pileup files compared above."""
    rows = read_gz_rows(test_dir / infile)
    ints = [(int(s), int(s) + 200, 1) if st == "+" else
            (max(int(e) - 200, 0), int(e), 1)
            for _, s, e, _, _, st in rows]
    run_pileup(macs3_argparser, ["-f", "BED", "-i", test_dir / infile,
                                 "--extsize", "200", "--outdir", tmp_path,
                                 "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == sweep_bedgraph("chr22", ints)


@pytest.mark.parametrize("infile", ["CTCF_PE_ChIP_chr22_50k.bedpe.gz",
                                    "CTCF_PE_CTRL_chr22_50k.bedpe.gz"])
def test_upstream_bedpe_matches_reference(macs3_argparser, test_dir,
                                          tmp_path, infile):
    rows = read_gz_rows(test_dir / infile)
    ints = [(int(s), int(e), 1) for _, s, e in rows]
    run_pileup(macs3_argparser, ["-f", "BEDPE", "-i", test_dir / infile,
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == sweep_bedgraph("chr22", ints)


@pytest.mark.parametrize("args", [[], ["--max-count", "1"],
                                  ["--barcodes", "barcodes.txt"],
                                  ["--barcodes", "barcodes.txt",
                                   "--max-count", "2"]])
def test_upstream_fragments_match_reference(macs3_argparser, test_dir,
                                            tmp_path, args):
    """test/test.fragments.tsv.gz (57k fragments on chr22) with counts,
    --max-count and the 50 barcodes of test/barcodes.txt."""
    rows = read_gz_rows(test_dir / "test.fragments.tsv.gz")
    maxc = int(args[args.index("--max-count") + 1]) \
        if "--max-count" in args else None
    keep = None
    if "--barcodes" in args:
        keep = set((test_dir / "barcodes.txt").read_text().split())
        args = [str(test_dir / a) if a == "barcodes.txt" else a
                for a in args]
    ints = [(int(s), int(e), min(int(n), maxc) if maxc else int(n))
            for _, s, e, bc, n in rows if keep is None or bc in keep]
    run_pileup(macs3_argparser, ["-f", "FRAG", "-i",
                                 test_dir / "test.fragments.tsv.gz"] + args +
               ["--outdir", tmp_path, "-o", "out.bdg"])
    assert read_lines(tmp_path / "out.bdg") == sweep_bedgraph("chr22", ints)


def test_contigs50k_buffer_size(run_macs3, test_dir, tmp_path):
    """Pins the current output. 50k contigs with --buffer-size 1000
    compared with upstream's standard_results_50kcontigs."""
    proc = run_macs3(["pileup", "-f", "BED", "-i",
                      test_dir / "contigs50k.bed.gz", "--extsize", "200",
                      "--outdir", tmp_path, "-o", "out.bdg",
                      "--buffer-size", "1000"], timeout=300)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out.bdg").read_bytes() == \
        (test_dir / "standard_results_50kcontigs" /
         "run_pileup_ChIP.bed.bdg").read_bytes()


# ------------------------------------
# load_tag_files_options / load_frag_files_options
# ------------------------------------

class _Opts:
    def __init__(self, parser, ifile, fmt="BED", maxcount=None,
                 buffer_size=100000):
        self.parser = parser
        self.ifile = ifile
        self.format = fmt
        self.maxcount = maxcount
        self.buffer_size = buffer_size
        self.messages = []
        self.info = self.messages.append


def test_load_tag_files_options_returns_tsize_and_track(write_bed):
    a = write_bed(bed_lines(SE_READS[:9]), name="a.bed")
    b = write_bed(bed_lines(SE_READS[9:], readlen=50), name="b.bed")
    opts = _Opts(BEDParser, [a, b])
    tsize, track = load_tag_files_options(opts)
    # the tag size comes from the first file only
    assert tsize == 36
    assert track.total == len(SE_READS)
    assert opts.messages == ["# read tags...",
                             "tag size is determined as 36 bps"]


def test_load_frag_files_options_bedpe(write_bedpe):
    opts = _Opts(BEDPEParser, [write_bedpe(PE_FRAGS)], fmt="BEDPE",
                 maxcount=1)
    track = load_frag_files_options(opts)
    assert track.total == len(PE_FRAGS)
    assert opts.messages == ["# read fragments..."]


@pytest.mark.parametrize("maxcount, total", [(None, 13), (0, 13), (1, 6),
                                             (2, 10)])
def test_load_frag_files_options_frag(write_frag, maxcount, total):
    a = write_frag(FRAG_ROWS[:2], name="a.tsv")
    b = write_frag(FRAG_ROWS[2:], name="b.tsv")
    opts = _Opts(FragParser, [a, b], fmt="FRAG", maxcount=maxcount)
    track = load_frag_files_options(opts)
    # the total of a fragments track is the sum of (capped) counts
    assert sum(min(r[4], maxcount) if maxcount else r[4]
               for r in FRAG_ROWS) == total
    assert track.total == total


# ------------------------------------
# opt_validate_pileup checks that argparse makes unreachable
# ------------------------------------

def test_validator_rejects_unknown_format(macs3_argparser, caplog, se_bed,
                                          tmp_path):
    """-f has fixed choices, so only a direct run() call reaches this."""
    options = macs3_argparser.parse_args(["pileup", "-i", se_bed,
                                          "--outdir", str(tmp_path),
                                          "-o", "out.bdg"])
    options.format = "bogus"
    with pytest.raises(SystemExit) as exc:
        run(options)
    assert exc.value.code == 1
    assert cmd_messages(caplog) == [
        ("ERROR", "Format \"BOGUS\" cannot be recognized!")]
    assert not (tmp_path / "out.bdg").exists()


@pytest.mark.parametrize("fmt", ["bed", "frag"])
def test_validator_uppercases_format(macs3_argparser, se_bed, frag_file,
                                     tmp_path, fmt):
    path = frag_file if fmt == "frag" else se_bed
    options = macs3_argparser.parse_args(["pileup", "-i", path,
                                          "--outdir", str(tmp_path),
                                          "-o", "out.bdg"])
    options.format = fmt
    run(options)
    assert options.format == fmt.upper()
    expect = frag_expected(FRAG_ROWS) if fmt == "frag" else \
        ref_bedgraph("chr1", se_intervals(SE_READS, "chr1", 200))
    assert read_lines(tmp_path / "out.bdg") == expect


# ------------------------------------
# boundaries: empty input
# ------------------------------------

def test_empty_bed_gives_empty_bedgraph(macs3_argparser, caplog, tmp_path):
    path = tmp_path / "empty.bed"
    path.write_text("")
    run_pileup(macs3_argparser, ["-i", path, "-f", "BED",
                                 "--outdir", tmp_path, "-o", "out.bdg"])
    assert (tmp_path / "out.bdg").read_text() == ""
    assert ("INFO", "# total tags in alignment file: 0") in \
        cmd_messages(caplog)


