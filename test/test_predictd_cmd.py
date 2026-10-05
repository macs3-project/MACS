#!/usr/bin/env python
"""Module Description: Test functions of the predictd subcommand
(MACS3/Commands/predictd_cmd.py), called directly and through
``macs3 predictd``.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import logging
import re

import numpy as np
import pytest

from MACS3.Commands.predictd_cmd import (run,
                                         load_tag_files_options,
                                         load_frag_files_options)
from MACS3.IO.Parser import (BEDParser,
                             BEDPEParser)

# ------------------------------------
# local helpers
# ------------------------------------

READLEN = 36
REFS = (("chr1", 2000000),)
SE_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BOWTIE", "BAM"]
AUTO_FORMATS = ["BED", "ELAND", "ELANDEXPORT", "BAM"]
MEM_PREFIX = re.compile(r"^\[\d+ MB\] ")
GSIZE = "900000"


def model_reads(n_sites=150, frag=200, ntags=10, first=10000, spacing=5000,
                chrom="chr1"):
    """Synthetic ChIP reads with a known strand shift.

    At each of ``n_sites`` sites, ``ntags`` plus-strand reads start at
    base, base+1, ... and ``ntags`` minus-strand reads end ``frag`` bp
    further (5' ends base+frag, base+frag+1, ...), i.e. every read pair
    brackets a ``frag``-bp fragment.
    """
    reads = []
    for i in range(n_sites):
        base = first + i * spacing
        for j in range(ntags):
            reads.append((chrom, base + j, "+"))
            reads.append((chrom, base + frag + j - READLEN, "-"))
    return reads


def tags_bound(total, fold, bw, gsize):
    """min_tags/max_tags of the model: total * fold * peaksize / gsize / 2
    with peaksize = 2 * bw, rounded."""
    return int(round(float(total) * fold * 2 * bw / gsize / 2))


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


def cmd_messages(caplog):
    return [(r.levelname, MEM_PREFIX.sub("", r.getMessage()))
            for r in caplog.records
            if r.name == "MACS3.Utilities.OptValidator"]


def predicted_d(messages):
    for _, m in messages:
        hit = re.match(r"# predicted fragment length is (-?\d+) bps$", m)
        if hit:
            return int(hit.group(1))
    return None


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


def run_predictd(macs3_argparser, args):
    options = macs3_argparser.parse_args(["predictd"] +
                                         [str(a) for a in args])
    run(options)
    return options


def not_enough_pairs(n):
    return [
        ("INFO", "#2 Total number of paired peaks: %d" % n),
        ("WARNING", "#2 MACS3 needs at least 100 paired peaks at + and - "
                    "strand to build the model, but can only find %d! Please "
                    "make your MFOLD range broader and try again. If MACS3 "
                    "still can't build the model, we suggest to use "
                    "--nomodel and --extsize 147 or other fixed number "
                    "instead." % n),
        ("WARNING", "#2 Process for pairing-model is terminated!"),
        ("WARNING", "# Can't find enough pairs of symmetric peaks to build "
                    "model!")]


@pytest.fixture
def model_bed(write_bed):
    return write_bed(bed_lines(model_reads()), name="model.bed")


PE_FRAGS = [("chr1", 100, 300), ("chr1", 1000, 1150), ("chr2", 50, 301)]


# ------------------------------------
# argument parser of `macs3 predictd`
# ------------------------------------

@pytest.mark.parametrize("dest, value", [
    ("format", "AUTO"), ("gsize", "hs"), ("tsize", None), ("bw", 300),
    ("d_min", 20), ("mfold", [5, 50]), ("outdir", ""),
    ("rfile", "predictd_model.R"), ("buffer_size", 100000), ("verbose", 2)])
def test_parser_defaults(macs3_argparser, dest, value):
    options = macs3_argparser.parse_args(["predictd", "-i", "x"])
    assert getattr(options, dest) == value


@pytest.mark.parametrize("args, message", [
    ([], "the following arguments are required: -i/--ifile"),
    (["-i", "x", "-m", "5"], "argument -m/--mfold: expected 2 arguments"),
    (["-i", "x", "-m", "5", "x"],
     "argument -m/--mfold: invalid int value: 'x'"),
    (["-i", "x", "--bw", "1.5"], "argument --bw: invalid int value: '1.5'"),
    (["-i", "x", "--d-min", "a"], "argument --d-min: invalid int value: 'a'"),
    (["-i", "x", "-f", "FRAG"],
     "argument -f/--format: invalid choice: 'FRAG'"),
])
def test_parser_errors(macs3_argparser, capsys, args, message):
    with pytest.raises(SystemExit) as exc:
        macs3_argparser.parse_args(["predictd"] + args)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_cli_argparse_error(run_macs3):
    proc = run_macs3(["predictd", "-i"], timeout=60)
    assert proc.returncode == 2
    assert proc.stderr.startswith("usage: macs3 predictd")


# ------------------------------------
# opt_validate_predictd reached through the CLI
# ------------------------------------

def test_cli_invalid_gsize(run_macs3, parse_log, model_bed, tmp_path):
    proc = run_macs3(["predictd", "-i", model_bed, "-g", "human",
                      "--outdir", tmp_path], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "Error when interpreting --gsize option: human"),
        ("ERROR", "Available shortcuts of effective genome sizes are "
                  "hs,mm,ce,dm")]


def test_cli_mfold_lower_above_upper(run_macs3, parse_log, model_bed,
                                     tmp_path):
    """The message is %-formatted with the mfold list, which Python
    treats as a mapping, so it prints unchanged."""
    proc = run_macs3(["predictd", "-i", model_bed, "-m", "50", "5",
                      "--outdir", tmp_path], timeout=60)
    assert proc.returncode == 1
    assert parse_log(proc.stderr) == [
        ("ERROR", "Upper limit of mfold should be greater than lower "
                  "limit!")]


@pytest.mark.parametrize("args", [["--d-min", "-1"], ["-m", "50", "5"]])
def test_cli_invalid_d_min_or_mfold_exit_code(run_macs3, model_bed, tmp_path,
                                              args):
    """The checks do stop the run (non-zero exit, no R script)."""
    proc = run_macs3(["predictd", "-i", model_bed, "--outdir", tmp_path] +
                     args, timeout=60)
    assert proc.returncode != 0
    assert not (tmp_path / "predictd_model.R").exists()


# ------------------------------------
# the predicted d on reads with a known strand shift
# ------------------------------------

def test_cli_predicted_d_log(run_macs3, parse_log, model_bed, tmp_path):
    proc = run_macs3(["predictd", "-i", model_bed, "-f", "BED", "-g", GSIZE,
                      "--outdir", tmp_path, "--rfile", "r.R"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    log = parse_log(proc.stderr)
    assert log[:9] == [
        ("INFO", "# read alignment files..."),
        ("INFO", "# read treatment tags..."),
        ("INFO", "tag size is determined as 36 bps"),
        ("INFO", "# tag size = 36"),
        ("INFO", "# total tags in alignment file: 3000"),
        ("INFO", "# Build Peak Model..."),
        ("INFO", "#2 looking for paired plus/minus strand peaks..."),
        ("INFO", "#2 Total number of paired peaks: 150"),
        ("INFO", "#2 Model building with cross-correlation: Done")]
    d = predicted_d(log)
    assert abs(d - 200) <= 1
    assert log[9:] == [
        ("INFO", "# finished!"),
        ("INFO", "# predicted fragment length is %d bps" % d),
        ("INFO", "# alternative fragment length(s) may be %d bps" % d),
        ("INFO", "# Generate R script for model : %s" % (tmp_path / "r.R"))]
    assert (tmp_path / "r.R").exists()


@pytest.mark.parametrize("frag", [150, 200, 250])
def test_predicted_d_within_one_bp(macs3_argparser, caplog, write_bed,
                                   tmp_path, frag):
    path = write_bed(bed_lines(model_reads(frag=frag)), name="m.bed")
    run_predictd(macs3_argparser, ["-i", path, "-g", GSIZE,
                                   "--outdir", tmp_path])
    d = predicted_d(cmd_messages(caplog))
    assert d is not None and abs(d - frag) <= 1


def test_debug_messages_show_model_bounds(macs3_argparser, caplog,
                                          model_bed, tmp_path):
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                   "--verbose", "3", "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    lo = tags_bound(3000, 5, 300, 900000)
    hi = tags_bound(3000, 50, 300, 900000)
    assert (lo, hi) == (5, 50)
    assert ("DEBUG", "#2 min_tags: %d; max_tags:%d; " % (lo, hi)) in msgs
    assert ("DEBUG", "Number of unique tags on + strand: 1500") in msgs
    assert ("DEBUG", "Number of peaks in + strand: 150") in msgs
    assert ("DEBUG", "Number of peaks in - strand: 150") in msgs
    assert ("DEBUG", "Number of paired peaks in this chromosome: 150") in msgs
    assert ("DEBUG", "#  Summary Model:") in msgs
    assert ("DEBUG", "#   min_tags: %d" % lo) in msgs
    d = predicted_d(msgs)
    assert ("DEBUG", "#   d: %d" % d) in msgs


@pytest.mark.parametrize("verbose, has_info, has_debug", [
    (0, False, False), (1, False, False), (2, True, False), (3, True, True)])
def test_verbose_levels(macs3_argparser, caplog, model_bed, tmp_path,
                        verbose, has_info, has_debug):
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                   "--verbose", verbose,
                                   "--outdir", tmp_path])
    levels = set(level for level, _ in cmd_messages(caplog))
    assert ("INFO" in levels) == has_info
    assert ("DEBUG" in levels) == has_debug
    # the R script is written whatever the verbosity
    assert (tmp_path / "predictd_model.R").exists()


# ------------------------------------
# -g, -m/--mfold, --bw and the minimum of 100 paired peaks
# ------------------------------------

def test_default_gsize_finds_no_peaks(macs3_argparser, caplog, model_bed,
                                      tmp_path):
    """With -g hs both bounds round to 0 for 3000 tags, so no strand peak
    passes 'height < max_tags'."""
    assert tags_bound(3000, 50, 300, 2913022398) == 0
    run_predictd(macs3_argparser, ["-i", model_bed, "--outdir", tmp_path])
    assert cmd_messages(caplog)[-4:] == not_enough_pairs(0)
    assert not (tmp_path / "predictd_model.R").exists()


@pytest.mark.parametrize("mfold, ok", [
    (["5", "50"], True),      # bounds 5, 50 around the peak height 10
    (["2", "11"], True),      # bounds 2, 11
    (["10", "50"], False),    # lower bound 10: height must exceed 10
    (["2", "10"], False),     # upper bound 10: height must be below 10
])
def test_mfold_bounds(macs3_argparser, caplog, model_bed, tmp_path, mfold,
                      ok):
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE, "-m"] +
                 mfold + ["--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    if ok:
        assert ("INFO", "#2 Total number of paired peaks: 150") in msgs
        assert abs(predicted_d(msgs) - 200) <= 1
    else:
        assert msgs[-4:] == not_enough_pairs(0)


@pytest.mark.parametrize("n_sites, ok", [(99, False), (100, True)])
def test_at_least_100_paired_peaks(macs3_argparser, caplog, write_bed,
                                   tmp_path, n_sites, ok):
    path = write_bed(bed_lines(model_reads(n_sites=n_sites)), name="m.bed")
    run_predictd(macs3_argparser, ["-i", path, "-g", GSIZE,
                                   "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    if ok:
        assert ("INFO", "#2 Total number of paired peaks: 100") in msgs
        assert abs(predicted_d(msgs) - 200) <= 1
    else:
        assert msgs[-4:] == not_enough_pairs(99)


def test_cli_not_enough_pairs_exit_code(run_macs3, parse_log, write_bed,
                                        tmp_path):
    path = write_bed(bed_lines(model_reads(n_sites=99)), name="m.bed")
    proc = run_macs3(["predictd", "-i", path, "-g", GSIZE,
                      "--outdir", tmp_path], timeout=60)
    assert proc.returncode == 0
    assert parse_log(proc.stderr)[-4:] == not_enough_pairs(99)
    assert not (tmp_path / "predictd_model.R").exists()


@pytest.mark.parametrize("bw, ok", [(50, False), (200, True), (300, True)])
def test_bw(macs3_argparser, caplog, model_bed, tmp_path, bw, ok):
    """--bw 50: strand peaks narrower than the 200-bp minimum length, so
    no pairs; --bw 200 and 300 find all 150 pairs."""
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE, "--bw", bw,
                                   "--verbose", "3", "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    lo = tags_bound(3000, 5, bw, 900000)
    hi = tags_bound(3000, 50, bw, 900000)
    assert ("DEBUG", "#2 min_tags: %d; max_tags:%d; " % (lo, hi)) in msgs
    if ok:
        assert abs(predicted_d(msgs) - 200) <= 1
        rscript = (tmp_path / "predictd_model.R").read_text()
        # the model window has 1 + 4 * bw + 10 points and the lags
        # cover 4 * bw points
        assert len(r_vector(rscript, "p")) == 1 + 4 * bw + 10
        assert len(r_vector(rscript, "ycorr")) == 4 * bw
    else:
        assert msgs[-4:] == not_enough_pairs(0)


@pytest.mark.parametrize("d_min", [0, 10, 150])
def test_d_min_below_d_keeps_d(macs3_argparser, caplog, model_bed, tmp_path,
                               d_min):
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                   "--d-min", d_min, "--outdir", tmp_path])
    assert abs(predicted_d(cmd_messages(caplog)) - 200) <= 1


def test_d_min_above_d_excludes_it(macs3_argparser, caplog, model_bed,
                                   tmp_path):
    """--d-min 300 excludes the lag of 200; whatever is reported must be
    above the minimum."""
    try:
        run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                       "--d-min", "300",
                                       "--outdir", tmp_path])
    except AssertionError as exc:
        assert "No proper d can be found! Tweak --mfold?" in str(exc)
        return
    d = predicted_d(cmd_messages(caplog))
    assert d > 300


def test_tsize_option(macs3_argparser, caplog, model_bed, tmp_path):
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE, "-s", "50",
                                   "--outdir", tmp_path])
    msgs = [m for _, m in cmd_messages(caplog)]
    assert "tag size is determined as 50 bps" in msgs
    assert "# tag size = 50" in msgs
    # the tag size does not enter the model
    assert abs(predicted_d(cmd_messages(caplog)) - 200) <= 1


# ------------------------------------
# --rfile and --outdir: the R script
# ------------------------------------

def r_vector(text, name):
    for line in text.splitlines():
        if line.startswith(name + " ") and "<- c(" in line:
            body = line.split("<- c(", 1)[1].rstrip(")")
            return [float(x) for x in body.split(",") if x.strip()]
    raise AssertionError("no vector %s in R script" % name)


def expected_strand_lines(frag, n_sites=150, ntags=10, bw=300):
    """Tag profiles around the paired centres, as projected by the model.

    Strand peak summits are base + (ntags-1)//2 and base + frag +
    (ntags-1)//2; the centre c is their integer mean.  A tag at t covers
    window indices [t - c + 2*bw, t - c + 2*bw + 10).
    """
    peaksize = 2 * bw
    w = 1 + 2 * peaksize + 10
    plus = np.zeros(w)
    minus = np.zeros(w)
    half = (ntags - 1) // 2
    c = (half + frag + half) // 2          # centre relative to base
    for j in range(ntags):
        s = j - c + peaksize
        plus[s:s + 10] += n_sites
        s = frag + j - c + peaksize
        minus[s:s + 10] += n_sites
    return plus, minus


def test_rfile_content(macs3_argparser, model_bed, tmp_path):
    out = tmp_path / "models"
    out.mkdir()
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                   "--rfile", "mymodel.R",
                                   "--outdir", out])
    text = (out / "mymodel.R").read_text()
    lines = text.splitlines()
    assert lines[:2] == ["# R script for Peak Model",
                         "#  -- generated by MACS"]
    plus, minus = expected_strand_lines(200)
    assert r_vector(text, "p") == pytest.approx(plus * 100 / plus.sum(),
                                                rel=1e-9)
    assert r_vector(text, "m") == pytest.approx(minus * 100 / minus.sum(),
                                                rel=1e-9)
    assert len(r_vector(text, "ycorr")) == 1200
    assert len(r_vector(text, "xcorr")) == 1200
    altd = r_vector(text, "altd")
    assert len(altd) == 1 and abs(altd[0] - 200) <= 1
    assert "pdf('mymodel.R_model.pdf',height=6,width=6)" in lines
    assert "legend('right','alt lag(s) : %d',bty='n')" % altd[0] in lines
    assert lines[-1] == "dev.off()"


def test_rfile_cross_correlation_peaks_at_true_lag(macs3_argparser,
                                                   model_bed, tmp_path):
    """The smoothed cross-correlation in the R script peaks at the index
    whose lag is exactly 200 (index = lag + peaksize - 1)."""
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                   "--outdir", tmp_path])
    ycorr = np.array(r_vector((tmp_path / "predictd_model.R").read_text(),
                              "ycorr"))
    assert int(np.argmax(ycorr)) - 600 + 1 == 200


def test_cli_rfile_in_outdir(run_macs3, model_bed, tmp_path):
    proc = run_macs3(["predictd", "-i", model_bed, "-g", GSIZE,
                      "--outdir", tmp_path / "new"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert [f.name for f in (tmp_path / "new").iterdir()] == \
        ["predictd_model.R"]


# ------------------------------------
# -f/--format
# ------------------------------------

@pytest.mark.parametrize("fmt", SE_FORMATS)
def test_formats_explicit(macs3_argparser, caplog, write_bed,
                          make_alignments, tmp_path, fmt):
    reads = model_reads(n_sites=100)
    path = write_se(fmt, reads, write_bed, make_alignments)
    run_predictd(macs3_argparser, ["-i", path, "-f", fmt, "-g", GSIZE,
                                   "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    assert ("INFO", "# total tags in alignment file: 2000") in msgs
    assert ("INFO", "#2 Total number of paired peaks: 100") in msgs
    assert abs(predicted_d(msgs) - 200) <= 1


@pytest.mark.parametrize("fmt", AUTO_FORMATS)
def test_formats_auto(macs3_argparser, caplog, write_bed, make_alignments,
                      tmp_path, fmt):
    path = write_se(fmt, model_reads(n_sites=100), write_bed,
                    make_alignments)
    run_predictd(macs3_argparser, ["-i", path, "-g", GSIZE,
                                   "--outdir", tmp_path])
    assert abs(predicted_d(cmd_messages(caplog)) - 200) <= 1


def test_multiple_input_files(macs3_argparser, caplog, write_bed, tmp_path):
    reads = model_reads()
    a = write_bed(bed_lines(reads[:1000]), name="a.bed")
    b = write_bed(bed_lines(reads[1000:]), name="b.bed")
    run_predictd(macs3_argparser, ["-i", a, b, "-g", GSIZE,
                                   "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    assert ("INFO", "# total tags in alignment file: 3000") in msgs
    assert abs(predicted_d(msgs) - 200) <= 1


@pytest.mark.parametrize("buffer_size", [1, 1000])
def test_buffer_size(macs3_argparser, caplog, model_bed, tmp_path,
                     buffer_size):
    run_predictd(macs3_argparser, ["-i", model_bed, "-g", GSIZE,
                                   "--buffer-size", buffer_size,
                                   "--outdir", tmp_path])
    assert abs(predicted_d(cmd_messages(caplog)) - 200) <= 1


# ------------------------------------
# paired-end formats: the average fragment length
# ------------------------------------

def write_bampe(frags, make_alignments, make_pe_pair):
    recs = []
    for i, (c, s, e) in enumerate(frags):
        recs += make_pe_pair("p%d" % i, c, s, e)
    return make_alignments(recs, refs=(("chr1", 10000), ("chr2", 10000)),
                           name="frags.bam")


@pytest.mark.parametrize("fmt", ["BEDPE", "BAMPE"])
def test_paired_end_average_length(macs3_argparser, caplog, write_bedpe,
                                   make_alignments, make_pe_pair, tmp_path,
                                   fmt):
    """Lengths 200, 150 and 251 average 200.33, reported as 200; the
    model options and --rfile have no effect."""
    if fmt == "BEDPE":
        path = write_bedpe(PE_FRAGS)
    else:
        path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair)
    run_predictd(macs3_argparser, ["-i", path, "-f", fmt, "--bw", "50",
                                   "-m", "10", "50", "--d-min", "500",
                                   "--rfile", "x.R", "--outdir", tmp_path])
    assert [m for _, m in cmd_messages(caplog)] == [
        "# read input file in Paired-end mode.",
        "# read treatment fragments...",
        "# total fragments/pairs in alignment file: 3",
        "# Build Peak Model...",
        "# Average insertion length of all pairs is 200 bps"]
    assert not (tmp_path / "x.R").exists()


def test_cli_paired_end_log(run_macs3, parse_log, make_alignments,
                            make_pe_pair, tmp_path):
    path = write_bampe(PE_FRAGS, make_alignments, make_pe_pair)
    proc = run_macs3(["predictd", "-i", path, "-f", "BAMPE",
                      "--outdir", tmp_path], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert parse_log(proc.stderr) == [
        ("INFO", "# read input file in Paired-end mode."),
        ("INFO", "# read treatment fragments..."),
        ("INFO", "3 fragments have been read."),
        ("INFO", "# total fragments/pairs in alignment file: 3"),
        ("INFO", "# Build Peak Model..."),
        ("INFO", "# Average insertion length of all pairs is 200 bps")]


@pytest.mark.parametrize("frags, expect", [
    ([("chr1", 0, 100)], 100),
    ([("chr1", 0, 100), ("chr1", 10, 111)], 100),     # 100.5 -> 100
    ([("chr1", 0, 99), ("chr1", 10, 112)], 100),      # 100.5 -> 100
    ([("chr1", 5, 6)] * 3, 1),
])
def test_paired_end_average_truncated(macs3_argparser, caplog, write_bedpe,
                                      tmp_path, frags, expect):
    path = write_bedpe(frags)
    run_predictd(macs3_argparser, ["-i", path, "-f", "BEDPE",
                                   "--outdir", tmp_path])
    assert cmd_messages(caplog)[-1] == \
        ("INFO", "# Average insertion length of all pairs is %d bps" % expect)


def test_paired_end_multiple_files(macs3_argparser, caplog, write_bedpe,
                                   tmp_path):
    a = write_bedpe(PE_FRAGS[:1], name="a.bedpe")
    b = write_bedpe(PE_FRAGS[1:], name="b.bedpe")
    run_predictd(macs3_argparser, ["-i", a, b, "-f", "BEDPE",
                                   "--outdir", tmp_path])
    assert cmd_messages(caplog)[-1] == \
        ("INFO", "# Average insertion length of all pairs is 200 bps")


# ------------------------------------
# upstream test data (test/cmdlinetest)
# ------------------------------------


@pytest.mark.parametrize("fmt, infile, std", [
    ("BAMPE", "CTCF_PE_ChIP_chr22_50k.bam", "run_predictd_bampe.txt"),
    ("BEDPE", "CTCF_PE_ChIP_chr22_50k.bedpe.gz", "run_predictd_bedpe.txt"),
])
def test_standard_result_pe(macs3_argparser, caplog, test_dir, tmp_path, fmt,
                            infile, std):
    """Pins the current output. The average fragment length comes from
    upstream's test/standard_results_predictd."""
    expect = (test_dir / "standard_results_predictd" / std).read_text().split()
    run_predictd(macs3_argparser, ["-i", test_dir / infile, "-f", fmt,
                                   "--outdir", tmp_path])
    assert cmd_messages(caplog)[-1] == \
        ("INFO", "# Average insertion length of all pairs is %s bps" %
         expect[0])


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


def test_load_tag_files_options(write_bed):
    reads = model_reads(n_sites=3)
    a = write_bed(bed_lines(reads[:20]), name="a.bed")
    b = write_bed(bed_lines(reads[20:]), name="b.bed")
    opts = _Opts(BEDParser, [a, b])
    track = load_tag_files_options(opts)
    assert track.total == 60
    assert opts.tsize == 36
    assert opts.messages == ["# read treatment tags...",
                             "tag size is determined as 36 bps"]


def test_load_tag_files_options_given_tsize(write_bed):
    opts = _Opts(BEDParser, [write_bed(bed_lines(model_reads(n_sites=1)))],
                 tsize=25)
    load_tag_files_options(opts)
    assert opts.tsize == 25
    assert opts.messages[-1] == "tag size is determined as 25 bps"


def test_load_frag_files_options(write_bedpe):
    a = write_bedpe(PE_FRAGS[:2], name="a.bedpe")
    b = write_bedpe(PE_FRAGS[2:], name="b.bedpe")
    opts = _Opts(BEDPEParser, [a, b])
    track = load_frag_files_options(opts)
    assert opts.messages == ["# read treatment fragments..."]
    assert track.total == 3
    assert track.average_template_length == pytest.approx(601 / 3, rel=1e-6)


def test_run_returns_none_and_sets_d(macs3_argparser, model_bed, tmp_path):
    options = macs3_argparser.parse_args(["predictd", "-i", str(model_bed),
                                          "-g", GSIZE,
                                          "--outdir", str(tmp_path)])
    assert run(options) is None
    assert options.PE_MODE is False
    assert abs(options.d - 200) <= 1
    assert options.modelR == str(tmp_path / "predictd_model.R")
    assert (options.lmfold, options.umfold) == (5, 50)


# ------------------------------------
# opt_validate_predictd checks that argparse makes unreachable
# ------------------------------------

def test_validator_rejects_unknown_format(macs3_argparser, caplog, model_bed,
                                          tmp_path):
    """-f has fixed choices, so only a direct run() call reaches this."""
    options = macs3_argparser.parse_args(["predictd", "-i", str(model_bed),
                                          "--outdir", str(tmp_path)])
    options.format = "frag"
    with pytest.raises(SystemExit) as exc:
        run(options)
    assert exc.value.code == 1
    assert cmd_messages(caplog) == [
        ("ERROR", "Format \"FRAG\" cannot be recognized!")]


@pytest.mark.parametrize("fmt, nomodel", [("bedpe", True), ("bed", None)])
def test_validator_uppercases_format_and_sets_nomodel(macs3_argparser,
                                                      write_bedpe, model_bed,
                                                      tmp_path, fmt, nomodel):
    path = write_bedpe(PE_FRAGS) if fmt == "bedpe" else model_bed
    options = macs3_argparser.parse_args(["predictd", "-i", str(path),
                                          "-g", GSIZE,
                                          "--outdir", str(tmp_path)])
    options.format = fmt
    run(options)
    assert options.format == fmt.upper()
    assert getattr(options, "nomodel", None) is nomodel


# ------------------------------------
# paired peaks from several chromosomes
# ------------------------------------

def test_paired_peaks_pooled_over_chromosomes(macs3_argparser, caplog,
                                              write_bed, tmp_path):
    """75 sites on each of two chromosomes: neither alone reaches 100
    pairs, together they give 150 and the model is built."""
    reads = model_reads(n_sites=75, chrom="chr1") + \
        model_reads(n_sites=75, chrom="chr2")
    path = write_bed(bed_lines(reads), name="two.bed")
    run_predictd(macs3_argparser, ["-i", path, "-g", GSIZE, "--verbose", "3",
                                   "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    assert [m for _, m in msgs].count(
        "Number of paired peaks in this chromosome: 75") == 2
    assert ("INFO", "#2 Total number of paired peaks: 150") in msgs
    assert abs(predicted_d(msgs) - 200) <= 1


def test_empty_input_finds_no_pairs(macs3_argparser, caplog, tmp_path):
    path = tmp_path / "empty.bed"
    path.write_text("")
    run_predictd(macs3_argparser, ["-i", path, "-f", "BED", "-g", GSIZE,
                                   "--outdir", tmp_path])
    msgs = cmd_messages(caplog)
    assert ("INFO", "# total tags in alignment file: 0") in msgs
    assert msgs[-4:] == not_enough_pairs(0)
    assert not (tmp_path / "predictd_model.R").exists()
