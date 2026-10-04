#!/usr/bin/env python

"""Module Description: Test the callpeak command (callpeak_cmd.py).

Most tests run ``macs3 callpeak`` in a subprocess (``run_macs3``) on
small synthetic inputs whose pileups, local lambdas, p-scores and peak
boundaries can be derived by hand, so the expected values are written
next to each test. Realistic runs use the CTCF chr22 subsets in
``test/``; their outputs are compared with upstream's
``test/standard_results_callpeak_*`` files (these pin the output),
except for the runs whose stored results carry a bug marked below.

Conventions used for the hand derivations below:

* single-end reads are extended from their 5' end towards 3' by
  ``d`` (``--extsize`` with ``--nomodel``): a ``+`` read starting at
  ``s`` covers ``[s, s + d)``; a ``-`` read ending at ``e`` covers
  ``[e - d, e)``. ``--shift S`` moves the 5' end by ``S`` towards 3'.
* without a control the local lambda is the treatment itself piled up
  in a window of ``llocal`` bp centred on each 5' end (clipped at 0),
  scaled by ``d / llocal``, and floored at the genome background
  ``lambda_bg = d * total / gsize``.
* with a control the lambda is the maximum over windows of ``d``,
  ``slocal`` and ``llocal`` bp centred on the control 5' ends, scaled by
  ``d / window`` (times the treatment/control depth ratio when the
  control is scaled to the treatment), floored at ``lambda_bg``.
* ``-log10 p = -log10 P(X > k; lambda)`` (Poisson upper tail without
  ``k``), stored as float32; q-scores follow MACS's ranking of p-scores
  weighted by the length they cover (``_ref_qscores``).
* the fold enrichment is ``(pileup + 1) / (lambda + 1)``.
* the output bedGraph and the peak scan stop at the last breakpoint of
  the shorter of the treatment and control tracks of a chromosome (the
  synthetic controls below reach past the treatment; a control that
  ends first truncates the treatment, which is marked as a bug).

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import filecmp
import logging
import math
import os
import re
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import binom, poisson

from MACS3.Commands.callpeak_cmd import (check_names,
                                         run,
                                         cal_max_dup_tags,
                                         load_frag_files_options,
                                         load_tag_files_options)
from MACS3.Utilities.OptValidator import (opt_validate_callpeak)
from MACS3.Utilities.Constants import (EFFECTIVEGS,
                                       MACS_VERSION)

# ------------------------------------
# shared constants and helpers
# ------------------------------------

@pytest.fixture(autouse=True)
def _restore_optvalidator_level():
    """opt_validate_callpeak sets the level of its module logger from
    --verbose; restore it after in-process calls."""
    lg = logging.getLogger("MACS3.Utilities.OptValidator")
    level = lg.level
    yield
    lg.setLevel(level)


# synthetic runs: no model, d = 100, small genome, keep all duplicates
SYN = ["--nomodel", "--extsize", "100", "-g", "100000", "--keep-dup", "all"]

NARROW_SUFFIXES = ("_peaks.xls", "_peaks.narrowPeak", "_summits.bed")
BDG_SUFFIXES = ("_treat_pileup.bdg", "_control_lambda.bdg")
COMPARE_SUFFIXES = ("_peaks.narrowPeak", "_summits.bed",
                    "_treat_pileup.bdg", "_control_lambda.bdg")

USAGE = "usage: macs3 callpeak"


def _block(n, start=1000, chrom="chr1", strand="+", readlen=50):
    """``n`` identical BED reads ``[start, start + readlen)``."""
    return [(chrom, start, start + readlen, "r%d" % i, 0, strand)
            for i in range(n)]


def _cp(run_macs3, args, outdir, name="t", timeout=60):
    """Run ``macs3 callpeak <args> -n name --outdir outdir``."""
    return run_macs3(["callpeak"] + [str(a) for a in args]
                     + ["-n", name, "--outdir", str(outdir)],
                     timeout=timeout)


def _ok(proc):
    assert proc.returncode == 0, proc.stderr[-3000:]
    return proc


def _out(outdir, name, suffix):
    return Path(outdir) / (name + suffix)


def _lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def _rows(path):
    """Non-comment, non-track, non-empty lines split on tabs."""
    return [ln.split("\t") for ln in _lines(path)
            if ln and not ln.startswith("#") and not ln.startswith("track")]


def _bdg(path):
    return [(r[0], int(r[1]), int(r[2]), r[3]) for r in _rows(path)]


def _bdg_value(rows, chrom, x):
    for c, s, e, v in rows:
        if c == chrom and s <= x < e:
            return float(v)
    raise KeyError((chrom, x))


def _xls_table(path):
    """Header row and data rows of a peaks.xls file."""
    rows = _rows(path)
    return rows[0], rows[1:]


def _f32(x):
    return float(np.float32(x))


def _pscore(k, lam):
    """-log10 P(X > k) for X ~ Poisson(lam), as MACS stores it (float32)."""
    return _f32(-poisson.logsf(k, _f32(lam)) / math.log(10))


def _ref_qscores(segments):
    """q-scores for p-scores covering given lengths.

    ``segments`` is an iterable of ``(pscore, length)``. P-scores are
    ranked from high to low; a p-score whose higher-scoring regions
    cover ``k - 1`` bp gets ``q = p + log10(k) - log10(N)`` with ``N``
    the total length, made monotonic, and set to 0 from the first
    non-positive value on (including the lowest p-score).
    """
    stat = {}
    for p, ln in segments:
        stat[p] = stat.get(p, 0) + ln
    n_total = sum(stat.values())
    f = _f32(-math.log10(n_total))
    values = sorted(stat, reverse=True)
    out = {}
    k = 1
    pre_q = 2147483647.0
    i = 0
    for i, v in enumerate(values):
        q = _f32(v + (math.log10(k) + f))
        if q > pre_q:
            q = pre_q
        if q <= 0:
            break
        out[v] = q
        pre_q = q
        k += stat[v]
    for v in values[i:]:
        out[v] = 0.0
    return out


def _assert_fields(row, expected, rel=1e-5):
    """Compare a split row with expected values: floats approximately
    (relative tolerance ``rel``), everything else as exact strings."""
    assert len(row) == len(expected), row
    for got, exp in zip(row, expected):
        if isinstance(exp, float):
            assert float(got) == pytest.approx(exp, rel=rel, abs=1e-6), row
        else:
            assert got == str(exp), row


def _coverage_segments(intervals, scale=1.0):
    """[(start, end, value)] of the weighted coverage of ``intervals``
    (``(s, e, w)``) from 0 to the last end, merging equal neighbours."""
    pts = sorted({0} | {s for s, _, _ in intervals}
                 | {e for _, e, _ in intervals})
    segs = []
    for a, b in zip(pts[:-1], pts[1:]):
        v = scale * sum(w for s, e, w in intervals if s <= a < e)
        if segs and segs[-1][2] == v:
            segs[-1] = (segs[-1][0], b, v)
        else:
            segs.append((a, b, v))
    return segs


def _bdg_lines(chrom, segments):
    return ["%s\t%d\t%d\t%.5f" % (chrom, s, e, v) for s, e, v in segments]


def _same_outputs(dir1, dir2, name, suffixes=COMPARE_SUFFIXES):
    for suffix in suffixes:
        f1 = _out(dir1, name, suffix)
        f2 = _out(dir2, name, suffix)
        assert f1.exists() == f2.exists(), suffix
        if f1.exists():
            assert filecmp.cmp(f1, f2, shallow=False), suffix


def _xls_without_cmdline(path):
    return [ln for ln in _lines(path) if not ln.startswith("# Command line:")]


def _messages(parse_log, proc, level=None):
    msgs = parse_log(proc.stderr)
    if level is None:
        return [m for _, m in msgs]
    return [m for lv, m in msgs if lv == level]


# reads used for the format-equivalence tests: 12 + and 12 - reads of
# 36 bp; (start0, strand)
SE_READS = ([(1000 + 3 * i, "+") for i in range(12)]
            + [(1100 + 3 * i, "-") for i in range(12)])
SE_READLEN = 36
SE_SEQ = "A" * SE_READLEN
SE_QUAL = "I" * SE_READLEN

# fragments used for the paired-end format tests: (left, right)
PE_FRAGS = [(5000 + 7 * i, 5000 + 7 * i + 150 + i) for i in range(15)]


def _se_expected_intervals(reads, d=100):
    out = []
    for s, strand in reads:
        if strand == "+":
            out.append((s, s + d, 1))
        else:
            e = s + SE_READLEN
            out.append((e - d, e, 1))
    return out


def _write_se_format(tmp_path, fmt, reads, make_alignments):
    """Write ``reads`` in one of the single-end formats; return the path."""
    L = SE_READLEN
    if fmt == "BED":
        path = tmp_path / "reads.bed"
        with open(path, "w") as fh:
            for i, (s, st) in enumerate(reads):
                fh.write("chr1\t%d\t%d\tr%d\t0\t%s\n" % (s, s + L, i, st))
    elif fmt == "BAM" or fmt == "SAM":
        recs = [dict(name="r%d" % i, ref="chr1", pos=s,
                     flag=0 if st == "+" else 16, cigar="%dM" % L)
                for i, (s, st) in enumerate(reads)]
        return make_alignments(recs, name="reads." + fmt.lower(),
                               fmt=fmt.lower())
    elif fmt == "ELAND":
        # name, sequence, match code, 3 counts, chrom(.fa), 1-based pos, F/R
        path = tmp_path / "reads.eland"
        with open(path, "w") as fh:
            for i, (s, st) in enumerate(reads):
                fh.write(">r%d\t%s\tU0\t1\t0\t0\tchr1.fa\t%d\t%s\t..\n"
                         % (i, SE_SEQ, s + 1, "F" if st == "+" else "R"))
    elif fmt == "ELANDMULTI":
        # name, sequence, hit counts, chrom.fa:<1-based pos><F/R><mismatches>
        path = tmp_path / "reads.elandmulti"
        with open(path, "w") as fh:
            for i, (s, st) in enumerate(reads):
                fh.write(">r%d\t%s\t1:0:0\tchr1.fa:%d%s0\n"
                         % (i, SE_SEQ, s + 1, "F" if st == "+" else "R"))
    elif fmt == "ELANDEXPORT":
        # 22 columns: [8] sequence, [10] chrom, [12] 1-based pos, [13] F/R
        path = tmp_path / "reads.elandexport"
        with open(path, "w") as fh:
            for i, (s, st) in enumerate(reads):
                fields = ["M1", "1", "1", "1", str(i), str(i), "0", "1",
                          SE_SEQ, SE_QUAL, "chr1", "", str(s + 1),
                          "F" if st == "+" else "R", str(L), "100",
                          "", "", "", "", "", "Y"]
                fh.write("\t".join(fields) + "\n")
    elif fmt == "BOWTIE":
        # name, strand, reference, 0-based offset, sequence, quals, count, mm
        path = tmp_path / "reads.bowtie"
        with open(path, "w") as fh:
            for i, (s, st) in enumerate(reads):
                fh.write("r%d\t%s\tchr1\t%d\t%s\t%s\t0\t\n"
                         % (i, st, s, SE_SEQ, SE_QUAL))
    else:
        raise ValueError(fmt)
    return str(path)


def _write_pe_format(tmp_path, fmt, frags, make_alignments, make_pe_pair):
    if fmt == "BAMPE":
        recs = []
        for i, (l, r) in enumerate(frags):
            recs += make_pe_pair("f%d" % i, "chr1", l, r, readlen=36)
        return make_alignments(recs, name="frags.bam")
    path = tmp_path / ("frags." + fmt.lower())
    with open(path, "w") as fh:
        for l, r in frags:
            if fmt == "BEDPE":
                fh.write("chr1\t%d\t%d\n" % (l, r))
            else:
                fh.write("chr1\t%d\t%d\tBC1\t1\n" % (l, r))
    return str(path)


# ------------------------------------
# realistic runs compared with upstream's standard results
# ------------------------------------

# (name, args, standard results folder) exactly as in test/cmdlinetest
STANDARD_RUNS = [
    ("run_callpeak_narrow0",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B",
      "--cutoff-analysis"], "callpeak_narrow"),
    ("run_callpeak_narrow1",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B", "--d-min", "15",
      "--call-summits"], "callpeak_narrow"),
    ("run_callpeak_narrow2",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B", "--nomodel",
      "--extsize", "100"], "callpeak_narrow"),
    ("run_callpeak_narrow3",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B", "--nomodel",
      "--extsize", "100", "--shift", "-50"], "callpeak_narrow"),
    ("run_callpeak_narrow4",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B", "--nomodel",
      "--nolambda", "--extsize", "100", "--shift", "-50"],
     "callpeak_narrow"),
    ("run_callpeak_narrow5",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B", "--scale-to",
      "large"], "callpeak_narrow"),
    ("run_callpeak_broad",
     ["-g", "52000000", "-t", "CHIP", "-c", "CTRL", "-B", "--broad"],
     "callpeak_broad"),
    ("run_callpeak_bampe_narrow",
     ["-g", "52000000", "-f", "BAMPE", "-t", "CHIPPE", "-c", "CTRLPE", "-B",
      "--call-summits"], "callpeak_pe_narrow"),
    ("run_callpeak_bedpe_narrow",
     ["-g", "52000000", "-f", "BEDPE", "-t", "CHIPBEDPE", "-c", "CTRLBEDPE",
      "-B", "--call-summits"], "callpeak_pe_narrow"),
    ("run_callpeak_pe_narrow_onlychip",
     ["-g", "52000000", "-f", "BEDPE", "-t", "CHIPBEDPE", "-B"],
     "callpeak_pe_narrow"),
    ("run_callpeak_bampe_broad",
     ["-g", "52000000", "-f", "BAMPE", "-t", "CHIPPE", "-c", "CTRLPE", "-B",
      "--broad"], "callpeak_pe_broad"),
    ("run_callpeak_bedpe_broad",
     ["-g", "52000000", "-f", "BEDPE", "-t", "CHIPBEDPE", "-c", "CTRLBEDPE",
      "-B", "--broad"], "callpeak_pe_broad"),
    ("run_callpeak_frag",
     ["-f", "FRAG", "-t", "FRAGFILE", "-B"], "callpeak_frag"),
    ("run_callpeak_frag_barcode",
     ["-f", "FRAG", "-t", "FRAGFILE", "--barcodes", "BARCODESFILE", "-B"],
     "callpeak_frag"),
    ("run_callpeak_narrow_revert",
     ["-g", "10000000", "--nomodel", "--extsize", "250", "-c", "CHIP",
      "-t", "CTRL", "-B"], "callpeak_narrow_revert"),
]

DATA_FILES = {"CHIP": "CTCF_SE_ChIP_chr22_50k.bed.gz",
              "CTRL": "CTCF_SE_CTRL_chr22_50k.bed.gz",
              "CHIPPE": "CTCF_PE_ChIP_chr22_50k.bam",
              "CTRLPE": "CTCF_PE_CTRL_chr22_50k.bam",
              "CHIPBEDPE": "CTCF_PE_ChIP_chr22_50k.bedpe.gz",
              "CTRLBEDPE": "CTCF_PE_CTRL_chr22_50k.bedpe.gz",
              "FRAGFILE": "test.fragments.tsv.gz",
              "BARCODESFILE": "barcodes.txt"}


def _ctcf_args(test_dir, args):
    return [str(test_dir / DATA_FILES[a]) if a in DATA_FILES else a
            for a in args]


# Runs whose standard results carry a bug marked below, so they are not
# pinned: the model-based runs use the 1 bp short d of PeakModel, and in
# the reverted run the control (the ChIP) ends before the treatment.
MODEL_RUN_NAMES = ("run_callpeak_narrow0", "run_callpeak_narrow1",
                   "run_callpeak_narrow5", "run_callpeak_broad")
REVERT_RUN = [r for r in STANDARD_RUNS
              if r[0] == "run_callpeak_narrow_revert"][0]
PINNED_RUNS = [r for r in STANDARD_RUNS
               if r[0] not in MODEL_RUN_NAMES and r is not REVERT_RUN]
OUTPUT_KINDS = ("peaks", "summits", "treat_pileup", "control_lambda")


def _output_files(folder, name):
    """Peak, summit and bedGraph files of run ``name`` in ``folder``."""
    return sorted(p for p in Path(folder).glob(name + "_*")
                  if p.name[len(name) + 1:].split(".")[0] in OUTPUT_KINDS
                  and not p.name.endswith(".xls"))


@pytest.mark.parametrize("name,args,folder", PINNED_RUNS,
                         ids=[r[0] for r in PINNED_RUNS])
def test_callpeak_matches_standard_results(run_macs3, tmp_path, test_dir,
                                           name, args, folder):
    """Every file upstream keeps for this cmdlinetest run is reproduced
    byte for byte.

    Pins the current output. The expected files are upstream's
    standard results for test/cmdlinetest; a real-data peak call has no
    practical independent derivation. The model-based runs and the
    reverted run are checked by the tests below instead.
    """
    out = tmp_path / "out"
    _ok(_cp(run_macs3, _ctcf_args(test_dir, args), out, name=name))
    expected = _output_files(test_dir / ("standard_results_" + folder), name)
    assert expected
    for std in expected:
        got = out / std.name
        assert got.exists(), std.name
        assert got.read_bytes() == std.read_bytes(), std.name


def _revert_run(run_macs3, tmp_path, test_dir):
    name, args, _ = REVERT_RUN
    out = tmp_path / "out"
    _ok(_cp(run_macs3, _ctcf_args(test_dir, args), out, name=name))
    return name, out


def test_callpeak_revert_bedgraphs_match_standard(run_macs3, tmp_path,
                                                  test_dir):
    """The reverted cmdlinetest run (control reads as treatment, ChIP
    reads as control) reproduces upstream's bedGraphs up to their last
    line.

    Pins the current output. The last line of each stored bedGraph ends
    where the shorter control ends (the treatment is truncated to the
    control's length), and the peak
    files depend on the truncated length through the q-values, so only
    the lines before it are compared.
    """
    name, out = _revert_run(run_macs3, tmp_path, test_dir)
    folder = test_dir / ("standard_results_" + REVERT_RUN[2])
    for suffix in BDG_SUFFIXES:
        std = _lines(folder / (name + suffix))
        assert len(std) > 1000
        assert _lines(_out(out, name, suffix))[:len(std) - 1] == std[:-1]


# ------------------------------------
# which files are written
# ------------------------------------

@pytest.mark.parametrize("extra,suffixes", [
    ([], NARROW_SUFFIXES),
    (["-B"], NARROW_SUFFIXES + BDG_SUFFIXES),
    (["--broad"], ("_peaks.xls", "_peaks.broadPeak", "_peaks.gappedPeak")),
    (["--broad", "-B"], ("_peaks.xls", "_peaks.broadPeak",
                         "_peaks.gappedPeak") + BDG_SUFFIXES),
    (["--cutoff-analysis"], NARROW_SUFFIXES + ("_cutoff_analysis.txt",)),
    (["-B", "--SPMR", "--cutoff-analysis", "--call-summits"],
     NARROW_SUFFIXES + BDG_SUFFIXES + ("_cutoff_analysis.txt",)),
    (["--SPMR"], NARROW_SUFFIXES),
], ids=["default", "bdg", "broad", "broad-bdg", "cutoff-analysis",
        "bdg-spmr-cutoff-summits", "spmr-without-bdg"])
def test_output_files_created(run_macs3, tmp_path, write_bed, extra,
                              suffixes):
    # --nomodel never writes _model.r; only the listed files appear
    bed = write_bed(_block(20))
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", bed] + SYN + extra, out))
    assert sorted(p.name for p in out.iterdir()) == \
        sorted("t" + s for s in suffixes)


def test_output_files_with_model(run_macs3, tmp_path, test_dir):
    # building the model adds NAME_model.r
    out = tmp_path / "out"
    _ok(_cp(run_macs3, _ctcf_args(test_dir, ["-g", "52000000", "-t", "CHIP",
                                             "-c", "CTRL"]), out))
    assert sorted(p.name for p in out.iterdir()) == \
        sorted("t" + s for s in NARROW_SUFFIXES + ("_model.r",))


def test_tmpdir_is_empty_after_success(run_macs3, tmp_path, write_bed,
                                       macs3_tmpdir):
    # the per-chromosome pileup files written to TMPDIR are removed
    bed = write_bed(_block(20) + _block(5, chrom="chr2"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN, tmp_path / "out"))
    assert list(macs3_tmpdir.iterdir()) == []


# ------------------------------------
# a hand-checkable block: 20 identical + reads at chr1:1000, no control
# ------------------------------------
# d = 100, so the treatment pileup is 20 on [1000, 1100).
# lambda_bg = 100 * 20 / 1e5 = 0.02; the llocal window (10000 bp) around
# position 1000 is [0, 6000) with 20 * 100 / 10000 = 0.2 > lambda_bg.
# p-score at the block: P20 = -log10 P(X > 20; 0.2) = 34.4696;
# outside the block (pileup 0, lambda 0.2) P0 = 0.741669.
# q at the block = P20 + log10(1) - log10(1100) = 31.4282.
# fold enrichment = (20 + 1) / (0.2 + 1) = 17.5; summit at 1050.

P20 = _pscore(20, 0.2)
Q20 = _ref_qscores([(P20, 100), (_pscore(0, 0.2), 1000)])[P20]


def test_block_treat_pileup(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN + ["-B"], out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000",
        "chr1\t1000\t1100\t20.00000"]


def test_block_control_lambda(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN + ["-B"], out))
    # truncated at the treatment's last breakpoint (1100)
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t1100\t0.20000"]


def test_block_narrowpeak(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN, out))
    rows = _rows(_out(out, "t", "_peaks.narrowPeak"))
    assert len(rows) == 1
    # score column = int(10 * q); summit column = offset from start
    _assert_fields(rows[0], ["chr1", "1000", "1100", "t_peak_1",
                             str(int(10 * Q20)), ".", 17.5, P20, Q20, "50"])


def test_block_summits(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN, out))
    rows = _rows(_out(out, "t", "_summits.bed"))
    assert len(rows) == 1
    _assert_fields(rows[0], ["chr1", "1050", "1051", "t_peak_1", Q20])


def test_block_xls(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    bed = write_bed(_block(20))
    args = ["-t", bed] + SYN
    _ok(_cp(run_macs3, args, out))
    lines = _lines(_out(out, "t", "_peaks.xls"))
    cmd = " ".join(["callpeak"] + args + ["-n", "t", "--outdir", str(out)])
    assert lines[:-2] == [
        "# This file is generated by MACS version %s" % MACS_VERSION,
        "# Command line: %s" % cmd,
        "# ARGUMENTS LIST:",
        "# name = t",
        "# format = AUTO",
        "# ChIP-seq file = ['%s']" % bed,
        "# control file = None",
        "# effective genome size = 1.00e+05",
        "# band width = 300",
        "# model fold = [5, 50]",
        "# qvalue cutoff = 5.00e-02",
        "# The maximum gap between significant sites is assigned as the "
        "read length/tag size.",
        "# The minimum length of peaks is assigned as the predicted "
        "fragment length \"d\".",
        "# Larger dataset will be scaled towards smaller dataset.",
        "# Range for calculating regional lambda is: 10000 bps",
        "# Broad region calling is off",
        "# Paired-End mode is off",
        "",
        "# tag size is determined as 50 bps",
        "# total tags in treatment: 20",
        "# d = 100"]
    assert lines[-2].split("\t") == [
        "chr", "start", "end", "length", "abs_summit", "pileup",
        "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"]
    # start is 1-based in the xls; abs_summit is 1-based too
    _assert_fields(lines[-1].split("\t"),
                   ["chr1", "1001", "1100", "100", "1051", "20", P20, 17.5,
                    Q20, "t_peak_1"])


def test_block_log_messages(run_macs3, tmp_path, write_bed, parse_log):
    out = tmp_path / "out"
    bed = write_bed(_block(20))
    proc = _ok(_cp(run_macs3, ["-t", bed, "-f", "BED"] + SYN + ["-B"], out))
    msgs = [m for lv, m in parse_log(proc.stderr) if lv]
    assert msgs == [
        "",
        "#1 read tag files...",
        "#1 read treatment tags...",
        "#1 tag size is determined as 50 bps",
        "#1 tag size = 50.0",
        "#1  total tags in treatment: 20",
        "#1 finished!",
        "#2 Build Peak Model...",
        "#2 Skipped...",
        "#2 Use 100 as fragment length",
        "#3 Call peaks...",
        "#3 Pre-compute pvalue-qvalue table...",
        "#3 In the peak calling step, the following will be performed "
        "simultaneously:",
        "#3   Write bedGraph files for treatment pileup (after scaling if "
        "necessary)... t_treat_pileup.bdg",
        "#3   Write bedGraph files for control lambda (after scaling if "
        "necessary)... t_control_lambda.bdg",
        "#3   Pileup will be based on sequencing depth in treatment.",
        "#3 Call peaks for each chromosome...",
        "#4 Write output xls file... %s" % _out(out, "t", "_peaks.xls"),
        "#4 Write peak in narrowPeak format file... %s"
        % _out(out, "t", "_peaks.narrowPeak"),
        "#4 Write summits bed file... %s" % _out(out, "t", "_summits.bed"),
        "Done!"]


def test_no_peaks_gives_empty_files(run_macs3, tmp_path, write_bed):
    # q at the block is 31.43; -q 1e-40 asks for 40, so nothing passes
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["-q", "1e-40"], out))
    assert _lines(_out(out, "t", "_peaks.narrowPeak")) == []
    assert _lines(_out(out, "t", "_summits.bed")) == []
    assert _xls_table(_out(out, "t", "_peaks.xls"))[1] == []


# ------------------------------------
# options acting on the treatment pileup
# ------------------------------------

@pytest.mark.parametrize("extsize", [50, 100, 150])
def test_extsize(run_macs3, tmp_path, write_bed, extsize):
    # block [1000, 1000 + E); lambda = max(20 E / 1e4, 20 E / 1e5) = E / 500
    out = tmp_path / "out"
    args = ["-t", write_bed(_block(20)), "--nomodel", "--extsize",
            str(extsize), "-g", "100000", "--keep-dup", "all", "-B"]
    proc = _ok(_cp(run_macs3, args, out))
    end = 1000 + extsize
    lam = extsize / 500.0
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000", "chr1\t1000\t%d\t20.00000" % end]
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t%d\t%.5f" % (end, lam)]
    p = _pscore(20, lam)
    q = _ref_qscores([(p, extsize), (_pscore(0, lam), 1000)])[p]
    _assert_fields(_rows(_out(out, "t", "_peaks.narrowPeak"))[0],
                   ["chr1", "1000", str(end), "t_peak_1", str(int(10 * q)),
                    ".", 21 / (1 + _f32(lam)), p, q, str(extsize // 2)])
    assert "# d = %d" % extsize in _lines(_out(out, "t", "_peaks.xls"))
    assert "#2 Use %d as fragment length" % extsize in proc.stderr


@pytest.mark.parametrize("shift,strand,block", [
    (-50, "+", (950, 1050)),
    (0, "+", (1000, 1100)),
    (50, "+", (1050, 1150)),
    # a - read 1050-1100 has its 5' end at 1100; a positive shift moves it
    # towards its 3' end, i.e. to the left
    (50, "-", (950, 1050)),
    (-50, "-", (1050, 1150)),
], ids=["minus50-plus", "zero-plus", "plus50-plus", "plus50-minus",
        "minus50-minus"])
def test_shift_moves_block(run_macs3, tmp_path, write_bed, shift, strand,
                           block):
    reads = _block(20, start=1000 if strand == "+" else 1050, strand=strand)
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(reads)] + SYN
            + ["--shift", str(shift), "-B"], out))
    s, e = block
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t%d\t0.00000" % s, "chr1\t%d\t%d\t20.00000" % (s, e)]
    row = _rows(_out(out, "t", "_peaks.narrowPeak"))[0]
    assert row[1:3] == [str(s), str(e)]


@pytest.mark.parametrize("shift,line", [
    (25, "# Sequencing ends will be shifted towards 3' by 25 bp(s)"),
    (-25, "# Sequencing ends will be shifted towards 5' by 25 bp(s)"),
    (0, None),
])
def test_shift_reported(run_macs3, tmp_path, write_bed, parse_log, shift,
                        line):
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
                   + ["--shift", str(shift)], out))
    xls = _lines(_out(out, "t", "_peaks.xls"))
    shift_lines = [ln for ln in xls if "shifted" in ln]
    msgs = [m for m in _messages(parse_log, proc, "INFO") if "shifted" in m]
    if line is None:
        assert shift_lines == [] and msgs == []
    else:
        assert shift_lines == [line]
        assert msgs == ["#2" + line[1:]]


# two blocks of 20 reads, [1000, 1100) and [1130, 1230), 30 bp apart;
# lambda = 0.4 (two llocal windows of 0.2), so both blocks pass the
# cutoff and the 30 bp gap (pileup 0) does not
def _two_blocks():
    return _block(20, start=1000) + _block(20, start=1130)


@pytest.mark.parametrize("extra,peaks", [
    ([], [(1000, 1230, 50)]),                       # max gap = tag size 50
    (["-s", "29"], [(1000, 1100, 50), (1130, 1230, 50)]),
    (["-s", "30"], [(1000, 1230, 50)]),
    (["--max-gap", "29"], [(1000, 1100, 50), (1130, 1230, 50)]),
    (["--max-gap", "30"], [(1000, 1230, 50)]),
    (["-s", "10", "--max-gap", "30"], [(1000, 1230, 50)]),
], ids=["default", "tsize29", "tsize30", "maxgap29", "maxgap30",
        "maxgap-overrides-tsize"])
def test_max_gap_and_tsize_merge_blocks(run_macs3, tmp_path, write_bed,
                                        extra, peaks):
    # merged peak: the two highest segments tie, the first is the summit
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_two_blocks())] + SYN + extra, out))
    rows = _rows(_out(out, "t", "_peaks.narrowPeak"))
    assert [(int(r[1]), int(r[2]), int(r[9])) for r in rows] == peaks
    assert [r[3] for r in rows] == ["t_peak_%d" % (i + 1)
                                    for i in range(len(peaks))]


@pytest.mark.parametrize("tsize", [29, 75])
def test_tsize_reported(run_macs3, tmp_path, write_bed, parse_log, tsize):
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
                   + ["-s", str(tsize)], out))
    assert "#1 tag size is determined as %d bps" % tsize in \
        _messages(parse_log, proc, "INFO")
    assert "# tag size is determined as %d bps" % tsize in \
        _lines(_out(out, "t", "_peaks.xls"))


@pytest.mark.parametrize("minlen,npeaks", [(100, 1), (101, 0)])
def test_min_length(run_macs3, tmp_path, write_bed, minlen, npeaks):
    # the block is 100 bp long
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["--min-length", str(minlen)], out))
    assert len(_rows(_out(out, "t", "_peaks.narrowPeak"))) == npeaks
    assert "# The minimum length of peaks = %d" % minlen in \
        _lines(_out(out, "t", "_peaks.xls"))


@pytest.mark.parametrize("minlen,npeaks", [(230, 1), (231, 0)])
def test_min_length_of_merged_peak(run_macs3, tmp_path, write_bed, minlen,
                                   npeaks):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_two_blocks())] + SYN
            + ["--min-length", str(minlen)], out))
    assert len(_rows(_out(out, "t", "_peaks.narrowPeak"))) == npeaks


# ------------------------------------
# cutoffs
# ------------------------------------

@pytest.mark.parametrize("qvalue,npeaks", [("1e-31", 1), ("1e-32", 0)])
def test_qvalue_cutoff(run_macs3, tmp_path, write_bed, qvalue, npeaks):
    # Q20 = 31.43 > 31 but < 32
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["-q", qvalue], out))
    assert len(_rows(_out(out, "t", "_peaks.narrowPeak"))) == npeaks
    assert "# qvalue cutoff = %.2e" % float(qvalue) in \
        _lines(_out(out, "t", "_peaks.xls"))


@pytest.mark.parametrize("pvalue,npeaks", [("1e-34", 1), ("1e-35", 0)])
def test_pvalue_cutoff(run_macs3, tmp_path, write_bed, parse_log, pvalue,
                       npeaks):
    # P20 = 34.47 > 34 but < 35
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
                   + ["-p", pvalue], out))
    assert len(_rows(_out(out, "t", "_peaks.narrowPeak"))) == npeaks
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# pvalue cutoff = %.2e" % float(pvalue) in xls
    assert "# qvalue cutoff" not in "\n".join(xls)
    assert "#3 Call peaks with given -log10pvalue cutoff: %.5f ..." % \
        (-math.log10(float(pvalue))) in _messages(parse_log, proc, "INFO")


def test_pvalue_score_columns(run_macs3, tmp_path, write_bed):
    # with -p, the narrowPeak score is int(10 * pscore) and the summit
    # BED score is the p-score
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["-p", "1e-5"], out))
    row = _rows(_out(out, "t", "_peaks.narrowPeak"))[0]
    assert row[4] == str(int(10 * P20))
    assert float(row[7]) == pytest.approx(P20, rel=1e-5)
    summit = _rows(_out(out, "t", "_summits.bed"))[0]
    assert float(summit[4]) == pytest.approx(P20, rel=1e-5)


@pytest.mark.parametrize("fe,npeaks", [("17.4", 1), ("17.6", 0)])
def test_fe_cutoff(run_macs3, tmp_path, write_bed, fe, npeaks):
    # the block's fold enrichment is 17.5
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["--fe-cutoff", fe], out))
    assert len(_rows(_out(out, "t", "_peaks.narrowPeak"))) == npeaks
    assert "# Additional cutoff on fold-enrichment is: %.2f" % float(fe) in \
        _lines(_out(out, "t", "_peaks.xls"))


# ------------------------------------
# lambda options without a control
# ------------------------------------

def test_nolambda(run_macs3, tmp_path, write_bed, parse_log):
    # the lambda is lambda_bg = 0.02 everywhere
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
                   + ["--nolambda", "-B"], out))
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t1100\t0.02000"]
    p = _pscore(20, 0.02)
    q = _ref_qscores([(p, 100), (_pscore(0, 0.02), 1000)])[p]
    _assert_fields(_rows(_out(out, "t", "_peaks.narrowPeak"))[0],
                   ["chr1", "1000", "1100", "t_peak_1", str(int(10 * q)),
                    ".", 21 / (1 + _f32(0.02)), p, q, "50"])
    assert "# local lambda is disabled!" in _lines(_out(out, "t",
                                                        "_peaks.xls"))
    info = _messages(parse_log, proc, "INFO")
    assert "# local lambda is disabled!" in info
    assert "#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!" in info


@pytest.mark.parametrize("llocal,lam_lines", [
    # window [1000 - L/2, 1000 + L/2) clipped at 0, value 20 * 100 / L
    (1000, ["chr1\t0\t500\t0.02000", "chr1\t500\t1100\t2.00000"]),
    (5000, ["chr1\t0\t1100\t0.40000"]),
    (20000, ["chr1\t0\t1100\t0.10000"]),
    (200000, ["chr1\t0\t1100\t0.02000"]),            # below lambda_bg
])
def test_llocal_without_control(run_macs3, tmp_path, write_bed, llocal,
                                lam_lines):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["--llocal", str(llocal), "-B"], out))
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == lam_lines
    assert "# Range for calculating regional lambda is: %d bps" % llocal in \
        _lines(_out(out, "t", "_peaks.xls"))


@pytest.mark.parametrize("gsize,lam", [
    ("10000", 0.2), ("100000", 0.02), ("1e6", 0.002),
    ("hs", 2000.0 / EFFECTIVEGS["hs"]), ("ce", 2000.0 / EFFECTIVEGS["ce"]),
])
def test_gsize_sets_background_lambda(run_macs3, tmp_path, write_bed, gsize,
                                      lam):
    # with --nolambda the lambda is d * total / gsize = 2000 / gsize
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20)), "--nomodel",
                        "--extsize", "100", "--keep-dup", "all",
                        "--nolambda", "-B", "-g", gsize], out))
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t1100\t%.5f" % _f32(lam)]


@pytest.mark.parametrize("gsize,text", [
    ("hs", "2.91e+09"), ("mm", "2.65e+09"), ("ce", "1.00e+08"),
    ("dm", "1.43e+08"), ("1000000000", "1.00e+09"), ("2.5e7", "2.50e+07"),
])
def test_gsize_header(run_macs3, tmp_path, write_bed, gsize, text):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20)), "--nomodel",
                        "-g", gsize], out))
    assert "# effective genome size = %s" % text in \
        _lines(_out(out, "t", "_peaks.xls"))


# ------------------------------------
# duplicates
# ------------------------------------

@pytest.mark.parametrize("keepdup,kept,gsize", [
    ("1", 1, "100000"),
    ("5", 5, "100000"),
    ("all", 20, "100000"),
    # auto: binomial 1 - 1e-5 quantile of Bin(20, 1/100) = 4
    ("auto", int(binom.ppf(1 - 1e-5, 20, 1 / 100.0)), "100"),
])
def test_keep_dup(run_macs3, tmp_path, write_bed, parse_log, keepdup, kept,
                  gsize):
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20)), "--nomodel",
                               "--extsize", "100", "-g", gsize, "-B",
                               "--keep-dup", keepdup], out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t1000\t1100\t%d.00000" % kept
    xls = _lines(_out(out, "t", "_peaks.xls"))
    info = _messages(parse_log, proc, "INFO")
    if keepdup == "all":
        assert not [ln for ln in xls if "after filtering" in ln]
        assert not [m for m in info if "filter out redundant" in m]
    else:
        assert "# tags after filtering in treatment: %d" % kept in xls
        assert ("# maximum duplicate tags at the same position in "
                "treatment = %d" % kept) in xls
        assert "# Redundant rate in treatment: %.2f" % ((20 - kept) / 20.0) \
            in xls
        assert ("#1 filter out redundant tags at the same location and the "
                "same strand by allowing at most %d tag(s)" % kept) in info
    if keepdup == "auto":
        assert "#1  max_dup_tags based on binomial = %d" % kept in info


def test_keep_dup_is_per_strand(run_macs3, tmp_path, write_bed):
    # 3 + and 3 - reads with the same coordinates: --keep-dup 1 keeps one
    # per strand; the + read covers [1000, 1100), the - read [950, 1050)
    reads = _block(3) + _block(3, strand="-")
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(reads), "--nomodel", "--extsize",
                        "100", "-g", "100000", "--keep-dup", "1", "-B"],
            out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t950\t0.00000", "chr1\t950\t1000\t1.00000",
        "chr1\t1000\t1050\t2.00000", "chr1\t1050\t1100\t1.00000"]
    assert "# tags after filtering in treatment: 2" in \
        _lines(_out(out, "t", "_peaks.xls"))


# ------------------------------------
# bedGraph options
# ------------------------------------

def test_spmr(run_macs3, tmp_path, write_bed, parse_log):
    # values divided by the treatment depth in millions (20 / 1e6)
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
                   + ["-B", "--SPMR"], out))
    denom = np.float32(20 / 1e6)
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000",
        "chr1\t1000\t1100\t%.5f" % np.float32(np.float32(20) / denom)]
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t1100\t%.5f" % np.float32(np.float32(0.2) / denom)]
    assert "# MACS will save fragment pileup signal per million reads" in \
        _lines(_out(out, "t", "_peaks.xls"))
    assert ("#3   --SPMR is requested, so pileup will be normalized by "
            "sequencing depth in million reads.") in \
        _messages(parse_log, proc, "INFO")


def test_spmr_does_not_change_peaks(run_macs3, tmp_path, write_bed):
    bed = write_bed(_block(20))
    _ok(_cp(run_macs3, ["-t", bed] + SYN, tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B", "--SPMR"], tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t",
                  ("_peaks.narrowPeak", "_summits.bed"))


def test_spmr_without_bdg_has_no_effect(run_macs3, tmp_path, write_bed):
    bed = write_bed(_block(20))
    _ok(_cp(run_macs3, ["-t", bed] + SYN, tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["--SPMR"], tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")
    assert _xls_without_cmdline(_out(tmp_path / "a", "t", "_peaks.xls")) == \
        _xls_without_cmdline(_out(tmp_path / "b", "t", "_peaks.xls"))


def test_trackline_narrow(run_macs3, tmp_path, write_bed):
    bed = write_bed(_block(20))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B", "--trackline"],
            tmp_path / "b"))
    a, b = tmp_path / "a", tmp_path / "b"
    narrow = _lines(_out(b, "t", "_peaks.narrowPeak"))
    assert narrow[0] == ('track type=narrowPeak name="t" description="t" '
                         'nextItemButton=on')
    assert narrow[1:] == _lines(_out(a, "t", "_peaks.narrowPeak"))
    summits = _lines(_out(b, "t", "_summits.bed"))
    assert re.fullmatch(r'track name="t \(summits\)" description="Summits '
                        r'for t \(Made with MACS v3, [^)]+\)" visibility=1',
                        summits[0])
    assert summits[1:] == _lines(_out(a, "t", "_summits.bed"))
    for suffix in BDG_SUFFIXES:
        bdg = _lines(_out(b, "t", suffix))
        assert bdg[0].startswith("track type=bedGraph name=")
        assert bdg[1:] == _lines(_out(a, "t", suffix))
    # the xls never gets a track line
    assert _xls_without_cmdline(_out(a, "t", "_peaks.xls")) == \
        _xls_without_cmdline(_out(b, "t", "_peaks.xls"))


def test_trackline_broad(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["--broad", "--trackline"], out))
    assert _lines(_out(out, "t", "_peaks.broadPeak"))[0] == (
        'track type=broadPeak name="t" description="t" nextItemButton=on')
    assert _lines(_out(out, "t", "_peaks.gappedPeak"))[0] == (
        'track name="t" description="t" type=gappedPeak nextItemButton=on')


# ------------------------------------
# names and output folders
# ------------------------------------

@pytest.mark.parametrize("name", ["x", "my.sample", "s_1-2"])
def test_name(run_macs3, tmp_path, write_bed, name):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN + ["-B"], out,
            name=name))
    assert sorted(p.name for p in out.iterdir()) == \
        sorted(name + s for s in NARROW_SUFFIXES + BDG_SUFFIXES)
    assert _rows(_out(out, name, "_peaks.narrowPeak"))[0][3] == \
        name + "_peak_1"
    assert _rows(_out(out, name, "_summits.bed"))[0][3] == name + "_peak_1"
    assert _xls_table(_out(out, name, "_peaks.xls"))[1][0][-1] == \
        name + "_peak_1"
    assert "# name = %s" % name in _lines(_out(out, name, "_peaks.xls"))


def test_default_name_is_NA(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(run_macs3(["callpeak", "-t", write_bed(_block(20))] + SYN
                  + ["--outdir", str(out)]))
    assert sorted(p.name for p in out.iterdir()) == \
        sorted("NA" + s for s in NARROW_SUFFIXES)
    assert _rows(_out(out, "NA", "_peaks.narrowPeak"))[0][3] == "NA_peak_1"


def test_outdir_nested_is_created(run_macs3, tmp_path, write_bed):
    out = tmp_path / "a" / "b" / "c"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN, out))
    assert sorted(p.name for p in out.iterdir()) == \
        sorted("t" + s for s in NARROW_SUFFIXES)


def test_without_outdir_writes_to_cwd(run_macs3, tmp_path, write_bed):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    bed = write_bed(_block(20))
    _ok(run_macs3(["callpeak", "-t", bed] + SYN + ["-n", "t"], cwd=cwd))
    assert sorted(p.name for p in cwd.iterdir()) == \
        sorted("t" + s for s in NARROW_SUFFIXES)


# ------------------------------------
# options with no effect on the outputs
# ------------------------------------

@pytest.mark.parametrize("extra", [
    ["--slocal", "300"],                 # slocal is used only with control
    ["--seed", "5"],                     # only with --down-sample
    ["--down-sample"],                   # only with a control
    ["--to-large"],                      # obsolete
    ["--ratio", "3.0"],                  # only with a control
    ["--scale-to", "large"],             # only with a control
    ["--fix-bimodal"],                   # only when building the model
    ["--bw", "100"],                     # only when building the model
    ["--d-min", "50"],
    ["-m", "10", "30"],
    ["--max-count", "2"],                # only for FRAG
    ["--buffer-size", "1"],
    ["--buffer-size", "3"],
    ["--verbose", "1"],
], ids=lambda x: "_".join(x).lstrip("-"))
def test_options_without_effect(run_macs3, tmp_path, write_bed, extra):
    reads = (_block(20) + _block(7, start=1500, strand="-")
             + _block(5, chrom="chr2", start=300))
    bed = write_bed(reads)
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B"] + extra, tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_barcodes_ignored_for_bed(run_macs3, tmp_path, write_bed):
    bed = write_bed(_block(20))
    bc = tmp_path / "bc.txt"
    bc.write_text("AAA\n")
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B", "--barcodes", bc],
            tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


@pytest.mark.parametrize("extra,line", [
    (["--ratio", "3.0"], "# Using a custom scaling factor: 3.00e+00"),
    (["--down-sample"],
     "# Larger dataset will be randomly sampled towards smaller dataset."),
    (["--scale-to", "large"],
     "# Smaller dataset will be scaled towards larger dataset."),
    (["--call-summits"], "# Searching for subpeak summits is on"),
    (["--max-gap", "77"],
     "# The maximum gap between significant sites = 77"),
], ids=["ratio", "down-sample", "scale-to-large", "call-summits", "max-gap"])
def test_header_lines(run_macs3, tmp_path, write_bed, extra, line):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN + extra, out))
    assert line in _lines(_out(out, "t", "_peaks.xls"))


def test_down_sample_seed_header(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
            + ["--down-sample", "--seed", "7"], out))
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# Random seed has been set as: 7" in xls


# ------------------------------------
# multiple chromosomes and multiple files
# ------------------------------------

def test_multiple_chromosomes(run_macs3, tmp_path, write_bed):
    # chr1: 20 reads at 1000, chr10: 10 at 2000, chr2: 15 at 3000;
    # lambda_bg = 100 * 45 / 1e5 = 0.045; local lambdas n / 100 per chrom
    reads = (_block(20, chrom="chr1", start=1000)
             + _block(10, chrom="chr10", start=2000)
             + _block(15, chrom="chr2", start=3000))
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(reads)] + SYN + ["-B"], out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000", "chr1\t1000\t1100\t20.00000",
        "chr10\t0\t2000\t0.00000", "chr10\t2000\t2100\t10.00000",
        "chr2\t0\t3000\t0.00000", "chr2\t3000\t3100\t15.00000"]
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t1100\t0.20000", "chr10\t0\t2100\t0.10000",
        "chr2\t0\t3100\t0.15000"]
    # q-scores are ranked over the whole genome
    blocks = [("chr1", 1000, 20, 0.2), ("chr10", 2000, 10, 0.1),
              ("chr2", 3000, 15, 0.15)]
    segs = []
    for _, s, k, lam in blocks:
        segs += [(_pscore(k, lam), 100), (_pscore(0, lam), s)]
    qtab = _ref_qscores(segs)
    rows = _rows(_out(out, "t", "_peaks.narrowPeak"))
    assert len(rows) == 3
    for i, (row, (chrom, s, k, lam)) in enumerate(zip(rows, blocks)):
        p = _pscore(k, lam)
        _assert_fields(row, [chrom, str(s), str(s + 100),
                             "t_peak_%d" % (i + 1), str(int(10 * qtab[p])),
                             ".", (k + 1) / (1 + _f32(lam)), p, qtab[p],
                             "50"], rel=2e-5)


def test_multiple_treatment_files_are_pooled(run_macs3, tmp_path, write_bed):
    a = write_bed(_block(12), name="a.bed")
    b = write_bed(_block(8) + _block(5, chrom="chr2"), name="b.bed")
    ab = write_bed(_block(12) + _block(8) + _block(5, chrom="chr2"),
                   name="ab.bed")
    _ok(_cp(run_macs3, ["-t", a, b] + SYN + ["-B"], tmp_path / "x"))
    _ok(_cp(run_macs3, ["-t", ab] + SYN + ["-B"], tmp_path / "y"))
    _same_outputs(tmp_path / "x", tmp_path / "y", "t")
    assert "# ChIP-seq file = ['%s', '%s']" % (a, b) in \
        _lines(_out(tmp_path / "x", "t", "_peaks.xls"))


# ------------------------------------
# with a control
# ------------------------------------
# treatment: 20 reads at chr1:1000. Control: 10 reads at 1300 and 10 at
# 5000 (+ strand), so both depths are 20 and the control is scaled to
# the treatment with ratio 1: control windows centred on 1300 and 5000:
#   d (100):      [1250, 1350) and [4950, 5050), 10 * 1 = 10
#   slocal 1000:  [800, 1800) and [4500, 5500), 10 * 100 / 1000 = 1
#   llocal 10000: [0, 6300) and [0, 10000), 10 * 100 / 10000 each = 0.2
#   lambda_bg = 100 * 20 / 1e5 = 0.02
# so the lambda is 0.2 on [0, 800) and 1 on [800, 1100) (the output stops
# at the treatment's last breakpoint, 1100).

def _control_reads(n_near=10, near=1300, n_far=10, far=5000):
    return _block(n_near, start=near) + _block(n_far, start=far)


def _run_with_control(run_macs3, tmp_path, write_bed, extra=(),
                      treat=None, ctrl=None, outname="out"):
    t = write_bed(treat if treat is not None else _block(20), name="t.bed")
    c = write_bed(ctrl if ctrl is not None else _control_reads(),
                  name="c.bed")
    out = tmp_path / outname
    proc = _ok(_cp(run_macs3, ["-t", t, "-c", c] + SYN + ["-B"]
                   + list(extra), out))
    return proc, out


def test_control_lambda_default(run_macs3, tmp_path, write_bed):
    _, out = _run_with_control(run_macs3, tmp_path, write_bed)
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000", "chr1\t1000\t1100\t20.00000"]
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t800\t0.20000", "chr1\t800\t1100\t1.00000"]
    # peak: pileup 20 against lambda 1
    p = _pscore(20, 1.0)
    q = _ref_qscores([(_pscore(0, 0.2), 800), (_pscore(0, 1.0), 200),
                      (p, 100)])[p]
    _assert_fields(_rows(_out(out, "t", "_peaks.narrowPeak"))[0],
                   ["chr1", "1000", "1100", "t_peak_1", str(int(10 * q)),
                    ".", 10.5, p, q, "50"])


def test_control_xls_header(run_macs3, tmp_path, write_bed):
    _, out = _run_with_control(run_macs3, tmp_path, write_bed)
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# control file = ['%s']" % (tmp_path / "c.bed") in xls
    assert ("# Range for calculating regional lambda is: 1000 bps and "
            "10000 bps") in xls
    assert "# total tags in control: 20" in xls


@pytest.mark.parametrize("extra,lam_lines", [
    (["--slocal", "2000"],          # [300, 2300), 10 * 100 / 2000 = 0.5
     ["chr1\t0\t300\t0.20000", "chr1\t300\t1100\t0.50000"]),
    (["--slocal", "0"],             # slocal skipped
     ["chr1\t0\t1100\t0.20000"]),
    (["--llocal", "20000"],         # 0.05 per window, 0.1 in total
     ["chr1\t0\t800\t0.10000", "chr1\t800\t1100\t1.00000"]),
    (["--llocal", "1000"],          # llocal not larger than slocal: skipped
     ["chr1\t0\t800\t0.02000", "chr1\t800\t1100\t1.00000"]),
    (["--slocal", "0", "--llocal", "0"],   # only the d window
     ["chr1\t0\t1100\t0.02000"]),
    (["--nolambda"],                # lambda_bg only
     ["chr1\t0\t1100\t0.02000"]),
    (["--ratio", "0.5"],            # every control window scaled by 0.5
     ["chr1\t0\t800\t0.10000", "chr1\t800\t1100\t0.50000"]),
], ids=["slocal2000", "slocal0", "llocal20000", "llocal1000",
        "slocal0-llocal0", "nolambda", "ratio0.5"])
def test_control_lambda_windows(run_macs3, tmp_path, write_bed, extra,
                                lam_lines):
    _, out = _run_with_control(run_macs3, tmp_path, write_bed, extra)
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == lam_lines
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000", "chr1\t1000\t1100\t20.00000"]


@pytest.mark.parametrize("extra,message", [
    (["--slocal", "50"], "AssertionError: 50 can't be smaller than 100!"),
    (["--llocal", "50"], "AssertionError: 50 can't be smaller than 100!"),
    (["--llocal", "500"], "AssertionError: 500 can't be smaller than 1000!"),
])
def test_control_window_smaller_than_d(run_macs3, tmp_path, write_bed,
                                       extra, message):
    t = write_bed(_block(20), name="t.bed")
    c = write_bed(_control_reads(), name="c.bed")
    proc = _cp(run_macs3, ["-t", t, "-c", c] + SYN + extra, tmp_path / "out")
    assert proc.returncode == 1
    assert proc.stderr.rstrip().splitlines()[-1] == message


# control of 40 reads (20 at 1300, 20 at 5000) against 20 treatment reads
@pytest.mark.parametrize("extra,treat_v,lam_lines,which", [
    # default: scale the control down by 20 / 40
    ([], 20, ["chr1\t0\t800\t0.20000", "chr1\t800\t1100\t1.00000"],
     "treatment"),
    # large: scale the treatment up by 2, control unscaled,
    # lambda_bg = 100 * 40 / 1e5 = 0.04
    (["--scale-to", "large"], 40,
     ["chr1\t0\t800\t0.40000", "chr1\t800\t1100\t2.00000"], "control"),
    # --to-large is obsolete and behaves like the default
    (["--to-large"], 20,
     ["chr1\t0\t800\t0.20000", "chr1\t800\t1100\t1.00000"], "treatment"),
], ids=["small", "large", "to-large"])
def test_scale_to(run_macs3, tmp_path, write_bed, parse_log, extra, treat_v,
                  lam_lines, which):
    proc, out = _run_with_control(run_macs3, tmp_path, write_bed, extra,
                                  ctrl=_control_reads(20, 1300, 20, 5000))
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t1000\t1100\t%d.00000" % treat_v
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == lam_lines
    assert "#3   Pileup will be based on sequencing depth in %s." % which in \
        _messages(parse_log, proc, "INFO")


def test_scale_to_small_when_treatment_is_larger(run_macs3, tmp_path,
                                                 write_bed, parse_log):
    # 40 treatment reads, 20 control reads: the treatment is scaled down
    # by 20 / 40 and the control is unscaled (lambda_bg from the control)
    proc, out = _run_with_control(run_macs3, tmp_path, write_bed,
                                  treat=_block(40))
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t1000\t1100\t20.00000"
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t800\t0.20000", "chr1\t800\t1100\t1.00000"]
    assert "#3   Pileup will be based on sequencing depth in control." in \
        _messages(parse_log, proc, "INFO")


def test_down_sample_control(run_macs3, tmp_path, write_bed, parse_log):
    # the 40 control reads at 1300 are sampled down to 20; any sample
    # gives the same track: slocal 20 * 100 / 1000 = 2, llocal 0.2
    proc, out = _run_with_control(run_macs3, tmp_path, write_bed,
                                  ["--down-sample", "--seed", "3"],
                                  ctrl=_block(40, start=1300))
    info = _messages(parse_log, proc, "INFO")
    assert "#3 User prefers to use random sampling instead of linear " \
        "scaling." in info
    assert "#3 MACS is random sampling control tags..." in info
    assert "#3 Random seed (3) is used." in info
    assert "#3 20 tags from control are kept" in info
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t800\t0.20000", "chr1\t800\t1100\t2.00000"]
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t1000\t1100\t20.00000"


def test_down_sample_treatment(run_macs3, tmp_path, write_bed, parse_log):
    proc, out = _run_with_control(run_macs3, tmp_path, write_bed,
                                  ["--down-sample", "--seed", "3"],
                                  treat=_block(40))
    info = _messages(parse_log, proc, "INFO")
    assert "#3 MACS is random sampling treatment tags..." in info
    assert "#3 20 Tags from treatment are kept" in info
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t1000\t1100\t20.00000"


def test_down_sample_without_seed_warns(run_macs3, tmp_path, write_bed,
                                        parse_log):
    proc, _ = _run_with_control(run_macs3, tmp_path, write_bed,
                                ["--down-sample"],
                                ctrl=_block(40, start=1300))
    assert "#3 Your results may not be reproducible due to the random " \
        "sampling!" in _messages(parse_log, proc, "WARNING")


def test_down_sample_seed_is_reproducible(run_macs3, tmp_path, write_bed):
    ctrl = [("chr1", 1300 + 7 * i, 1350 + 7 * i, "c", 0, "+")
            for i in range(40)] + _block(1, start=5000)
    _run_with_control(run_macs3, tmp_path, write_bed,
                      ["--down-sample", "--seed", "11"], ctrl=ctrl,
                      outname="a")
    _run_with_control(run_macs3, tmp_path, write_bed,
                      ["--down-sample", "--seed", "11"], ctrl=ctrl,
                      outname="b")
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_multiple_control_files_are_pooled(run_macs3, tmp_path, write_bed):
    t = write_bed(_block(20), name="t.bed")
    c1 = write_bed(_block(10, start=1300), name="c1.bed")
    c2 = write_bed(_block(10, start=5000), name="c2.bed")
    c12 = write_bed(_control_reads(), name="c12.bed")
    _ok(_cp(run_macs3, ["-t", t, "-c", c1, c2] + SYN + ["-B"], tmp_path / "x"))
    _ok(_cp(run_macs3, ["-t", t, "-c", c12] + SYN + ["-B"], tmp_path / "y"))
    _same_outputs(tmp_path / "x", tmp_path / "y", "t")


def test_control_only_common_chromosomes(run_macs3, tmp_path, write_bed):
    # chr2 is absent from the control, so it is not called
    treat = _block(20) + _block(20, chrom="chr2")
    _, out = _run_with_control(run_macs3, tmp_path, write_bed, treat=treat)
    assert {r[0] for r in _bdg(_out(out, "t", "_treat_pileup.bdg"))} == \
        {"chr1"}
    assert [r[0] for r in _rows(_out(out, "t", "_peaks.narrowPeak"))] == \
        ["chr1"]


def test_keep_dup_auto_uses_control_maximum(run_macs3, tmp_path, write_bed):
    """Regression test: with --keep-dup auto the control's own maximum
    duplicate count was computed and logged but the treatment's value was
    used to filter it and in the xls header.

    Fixed upstream in ed436c1 (#753).
    """
    # -g 1000: treatment 20 tags -> max 2; control 200 tags -> max 4
    tmax = int(binom.ppf(1 - 1e-5, 20, 1e-3))
    cmax = int(binom.ppf(1 - 1e-5, 200, 1e-3))
    assert (tmax, cmax) == (2, 4)
    treat = _block(10, start=1000) + _block(10, start=1500)
    ctrl = []
    for j in range(10):
        ctrl += _block(20, start=2000 + 100 * j)
    t = write_bed(treat, name="t.bed")
    c = write_bed(ctrl, name="c.bed")
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", t, "-c", c, "--nomodel", "-g", "1000",
                        "--keep-dup", "auto"], out))
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# tags after filtering in control: %d" % (10 * cmax) in xls
    assert ("# maximum duplicate tags at the same position in control = %d"
            % cmax) in xls


# ------------------------------------
# --cutoff-analysis
# ------------------------------------
# On the block (no control, lambda 0.2 over [0, 1100)) every cutoff c in
# 0.3, 0.6, ..., 9.9 lies below P20, so [1000, 1100) is one peak of 100
# bp. Cutoffs below P0 = 0.742 also include [0, 1000): the whole
# [0, 1100) becomes one peak of 1100 bp. The q-score of a cutoff is
# ranked together with the observed p-scores (the cutoffs cover 0 bp):
# q(c) = c + log10(101) - log10(1100) while positive.

def _cutoff_rows():
    cutoffs = [round(x, 5) for x in sorted(np.arange(0.3, 10.0, 0.3),
                                           reverse=True)]
    p0 = _pscore(0, 0.2)
    qtab = _ref_qscores([(P20, 100), (p0, 1000)] + [(c, 0) for c in cutoffs])
    rows = []
    for c in cutoffs:
        if c > p0:
            rows.append("%.2f\t%.2f\t1\t100\t100.00" % (c, qtab[c]))
        else:
            rows.append("%.2f\t%.2f\t1\t1100\t1100.00" % (c, qtab[c]))
    return rows


def test_cutoff_analysis_high_cutoffs(run_macs3, tmp_path, write_bed,
                                      parse_log):
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
                   + ["--cutoff-analysis"], out))
    lines = _lines(_out(out, "t", "_cutoff_analysis.txt"))
    assert lines[0] == "pscore\tqscore\tnpeaks\tlpeaks\tavelpeak"
    expected = [r for r in _cutoff_rows() if "\t100\t" in r]
    assert lines[1:len(expected) + 1] == expected
    assert "#3 Cutoff vs peaks called will be analyzed!" in \
        _messages(parse_log, proc, "INFO")


def test_cutoff_analysis_does_not_change_peaks(run_macs3, tmp_path,
                                               write_bed):
    bed = write_bed(_two_blocks() + _block(9, chrom="chr2", start=4000))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bed] + SYN + ["-B", "--cutoff-analysis"],
            tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_cutoff_analysis_failure_exit_code(run_macs3, tmp_path, write_bed):
    # the cutoff-analysis file cannot be opened (a directory is in the way)
    out = tmp_path / "out"
    (out / "t_cutoff_analysis.txt").mkdir(parents=True)
    proc = _cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
               + ["--cutoff-analysis"], out)
    assert proc.returncode == 1
    assert "IsADirectoryError" in proc.stderr.splitlines()[-1]


def test_tempdir_missing_fails(run_macs3, tmp_path, write_bed):
    # the pileup files go to --tempdir; a missing folder makes mkstemp fail
    proc = _cp(run_macs3, ["-t", write_bed(_block(20))] + SYN
               + ["--tempdir", tmp_path / "missing"], tmp_path / "out")
    assert proc.returncode == 1
    assert proc.stderr.splitlines()[-1].startswith("FileNotFoundError")


# ------------------------------------
# --call-summits
# ------------------------------------
# three ramps: 20 + reads at 1000, 1003, ..., 1057; 6 at 1150, 1160, ...,
# 1200; 20 at 1300, ..., 1357 (extsize 100): plateaus of 20 on
# [1057, 1100), of 6 on [1200, 1250) and of 20 on [1357, 1400); the whole
# [1003, 1454) is above the cutoff.

def _ramps():
    reads = [("chr1", 1000 + 3 * i, 1050 + 3 * i, "a", 0, "+")
             for i in range(20)]
    reads += [("chr1", 1150 + 10 * i, 1200 + 10 * i, "b", 0, "+")
              for i in range(6)]
    reads += [("chr1", 1300 + 3 * i, 1350 + 3 * i, "c", 0, "+")
              for i in range(20)]
    return reads


def test_call_summits_subpeaks(run_macs3, tmp_path, write_bed):
    """Re-derived for upstream 9597df4 (#750, issue #748), which fixed
    enforce_peakyness (the first peak no longer subtracts sqrt(threshold)
    twice, and hard_clip now clips on the left). Before it, the three
    ramps gave three sub-peaks (1a, 1b, 1c); now only the middle one passes.

    By hand, with pileup 1 below the cutoff (the peak starts at 1003, where
    the pileup reaches 2), the summit-search signal is 0 in the gaps
    [1157, 1160) and [1290, 1303). maxima() (smoothing d = 100) finds the
    two plateaus of 20 as adjacent pairs (1077/1078 and 1377/1378) and the
    middle hump at 1224, so enforce_peakyness sees minima at 1077 (20), in
    the first gap (0), in the second gap (0) and at 1377 (20):

    - the four maxima on the plateaus of 20 each have an adjacent minimum
      of 20, so the threshold is 20 + sqrt(20) and their whole region is
      negative: rejected;
    - the middle maximum has threshold 0 + sqrt(0) = 0 between the two
      gaps, so its region [1157, 1290) is nonnegative, 133 bp wide, with
      the 6 distinct values 0, 2, 3, 4, 5, 6: kept.

    A single sub-peak is reported without a letter suffix.
    """
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_ramps()), "--nomodel", "--extsize",
                        "100", "-g", "100000", "--call-summits"], out))
    rows = _rows(_out(out, "t", "_peaks.narrowPeak"))
    assert [r[3] for r in rows] == ["t_peak_1"]
    assert {(r[1], r[2]) for r in rows} == {("1003", "1454")}
    plateaus = [(1200, 1250)]
    summits = [1003 + int(r[9]) for r in rows]
    for summit, (s, e) in zip(summits, plateaus):
        assert s <= summit < e
    srows = _rows(_out(out, "t", "_summits.bed"))
    assert [int(r[1]) for r in srows] == summits
    assert [r[3] for r in srows] == ["t_peak_1"]
    _, xrows = _xls_table(_out(out, "t", "_peaks.xls"))
    assert [r[5] for r in xrows] == ["6"]
    assert [int(r[4]) for r in xrows] == [x + 1 for x in summits]


def test_without_call_summits_one_summit(run_macs3, tmp_path, write_bed):
    # the two plateaus of 20 tie; the first one's middle is the summit
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_ramps()), "--nomodel", "--extsize",
                        "100", "-g", "100000"], out))
    rows = _rows(_out(out, "t", "_peaks.narrowPeak"))
    assert [(r[1], r[2], r[3], r[9]) for r in rows] == \
        [("1003", "1454", "t_peak_1", str((1057 + 1100) // 2 - 1003))]


# ------------------------------------
# --broad and --broad-cutoff
# ------------------------------------
# 20 + reads at 1000, 2 at 1100 and 2 at 3000 (extsize 100); lambda is
# 24 * 100 / 10000 = 0.24 on [0, 3100); pileup 20 on [1000, 1100) and 2
# on [1100, 1200) and [3000, 3100). With -q 0.01 the strong level is
# q > 2 and with the default --broad-cutoff 0.1 the weak level is q > 1:
# q(20) = 30.1, q(2) = 1.23. Broad peaks are the weak regions; their
# values are length-weighted means; strong regions inside are blocks.

def _broad_reads():
    return _block(20) + _block(2, start=1100) + _block(2, start=3000)


def _broad_expect():
    lam = 0.24
    p20, p2, p0 = _pscore(20, lam), _pscore(2, lam), _pscore(0, lam)
    qtab = _ref_qscores([(p20, 100), (p2, 200), (p0, 2800)])
    fc20, fc2 = 21 / (1 + _f32(lam)), 3 / (1 + _f32(lam))
    peak1 = dict(q=(qtab[p20] + qtab[p2]) / 2, p=(p20 + p2) / 2,
                 fc=(fc20 + fc2) / 2, pileup=11)
    peak2 = dict(q=qtab[p2], p=p2, fc=fc2, pileup=2)
    return peak1, peak2


def test_broad_peaks(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_broad_reads())] + SYN
            + ["-q", "0.01", "--broad"], out))
    pk1, pk2 = _broad_expect()
    rows = _rows(_out(out, "t", "_peaks.broadPeak"))
    assert len(rows) == 2
    _assert_fields(rows[0], ["chr1", "1000", "1200", "t_peak_1",
                             str(int(10 * pk1["q"])), ".", pk1["fc"],
                             pk1["p"], pk1["q"]])
    _assert_fields(rows[1], ["chr1", "3000", "3100", "t_peak_2",
                             str(int(10 * pk2["q"])), ".", pk2["fc"],
                             pk2["p"], pk2["q"]])


def test_broad_gapped_peaks(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_broad_reads())] + SYN
            + ["-q", "0.01", "--broad"], out))
    pk1, pk2 = _broad_expect()
    rows = _rows(_out(out, "t", "_peaks.gappedPeak"))
    # peak 1: strong block [1000, 1100) plus a 1 bp block at the end;
    # peak 2 has no strong block: two 1 bp blocks at its ends
    _assert_fields(rows[0], ["chr1", "1000", "1200", "t_peak_1",
                             str(int(10 * pk1["q"])), ".", "0", "0", "0",
                             "2", "100,1", "0,199", pk1["fc"], pk1["p"],
                             pk1["q"]])
    _assert_fields(rows[1], ["chr1", "3000", "3100", "t_peak_2",
                             str(int(10 * pk2["q"])), ".", "0", "0", "0",
                             "2", "1,1", "0,99", pk2["fc"], pk2["p"],
                             pk2["q"]])


def test_broad_xls(run_macs3, tmp_path, write_bed):
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_broad_reads())] + SYN
            + ["-q", "0.01", "--broad"], out))
    pk1, pk2 = _broad_expect()
    header, rows = _xls_table(_out(out, "t", "_peaks.xls"))
    assert header == ["chr", "start", "end", "length", "pileup",
                      "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)",
                      "name"]
    _assert_fields(rows[0], ["chr1", "1001", "1200", "200", "11", pk1["p"],
                             pk1["fc"], pk1["q"], "t_peak_1"])
    _assert_fields(rows[1], ["chr1", "3001", "3100", "100", "2", pk2["p"],
                             pk2["fc"], pk2["q"], "t_peak_2"])
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# qvalue cutoff for narrow/strong regions = 1.00e-02" in xls
    assert "# qvalue cutoff for broad/weak regions = 1.00e-01" in xls
    assert "# Broad region calling is on" in xls


def test_broad_cutoff_equal_to_strong(run_macs3, tmp_path, write_bed):
    # --broad-cutoff 0.01 equals -q: the broad peak is the strong region
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", write_bed(_broad_reads())] + SYN
            + ["-q", "0.01", "--broad", "--broad-cutoff", "0.01"], out))
    lam = 0.24
    p20 = _pscore(20, lam)
    q20 = _ref_qscores([(p20, 100), (_pscore(2, lam), 200),
                        (_pscore(0, lam), 2800)])[p20]
    rows = _rows(_out(out, "t", "_peaks.gappedPeak"))
    _assert_fields(rows[0], ["chr1", "1000", "1100", "t_peak_1",
                             str(int(10 * q20)), ".", "0", "0", "0", "1",
                             "100", "0", 21 / (1 + _f32(lam)), p20, q20])
    assert len(rows) == 1


def test_broad_with_pvalue(run_macs3, tmp_path, write_bed, parse_log):
    # -p 1e-3 (strong: p > 3) and --broad-cutoff 0.5 (weak: p > 0.30103)
    # are p-value cutoffs. Pileup 0 against lambda 0.24 has p = 0.67, so
    # the whole [0, 3100) is one weak region; only [1000, 1100) is strong
    # (p(2) = 2.82). Values are means weighted by length over [0, 3100);
    # with -p the score column is int(10 * mean p-score).
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_broad_reads())] + SYN
                   + ["-p", "1e-3", "--broad", "--broad-cutoff", "0.5"],
                   out))
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# pvalue cutoff for narrow/strong regions = 1.00e-03" in xls
    assert "# pvalue cutoff for broad/weak regions = 5.00e-01" in xls
    assert ("#3 Call broad peaks with given level1 -log10pvalue cutoff and "
            "level2: 3.00000, 0.30103...") in \
        _messages(parse_log, proc, "INFO")
    lam = _f32(0.24)
    length = {20: 100, 2: 200, 0: 2800}
    pscores = {k: _pscore(k, lam) for k in length}
    qtab = _ref_qscores([(pscores[k], ln) for k, ln in length.items()])
    pmean = sum(pscores[k] * ln for k, ln in length.items()) / 3100
    qmean = sum(qtab[pscores[k]] * ln for k, ln in length.items()) / 3100
    fcmean = sum((k + 1) / (1 + lam) * ln for k, ln in length.items()) / 3100
    _assert_fields(_rows(_out(out, "t", "_peaks.broadPeak"))[0],
                   ["chr1", "0", "3100", "t_peak_1", str(int(10 * pmean)),
                    ".", fcmean, pmean, qmean])
    # the strong block sits inside 1 bp blocks at both ends
    _assert_fields(_rows(_out(out, "t", "_peaks.gappedPeak"))[0],
                   ["chr1", "0", "3100", "t_peak_1", str(int(10 * pmean)),
                    ".", "0", "0", "0", "3", "1,100,1", "0,1000,3099",
                    fcmean, pmean, qmean])
    # mean pileup (20 * 100 + 2 * 200) / 3100 rounded to 2 decimals
    _assert_fields(_xls_table(_out(out, "t", "_peaks.xls"))[1][0],
                   ["chr1", "1", "3100", "3100", "0.77", pmean, fcmean,
                    qmean, "t_peak_1"])


# ------------------------------------
# shifting model (CTCF data)
# ------------------------------------

CTCF_SE = ["-g", "52000000", "-t", "CHIP", "-c", "CTRL"]


MODEL_IDS = ["default", "bw200", "bw150", "dmin15", "mfold2-50"]


@pytest.mark.parametrize("extra,npairs", [
    ([], 469), (["--bw", "200"], 294), (["--bw", "150"], 162),
    (["--d-min", "15"], 469), (["-m", "2", "50"], 469),
], ids=MODEL_IDS)
def test_model_paired_peaks(run_macs3, tmp_path, test_dir, parse_log, extra,
                            npairs):
    """Number of paired peaks; the log and the xls report the same d
    and alternative d.

    Pins the current output. The paired-peak search on real reads has
    no practical hand derivation.
    """
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE) + extra, out))
    info = _messages(parse_log, proc, "INFO")
    assert "#2 Total number of paired peaks: %d" % npairs in info
    d = [m for m in info if m.startswith("#2 predicted fragment length is ")]
    altd = [m for m in info
            if m.startswith("#2 alternative fragment length(s) may be ")]
    assert len(d) == len(altd) == 1
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# d = %s" % d[0].split()[-2] in xls
    assert "# " + altd[0][3:] in xls


@pytest.mark.parametrize("bw", [300, 200])
def test_model_r_script(run_macs3, tmp_path, test_dir, bw):
    # the R script holds the +/- tag profiles (percentages summing to 100,
    # 4 * bw + 11 points), the cross-correlation at 4 * bw lags, and the
    # alternative d values. The xcorr labels are not checked: they run
    # from -2 bw to 2 bw (np.linspace) while the correlation lags run
    # from -2 bw + 1 to 2 bw.
    out = tmp_path / "out"
    _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE)
            + ["--bw", str(bw)], out))
    text = _out(out, "t", "_model.r").read_text()
    vec = {}
    for key in ("p", "m", "ycorr", "xcorr", "altd"):
        m = re.search(r"^%s\s*<- c\((.*)\)$" % key, text, re.M)
        vec[key] = [float(x) for x in m.group(1).split(",")]
    assert len(vec["p"]) == len(vec["m"]) == 4 * bw + 11
    assert sum(vec["p"]) == pytest.approx(100.0)
    assert sum(vec["m"]) == pytest.approx(100.0)
    assert len(vec["ycorr"]) == len(vec["xcorr"]) == 4 * bw
    lines = text.splitlines()
    assert lines[:2] == ["# R script for Peak Model",
                         "#  -- generated by MACS"]
    assert "pdf('t_model.pdf',height=6,width=6)" in lines
    assert lines[-1] == "dev.off()"
    d = int(re.search(r"^# d = (\d+)$", _out(out, "t", "_peaks.xls")
                      .read_text(), re.M).group(1))
    assert d in [int(x) for x in vec["altd"]]


@pytest.mark.parametrize("mfold,npairs", [(["10", "30"], 95),
                                          (["40", "50"], 5)])
def test_model_not_enough_pairs(run_macs3, tmp_path, test_dir, parse_log,
                                mfold, npairs):
    """Too few paired peaks stops the run with exit code 1.

    Pins the current output. The number of paired peaks found in real
    data within the MFOLD range has no practical hand derivation.
    """
    out = tmp_path / "out"
    proc = _cp(run_macs3, _ctcf_args(test_dir, CTCF_SE) + ["-m"] + mfold, out)
    assert proc.returncode == 1
    warn = _messages(parse_log, proc, "WARNING")
    assert warn == [
        "#2 MACS3 needs at least 100 paired peaks at + and - strand to build "
        "the model, but can only find %d! Please make your MFOLD range "
        "broader and try again. If MACS3 still can't build the model, we "
        "suggest to use --nomodel and --extsize 147 or other fixed number "
        "instead." % npairs,
        "#2 Process for pairing-model is terminated!"]
    assert not _out(out, "t", "_peaks.xls").exists()


def test_model_not_enough_pairs_synthetic(run_macs3, tmp_path, write_bed,
                                          parse_log):
    # one block has no + / - peak pair at all
    proc = _cp(run_macs3, ["-t", write_bed(_block(20)), "-g", "100000"],
               tmp_path / "out")
    assert proc.returncode == 1
    assert "#2 Total number of paired peaks: 0" in \
        _messages(parse_log, proc, "INFO")


def test_fix_bimodal_falls_back_to_extsize(run_macs3, tmp_path, write_bed,
                                           parse_log):
    # same input as above, but --fix-bimodal continues with d = extsize
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", write_bed(_block(20)), "-g", "100000",
                               "--keep-dup", "all", "--fix-bimodal",
                               "--extsize", "100"], out))
    warn = _messages(parse_log, proc, "WARNING")
    assert warn[-2:] == ["#2 Skipped...",
                         "#2 Since --fix-bimodal is set, MACS will use 100 "
                         "as fragment length"]
    assert "# d = 100" in _lines(_out(out, "t", "_peaks.xls"))
    assert not _out(out, "t", "_model.r").exists()
    assert _rows(_out(out, "t", "_peaks.narrowPeak"))[0][1:3] == \
        ["1000", "1100"]


@pytest.mark.parametrize("fix", [False, True])
def test_model_d_below_twice_tag_size(run_macs3, tmp_path, test_dir,
                                      parse_log, fix):
    """With -s 200 the predicted d (about 229) is below 2 * tag size:
    the warning names the predicted d, and --fix-bimodal replaces it by
    the tag size.

    Pins the current output. The warning branch is chosen from the
    model's prediction on real data.
    """
    out = tmp_path / "out"
    extra = ["-s", "200"] + (["--fix-bimodal"] if fix else [])
    proc = _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE) + extra, out))
    info = _messages(parse_log, proc, "INFO")
    d = int(re.fullmatch(r"#2 predicted fragment length is (\d+) bps",
                         [m for m in info if m.startswith(
                             "#2 predicted fragment length")][0]).group(1))
    assert d < 2 * 200
    warn = _messages(parse_log, proc, "WARNING")
    assert warn[0] == ("#2 Since the d (%d) calculated from paired-peaks are "
                       "smaller than 2*tag length, it may be influenced by "
                       "unknown sequencing problem!" % d)
    if fix:
        assert warn[1:] == [
            "#2 MACS will use 200 as EXTSIZE/fragment length d. NOTE: if the "
            "d calculated is still acceptable, please do not use "
            "--fix-bimodal option!"]
        assert "# d = 200" in _lines(_out(out, "t", "_peaks.xls"))
    else:
        assert warn[1:] == [
            "#2 You may need to consider one of the other alternative d(s): "
            "%d" % d,
            "#2 You can restart the process with --nomodel --extsize XXX "
            "with your choice or an arbitrary number. Nontheless, MACS will "
            "continute computing."]
        assert "# d = %d" % d in _lines(_out(out, "t", "_peaks.xls"))


def test_model_verbose3_debug_lines(run_macs3, tmp_path, test_dir,
                                    parse_log):
    """--verbose 3 adds the model summary at DEBUG level.

    Pins the current output. min_tags, max_tags and the peak counts come
    from the model built on real data. The chromosome line (a bytes
    repr) and the d and scan_window lines are checked only for their
    form and consistency.
    """
    proc = _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE)
                   + ["--verbose", "3"], tmp_path / "out"))
    debug = _messages(parse_log, proc, "DEBUG")
    assert debug[1].startswith("Chromosome: ") and "chr22" in debug[1]
    d = [m for m in _messages(parse_log, proc, "INFO")
         if m.startswith("#2 predicted fragment length is ")][0].split()[-2]
    assert debug[-2:] == ["#2   d: %s" % d,
                          "#2   scan_window: %d" % (2 * int(d))]
    debug = debug[:1] + debug[2:-2]
    assert debug == [
        "#2 min_tags: 1; max_tags:14; ",
        "Number of unique tags on + strand: 24079",
        "Number of peaks in + strand: 1726",
        "plus peaks: first - (16977271, 2.0) ... last - (51222024, 4.0)",
        "Number of unique tags on - strand: 23968",
        "Number of peaks in - strand: 1687",
        "minus peaks: first - (17255609, 7.0) ... last - (51213904, 8.0)",
        "ip_max: 1726; im_max: 1687",
        "Paired centers: first - 17255565 ... second - 51213799 ",
        "Number of paired peaks in this chromosome: 469",
        "start model_add_line...",
        "start X-correlation...",
        "#2  Summary Model:",
        "#2   min_tags: 1"]
    # the default level shows none of them
    proc2 = _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE),
                    tmp_path / "out2"))
    assert _messages(parse_log, proc2, "DEBUG") == []


def test_ctcf_outputs_consistent(run_macs3, tmp_path, test_dir):
    # narrowPeak, xls and summits describe the same peaks
    out = tmp_path / "out"
    _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE), out))
    narrow = _rows(_out(out, "t", "_peaks.narrowPeak"))
    summits = _rows(_out(out, "t", "_summits.bed"))
    _, xls = _xls_table(_out(out, "t", "_peaks.xls"))
    assert len(narrow) == len(summits) == len(xls) > 0
    for i, (n, s, x) in enumerate(zip(narrow, summits, xls)):
        name = "t_peak_%d" % (i + 1)
        start, end, off = int(n[1]), int(n[2]), int(n[9])
        assert n[3] == s[3] == x[9] == name
        assert 0 <= off < end - start
        assert int(s[1]) == start + off and int(s[2]) == start + off + 1
        assert x[:4] == [n[0], str(start + 1), str(end), str(end - start)]
        assert int(x[4]) == start + off + 1
        assert (x[6], x[7], x[8]) == (n[7], n[6], n[8])
        assert s[4] == n[8]
        q = float(n[8])
        assert q > -math.log10(0.05)
        assert abs(int(n[4]) - int(10 * q)) <= 1


def test_ctcf_verbose0_quiet_optvalidator_logger(run_macs3, tmp_path,
                                                 test_dir, parse_log):
    # --verbose 0 silences the messages of callpeak itself
    proc = _ok(_cp(run_macs3, _ctcf_args(test_dir, CTCF_SE)
                   + ["--verbose", "0"], tmp_path / "out"))
    msgs = _messages(parse_log, proc)
    assert not [m for m in msgs if m.startswith("#1") or m.startswith("#2")
                or m.startswith("#4") or m == "Done!"]


# ------------------------------------
# --verbose with FRAG input (FRAG always warns about --keep-dup)
# ------------------------------------


# ------------------------------------
# input formats: single end
# ------------------------------------

SE_FORMATS = ["BED", "BAM", "ELAND", "ELANDEXPORT", "BOWTIE"]


@pytest.mark.parametrize("fmt", SE_FORMATS)
def test_se_formats_give_same_pileup(run_macs3, tmp_path, make_alignments,
                                     fmt):
    path = _write_se_format(tmp_path, fmt, SE_READS, make_alignments)
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", path, "-f", fmt] + SYN + ["-B"], out))
    expected = _coverage_segments(_se_expected_intervals(SE_READS))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == \
        _bdg_lines("chr1", expected)
    assert "# tag size is determined as 36 bps" in \
        _lines(_out(out, "t", "_peaks.xls"))
    assert "# total tags in treatment: 24" in \
        _lines(_out(out, "t", "_peaks.xls"))


@pytest.mark.parametrize("fmt", SE_FORMATS[1:])
def test_se_formats_match_bed(run_macs3, tmp_path, make_alignments, fmt):
    bed = _write_se_format(tmp_path, "BED", SE_READS, make_alignments)
    path = _write_se_format(tmp_path, fmt, SE_READS, make_alignments)
    _ok(_cp(run_macs3, ["-t", bed, "-f", "BED"] + SYN + ["-B"],
            tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", path, "-f", fmt] + SYN + ["-B"],
            tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_sam_plus_strand_matches_bed(run_macs3, tmp_path, make_alignments):
    reads = [r for r in SE_READS if r[1] == "+"]
    bed = _write_se_format(tmp_path, "BED", reads, make_alignments)
    sam = _write_se_format(tmp_path, "SAM", reads, make_alignments)
    _ok(_cp(run_macs3, ["-t", bed, "-f", "BED"] + SYN + ["-B"],
            tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", sam, "-f", "SAM"] + SYN + ["-B"],
            tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


@pytest.mark.parametrize("fmt", ["BED", "BAM", "SAM", "ELAND", "ELANDEXPORT"])
def test_auto_detects_format(run_macs3, tmp_path, make_alignments, parse_log,
                             fmt):
    # SAM with + strand reads only: SAMParser raises TypeError on
    # minus-strand reads in this version
    reads = SE_READS if fmt != "SAM" else [r for r in SE_READS
                                           if r[1] == "+"]
    path = _write_se_format(tmp_path, fmt, reads, make_alignments)
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", path] + SYN + ["-B"], out))
    info = _messages(parse_log, proc, "INFO")
    assert "Detected format is: %s" % fmt in info
    assert ("* Input file is gzipped." in info) == (fmt == "BAM")
    expected = _coverage_segments(_se_expected_intervals(reads))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == \
        _bdg_lines("chr1", expected)


def test_auto_gzipped_bed(run_macs3, tmp_path, write_bed, parse_log):
    gz = write_bed(_block(20), name="reads.bed.gz")
    plain = write_bed(_block(20), name="reads.bed")
    proc = _ok(_cp(run_macs3, ["-t", gz] + SYN + ["-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", plain] + SYN + ["-B"], tmp_path / "b"))
    info = _messages(parse_log, proc, "INFO")
    assert "Detected format is: BED" in info
    assert "* Input file is gzipped." in info
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_bed_comment_lines_skipped(run_macs3, tmp_path, write_bed):
    rows = (["track name=x", "browser position chr1:1-100", "# comment"]
            + _block(20))
    _ok(_cp(run_macs3, ["-t", write_bed(rows, name="c.bed"), "-f", "BED"]
            + SYN + ["-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", write_bed(_block(20)), "-f", "BED"] + SYN
            + ["-B"], tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_bam_chromosome_length_clips_pileup(run_macs3, tmp_path,
                                            make_alignments):
    # chr1 is 1050 bp long in the BAM header; reads at 1000 extended by
    # 100 are clipped at 1050
    recs = [dict(name="r%d" % i, ref="chr1", pos=1000, flag=0, cigar="36M")
            for i in range(20)]
    bam = make_alignments(recs, refs=(("chr1", 1050),))
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", bam, "-f", "BAM"] + SYN + ["-B"], out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t1000\t0.00000", "chr1\t1000\t1050\t20.00000"]


def test_bam_skips_unmapped_and_secondary(run_macs3, tmp_path,
                                          make_alignments):
    # flags 4 (unmapped), 256 (secondary), 1024 is kept (duplicates are
    # MACS3's own business), 2048 (supplementary)
    recs = [dict(name="r%d" % i, ref="chr1", pos=1000, flag=0, cigar="36M")
            for i in range(10)]
    recs += [dict(name="u", ref="chr1", pos=1000, flag=4, cigar="36M"),
             dict(name="s", ref="chr1", pos=1000, flag=256, cigar="36M"),
             dict(name="x", ref="chr1", pos=1000, flag=2048, cigar="36M"),
             dict(name="d", ref="chr1", pos=1000, flag=1024, cigar="36M")]
    bam = make_alignments(recs)
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", bam, "-f", "BAM"] + SYN + ["-B"], out))
    assert "# total tags in treatment: 11" in \
        _lines(_out(out, "t", "_peaks.xls"))


# ------------------------------------
# input formats: paired end
# ------------------------------------

PE_FORMATS = ["BAMPE", "BEDPE", "FRAG"]


@pytest.mark.parametrize("fmt", PE_FORMATS)
def test_pe_formats_give_same_pileup(run_macs3, tmp_path, make_alignments,
                                     make_pe_pair, parse_log, fmt):
    path = _write_pe_format(tmp_path, fmt, PE_FRAGS, make_alignments,
                            make_pe_pair)
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", path, "-f", fmt, "-g", "10000", "-B"],
                   out))
    expected = _coverage_segments([(l, r, 1) for l, r in PE_FRAGS])
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == \
        _bdg_lines("chr1", expected)
    mean = sum(r - l for l, r in PE_FRAGS) / len(PE_FRAGS)    # 157.0
    info = _messages(parse_log, proc, "INFO")
    assert "#1 mean fragment size is determined as %.1f bp from treatment" \
        % mean in info
    assert "#1  total fragments in treatment: 15" in info
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# fragment size is determined as %d bps" % mean in xls
    assert "# d = %d" % mean in xls
    assert "# Paired-End mode is on" in xls


@pytest.mark.parametrize("fmt", ["BEDPE", "FRAG"])
def test_pe_formats_match_bampe(run_macs3, tmp_path, make_alignments,
                                make_pe_pair, fmt):
    bam = _write_pe_format(tmp_path, "BAMPE", PE_FRAGS, make_alignments,
                           make_pe_pair)
    path = _write_pe_format(tmp_path, fmt, PE_FRAGS, make_alignments,
                            make_pe_pair)
    _ok(_cp(run_macs3, ["-t", bam, "-f", "BAMPE", "-g", "10000", "-B"],
            tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", path, "-f", fmt, "-g", "10000", "-B"],
            tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_bampe_skips_improper_pairs(run_macs3, tmp_path, make_alignments,
                                    make_pe_pair):
    recs = []
    for i, (l, r) in enumerate(PE_FRAGS):
        recs += make_pe_pair("f%d" % i, "chr1", l, r, readlen=36)
    # a pair without the proper-pair flag (99 -> 97, 147 -> 145)
    bad = make_pe_pair("bad", "chr1", 5003, 5303, readlen=36)
    bad[0]["flag"], bad[1]["flag"] = 97, 145
    bam = make_alignments(recs + bad)
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", bam, "-f", "BAMPE", "-g", "10000", "-B"], out))
    assert "# total fragments in treatment: 15" in \
        _lines(_out(out, "t", "_peaks.xls"))


def test_bedpe_keep_dup_default_filters(run_macs3, tmp_path, write_bedpe,
                                        parse_log):
    # 5 identical fragments: one is kept by default
    path = write_bedpe([("chr1", 10000, 10200)] * 5)
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", path, "-f", "BEDPE", "-g", "100000",
                               "-B"], out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t10000\t10200\t1.00000"
    assert "#1 filter out redundant fragments by allowing at most 1 " \
        "identical fragment(s)" in _messages(parse_log, proc, "INFO")
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# fragments after filtering in treatment: 1" in xls
    assert "# maximum duplicate fragments in treatment = 1" in xls
    assert "# Redundant rate in treatment: 0.80" in xls


def test_frag_counts_equal_duplicate_bedpe_lines(run_macs3, tmp_path,
                                                 write_frag, write_bedpe):
    # a FRAG count of 5 is 5 identical BEDPE lines kept with --keep-dup all
    frag = write_frag([("chr1", 10000, 10200, "A", 5),
                       ("chr1", 20000, 20200, "B", 3)])
    bedpe = write_bedpe([("chr1", 10000, 10200)] * 5
                        + [("chr1", 20000, 20200)] * 3)
    _ok(_cp(run_macs3, ["-t", frag, "-f", "FRAG", "-g", "100000", "-B"],
            tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bedpe, "-f", "BEDPE", "--keep-dup", "all",
                        "-g", "100000", "-B"], tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_frag_max_count_one_equals_bedpe(run_macs3, tmp_path, write_frag,
                                         write_bedpe):
    # documented: "-f FRAG --max-count 1" == "-f BEDPE --keep-dup all"
    frag = write_frag([("chr1", 10000, 10200, "A", 5),
                       ("chr1", 20000, 20200, "B", 3)])
    bedpe = write_bedpe([("chr1", 10000, 10200), ("chr1", 20000, 20200)])
    _ok(_cp(run_macs3, ["-t", frag, "-f", "FRAG", "--max-count", "1",
                        "-g", "100000", "-B"], tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", bedpe, "-f", "BEDPE", "--keep-dup", "all",
                        "-g", "100000", "-B"], tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


@pytest.mark.parametrize("maxcount,heights,total", [
    (None, (5, 3), 8), ("0", (5, 3), 8), ("2", (2, 2), 4), ("4", (4, 3), 7),
])
def test_frag_max_count(run_macs3, tmp_path, write_frag, maxcount, heights,
                        total):
    frag = write_frag([("chr1", 10000, 10200, "A", 5),
                       ("chr1", 20000, 20200, "B", 3)])
    out = tmp_path / "out"
    extra = ["--max-count", maxcount] if maxcount else []
    _ok(_cp(run_macs3, ["-t", frag, "-f", "FRAG", "-g", "100000", "-B"]
            + extra, out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t10000\t0.00000",
        "chr1\t10000\t10200\t%d.00000" % heights[0],
        "chr1\t10200\t20000\t0.00000",
        "chr1\t20000\t20200\t%d.00000" % heights[1]]
    xls = _lines(_out(out, "t", "_peaks.xls"))
    assert "# total fragments in treatment: %d" % total in xls
    assert ("# Maximum count in fragment file is set as %s" % maxcount
            in xls) == (maxcount not in (None, "0"))


def test_frag_barcodes(run_macs3, tmp_path, write_frag, parse_log):
    frag = write_frag([("chr1", 10000, 10200, "A", 5),
                       ("chr1", 20000, 20200, "B", 3),
                       ("chr1", 30000, 30200, "A", 2)])
    bc = tmp_path / "bc.txt"
    bc.write_text("A\nZZZ\n")
    out = tmp_path / "out"
    proc = _ok(_cp(run_macs3, ["-t", frag, "-f", "FRAG", "-g", "100000",
                               "-B", "--barcodes", bc], out))
    info = _messages(parse_log, proc, "INFO")
    assert "#1 extract fragments with given barcodes" in info
    assert "#   extracted 7 fragments in treatment" in info
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t10000\t0.00000", "chr1\t10000\t10200\t5.00000",
        "chr1\t10200\t30000\t0.00000", "chr1\t30000\t30200\t2.00000"]
    # the xls reports the total before extraction
    assert "# total fragments in treatment: 10" in \
        _lines(_out(out, "t", "_peaks.xls"))


@pytest.mark.parametrize("extra", [["--shift", "30"], ["--extsize", "50"],
                                   ["-s", "20"], ["--bw", "100"]],
                         ids=lambda x: "_".join(x).lstrip("-"))
def test_pe_ignores_se_options(run_macs3, tmp_path, write_bedpe, extra):
    # paired-end mode forces --nomodel, neutralises --shift and takes d
    # and the tag size from the mean fragment length
    path = write_bedpe([(c, l, r) for c, (l, r) in
                        zip(["chr1"] * len(PE_FRAGS), PE_FRAGS)])
    _ok(_cp(run_macs3, ["-t", path, "-f", "BEDPE", "-g", "10000", "-B"],
            tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", path, "-f", "BEDPE", "-g", "10000", "-B"]
            + extra, tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")
    assert not _out(tmp_path / "b", "t", "_model.r").exists()


@pytest.mark.parametrize("fmt", PE_FORMATS)
def test_buffer_size_paired_end(run_macs3, tmp_path, make_alignments,
                                make_pe_pair, fmt):
    # growing the fragment arrays one entry at a time changes nothing
    frags = PE_FRAGS + [(l + 30000, r + 30000) for l, r in PE_FRAGS]
    path = _write_pe_format(tmp_path, fmt, frags, make_alignments,
                            make_pe_pair)
    _ok(_cp(run_macs3, ["-t", path, "-f", fmt, "-g", "10000", "-B"],
            tmp_path / "a"))
    _ok(_cp(run_macs3, ["-t", path, "-f", fmt, "-g", "10000", "-B",
                        "--buffer-size", "1"], tmp_path / "b"))
    _same_outputs(tmp_path / "a", tmp_path / "b", "t")


def test_pe_without_control_lambda(run_macs3, tmp_path, write_bedpe):
    # 10 fragments [10000, 10200): length 2000, lambda_bg = 2000 / 1e5 =
    # 0.02; local lambda: windows of 10000 bp around both ends, scaled by
    # 2000 / (10000 * 10 * 2) = 0.01, i.e. 20 * 0.01 = 0.2 near the
    # fragments
    path = write_bedpe([("chr1", 10000, 10200)] * 10)
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", path, "-f", "BEDPE", "--keep-dup", "all",
                        "-g", "100000", "-B"], out))
    assert _lines(_out(out, "t", "_treat_pileup.bdg")) == [
        "chr1\t0\t10000\t0.00000", "chr1\t10000\t10200\t10.00000"]
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t5000\t0.02000", "chr1\t5000\t5200\t0.10000",
        "chr1\t5200\t10200\t0.20000"]
    # pileup 10 against lambda 0.2; minimum length = d = 200
    p = _pscore(10, 0.2)
    q = _ref_qscores([(_pscore(0, 0.02), 5000), (_pscore(0, 0.1), 200),
                      (_pscore(0, 0.2), 4800), (p, 200)])[p]
    _assert_fields(_rows(_out(out, "t", "_peaks.narrowPeak"))[0],
                   ["chr1", "10000", "10200", "t_peak_1", str(int(10 * q)),
                    ".", 11 / (1 + _f32(0.2)), p, q, "100"])


# FRAG with a control: treatment (100000, 100200) x10 and (100500,
# 100700) x1; control (100000, 100200) x1 and (103000, 103200) x1.
# treatment total 11 > 2 * control total 4, so the treatment is scaled to
# the control; windows d = 200 (mean treatment fragment length), 1000,
# 10000 around each control fragment end, scaled by 1, 0.2, 0.02;
# lambda_bg = 4 * 200 / 1e5 = 0.008.

def _frag_control_reference(x):
    ends = [100000, 100200, 103000, 103200]
    lam = 0.008
    for d, scale in ((200, 1.0), (1000, 0.2), (10000, 0.02)):
        cov = sum(1 for e in ends if e - d // 2 <= x < e - d // 2 + d)
        lam = max(lam, scale * cov)
    return lam


def test_frag_control_lambda(run_macs3, tmp_path, write_frag):
    """Regression test: with FRAG input and a control, merging the
    d/slocal/llocal lambdas read positions through a raw pointer into a
    structured array, mixing position and value bits.

    Fixed upstream in 5456b02 (#737): PETrackII.pileup_a_chromosome_c now
    passes contiguous position arrays to over_two_pv_array.
    """
    t = write_frag([("chr1", 100000, 100200, "A", 10),
                    ("chr1", 100500, 100700, "B", 1)], name="t.tsv")
    c = write_frag([("chr1", 100000, 100200, "A", 1),
                    ("chr1", 103000, 103200, "A", 1)], name="c.tsv")
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", t, "-c", c, "-f", "FRAG", "-g", "100000",
                        "-B"], out))
    rows = _bdg(_out(out, "t", "_control_lambda.bdg"))
    # contiguous, increasing coordinates
    assert all(s < e for _, s, e, _ in rows)
    assert all(a[2] == b[1] for a, b in zip(rows[:-1], rows[1:]))
    # value at the middle of every segment of the reference step function
    breaks = [0, 95000, 95200, 98000, 98200, 99500, 99700, 99900, 100300,
              100500, 100700]
    for a, b in zip(breaks[:-1], breaks[1:]):
        x = (a + b) // 2
        assert _bdg_value(rows, "chr1", x) == \
            pytest.approx(_frag_control_reference(x), abs=1e-5), x


def test_frag_control_lambda_single_window(run_macs3, tmp_path, write_frag):
    # with --slocal 0 --llocal 0 only the d window is used, so no merge
    # happens: lambda = max(0.008, coverage of [e - 100, e + 100))
    t = write_frag([("chr1", 100000, 100200, "A", 10),
                    ("chr1", 100500, 100700, "B", 1)], name="t.tsv")
    c = write_frag([("chr1", 100000, 100200, "A", 1),
                    ("chr1", 103000, 103200, "A", 1)], name="c.tsv")
    out = tmp_path / "out"
    _ok(_cp(run_macs3, ["-t", t, "-c", c, "-f", "FRAG", "-g", "100000", "-B",
                        "--slocal", "0", "--llocal", "0"], out))
    assert _lines(_out(out, "t", "_control_lambda.bdg")) == [
        "chr1\t0\t99900\t0.00800", "chr1\t99900\t100300\t1.00000",
        "chr1\t100300\t100700\t0.00800"]
    # treatment scaled by control/treatment depth: 800 / 2200
    assert _lines(_out(out, "t", "_treat_pileup.bdg"))[1] == \
        "chr1\t100000\t100200\t%.5f" % _f32(10 * _f32(1 / (2200 / 800)))


# ------------------------------------
# option validation through the command line
# ------------------------------------

@pytest.mark.parametrize("extra,errors", [
    (["-g", "foo"], ["Error when interpreting --gsize option: foo",
                     "Available shortcuts of effective genome sizes are "
                     "hs,mm,ce,dm"]),
    (["--keep-dup", "foo"], ["--keep-dup should be 'auto', 'all' or an "
                             "integer!"]),
    (["--keep-dup", "-1"], ["--keep-dup should be 'auto', 'all' or an "
                            "integer!"]),
    (["--keep-dup", "1.5"], ["--keep-dup should be 'auto', 'all' or an "
                             "integer!"]),
    (["--nomodel", "--extsize", "0"], ["--extsize must >= 1!"]),
    (["--extsize", "-5"], ["--extsize must >= 1!"]),
    (["--broad", "--call-summits"], ["--broad can't be combined with "
                                     "--call-summits!"]),
    (["-m", "50", "5"], ["Upper limit of mfold should be greater than lower "
                         "limit!"]),
    (["-f", "FRAG", "--max-count", "-1"], ["--max-count can't be a negative "
                                           "value"]),
], ids=["gsize", "keepdup-word", "keepdup-negative", "keepdup-float",
        "extsize0", "extsize-negative", "broad-summits", "mfold",
        "maxcount-negative"])
def test_validation_errors(run_macs3, tmp_path, write_bed, parse_log, extra,
                           errors):
    proc = _cp(run_macs3, ["-t", write_bed(_block(20))] + extra,
               tmp_path / "out")
    assert proc.returncode == 1
    assert [(lv, m) for lv, m in parse_log(proc.stderr)] == \
        [("ERROR", e) for e in errors]
    assert list((tmp_path / "out").iterdir()) == []


@pytest.mark.parametrize("args,error", [
    ([], "the following arguments are required: -t/--treatment"),
    (["-t"], "argument -t/--treatment: expected at least one argument"),
    (["-t", "X", "-f", "bed"],
     "argument -f/--format: invalid choice: 'bed' (choose from AUTO, BAM, "
     "SAM, BED, ELAND, ELANDMULTI, ELANDEXPORT, BOWTIE, BAMPE, BEDPE, FRAG)"),
    (["-t", "X", "-q", "0.1", "-p", "0.1"],
     "argument -p/--pvalue: not allowed with argument -q/--qvalue"),
    (["-t", "X", "--scale-to", "medium"],
     "argument --scale-to: invalid choice: 'medium' (choose from large, "
     "small)"),
    (["-t", "X", "--extsize", "abc"],
     "argument --extsize: invalid int value: 'abc'"),
    (["-t", "X", "-m", "5"], "argument -m/--mfold: expected 2 arguments"),
    (["-t", "X", "--broad-cutoff", "x"],
     "argument --broad-cutoff: invalid float value: 'x'"),
    (["-t", "X", "-s", "1.5"], "argument -s/--tsize: invalid int value: "
     "'1.5'"),
    (["-t", "X", "--verbose", "high"],
     "argument --verbose: invalid int value: 'high'"),
    (["-t", "X", "--buffer-size", "1e3"],
     "argument --buffer-size: invalid int value: '1e3'"),
], ids=["missing-t", "empty-t", "bad-format", "q-and-p", "bad-scale-to",
        "extsize-type", "mfold-one", "broad-cutoff-type", "tsize-type",
        "verbose-type", "buffer-size-type"])
def test_argparse_errors(run_macs3, tmp_path, args, error):
    proc = run_macs3(["callpeak", "--outdir", str(tmp_path / "o")] + args)
    assert proc.returncode == 2
    assert USAGE in proc.stderr
    last = proc.stderr.rstrip().splitlines()[-1]
    # Python 3.12 quotes each choice in "(choose from ...)"; 3.13+ does not
    last = re.sub(r"\(choose from [^)]*\)",
                  lambda m: m.group(0).replace("'", ""), last)
    assert last == "macs3 callpeak: error: " + error
    assert not (tmp_path / "o").exists()


def test_argparse_unknown_option(run_macs3, tmp_path):
    # unknown options are reported by the top-level macs3 parser
    proc = run_macs3(["callpeak", "-t", "X", "--no-such-option", "--outdir",
                      str(tmp_path / "o")])
    assert proc.returncode == 2
    assert proc.stderr.startswith("usage: macs3 ")
    assert proc.stderr.rstrip().splitlines()[-1] == \
        "macs3: error: unrecognized arguments: --no-such-option"


def test_help_lists_every_option(run_macs3):
    proc = run_macs3(["callpeak", "-h"])
    assert proc.returncode == 0
    for opt in ["--treatment", "--control", "--format", "--gsize", "--tsize",
                "--keep-dup", "--barcodes", "--max-count", "--outdir",
                "--name", "--bdg", "--verbose", "--trackline", "--SPMR",
                "--nomodel", "--shift", "--extsize", "--bw", "--d-min",
                "--mfold", "--fix-bimodal", "--qvalue", "--pvalue",
                "--scale-to", "--down-sample", "--seed", "--tempdir",
                "--nolambda", "--slocal", "--llocal", "--max-gap",
                "--min-length", "--broad", "--broad-cutoff",
                "--cutoff-analysis", "--call-summits", "--fe-cutoff",
                "--to-large", "--ratio", "--buffer-size"]:
        assert opt in proc.stdout, opt


# ------------------------------------
# check_names
# ------------------------------------

class _Named:
    """Minimal stand-in for a track: only get_chr_names() is used."""

    def __init__(self, names):
        self._names = names

    def get_chr_names(self):
        return set(self._names)


@pytest.mark.parametrize("tnames,cnames", [
    ({"chr1"}, {"chr1"}),
    ({"chr1", "chr2"}, {"chr2", "chr3"}),
    ({b"chr1", b"chrX"}, {b"chrX"}),
])
def test_check_names_common(tnames, cnames):
    messages = []
    assert check_names(_Named(tnames), _Named(cnames),
                       messages.append) is None
    assert messages == []


def test_check_names_none_common_exits():
    messages = []
    with pytest.raises(SystemExit):
        check_names(_Named({"chr2", "chr1"}), _Named({"chrX"}),
                    messages.append)
    assert messages == [
        "No common chromosome names can be found from treatment and "
        "control!",
        "Please make sure that the treatment and control alignment files "
        "were generated by using the same genome assembly!",
        "Chromosome names in treatment: chr1,chr2",
        "Chromosome names in control: chrX"]


def test_check_names_empty_control_exits():
    messages = []
    with pytest.raises(SystemExit):
        check_names(_Named({"chr1"}), _Named(set()), messages.append)
    assert messages[-1] == "Chromosome names in control: "


# ------------------------------------
# cal_max_dup_tags
# ------------------------------------

@pytest.mark.parametrize("gsize,n", [
    (1000, 20), (1000, 200), (100, 20), (1e9, 10 ** 6),
    (EFFECTIVEGS["hs"], 10 ** 7), (10, 5), (2, 10), (1e5, 0),
])
def test_cal_max_dup_tags(gsize, n):
    # the smallest x with P(X <= x) >= 1 - 1e-5 for X ~ Bin(n, 1/gsize)
    assert cal_max_dup_tags(gsize, n) == int(binom.ppf(1 - 1e-5, n,
                                                        1.0 / gsize))


@pytest.mark.parametrize("p", [0.5, 1e-2, 1e-8])
def test_cal_max_dup_tags_pvalue(p):
    assert cal_max_dup_tags(1000, 500, p=p) == int(binom.ppf(1 - p, 500,
                                                             1e-3))


def test_cal_max_dup_tags_all_tags():
    # gsize 1: every tag sits at the one position
    assert cal_max_dup_tags(1, 7) == 7


# ------------------------------------
# load_tag_files_options / load_frag_files_options
# ------------------------------------

def _options(macs3_argparser, tmp_path, args):
    out = tmp_path / "opt_out"
    out.mkdir(exist_ok=True)
    ns = macs3_argparser.parse_args(["callpeak"] + [str(a) for a in args]
                                    + ["--outdir", str(out), "-n", "t"])
    return opt_validate_callpeak(ns)


def test_load_tag_files_options_treatment_only(macs3_argparser, tmp_path,
                                               make_alignments):
    bed = _write_se_format(tmp_path, "BED", SE_READS, make_alignments)
    opts = _options(macs3_argparser, tmp_path, ["-t", bed, "-f", "BED"])
    treat, control = load_tag_files_options(opts)
    assert control is None
    assert opts.tsize == 36
    assert treat.total == 24
    plus, minus = treat.get_locations_by_chr(b"chr1")
    assert list(plus) == [s for s, st in SE_READS if st == "+"]
    assert list(minus) == [s + 36 for s, st in SE_READS if st == "-"]


def test_load_tag_files_options_keeps_given_tsize(macs3_argparser, tmp_path,
                                                  write_bed):
    bed = write_bed(_block(20))
    opts = _options(macs3_argparser, tmp_path, ["-t", bed, "-s", "77"])
    load_tag_files_options(opts)
    assert opts.tsize == 77


def test_load_tag_files_options_multiple_files(macs3_argparser, tmp_path,
                                               write_bed):
    t1 = write_bed(_block(3), name="t1.bed")
    t2 = write_bed(_block(4, chrom="chr2"), name="t2.bed")
    c1 = write_bed(_block(5), name="c1.bed")
    c2 = write_bed(_block(6, chrom="chr3"), name="c2.bed")
    opts = _options(macs3_argparser, tmp_path,
                    ["-t", t1, t2, "-c", c1, c2, "-f", "BED"])
    treat, control = load_tag_files_options(opts)
    assert treat.total == 7 and control.total == 11
    assert treat.get_chr_names() == {b"chr1", b"chr2"}
    assert control.get_chr_names() == {b"chr1", b"chr3"}


def test_load_tag_files_options_bam_lengths(macs3_argparser, tmp_path,
                                            make_alignments):
    bam = _write_se_format(tmp_path, "BAM", SE_READS, make_alignments)
    opts = _options(macs3_argparser, tmp_path, ["-t", bam, "-f", "BAM"])
    treat, _ = load_tag_files_options(opts)
    assert treat.get_rlengths()[b"chr1"] == 100000


def test_load_frag_files_options_bedpe(macs3_argparser, tmp_path,
                                       write_bedpe):
    path = write_bedpe([("chr1", l, r) for l, r in PE_FRAGS])
    opts = _options(macs3_argparser, tmp_path, ["-t", path, "-f", "BEDPE"])
    treat, control = load_frag_files_options(opts)
    assert control is None
    assert treat.total == 15
    assert opts.tsize == pytest.approx(157.0)
    locs = treat.get_locations_by_chr(b"chr1")
    assert [(int(a), int(b)) for a, b in zip(locs["l"], locs["r"])] == \
        PE_FRAGS


def test_load_frag_files_options_frag_max_count(macs3_argparser, tmp_path,
                                                write_frag):
    t1 = write_frag([("chr1", 100, 300, "A", 5)], name="t1.tsv")
    t2 = write_frag([("chr1", 500, 600, "B", 9)], name="t2.tsv")
    c1 = write_frag([("chr1", 100, 300, "A", 4)], name="c1.tsv")
    opts = _options(macs3_argparser, tmp_path,
                    ["-t", t1, t2, "-c", c1, "-f", "FRAG", "--max-count",
                     "3"])
    treat, control = load_frag_files_options(opts)
    # counts capped at 3
    assert treat.total == 6 and control.total == 3
    assert treat.length == 3 * 200 + 3 * 100


def test_load_frag_files_options_bampe(macs3_argparser, tmp_path,
                                       make_alignments, make_pe_pair):
    bam = _write_pe_format(tmp_path, "BAMPE", PE_FRAGS, make_alignments,
                           make_pe_pair)
    opts = _options(macs3_argparser, tmp_path, ["-t", bam, "-c", bam,
                                                "-f", "BAMPE"])
    treat, control = load_frag_files_options(opts)
    assert treat.total == control.total == 15
    assert opts.tsize == pytest.approx(157.0)


# ------------------------------------
# opt_validate_callpeak as used by run()
# ------------------------------------

def test_validated_defaults(macs3_argparser, tmp_path, write_bed):
    bed = write_bed(_block(1))
    opts = _options(macs3_argparser, tmp_path, ["-t", bed])
    out = str(tmp_path / "opt_out")
    assert opts.gsize == EFFECTIVEGS["hs"]
    assert opts.parser.__name__ == "guess_parser"
    assert opts.log_qvalue == pytest.approx(-math.log10(0.05))
    assert opts.log_pvalue is None
    assert (opts.lmfold, opts.umfold) == (5, 50)
    assert opts.peakxls == os.path.join(out, "t_peaks.xls")
    assert opts.peakNarrowPeak == os.path.join(out, "t_peaks.narrowPeak")
    assert opts.summitbed == os.path.join(out, "t_summits.bed")
    assert opts.bdg_treat == os.path.join(out, "t_treat_pileup.bdg")
    assert opts.bdg_control == os.path.join(out, "t_control_lambda.bdg")
    assert opts.modelR == os.path.join(out, "t_model.r")
    assert opts.cutoff_analysis_file == "None"


@pytest.mark.parametrize("args,attrs", [
    (["-p", "0.001"], dict(log_pvalue=3.0, log_qvalue=None)),
    (["--broad", "--broad-cutoff", "0.01"], dict(log_broadcutoff=2.0)),
    (["-m", "10", "30"], dict(lmfold=10, umfold=30)),
    (["-g", "mm"], dict(gsize=EFFECTIVEGS["mm"])),
    (["-g", "1.5e6"], dict(gsize=1.5e6)),
    (["-f", "BEDPE", "--shift", "40"], dict(nomodel=True, shift=0)),
    (["-f", "BAMPE"], dict(nomodel=True, gzip_flag=True)),
    (["-f", "FRAG", "--keep-dup", "3"], dict(nomodel=True,
                                             keepduplicates="all")),
    (["-f", "BAM"], dict(gzip_flag=True, nomodel=False)),
    (["--shift", "-40"], dict(shift=-40)),
], ids=["pvalue", "broad-cutoff", "mfold", "gsize-mm", "gsize-float",
        "bedpe", "bampe", "frag", "bam", "shift"])
def test_validated_attributes(macs3_argparser, tmp_path, write_bed, args,
                              attrs):
    bed = write_bed(_block(1))
    opts = _options(macs3_argparser, tmp_path, ["-t", bed] + args)
    for k, v in attrs.items():
        if isinstance(v, float):
            assert getattr(opts, k) == pytest.approx(v), k
        else:
            assert getattr(opts, k) == v, k


def test_validated_cutoff_analysis_path(macs3_argparser, tmp_path, write_bed):
    opts = _options(macs3_argparser, tmp_path,
                    ["-t", write_bed(_block(1)), "--cutoff-analysis"])
    assert opts.cutoff_analysis_file == \
        os.path.join(str(tmp_path / "opt_out"), "t_cutoff_analysis.txt")


def test_validated_unknown_format_exits(macs3_argparser, tmp_path, write_bed):
    # not reachable through argparse's choices; checked when the format
    # attribute is set directly
    ns = macs3_argparser.parse_args(["callpeak", "-t", write_bed(_block(1))])
    ns.format = "XYZ"
    with pytest.raises(SystemExit) as excinfo:
        opt_validate_callpeak(ns)
    assert excinfo.value.code == 1


# ------------------------------------
# run() called in-process
# ------------------------------------

def test_run_in_process(macs3_argparser, tmp_path, write_bed, monkeypatch):
    monkeypatch.setattr(tempfile, "tempdir", tempfile.tempdir)
    out = tmp_path / "out"
    out.mkdir()
    ns = macs3_argparser.parse_args(["callpeak", "-t", write_bed(_block(20))]
                                    + SYN + ["-n", "t", "--outdir", str(out),
                                             "--tempdir", str(tmp_path)])
    assert run(ns) is None
    _assert_fields(_rows(_out(out, "t", "_peaks.narrowPeak"))[0],
                   ["chr1", "1000", "1100", "t_peak_1", str(int(10 * Q20)),
                    ".", 17.5, P20, Q20, "50"])
    assert ns.PE_MODE is False and ns.d == 100


def test_run_in_process_validation_exit(macs3_argparser, tmp_path,
                                        write_bed):
    ns = macs3_argparser.parse_args(["callpeak", "-t", write_bed(_block(20)),
                                     "-g", "nope", "--outdir",
                                     str(tmp_path)])
    with pytest.raises(SystemExit) as excinfo:
        run(ns)
    assert excinfo.value.code == 1


def test_run_in_process_not_enough_pairs_exit(macs3_argparser, tmp_path,
                                              write_bed, monkeypatch):
    monkeypatch.setattr(tempfile, "tempdir", tempfile.tempdir)
    ns = macs3_argparser.parse_args(["callpeak", "-t", write_bed(_block(20)),
                                     "-g", "100000", "--outdir",
                                     str(tmp_path), "--tempdir",
                                     str(tmp_path)])
    with pytest.raises(SystemExit) as excinfo:
        run(ns)
    assert excinfo.value.code == 1


def test_run_uses_tempdir(macs3_argparser, tmp_path, write_bed, monkeypatch):
    # the per-chromosome pileup files are created in --tempdir and removed
    monkeypatch.setattr(tempfile, "tempdir", tempfile.tempdir)
    td = tmp_path / "td"
    td.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    removed = []
    real_unlink = os.unlink

    def spy(path, *a, **kw):
        removed.append(Path(os.fsdecode(path)))
        return real_unlink(path, *a, **kw)

    monkeypatch.setattr(os, "unlink", spy)
    reads = _block(20) + _block(10, chrom="chr2")
    ns = macs3_argparser.parse_args(["callpeak", "-t", write_bed(reads)]
                                    + SYN + ["-n", "t", "--outdir", str(out),
                                             "--tempdir", str(td)])
    run(ns)
    assert len([p for p in removed if p.parent == td]) == 2
    assert list(td.iterdir()) == []
