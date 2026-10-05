#!/usr/bin/env python

"""Module Description: Test the `macs3 hmmratac` command and the public
output functions in MACS3.Commands.hmmratac_cmd.

Most command-line runs use a small synthetic ATAC-seq-like BEDPE file
(open regions full of short fragments, flanked by positioned
nucleosomes, over a sparse background of mono-, di- and
tri-nucleosomal fragments) on which `hmmratac` finishes in well under a
second. They call `run()` in the test process with arguments parsed by
the real `macs3` argument parser, which is what `bin/macs3` does; the
exit codes, log messages and argument errors are also checked through
`bin/macs3` in a subprocess. Runs on the real yeast data take about
11 s each and are marked slow.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import gzip
import io
import json
import logging
import os
import re
import sys
import tempfile
import traceback
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import norm

from MACS3.Commands.hmmratac_cmd import (run as hmmratac_run,
                                         save_proba_to_bedGraph,
                                         save_states_bed,
                                         generate_states_path,
                                         save_accessible_regions)
from MACS3.Signal.BedGraph import bedGraphTrackI

# ------------------------------------
# Helpers: synthetic data
# ------------------------------------

# cutoffs suited to the synthetic data (fold changes up to ~10)
TUNE = ["-l", "4", "-u", "20", "-c", "2"]

LAYOUT = {b"chr1": (60000, [8000, 16000, 24000, 32000, 40000, 48000]),
          b"chr2": (30000, [7000, 15000, 23000])}


def make_atac_fragments(seed=2026):
    """Synthetic ATAC-seq fragments as sorted (chrom, left, right).

    Each open region centred at c has short fragments in [c-120,
    c+120), positioned mono-nucleosome fragments centred at c-230 and
    c+230 and a few di-nucleosome fragments spanning the centre. The
    background is sparse and mostly nucleosomal. A few fragments are
    repeated exactly so that duplicates exist.
    """
    rs = np.random.RandomState(seed)
    frags = []
    for chrom, (size, centers) in LAYOUT.items():
        for _ in range(size // 100):
            kind = rs.rand()
            if kind < 0.15:
                ln = rs.normal(55, 10)
            elif kind < 0.70:
                ln = rs.normal(200, 20)
            elif kind < 0.92:
                ln = rs.normal(400, 25)
            else:
                ln = rs.normal(600, 30)
            ln = int(max(25, ln))
            s = int(rs.randint(0, size - ln))
            frags.append((chrom, s, s + ln))
        for c in centers:
            for _ in range(200):
                ln = int(np.clip(rs.normal(55, 12), 25, 100))
                s = int(rs.randint(c - 120, c + 120 - ln))
                frags.append((chrom, s, s + ln))
            for side in (-1, 1):
                for _ in range(45):
                    ln = int(np.clip(rs.normal(190, 12), 150, 240))
                    mid = c + side * 230 + int(rs.normal(0, 15))
                    frags.append((chrom, mid - ln // 2, mid - ln // 2 + ln))
            for _ in range(12):
                ln = int(np.clip(rs.normal(400, 20), 330, 470))
                mid = c + int(rs.normal(0, 40))
                frags.append((chrom, mid - ln // 2, mid - ln // 2 + ln))
    dup = rs.choice(len(frags), 25, replace=False)
    frags.extend([frags[i] for i in dup])
    frags.sort()
    return frags


def write_rows(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(x.decode() if isinstance(x, bytes) else str(x)
                               for x in row) + "\n")
    return str(path)


def coverage(frags):
    """Per-chromosome fragment coverage from 0 to the last fragment end."""
    out = {}
    for chrom in sorted(set(c for c, _, _ in frags)):
        fr = [(l, r) for c, l, r in frags if c == chrom]
        end = max(r for _, r in fr)
        d = np.zeros(end + 1, dtype=np.int64)
        for l, r in fr:
            d[l] += 1
            d[r] -= 1
        out[chrom] = np.cumsum(d)[:end]
    return out


def fold_change(frags):
    """Fold change over the average pileup, as hmmratac computes it.

    The average is the total fragment length over the summed lengths of
    the per-chromosome pileups (from 0 to the last fragment end), in
    float32; each value is divided in float64 and stored as float32.
    """
    cov = coverage(frags)
    total = sum(int(a.sum()) for a in cov.values())
    n = sum(len(a) for a in cov.values())
    mean = np.float32(total) / np.float32(n)
    fc = {c: (a.astype(np.float64) / np.float64(mean)).astype(np.float32)
          for c, a in cov.items()}
    return fc, mean


def ref_call_peaks(fc, cutoff, min_length, max_gap):
    """Reference for bedGraphTrackI.call_peaks on per-base arrays.

    Runs of bases with fc >= cutoff are merged when the gap between
    them is at most ``max_gap``; merged peaks shorter than
    ``min_length`` are dropped. Returns (chrom, start, end, max fc).
    """
    out = []
    cutoff = np.float32(cutoff)
    for chrom in sorted(fc):
        a = fc[chrom]
        above = np.concatenate([[False], a >= cutoff, [False]])
        d = np.diff(above.astype(np.int8))
        merged = []
        for s, e in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
            if merged and s - merged[-1][1] <= max_gap:
                merged[-1][1] = int(e)
            else:
                merged.append([int(s), int(e)])
        for s, e in merged:
            if e - s >= min_length:
                out.append((chrom, s, e, float(a[s:e].max())))
    return out


def ref_cutoff_analysis(fc, min_length, max_gap, max_score=100, steps=100):
    """Reference for the cutoff analysis report on per-base arrays.

    Cutoffs run from the smallest fold change (at least 0) towards
    min(largest fold change, max_score) in ``steps`` float32 steps,
    rounded to 3 decimals. For each cutoff, runs of bases with fc above
    it are merged across gaps up to ``max_gap``; merged peaks of at
    least ``min_length`` are counted. Rows with peaks are reported from
    the highest cutoff down.
    """
    minv = max(0.0, float(min(a.min() for a in fc.values())))
    maxv = min(float(max(a.max() for a in fc.values())),
               float(np.float32(max_score)))
    step = float(np.float32((maxv - minv) / steps))
    rows = []
    for cutoff in [round(v, 3) for v in np.arange(minv, maxv, step)]:
        c32 = np.float32(cutoff)
        n = total = 0
        for chrom in sorted(fc):
            above = np.concatenate([[False], fc[chrom] > c32, [False]])
            d = np.diff(above.astype(np.int8))
            merged = []
            for s, e in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
                if merged and s - merged[-1][1] <= max_gap:
                    merged[-1][1] = int(e)
                else:
                    merged.append([int(s), int(e)])
            for s, e in merged:
                if e - s >= min_length:
                    n += 1
                    total += e - s
        if n:
            rows.append("%.2f\t%d\t%d\t%.2f" % (float(c32), n, total,
                                                total / n))
    return ["score\tnpeaks\tlpeaks\tavelpeak"] + rows[::-1]


def ref_training_regions(peaks, flanking):
    """Expand peaks by ``flanking`` on both sides and merge overlaps."""
    out = []
    for chrom in sorted(set(p[0] for p in peaks)):
        regs = sorted((max(0, s - flanking), e + flanking)
                      for c, s, e, _ in peaks if c == chrom)
        merged = []
        for s, e in regs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = e
            else:
                merged.append([s, e])
        out.extend((chrom, s, e) for s, e in merged)
    return out


def ref_weights(fl, means, sds, min_frag_p=0.001):
    p = [norm.pdf(fl, m, s) for m, s in zip(means, sds)]
    if all(x < min_frag_p for x in p):
        return [0.0, 0.0, 0.0, 0.0]
    return [x / sum(p) for x in p]


def ref_digested(frags, means, sds, min_frag_p=0.001):
    """Per-chromosome, per-track (short, mono, di, tri) weighted pileups."""
    out = {}
    lens = sorted(set(r - l for _, l, r in frags))
    w = {fl: ref_weights(fl, means, sds, min_frag_p) for fl in lens}
    for chrom in sorted(set(c for c, _, _ in frags)):
        fr = [(l, r) for c, l, r in frags if c == chrom]
        end = max(r for _, r in fr)
        tracks = np.zeros((4, end + 1))
        for l, r in fr:
            for k in range(4):
                tracks[k, l] += w[r - l][k]
                tracks[k, r] -= w[r - l][k]
        out[chrom] = np.cumsum(tracks, axis=1)[:, :end]
    return out


def read_bdg(path, skip_track=True):
    """bedGraph rows as (chrom, start, end, value)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if skip_track and line.startswith("track"):
                continue
            c, s, e, v = line.rstrip("\n").split("\t")
            rows.append((c, int(s), int(e), float(v)))
    return rows


def bdg_to_arrays(rows):
    """Per-chromosome per-base arrays from contiguous bedGraph rows."""
    out = {}
    for c in sorted(set(r[0] for r in rows)):
        rr = [r for r in rows if r[0] == c]
        a = np.zeros(rr[-1][2])
        prev = 0
        for _, s, e, v in rr:
            assert s == prev, "bedGraph rows are not contiguous"
            a[s:e] = v
            prev = e
        out[c] = a
    return out


def read_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def narrowpeak_regions(path):
    return [(f[0], int(f[1]), int(f[2]))
            for f in (line.split("\t") for line in read_lines(path))]


def triplets_from_states(rows, minlen):
    """Open regions flanked by adjacent nuc regions (states.bed rows).

    A region qualifies when nuc-open-nuc are contiguous on one
    chromosome and the three together span more than ``minlen``.
    """
    out = []
    for a, b, c in zip(rows, rows[1:], rows[2:]):
        if (a[3], b[3], c[3]) == ("nuc", "open", "nuc") and \
                a[0] == b[0] == c[0] and a[2] == b[1] and b[2] == c[1] and \
                c[2] - a[1] > minlen:
            out.append((b[0], b[1], b[2]))
    return out


def train_chrom(name):
    """Chromosome of a training data row. The original writes the bytes
    repr (b'chr1'); this accepts both spellings."""
    return re.sub(r"^b'(.*)'$", r"\1", name)


def read_states(path):
    return [(f[0], int(f[1]), int(f[2]), f[3])
            for f in (line.split("\t") for line in read_lines(path))]


# ------------------------------------
# Helpers: running hmmratac
# ------------------------------------

class _Collect(logging.Handler):
    def __init__(self):
        super().__init__(level=1)
        self.records = []

    def emit(self, record):
        self.records.append(record)


def split_messages(records):
    """(LEVEL, message) pairs of MACS3's log records, like conftest's
    parse_log: the memory prefix is removed and multi-line messages are
    split, continuation lines getting level ''."""
    out = []
    for r in records:
        if not r.name.startswith("MACS3"):
            continue            # e.g. hmmlearn's convergence messages
        msg = re.sub(r"^\[\d+ MB\] ", "", r.getMessage())
        lines = msg.split("\n")
        out.append((r.levelname, lines[0]))
        out.extend(("", x) for x in lines[1:])
    return out


def run_inprocess(argparser, args, workdir):
    """Run `macs3 hmmratac <args>` in this process.

    Parses ``args`` with the real parser, creates --outdir as bin/macs3
    does, and calls hmmratac_cmd.run. The temporary file hmmratac uses
    goes under ``workdir``. Returns a namespace with returncode (the
    SystemExit code, 1 for an exception), messages, exc and the paths.
    """
    full = ["hmmratac"] + [str(a) for a in args]
    vlogger = logging.getLogger("MACS3.Utilities.OptValidator")
    old_level = vlogger.level
    old_argv = sys.argv
    old_tempdir = tempfile.tempdir
    tmpd = Path(workdir) / "TMPDIR"
    tmpd.mkdir(parents=True, exist_ok=True)
    handler = _Collect()
    root = logging.getLogger()
    root.addHandler(handler)
    # the macs3 command sets the root level to INFO (logging.basicConfig
    # in MACS3.Utilities.Logger, a no-op under pytest)
    old_root_level = root.level
    root.setLevel(logging.INFO)
    code, exc, tb = 0, None, ""
    ns = argparser.parse_args(full)
    try:
        sys.argv = ["macs3"] + full
        tempfile.tempdir = str(tmpd)
        if ns.outdir:
            os.makedirs(ns.outdir, exist_ok=True)
        hmmratac_run(ns)
    except SystemExit as e:
        code = 0 if e.code is None else e.code
    except Exception as e:      # the process would exit with status 1
        code, exc, tb = 1, e, traceback.format_exc()
    finally:
        sys.argv = old_argv
        tempfile.tempdir = old_tempdir
        vlogger.setLevel(old_level)
        root.setLevel(old_root_level)
        root.removeHandler(handler)
    outdir = Path(ns.outdir) if ns.outdir else Path.cwd()
    return SimpleNamespace(returncode=code, exc=exc, traceback=tb,
                           messages=split_messages(handler.records),
                           outdir=outdir, name=ns.name,
                           path=lambda suffix: outdir / (ns.name + suffix),
                           files=sorted(os.listdir(outdir)))


def msgs(result, level=None):
    return [m for lv, m in result.messages if level is None or lv == level]


def log_value(result, prefix):
    """The message starting with ``prefix`` (exactly one must exist)."""
    found = [m for m in msgs(result) if m.startswith(prefix)]
    assert len(found) == 1, (prefix, found)
    return found[0]


def em_params(result):
    """(means, stddevs) printed in the '#  The means and stddevs' block."""
    means = log_value(result, "#             means:").split()[2:]
    sds = log_value(result, "#           stddevs:").split()[2:]
    return [float(x) for x in means], [float(x) for x in sds]


OUTPUT_SUFFIXES = {
    "narrowpeak": "_accessible_regions.narrowPeak",
    "cutoff": "_cutoff_analysis.tsv",
    "short": "_digested_short.bdg",
    "mono": "_digested_mono.bdg",
    "di": "_digested_di.bdg",
    "tri": "_digested_tri.bdg",
    "states": "_states.bed",
    "open": "_open.bdg",
    "nuc": "_nuc.bdg",
    "bg": "_bg.bdg",
    "train_bed": "_training_regions.bed",
    "train_data": "_training_data.txt",
    "train_len": "_training_lengths.txt",
    "model": "_model.json",
}


def expected_files(name, keys):
    return sorted(name + OUTPUT_SUFFIXES[k] for k in keys)


# ------------------------------------
# Fixtures
# ------------------------------------

@pytest.fixture(scope="module")
def atac(tmp_path_factory):
    d = tmp_path_factory.mktemp("hmmratac_data")
    frags = make_atac_fragments()
    bedpe = write_rows(d / "atac.bedpe", frags)
    fc, mean = fold_change(frags)
    total_len = sum(r - l for _, l, r in frags)
    minlen = int(np.float32(total_len) / np.float32(len(frags)))
    return SimpleNamespace(dir=d, frags=frags, bedpe=bedpe, fc=fc, mean=mean,
                           minlen=minlen)


@pytest.fixture(scope="module")
def full(atac, macs3_argparser, tmp_path_factory):
    """A run on the synthetic BEDPE saving every optional output."""
    out = tmp_path_factory.mktemp("full")
    return run_inprocess(macs3_argparser,
                         ["-i", atac.bedpe, "-f", "BEDPE", "-n", "full",
                          "--outdir", out, "--save-digested", "--save-states",
                          "--save-likelihoods", "--save-training-data"]
                         + TUNE, tmp_path_factory.mktemp("full_work"))


@pytest.fixture(scope="module")
def model_file(full):
    assert full.returncode == 0, full.traceback
    return str(full.path("_model.json"))


@pytest.fixture(scope="module")
def base(atac, model_file, macs3_argparser, tmp_path_factory):
    """The synthetic BEDPE decoded with the model of the full run."""
    out = tmp_path_factory.mktemp("base")
    return run_inprocess(macs3_argparser,
                         ["-i", atac.bedpe, "-f", "BEDPE", "-n", "base",
                          "--outdir", out, "--model", model_file] + TUNE,
                         tmp_path_factory.mktemp("base_work"))


@pytest.fixture(scope="module")
def base_noem(atac, model_file, macs3_argparser, tmp_path_factory):
    """As `base`, with EM skipped (default means and stddevs)."""
    out = tmp_path_factory.mktemp("base_noem")
    return run_inprocess(macs3_argparser,
                         ["-i", atac.bedpe, "-f", "BEDPE", "-n", "base",
                          "--outdir", out, "--model", model_file,
                          "--no-fragem"] + TUNE,
                         tmp_path_factory.mktemp("base_noem_work"))


@pytest.fixture
def quick(macs3_argparser, tmp_path):
    """``quick(args)`` runs hmmratac in-process with --outdir tmp_path/out."""
    def _q(args, name="t"):
        out = tmp_path / "out"
        return run_inprocess(macs3_argparser,
                             ["-n", name, "--outdir", out] + list(args),
                             tmp_path)
    return _q


# ------------------------------------
# save_proba_to_bedGraph
# ------------------------------------

PROBA = (b"chr1,30,0.700000,0.200000,0.100000\n"
         b"chr1,40,0.100000,0.300000,0.600000\n"
         b"chr1,70,0.250000,0.250000,0.500000\n"
         b"chr2,10,0.000000,0.500000,0.500000\n")


def test_save_proba_to_bedGraph_hand_example(tmp_path):
    # columns after chrom and bin end are states 0, 1, 2; here state 2
    # is open, 1 nuc and 0 bg. The region before the first bin of a
    # chromosome and gaps between bins are background (bg=1, others 0).
    f_open, f_nuc, f_bg = (str(tmp_path / x) for x in ("o.bdg", "n.bdg",
                                                       "b.bdg"))
    save_proba_to_bedGraph(io.BytesIO(PROBA), 10, f_open, f_nuc, f_bg,
                           2, 1, 0)
    assert read_lines(f_open) == ["chr1\t0\t20\t0.00000",
                                  "chr1\t20\t30\t0.10000",
                                  "chr1\t30\t40\t0.60000",
                                  "chr1\t40\t60\t0.00000",
                                  "chr1\t60\t70\t0.50000",
                                  "chr2\t0\t10\t0.50000"]
    assert read_lines(f_nuc) == ["chr1\t0\t20\t0.00000",
                                 "chr1\t20\t30\t0.20000",
                                 "chr1\t30\t40\t0.30000",
                                 "chr1\t40\t60\t0.00000",
                                 "chr1\t60\t70\t0.25000",
                                 "chr2\t0\t10\t0.50000"]
    assert read_lines(f_bg) == ["chr1\t0\t20\t1.00000",
                                "chr1\t20\t30\t0.70000",
                                "chr1\t30\t40\t0.10000",
                                "chr1\t40\t60\t1.00000",
                                "chr1\t60\t70\t0.25000",
                                "chr2\t0\t10\t0.00000"]


def test_save_proba_to_bedGraph_merges_equal_neighbours(tmp_path):
    proba = (b"chr1,10,0.500000,0.250000,0.250000\n"
             b"chr1,20,0.500000,0.250000,0.250000\n"
             b"chr1,30,0.000000,0.000000,1.000000\n"
             b"chr1,50,0.000000,0.100000,0.900000\n")
    f = [str(tmp_path / x) for x in ("o.bdg", "n.bdg", "b.bdg")]
    # state 0 is open, 1 nuc, 2 bg
    save_proba_to_bedGraph(io.BytesIO(proba), 10, f[0], f[1], f[2], 0, 1, 2)
    # open: 0.5 on [0,20), 0 on [20,30), the gap [30,40) and [40,50)
    assert read_lines(f[0]) == ["chr1\t0\t20\t0.50000",
                                "chr1\t20\t50\t0.00000"]
    assert read_lines(f[1]) == ["chr1\t0\t20\t0.25000",
                                "chr1\t20\t40\t0.00000",
                                "chr1\t40\t50\t0.10000"]
    assert read_lines(f[2]) == ["chr1\t0\t20\t0.25000",
                                "chr1\t20\t40\t1.00000",
                                "chr1\t40\t50\t0.90000"]


@pytest.mark.parametrize("binsize", [5, 10, 25])
def test_save_proba_to_bedGraph_binsize(tmp_path, binsize):
    proba = b"chrA,100,0.1,0.2,0.7\nchrA,%d,0.3,0.3,0.4\n" % (100 + binsize)
    f = [str(tmp_path / x) for x in ("o.bdg", "n.bdg", "b.bdg")]
    save_proba_to_bedGraph(io.BytesIO(proba), binsize, f[0], f[1], f[2],
                           0, 1, 2)
    s = 100 - binsize
    assert read_lines(f[2]) == ["chrA\t0\t%d\t1.00000" % s,
                                "chrA\t%d\t100\t0.70000" % s,
                                "chrA\t100\t%d\t0.40000" % (100 + binsize)]


def test_save_proba_to_bedGraph_empty(tmp_path):
    f = [str(tmp_path / x) for x in ("o.bdg", "n.bdg", "b.bdg")]
    save_proba_to_bedGraph(io.BytesIO(b""), 10, f[0], f[1], f[2], 0, 1, 2)
    for x in f:
        assert read_lines(x) == []


# ------------------------------------
# generate_states_path
# ------------------------------------

# probabilities per bin as (open, nuc, bg)
STATE_BINS = [(b"chr1", 20, (0.1, 0.8, 0.1)),
              (b"chr1", 30, (0.1, 0.7, 0.2)),
              (b"chr1", 40, (0.9, 0.05, 0.05)),
              (b"chr1", 50, (0.2, 0.6, 0.2)),
              (b"chr1", 80, (0.1, 0.1, 0.8)),
              (b"chr1", 90, (0.5, 0.5, 0.0)),
              (b"chr2", 10, (0.2, 0.2, 0.6))]

# by hand: chr1 starts at 10, so [0,10) is bg; two nuc bins merge; the
# gap [50,70) is bg and the bg bin [70,80) extends it; the tie at
# [80,90) goes to open (open, nuc, bg are tried in this order); chr2
# starts at 0 so no leading bg is added
STATE_PATH = [(b"chr1", 0, 10, "bg"), (b"chr1", 10, 30, "nuc"),
              (b"chr1", 30, 40, "open"), (b"chr1", 40, 50, "nuc"),
              (b"chr1", 50, 80, "bg"), (b"chr1", 80, 90, "open"),
              (b"chr2", 0, 10, "bg")]


def proba_file(bins, i_open, i_nuc, i_bg):
    lines = []
    for chrom, end, (po, pn, pb) in bins:
        cols = [0.0, 0.0, 0.0]
        cols[i_open], cols[i_nuc], cols[i_bg] = po, pn, pb
        lines.append(b"%s,%d,%f,%f,%f\n" % (chrom, end, cols[0], cols[1],
                                           cols[2]))
    return io.BytesIO(b"".join(lines))


@pytest.mark.parametrize("i_open, i_nuc, i_bg", [
    (0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)])
def test_generate_states_path_hand_example(i_open, i_nuc, i_bg):
    path = generate_states_path(proba_file(STATE_BINS, i_open, i_nuc, i_bg),
                                10, i_open, i_nuc, i_bg)
    assert path == STATE_PATH


def test_generate_states_path_tie_nuc_bg_goes_to_nuc():
    bins = [(b"chr1", 10, (0.0, 0.5, 0.5))]
    assert generate_states_path(proba_file(bins, 0, 1, 2), 10, 0, 1, 2) == \
        [(b"chr1", 0, 10, "nuc")]


def test_generate_states_path_binsize():
    bins = [(b"chr1", 100, (0.9, 0.05, 0.05)),
            (b"chr1", 125, (0.9, 0.05, 0.05))]
    assert generate_states_path(proba_file(bins, 0, 1, 2), 25, 0, 1, 2) == \
        [(b"chr1", 0, 75, "bg"), (b"chr1", 75, 125, "open")]


def collapse(path):
    """Merge adjacent segments with the same label."""
    out = []
    for c, s, e, lab in path:
        if out and out[-1][0] == c and out[-1][2] == s and out[-1][3] == lab:
            out[-1] = (c, out[-1][1], e, lab)
        else:
            out.append((c, s, e, lab))
    return out


def test_generate_states_path_gaps_are_background():
    # bins [0,10) open, [30,40) bg, [50,60) bg; the gaps [10,30) and
    # [40,50) are background
    bins = [(b"chr1", 10, (0.9, 0.05, 0.05)),
            (b"chr1", 40, (0.0, 0.1, 0.9)),
            (b"chr1", 60, (0.0, 0.1, 0.9))]
    path = generate_states_path(proba_file(bins, 0, 1, 2), 10, 0, 1, 2)
    assert collapse(path) == [(b"chr1", 0, 10, "open"),
                              (b"chr1", 10, 60, "bg")]
    # segments tile the chromosome without overlap
    for a, b in zip(path, path[1:]):
        assert a[2] == b[1]


def test_generate_states_path_empty():
    assert generate_states_path(io.BytesIO(b""), 10, 0, 1, 2) == []


def test_generate_states_path_single_bin_at_zero():
    bins = [(b"chrM", 10, (0.1, 0.1, 0.8))]
    assert generate_states_path(proba_file(bins, 2, 1, 0), 10, 2, 1, 0) == \
        [(b"chrM", 0, 10, "bg")]


# ------------------------------------
# save_states_bed
# ------------------------------------

def test_save_states_bed_skips_background():
    fh = io.StringIO()
    save_states_bed(STATE_PATH, fh)
    assert fh.getvalue() == ("chr1\t10\t30\tnuc\n"
                             "chr1\t30\t40\topen\n"
                             "chr1\t40\t50\tnuc\n"
                             "chr1\t80\t90\topen\n")


def test_save_states_bed_empty_and_all_background():
    fh = io.StringIO()
    save_states_bed([], fh)
    save_states_bed([(b"chr1", 0, 100, "bg"), (b"chr2", 0, 5, "bg")], fh)
    assert fh.getvalue() == ""


def test_save_states_bed_large_coordinates():
    fh = io.StringIO()
    save_states_bed([(b"chr1", 2**31 - 11, 2**31 - 1, "open")], fh)
    assert fh.getvalue() == "chr1\t2147483637\t2147483647\topen\n"


# ------------------------------------
# save_accessible_regions
# ------------------------------------

def score_track(blocks=((0, 250, 2.0), (250, 280, 5.5), (280, 320, 1.0),
                       (320, 560, 1.5), (560, 700, 1.0), (700, 1000, 0.5))):
    """Fold-change track on chr1 from (start, end, value) blocks."""
    bdg = bedGraphTrackI()
    for s, e, v in blocks:
        bdg.add_loc(b"chr1", s, e, v)
    return bdg


ACC_PATH = [(b"chr1", 0, 100, "bg"), (b"chr1", 100, 200, "nuc"),
            (b"chr1", 200, 300, "open"), (b"chr1", 300, 400, "nuc"),
            (b"chr1", 400, 500, "bg"), (b"chr1", 500, 550, "nuc"),
            (b"chr1", 550, 600, "open"), (b"chr1", 600, 650, "nuc"),
            (b"chr1", 650, 700, "bg"), (b"chr1", 700, 720, "nuc"),
            (b"chr1", 720, 730, "open"), (b"chr1", 730, 740, "nuc")]

# nuc-open-nuc at 100-400 (span 300) and 500-650 (span 150) qualify;
# 700-740 spans only 40. Region 1 [200,300): blocks 2.0, 5.5, 1.0; the
# summit is the middle of the 5.5 block, 265 (offset 65), score
# int(10 * 5.5) = 55. Region 2 [550,600): blocks [550,560)=1.5 and
# [560,600)=1.0; summit 555 (offset 5), score 15.
ACC_EXPECTED = ["chr1\t200\t300\tMACS_peak_1\t55\t.\t0\t0\t0\t65",
                "chr1\t550\t600\tMACS_peak_2\t15\t.\t0\t0\t0\t5"]


def accessible(path, minlen=100, bdg=None):
    fh = io.StringIO()
    save_accessible_regions(path, fh, minlen, bdg or score_track())
    return fh.getvalue().splitlines()


def test_save_accessible_regions_first_region_exact():
    assert accessible(ACC_PATH)[0] == ACC_EXPECTED[0]


def test_save_accessible_regions_reports_only_qualifying_regions():
    lines = accessible(ACC_PATH)
    regions = [tuple(x.split("\t")[1:3]) for x in lines]
    assert set(regions) <= {("200", "300"), ("550", "600")}


@pytest.mark.parametrize("minlen, reported", [
    (299, True), (300, False), (0, True)])
def test_save_accessible_regions_minlen_is_strict(minlen, reported):
    # the first triplet spans exactly 300; a longer second triplet
    # follows so that the first one is not the last region
    path = [(b"chr1", 100, 200, "nuc"), (b"chr1", 200, 300, "open"),
            (b"chr1", 300, 400, "nuc"), (b"chr1", 400, 500, "bg"),
            (b"chr1", 500, 600, "nuc"), (b"chr1", 600, 700, "open"),
            (b"chr1", 700, 900, "nuc")]
    lines = accessible(path, minlen=minlen)
    assert (ACC_EXPECTED[0] in lines) == reported


def test_save_accessible_regions_shared_nucleosome():
    # nuc-open-nuc-open-nuc: the middle nuc closes the first region and
    # opens the second; both open regions are reported
    path = [(b"chr1", 100, 200, "nuc"), (b"chr1", 200, 300, "open"),
            (b"chr1", 300, 400, "nuc"), (b"chr1", 400, 500, "open"),
            (b"chr1", 500, 600, "nuc"), (b"chr1", 600, 700, "bg"),
            (b"chr1", 700, 800, "nuc"), (b"chr1", 800, 850, "open"),
            (b"chr1", 850, 950, "nuc")]
    lines = accessible(path)
    # region 2 [400,500) lies in the 1.5 block: summit 450, score 15
    assert lines[:2] == ["chr1\t200\t300\tMACS_peak_1\t55\t.\t0\t0\t0\t65",
                         "chr1\t400\t500\tMACS_peak_2\t15\t.\t0\t0\t0\t50"]


# a qualifying region after the one under test, so that a wrongly
# accepted region would not be the (dropped) last one
TRAIL = [(b"chr1", 450, 500, "bg")] + ACC_PATH[5:8]


def test_save_accessible_regions_requires_contiguous_triplet():
    # a gap between nuc and open breaks the pattern
    path = [(b"chr1", 100, 200, "nuc"), (b"chr1", 210, 300, "open"),
            (b"chr1", 300, 400, "nuc")] + TRAIL
    lines = accessible(path)
    assert not any(x.startswith("chr1\t210\t300") for x in lines)


@pytest.mark.parametrize("pattern", [
    ("open", "open", "nuc"), ("nuc", "open", "bg"), ("bg", "open", "nuc"),
    ("nuc", "nuc", "nuc"), ("nuc", "bg", "nuc")])
def test_save_accessible_regions_other_patterns_not_reported(pattern):
    path = [(b"chr1", 100, 200, pattern[0]), (b"chr1", 200, 300, pattern[1]),
            (b"chr1", 300, 400, pattern[2])] + TRAIL
    lines = accessible(path)
    assert not any(x.startswith("chr1\t200\t300") for x in lines)


def test_save_accessible_regions_trailing_region_control():
    # control for the two tests above: the same layout with nuc-open-nuc
    # reports the first region
    path = [(b"chr1", 100, 200, "nuc"), (b"chr1", 200, 300, "open"),
            (b"chr1", 300, 400, "nuc")] + TRAIL
    assert accessible(path)[0] == ACC_EXPECTED[0]


def test_save_accessible_regions_summit_ties_take_middle():
    # three equal maxima in the first region: the summit is the middle
    # one (index (3+1)//2 - 1 = 1 of the tied blocks)
    bdg = bedGraphTrackI()
    for s, e, v in [(0, 210, 1.0), (210, 220, 3.0), (220, 240, 1.0),
                    (240, 250, 3.0), (250, 270, 1.0), (270, 290, 3.0),
                    (290, 1000, 1.0)]:
        bdg.add_loc(b"chr1", s, e, v)
    lines = accessible(ACC_PATH, bdg=bdg)
    assert lines[0] == "chr1\t200\t300\tMACS_peak_1\t30\t.\t0\t0\t0\t45"


def test_save_accessible_regions_empty_inputs():
    assert accessible([]) == []
    assert accessible(ACC_PATH[:2]) == []
    assert accessible([(b"chr1", 0, 1000, "bg")]) == []


# ------------------------------------
# macs3 hmmratac: the full synthetic run
# ------------------------------------

def test_full_run_succeeds_and_creates_all_outputs(full):
    assert full.returncode == 0, full.traceback
    assert full.files == sorted(expected_files(
        "full", OUTPUT_SUFFIXES.keys()))


def test_full_run_log_milestones(full, atac):
    m = msgs(full, "INFO")
    total = len(atac.frags)
    for line in ["#1 Read fragments from BEDPE file...",
                 "#  Read %d fragments." % total,
                 "#2 Use EM algorithm to estimate means and stddevs of "
                 "fragment lengths",
                 "#  for mono-, di-, and tri-nucleosomal signals...",
                 "#  The means and stddevs after EM:",
                 "#  Compute the weights for each fragment length for each "
                 "of the four signal types",
                 "#  Generate short, mono-, di-, and tri-nucleosomal signals",
                 "#  Pile up all fragments",
                 "#  Convert pileup to fold-change over average signal",
                 "#3 Look for training set from %d fragments" % total,
                 "#  Call peak above within fold-change range of 4 and 20.",
                 "#   The minimum length of the region is set as the average "
                 "template/fragment length in the dataset: %d" % atac.minlen,
                 "#   The maximum gap to merge nearby significant regions is "
                 "set as the flanking size to extend training regions: 1000",
                 "#  We expand the training regions with 1000 basepairs and "
                 "merge overlap",
                 "#  Training regions have been saved to "
                 "`full_training_regions.bed` ",
                 "#4 Train Hidden Markov Model with Multivariate Gaussian "
                 "Emission",
                 "#  Extract signals in training regions with bin size of 10",
                 "#  Use Baum-Welch algorithm to train the HMM",
                 "#  Write HMM parameters into JSON: %s"
                 % full.path("_model.json"),
                 "#  The Hidden Markov Model for signals of binsize of 10 "
                 "basepairs:",
                 "#   HMM Emissions (means): ",
                 "#5 Decode with Viterbi to predict states",
                 "# Write the likelihoods for each states into three bedGraph "
                 "files full_open.bdg, full_nuc.bdg, and full_bg.bdg",
                 "# Write states assignments in a BED file: full_states.bed",
                 "# Write accessible regions in a narrowPeak file: "
                 "full_accessible_regions.narrowPeak"]:
        assert line in m, line
    assert m[-1] == "# Finished"
    assert msgs(full, "WARNING") == [] and msgs(full, "ERROR") == []


def test_full_run_training_counts_match_reference(full, atac):
    peaks = ref_call_peaks(atac.fc, 4, atac.minlen, 1000)
    kept = [p for p in peaks if 4 <= p[3] < 20]
    assert log_value(full, "#  Total training regions called after applying "
                     "the lower cutoff 4:") == \
        "#  Total training regions called after applying the lower " \
        "cutoff 4: %d" % len(peaks)
    assert log_value(full, "#  Total training regions after filtering with "
                     "upper cutoff 20:") == \
        "#  Total training regions after filtering with upper cutoff " \
        "20: %d" % len(kept)
    assert len(kept) >= 5


def test_full_run_training_regions_bed_matches_reference(full, atac):
    peaks = [p for p in ref_call_peaks(atac.fc, 4, atac.minlen, 1000)
             if 4 <= p[3] < 20]
    expected = ["%s\t%d\t%d" % (c.decode(), s, e)
                for c, s, e in ref_training_regions(peaks, 1000)]
    assert read_lines(full.path("_training_regions.bed")) == expected


def test_full_run_candidate_peaks_match_reference(full, atac):
    n = len(ref_call_peaks(atac.fc, 2, atac.minlen, 1000))
    assert "#5  Total candidate peaks : %d" % n in msgs(full)


def test_full_run_training_data_consistent(full):
    data = [line.split("\t") for line in
            read_lines(full.path("_training_data.txt"))]
    lengths = [int(x) for x in read_lines(full.path("_training_lengths.txt"))]
    regions = [line.split("\t") for line in
               read_lines(full.path("_training_regions.bed"))]
    assert sum(lengths) == len(data)
    assert len(lengths) == len(regions)
    assert all(len(row) == 6 for row in data)
    for c, pos, *vals in data:
        c = train_chrom(c)
        assert int(pos) % 10 == 0
        assert all(float(v) >= 0.0001 for v in vals)
        assert any(c == rc and int(rs) // 10 * 10 < int(pos) <= int(re_)
                   for rc, rs, re_ in regions)


def test_full_run_model_json_structure(full):
    with open(full.path("_model.json")) as fh:
        m = json.load(fh)
    assert list(m.keys()) == ["startprob", "transmat", "means", "covars",
                              "covariance_type", "n_features",
                              "i_open_region", "i_background_region",
                              "i_nucleosomal_region", "hmm_binsize",
                              "hmm_type"]
    assert np.array(m["startprob"]).shape == (3,)
    assert np.array(m["transmat"]).shape == (3, 3)
    assert np.array(m["means"]).shape == (3, 4)
    assert np.array(m["covars"]).shape == (3, 4, 4)
    assert m["covariance_type"] == "full"
    assert m["n_features"] == 4
    assert m["hmm_binsize"] == 10
    assert m["hmm_type"] == "gaussian"
    assert sorted([m["i_open_region"], m["i_background_region"],
                   m["i_nucleosomal_region"]]) == [0, 1, 2]
    sums = np.array(m["means"]).sum(axis=1)
    # open = largest total emission, background = smallest
    assert m["i_open_region"] == int(np.argmax(sums))
    assert m["i_background_region"] == int(np.argmin(sums))


def test_full_run_state_assignment_log(full):
    with open(full.path("_model.json")) as fh:
        m = json.load(fh)
    assert "#   open state index: state%d" % m["i_open_region"] in msgs(full)
    assert "#   nucleosomal state index: state%d" % \
        m["i_nucleosomal_region"] in msgs(full)
    assert "#   background state index: state%d" % \
        m["i_background_region"] in msgs(full)


def test_full_run_narrowpeak_format(full):
    lines = read_lines(full.path("_accessible_regions.narrowPeak"))
    assert len(lines) >= 5
    prev = None
    for i, line in enumerate(lines):
        f = line.split("\t")
        assert len(f) == 10
        assert f[3] == "MACS_peak_%d" % (i + 1)
        assert f[5] == "." and f[6:9] == ["0", "0", "0"]
        s, e, summit = int(f[1]), int(f[2]), int(f[9])
        assert 0 <= s < e and 0 <= summit < e - s
        assert s % 10 == 0 and e % 10 == 0
        if prev and prev[0] == f[0]:
            assert s >= prev[1]
        prev = (f[0], e)


def test_full_run_narrowpeak_scores_and_summits(full, atac):
    # score: int(10 * max fold change in the region); summit: middle of
    # the highest block (the middle one of tied blocks)
    for line in read_lines(full.path("_accessible_regions.narrowPeak")):
        f = line.split("\t")
        chrom, s, e = f[0].encode(), int(f[1]), int(f[2])
        a = atac.fc[chrom][s:e]
        assert int(f[4]) == int(10 * float(a.max()))
        top = np.concatenate([[False], a == a.max(), [False]])
        d = np.diff(top.astype(np.int8))
        blocks = list(zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]))
        bs, be = blocks[(len(blocks) + 1) // 2 - 1]
        assert int(f[9]) == int((s + bs + s + be) / 2) - s


def test_full_run_finds_designed_open_regions(full):
    regions = narrowpeak_regions(full.path("_accessible_regions.narrowPeak"))
    found = 0
    for chrom, (_, centers) in LAYOUT.items():
        for c in centers:
            if any(rc == chrom.decode() and rs - 300 <= c <= re_ + 300
                   for rc, rs, re_ in regions):
                found += 1
    assert found >= 7


def test_full_run_accessible_regions_are_states_triplets(full):
    states = read_states(full.path("_states.bed"))
    assert set(r[3] for r in states) <= {"open", "nuc"}
    qualifying = triplets_from_states(states, 100)
    regions = narrowpeak_regions(full.path("_accessible_regions.narrowPeak"))
    assert set(regions) <= set(qualifying)
    assert len(regions) >= len(qualifying) - 1


def test_full_run_likelihoods_sum_to_one(full):
    arrays = [bdg_to_arrays(read_bdg(full.path(s))) for s in
              ("_open.bdg", "_nuc.bdg", "_bg.bdg")]
    for x in ("_open.bdg", "_nuc.bdg", "_bg.bdg"):
        assert not read_lines(full.path(x))[0].startswith("track")
    for chrom in arrays[0]:
        n = len(arrays[0][chrom])
        assert len(arrays[1][chrom]) == n and len(arrays[2][chrom]) == n
        total = arrays[0][chrom] + arrays[1][chrom] + arrays[2][chrom]
        assert total == pytest.approx(np.ones(n), abs=2e-5)


def test_full_run_likelihoods_agree_with_states(full):
    # the state of every bin is the one with the largest likelihood
    arrays = {k: bdg_to_arrays(read_bdg(full.path("_%s.bdg" % k)))
              for k in ("open", "nuc", "bg")}
    for c, s, e, label in read_states(full.path("_states.bed")):
        mid = (s + e) // 2
        v = {k: arrays[k][c][mid] for k in arrays}
        assert v[label] >= max(v.values()) - 1e-5


def test_full_run_digested_track_lines(full):
    for k in ("short", "mono", "di", "tri"):
        first = read_lines(full.path("_digested_%s.bdg" % k))[0]
        assert first == ('track type=bedGraph name="%s" description="%s" '
                         'visibility=2 alwaysZero=on' % (k, k))


def test_full_run_digested_signals_match_reference(full, atac):
    means, sds = em_params(full)
    ref = ref_digested(atac.frags, means, sds)
    for k, key in enumerate(("short", "mono", "di", "tri")):
        got = bdg_to_arrays(read_bdg(full.path("_digested_%s.bdg" % key)))
        assert sorted(got) == ["chr1", "chr2"]
        for chrom, a in got.items():
            r = ref[chrom.encode()][k]
            assert len(a) == len(r)
            assert a == pytest.approx(r, abs=5e-3)


def test_full_run_em_means_near_synthetic_modes(full):
    means, sds = em_params(full)
    # the short mean is not fitted by EM; mono/di/tri come from EM
    assert means[0] == 50
    assert sds[0] == 20
    assert 160 < means[1] < 230
    assert 340 < means[2] < 460
    assert 520 < means[3] < 680


def test_full_run_cutoff_analysis_matches_reference(full, atac):
    assert read_lines(full.path("_cutoff_analysis.tsv")) == \
        ref_cutoff_analysis(atac.fc, atac.minlen, 1000)


def test_full_run_cutoff_analysis_report(full):
    lines = read_lines(full.path("_cutoff_analysis.tsv"))
    assert lines[0] == "score\tnpeaks\tlpeaks\tavelpeak"
    rows = [line.split("\t") for line in lines[1:]]
    assert 1 <= len(rows) <= 100
    scores = [float(r[0]) for r in rows]
    assert scores == sorted(scores, reverse=True)
    lpeaks = [int(r[2]) for r in rows]
    # a lower cutoff never gives less total peak length
    assert lpeaks == sorted(lpeaks)
    for r in rows:
        assert int(r[1]) > 0
        assert r[3] == "%.2f" % (int(r[2]) / int(r[1]))


# ------------------------------------
# macs3 hmmratac: options, run in this process
# ------------------------------------

def test_model_option_reproduces_training_run(full, base):
    # -- model: the saved model gives the same accessible regions as the
    # run that trained it
    assert base.returncode == 0, base.traceback
    assert read_lines(base.path("_accessible_regions.narrowPeak")) == \
        read_lines(full.path("_accessible_regions.narrowPeak"))


def test_model_option_skips_training(base, model_file):
    assert base.files == ["base_accessible_regions.narrowPeak"]
    m = msgs(base)
    assert "#3 Skip this step of looking for training set since a Hidden " \
        "Markov Model file has been provided!" in m
    assert "#4 Load Hidden Markov Model from given model file" in m
    assert not any(x.startswith("#4 Train") for x in m)


def test_name_default_is_NA(macs3_argparser, atac, tmp_path):
    out = tmp_path / "o"
    r = run_inprocess(macs3_argparser,
                      ["-i", atac.bedpe, "-f", "BEDPE", "--outdir", out,
                       "--cutoff-analysis-only"], tmp_path)
    assert r.files == ["NA_cutoff_analysis.tsv"]


def test_name_option_prefixes_files(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only"],
              name="my.sample")
    assert r.files == ["my.sample_cutoff_analysis.tsv"]


def test_input_multiple_files_are_pooled(quick, atac, base, model_file,
                                         tmp_path):
    a = write_rows(tmp_path / "a.bedpe", atac.frags[0::2])
    b = write_rows(tmp_path / "b.bedpe", atac.frags[1::2])
    r = quick(["-i", a, b, "-f", "BEDPE", "--model", model_file] + TUNE)
    assert r.returncode == 0, r.traceback
    assert "#  Read %d fragments." % len(atac.frags) in msgs(r)
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


def test_format_bampe_same_as_bedpe(quick, atac, base, model_file,
                                    make_alignments, make_pe_pair):
    reads = []
    for i, (c, l, r) in enumerate(atac.frags):
        reads.extend(make_pe_pair("f%d" % i, c.decode(), l, r,
                                  readlen=min(36, r - l)))
    bam = make_alignments(reads, refs=(("chr1", 60000), ("chr2", 30000)),
                          name="atac.bam")
    r = quick(["-i", bam, "--model", model_file] + TUNE)
    assert r.returncode == 0, r.traceback
    assert "#1 Read fragments from BAMPE file..." in msgs(r)
    assert "#  Read %d fragments." % len(atac.frags) in msgs(r)
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


def frag_rows(frags, counts=None, barcodes=4):
    return [(c, l, r, "BC%d" % (i % barcodes),
             1 if counts is None else counts[i])
            for i, (c, l, r) in enumerate(frags)]


def test_format_frag_same_as_bedpe_without_em(quick, atac, base_noem,
                                              model_file, tmp_path):
    frag = write_rows(tmp_path / "atac.frag.tsv", frag_rows(atac.frags))
    r = quick(["-i", frag, "-f", "FRAG", "--model", model_file,
               "--no-fragem"] + TUNE)
    assert r.returncode == 0, r.traceback
    assert "#1 Read fragments from FRAG file..." in msgs(r)
    assert "#  Read %d fragments." % len(atac.frags) in msgs(r)
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base_noem.path("_accessible_regions.narrowPeak"))


def test_format_frag_with_em_negative_variance_exits(quick, atac, tmp_path,
                                                    model_file, capsys):
    """Upstream 248fd6d added --jump with default 0.5 (HMMR_EM's default
    was 1.5), with which the variance update cannot go negative; this
    test now passes --jump 1.5. The expectation is unchanged."""
    # The fragments of the BEDPE runs, as FRAG with EM. PETrackII draws
    # its own 10% EM subsample (seed 10151), which holds only two
    # tri-nucleosome lengths, 584 and 607 (population variance 132.25),
    # so HMMR_EM's over-relaxed update 400 + 1.5 * (132.25 - 400) is
    # negative (PETrackI's subsample has seven lengths above 500 and EM
    # succeeds). HMMR_EM prints its advice to change --means and
    # --stddevs, and the next density evaluation exits with status 1
    # before any output is written; test_HMMR_EM.py covers the same path
    # on small inputs.
    frag = write_rows(tmp_path / "atac.frag.tsv", frag_rows(atac.frags))
    r = quick(["-i", frag, "-f", "FRAG", "--model", model_file,
               "--jump", "1.5"] + TUNE)
    assert r.returncode == 1 and r.exc is None
    assert capsys.readouterr().out == \
        " ValueError:  Adjust --means and --stddevs options and re-run " \
        "command\n"
    assert "# Downsampled 167 fragments will be used for EM training..." \
        in msgs(r)
    assert r.files == []


def test_barcodes_select_fragments(quick, atac, model_file, tmp_path):
    rows = frag_rows(atac.frags)
    frag = write_rows(tmp_path / "all.frag.tsv", rows)
    bc = tmp_path / "bc.txt"
    bc.write_text("BC0\nBC2\n")
    keep = [x for x in rows if x[3] in ("BC0", "BC2")]
    sub = write_rows(tmp_path / "sub.frag.tsv", keep)
    r1 = quick(["-i", frag, "-f", "FRAG", "--barcodes", bc, "--model",
                model_file, "--no-fragem"] + TUNE, name="bc")
    r2 = quick(["-i", sub, "-f", "FRAG", "--model", model_file,
                "--no-fragem"] + TUNE, name="sub")
    assert r1.returncode == 0 and r2.returncode == 0
    assert "#1 extract fragments with given barcodes" in msgs(r1)
    assert "#   extracted %d fragments in treatment" % len(keep) in msgs(r1)
    assert "#  Read %d fragments." % len(keep) in msgs(r1)
    assert read_lines(r1.path("_accessible_regions.narrowPeak")) == \
        read_lines(r2.path("_accessible_regions.narrowPeak"))


def test_barcodes_ignored_for_bedpe(quick, atac, base, model_file, tmp_path):
    bc = tmp_path / "bc.txt"
    bc.write_text("BC0\n")
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--barcodes", bc,
               "--model", model_file] + TUNE)
    assert "#1 extract fragments with given barcodes" not in msgs(r)
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


def test_barcodes_none_matching_fails(quick, atac, model_file, tmp_path):
    frag = write_rows(tmp_path / "all.frag.tsv", frag_rows(atac.frags))
    bc = tmp_path / "bc.txt"
    bc.write_text("NOSUCHBARCODE\n")
    r = quick(["-i", frag, "-f", "FRAG", "--barcodes", bc, "--model",
               model_file] + TUNE)
    assert r.returncode == 1
    assert isinstance(r.exc, AssertionError)
    assert "no fragments in PETrackII" in str(r.exc)


@pytest.mark.parametrize("max_count", [None, 0, 1, 2])
def test_max_count_caps_frag_counts(quick, atac, model_file, tmp_path,
                                    max_count):
    uniq = sorted(set(atac.frags))
    counts = [3 if i % 7 == 0 else 1 for i in range(len(uniq))]
    frag = write_rows(tmp_path / "c.frag.tsv", frag_rows(uniq, counts))
    args = ["-i", frag, "-f", "FRAG", "--model", model_file, "--no-fragem",
            "--cutoff-analysis-only"]
    if max_count is not None:
        args += ["--max-count", max_count]
    r = quick(args)
    cap = max_count or 10**9
    assert "#  Read %d fragments." % sum(min(c, cap) for c in counts) in \
        msgs(r)


def test_max_count_one_equals_unit_counts(quick, atac, model_file, tmp_path):
    uniq = sorted(set(atac.frags))
    counts = [3 if i % 7 == 0 else 1 for i in range(len(uniq))]
    f3 = write_rows(tmp_path / "c3.frag.tsv", frag_rows(uniq, counts))
    f1 = write_rows(tmp_path / "c1.frag.tsv", frag_rows(uniq))
    r3 = quick(["-i", f3, "-f", "FRAG", "--max-count", "1", "--model",
                model_file, "--no-fragem"] + TUNE, name="c3")
    r1 = quick(["-i", f1, "-f", "FRAG", "--model", model_file,
                "--no-fragem"] + TUNE, name="c1")
    assert r3.returncode == 0 and r1.returncode == 0
    assert read_lines(r3.path("_accessible_regions.narrowPeak")) == \
        read_lines(r1.path("_accessible_regions.narrowPeak"))


def test_max_count_ignored_for_bedpe(quick, atac, base, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--max-count", "1",
               "--model", model_file] + TUNE)
    assert "#  Read %d fragments." % len(atac.frags) in msgs(r)
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


def test_outdir_is_created(macs3_argparser, atac, tmp_path):
    out = tmp_path / "a" / "b"
    r = run_inprocess(macs3_argparser,
                      ["-i", atac.bedpe, "-f", "BEDPE", "-n", "x",
                       "--outdir", out, "--cutoff-analysis-only"], tmp_path)
    assert r.files == ["x_cutoff_analysis.tsv"]
    assert (out / "x_cutoff_analysis.tsv").is_file()


def test_cutoff_analysis_only_outputs(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only"])
    assert r.files == ["t_cutoff_analysis.tsv"]
    m = msgs(r)
    assert "#3 Generate cutoff analysis report from %d fragments" % \
        len(atac.frags) in m
    assert "#   Please review the cutoff analysis result in %s" % \
        r.path("_cutoff_analysis.tsv") in m
    assert not any(x.startswith("#4") for x in m)
    lines = read_lines(r.path("_cutoff_analysis.tsv"))
    assert lines[0] == "score\tnpeaks\tlpeaks\tavelpeak"
    assert 1 < len(lines) <= 101


def cutoff_rows(r):
    return [line.split("\t")
            for line in read_lines(r.path("_cutoff_analysis.tsv"))[1:]]


@pytest.mark.parametrize("cmax", [3, 5, 100])
def test_cutoff_analysis_max(quick, atac, cmax):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--cutoff-analysis-max", cmax])
    scores = [float(x[0]) for x in cutoff_rows(r)]
    assert scores and max(scores) <= cmax
    # with the maximum below the largest fold change, the cutoffs form
    # a grid from 0 with step cmax/100
    if cmax <= 5:
        step = cmax / 100
        for s in scores:
            assert abs(s - round(s / step) * step) <= 0.0056


@pytest.mark.parametrize("cmax, steps, flank", [(3, 10, 1000),
                                                (100, 100, 1000),
                                                (8, 37, 300)])
def test_cutoff_analysis_only_matches_reference(quick, atac, cmax, steps,
                                                flank):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--cutoff-analysis-max", cmax, "--cutoff-analysis-steps",
               steps, "--training-flanking", flank])
    assert read_lines(r.path("_cutoff_analysis.tsv")) == \
        ref_cutoff_analysis(atac.fc, atac.minlen, flank, cmax, steps)


@pytest.mark.parametrize("steps", [5, 10, 37])
def test_cutoff_analysis_steps(quick, atac, steps):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--cutoff-analysis-max", "3", "--cutoff-analysis-steps",
               steps])
    rows = cutoff_rows(r)
    # np.arange(0, max, max/steps) can give one extra cutoff equal to
    # max through floating point rounding
    assert 1 <= len(rows) <= steps + 1
    # cutoffs k * 3 / steps, rounded to 3 decimals, printed with 2
    step = 3 / steps
    for x in rows:
        s = float(x[0])
        assert abs(s - round(s / step) * step) <= 0.0056


def test_cutoff_analysis_steps_more_rows_with_more_steps(quick, atac):
    n = []
    for steps in (10, 50):
        r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
                   "--cutoff-analysis-max", "3", "--cutoff-analysis-steps",
                   steps], name="s%d" % steps)
        n.append(len(cutoff_rows(r)))
    assert n[0] < n[1]


def test_save_flags_absent_by_default(base):
    for k in ("short", "mono", "di", "tri", "states", "open", "nuc", "bg",
              "train_bed", "train_data", "train_len"):
        assert not base.path(OUTPUT_SUFFIXES[k]).exists()


def test_save_digested_with_cutoff_analysis_only(quick, atac):
    # the digested signals are written before the cutoff analysis exit
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--save-digested", "--no-fragem"])
    assert r.files == expected_files(
        "t", ["cutoff", "short", "mono", "di", "tri"])
    ref = ref_digested(atac.frags, [50, 200, 400, 600], [20, 20, 20, 20])
    got = bdg_to_arrays(read_bdg(r.path("_digested_short.bdg")))
    for chrom, a in got.items():
        assert a == pytest.approx(ref[chrom.encode()][0], abs=5e-3)


def test_save_states_only(quick, atac, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--save-states"] + TUNE)
    assert r.files == expected_files("t", ["narrowpeak", "states"])
    assert "# Write states assignments in a BED file: t_states.bed" in \
        msgs(r)


def test_save_likelihoods_only(quick, atac, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--save-likelihoods"] + TUNE)
    assert r.files == expected_files(
        "t", ["narrowpeak", "open", "nuc", "bg"])
    assert "# finished writing proba_to_bedgraph" in msgs(r)


def test_save_training_data_with_training_bed(quick, atac, tmp_path):
    # -t: training data and lengths are saved, the training regions BED
    # is not (the regions came from the user)
    bed = write_rows(tmp_path / "train.bed", [("chr1", 15000, 17000),
                                              ("chr2", 14500, 15500)])
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "-t", bed, "--modelonly",
               "--save-training-data"] + TUNE)
    assert r.returncode == 0, r.traceback
    assert r.files == expected_files(
        "t", ["train_data", "train_len", "model"])
    m = msgs(r)
    assert "#3 Read training regions from BED file: %s" % bed in m
    assert "#  Training regions have been read from bedfile" in m
    lengths = [int(x) for x in read_lines(r.path("_training_lengths.txt"))]
    # 200 bins of 10 bp in chr1:15000-17000 and 100 in chr2:14500-15500
    assert sorted(lengths) == [100, 200]
    data = [x.split("\t") for x in read_lines(r.path("_training_data.txt"))]
    assert len(data) == 300
    pos = sorted((train_chrom(c), int(p)) for c, p, *_ in data)
    assert pos == sorted([("chr1", p) for p in range(15010, 17001, 10)] +
                         [("chr2", p) for p in range(14510, 15501, 10)])


def test_no_fragem_uses_given_means(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--no-fragem", "--means", "40", "180", "380", "580",
               "--stddevs", "15", "25", "30", "35"])
    m = msgs(r)
    assert "#2 EM is skipped. The following means and stddevs will be used:"\
        in m
    i = m.index("#2 EM is skipped. The following means and stddevs will be "
                "used:")
    # labels right-aligned in 10 characters, numbers as %10.4g
    assert m[i + 1:i + 4] == [
        "#" + " " * 20 + "%10s %10s %10s %10s" % ("short", "mono", "di",
                                                 "tri"),
        "#             means: " + "%10s %10s %10s %10s" % ("40", "180", "380",
                                                         "580"),
        "#           stddevs: " + "%10s %10s %10s %10s" % ("15", "25", "30",
                                                         "35")]
    assert "# EM training not performed on fragment distribution. " in m


def test_no_fragem_default_means(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--no-fragem"])
    assert em_params(r) == ([50, 200, 400, 600], [20, 20, 20, 20])


def test_means_with_em_keep_short_parameters(quick, atac):
    # EM fits mono, di and tri only; the short mean and sd are kept
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--means", "45", "210", "410", "610", "--stddevs", "12", "20",
               "20", "20"])
    means, sds = em_params(r)
    assert means[0] == 45 and sds[0] == 12
    assert means[1:] != [210, 410, 610]


@pytest.mark.parametrize("p", [0.001, 0.01])
def test_min_frag_p_controls_digested_signals(quick, atac, p):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--save-digested", "--no-fragem", "--min-frag-p", p])
    ref = ref_digested(atac.frags, [50, 200, 400, 600], [20, 20, 20, 20],
                       min_frag_p=p)
    for k, key in enumerate(("short", "mono", "di", "tri")):
        got = bdg_to_arrays(read_bdg(r.path("_digested_%s.bdg" % key)))
        for chrom, a in got.items():
            assert a == pytest.approx(ref[chrom.encode()][k], abs=5e-3)


def test_binsize_option(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--modelonly",
               "--save-training-data", "--binsize", "20"] + TUNE)
    assert r.returncode == 0, r.traceback
    with open(r.path("_model.json")) as fh:
        assert json.load(fh)["hmm_binsize"] == 20
    assert "#  Extract signals in training regions with bin size of 20" in \
        msgs(r)
    for line in read_lines(r.path("_training_data.txt")):
        assert int(line.split("\t")[1]) % 20 == 0


def test_binsize_taken_from_model(quick, atac, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--binsize", "50", "--save-states"] + TUNE)
    assert "#  The Hidden Markov Model for signals of binsize of 10 " \
        "basepairs:" in msgs(r)
    # bins of 10 bp: some state boundaries are not multiples of 50
    states = read_states(r.path("_states.bed"))
    assert all(s % 10 == 0 and e % 10 == 0 for c, s, e, lab in states)
    assert any(s % 50 or e % 50 for c, s, e, lab in states)


@pytest.mark.parametrize("lower, upper", [(4, 20), (4, 8), (6, 20)])
def test_lower_upper_select_training_regions(quick, atac, lower, upper):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--modelonly",
               "--save-training-data", "-l", lower, "-u", upper, "-c", "2"])
    peaks = ref_call_peaks(atac.fc, lower, atac.minlen, 1000)
    kept = [p for p in peaks if lower <= p[3] < upper]
    m = msgs(r)
    assert "#  Call peak above within fold-change range of %d and %d." % (
        lower, upper) in m
    assert "#  Total training regions called after applying the lower " \
        "cutoff %d: %d" % (lower, len(peaks)) in m
    assert "#  Total training regions after filtering with upper cutoff " \
        "%d: %d" % (upper, len(kept)) in m
    if kept:
        assert read_lines(r.path("_training_regions.bed")) == [
            "%s\t%d\t%d" % (c.decode(), s, e)
            for c, s, e in ref_training_regions(kept, 1000)]


def test_no_training_regions_found(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "-l", "50", "-u", "60",
               "-c", "2"])
    assert r.returncode == 1
    assert ("CRITICAL", "# No training regions found. Please adjust the "
            "lower or upper cutoff.") in r.messages
    assert str(r.exc) == "Not enough training regions!"
    # the cutoff analysis was written before the failure
    assert r.files == ["t_cutoff_analysis.tsv"]


def test_maxtrain_limits_training_regions(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--modelonly",
               "--save-training-data", "--maxTrain", "3"] + TUNE)
    assert r.returncode == 0, r.traceback
    assert "#  We randomly pick 3 regions for training" in msgs(r)
    peaks = [p for p in ref_call_peaks(atac.fc, 4, atac.minlen, 1000)
             if 4 <= p[3] < 20]
    possible = set(ref_training_regions(peaks, 1000))
    got = [tuple(x.split("\t")) for x in
           read_lines(r.path("_training_regions.bed"))]
    assert 1 <= len(got) <= 3
    for c, s, e in got:
        assert (c.encode(), int(s), int(e)) in possible


@pytest.mark.parametrize("flank", [300, 2000])
def test_training_flanking(quick, atac, flank):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--modelonly",
               "--save-training-data", "--training-flanking", flank] + TUNE)
    assert r.returncode == 0, r.traceback
    peaks = [p for p in ref_call_peaks(atac.fc, 4, atac.minlen, flank)
             if 4 <= p[3] < 20]
    assert read_lines(r.path("_training_regions.bed")) == [
        "%s\t%d\t%d" % (c.decode(), s, e)
        for c, s, e in ref_training_regions(peaks, flank)]
    assert "#  We expand the training regions with %d basepairs and merge " \
        "overlap" % flank in msgs(r)


def test_modelonly(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--modelonly"] + TUNE)
    assert r.returncode == 0
    assert r.files == expected_files("t", ["cutoff", "model"])
    m = msgs(r)
    assert m[-1] == "#  Complete - HMM model was saved, program exited " \
        "(--modelonly option was provided) "
    assert "# Program will stop after generating model, which can be later " \
        "applied with '--model'. " in m


def test_modelonly_ignored_with_model(quick, atac, base, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--modelonly"] + TUNE)
    assert r.returncode == 0
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


def test_hmm_type_poisson_model(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--modelonly",
               "--hmm-type", "poisson", "--save-training-data"] + TUNE)
    assert r.returncode == 0, r.traceback
    with open(r.path("_model.json")) as fh:
        m = json.load(fh)
    assert list(m.keys()) == ["startprob", "transmat", "lambdas",
                              "n_features", "i_open_region",
                              "i_background_region", "i_nucleosomal_region",
                              "hmm_binsize", "hmm_type"]
    assert m["hmm_type"] == "poisson"
    assert np.array(m["lambdas"]).shape == (3, 4)
    sums = np.array(m["lambdas"]).sum(axis=1)
    assert m["i_open_region"] == int(np.argmax(sums))
    # poisson training data are integers
    for line in read_lines(r.path("_training_data.txt")):
        for v in line.split("\t")[2:]:
            assert re.fullmatch(r"\d+", v)


def test_hmm_type_poisson_full_run(quick, atac):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--hmm-type", "poisson"]
              + TUNE)
    assert r.returncode == 0, r.traceback
    assert "#   HMM Emissions (lambdas): " in msgs(r)
    assert "#   HMM Emissions (means): " not in msgs(r)
    assert r.path("_accessible_regions.narrowPeak").exists()


def test_hmm_type_taken_from_model(quick, atac, base, model_file):
    # the gaussian model file wins over --hmm-type poisson
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--hmm-type", "poisson"] + TUNE)
    assert "#   HMM Emissions (means): " in msgs(r)
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


@pytest.mark.parametrize("c", [1.5, 2, 3])
def test_prescan_cutoff_candidate_peaks(quick, atac, model_file, c):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "-l", "4", "-u", "20", "-c", c])
    n = len(ref_call_peaks(atac.fc, c, atac.minlen, 1000))
    assert "#5  Total candidate peaks : %d" % n in msgs(r)


@pytest.mark.parametrize("minlen", [0, 100, 400, 100000])
def test_minlen_filters_accessible_regions(quick, atac, model_file, minlen):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--save-states", "--minlen", minlen] + TUNE)
    assert r.returncode == 0, r.traceback
    regions = narrowpeak_regions(r.path("_accessible_regions.narrowPeak"))
    qualifying = triplets_from_states(read_states(r.path("_states.bed")),
                                      minlen)
    assert set(regions) <= set(qualifying)
    assert len(regions) >= len(qualifying) - 1
    if minlen == 100000:
        assert regions == []


def test_pileup_short(quick, atac):
    r1 = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
                "--pileup-short"], name="short")
    r2 = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only"],
               name="all")
    assert "#  Pile up ONLY short fragments" in msgs(r1)
    assert "#  Pile up all fragments" not in msgs(r1)
    assert "#  Pile up all fragments" in msgs(r2)
    assert read_lines(r1.path("_cutoff_analysis.tsv")) != \
        read_lines(r2.path("_cutoff_analysis.tsv"))


@pytest.mark.parametrize("seed", [10151, 7, 0])
def test_random_seed(quick, atac, seed):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--randomSeed", seed])
    m = msgs(r)
    # a seed of 0 is used but not announced in the command summary
    assert ("# Random seed selected as: %d" % seed in m) == (seed != 0)
    assert "# A random seed %d has been used in the sampling function" % \
        seed in m


@pytest.mark.parametrize("steps", [1, 2])
def test_decoding_steps(quick, atac, model_file, steps):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--decoding-steps", steps] + TUNE)
    assert r.returncode == 0, r.traceback
    m = msgs(r)
    n = int(log_value(r, "#   after expanding and merging, we have")
            .split()[7])
    decoding = [x for x in m if x.startswith("#    decoding ")]
    assert decoding == ["#    decoding %d..." % (k * steps)
                        for k in range(1, -(-n // steps) + 1)]


def test_decoding_steps_default_single_batch(base):
    assert [x for x in msgs(base) if x.startswith("#    decoding ")] == \
        ["#    decoding 5000..."]


def test_decoding_steps_zero_fails_assertion(quick, atac, model_file):
    """Pins the current output.

    --decoding-steps is not validated and its help gives no lower bound.
    With 0 regions per batch the first batch holds no bins, so
    extract_signals_from_regions fails its `assert nn > 0` and nothing is
    written.
    """
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--decoding-steps", "0"] + TUNE)
    assert r.returncode == 1
    assert isinstance(r.exc, AssertionError)
    assert r.files == []


def test_blacklist_removes_fragments(quick, atac, model_file, tmp_path):
    bl = write_rows(tmp_path / "bl.bed", [("chr1", 15500, 16500),
                                          ("chr2", 22800, 23200)])
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "-e", bl] + TUNE)
    assert r.returncode == 0, r.traceback
    n = sum(1 for c, l, rr in atac.frags
            if (c == b"chr1" and l < 16500 and 15500 < rr) or
            (c == b"chr2" and l < 23200 and 22800 < rr))
    m = msgs(r)
    assert "#  Read blacklist file..." in m
    assert "#  We removed %d fragments overlapping with blacklisted " \
        "regions." % n in m
    assert "#  There are %d fragments left." % (len(atac.frags) - n) in m
    regions = narrowpeak_regions(r.path("_accessible_regions.narrowPeak"))
    for c, s, e in regions:
        assert not (c == "chr1" and s < 16500 and 15500 < e)
        assert not (c == "chr2" and s < 23200 and 22800 < e)


def test_remove_dup(quick, atac, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--remove-dup"] + TUNE)
    assert r.returncode == 0, r.traceback
    n_dup = len(atac.frags) - len(set(atac.frags))
    assert n_dup >= 25
    m = msgs(r)
    assert "#  Removing duplicated fragments..." in m
    assert "#  We removed %d duplicated fragments." % n_dup in m
    assert "#  There are %d fragments left." % len(set(atac.frags)) in m


def test_no_remove_dup_by_default(base):
    assert "#  Removing duplicated fragments..." not in msgs(base)


@pytest.mark.parametrize("verbose, shown", [(0, False), (1, False),
                                            (2, True), (3, True)])
def test_verbose(quick, atac, verbose, shown):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
               "--verbose", verbose])
    m = msgs(r)
    assert ("#1 Read fragments from BEDPE file..." in m) == shown
    assert ("#3 Generate cutoff analysis report from %d fragments"
            % len(atac.frags) in m) == shown
    assert r.path("_cutoff_analysis.tsv").exists()


def test_verbose_does_not_change_outputs(quick, atac):
    r0 = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only",
                "--verbose", "0"], name="v0")
    r2 = quick(["-i", atac.bedpe, "-f", "BEDPE", "--cutoff-analysis-only"],
               name="v2")
    assert read_lines(r0.path("_cutoff_analysis.tsv")) == \
        read_lines(r2.path("_cutoff_analysis.tsv"))


def test_buffer_size(quick, atac, base, model_file):
    r = quick(["-i", atac.bedpe, "-f", "BEDPE", "--model", model_file,
               "--buffer-size", "7"] + TUNE)
    assert r.returncode == 0, r.traceback
    assert read_lines(r.path("_accessible_regions.narrowPeak")) == \
        read_lines(base.path("_accessible_regions.narrowPeak"))


# ------------------------------------
# macs3 hmmratac: through bin/macs3
# ------------------------------------

def test_cli_run_end_to_end(run_macs3, parse_log, atac, tmp_path, full):
    out = tmp_path / "new" / "dir"
    p = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE", "-n", "cli",
                   "--outdir", out, "--save-states"] + TUNE, timeout=60)
    assert p.returncode == 0, p.stderr
    assert sorted(os.listdir(out)) == expected_files(
        "cli", ["narrowpeak", "cutoff", "states", "model"])
    log = parse_log(p.stderr)
    assert log[-1] == ("INFO", "# Finished")
    assert ("", "# Command line: hmmratac -i %s -f BEDPE -n cli --outdir %s "
            "--save-states %s" % (atac.bedpe, out, " ".join(TUNE))) in log
    # the cutoff analysis does not depend on the HMM
    assert read_lines(out / "cli_cutoff_analysis.tsv") == \
        read_lines(full.path("_cutoff_analysis.tsv"))


def test_cli_model_json_reproducible_same_threads(run_macs3, atac,
                                                 tmp_path):
    texts = []
    for k in range(2):
        out = tmp_path / ("r%d" % k)
        p = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE",
                       "--outdir", out, "--modelonly"] + TUNE,
                      timeout=60, env={"OMP_NUM_THREADS": "1"})
        assert p.returncode == 0, p.stderr
        texts.append((out / "NA_model.json").read_text())
    assert texts[0] == texts[1]


def test_cli_model_json_thread_count_changes_only_rounding(run_macs3, atac,
                                                          tmp_path):
    # The HMM is initialised by scikit-learn's KMeans (through hmmlearn),
    # which splits its sums over OpenMP threads, so OMP_NUM_THREADS=1 and
    # 2 give model files that differ in the last digits (below 1e-12
    # relative on this input). Otherwise the model is the same.
    models = []
    for k, threads in enumerate(("1", "2")):
        out = tmp_path / ("r%d" % k)
        p = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE",
                       "--outdir", out, "--modelonly"] + TUNE,
                      timeout=60, env={"OMP_NUM_THREADS": threads})
        assert p.returncode == 0, p.stderr
        with open(out / "NA_model.json") as fh:
            models.append(json.load(fh))
    a, b = models
    assert list(a) == list(b)
    for key in a:
        if isinstance(a[key], list):
            np.testing.assert_allclose(np.array(a[key]), np.array(b[key]),
                                       rtol=1e-9, atol=1e-12)
        else:
            assert a[key] == b[key]


def test_cli_cutoff_analysis_only_exit_status(run_macs3, atac, tmp_path):
    """Regression test: --cutoff-analysis-only exited with status 1 after
    writing its report (--modelonly, the other documented early stop,
    exits 0).

    Fixed upstream in 4a2f0b5 (#731, issue #704).
    """
    p = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE",
                   "--outdir", tmp_path / "o", "--cutoff-analysis-only"],
                  timeout=60)
    assert (tmp_path / "o" / "NA_cutoff_analysis.tsv").is_file()
    assert p.returncode == 0


def test_cli_modelonly_exit_status(run_macs3, atac, tmp_path):
    p = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE",
                   "--outdir", tmp_path / "o", "--modelonly"] + TUNE,
                  timeout=60)
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "o" / "NA_model.json").is_file()


VALIDATION = [
    (["--means", "-1", "200", "400", "600"],
     " `--means` should not be negative! "),
    # the message spells the option `--stddev`, so its wording is not
    # checked
    (["--stddevs", "20", "-5", "20", "20"], None),
    (["--min-frag-p", "0"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--min-frag-p", "1"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--min-frag-p", "-0.5"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--binsize", "0"], " `--binsize` must be larger than 0."),
    (["-l", "-1"], " `-l` or `--lower` should not be negative! "),
    (["-u", "-1", "-l", "0"], " `-u` or `--upper` should not be negative! "),
    (["--maxTrain", "0"], " `--maxTrain` should be larger than 0!"),
    (["-c", "1"], " In order to use -c or --prescan-cutoff, the cutoff "
     "must be larger than 1."),
    (["--minlen", "-1"], " In order to use --minlen, the length should not "
     "be negative."),
]


@pytest.mark.parametrize("args, message", VALIDATION)
def test_cli_option_validation(run_macs3, parse_log, atac, tmp_path, args,
                               message):
    out = tmp_path / "o"
    p = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE", "--outdir",
                   out] + args, timeout=60)
    assert p.returncode == 1
    log = parse_log(p.stderr)
    if message is None:
        assert [level for level, _ in log] == ["ERROR"]
    else:
        assert log == [("ERROR", message)]
    assert os.listdir(out) == []


ARGPARSE_ERRORS = [
    ([], "the following arguments are required: -i/--input"),
    (["-f", "BAM"], "argument -f/--format: invalid choice: 'BAM'"),
    (["--hmm-type", "gamma"], "argument --hmm-type: invalid choice: 'gamma'"),
    (["--means", "50", "200", "400"], "argument --means: expected 4 "
     "arguments"),
    (["--stddevs", "20"], "argument --stddevs: expected 4 arguments"),
    (["--binsize", "ten"], "argument --binsize: invalid int value: 'ten'"),
    (["--min-frag-p", "small"], "argument --min-frag-p: invalid float value: "
     "'small'"),
    (["-c", "x"], "argument -c/--prescan-cutoff: invalid float value: 'x'"),
    (["--max-count", "1.5"], "argument --max-count: invalid int value: "
     "'1.5'"),
    (["--no-such-option"], "unrecognized arguments: --no-such-option"),
]


@pytest.mark.parametrize("args, message", ARGPARSE_ERRORS)
def test_cli_argparse_errors(run_macs3, atac, args, message):
    # every case gives -i except the one about the missing -i
    base_args = ["-i", atac.bedpe] if args else []
    p = run_macs3(["hmmratac"] + base_args + args, timeout=60)
    assert p.returncode == 2
    # unknown options are reported by the top-level parser
    assert "usage: macs3" in p.stderr
    assert message in p.stderr


def test_cli_help_lists_every_option(run_macs3):
    p = run_macs3(["hmmratac", "-h"], timeout=60)
    assert p.returncode == 0
    for opt in ["-i", "--input", "-f", "--format", "--barcodes",
                "--max-count", "--outdir", "-n", "--name",
                "--cutoff-analysis-only", "--cutoff-analysis-max",
                "--cutoff-analysis-steps", "--save-digested", "--save-states",
                "--save-likelihoods", "--save-training-data", "--no-fragem",
                "--means", "--stddevs", "--min-frag-p", "--binsize", "-u",
                "--upper", "-l", "--lower", "--maxTrain",
                "--training-flanking", "-t", "--training", "--model",
                "--modelonly", "--hmm-type", "-c", "--prescan-cutoff",
                "--minlen", "--pileup-short", "--randomSeed",
                "--decoding-steps", "-e", "--blacklist", "--remove-dup",
                "--verbose", "--buffer-size"]:
        assert re.search(r"(^|[\s,\[])%s([\s,\]]|$)" % re.escape(opt),
                         p.stdout, re.M), opt


def test_parser_defaults(macs3_argparser):
    ns = macs3_argparser.parse_args(["hmmratac", "-i", "x.bam"])
    assert ns.input_file == ["x.bam"]
    assert (ns.format, ns.barcodefile, ns.maxcount, ns.outdir, ns.name) == \
        ("BAMPE", "", None, "", "NA")
    assert (ns.cutoff_analysis_only, ns.cutoff_analysis_max,
            ns.cutoff_analysis_steps) == (False, 100, 100)
    assert (ns.save_digested, ns.save_states, ns.save_likelihoods,
            ns.save_train) == (False, False, False, False)
    assert (ns.em_skip, ns.em_means, ns.em_stddevs, ns.min_frag_p) == \
        (False, [50, 200, 400, 600], [20, 20, 20, 20], 0.001)
    assert (ns.hmm_binsize, ns.hmm_upper, ns.hmm_lower, ns.hmm_maxTrain,
            ns.hmm_training_flanking, ns.hmm_training_regions, ns.hmm_file,
            ns.hmm_modelonly, ns.hmm_type) == \
        (10, 20, 10, 1000, 1000, None, None, False, "gaussian")
    assert (ns.prescan_cutoff, ns.openregion_minlen, ns.pileup_short,
            ns.hmm_randomSeed, ns.decoding_steps, ns.blacklist,
            ns.misc_remove_duplicates, ns.verbose, ns.buffer_size) == \
        (1.2, 100, False, 10151, 5000, None, False, 2, 100000)


# ------------------------------------
# macs3 hmmratac on the yeast test data (slow)
# ------------------------------------

def jaccard(a, b):
    """Base-pair Jaccard index of two lists of (chrom, start, end)."""
    inter = union = 0
    for chrom in set(r[0] for r in a) | set(r[0] for r in b):
        end = max(r[2] for r in a + b if r[0] == chrom)
        ma = np.zeros(end, dtype=bool)
        mb = np.zeros(end, dtype=bool)
        for c, s, e in a:
            if c == chrom:
                ma[s:e] = True
        for c, s, e in b:
            if c == chrom:
                mb[s:e] = True
        inter += int((ma & mb).sum())
        union += int((ma | mb).sum())
    return inter / union


def standard_regions(test_dir, name):
    rows = []
    with open(test_dir / "standard_results_hmmratac" /
              ("%s_accessible_regions.narrowPeak" % name)) as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            f = line.split("\t")
            rows.append((f[0], int(f[1]), int(f[2])))
    return rows


@pytest.mark.slow
@pytest.mark.parametrize("inp, fmt, extra, std", [
    ("yeast_500k_SRR1822137.bam", "BAMPE", [], "hmmratac_yeast500k"),
    ("yeast_500k_SRR1822137.bedpe.gz", "BEDPE", [],
     "hmmratac_yeast500k_bedpe"),
    ("yeast_500k_SRR1822137.bam", "BAMPE", ["--hmm-type", "poisson"],
     "hmmratac_yeast500k_poisson"),
])
def test_yeast_matches_standard_results(run_macs3, test_dir, tmp_path, inp,
                                        fmt, extra, std):
    """Pins the current output.

    The reference is MACS3's own test/standard_results_hmmratac (there is
    no independent expectation for an HMM fit on real data). Those files
    use an older peak naming, so the regions are compared with the
    Jaccard index > 0.99 that test/cmdlinetest uses.

    Upstream 248fd6d added --jump (default 0.5; HMMR_EM's was 1.5) and
    made test/cmdlinetest pass --jump 1.5 to the runs that produce these
    standard results, so this test passes --jump 1.5 too.
    """
    p = run_macs3(["hmmratac", "-i", test_dir / inp, "-f", fmt, "-n", "y",
                   "--jump", "1.5", "--outdir", tmp_path] + extra,
                  timeout=300)
    assert p.returncode == 0, p.stderr
    got = narrowpeak_regions(tmp_path / "y_accessible_regions.narrowPeak")
    assert jaccard(got, standard_regions(test_dir, std)) > 0.99


@pytest.mark.slow
def test_yeast_frag_with_barcodes(run_macs3, test_dir, tmp_path):
    """Pins the current output.

    As test_yeast_matches_standard_results, for the scATAC fragment file
    run of test/cmdlinetest, which passes --jump 1.5 since upstream
    248fd6d.
    """
    p = run_macs3(["hmmratac", "-i", test_dir / "test.fragments.tsv.gz",
                   "-f", "FRAG", "--barcodes", test_dir / "barcodes.txt",
                   "--hmm-type", "poisson", "--jump", "1.5", "-n", "sc",
                   "--outdir", tmp_path],
                  timeout=300)
    assert p.returncode == 0, p.stderr
    got = narrowpeak_regions(tmp_path / "sc_accessible_regions.narrowPeak")
    assert jaccard(got, standard_regions(test_dir,
                                         "hmmratac_scatac_test")) > 0.99


@pytest.mark.slow
def test_yeast_frag_file(run_macs3, parse_log, test_dir, tmp_path):
    # the yeast fragments in FRAG format (barcodes AAAA... and BBBB...),
    # keeping barcode BBBB... with --barcodes and capping counts at 1
    frag = test_dir / "yeast_500k_SRR1822137.frag.gz"
    bc = test_dir / "yeast_500k_SRR1822137.barcode.txt"
    keep = bc.read_text().split()
    with gzip.open(frag, "rt") as fh:
        n = sum(1 for line in fh if line.split("\t")[3] in keep)
    p = run_macs3(["hmmratac", "-i", frag, "-f", "FRAG", "--barcodes", bc,
                   "--max-count", "1", "-n", "yf", "--outdir", tmp_path],
                  timeout=300)
    assert p.returncode == 0, p.stderr
    log = parse_log(p.stderr)
    assert ("INFO", "#   extracted %d fragments in treatment" % n) in log
    assert ("INFO", "#  Read %d fragments." % n) in log
    lines = read_lines(tmp_path / "yf_accessible_regions.narrowPeak")
    assert len(lines) > 50
    assert all(len(x.split("\t")) == 10 for x in lines)


