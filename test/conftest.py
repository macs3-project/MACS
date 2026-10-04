"""Shared fixtures for the MACS3 test suite.

The existing test files open data with paths relative to the repository
root (``test/tiny.bed.gz``), so ``pytest`` is expected to run from the
repository root. New tests use the absolute paths given by the
``test_dir`` and ``data_dir`` fixtures and run from anywhere.

Command-line tests run ``bin/macs3`` in a subprocess with the same
interpreter and the same ``MACS3`` package that the test process
imported, so a suite run against an in-place build exercises that
build. Each subprocess gets ``PYTHONHASHSEED=0`` (``filterdup`` and
``randsample`` iterate a set in hash order), capped thread pools, and
``TMPDIR`` inside the test's ``tmp_path``.

Tests slower than about 2 s are marked ``slow`` and skipped unless
``--run-slow`` is given.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

# Cap thread pools before numpy/scipy/scikit-learn are imported by any test
# module. Two threads keep the shared machine responsive, and scikit-learn's
# KMeans (hmmratac's HMM initialisation) gives results that depend on thread
# arrival order with 3 or more OpenMP threads, which would make in-process
# hmmratac tests nondeterministic. Set MACS3_TEST_THREADS to change it.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_var] = os.environ.get("MACS3_TEST_THREADS", "2")

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent
DATA_DIR = TEST_DIR / "data"
MACS3_SCRIPT = REPO_ROOT / "bin" / "macs3"


# ------------------------------------
# markers and options
# ------------------------------------

def pytest_addoption(parser):
    parser.addoption("--run-slow", action="store_true", default=False,
                     help="also run tests marked 'slow' (> ~2 s each)")


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: takes more than about 2 s; skipped unless --run-slow is given")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="slow test; use --run-slow to run it")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


# ------------------------------------
# paths
# ------------------------------------

@pytest.fixture(scope="session")
def test_dir():
    """Absolute path of the ``test/`` directory."""
    return TEST_DIR


@pytest.fixture(scope="session")
def data_dir():
    """Absolute path of ``test/data/`` (small inputs added for unit tests)."""
    return DATA_DIR


# ------------------------------------
# running the macs3 command
# ------------------------------------

def macs3_env(tmpdir=None, extra=None):
    """Environment for a ``macs3`` subprocess.

    ``PYTHONPATH`` points at the directory holding the imported
    ``MACS3`` package, so the subprocess runs the same build.
    """
    import MACS3
    pkg_parent = str(Path(MACS3.__file__).resolve().parent.parent)
    env = os.environ.copy()
    env["PYTHONPATH"] = pkg_parent + (os.pathsep + env["PYTHONPATH"]
                                      if env.get("PYTHONPATH") else "")
    env["PYTHONHASHSEED"] = "0"
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS"):
        env[var] = os.environ.get("MACS3_TEST_THREADS", "2")
    if tmpdir is not None:
        env["TMPDIR"] = str(tmpdir)
    if extra:
        env.update(extra)
    return env


def no_core_dumps():
    """``preexec_fn`` for subprocesses: some tests drive MACS3 into a
    crash (segfault, C assertion abort); do not leave core files."""
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    except (ImportError, ValueError, OSError):
        pass


@pytest.fixture
def macs3_tmpdir(tmp_path):
    """The TMPDIR given to ``macs3`` subprocesses of this test."""
    d = tmp_path / "TMPDIR"
    d.mkdir(exist_ok=True)
    return d


@pytest.fixture
def run_macs3(tmp_path, macs3_tmpdir):
    """Return a function that runs ``macs3 <args>`` in a subprocess.

    ``run_macs3(["callpeak", "-t", f, ...], timeout=120, cwd=None,
    env=None, stdin=None)`` returns the ``subprocess.CompletedProcess``
    with text ``stdout``/``stderr``. ``cwd`` defaults to ``tmp_path``.
    It never raises on a non-zero exit code; tests assert on
    ``returncode`` themselves. ``subprocess.TimeoutExpired`` propagates.
    Core dumps are disabled in the child.
    """
    def _run(args, timeout=120, cwd=None, env=None, stdin=None):
        cmd = [sys.executable, str(MACS3_SCRIPT)] + [str(a) for a in args]
        return subprocess.run(cmd, cwd=str(cwd or tmp_path),
                              env=macs3_env(macs3_tmpdir, env),
                              input=stdin, capture_output=True, text=True,
                              timeout=timeout, preexec_fn=no_core_dumps)
    return _run


@pytest.fixture(scope="session")
def macs3_argparser():
    """The ``argparse`` parser built by ``bin/macs3``'s ``prepare_argparser``."""
    import importlib.machinery
    import importlib.util
    loader = importlib.machinery.SourceFileLoader("macs3_main",
                                                  str(MACS3_SCRIPT))
    spec = importlib.util.spec_from_loader("macs3_main", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod.prepare_argparser()


# ------------------------------------
# writing small synthetic inputs
# ------------------------------------

def _open_w(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return open(path, "w")


def write_lines(path, lines, trailing_newline=True):
    """Write ``lines`` (iterables of fields or strings) as a TSV file.

    A ``.gz`` suffix gzips the file. Returns the path as ``str``.
    """
    rows = ["\t".join(map(str, x)) if not isinstance(x, str) else x
            for x in lines]
    text = "\n".join(rows)
    if trailing_newline and rows:
        text += "\n"
    with _open_w(path) as fh:
        fh.write(text)
    return str(path)


@pytest.fixture
def write_bed(tmp_path):
    """``write_bed(rows, name="reads.bed")``: rows are
    (chrom, start, end, name, score, strand) or shorter tuples."""
    def _w(rows, name="reads.bed", trailing_newline=True):
        return write_lines(tmp_path / name, rows, trailing_newline)
    return _w


@pytest.fixture
def write_bedpe(tmp_path):
    """``write_bedpe(rows, name="frags.bedpe")``: rows are (chrom, start, end)."""
    def _w(rows, name="frags.bedpe", trailing_newline=True):
        return write_lines(tmp_path / name, rows, trailing_newline)
    return _w


@pytest.fixture
def write_frag(tmp_path):
    """``write_frag(rows, name="frags.tsv")``: rows are
    (chrom, start, end, barcode, count)."""
    def _w(rows, name="frags.tsv", trailing_newline=True):
        return write_lines(tmp_path / name, rows, trailing_newline)
    return _w


@pytest.fixture
def write_bedgraph(tmp_path):
    """``write_bedgraph(rows, name="x.bdg", header=None)``: rows are
    (chrom, start, end, value); ``header`` is an optional first line."""
    def _w(rows, name="x.bdg", header=None):
        lines = ([header] if header else []) + list(rows)
        return write_lines(tmp_path / name, lines)
    return _w


@pytest.fixture
def make_alignments(tmp_path):
    """Write a SAM or BAM file with pysam.

    ``make_alignments(reads, refs=(("chr1", 100000),), name="reads.bam",
    fmt="bam", index=False)``. Each read is a dict with keys
    ``name, ref, pos (0-based), flag, cigar ("36M"), mapq (default 30),
    seq (default "A" * read length), next_ref, next_pos, tlen`` and an
    optional ``tags`` list of (tag, value). Reads are coordinate sorted
    before writing; ``index=True`` also writes a ``.bai`` (BAM only).
    Returns the path as ``str``.
    """
    pysam = pytest.importorskip("pysam")

    def _make(reads, refs=(("chr1", 100000),), name="reads.bam",
              fmt="bam", index=False):
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": r, "LN": ln} for r, ln in refs]}
        refidx = {r: i for i, (r, _) in enumerate(refs)}
        path = str(tmp_path / name)
        mode = "wb" if fmt == "bam" else "w"
        recs = sorted(reads, key=lambda r: (refidx.get(r.get("ref"), 1 << 30),
                                            r.get("pos", -1)))
        with pysam.AlignmentFile(path, mode, header=header) as out:
            for r in recs:
                a = pysam.AlignedSegment(out.header)
                a.query_name = r["name"]
                a.flag = r.get("flag", 0)
                cigar = r.get("cigar", "36M")
                if r.get("ref") is not None and r["ref"] != "*":
                    a.reference_id = refidx[r["ref"]]
                    a.reference_start = r["pos"]
                else:
                    a.reference_id = -1
                    a.reference_start = -1
                a.mapping_quality = r.get("mapq", 30)
                a.cigarstring = cigar if cigar != "*" else None
                if "seq" in r:
                    seq = r["seq"]
                else:
                    qlen = a.query_length if cigar != "*" else 36
                    seq = "A" * (qlen or 36)
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array(
                    r.get("qual", "I" * len(seq)))
                nr = r.get("next_ref")
                a.next_reference_id = refidx[nr] if nr not in (None, "*") else -1
                a.next_reference_start = r.get("next_pos", -1)
                a.template_length = r.get("tlen", 0)
                for tag, val in r.get("tags", []):
                    a.set_tag(tag, val)
                out.write(a)
        if index and fmt == "bam":
            pysam.index(path)
        return path
    return _make


def pe_pair(name, ref, start, end, readlen=36, mapq=30):
    """Two pysam read dicts for a properly paired fragment [start, end)."""
    tlen = end - start
    r1 = dict(name=name, ref=ref, pos=start, flag=99,
              cigar="%dM" % readlen, mapq=mapq, next_ref=ref,
              next_pos=end - readlen, tlen=tlen)
    r2 = dict(name=name, ref=ref, pos=end - readlen, flag=147,
              cigar="%dM" % readlen, mapq=mapq, next_ref=ref,
              next_pos=start, tlen=-tlen)
    return [r1, r2]


@pytest.fixture
def make_pe_pair():
    """``make_pe_pair(name, ref, start, end, readlen=36)`` -> [r1, r2] dicts
    for ``make_alignments`` (flags 99/147, proper pair)."""
    return pe_pair


def read_table(path, comment="#"):
    """Read a (possibly gzipped) TSV into a list of field lists, skipping
    blank lines, lines starting with ``comment`` and ``track`` lines."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(comment) or line.startswith("track"):
                continue
            rows.append(line.split("\t"))
    return rows


@pytest.fixture
def read_tsv():
    """``read_tsv(path)`` -> list of field lists (skips comments, track lines)."""
    return read_table


def log_messages(stderr):
    """Split ``macs3`` log output into (LEVEL, message) pairs.

    Drops the timestamp and the ``[N MB]`` memory prefix that MACS3's
    logger adds, so messages can be compared across runs. Lines that
    are not log records (tracebacks, argparse usage) are returned with
    level ``""``.
    """
    import re
    pat = re.compile(r"^(\w+)\s+@ \d\d \w\w\w \d{4} \d\d:\d\d:\d\d: "
                     r"(?:\[\d+ MB\] )?(.*?) ?$")
    out = []
    for line in stderr.splitlines():
        m = pat.match(line)
        if m:
            out.append((m.group(1), m.group(2)))
        elif line.strip():
            out.append(("", line))
    return out


@pytest.fixture
def parse_log():
    """``parse_log(stderr)`` -> list of (LEVEL, message) without timestamps."""
    return log_messages
