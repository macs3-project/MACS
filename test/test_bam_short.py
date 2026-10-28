"""Public API and CLI regressions for #801, part of #757."""

import gzip
from pathlib import Path
import struct
import subprocess
import sys

import pytest

from conftest import macs3_env, no_core_dumps


def run_probe(tmp_path, path, source):
    script = tmp_path / "probe.py"
    script.write_text('''
from pathlib import Path
import sys
from MACS3.IO import Parser
assert Path(Parser.__file__).resolve() == Path(sys.argv[2]).resolve()
assert Path(Parser.__file__).suffix != '.py'
path = sys.argv[1]
''' + source)
    from MACS3.IO import Parser
    result = subprocess.run(
        [sys.executable, str(script), str(path), Parser.__file__],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True,
        text=True, timeout=30, preexec_fn=no_core_dumps)
    assert result.returncode == 0, result.stdout + result.stderr


def bam_reads(n):
    return [dict(name=f"r{i}", ref="chr1", pos=1000 + 100*i,
                 cigar=f"{50 + i}M", seq="A" * (50 + i),
                 qual="I" * (50 + i)) for i in range(n)]


@pytest.mark.parametrize("n", range(1, 13))
@pytest.mark.parametrize("compressed", [False, True], ids=["raw", "bgzf"])
def test_short_bam_public_api(make_alignments, tmp_path, n, compressed):
    path = Path(make_alignments(bam_reads(n)))
    if not compressed:
        raw = tmp_path / "raw.bam"
        raw.write_bytes(gzip.decompress(path.read_bytes()))
        path = raw
    # Retain the first-ten sampling rule and integer truncation of the mean.
    expected = sum(50 + i for i in range(min(n, 10))) // min(n, 10)
    run_probe(tmp_path, path, f'''
for factory in (Parser.BAMParser, Parser.guess_parser):
    p = factory(path)
    assert p.sniff() is True
    assert p.tsize() == {expected}
    track = p.build_fwtrack()
    track.finalize()
    plus, minus = track.get_locations_by_chr(b'chr1')
    assert plus.tolist() == {[1000 + 100*i for i in range(n)]!r}
    assert minus.tolist() == []
    assert p.tsize() == {expected}  # cached after the input stream closes
''')


@pytest.mark.parametrize("fmt", ["BAM", "AUTO"])
@pytest.mark.parametrize("explicit_size", [False, True])
@pytest.mark.parametrize("n", [1, 9, 10])
def test_short_bam_filterdup(make_alignments, run_macs3, tmp_path,
                             fmt, explicit_size, n):
    # Original reproduction: all reads have length 50.
    reads = bam_reads(n)
    for read in reads:
        read.update(cigar="50M", seq="A"*50, qual="I"*50)
    path = make_alignments(reads)
    out = tmp_path / "out.bed"
    args = ["filterdup", "-i", path, "-f", fmt, "-g", "1e6", "-o", out,
            "--keep-dup", "1"]
    if explicit_size:
        args += ["-s", "50"]
    result = run_macs3(args, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert out.read_text().splitlines() == [
        f"chr1\t{1000 + 100*i}\t{1050 + 100*i}\t.\t.\t+"
        for i in range(n)]


def test_original_filterdup_defaults(make_alignments, run_macs3, tmp_path):
    reads = bam_reads(9)
    for read in reads:
        read.update(cigar="50M", seq="A"*50, qual="I"*50)
    path = make_alignments(reads)
    outputs = []
    for extra in ([], ["-s", "50"]):
        out = tmp_path / "original.bed"
        result = run_macs3(["filterdup", "-i", path, "-f", "BAM", "-g", "1e6",
                            "-o", out] + extra, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append(out.read_text())
    # Preserve automatic duplicate filtering, including its small-sample cutoff.
    assert outputs[0] == outputs[1]
    assert outputs[0]


@pytest.mark.parametrize("compressed", [False, True], ids=["raw", "bgzf"])
def test_header_only_bam(make_alignments, tmp_path, compressed):
    path = Path(make_alignments([]))
    if not compressed:
        raw = tmp_path / "raw.bam"
        raw.write_bytes(gzip.decompress(path.read_bytes()))
        path = raw
    run_probe(tmp_path, path, '''
for method in ('tsize', 'sniff'):
    p = Parser.BAMParser(path)
    try:
        getattr(p, method)()
    except ValueError as error:
        assert 'No alignments' in str(error)
    else:
        raise AssertionError('Header-only BAM must report no alignments')
    finally:
        p.close()
''')


@pytest.mark.parametrize("fmt", ["BAM", "AUTO"])
@pytest.mark.parametrize("explicit_size", [False, True])
def test_header_only_bam_filterdup(make_alignments, run_macs3, tmp_path,
                                   fmt, explicit_size):
    args = ["filterdup", "-i", make_alignments([]), "-f", fmt, "-g", "1e6",
            "-o", tmp_path / "empty.bed"]
    if explicit_size:
        args += ["-s", "50"]
    result = run_macs3(args, timeout=30)
    assert result.returncode > 0
    assert "No alignments" in result.stderr
    assert "ZeroDivisionError" not in result.stderr
    assert "struct.error" not in result.stderr


@pytest.mark.parametrize("compressed", [False, True], ids=["raw", "gzip"])
@pytest.mark.parametrize("case", ["empty", "header", "size", "body",
                                 "short-core", "negative-size"])
def test_malformed_bam(make_alignments, tmp_path, case, compressed):
    header = gzip.decompress(Path(make_alignments([], name="empty.bam")).read_bytes())
    record = gzip.decompress(Path(make_alignments(bam_reads(1))).read_bytes())
    data = {"empty": b"", "header": b"BAM\1\0", "size": header + b"\1\2",
            "body": record[:-10], "short-core": header + struct.pack('<i', 20) + b'\0'*20,
            "negative-size": header + struct.pack('<i', -1)}[case]
    path = tmp_path / "malformed.bam"
    path.write_bytes(gzip.compress(data) if compressed else data)
    error = "struct.error" if case in ("empty", "header") else (
        "ValueError" if case in ("short-core", "negative-size") else "EOFError")
    run_probe(tmp_path, path, f'''
import struct
for method in ('tsize', 'sniff'):
    # An empty stream has no BAM magic, so sniff correctly returns False.
    if {case!r} == 'empty' and method == 'sniff':
        continue
    p = Parser.BAMParser(path)
    try:
        getattr(p, method)()
    except {error}:
        pass
    else:
        raise AssertionError('Malformed BAM must be rejected')
    finally:
        p.close()
''')
