"""Bounded public FRAG EOF regressions for #765 / #757 item 7."""
import gzip
import struct
import subprocess
import sys
import zlib

import pytest

from conftest import macs3_env


def write_frag(path, data, encoding):
    if encoding == "plain":
        path.write_bytes(data)
    elif encoding == "gzip":
        path.write_bytes(gzip.compress(data))
    else:
        # BGZF members have a 64 KiB size limit. No optional pysam dependency.
        blocks = []
        for offset in range(0, len(data), 60000):
            chunk = data[offset:offset + 60000]
            compressor = zlib.compressobj(wbits=-15)
            payload = compressor.compress(chunk) + compressor.flush()
            blocks.append(bytes.fromhex("1f8b08040000000000ff060042430200")
                          + struct.pack("<H", len(payload) + 25) + payload
                          + struct.pack("<II", zlib.crc32(chunk), len(chunk)))
        blocks.append(bytes.fromhex(
            "1f8b08040000000000ff0600424302001b0003000000000000000000"))
        path.write_bytes(b"".join(blocks))
    return path


SCRIPT = r'''
from pathlib import Path
import sys
from MACS3.IO.Parser import FragParser
from MACS3.Signal.PairedEndTrack import PETrackII

mode, operation, cap = sys.argv[1:4]
cap = int(cap)
paths = sys.argv[4:]


def snapshot(track):
    track.finalize()
    return (dict(track.barcode_dict),
            {chrom: (track.get_locations_by_chr(chrom).tolist(),
                     track.barcodes[chrom].tolist())
             for chrom in track.get_chr_names()})


def load(path):
    parser = FragParser(path)
    if operation == 'build':
        track = parser.build_petrack(max_count=cap)
    else:
        track = PETrackII(buffer_size=2)
        track.add_loc(b'chr0', 1, 11, b'A', 1)
        result = parser.append_petrack(track, max_count=cap)
        assert result is track
    return parser, track


if mode == 'valid':
    results = []
    for path in paths:
        p, track = load(path)
        results.append((p.n, p.d, snapshot(track)))
    assert results[0] == results[1]
    n, d, (barcodes, records) = results[0]
    assert n == 3 and d == 40.0
    assert barcodes == {b'A': 0, b'B': 1}
    assert records[b'chr1'] == ([(10, 50, min(7, cap) if cap else 7),
                                (60, 110, 1)], [0, 1])
    assert records[b'chr2'] == ([(0, 30, min(3, cap) if cap else 3)], [0])

elif mode == 'edges':
    results = []
    for path in paths:
        p, track = load(path)
        results.append((p.n, p.d, snapshot(track)))
    assert results[0] == results[1]
    n, d, (barcodes, records) = results[0]
    assert n == 1 and d == 40.0
    assert records[b'chr1'] == ([(10, 50, 7)], [0])
    assert barcodes == {b'A': 0}

elif mode == 'malformed':
    errors = []
    for path in paths:
        try:
            load(path)
        except Exception as exc:
            errors.append((type(exc).__name__, str(exc)))
        else:
            raise AssertionError('malformed record accepted')
    assert errors[0] == errors[1], errors
    assert errors[0][0] == ('Exception' if cap == 0 else 'AssertionError')
    assert ('Less than 5 columns' if cap == 0 else
            'Right position must be larger') in errors[0][1]

elif mode == 'empty':
    for path in paths:
        try:
            load(path)
        except ValueError as exc:
            assert str(exc) == 'No valid fragments found in FRAG input'
        else:
            raise AssertionError('empty input accepted')

elif mode == 'empty_append_with_history':
    parser, track = load(paths[0])
    before = (parser.n, parser.d, snapshot(track))
    for path in paths[1:]:
        # Reopen the same public parser instance; n/d retain its history.
        parser.__init__(path)
        assert parser.append_petrack(track) is track
        assert (parser.n, parser.d, snapshot(track)) == before

elif mode == 'fallback':
    class Receiver:
        def __init__(self):
            self.records = []
        def add_loc(self, *record):
            self.records.append(record)
        def set_rlengths(self, lengths):
            pass
    results = []
    for path in paths:
        track = Receiver()
        p = FragParser(path)
        assert p.append_petrack(track, max_count=2) is track
        results.append((p.n, p.d, track.records))
    assert results[0] == results[1]
    assert results[0][2][-1] == (b'chr2', 0, 30, b'A', 2)

elif mode in ('parallel', 'boundary'):
    results = []
    for path in paths:
        p, track = load(path)
        results.append((p.n, p.d, snapshot(track)))
    assert results[0] == results[1]
    n, d, (barcodes, records) = results[0]
    if mode == 'parallel':
        assert n == 10001 and d == 40.0
        assert barcodes == {b'A': 0, b'LAST': 1}
        assert records[b'chr2'] == ([(2, 42, 9)], [1])
    else:
        assert n == 2 and d == 40.0
        assert len(barcodes) == 2
        assert records[b'chr2'][0] == [(2, 42, 9)]
'''


def probe(tmp_path, mode, operation, cap, paths):
    script = tmp_path / "probe.py"
    script.write_text(SCRIPT)
    result = subprocess.run(
        [sys.executable, str(script), mode, operation, str(cap),
         *map(str, paths)], cwd=tmp_path, env=macs3_env(),
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("operation", ["build", "append"])
@pytest.mark.parametrize("cap", [0, 2])
def test_frag_final_record(tmp_path, encoding, operation, cap):
    data = (b'# sample\ntrack name=sample\nbrowser position chr1\n'
            b'chr1\t10\t50\tA\t7\nchr1\t60\t110\tB\t1\n'
            b'chr2\t0\t30\tA\t3\textra\tfields   ')
    paths = [write_frag(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, "valid", operation, cap, paths)


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("operation", ["build", "append"])
@pytest.mark.parametrize("bad,cap", [(b'chr1\t10\t50\tA', 0),
                                     (b'chr1\t50\t10\tA\t1', 2)])
def test_frag_malformed_final_record(tmp_path, encoding, operation, bad, cap):
    data = b'chr1\t1\t41\tA\t1\n' + bad
    paths = [write_frag(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, "malformed", operation, cap, paths)


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("operation", ["build", "append"])
def test_frag_no_accepted_records(tmp_path, encoding, operation):
    inputs = [b'', b'# header', b'# header\n',
              b'track name=x\nbrowser position chr1',
              b'chr1\t-1\t40\tA\t1', b'\t0\t40\tA\t1\n']
    paths = [write_frag(tmp_path / str(i), data, encoding)
             for i, data in enumerate(inputs)]
    probe(tmp_path, "empty", operation, 0, paths)


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
def test_frag_final_record_fallback_receiver(tmp_path, encoding):
    data = b'chr1\t10\t50\tA\t7\nchr2\t0\t30\tA\t3'
    paths = [write_frag(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, "fallback", "append", 0, paths)


@pytest.mark.parametrize("mode,encoding", [("parallel", "bgzf"),
                                             ("boundary", "plain"),
                                             ("boundary", "gzip"),
                                             ("boundary", "bgzf")])
def test_frag_final_record_across_blocks(tmp_path, mode, encoding):
    # ~0.1 s per probe. The first fixture exceeds the worker chunk threshold;
    # the second exceeds the 4 MiB refill window in one unterminated record.
    if mode == "parallel":
        data = b'chr1\t1\t41\tA\t1\n' * 10000 + b'chr2\t2\t42\tLAST\t9'
    else:
        data = b'chr1\t1\t41\tA\t1\nchr2\t2\t42\t' + b'B' * (4194304 + 1) + b'\t9'
    paths = [write_frag(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, mode, "build", 0, paths)


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
def test_frag_pileup_final_record(tmp_path, run_macs3, encoding):
    data = b'chr1\t10\t50\tA\t2\nchr1\t30\t70\tB\t3'
    outputs = []
    for i, ending in enumerate([b'\n', b'']):
        path = write_frag(tmp_path / str(i), data + ending, encoding)
        output = tmp_path / (str(i) + '.bdg')
        result = run_macs3(['pileup', '-f', 'FRAG', '-i', path, '-o', output],
                          timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append(output.read_text())
    assert outputs[0] == outputs[1]
    assert outputs[0].strip()


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("tail", [b'', b'\nchr2\t-1\t30\tB\t3',
                                  b'\n\t0\t30\tB\t3'])
def test_frag_single_record_and_skipped_tail(tmp_path, encoding, tail):
    data = b'chr1\t10\t50\tA\t7' + tail
    paths = [write_frag(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, "edges", "build", 0, paths)


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
def test_frag_empty_append_preserves_parser_history(tmp_path, encoding):
    inputs = [b'chr1\t10\t50\tA\t7', b'', b'# header', b'# header\n']
    paths = [write_frag(tmp_path / str(i), data, encoding)
             for i, data in enumerate(inputs)]
    probe(tmp_path, "empty_append_with_history", "build", 0, paths)
