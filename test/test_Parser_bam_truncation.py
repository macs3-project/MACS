"""Bounded public BAM framing regressions for #766; no pysam required."""
import gzip
import json
from pathlib import Path
import random
import struct
import subprocess
import sys
import zlib

import pytest

from conftest import macs3_env


HEADER = (b'BAM\1' + struct.pack('<ii', 0, 1) + struct.pack('<i', 5)
          + b'chr1\0' + struct.pack('<i', 100000))


def record(flag=99, aux=b''):
    body = (struct.pack('<iiIIiiii', 0, 10, (30 << 8) | 2,
                        (flag << 16) | 1, 1, 0, 40, 40)
            + b'r\0' + struct.pack('<I', 16) + b'\x10\x1e' + aux)
    return struct.pack('<i', len(body)) + body


def byte_array(data):
    return b'ZZBC' + struct.pack('<I', len(data)) + data


def write_bam(path, data, encoding, cuts=()):
    if encoding == 'raw':
        path.write_bytes(data)
    elif encoding == 'gzip':
        path.write_bytes(gzip.compress(data))
    else:
        blocks = []
        # Explicit cuts allow size fields and bodies to cross members.
        positions = sorted({0, len(data), *cuts})
        for left, right in zip(positions, positions[1:]):
            for offset in range(left, right, 60000):
                chunk = data[offset:min(offset + 60000, right)]
                compressor = zlib.compressobj(wbits=-15)
                payload = compressor.compress(chunk) + compressor.flush()
                blocks.append(bytes.fromhex('1f8b08040000000000ff060042430200')
                              + struct.pack('<H', len(payload) + 25) + payload
                              + struct.pack('<II', zlib.crc32(chunk), len(chunk)))
        blocks.append(bytes.fromhex(
            '1f8b08040000000000ff0600424302001b0003000000000000000000'))
        path.write_bytes(b''.join(blocks))
    if encoding != 'raw':
        assert gzip.decompress(path.read_bytes()) == data
    return str(path)


SCRIPT = r'''
import json
import sys
from MACS3.IO.Parser import BAMParser, BAMPEParser
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI

mode, operation, manifest = sys.argv[1:]
paired = mode == 'pe'
for path, error, count in json.loads(open(manifest).read()):
    parser = (BAMPEParser if paired else BAMParser)(path, buffer_size=2)
    track = None
    if operation == 'append':
        track = PETrackI(buffer_size=2) if paired else FWTrack(buffer_size=2)
        if paired:
            track.add_loc(b'chr0', 1, 11)
        else:
            track.add_loc(b'chr0', 1, 0)
    try:
        if paired:
            result = (parser.build_petrack() if track is None else
                      parser.append_petrack(track))
        else:
            result = (parser.build_fwtrack() if track is None else
                      parser.append_fwtrack(track))
    except EOFError as exc:
        assert error is not None, (path, str(exc))
        assert str(exc) == error, (path, str(exc), error)
        if paired:
            assert parser.n == 0 and parser.d == 0
    else:
        assert error is None, 'truncated input accepted: ' + path
        if track is not None:
            assert result is track
        result.finalize()
        assert result.total == count + (operation == 'append')
        if count:
            locations = result.get_locations_by_chr(b'chr1')
            if paired:
                assert locations.tolist() == [(10, 50)] * count
                assert parser.n == count and parser.d == 40
            else:
                assert locations[0].tolist() == [10] * count
                assert locations[1].tolist() == []
'''


def probe(tmp_path, mode, operation, cases):
    script = tmp_path / 'probe.py'
    script.write_text(SCRIPT)
    manifest = tmp_path / 'cases.json'
    manifest.write_text(json.dumps(cases))
    result = subprocess.run(
        [sys.executable, str(script), mode, operation, str(manifest)],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True,
        text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
@pytest.mark.parametrize('mode', ['se', 'pe'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_bam_truncated_records(tmp_path, encoding, mode, operation):
    complete = record()
    tails = [(complete[:n], 'Truncated BAM block_size field: '
              f'expected 4 bytes, found {n}') for n in (1, 2, 3)]
    tails += [(complete[:4 + n], 'Truncated BAM record body: '
               f'expected {len(complete) - 4} bytes, found {n}')
              for n in (0, 1, 10, 31, 32, len(complete) - 5)]
    cases = []
    # A kept prefix, no prefix, or only a filtered prefix must all fail.
    for prefix in (complete, b'', record(flag=4)):
        for tail, error in tails:
            path = write_bam(tmp_path / f'bad{len(cases)}.bam',
                             HEADER + prefix + tail, encoding)
            cases.append((path, error, None))
    probe(tmp_path, mode, operation, cases)


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
@pytest.mark.parametrize('mode', ['se', 'pe'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_bam_complete_records(tmp_path, encoding, mode, operation):
    cases = []
    for count in (1, 3):
        data = HEADER + record() * count + record(flag=4)
        start = len(HEADER)
        path = write_bam(tmp_path / f'good{count}.bam', data, encoding,
                         cuts=(start + 1, start + 3, start + 8,
                               start + len(record()) - 1))
        cases.append((path, None, count))
    if mode == 'se':
        cases.append((write_bam(tmp_path / 'empty.bam', HEADER, encoding),
                      None, 0))
    probe(tmp_path, mode, operation, cases)


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
@pytest.mark.parametrize('mode', ['se', 'pe'])
def test_bam_serial_window_boundary(tmp_path, encoding, mode):
    cases = []
    # Put the next size field at/before the 4 MiB decompressed boundary.
    for delta in (0, 1, 3, 20):
        target = 4 * 1024 * 1024 - delta
        filler = record(aux=byte_array(b'\0' * (
            target - len(HEADER) - len(record()) - 8)))
        assert len(HEADER + filler) == target
        prefix = HEADER + filler
        good = write_bam(tmp_path / f'good{delta}.bam',
                         prefix + record(), encoding)
        bad = write_bam(tmp_path / f'bad{delta}.bam',
                        prefix + record()[:-1], encoding)
        cases.extend([(good, None, 2),
                      (bad, 'Truncated BAM record body: expected 40 bytes, '
                       'found 39', None)])
    probe(tmp_path, mode, 'build', cases)


@pytest.mark.parametrize('mode', ['se', 'pe'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_bam_parallel_batch_boundary(tmp_path, mode, operation):
    # >4 MiB compressed and >16 MiB decompressed: exceeds both batch limits.
    # Also exercises a single record larger than the stream's initial buffer.
    payload = random.Random(766).randbytes(5 * 1024 * 1024)
    payload += b'\0' * (12 * 1024 * 1024)
    prefix = HEADER + record(aux=byte_array(payload))
    good = write_bam(tmp_path / 'good.bam', prefix + record(), 'bgzf')
    assert Path(good).stat().st_size > 4 * 1024 * 1024
    bad = write_bam(tmp_path / 'bad.bam', prefix + record()[:-1], 'bgzf')
    probe(tmp_path, mode, operation,
          [(good, None, 2), (bad, 'Truncated BAM record body: '
                                'expected 40 bytes, found 39', None)])


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
@pytest.mark.parametrize('fmt', ['BAM', 'BAMPE'])
@pytest.mark.parametrize('tail', [b'\x28\0', record()[:-1]],
                         ids=['size-field', 'body'])
def test_bam_truncation_cli(tmp_path, run_macs3, encoding, fmt, tail):
    data = HEADER + record() * 10  # BAM tsize samples ten complete records.
    good = write_bam(tmp_path / 'good.bam', data, encoding)
    bad = write_bam(tmp_path / 'bad.bam', data + tail, encoding)
    output = tmp_path / 'pileup.bdg'
    args = ['pileup', '-f', fmt, '-o', str(output)]
    result = run_macs3([*args, '-i', good], timeout=30)
    assert result.returncode == 0, result.stderr
    assert output.read_text().strip()
    expected = ('chr1\t0\t10\t0.00000\nchr1\t10\t50\t10.00000' if fmt == 'BAMPE'
                else 'chr1\t0\t10\t0.00000\nchr1\t10\t210\t10.00000')
    assert output.read_text().strip() == expected
    output.unlink()
    result = run_macs3([*args, '-i', bad], timeout=30)
    assert result.returncode != 0
    assert 'EOFError: Truncated BAM' in result.stderr
    assert not output.exists()
