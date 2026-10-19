"""Compressed framing/fallback matrix; each bounded child takes about 0.5 s."""

import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import random
import signal
import subprocess
import sys

import pytest

from conftest import macs3_env
from MACS3.Utilities.Constants import READ_BUFFER_SIZE
from test_Parser_bam_truncation import HEADER, byte_array, record, write_bam
from test_Parser_frag_eof import write_frag


def native_probe(tmp_path, scenario, *args, extra_env=None):
    script = tmp_path / 'native_probe.py'
    script.write_bytes(Path(__file__).with_name('native_probe.py').read_bytes())
    (tmp_path / 'execution_probe.py').write_bytes(
        Path(__file__).with_name('execution_probe.py').read_bytes())
    temp = tmp_path / 'tmp'
    temp.mkdir(exist_ok=True)
    proc = subprocess.Popen([sys.executable, str(script), scenario, *map(str, args)],
                            cwd=tmp_path, env=macs3_env(temp, extra_env),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, start_new_session=True)
    try:
        out, err = proc.communicate(timeout=45)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        out, err = proc.communicate()
        pytest.fail('native probe timed out\n' + out + err)
    assert proc.returncode == 0, out + err
    return json.loads(out)


def encoded(tmp_path, kind, data, encoding, name):
    path = tmp_path / name
    if kind == 'frag':
        write_frag(path, data, encoding)
    else:
        write_bam(path, data, encoding)
    return path


@pytest.mark.parametrize('kind', ['bam', 'bampe', 'frag'])
@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('damage', ['crc', 'size', 'deflate', 'truncated'])
@pytest.mark.parametrize('native', [False, pytest.param(True, marks=pytest.mark.slow)])
def test_bad_compressed_member_across_budgets(tmp_path, kind, operation, damage, native):
    prefix = (b'chr1\t10\t50\tA\t1\n' * 5000 if kind == 'frag' else
              HEADER + record() * 5000)
    suffix = b'chr1\t20\t60\tB\t2\n' * 5000 if kind == 'frag' else record() * 5000
    if native:
        # Past the Python header reader's 10 MB lookahead, so corruption is
        # reached by the native stream itself (~1 s, <100 MiB per case).
        if kind == 'frag':
            row = b'chr1\t10\t50\tA\t1\n'
            prefix = row * ((READ_BUFFER_SIZE + 65536) // len(row) + 1)
        else:
            prefix = HEADER + record(aux=byte_array(b'\0' * (READ_BUFFER_SIZE + 65536)))
    first = encoded(tmp_path, kind, prefix, 'bgzf', 'first').read_bytes()[:-28]
    last = encoded(tmp_path, kind, suffix, 'bgzf', 'last').read_bytes()
    # Corrupt the first member of the suffix after a healthy multi-block prefix.
    length = int.from_bytes(last[16:18], 'little') + 1
    bad = bytearray(last)
    if damage == 'truncated':
        bad = bad[:length - 4]
    elif damage == 'deflate':
        bad[18] = (bad[18] & ~6) | 6  # reserved DEFLATE block type
    else:
        bad[length - (8 if damage == 'crc' else 4)] ^= 1
    path = tmp_path / 'bad.bgzf'
    path.write_bytes(first + bad)
    manifest = tmp_path / 'cases.json'
    manifest.write_text(json.dumps([[str(path), kind, operation, budget, 0]
                                    for budget in (1, 4)]))
    results = native_probe(tmp_path, 'streams', manifest)
    expected = ('EOFError' if damage == 'truncated' else
                'error' if damage == 'deflate' and not native else 'BadGzipFile')
    assert [r[0] for r in results] == [expected, expected], results
    assert all(r[1] for r in results)


@pytest.mark.parametrize('kind', ['bam', 'bampe', 'frag'])
@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('native', [False, pytest.param(True, marks=pytest.mark.slow)])
def test_mixed_members_and_padding_fall_back_without_data_loss(tmp_path, kind, operation, native):
    row = b'chr1\t10\t50\tA\t2\n' if kind == 'frag' else record()
    prefix = b'' if kind == 'frag' else HEADER
    prefix += row * 5000
    count = 15000
    if native:
        # Exceed compressed batch input and decompressed header lookahead.
        # A staged second input and a large record survive the serial handoff.
        # Generated only for slow tests; ~2–4 s and <350 MiB per child.
        noise = random.Random(771).randbytes(5 * 1024 * 1024)
        if kind == 'frag':
            prefix += row * ((READ_BUFFER_SIZE + 65536) // len(row) + 1)
            count += (READ_BUFFER_SIZE + 65536) // len(row) + 1
            prefix += b'chr2\t20\t60\t' + noise.hex().encode() + b'\t2\n'
        else:
            prefix += record(aux=byte_array(noise + b'\0' * (12 * 1024 * 1024)))
        count += 1
    reference = encoded(tmp_path, kind, prefix + row * 10000, 'bgzf', 'reference')
    first = encoded(tmp_path, kind, prefix, 'bgzf', 'first').read_bytes()[:-28]
    if native:
        assert len(first) > 4 * 1024 * 1024
    tail = encoded(tmp_path, kind, row * 5000, 'bgzf', 'tail').read_bytes()
    mixed = tmp_path / 'mixed'
    mixed.write_bytes(first + b'\0' * 17 + gzip.compress(row * 5000) + tail)
    manifest = tmp_path / 'cases.json'
    manifest.write_text(json.dumps([[str(path), kind, operation, budget, 0]
                                    for path in (reference, mixed) for budget in (1, 4)]))
    results = native_probe(tmp_path, 'streams', manifest)
    assert results[0][0] == 'ok', results
    assert results[0][1] == count * (2 if kind == 'frag' else 1)
    assert all(result == results[0] for result in results), results


@pytest.mark.parametrize('kind', ['bam', 'bampe', 'frag'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_native_stream_failure_releases_workers_before_next_load(tmp_path, kind, operation):
    # Beyond header lookahead so failure occurs in the compiled stream.
    # ~1–2 s per bounded child. Repeated failures must not retain worker pools.
    if kind == 'frag':
        row = b'chr1\t10\t50\tA\t2\n'
        prefix = row * ((READ_BUFFER_SIZE + 65536) // len(row) + 1)
    else:
        row = record()
        prefix = HEADER + record(aux=byte_array(b'\0' * (READ_BUFFER_SIZE + 65536)))
    first = encoded(tmp_path, kind, prefix, 'bgzf', 'prefix').read_bytes()[:-28]
    good = encoded(tmp_path, kind, prefix + row, 'bgzf', 'good')
    damaged = bytearray(encoded(tmp_path, kind, row, 'bgzf', 'last').read_bytes())
    size = int.from_bytes(damaged[16:18], 'little') + 1
    damaged[size - 8] ^= 1
    bad = tmp_path / 'bad'
    bad.write_bytes(first + damaged)
    cases = [[str(path), kind, operation, budget, 0]
             for _ in range(3) for budget in (1, 4) for path in (good, bad, good)]
    manifest = tmp_path / 'cases.json'
    manifest.write_text(json.dumps(cases))
    results = native_probe(tmp_path, 'streams', manifest, 'cleanup')
    for index in range(0, len(results), 3):
        assert results[index] == results[index + 2] == results[0]
        assert results[index + 1][0] == 'BadGzipFile', results


@pytest.mark.parametrize('encoding', ['plain', 'gzip', 'bgzf'])
@pytest.mark.parametrize('cap', [0, 2, 65535])
@pytest.mark.slow
def test_frag_long_unusual_barcodes_and_count_limits(tmp_path, encoding, cap):
    # 4 MiB line crosses the serial buffer and many BGZF blocks; the final
    # record has no newline. Under 1 s per child, ~20 MiB peak fixture memory.
    bars = [b'A\0B', b'\xff\xfe', '条碼'.encode(), b'X' * ((1 << 22) + 17),
            b'OVER', b'ZERO']
    rows = [(b'chr1', i * 50, i * 50 + 40, bar, count)
            for i, (bar, count) in enumerate(zip(bars, (1, 2, 65535, 65535, 65536, 0)))]
    data = b'\n'.join(b'%s\t%d\t%d\t%s\t%d' % row for row in rows)
    path = encoded(tmp_path, 'frag', data, encoding, 'long')
    manifest = tmp_path / 'cases.json'
    manifest.write_text(json.dumps([[str(path), 'frag', operation, budget, cap]
                                    for operation in ('build', 'append', 'fallback')
                                    for budget in (1, 4)]))
    results = native_probe(tmp_path, 'streams', manifest)
    assert results[0][0] == 'ok', results
    # Compatibility with pe_parse_line's existing C unsigned-short conversion:
    # conversion precedes max_count, so 65536 becomes zero. The scalar public
    # append receiver above independently checks this legacy behavior.
    counts = [min(row[4] & 65535, cap or 65535) for row in rows]
    assert results[0][1] == sum(counts)
    assert results[0][3] == [[hashlib.sha256(bar).hexdigest(), code]
                            for code, bar in enumerate(bars)]
    assert results[0][4] == {b'chr1'.hex(): [[row[1], row[2], count]
                                          for row, count in zip(rows, counts)]}
    assert all(result == results[0] for result in results), results


@pytest.mark.skipif('fork' not in multiprocessing.get_all_start_methods(), reason='requires fork')
@pytest.mark.parametrize('eligibility', ['host', 'controlled'])
@pytest.mark.parametrize('phase', ['histogram', 'calling'])
@pytest.mark.parametrize('fail_at', [1, 2])
def test_partial_process_startup_joins_and_cleans(tmp_path, eligibility, phase, fail_at):
    if eligibility == 'host' and sys.platform != 'linux':
        pytest.skip('production fork backend requires Linux')
    assert native_probe(tmp_path, 'process_start_failure', eligibility, phase, fail_at) == {
        'started': fail_at - 1, 'joined': True}
