"""Bounded public BEDPE/FRAG coordinate regressions for #815 / #757 item 18."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import macs3_env
from test_Parser_frag_eof import write_frag


SCRIPT = r'''
import json
from pathlib import Path
import sys
import MACS3.IO.Parser as parser_module
from MACS3.IO.Parser import BEDPEParser, FragParser
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII

fmt, operation, threads, manifest = sys.argv[1:]
cls = BEDPEParser if fmt == 'BEDPE' else FragParser
track_cls = PETrackI if fmt == 'BEDPE' else PETrackII
results = []
for path in json.loads(Path(manifest).read_text()):
    parser = cls(path, buffer_size=2, threads=int(threads))
    try:
        if operation == 'build':
            track = parser.build_petrack()
        else:
            track = track_cls(buffer_size=2)
            if fmt == 'BEDPE':
                track.add_loc(b'chr0', 1, 11)
            else:
                track.add_loc(b'chr0', 1, 11, b'OLD', 1)
            returned = parser.append_petrack(track)
            if returned is not track:
                raise RuntimeError('append did not return the supplied track')
        track.finalize()
        records = {chrom.decode(): track.get_locations_by_chr(chrom).tolist()
                   for chrom in track.get_chr_names()}
        barcodes = ({bc.decode(): i for bc, i in track.barcode_dict.items()}
                    if fmt == 'FRAG' else None)
        results.append(dict(n=parser.n, d=parser.d, records=records,
                            barcodes=barcodes))
    except Exception as exc:
        results.append(dict(error=type(exc).__name__, message=str(exc)))
    finally:
        parser.close()
print(json.dumps(dict(extension=parser_module.__file__, results=results)))
'''


def probe(tmp_path, fmt, operation, optimized, encoding, rows, threads=1):
    paths = [str(write_frag(tmp_path / f'input-{i}', data, encoding))
             for i, data in enumerate(rows)]
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(paths))
    script = tmp_path / 'probe.py'
    script.write_text(SCRIPT)
    result = subprocess.run(
        [sys.executable, *(['-O'] if optimized else []), str(script), fmt,
         operation, str(threads), str(manifest)], cwd=tmp_path,
        env=macs3_env(tmp_path), capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    output = json.loads(result.stdout)
    import MACS3.IO.Parser as parser_module
    assert Path(output['extension']).resolve() == Path(parser_module.__file__).resolve()
    assert Path(output['extension']).suffix != '.py', 'a compiled build is required'
    assert len(output['results']) == len(paths)
    return output['results']


# No assertions live inside the probe: -O must exercise the same checks.
@pytest.mark.parametrize('fmt', ['BEDPE', 'FRAG'])
@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('optimized', [False, True], ids=['normal', 'optimized'])
@pytest.mark.parametrize('encoding', ['plain', 'gzip', 'bgzf'])
def test_fragment_coordinates(tmp_path, fmt, operation, optimized, encoding):
    bad_fields = [b'start', b'', b'10junk', b'1.5', b'1e2', b'+', b'--1',
                  b'1_0', b'\x0010', b'2147483648', b'-2147483649',
                  b'9' * 100]
    cases = []
    for field in ('start', 'end'):
        for value in bad_fields:
            left, right = (value, b'20') if field == 'start' else (b'10', value)
            line = b'chr1\t' + left + b'\t' + right + b'\tBC\t1'
            cases.append((line, field, value))
    cases.extend([(b'chrom\tstart\tend\tbarcode\tcount', 'start', b'start'),
                  (b'chr1\t10\t10\tBC\t1', 'interval', None),
                  (b'chr1\t20\t10\tBC\t1', 'interval', None)])
    inputs = []
    for line, _, _ in cases:
        # First row and later row; complete line and EOF tail use distinct paths.
        for prefix in (b'', b'chr1\t1\t5\tBC\t1\n'):
            for ending in (b'\n', b''):
                inputs.append(prefix + line + ending)
    actual = probe(tmp_path, fmt, operation, optimized, encoding, inputs)
    for (line, field, value), group in zip(cases, range(0, len(actual), 4)):
        for result in actual[group:group + 4]:
            assert result.get('error') == 'ValueError', result
            if field == 'interval':
                assert 'Right position must be larger' in result['message']
            else:
                assert f'Invalid {field} coordinate {value!r}' in result['message']
            assert repr(line) in result['message']


@pytest.mark.parametrize('fmt', ['BEDPE', 'FRAG'])
@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('optimized', [False, True], ids=['normal', 'optimized'])
def test_fragment_coordinate_controls(tmp_path, fmt, operation, optimized):
    valid = (b'#chrom\tstart\tend\tbarcode\tcount\n'
             b'chr1\t  +00010 \t00020\tBC\t2\n'
             b'chr1\t2147483600\t2147483647\tBC\t3\n'
             b'chr1\t-2147483648\t20\tSKIP\t1\n'
             b'\t10\t20\tSKIP\t1\n')
    data = [valid, valid.rstrip(b'\n'), b'', b'#comment', b'#comment\n',
            b'chr1\t10', b'chr1\t10\t20\tBC' if fmt == 'FRAG' else b'chr1']
    results = probe(tmp_path, fmt, operation, optimized, 'plain', data)
    expected = {'chr1': [[10, 20], [2147483600, 2147483647]]}
    if fmt == 'FRAG':
        expected['chr1'] = [[10, 20, 2], [2147483600, 2147483647, 3]]
    if operation == 'append':
        expected['chr0'] = [[1, 11, 1]] if fmt == 'FRAG' else [[1, 11]]
    for result in results[:2]:
        assert result['n'] == 2 and result['d'] == 28.5
        assert result['records'] == expected
        if fmt == 'FRAG':
            assert result['barcodes'] == ({'OLD': 0, 'BC': 1} if operation == 'append'
                                          else {'BC': 0})
    for result in results[2:5]:
        assert result == dict(error='ValueError', message=f'No valid fragments found in {fmt} input')
    for result in results[5:]:
        assert result['error'] == 'Exception'
        assert f'Less than {5 if fmt == "FRAG" else 3} columns' in result['message']


@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('optimized', [False, True], ids=['normal', 'optimized'])
@pytest.mark.parametrize('threads', [1, 2])
def test_frag_coordinate_after_large_bgzf_prefix(tmp_path, operation, optimized, threads):
    # >64 KiB triggers the threaded FRAG line walk with two threads. Include
    # a malformed row before more valid rows as well as an EOF-tail variant.
    prefix = b'chr1\t10\t20\tBC\t1\n' * 10000
    bad = b'chr1\t10junk\t20\tBC\t1'
    inputs = [prefix + bad + b'\n' + prefix, prefix + bad, prefix]
    results = probe(tmp_path, 'FRAG', operation, optimized, 'bgzf', inputs, threads)
    for result in results[:2]:
        assert result.get('error') == 'ValueError', result
        assert "Invalid start coordinate b'10junk'" in result['message']
        assert repr(bad) in result['message']
    assert results[2]['n'] == 10000
    assert results[2]['d'] == 10.0
    assert len(results[2]['records']['chr1']) == 10000
