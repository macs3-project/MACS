"""Bounded BED/BEDPE EOF regressions for #804 / #757 item 7."""

import gzip
import subprocess
import sys

import pytest

from conftest import macs3_env
from MACS3.Utilities.Constants import READ_BUFFER_SIZE


SCRIPT = r'''
import sys
from MACS3.IO.Parser import BEDParser, BEDPEParser, StrandFormatError
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI

fmt, operation, mode = sys.argv[1:4]
paths = sys.argv[4:]
parser_class, track_class, suffix = (
    (BEDParser, FWTrack, 'fwtrack') if fmt == 'BED' else
    (BEDPEParser, PETrackI, 'petrack'))


def snapshot(track):
    track.finalize()
    result = {}
    for chrom in track.get_chr_names():
        locs = track.get_locations_by_chr(chrom)
        result[chrom] = ([x.tolist() for x in locs] if fmt == 'BED' else
                         locs.tolist())
    return result


def load(path):
    parser = parser_class(path, buffer_size=2)
    if operation == 'build':
        track = getattr(parser, 'build_' + suffix)()
    else:
        track = track_class(buffer_size=2)
        if fmt == 'BED':
            track.add_loc(b'chr0', 7, 0)
        else:
            track.add_loc(b'chr0', 7, 17)
        assert getattr(parser, 'append_' + suffix)(track) is track
    return parser, track


if mode in ('valid', 'single', 'skipped', 'boundary'):
    results = []
    for path in paths:
        parser, track = load(path)
        stats = None if fmt == 'BED' else (parser.n, parser.d)
        results.append((snapshot(track), stats))
    assert results[0] == results[1], results
    expected = {}
    if operation == 'append':
        expected[b'chr0'] = [[7], []] if fmt == 'BED' else [(7, 17)]
    if mode != 'single':
        expected[b'chr1'] = [[10], []] if fmt == 'BED' else [(10, 50)]
    if mode != 'skipped':
        expected[b'chr2'] = [[], [130]] if fmt == 'BED' else [(60, 130)]
    assert results[0][0] == expected, results[0]
    if fmt == 'BEDPE':
        assert results[0][1] == ((1, 70.0) if mode == 'single' else
                                 (1, 40.0) if mode == 'skipped' else
                                 (2, 55.0)), results[0]

elif mode in ('malformed', 'invalid'):
    errors = []
    for path in paths:
        try:
            load(path)
        except (Exception, StrandFormatError) as exc:
            errors.append((type(exc).__name__, str(exc)))
        else:
            raise AssertionError('malformed record accepted')
    assert errors[0] == errors[1], errors
    expected = (('IndexError' if fmt == 'BED' else 'Exception')
                if mode == 'malformed' else
                ('StrandFormatError' if fmt == 'BED' else 'AssertionError'))
    assert errors[0][0] == expected, errors
    if mode == 'malformed' and fmt == 'BEDPE':
        assert 'broken' in errors[0][1], errors

elif mode == 'empty':
    for path in paths:
        if fmt == 'BEDPE':
            try:
                load(path)
            except ValueError as exc:
                assert str(exc) == 'No valid fragments found in BEDPE input'
            else:
                raise AssertionError('empty BEDPE input accepted')
        else:
            parser, track = load(path)
            assert snapshot(track) == (
                {b'chr0': [[7], []]} if operation == 'append' else {})

elif mode == 'history':
    parser, track = load(paths[0])
    before = (parser.n, parser.d, snapshot(track))
    for path in paths[1:]:
        parser.__init__(path, buffer_size=2)
        assert parser.append_petrack(track) is track
        assert (parser.n, parser.d, snapshot(track)) == before
'''


def write_input(path, data, encoding):
    path.write_bytes(gzip.compress(data) if encoding == "gzip" else data)
    return path


def probe(tmp_path, fmt, operation, mode, paths):
    script = tmp_path / "probe.py"
    script.write_text(SCRIPT)
    result = subprocess.run(
        [sys.executable, str(script), fmt, operation, mode, *map(str, paths)],
        cwd=tmp_path, env=macs3_env(), capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def records(fmt):
    if fmt == "BED":
        return b'chr1\t10\t50\ta\t0\t+', b'chr2\t60\t130\tb\t0\t-'
    return b'chr1\t10\t50', b'chr2\t60\t130'


@pytest.mark.parametrize("fmt", ["BED", "BEDPE"])
@pytest.mark.parametrize("operation", ["build", "append"])
@pytest.mark.parametrize("encoding", ["plain", "gzip"])
@pytest.mark.parametrize("mode", ["valid", "single", "skipped",
                                  "malformed", "invalid"])
def test_text_final_record(tmp_path, fmt, operation, encoding, mode):
    first, last = records(fmt)
    if mode == "single":
        data = last
    else:
        if mode == "malformed":
            last = b'broken'
        elif mode == "invalid":
            last = (b'chr2\t60\t130\tb\t0\t?' if fmt == "BED" else
                    b'chr2\t130\t60')
        elif mode == "skipped":
            last = b'chr2\t-1\t130'
        headers = b'# sample\ntrack name=sample\nbrowser position chr1\n'
        data = headers + first + b'\n' + last
    paths = [write_input(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, fmt, operation, mode, paths)


@pytest.mark.parametrize("fmt", ["BED", "BEDPE"])
@pytest.mark.parametrize("operation", ["build", "append"])
@pytest.mark.parametrize("encoding", ["plain", "gzip"])
def test_text_no_accepted_records(tmp_path, fmt, operation, encoding):
    inputs = [b'', b'# header', b'# header\n',
              b'track name=x\nbrowser position chr1',
              b'track name=x\nbrowser position chr1\n',
              b'chr1\t-1\t40', b'chr1\t-1\t40\n',
              b'\t0\t40', b'\t0\t40\n']
    paths = [write_input(tmp_path / str(i), data, encoding)
             for i, data in enumerate(inputs)]
    probe(tmp_path, fmt, operation, "empty", paths)


@pytest.mark.parametrize("fmt", ["BED", "BEDPE"])
@pytest.mark.parametrize("operation", ["build", "append"])
@pytest.mark.parametrize("encoding", ["plain", "gzip"])
def test_text_final_record_across_blocks(tmp_path, fmt, operation, encoding):
    # The final record exceeds the 10 MB read buffer and spans refills at EOF.
    # The name/extra column is ignored by the public parser. ~0.4 s per probe.
    last = b'chr2\t60\t130\t' + b'X' * READ_BUFFER_SIZE
    if fmt == "BED":
        last += b'\t0\t-'
    data = records(fmt)[0] + b'\n' + last
    paths = [write_input(tmp_path / str(i), data + ending, encoding)
             for i, ending in enumerate([b'\n', b''])]
    probe(tmp_path, fmt, operation, "boundary", paths)


@pytest.mark.parametrize("encoding", ["plain", "gzip"])
def test_bedpe_empty_append_preserves_parser_history(tmp_path, encoding):
    inputs = [b'chr1\t10\t50', b'', b'# header', b'# header\n',
              b'chr2\t-1\t40', b'chr2\t-1\t40\n']
    paths = [write_input(tmp_path / str(i), data, encoding)
             for i, data in enumerate(inputs)]
    probe(tmp_path, "BEDPE", "build", "history", paths)


@pytest.mark.parametrize("fmt", ["BED", "BEDPE"])
@pytest.mark.parametrize("encoding", ["plain", "gzip"])
def test_text_pileup_final_record(tmp_path, run_macs3, fmt, encoding):
    data = b'\n'.join(records(fmt))
    outputs = []
    for i, ending in enumerate([b'\n', b'']):
        path = write_input(tmp_path / str(i), data + ending, encoding)
        output = tmp_path / (str(i) + '.bdg')
        result = run_macs3(['pileup', '-f', fmt, '-i', path, '-o', output],
                          timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append(output.read_text())
    assert outputs[0] == outputs[1]
    assert outputs[0].strip()
