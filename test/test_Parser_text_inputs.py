"""Public text-parser and bounded CLI regressions for #811 / #757 item 14."""

import gzip
import subprocess
import sys

import pytest

from conftest import macs3_env


HEADER = b'@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n@PG\tID:bwa\tPN:bwa\n'
SAM = (b'r1\t0\tchr1\t101\t60\t25M\t*\t0\t0\t' + b'A' * 25 +
       b'\t' + b'I' * 25 + b'\n' +
       b'r2\t16\tchr1\t201\t60\t25M\t*\t0\t0\t' + b'A' * 25 +
       b'\t' + b'I' * 25 + b'\n')
BED = b'chr1\t100\t125\tr1\t0\t+\nchr1\t200\t225\tr2\t0\t-\n'
SHORT_BED = b'chr1\t100\t101\tr1\t0\t+\nchr1\t200\t201\tr2\t0\t-\n'


SCRIPT = r'''
import sys
from MACS3.IO.Parser import (
    BEDParser, SAMParser, BowtieParser, ELANDResultParser,
    ELANDMultiParser, ELANDExportParser, guess_parser)
from MACS3.Signal.FixWidthTrack import FWTrack

mode, fmt, operation, *paths = sys.argv[1:]
classes = dict(BED=BEDParser, SAM=SAMParser, BOWTIE=BowtieParser,
               ELAND=ELANDResultParser, ELANDMULTI=ELANDMultiParser,
               ELANDEXPORT=ELANDExportParser)
for path in paths:
    if mode == 'unknown':
        try:
            guess_parser(path)
        except ValueError as exc:
            assert str(exc) == "Can't detect format!", str(exc)
        else:
            raise AssertionError('unrecognised input accepted')
        continue
    parser = guess_parser(path) if mode == 'auto' else classes[fmt](path)
    if mode == 'auto':
        assert type(parser) is classes[fmt]
    if mode == 'empty':
        assert parser.tsize() == -1
        assert parser.sniff() is False
    elif mode == 'sniff-miss':
        assert parser.sniff() is False
        # A failed sniff must rewind even when field parsing raises.
        # BED coordinates below remain usable with an explicit format.
    else:
        assert parser.tsize() == (1 if mode == 'short' else 25)
        if mode != 'short':
            assert parser.sniff() is True
    if operation == 'append':
        track = FWTrack(buffer_size=2)
        track.add_loc(b'chr0', 7, 0)
        assert parser.append_fwtrack(track) is track
    else:
        track = parser.build_fwtrack()
    track.finalize()
    actual = {chrom: [a.tolist() for a in track.get_locations_by_chr(chrom)]
              for chrom in track.get_chr_names()}
    expected = {b'chr0': [[7], []]} if operation == 'append' else {}
    if mode != 'empty':
        expected[b'chr1'] = [[100], [201 if mode in ('short', 'sniff-miss') else 225]]
    assert actual == expected, (actual, expected)
'''


def write_input(path, data, encoding):
    path.write_bytes(gzip.compress(data) if encoding == 'gzip' else data)
    return path


def probe(tmp_path, mode, fmt, operation, paths):
    script = tmp_path / 'probe.py'
    script.write_text(SCRIPT)
    result = subprocess.run(
        [sys.executable, str(script), mode, fmt, operation, *map(str, paths)],
        cwd=tmp_path, env=macs3_env(), capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('fmt', ['BED', 'SAM'])
def test_blank_and_incomplete_rows(tmp_path, encoding, operation, fmt):
    # Include each truncated SAM field count, including 10 fields (SEQ but
    # missing QUAL), and an embedded header after alignments. Compare to
    # clean input to pin valid coordinates, tag sampling and stream rewind.
    bad_sam = b'\n'.join(b'\t'.join(SAM.splitlines()[0].split(b'\t')[:n])
                         for n in range(1, 11)) + b'\n'
    clean = BED if fmt == 'BED' else HEADER + SAM
    dirty = (b'\n \t\n' + BED + b'\n\t\n' if fmt == 'BED' else
             HEADER + bad_sam + SAM + b'@CO\tcomment\n\n')
    paths = [write_input(tmp_path / str(i), data, encoding)
             for i, data in enumerate([clean, dirty, dirty.rstrip(b'\n')])]
    probe(tmp_path, 'valid', fmt, operation, paths)


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
@pytest.mark.parametrize('operation', ['build', 'append'])
@pytest.mark.parametrize('fmt', ['BED', 'SAM', 'BOWTIE', 'ELAND',
                                 'ELANDMULTI', 'ELANDEXPORT'])
def test_empty_and_header_only(tmp_path, encoding, operation, fmt):
    header = HEADER if fmt == 'SAM' else (
        b'track name=x\nbrowser position chr1\n# comment\n' if fmt == 'BED'
        else b'# comment\n')
    inputs = [b'', b'\n \t\n', header, header.rstrip(b'\n')]
    # Short SAM rows are skipped even when no usable alignment follows.
    if fmt == 'SAM':
        inputs += [b'r1\t0\tchr1\t101\n', header + b'broken']
    paths = [write_input(tmp_path / str(i), data, encoding)
             for i, data in enumerate(inputs)]
    probe(tmp_path, 'empty', fmt, operation, paths)


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
def test_auto_rejects_unknown_inputs(tmp_path, encoding):
    inputs = [b'', b'\n \t\n', b'broken\n', HEADER,
              b'# comment\n', SHORT_BED, b'chr1\t100\n']
    paths = [write_input(tmp_path / str(i), data, encoding)
             for i, data in enumerate(inputs)]
    probe(tmp_path, 'unknown', 'AUTO', 'build', paths)


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
@pytest.mark.parametrize('fmt', ['BED', 'SAM'])
def test_auto_preserves_valid_reads(tmp_path, encoding, fmt):
    path = write_input(tmp_path / 'input', BED if fmt == 'BED' else HEADER + SAM,
                       encoding)
    probe(tmp_path, 'auto', fmt, 'build', [path])


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_short_bed_explicit_format(tmp_path, encoding, operation):
    path = write_input(tmp_path / 'input', SHORT_BED, encoding)
    probe(tmp_path, 'short', 'BED', operation, [path])


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
def test_failed_sniff_rewinds(tmp_path, encoding):
    # The third column is absent, raising IndexError during BED sampling;
    # explicit BED loading still supports a plus-strand start-only record.
    path = write_input(tmp_path / 'input',
                       b'chr1\t100\nchr1\t200\t201\tx\t0\t-\n', encoding)
    probe(tmp_path, 'sniff-miss', 'BED', 'build', [path])


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
@pytest.mark.parametrize('fmt,data', [('BED', SHORT_BED), ('BED', BED + b'\n'),
                                     ('SAM', HEADER + SAM + b'@CO\tx\n')])
def test_filterdup_text_inputs(tmp_path, run_macs3, encoding, fmt, data):
    path = write_input(tmp_path / 'input', data, encoding)
    output = tmp_path / 'out.bed'
    result = run_macs3(['filterdup', '-i', path, '-f', fmt, '-g', '1e6',
                        '--keep-dup', 'all', '-o', output], timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = [line.split('\t') for line in output.read_text().splitlines()]
    assert len(rows) == 2
    assert sorted((r[0], int(r[1]), int(r[2]), r[5]) for r in rows) == [
        ('chr1', 100, 101 if data == SHORT_BED else 125, '+'),
        ('chr1', 200, 201 if data == SHORT_BED else 225, '-')]


@pytest.mark.parametrize('encoding', ['plain', 'gzip'])
@pytest.mark.parametrize('fmt,data', [('BED', b''), ('BED', b'\n \t\n'),
                                     ('SAM', HEADER), ('SAM', b'r1\t0\tchr1\n'),
                                     ('AUTO', b''), ('AUTO', SHORT_BED)])
def test_cli_rejects_unusable_inputs(tmp_path, run_macs3, encoding, fmt, data):
    path = write_input(tmp_path / 'input', data, encoding)
    for command in ('callpeak', 'filterdup'):
        args = (['callpeak', '-t', path, '-n', 'probe'] if command == 'callpeak'
                else ['filterdup', '-i', path, '-o', tmp_path / 'out.bed'])
        result = run_macs3([*args, '-f', fmt, '-g', '1e6'], timeout=30)
        assert result.returncode == 1, result.stdout + result.stderr
        expected = ("Can't detect format!" if fmt == 'AUTO' else
                    'No usable alignments remain')
        assert expected in result.stderr
        assert 'Traceback' not in result.stderr
