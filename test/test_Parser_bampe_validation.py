"""Public BAMPE regressions for #806, including bounded native probes."""

import json
import struct
import subprocess
import sys

import pytest

from conftest import macs3_env
from test_Parser_bam_truncation import write_bam


HEADER = (b'BAM\1' + struct.pack('<ii', 0, 2)
          + struct.pack('<i', 5) + b'chr1\0' + struct.pack('<i', 1000000)
          + struct.pack('<i', 5) + b'chr2\0' + struct.pack('<i', 1000000))


def alignment(tlen=200, ref=0, mate_ref=0, pos=100, mate_pos=250,
              flag=99):
    """One 50M alignment; BAM core fields are little-endian int32."""
    body = (struct.pack('<iiIIiiii', ref, pos, (60 << 8) | 2,
                        (flag << 16) | 1, 50, mate_ref, mate_pos, tlen)
            + b'r\0' + struct.pack('<I', 50 << 4)
            + b'\x11' * 25 + b'\x28' * 50)
    return struct.pack('<i', len(body)) + body


INVALID = {
    'zero_tlen': dict(tlen=0),
    'other_chromosome': dict(mate_ref=1),
    'missing_mate_reference': dict(mate_ref=-1),
    'out_of_range_mate_reference': dict(mate_ref=2),
    'other_chromosome_zero_tlen': dict(mate_ref=1, tlen=0),
    'missing_reference': dict(ref=-1, mate_ref=-1),
    'negative_reference': dict(ref=-2, mate_ref=-2),
    'out_of_range_reference': dict(ref=2, mate_ref=2),
    'int_min_tlen': dict(tlen=-(2**31)),
    'overflow_end': dict(tlen=2**31 - 1),
    'negative_position': dict(pos=-1),
    'missing_mate_position': dict(mate_pos=-1),
}


PROBE = '''
import json
import math
import sys
from MACS3.IO.Parser import BAMPEParser
from MACS3.Signal.PairedEndTrack import PETrackI

operation, manifest = sys.argv[1:]
for path, expected in json.loads(open(manifest).read()):
    parser = BAMPEParser(path, buffer_size=2)
    track = PETrackI(buffer_size=2)
    if operation == 'append':
        track.add_loc(b'chr2', 10, 60)
    try:
        result = (parser.append_petrack(track) if operation == 'append'
                  else parser.build_petrack())
    except ValueError as exc:
        assert not expected, (path, str(exc))
        assert 'No valid fragments found in BAMPE input' in str(exc)
        assert parser.n == 0 and parser.d == 0
        if operation == 'append':
            track.finalize()
            assert track.total == 1
            assert track.get_locations_by_chr(b'chr2').tolist() == [(10, 60)]
        continue
    assert expected, 'invalid input accepted: ' + path
    if operation == 'append':
        assert result is track
    assert parser.n == len(expected), (path, parser.n)
    mean = sum(r - l for l, r in expected) / len(expected) if expected else 0
    assert math.isclose(parser.d, mean, rel_tol=1e-6), (path, parser.d, mean)
    result.finalize()
    assert result.total == len(expected) + (operation == 'append')
    observed = (result.get_locations_by_chr(b'chr1').tolist()
                if b'chr1' in result.get_chr_names() else [])
    assert observed == [tuple(pair) for pair in expected], (path, observed)
    if operation == 'append':
        assert result.get_locations_by_chr(b'chr2').tolist() == [(10, 60)]
'''


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
@pytest.mark.parametrize('operation', ['build', 'append'])
def test_bampe_fragment_validation(tmp_path, encoding, operation):
    valid = alignment() + alignment(tlen=-300, pos=450, mate_pos=200)
    cases = [('valid', valid, [(100, 300), (200, 500)]),
             ('max_endpoint', alignment(tlen=2**31 - 1, pos=0, mate_pos=0),
              [(0, 2**31 - 1)]),
             ('negative_max_tlen', alignment(tlen=-(2**31 - 1), pos=0,
                                             mate_pos=0), [(0, 2**31 - 1)]),
             ('second_mate', valid + alignment(flag=147),
              [(100, 300), (200, 500)]),
             ('empty', b'', []),
             ('filtered', alignment(flag=4), [])]
    for name, fields in INVALID.items():
        invalid = alignment(**fields)
        cases.extend([(name, invalid, []),
                      (name + '_mixed', invalid + valid,
                       [(100, 300), (200, 500)])])
    cases.append(('all_invalid', b''.join(alignment(**f)
                                         for f in INVALID.values()), []))
    manifest = []
    for name, data, expected in cases:
        path = write_bam(tmp_path / (name + '.bam'), HEADER + data, encoding)
        manifest.append((path, expected))
    script = tmp_path / 'probe.py'
    script.write_text(PROBE)
    manifest_path = tmp_path / 'cases.json'
    manifest_path.write_text(json.dumps(manifest))
    result = subprocess.run(
        [sys.executable, str(script), operation, str(manifest_path)],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True,
        text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('sample', ['empty', 'zero_tlen', 'other_chromosome'])
@pytest.mark.parametrize('role', ['treatment', 'control'])
def test_callpeak_rejects_no_usable_bampe_fragments(tmp_path, run_macs3,
                                                 sample, role):
    data = b'' if sample == 'empty' else alignment(**INVALID[sample]) * 20
    bad = write_bam(tmp_path / 'bad.bam', HEADER + data, 'bgzf')
    good = write_bam(tmp_path / 'good.bam', HEADER + alignment() * 20, 'bgzf')
    args = ['callpeak', '-f', 'BAMPE', '-g', '1e6', '-n', 'x',
            '-t', bad if role == 'treatment' else good]
    if role == 'control':
        args += ['-c', bad]
    result = run_macs3(args, timeout=30)
    assert result.returncode != 0
    assert 'No valid fragments found in BAMPE input' in result.stderr
    assert 'bad.bam' in result.stderr
    assert 'Traceback' not in result.stderr
    assert 'Lambda must' not in result.stderr


@pytest.mark.parametrize('with_control', [False, True])
def test_callpeak_bampe_mixed_matches_valid(tmp_path, run_macs3, with_control):
    valid = alignment() * 20
    good = write_bam(tmp_path / 'good.bam', HEADER + valid, 'bgzf')
    mixed = write_bam(tmp_path / 'mixed.bam', HEADER + valid +
                      b''.join(alignment(**f) for f in INVALID.values()), 'bgzf')
    for name, path in [('good', good), ('mixed', mixed)]:
        args = ['callpeak', '-f', 'BAMPE', '-g', '1e6', '--keep-dup', 'all',
                '-B', '-n', name, '-t', path, '--outdir', str(tmp_path)]
        if with_control:
            args += ['-c', good]
        result = run_macs3(args, timeout=30)
        assert result.returncode == 0, result.stderr
        assert 'fragment size = 200.0' in result.stderr
    for suffix in ['_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg',
                   '_control_lambda.bdg']:
        # Names are the only expected difference in the output records.
        expected = (tmp_path / ('good' + suffix)).read_text()
        observed = (tmp_path / ('mixed' + suffix)).read_text()
        assert observed.replace('mixed', 'good') == expected
