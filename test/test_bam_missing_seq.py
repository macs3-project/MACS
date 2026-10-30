"""Public API/CLI regressions for #810 (part of #757)."""

import json
from pathlib import Path
import struct
import subprocess
import sys

import pytest

from conftest import macs3_env, no_core_dumps
from test_Parser_bam_truncation import write_bam


HEADER = (b'BAM\1' + struct.pack('<ii', 0, 1) + struct.pack('<i', 5)
          + b'chr1\0' + struct.pack('<i', 1000000))


def record(cigar=((50, 0),), seq_length=0, pos=1000, flag=0):
    core = struct.pack('<iiIIiiii', 0, pos, (60 << 8) | 2,
                       (flag << 16) | len(cigar), seq_length, -1, -1, 0)
    body = core + b'r\0' + b''.join(struct.pack('<I', n << 4 | op)
                                  for n, op in cigar)
    if seq_length > 0:
        body += b'\x11' * ((seq_length + 1)//2) + b'\x1e' * seq_length
    return struct.pack('<i', len(body)) + body


def probe(tmp_path, source, args=()):
    from MACS3.IO import Parser
    script = tmp_path / 'probe.py'
    script.write_text('''
from pathlib import Path
import sys
from MACS3.IO import Parser
assert Path(Parser.__file__).resolve() == Path(sys.argv[1]).resolve()
assert Path(Parser.__file__).suffix != '.py'
''' + source)
    result = subprocess.run(
        [sys.executable, str(script), Parser.__file__, *map(str, args)],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True,
        text=True, timeout=30, preexec_fn=no_core_dumps)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
def test_missing_seq_inference(tmp_path, encoding):
    complex_cigar = ((5, 4), (20, 0), (3, 1), (10, 7), (7, 2),
                     (11000, 3), (15, 8), (2, 5), (1, 6))
    cases = [
        ('one', record(), 50),
        ('short', record() * 9, 50),
        ('ten', record() * 10, 50),
        ('complex', record(complex_cigar), 53),
        ('strands', record(complex_cigar) + record(complex_cigar, flag=16), 53),
        ('mixed', record() + record(seq_length=60), 55),
        ('mean', record() + record(((51, 0),)), 50),
        ('first-ten', record()*10 + record(((100, 0),)), 50),
        # Positive l_seq retains precedence even when CIGAR query length differs.
        ('sequence-precedence', record(seq_length=60), 60),
    ]
    manifest = []
    for name, data, expected in cases:
        path = write_bam(tmp_path/(name+'.bam'), HEADER+data, encoding)
        manifest.append((path, expected))
    probe(tmp_path, '''
import json
for path, expected in json.loads(sys.argv[2]):
    for factory in (Parser.BAMParser, Parser.guess_parser):
        p = factory(path)
        try:
            assert p.sniff() is True
            assert p.tsize() == expected, (path, p.tsize(), expected)
            track = p.build_fwtrack()
            track.finalize()
            assert track.total > 0  # inference/sniff rewound the stream
            if Path(path).stem == 'strands':
                plus, minus = track.get_locations_by_chr(b'chr1')
                assert plus.tolist() == [1000]
                assert minus.tolist() == [12052]
            assert p.tsize() == expected  # cached after loader closes stream
        finally:
            p.close()
''', [json.dumps(manifest)])


@pytest.mark.parametrize('encoding', ['raw', 'gzip', 'bgzf'])
def test_unusable_or_malformed_lengths(tmp_path, encoding):
    truncated_cigar = bytearray(record())
    struct.pack_into('<H', truncated_cigar, 4+12, 2)
    missing_name = bytearray(record())
    missing_name[4+8] = 0
    cases = [
        ('no-cigar', record(()), 'Cannot infer BAM tag size'),
        ('no-query', record(((10, 2), (20, 3), (5, 5), (1, 6))),
         'Cannot infer BAM tag size'),
        ('mixed-unusable', record()+record(()), 'Cannot infer BAM tag size'),
        ('zero-op', record(((0, 0),)), 'invalid CIGAR operation'),
        ('unknown-op', record(((50, 15),)), 'invalid CIGAR operation'),
        ('negative-seq', record(seq_length=-1), 'negative l_seq'),
        ('short-cigar', bytes(truncated_cigar), 'CIGAR runs past'),
        ('missing-name', bytes(missing_name), 'read name is missing'),
        ('overflow', record(((2**28-1, 0),)*9), 'exceeds INT_MAX'),
        ('header-only', b'', 'No alignments'),
    ]
    manifest = [(write_bam(tmp_path/(name+'.bam'), HEADER+data, encoding), message)
                for name, data, message in cases]
    probe(tmp_path, '''
import json
for path, message in json.loads(sys.argv[2]):
    for factory in (Parser.BAMParser, Parser.guess_parser):
        p = None
        try:
            p = factory(path)
            p.tsize()
        except ValueError as error:
            assert message in str(error), (path, str(error))
        else:
            raise AssertionError('Unusable BAM accepted: ' + path)
        finally:
            if p is not None:
                p.close()
''', [json.dumps(manifest)])


@pytest.mark.parametrize('fmt', ['BAM', 'AUTO'])
@pytest.mark.parametrize('command', ['filterdup', 'randsample', 'pileup', 'callpeak'])
def test_missing_seq_commands(tmp_path, run_macs3, fmt, command):
    # Original reproduction, including both strands and the original commands.
    reads = b''.join(record(pos=1000+400*i, flag=16*(i % 2))
                     for i in range(2000))
    path = write_bam(tmp_path/'noseq.bam', HEADER+reads, 'bgzf')
    outputs = []
    for explicit in ([False] if command == 'pileup' else [False, True]):
        out = tmp_path/('explicit' if explicit else 'inferred')
        if command == 'callpeak':
            args = [command, '-t', path, '-f', fmt, '-g', '1e6',
                    '--nomodel', '-n', 'test', '--outdir', out]
        else:
            args = [command, '-i', path, '-f', fmt, '-o', out]
            if command == 'filterdup':
                args += ['-g', '1e6']
            elif command == 'randsample':
                args += ['-p', '50', '--seed', '1']
        if explicit:
            args += ['-s', '50']
        result = run_macs3(args, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        assert ('size = 50' in result.stderr
                or 'size is determined as 50' in result.stderr)
        assert 'size = 0' not in result.stderr
        assert 'size is determined as 0' not in result.stderr
        if command == 'callpeak':
            xls = (out/'test_peaks.xls').read_text()
            assert '# tag size is determined as 50 bps' in xls
            outputs.append('\n'.join(line for line in xls.splitlines()
                                     if not line.startswith('#')))
        else:
            outputs.append(out.read_text())
            assert outputs[-1]
    if len(outputs) == 2:
        assert outputs[0] == outputs[1]
    if command in ('filterdup', 'randsample'):
        lines = outputs[0].splitlines()
        assert len(lines) == (2000 if command == 'filterdup' else 1000)
        expected = {f'chr1\t{1000+400*i}\t{1050+400*i}\t.\t.\t'
                    + ('-' if i % 2 else '+') for i in range(2000)}
        assert set(lines) <= expected


@pytest.mark.parametrize('command', ['filterdup', 'randsample', 'pileup', 'callpeak'])
def test_uninferable_cli_diagnostic(tmp_path, run_macs3, command):
    path = write_bam(tmp_path/'unknown.bam', HEADER+record(()), 'bgzf')
    if command == 'callpeak':
        args = [command, '-t', path, '-f', 'BAM', '-g', '1e6', '--nomodel',
                '--outdir', tmp_path, '-n', 'unknown']
    else:
        args = [command, '-i', path, '-f', 'BAM', '-o', tmp_path/'out']
        if command == 'randsample':
            args += ['-p', '50']
    result = run_macs3(args, timeout=30)
    assert result.returncode > 0
    assert 'Cannot infer BAM tag size' in result.stderr
    assert 'for commands supporting it, specify -s/--tsize' in result.stderr
    assert 'AssertionError' not in result.stderr


@pytest.mark.parametrize('command', ['filterdup', 'randsample', 'callpeak'])
def test_explicit_size_bypasses_inference(tmp_path, run_macs3, command):
    path = write_bam(tmp_path/'unknown.bam', HEADER+record(()), 'bgzf')
    if command == 'callpeak':
        args = [command, '-t', path, '-f', 'BAM', '--nomodel', '-g', '1e6',
                '--outdir', tmp_path, '-n', 'explicit']
    else:
        args = [command, '-i', path, '-f', 'BAM', '-o', tmp_path/'out']
        args += (['--keep-dup', 'all'] if command == 'filterdup' else ['-p', '100'])
    result = run_macs3(args+['-s', '50'], timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
