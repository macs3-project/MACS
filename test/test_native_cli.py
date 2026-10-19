"""End-to-end native parser/caller combinations beyond individual fix tests."""

import struct
import sys

import pytest

from MACS3.Utilities.Performance import DEFAULTS
from MACS3.Utilities.Execution import resolve_execution_policy
from test_native_streams import encoded


def bam_record(chrom, left):
    # One mapped first mate: the same interval is used for BAMPE and FRAG.
    body = (struct.pack('<iiIIiiii', chrom, left, (30 << 8) | 2,
                        (99 << 16) | 1, 36, chrom, left + 114, 150)
            + b'r\0' + struct.pack('<I', 36 << 4) + b'\x11' * 18 + b'\x1e' * 36)
    return struct.pack('<i', len(body)) + body


def input_data(kind, rows):
    if kind == 'frag':
        return b''.join(b'chr%d\t%d\t%d\tBC%d\t1\n' %
                        (chrom, left, left + 150, left % 3) for chrom, left in rows)
    header = b'BAM\1' + struct.pack('<ii', 0, 3)
    for chrom in range(3):
        name = b'chr%d\0' % chrom
        header += struct.pack('<i', len(name)) + name + struct.pack('<i', 100000)
    return header + b''.join(bam_record(chrom, left) for chrom, left in rows)


def artifacts(outdir, mode, extra):
    suffix = '.broadPeak' if mode == 'broad' else '.narrowPeak'
    assert (outdir / ('peaks_peaks' + suffix)).read_bytes().strip()
    if extra:
        assert (outdir / 'peaks_treat_pileup.bdg').read_bytes().strip()
        assert (outdir / 'peaks_control_lambda.bdg').read_bytes().strip()
        assert any('cutoff' in path.name and path.stat().st_size for path in outdir.iterdir())
    stable = {p.name: p.read_bytes() for p in outdir.iterdir() if p.suffix != '.xls'}
    stable['xls-data'] = b'\n'.join(
        line for line in (outdir / 'peaks_peaks.xls').read_bytes().splitlines()
        if line and not line.startswith(b'#'))
    return stable


@pytest.mark.parametrize('kind', ['bam', 'bampe', 'frag'])
@pytest.mark.parametrize('mode', ['narrow', 'summits', 'broad'])
def test_multi_input_native_cli_artifacts(tmp_path, run_macs3, kind, mode):
    # Three invocations take ~1–2 s total; tiny generated inputs only.
    rows = [(chrom, left) for chrom in range(3)
            for left in [1000 + chrom * 100] * 40 + [9000]]
    outputs = []
    for index, (budget, backend, encoding, split) in enumerate([
            (1, 'serial', 'plain' if kind == 'frag' else 'raw', False),
            (4, 'serial', 'gzip', True), (4, 'auto', 'bgzf', True)]):
        parts = (rows[::2], rows[1::2]) if split else (rows,)
        paths = [encoded(tmp_path, kind, input_data(kind, part), encoding,
                          f'input-{index}-{part_index}')
                 for part_index, part in enumerate(parts)]
        outdir = tmp_path / str(index)
        args = ['callpeak', '-t', *paths, '-f', kind.upper(), '-g', '100000',
                '-n', 'peaks', '--outdir', outdir, '--threads', budget,
                '--process-backend', backend, '--keep-dup', 'all', '-B', '--cutoff-analysis']
        if kind == 'bam':
            args += ['--nomodel', '--extsize', 150]
        if mode != 'narrow':
            args += ['--broad' if mode == 'broad' else '--call-summits']
        result = run_macs3(args, timeout=45)
        assert result.returncode == 0, result.stderr
        outputs.append(artifacts(outdir, mode, True))
    assert outputs[0] == outputs[1] == outputs[2]


@pytest.mark.slow
@pytest.mark.parametrize('kind', ['bam', 'bampe', 'frag'])
def test_native_cli_reaches_production_backend(tmp_path, run_macs3, kind):
    # ~3–5 s per case; 262k stored records exceed the actual process gate.
    # Serial/native-output equality on macOS; real fork workers on Linux.
    rows = [(chrom, 1000 + chrom * 100) for chrom in range(3)
            for _ in range(DEFAULTS.min_parallel_peak_locations // 3 + 1)]
    path = encoded(tmp_path, kind, input_data(kind, rows), 'bgzf', 'large')
    outputs = []
    for budget, backend in [(1, 'serial'), (4, 'auto')]:
        outdir = tmp_path / backend
        args = ['callpeak', '-t', path, '-f', kind.upper(), '-g', '100000',
                '-n', 'peaks', '--outdir', outdir, '--threads', budget,
                '--process-backend', backend, '--keep-dup', 'all']
        if kind == 'bam':
            args += ['--nomodel', '--extsize', 150]
        result = run_macs3(args, timeout=45)
        assert result.returncode == 0, result.stderr
        if (backend == 'auto' and sys.platform == 'linux' and
                resolve_execution_policy(4).threads > 1):
            assert 'Peak process backend: fork, workers=' in result.stderr
        else:
            assert 'Peak process backend: fork' not in result.stderr
        outputs.append(artifacts(outdir, 'narrow', False))
    assert outputs[0] == outputs[1]
