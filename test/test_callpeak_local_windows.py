"""Public API and CLI regressions for long fragments (#807)."""

import math
from types import SimpleNamespace

import pytest

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.Signal.PeakDetect import PeakDetect


def inputs(tmp_path, fmt, length):
    paths = []
    for sample, n, size in [('t', 40, length), ('c', 30, 200)]:
        rows = []
        alignments = []
        for i in range(n):
            start = 30000 + i * 10 if sample == 't' else 10000 + i * 2000
            end = start + size + (i % 2 if sample == 't' and fmt != 'BED' else 0)
            row = f'chr1\t{start}\t{end}'
            if fmt == 'BED':
                row = f'chr1\t{start}\t{start+50}\tr{i}\t0\t+'
            elif fmt == 'FRAG':
                row += '\tbarcode\t1'
            elif fmt == 'BAMPE':
                from test_Parser_bampe_validation import alignment
                alignments.append(alignment(tlen=end-start, pos=start,
                                            mate_pos=end-50))
            rows.append(row)
        path = tmp_path / f'{sample}.{fmt.lower()}'
        if fmt == 'BAMPE':
            from test_Parser_bampe_validation import HEADER
            from test_Parser_bam_truncation import write_bam
            write_bam(path, HEADER + b''.join(alignments), 'bgzf')
        else:
            path.write_text('\n'.join(rows) + '\n')
        paths.append(path)
    return paths


def call(run_macs3, tmp_path, fmt, length, windows, extra=(), name='auto'):
    treatment, control = inputs(tmp_path, fmt, length)
    outdir = tmp_path / name
    args = ['callpeak', '-f', fmt, '-t', treatment, '-c', control,
            '-g', '100000', '--keep-dup', 'all', '--nomodel',
            '--extsize', str(length), '-p', '0.1', '-B',
            '--slocal', str(windows[0]), '--llocal', str(windows[1]),
            '-n', 'result', '--outdir', outdir, *extra]
    return run_macs3(args, timeout=30), outdir


@pytest.mark.parametrize('fmt', ['BEDPE', 'BED'])
@pytest.mark.parametrize('length,windows', [
    (1500, (1000, 10000)), (12000, (1000, 10000)),
    (1500, (0, 1000)), (1500, (1000, 0)),
    (1500, (0, 0)), (200, (1000, 10000)),
    (1500, (1500, 10000)), (1500, (1501, 10000)),
])
@pytest.mark.parametrize('extra', [(), ('--scale-to', 'large')])
def test_cli_matches_explicit_windows(run_macs3, tmp_path, fmt, length,
                                      windows, extra):
    # PE mean is length + 0.5, exercising upward rounding rather than truncation.
    minimum = math.ceil(length + (0.5 if fmt != 'BED' else 0))
    effective = tuple(max(w, minimum) if w else 0 for w in windows)
    result, outdir = call(run_macs3, tmp_path, fmt, length, windows, extra)
    assert result.returncode == 0, result.stderr
    reference, refdir = call(run_macs3, tmp_path, fmt, length, effective,
                            extra, name='explicit')
    assert reference.returncode == 0, reference.stderr
    for suffix in ('_treat_pileup.bdg', '_control_lambda.bdg',
                   '_peaks.narrowPeak', '_summits.bed'):
        actual = (outdir / ('result' + suffix)).read_text()
        assert actual == (refdir / ('result' + suffix)).read_text()
        if suffix.endswith('.bdg'):
            assert actual.strip()
    for option, old, new in zip(('slocal', 'llocal'), windows, effective):
        assert (f'Increasing --{option} from {old} to {new}' in result.stderr) == (old != new)
    metadata = (outdir / 'result_peaks.xls').read_text()
    if effective != windows:
        assert f'# Effective local windows: slocal={effective[0]}, llocal={effective[1]}' in metadata


@pytest.mark.parametrize('fmt', ['BAMPE', 'FRAG'])
@pytest.mark.parametrize('length', [1500, 12000])
def test_other_paired_formats(run_macs3, tmp_path, fmt, length):
    test_cli_matches_explicit_windows(run_macs3, tmp_path, fmt, length,
                                      (1000, 10000), ())


@pytest.mark.parametrize('fmt', ['BEDPE', 'BED'])
@pytest.mark.parametrize('windows', [(1000, 10000), (-1, -2), (2000, 1000)])
def test_nolambda_ignores_local_windows(run_macs3, tmp_path, fmt, windows):
    result, outdir = call(run_macs3, tmp_path, fmt, 1500, windows,
                          extra=('--nolambda',))
    assert result.returncode == 0, result.stderr
    assert 'Increasing --' not in result.stderr
    lines = (outdir / 'result_control_lambda.bdg').read_text().splitlines()
    assert lines
    # Default scaling targets the smaller depth: 40 PE fragments versus
    # 60 control ends, or 30 BED control reads versus 40 treatment reads.
    expected = (40 * 1500.5 if fmt == 'BEDPE' else 30 * 1500) / 100000
    assert all(float(line.split()[3]) == pytest.approx(expected, abs=1e-5)
               for line in lines)


@pytest.mark.parametrize('windows,message', [
    ((-1, 10000), '--slocal and --llocal must be non-negative'),
    ((1000, -1), '--slocal and --llocal must be non-negative'),
    ((2000, 1000), '--llocal must be at least --slocal'),
])
def test_invalid_windows(run_macs3, tmp_path, windows, message):
    result, _ = call(run_macs3, tmp_path, 'BEDPE', 1500, windows)
    assert result.returncode > 0, result.stderr
    assert message in result.stderr


@pytest.mark.parametrize('content', ['', 'chr1\tbad\tend\n'])
def test_invalid_fragments_terminate(run_macs3, tmp_path, content):
    treatment = tmp_path / 'invalid.bedpe'
    treatment.write_text(content)
    control = tmp_path / 'control.bedpe'
    control.write_text('chr1\t100\t300\n')
    result = run_macs3(['callpeak', '-f', 'BEDPE', '-t', treatment,
                        '-c', control, '-g', '100000', '-n', 'invalid'], timeout=30)
    assert result.returncode > 0, result.stderr
    assert 'No valid fragments' in result.stderr or 'Right position' in result.stderr


@pytest.mark.parametrize('paired', [False, True])
def test_public_api_matches_explicit_windows(tmp_path, paired):
    def tracks():
        t, c = (PETrackI(), PETrackI()) if paired else (FWTrack(fw=50), FWTrack(fw=50))
        for track, n, size in [(t, 40, 1500), (c, 20, 200)]:
            for i in range(n):
                start = 30000 + i * 10 if track is t else 10000 + i * 2000
                if paired:
                    track.add_loc(b'chr1', start, start + size + (i % 2 if track is t else 0))
                else:
                    track.add_loc(b'chr1', start, 0)
            track.finalize()
        return t, c

    def detect(small, prefix):
        warnings = []
        opt = SimpleNamespace(
            info=lambda *a: None, debug=lambda *a: None, warn=warnings.append,
            PE_MODE=paired, d=1500.5 if paired else 1500, maxgap=50,
            tsize=50, minlen=50, shift=0, gsize=100000, nolambda=False,
            smalllocal=small, largelocal=10000, ratio=1., tocontrol=False,
            log_pvalue=1., log_qvalue=None, call_summits=False,
            store_bdg=True, name=str(tmp_path / prefix),
            bdg_treat=str(tmp_path / (prefix + '_t.bdg')),
            bdg_control=str(tmp_path / (prefix + '_c.bdg')),
            do_SPMR=False, cutoff_analysis_file='', trackline=False,
            broad=False, cutoff_analysis=False)
        t, c = tracks()
        detector = PeakDetect(opt=opt, treat=t, control=c)
        peaks = detector.call_peaks()
        return detector, peaks, warnings

    effective = 1501 if paired else 1500
    auto, peaks, warnings = detect(1000, 'auto')
    explicit, reference, _ = detect(effective, 'explicit')
    assert auto.sregion == explicit.sregion == effective
    assert warnings and f'Increasing --slocal from 1000 to {effective}' in warnings[0]
    def rows(collection):
        return [tuple(p[key] for key in ('start', 'end', 'summit', 'pileup',
                                         'pscore', 'qscore', 'fc'))
                for p in collection.get_data_from_chrom(b'chr1')]
    assert rows(peaks) == rows(reference)
    assert peaks.get_data_from_chrom(b'chr1')
    assert (tmp_path / 'auto_c.bdg').read_text() == (tmp_path / 'explicit_c.bdg').read_text()
