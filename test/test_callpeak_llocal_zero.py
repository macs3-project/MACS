"""Public API and bounded CLI regressions for disabled llocal (#809)."""

from types import SimpleNamespace

import pytest

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.Signal.PeakDetect import PeakDetect
from test_callpeak_local_windows import inputs


@pytest.mark.parametrize('fmt', ['BED', 'BEDPE', 'BAMPE', 'FRAG'])
@pytest.mark.parametrize('small', [0, 1000])
@pytest.mark.parametrize('broad', [False, True])
def test_cli_llocal_zero_matches_global_background(
        run_macs3, tmp_path, fmt, small, broad):
    treatment, _ = inputs(tmp_path, fmt, 200)
    args = ['callpeak', '-f', fmt, '-t', treatment, '-g', '100000',
            '--keep-dup', 'all', '--nomodel', '--extsize', '200',
            '--slocal', str(small), '-B', '-n', 'result']
    if broad:
        args.append('--broad')
    directories = []
    for name, extra in [('zero', ['--llocal', '0']),
                        ('global', ['--nolambda'])]:
        outdir = tmp_path / name
        result = run_macs3([*args, '--outdir', outdir, *extra], timeout=30)
        assert result.returncode == 0, result.stderr
        directories.append(outdir)
    suffixes = ['_treat_pileup.bdg', '_control_lambda.bdg']
    suffixes += (['_peaks.broadPeak', '_peaks.gappedPeak'] if broad else
                 ['_peaks.narrowPeak', '_summits.bed'])
    for suffix in suffixes:
        actual = (directories[0] / ('result' + suffix)).read_text()
        assert actual.strip(), suffix
        assert actual == (directories[1] / ('result' + suffix)).read_text()
    background = (40 * 200 if fmt == 'BED' else 40 * 200.5) / 100000
    lines = (directories[0] / 'result_control_lambda.bdg').read_text().splitlines()
    assert all(float(line.split()[3]) == pytest.approx(background, abs=1e-5)
               for line in lines)


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('cutoff', ['p', 'q'])
def test_public_api_llocal_zero_matches_global_background(tmp_path, paired,
                                                        cutoff):
    def detect(name, local, nolambda):
        track = PETrackI() if paired else FWTrack(fw=50)
        for i in range(40):
            start = 30000 + i * 10
            if paired:
                track.add_loc(b'chr1', start, start + 200 + i % 2)
            else:
                track.add_loc(b'chr1', start, i % 2)
        track.finalize()
        opt = SimpleNamespace(
            info=lambda *a: None, debug=lambda *a: None, warn=lambda *a: None,
            PE_MODE=paired, d=200, maxgap=50, tsize=50, minlen=50,
            shift=0, gsize=100000, nolambda=nolambda, smalllocal=1000,
            largelocal=local, ratio=1., tocontrol=False,
            log_pvalue=1. if cutoff == 'p' else None,
            log_qvalue=1. if cutoff == 'q' else None,
            call_summits=False, store_bdg=True, name=str(tmp_path / name),
            bdg_treat=str(tmp_path / (name + '_t.bdg')),
            bdg_control=str(tmp_path / (name + '_c.bdg')),
            do_SPMR=False, cutoff_analysis_file='', trackline=False,
            broad=False, cutoff_analysis=False)
        return PeakDetect(opt=opt, treat=track).call_peaks()

    actual = detect('zero', 0, False)
    reference = detect('global', 10000, True)

    def rows(collection):
        return [tuple(peak[key] for key in ('start', 'end', 'summit', 'pileup',
                                            'pscore', 'qscore', 'fc'))
                for peak in collection.get_data_from_chrom(b'chr1')]

    assert rows(actual)
    assert rows(actual) == rows(reference)
    for suffix in ('_t.bdg', '_c.bdg'):
        assert (tmp_path / ('zero' + suffix)).read_text() == (
            tmp_path / ('global' + suffix)).read_text()


@pytest.mark.parametrize('fmt', ['BED', 'BEDPE', 'FRAG'])
@pytest.mark.parametrize('empty', [False, True])
def test_llocal_zero_invalid_input_terminates(run_macs3, tmp_path, fmt, empty):
    if empty:
        content = ''
        diagnostic = 'No usable alignments' if fmt == 'BED' else 'No valid fragments'
    elif fmt == 'BED':
        content = 'chr1\t100\t150\tr1\t0\tinvalid\n'
        diagnostic = 'Strand'
    else:
        content = 'chr1\tbad\tend' + ('\tbarcode\t1' if fmt == 'FRAG' else '') + '\n'
        diagnostic = 'Right position'
    treatment = tmp_path / ('invalid.' + fmt.lower())
    treatment.write_text(content)
    result = run_macs3(['callpeak', '-f', fmt, '-t', treatment, '-g', '100000',
                        '--nomodel', '--extsize', '200', '--llocal', '0',
                        '-n', 'invalid'], timeout=30)
    assert result.returncode > 0, result.stderr
    assert diagnostic in result.stderr
    assert 'ZeroDivisionError' not in result.stderr
