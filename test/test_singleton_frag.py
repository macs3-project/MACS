"""Regressions for #798 (part of #757): singleton FRAG chromosomes.

Keep native pileup probes in bounded subprocesses using the parent's
installation. The tiny inputs exercise both callpeak's chromosome pileup
and hmmratac's bedGraph pileup without needing HMM training regions.
"""

from pathlib import Path
import subprocess
import sys
import textwrap

import MACS3
import pytest

from conftest import macs3_env, no_core_dumps


def run_probe(tmp_path, source):
    script = tmp_path / "probe.py"
    script.write_text(textwrap.dedent(source))
    package_root = Path(MACS3.__file__).resolve().parent.parent
    result = subprocess.run(
        [sys.executable, str(script), str(package_root)],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True,
        text=True, timeout=30, preexec_fn=no_core_dumps)
    assert result.returncode == 0, result.stdout + result.stderr


def test_lrc_small_pileups_preserve_input(tmp_path):
    run_probe(tmp_path, """
        from pathlib import Path
        import sys
        import numpy as np
        import MACS3.Signal.PileupV2 as pileup

        assert Path(pileup.__file__).resolve().is_relative_to(Path(sys.argv[1]))
        assert Path(pileup.__file__).suffix != '.py'
        for n in (0, 1, 2):
            for count in (1, 7, 65535):
                rows = [(10, 30, count), (20, 40, count)][:n]
                for left_sorted in (False, True):
                    for readonly in (False, True):
                        data = np.array(rows, dtype=pileup.LRC_DTYPE)
                        if not left_sorted:
                            data = data[::-1]
                        if readonly:
                            data.flags.writeable = False
                        before = data.copy()
                        for scale, baseline in ((1.0, 0.0), (0.5, 0.25)):
                            p, v = pileup.pileup_from_LRC_as_list(
                                data, scale, baseline, left_sorted)
                            positions = ([], [10, 30], [10, 20, 30, 40])[n]
                            values = ([], [baseline, count * scale],
                                      [baseline, count * scale,
                                       2 * count * scale, count * scale])[n]
                            np.testing.assert_array_equal(p, positions)
                            np.testing.assert_array_equal(v, values)
                            np.testing.assert_array_equal(data, before)
    """)


def test_frag_track_small_chromosome_pileups(tmp_path):
    run_probe(tmp_path, """
        from pathlib import Path
        import sys
        import numpy as np
        import MACS3.IO.Parser as parser

        assert Path(parser.__file__).resolve().is_relative_to(Path(sys.argv[1]))
        assert Path(parser.__file__).suffix != '.py'
        for n in (0, 1, 2):
            for count in (1, 7, 65535):
                path = Path('fragments.tsv')
                rows = ['chr1\\t100\\t150\\tBC0\\t1\\n']
                rows += [f'chr2\\t{10+10*i}\\t{30+10*i}\\tBC0\\t{count}\\n'
                         for i in range(n)]
                path.write_text(''.join(rows))
                track = parser.FragParser(str(path)).build_petrack()
                track.finalize()
                before = track.get_locations_by_chr(b'chr1').copy()
                if n:
                    before2 = track.get_locations_by_chr(b'chr2').copy()
                for repeat in range(2):
                    bdg = track.pileup_bdg()
                    assert bdg.get_chr_names() == ({b'chr1', b'chr2'} if n else {b'chr1'})
                    for chrom in track.get_chr_names():
                        if chrom == b'chr1':
                            expected = ([100, 150], [0, 1])
                        elif n == 1:
                            expected = ([10, 30], [0, count])
                        else:
                            expected = ([10, 20, 30, 40], [0, count, 2*count, count])
                        for p, v in (track.pileup_a_chromosome(chrom),
                                     bdg.get_data_by_chr(chrom)):
                            np.testing.assert_array_equal(p, expected[0])
                            np.testing.assert_array_equal(v, expected[1])
                np.testing.assert_array_equal(track.get_locations_by_chr(b'chr1'), before)
                if n:
                    np.testing.assert_array_equal(track.get_locations_by_chr(b'chr2'), before2)
    """)


def test_frag_empty_and_malformed_inputs(tmp_path):
    run_probe(tmp_path, """
        from pathlib import Path
        from MACS3.IO.Parser import FragParser

        cases = [('', ValueError, 'No valid fragments'),
                 ('chr1\\t10\\t30\\tBC0\\n', Exception, 'Less than 5 columns'),
                 ('chr1\\t30\\t10\\tBC0\\t1\\n', AssertionError,
                  'Right position must be larger'),
                 ('chr1\\t10\\t10\\tBC0\\t1\\n', AssertionError,
                  'Right position must be larger')]
        for text, error, message in cases:
            path = Path('invalid.tsv')
            path.write_text(text)
            try:
                track = FragParser(str(path)).build_petrack()
                track.finalize()
            except error as exc:
                assert message in str(exc), str(exc)
            else:
                raise AssertionError('invalid FRAG input was accepted')
    """)


@pytest.mark.parametrize("n", [0, 1, 2], ids=["absent", "singleton", "two"])
def test_callpeak_frag_small_chromosome(run_macs3, write_frag, tmp_path, n):
    rows = [("chr1", 1000 + 37*i, 1150 + 37*i + (i % 7)*10,
             f"BC{i % 5}", 1) for i in range(200)]
    rows += [("chr2", 5000 + 100*i, 5200 + 100*i, "BC0", 1)
             for i in range(n)]
    path = write_frag(rows)
    result = run_macs3(["callpeak", "-f", "FRAG", "-g", "1e6",
                        "-t", path, "-n", "singleton", "-B"], timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "singleton_peaks.xls").is_file()
    rows = [line.split() for line in
            (tmp_path / "singleton_treat_pileup.bdg").read_text().splitlines()]
    observed = [(int(start), int(end), float(value))
                for chrom, start, end, value in rows if chrom == "chr2"]
    expected = ([], [(0, 5000, 0.0), (5000, 5200, 1.0)],
                [(0, 5000, 0.0), (5000, 5100, 1.0),
                 (5100, 5200, 2.0), (5200, 5300, 1.0)])[n]
    assert observed == expected
