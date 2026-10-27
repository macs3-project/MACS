"""Regressions for #799 / #757: pileups after sampling empties a chromosome."""

from pathlib import Path
import random
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


def test_legacy_pileups_empty_and_one_sided(tmp_path):
    run_probe(tmp_path, """
        from pathlib import Path
        import sys
        import numpy as np
        import MACS3.Signal.PileupV2 as pileup

        assert Path(pileup.__file__).resolve().is_relative_to(Path(sys.argv[1]))
        assert Path(pileup.__file__).suffix != '.py'
        for scale, baseline in ((1.0, 0.0), (0.5, 0.25), (0.5, 2.0)):
            for plus, minus, positions in (([], [], []),
                                           ([100], [], [100, 150]),
                                           ([], [200], [150, 200])):
                p, n = np.array(plus, dtype='i4'), np.array(minus, dtype='i4')
                before = (p.copy(), n.copy())
                for repeat in range(2):
                    pos, val = pileup.pileup_from_PN_shifted(
                        p, n, 0, 50, 1000, scale, baseline)
                    assert pos.dtype == np.int32 and val.dtype == np.float32
                    np.testing.assert_array_equal(pos, positions)
                    expected = [baseline, max(scale, baseline)] if positions else []
                    np.testing.assert_array_equal(val, expected)
                    np.testing.assert_array_equal(p, before[0])
                    np.testing.assert_array_equal(n, before[1])
            for rows in ([], [(100, 150)]):
                data = np.array(rows, dtype=pileup.LR_DTYPE)
                data.flags.writeable = False
                for repeat in range(2):
                    pos, val = pileup.pileup_from_LR_as_list(data, scale, baseline)
                    assert pos.dtype == np.int32 and val.dtype == np.float32
                    np.testing.assert_array_equal(pos, [100, 150] if rows else [])
                    expected = [baseline, max(scale, baseline)] if rows else []
                    np.testing.assert_array_equal(val, expected)
                    np.testing.assert_array_equal(data, np.array(rows, dtype=pileup.LR_DTYPE))
    """)


def test_sampled_empty_track_pileups(tmp_path):
    run_probe(tmp_path, """
        from pathlib import Path
        import sys
        import numpy as np
        import MACS3.Signal.FixWidthTrack as fw
        import MACS3.Signal.PairedEndTrack as pe

        for module in (fw, pe):
            assert Path(module.__file__).resolve().is_relative_to(Path(sys.argv[1]))
            assert Path(module.__file__).suffix != '.py'
        for paired in (False, True):
            track = pe.PETrackI() if paired else fw.FWTrack(fw=50)
            for i in range(10):
                if paired:
                    track.add_loc(b'chr1', 100+100*i, 150+100*i)
                else:
                    track.add_loc(b'chr1', 100+100*i, i % 2)
            if paired:
                track.add_loc(b'chr2', 100, 150)
            else:
                track.add_loc(b'chr2', 100, 0)
                track.add_loc(b'chr2', 200, 1)
            track.finalize()
            track.sample_percent(0.2, seed=7)
            assert track.total == 2
            assert b'chr2' in track.get_chr_names()
            locs = track.get_locations_by_chr(b'chr2')
            assert len(locs) == 0 if paired else all(len(x) == 0 for x in locs)
            for repeat in range(2):
                for baseline in (0.0, 0.25):
                    if paired:
                        result = track.pileup_a_chromosome(b'chr2', 0.5, baseline)
                    else:
                        result = track.pileup_a_chromosome(b'chr2', 200, 0.5, baseline)
                    results = [result]
                    for ds in ([200], [200, 1000], [200, 1000, 10000]):
                        results.append(track.pileup_a_chromosome_c(
                            b'chr2', ds, [0.5]*len(ds), baseline))
                    for pos, val in results:
                        assert pos.dtype == np.int32 and val.dtype == np.float32
                        assert pos.size == val.size == 0
            # An empty chromosome must not discard retained reads elsewhere.
            retained = (track.pileup_a_chromosome(b'chr1') if paired else
                        track.pileup_a_chromosome(b'chr1', 50))
            assert retained[0].size and retained[1].max() > 0
            try:
                track.pileup_a_chromosome_c(b'chr2', [200], [])
            except AssertionError as exc:
                assert 'same length' in str(exc)
            else:
                raise AssertionError('mismatched control windows accepted')
    """)


def test_legacy_pileup_missing_fields_rejected(tmp_path):
    run_probe(tmp_path, """
        import numpy as np
        from MACS3.Signal.PileupV2 import pileup_from_LR_as_list

        for n in (0, 1):
            for dtype in ([('l', 'i4')], [('r', 'i4')]):
                try:
                    pileup_from_LR_as_list(np.zeros(n, dtype=dtype))
                except ValueError as exc:
                    assert 'no field' in str(exc)
                else:
                    raise AssertionError('missing interval coordinate accepted')
    """)


@pytest.mark.parametrize("fmt", ["BED", "BEDPE", "BAMPE"])
@pytest.mark.parametrize("swapped", [False, True])
@pytest.mark.parametrize("downsample", [False, True])
def test_callpeak_downsample_small_chromosome(
        run_macs3, tmp_path, request, make_pe_pair,
        fmt, swapped, downsample):
    if fmt == "BAMPE":
        make_alignments = request.getfixturevalue("make_alignments")
    rng = random.Random(1)
    paths = []
    for name, count in (("big", 1000), ("small", 200)):
        rows, reads = [], []
        for chrom, n, span in (("chr1", count, 200000), ("chr2", 1, 5000)):
            if fmt == "BED":
                for strand in "+-":
                    for _ in range(n):
                        start = rng.randrange(1000, span)
                        rows.append(f"{chrom}\t{start}\t{start+50}\tr\t0\t{strand}\n")
            else:
                for i in range(n):
                    start = rng.randrange(1000, span)
                    rows.append(f"{chrom}\t{start}\t{start+200}\n")
                    if fmt == "BAMPE":
                        reads.extend(make_pe_pair(f"{chrom}_{i}", chrom, start,
                                                  start+200, readlen=50))
        if fmt == "BAMPE":
            path = make_alignments(reads, refs=(("chr1", 1000000), ("chr2", 1000000)),
                                   name=f"{name}.bam")
        else:
            path = tmp_path / f"{name}.bed"
            path.write_text("".join(rows))
        paths.append(path)
    if swapped:
        paths.reverse()
    args = ["callpeak", "-f", fmt, "-t", paths[0], "-c", paths[1],
            "-g", "1e6", "--nomodel", "--extsize", "200", "--seed", "7",
            "-n", "sampled"]
    if downsample:
        args.append("--down-sample")
    result = run_macs3(args, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Done!" in result.stderr
    assert (tmp_path / "sampled_peaks.xls").is_file()
    assert (tmp_path / "sampled_peaks.narrowPeak").is_file()
