"""Public-API and CLI regressions for empty alignment tracks (#808)."""

import subprocess
import sys

import pytest

from conftest import macs3_env
from test_Parser_bam_truncation import write_bam
from test_Parser_bampe_validation import HEADER, alignment


TRACK_PROBE = '''
import numpy as np
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Signal.Region import Regions

for cls in (PETrackI, PETrackII):
    for operation in ('sample_percent', 'sample_percent_copy',
                      'sample_num', 'sample_num_copy'):
        for populated, fraction in ((False, 0.5), (True, 0.5), (True, 0.0)):
            track = cls(buffer_size=2)
            if populated:
                for chrom in (b'chr1', b'chr2'):
                    if cls is PETrackI:
                        track.add_loc(chrom, 100, 300)
                    else:
                        track.add_loc(chrom, 100, 300, b'BC', 1)
            track.finalize()
            before = (track.total, track.length, track.get_rlengths().copy())
            value = fraction if 'percent' in operation else 0
            result = getattr(track, operation)(value, seed=17)
            sampled = result if operation.endswith('_copy') else track
            assert sampled.total == sampled.length == sampled.average_template_length == 0
            assert sampled.fraglengths().size == 0
            assert sampled.fraglengths().dtype == np.int32
            if operation.endswith('_copy'):
                assert (track.total, track.length, track.get_rlengths()) == before
            sampled.finalize()
            if cls is PETrackI:
                sampled.filter_dup(1)
            assert sampled.total == sampled.length == sampled.average_template_length == 0

    track = cls(buffer_size=2)
    if cls is PETrackI:
        track.add_loc(b'chr1', 100, 300)
    else:
        track.add_loc(b'chr1', 100, 300, b'BC', 2)
    track.finalize()
    regions = Regions()
    regions.add_loc(b'chr1', 0, 1000)
    track.exclude(regions)
    assert track.total == track.length == track.average_template_length == 0
    assert track.fraglengths().size == 0
'''


def test_empty_track_operations(tmp_path):
    script = tmp_path / "track_probe.py"
    script.write_text(TRACK_PROBE)
    result = subprocess.run(
        [sys.executable, str(script)], cwd=tmp_path,
        env=macs3_env(tmp_path), capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def assert_empty_error(result, message):
    assert result.returncode == 1, result.stderr
    assert message in result.stderr
    assert "Traceback" not in result.stderr
    assert "ZeroDivisionError" not in result.stderr


@pytest.mark.parametrize("keep", ["1", "auto", "all"])
@pytest.mark.parametrize("role", ["treatment", "control"])
def test_callpeak_empty_bed(tmp_path, run_macs3, keep, role):
    empty = tmp_path / "empty.bed"
    empty.touch()
    valid = tmp_path / "valid.bed"
    valid.write_text("chr1\t100\t150\tr\t0\t+\n")
    args = ["callpeak", "-f", "BED", "-g", "1e5", "--nomodel",
            "--keep-dup", keep, "-n", "empty", "-t",
            empty if role == "treatment" else valid]
    if role == "control":
        args += ["-c", empty]
    result = run_macs3(args, timeout=30)
    assert_empty_error(result, f"No usable alignments remain in {role}")
    assert not (tmp_path / "empty_peaks.narrowPeak").exists()


@pytest.mark.parametrize("keep", ["1", "auto", "all"])
def test_filterdup_empty_bed(tmp_path, run_macs3, keep):
    empty = tmp_path / "empty.bed"
    empty.touch()
    result = run_macs3(["filterdup", "-f", "BED", "-g", "1e5",
                        "--keep-dup", keep, "-i", empty, "-o", "out.bed"],
                       timeout=30)
    assert_empty_error(result, "No usable alignments remain in input")


@pytest.mark.parametrize("flag", [97, 147, 73])
@pytest.mark.parametrize("command", ["callpeak", "filterdup", "randsample"])
@pytest.mark.parametrize("fmt", ["BAM", "BAMPE"])
def test_bam_with_no_accepted_reads(tmp_path, run_macs3, flag, command, fmt):
    path = write_bam(tmp_path / "filtered.bam", HEADER + alignment(flag=flag),
                     "bgzf")
    if command == "callpeak":
        args = ["-t", path, "-g", "1e5", "--nomodel", "-n", "empty"]
    else:
        args = ["-i", path, "-o", "out.bed"]
        args += ["-g", "1e5"] if command == "filterdup" else ["-p", "50"]
    result = run_macs3([command, "-f", fmt, *args], timeout=30)
    if fmt == "BAMPE":
        assert_empty_error(result, "No valid fragments found in BAMPE input")
    else:
        assert_empty_error(result, "No usable alignments remain")
        assert "proper pairs (flag 0x2)" in result.stderr


@pytest.mark.parametrize("command, fmt", [
    ("callpeak", "BEDPE"), ("callpeak", "FRAG"),
    ("hmmratac", "BEDPE"), ("hmmratac", "FRAG"),
    ("filterdup", "BEDPE"), ("randsample", "BEDPE"),
])
def test_empty_paired_input(tmp_path, run_macs3, command, fmt):
    path = tmp_path / "empty.tsv"
    path.touch()
    if command == "callpeak":
        args = ["-t", path, "-g", "1e5", "-n", "empty"]
    elif command == "hmmratac":
        args = ["-i", path, "-n", "empty"]
    else:
        args = ["-i", path, "-o", "out.bed"]
        args += ["-g", "1e5"] if command == "filterdup" else ["-p", "50"]
    result = run_macs3([command, "-f", fmt, *args], timeout=30)
    assert_empty_error(result, f"No valid fragments found in {fmt} input")


@pytest.mark.parametrize("percent", [0, 50, 100])
def test_randsample_sparse_succeeds(tmp_path, run_macs3, percent):
    data = alignment() + alignment(ref=1, mate_ref=1)
    path = write_bam(tmp_path / "sparse.bam", HEADER + data, "bgzf")
    result = run_macs3(["randsample", "-f", "BAMPE", "-i", path,
                        "-p", str(percent), "-o", "sample.bed"], timeout=30)
    assert result.returncode == 0, result.stderr
    expected = "chr1\t100\t300\nchr2\t100\t300\n" if percent == 100 else ""
    assert (tmp_path / "sample.bed").read_text() == expected


@pytest.mark.parametrize("kind", ["sparse", "out_of_range"])
def test_hmmratac_empty_em_sample(tmp_path, run_macs3, kind):
    if kind == "sparse":
        data = alignment() + alignment(ref=1, mate_ref=1)
        message = "No fragments remain after sampling"
    else:
        data = b"".join(alignment(tlen=2000, pos=100 + i * 3000,
                                  mate_pos=2050 + i * 3000)
                        for i in range(200))
        message = "No sampled fragments within the EM fragment-length range"
    path = write_bam(tmp_path / "input.bam", HEADER + data, "bgzf")
    result = run_macs3(["hmmratac", "-f", "BAMPE", "-i", path,
                        "-n", "empty"], timeout=30)
    assert_empty_error(result, message)
    assert "--no-fragem" in result.stderr


@pytest.mark.parametrize("command", ["callpeak", "hmmratac"])
def test_no_matching_barcodes(tmp_path, run_macs3, command):
    fragments = tmp_path / "input.frag"
    fragments.write_text("chr1\t100\t300\tBC\t2\n")
    barcodes = tmp_path / "barcodes.txt"
    barcodes.write_text("NO_MATCH\n")
    args = ["-t" if command == "callpeak" else "-i", fragments]
    if command == "callpeak":
        args += ["-g", "1e5", "--nomodel"]
    result = run_macs3([command, "-f", "FRAG", *args, "--barcodes", barcodes,
                        "-n", "empty"], timeout=30)
    assert_empty_error(result, "No usable alignments remain")


@pytest.mark.parametrize("fmt", ["BEDPE", "FRAG"])
def test_hmmratac_blacklist_removes_all(tmp_path, run_macs3, fmt):
    path = tmp_path / "input.tsv"
    path.write_text("chr1\t100\t300" + ("\tBC\t2" if fmt == "FRAG" else "")
                    + "\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t1000\n")
    result = run_macs3(["hmmratac", "-f", fmt, "-i", path, "--blacklist",
                        blacklist, "-n", "empty"], timeout=30)
    assert_empty_error(result, "No usable alignments remain in input after blacklist")
