"""Public compiled-track and CLI regressions for #802 (part of #757)."""

from pathlib import Path
import random

import pytest

from MACS3.Commands import callpeak_cmd, callvar_cmd
from MACS3.Signal import FixWidthTrack, PairedEndTrack


@pytest.mark.parametrize("command", [callpeak_cmd, callvar_cmd])
@pytest.mark.parametrize("track_module", [FixWidthTrack, PairedEndTrack])
@pytest.mark.parametrize("treatment,control", [
    ([b"chr2", b"chr1"], [b"1", b"2"]),
    ([], [b"chr1"]), ([b"chr1"], []), ([], []),
    ([b"chr\xff"], [b"chr1"]),
    (["chré".encode()], [b"chr1"]),
    ([b"chr1", b"chr2"], [b"chr2", b"chr3"]),
])
def test_compiled_track_names(command, track_module, treatment, control):
    assert Path(track_module.__file__).suffix != ".py"
    tracks = []
    for names in (treatment, control):
        if track_module is FixWidthTrack:
            track = track_module.FWTrack(buffer_size=10)
            for name in names:
                track.add_loc(name, 100, 0)
        else:
            track = track_module.PETrackI(buffer_size=10)
            for name in names:
                track.add_loc(name, 100, 150)
        track.finalize()
        tracks.append(track)
    messages = []
    if set(treatment).intersection(control):
        assert command.check_names(*tracks, messages.append) is None
        assert messages == []
        return
    with pytest.raises(SystemExit) as exc:
        command.check_names(*tracks, messages.append)
    assert exc.value.code == 1
    assert "No common chromosome names" in messages[0]
    assert messages[-2:] == [
        "Chromosome names in treatment: " + ",".join(sorted(
            n.decode("utf-8", errors="backslashreplace") for n in treatment)),
        "Chromosome names in control: " + ",".join(sorted(
            n.decode("utf-8", errors="backslashreplace") for n in control)),
    ]


@pytest.mark.parametrize("fmt,treatment,control", [
    (fmt, treatment, control)
    for fmt in ("BED", "BEDPE", "FRAG")
    for treatment, control in (
        ([b"chr2", b"chr1"], [b"1"]),
        ([], [b"chr1"]), ([b"chr1"], []), ([], []),
        ([b"chr\xff"], [b"chr1"]),
    )
    # Empty paired-end files fail during loading, before the name check (#808).
    if fmt == "BED" or (treatment and control)
])
def test_callpeak_disjoint_names(run_macs3, tmp_path, fmt, treatment, control):
    for filename, names in (("treatment", treatment), ("control", control)):
        with (tmp_path / filename).open("wb") as stream:
            for name in names:
                for i in range(20):
                    row = name + f"\t{100 + i*100}\t{150 + i*100}".encode()
                    if fmt == "BED":
                        row += b"\tr\t0\t+"
                    elif fmt == "FRAG":
                        row += b"\tBC\t1"
                    stream.write(row + b"\n")
    result = run_macs3([
        "callpeak", "-t", tmp_path / "treatment", "-c", tmp_path / "control",
        "-f", fmt, "-g", "1e6", "--nomodel", "--extsize", "100", "-n", "x",
    ], timeout=30)
    assert result.returncode == 1, result.stderr
    if not treatment or not control:
        role = "treatment" if not treatment else "control"
        assert f"No usable alignments remain in {role}" in result.stderr
        assert "No common chromosome names" not in result.stderr
    else:
        assert "No common chromosome names" in result.stderr
        for label, names in (("treatment", treatment), ("control", control)):
            expected = ",".join(sorted(n.decode("utf-8", errors="backslashreplace")
                                       for n in names))
            assert f"Chromosome names in {label}: {expected}" in result.stderr
    assert "Traceback" not in result.stderr
    assert not list(tmp_path.glob("x_peaks.*"))


def test_original_reproduction_and_overlap(run_macs3, tmp_path):
    random.seed(0)
    for filename, chrom in (("t.bed", "chr1"), ("c.bed", "chr2")):
        with (tmp_path / filename).open("w") as stream:
            for i in range(2000):
                p = random.randint(1000, 900000)
                stream.write(f'{chrom}\t{p}\t{p+50}\tr{i}\t0\t{"+-"[i%2]}\n')
    args = ["callpeak", "-t", tmp_path / "t.bed", "-f", "BED", "-g", "1e6",
            "-n", "x", "--nomodel", "--extsize", "100", "-c"]
    result = run_macs3([*args, tmp_path / "c.bed"], timeout=30)
    assert result.returncode == 1
    assert "Chromosome names in treatment: chr1" in result.stderr
    assert "Chromosome names in control: chr2" in result.stderr
    assert "Traceback" not in result.stderr
    result = run_macs3([*args, tmp_path / "t.bed"], timeout=30)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "x_peaks.narrowPeak").exists()


@pytest.mark.parametrize("treatment,control,overlap", [
    (("chr2", "chr1"), ("1", "2"), False),
    (("chr1", "chr2"), ("chr3", "chr2"), True),
    (("chr1",), ("chr1",), True),
    ((), ("chr1",), False), (("chr1",), (), False), ((), (), False),
])
@pytest.mark.parametrize("optimized", [False, True])
def test_callvar_header_names(make_alignments, run_macs3, tmp_path,
                             treatment, control, overlap, optimized):
    # Header-only BAMs isolate the name check from variant/assembly behavior.
    treat = make_alignments([], refs=tuple((n, 10000) for n in treatment),
                            name="t.bam", index=True)
    ctrl = make_alignments([], refs=tuple((n, 10000) for n in control),
                           name="c.bam", index=True)
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("")
    out = tmp_path / "out.vcf"
    result = run_macs3(["callvar", "-t", treat, "-c", ctrl, "-b", peaks,
                        "-o", out], timeout=30,
                       env={"PYTHONOPTIMIZE": "1" if optimized else "0"})
    assert result.returncode == (0 if overlap else 1), result.stderr
    assert "Traceback" not in result.stderr
    if overlap:
        assert out.read_text().startswith("##fileformat=VCF")
        assert "No common chromosome names" not in result.stderr
    else:
        assert "No common chromosome names" in result.stderr
        assert ("Chromosome names in treatment: " + ",".join(sorted(treatment))) in result.stderr
        assert ("Chromosome names in control: " + ",".join(sorted(control))) in result.stderr
        assert not out.exists()
