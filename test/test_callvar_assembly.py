"""Public API and CLI regressions for #813 / #757.

Native failures run in bounded subprocesses using the parent's installation.
"""

from pathlib import Path
import subprocess
import sys

import pytest

from conftest import macs3_env, no_core_dumps
from test_callvar_cmd import (REF1, GENOME, REFS, DEL,
                              het_reads, hom_reads, records, write_peaks)


@pytest.mark.parametrize("case", [
    "short", "equal", "filtered", "ambiguous", "mixed-treatment",
    "mixed-control", "minimum", "maximum", "invalid-0", "invalid--1",
    "invalid-64", "invalid-2147483648", "bad-span", "splice",
])
def test_assembly_api(tmp_path, case):
    from MACS3.Signal import RACollection
    probe = Path(__file__).with_name("callvar_assembly_probe.py")
    result = subprocess.run(
        [sys.executable, str(probe), RACollection.__file__, case],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True, text=True,
        timeout=30, preexec_fn=no_core_dumps)
    assert result.returncode == 0, result.stdout + result.stderr


def run(run_macs3, tmp_path, bam, peaks, mode, *extra, name="out"):
    output = tmp_path / (name + ".vcf")
    result = run_macs3(["callvar", "-t", bam, "-b", peaks, "-o", output,
                       "-F", mode, "-m", "1", *extra], timeout=30)
    return result, output


@pytest.mark.parametrize("mode", ["auto", "on", "off"])
@pytest.mark.parametrize("overlap", [-1, 0, 64, 2**31])
def test_invalid_overlap_rejected_before_inputs(tmp_path, run_macs3, mode,
                                               overlap):
    res, out = run(run_macs3, tmp_path, "missing.bam", "missing.bed", mode,
                   "--fermi-overlap", overlap)
    assert res.returncode > 0
    assert "--fermi-overlap must be between 1 and 63" in res.stderr
    assert "Traceback" not in res.stderr
    assert not out.exists()


def short_reads(indel=False):
    reads = []
    for i in range(20):
        start = DEL - 18 + i // 2
        k = DEL - start
        seq = (REF1[start:DEL] + REF1[DEL + 1:start + 26]
               if indel else REF1[start:start + 25])
        cigar = f"{k}M1D{25-k}M" if indel else "25M"
        md = f"{k}^{REF1[DEL]}{25-k}" if indel else "25"
        reads.append(dict(name=f"short{i}", ref="chr1", pos=start,
                          flag=16 * (i % 2), cigar=cigar, seq=seq,
                          tags=[("MD", md)]))
    return reads


@pytest.mark.parametrize("mode", ["auto", "on", "off"])
@pytest.mark.parametrize("indel", [False, True])
def test_short_reads_and_empty_unitigs(tmp_path, make_alignments, run_macs3,
                                      mode, indel):
    # Keep an ordinary SNV as a control alongside the short-read peak.
    bam = make_alignments(short_reads(indel) + hom_reads(), refs=REFS, index=True)
    peaks = write_peaks(tmp_path, [("chr1", 1000, 1200), ("chr1", 600, 800)])
    res, out = run(run_macs3, tmp_path, bam, peaks, mode)
    assert res.returncode == 0, res.stderr
    expected, expected_out = run(run_macs3, tmp_path, bam, peaks, "off",
                                 name="expected")
    assert expected.returncode == 0, expected.stderr
    if indel:
        assert [r for r in records(expected_out) if 1000 <= int(r[1])-1 < 1200]
    if mode != "on":
        assert records(out) == records(expected_out)
    else:
        assert not [r for r in records(out) if 1000 <= int(r[1])-1 < 1200]
    assert "Peak: chr1 1000 1200" in res.stderr


def spliced_read():
    ref = GENOME["chr1"]
    return dict(name="splice", ref="chr1", pos=200, cigar="25M2000N25M",
                seq=ref[200:225] + ref[2225:2250], tags=[("MD", "50")])


@pytest.mark.parametrize("mode", ["auto", "on", "off"])
@pytest.mark.parametrize("sample", ["treatment", "control", "only"])
def test_spliced_reads_skipped(tmp_path, make_alignments, run_macs3, mode,
                               sample):
    ordinary = het_reads()
    treatment = ([] if sample == "only" else ordinary)
    if sample != "control":
        treatment = treatment + [spliced_read()]
    bam = make_alignments(treatment, refs=REFS, index=True)
    extra = []
    if sample == "control":
        control = make_alignments([spliced_read()], refs=REFS,
                                  name="control.bam", index=True)
        extra = ["-c", control]
    peaks = write_peaks(tmp_path, [("chr1", 190, 2300)])
    res, out = run(run_macs3, tmp_path, bam, peaks, mode, *extra)
    assert res.returncode == 0, res.stderr
    assert "Skipped 1 spliced" in res.stderr
    if sample == "only":
        assert records(out) == []
        assert "No reads found in this peak. Skipped" in res.stderr
    else:
        clean = make_alignments(ordinary, refs=REFS, name="clean.bam", index=True)
        expected, expected_out = run(run_macs3, tmp_path, clean, peaks, mode,
                                     name="expected")
        assert expected.returncode == 0, expected.stderr
        assert records(out) == records(expected_out)


def test_empty_assembly_continues_to_later_peak(tmp_path, make_alignments,
                                               run_macs3):
    # Four 30 bp reads have no unitigs. The next peak tiles 300 bp, carries
    # a homozygous deletion, and has a known nonempty successful assembly.
    first = [dict(name=f"empty{i}", ref="chr1", pos=100+i*10, cigar="30M",
                  seq=REF1[100+i*10:130+i*10], tags=[("MD", "30")])
             for i in range(4)]
    tiled = []
    for i, start in enumerate(range(DEL-140, DEL+110, 4)):
        if start <= DEL < start+50:
            seq = REF1[start:DEL] + REF1[DEL+1:start+51]
            k = DEL-start
            cigar = (f"{k}M" if k else "") + f"1D{50-k}M"
            md = f"{k}^{REF1[DEL]}{50-k}"
        else:
            seq = REF1[start:start+50]
            cigar, md = "50M", "50"
        tiled.append(dict(name=f"tile{i}", ref="chr1", pos=start,
                          flag=16*(i%2), cigar=cigar, seq=seq,
                          tags=[("MD", md)]))
    bam = make_alignments(first + tiled, refs=REFS, index=True)
    peaks = write_peaks(tmp_path, [("chr1", 90, 180), ("chr1", 1000, 1200)])
    res, out = run(run_macs3, tmp_path, bam, peaks, "on")
    assert res.returncode == 0, res.stderr
    assert "Failed to assemble unitigs" in res.stderr
    expected_peaks = write_peaks(tmp_path, [("chr1", 1000, 1200)], name="later.bed")
    expected, expected_out = run(run_macs3, tmp_path, bam, expected_peaks, "on",
                                 name="expected")
    assert expected.returncode == 0, expected.stderr
    assert records(out) == records(expected_out)
    assert records(out)


def test_short_reads_with_smaller_overlap(tmp_path, make_alignments, run_macs3):
    bam = make_alignments(short_reads(indel=True), refs=REFS, index=True)
    peaks = write_peaks(tmp_path, [("chr1", 1000, 1200)])
    res, out = run(run_macs3, tmp_path, bam, peaks, "on", "--fermi-overlap", 20)
    assert res.returncode == 0, res.stderr
    assert "Failed to assemble unitigs" not in res.stderr
    expected, expected_out = run(run_macs3, tmp_path, bam, peaks, "off",
                                 name="expected")
    assert expected.returncode == 0, expected.stderr
    assert records(out)
    assert [r[:5] for r in records(out)] == [r[:5] for r in records(expected_out)]


@pytest.mark.parametrize("empty", ["bam", "peaks"])
@pytest.mark.parametrize("mode", ["auto", "on", "off"])
def test_empty_input(tmp_path, make_alignments, run_macs3, mode, empty):
    bam = make_alignments([] if empty == "bam" else het_reads(),
                          refs=REFS, index=True)
    peaks = write_peaks(tmp_path, [] if empty == "peaks" else
                        [("chr1", 200, 400)])
    res, out = run(run_macs3, tmp_path, bam, peaks, mode)
    assert res.returncode == 0, res.stderr
    assert out.read_text().startswith("##fileformat=VCFv4.1")
    assert records(out) == []
