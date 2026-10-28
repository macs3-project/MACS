"""Regression coverage for singleton-strand accounting (issues #805 and #757)."""

from collections import Counter

import pytest

from MACS3.Signal.FixWidthTrack import FWTrack


@pytest.mark.parametrize("maxnum", [-1, 1, 2])
@pytest.mark.parametrize("reads", [
    [],
    [(b"chr1", 100, 0)],
    [(b"chr1", 100, 1)],
    [(b"chr1", 100, 0), (b"chr1", 150, 1)],
    [(b"chr1", 100, 0), (b"chr2", 100, 1)],
    [(b"chr1", p, 0) for p in (300, 100, 200)] + [(b"chr2", 100, 1)],
    [(b"chr1", 100, 0)] * 3 + [(b"chr2", 100, 1)],
])
def test_filter_dup_singleton_totals(reads, maxnum):
    track = FWTrack(fw=25)
    for read in reads:
        track.add_loc(*read)
    track.finalize()
    expected = Counter(reads)
    if maxnum >= 0:
        expected = Counter({read: min(count, maxnum)
                            for read, count in expected.items()})
    # A second pass also covers strands reduced to singletons by the first.
    for _ in range(2):
        assert track.filter_dup(maxnum) == sum(expected.values())
        assert track.total == sum(expected.values())
        assert track.length == 25 * sum(expected.values())
        actual = Counter()
        for chrom in track.get_chr_names():
            for strand, positions in enumerate(track.get_locations_by_chr(chrom)):
                actual.update((chrom, int(pos), strand) for pos in positions)
        assert actual == expected


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("keepdup", ["1", "auto", "all"])
def test_callpeak_singleton_strands(tmp_path, write_bed, run_macs3,
                                   strand, keepdup):
    path = write_bed([
        ("chr1", 20051074, 20051098, "r1", 0, strand),
        ("chr2", 11326107, 11326132, "r2", 0, strand),
    ])
    result = run_macs3([
        "callpeak", "-f", "BED", "-t", path, "-g", "hs", "--nomodel",
        "--extsize", "200", "--keep-dup", keepdup, "-n", "singleton",
    ], timeout=30)
    assert result.returncode == 0, result.stderr
    assert "total tags in treatment: 2" in result.stderr
    if keepdup != "all":
        assert "tags after filtering in treatment: 2" in result.stderr
        assert "Redundant rate of treatment: 0.00" in result.stderr
    assert (tmp_path / "singleton_peaks.xls").is_file()
    assert (tmp_path / "singleton_peaks.narrowPeak").is_file()


@pytest.mark.parametrize("duplicate", [False, True])
def test_filterdup_mixed_singletons(tmp_path, write_bed, run_macs3, duplicate):
    rows = [("chr1", p, p + 25, "r", 0, "+") for p in (100, 200, 300)]
    rows += [("chr2", 100, 125, "r", 0, "-")]
    if duplicate:
        rows.append(rows[0])
    path = write_bed(rows)
    result = run_macs3([
        "filterdup", "-f", "BED", "-i", path, "-g", "hs",
        "--keep-dup", "1", "-o", "filtered.bed",
    ], timeout=30)
    assert result.returncode == 0, result.stderr
    assert "tags after filtering in alignment file: 4" in result.stderr
    rate = "0.20" if duplicate else "0.00"
    assert f"Redundant rate of alignment file: {rate}" in result.stderr
    output = [line.split("\t") for line in
              (tmp_path / "filtered.bed").read_text().splitlines()]
    assert sorted((r[0], int(r[1]), int(r[2]), r[5]) for r in output) == [
        ("chr1", 100, 125, "+"), ("chr1", 200, 225, "+"),
        ("chr1", 300, 325, "+"), ("chr2", 100, 125, "-"),
    ]


@pytest.mark.parametrize("command", ["callpeak", "filterdup"])
@pytest.mark.parametrize("content, error", [
    ("", "No alignments found"), ("chr1\n", "IndexError"),
])
def test_duplicate_filter_invalid_input(tmp_path, run_macs3,
                                       command, content, error):
    path = tmp_path / "invalid.bed"
    path.write_text(content)
    args = [command, "-f", "BED", "-g", "hs", "--keep-dup", "1"]
    if command == "callpeak":
        args += ["-t", path, "--nomodel", "--extsize", "200", "-n", "invalid"]
    else:
        args += ["-i", path, "-o", "invalid-out.bed"]
    result = run_macs3(args, timeout=30)
    assert result.returncode > 0, result.stderr
    if not content and command == "callpeak":
        # Existing empty-input handling fails while calculating redundancy.
        assert "ZeroDivisionError" in result.stderr
    else:
        assert error in result.stderr
