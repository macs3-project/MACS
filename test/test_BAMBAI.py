#!/usr/bin/env python
# Time-stamp: <2025-09-29 14:18:42 Tao Liu>

import re
import struct
import unittest
from collections import Counter

import pytest

from MACS3.IO.BAM import BAIFile, BAMaccessor
from MACS3.IO.BAM import (StrandFormatError,
                          MDTagMissingError)


class Test_BAIFile(unittest.TestCase):

    def setUp(self):
        self.bamfile = "test/tiny.bam"
        self.baifile = "test/tiny.bam.bai"
        self.bai = BAIFile(self.baifile)

    def test_load(self):
        expected = {'ref_beg': 50724864, 'ref_end': 2709258240,
                    'n_mapped': 971, 'n_unmapped': 0}
        meta = self.bai.get_metadata_by_refseq(0)
        self.assertDictEqual(expected, meta)

    def test_get_chunks_bin(self):
        expected = [(2178362419, 2178363422)]
        chunks = self.bai.get_chunks_by_bin(0, 4728)
        self.assertListEqual(chunks, expected)

    def test_get_chunks_list_bins(self):
        expected = [(50741975, 50751615),
                    (2178362419, 2178363422),
                    (2178363422, 2709258240)]
        chunks = self.bai.get_chunks_by_list_of_bins(0, [591, 4694, 4728])
        self.assertListEqual(chunks, expected)

    def test_get_chunks_region(self):
        expected = [(1130775424, 1130820164), (1130820164, 1130822629),
                    (1130822629, 2178352280), (2178352280, 2178353973),
                    (2178353973, 2178356807), (2178356807, 2178358139),
                    (2178358139, 2178359954), (2178359954, 2178362419),
                    (2178362419, 2178363422), (2178363422, 2709258240)]
        chunks = self.bai.get_chunks_by_region(0, 557000, 851600)
        self.assertListEqual(chunks, expected)

    def test_get_chunks_list_regions(self):
        expected = [(1130770109, 1130772126), (1130772126, 1130775424),
                    (1130775424, 1130820164), (2178363422, 2709258240)]
        chunks = self.bai.get_chunks_by_list_of_regions(0, [(500000, 600000),
                                                            (800000, 850000)])
        self.assertListEqual(chunks, expected)


class Test_BAM_w_BAI(unittest.TestCase):

    def setUp(self):
        self.bamfile = "test/tiny.bam"
        self.baifile = "test/tiny.bam.bai"
        self.bam = BAMaccessor(self.bamfile)

    def test_get_reads_large_dedup(self):
        c = b"chr10"
        s = 100
        e = 900000
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=1)
        self.assertEqual(len(a), 941)

    def test_get_reads_large(self):
        c = b"chr10"
        s = 100
        e = 900000
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=100)
        self.assertEqual(len(a), 971)

    def test_get_reads_mod(self):
        c = b"chr10"
        s = 100000
        e = 550000
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=100)
        self.assertEqual(len(a), 520)

    def test_get_reads_few(self):
        c = b"chr10"
        s = 142999
        e = 163000
        expected = """chr10	145930	145966	MAGNUM:8:93:1052:1153#0	36	+
chr10	148133	148169	MAGNUM:8:108:19082:3782#0	36	-
chr10	148930	148966	MAGNUM:8:34:18198:13544#0	36	+
chr10	152575	152611	ROCKFORD:1:114:13238:15292#0	36	+
chr10	152927	152963	ROCKFORD:1:21:4055:6592#0	36	-
chr10	153995	154031	MAGNUM:8:44:12716:12267#0	36	+
chr10	156153	156189	MAGNUM:8:74:16381:18650#0	36	-
chr10	159259	159295	ROCKFORD:1:112:5095:4754#0	36	+"""
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=100)
        result = "\n".join(map(str, a))
        self.assertEqual(len(a), 8)
        self.assertEqual(result, expected)

    def test_get_empty1(self):
        # region outside of bams
        c = b"chr10"
        s = 961000
        e = 963000
        a = self.bam.get_reads_in_region(c, s, e)
        self.assertEqual(len(a), 0)

    def test_get_empty2(self):
        # small region missed in bam
        c = b"chr10"
        s = 70000
        e = 90000
        a = self.bam.get_reads_in_region(c, s, e)
        self.assertEqual(len(a), 0)


# ====================================================================
# Comprehensive tests, added below the original ones.
#
# Expected values come from pysam (record coordinates, flags, BGZF
# virtual offsets, index statistics), from a BAI reader written below
# from the SAM specification (section 5.2), and from the spec's reg2bins.
#
# The C-only functions are covered through the Python-visible callers:
#   get_bins_by_region -> BAIFile.get_chunks_by_region,
#       get_chunks_by_list_of_regions, get_coffset_by_region and
#       BAMaccessor.get_reads_in_region (the bin sets are checked exactly
#       with a BAI in which every bin holds one distinct chunk);
#   reg2bins -> no caller in MACS3, so it cannot be reached from Python.
# The private BAMaccessor methods (__parse_header, __check_sorted,
# __decode_voffset, __seek, __retrieve_cdata_from_bgzf_block,
# __fw_binary_parse) are covered through the constructor and
# get_reads_in_region.
# ====================================================================

# ------------------------------------
# helpers
# ------------------------------------

def spec_reg2bins(beg, end):
    """Bins overlapping the 0-based half-open [beg, end) (SAM spec C code)."""
    end -= 1
    bins = [0]
    for first, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        bins.extend(range(first + (beg >> shift), first + (end >> shift) + 1))
    return bins


def read_bai(path):
    """Parse a BAI file following the SAM specification.

    Returns ``(bins, meta)``: ``bins[ref]`` maps bin -> list of
    (begin, end) virtual offsets (pseudo-bin 37450 excluded) and
    ``meta[ref]`` is the pseudo-bin as a dict, or None when the
    reference has no bins.
    """
    with open(path, "rb") as fh:
        data = fh.read()
    assert data[:4] == b"BAI\x01"
    (n_ref,) = struct.unpack_from("<i", data, 4)
    off = 8
    bins, meta = [], []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", data, off)
        off += 4
        ref_bins, ref_meta = {}, None
        for _ in range(n_bin):
            b, n_chunk = struct.unpack_from("<Ii", data, off)
            off += 8
            chunks = [struct.unpack_from("<QQ", data, off + 16 * k)
                      for k in range(n_chunk)]
            off += 16 * n_chunk
            if b == 37450:
                ref_meta = {"ref_beg": chunks[0][0], "ref_end": chunks[0][1],
                            "n_mapped": chunks[1][0],
                            "n_unmapped": chunks[1][1]}
            else:
                ref_bins[b] = chunks
        (n_intv,) = struct.unpack_from("<i", data, off)
        off += 4 + 8 * n_intv
        bins.append(ref_bins)
        meta.append(ref_meta)
    return bins, meta


def bam_records(path):
    """Every record of a BAM, in file order, read with pysam.

    Each item is a dict with the reference id, 0-based start, reference
    end, flag, MAPQ, name, query length, strand, CIGAR, and the BGZF
    virtual offsets where the record begins and ends.
    """
    pysam = pytest.importorskip("pysam")
    out = []
    with pysam.AlignmentFile(path, "rb") as fh:
        v_beg = fh.tell()
        for r in fh.fetch(until_eof=True):
            v_end = fh.tell()
            out.append(dict(ref=r.reference_id, start=r.reference_start,
                            end=r.reference_end, flag=r.flag,
                            mapq=r.mapping_quality, name=r.query_name,
                            qlen=r.query_length, reverse=r.is_reverse,
                            cigar=r.cigartuples, v_beg=v_beg, v_end=v_end))
            v_beg = v_end
    return out


def expected_reads(path, chrom, left, right, max_dup=1):
    """``str()`` of the reads ``get_reads_in_region`` should return.

    pysam's fetch gives the records overlapping [left, right) in file
    order. The accessor's documented filter is applied: unmapped,
    secondary, QC-fail and supplementary records, paired records that
    are not a proper pair or whose mate is unmapped, and MAPQ 0 or 255
    are dropped. At most ``max_dup`` identical alignments (same start,
    end, strand and CIGAR) are kept.
    """
    pysam = pytest.importorskip("pysam")
    seen = Counter()
    out = []
    with pysam.AlignmentFile(path, "rb") as fh:
        for r in fh.fetch(chrom, left, right):
            if r.flag & (4 | 256 | 512 | 2048):
                continue
            if r.flag & 1 and (not r.flag & 2 or r.flag & 8):
                continue
            if r.mapping_quality in (0, 255):
                continue
            key = (r.reference_start, r.reference_end, r.is_reverse,
                   tuple(r.cigartuples))
            seen[key] += 1
            if seen[key] > max_dup:
                continue
            out.append("%s\t%d\t%d\t%s\t%d\t%s" % (
                chrom, r.reference_start, r.reference_end, r.query_name,
                r.query_length, "-" if r.is_reverse else "+"))
    return out


def md_reads(reads):
    """Give every read dict an MD tag (36 matches) unless it has tags."""
    return [dict(r, tags=r.get("tags", [("MD", "36")])) for r in reads]


def names(reads):
    return [r["readname"].decode() for r in reads]


def write_bam_with_header(path, header, reads):
    """Write a BAM with an arbitrary pysam header dict (no sorting)."""
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, pos in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 30
            a.cigarstring = "36M"
            a.query_sequence = "A" * 36
            a.query_qualities = pysam.qualitystring_to_array("I" * 36)
            a.set_tag("MD", "36")
            out.write(a)
    return str(path)


@pytest.fixture(scope="module")
def all_bins_bai(tmp_path_factory):
    """A BAI with one reference in which each bin b of 0..37449 holds the
    single chunk (b << 16, (b << 16) + 1), so a chunk names its bin."""
    path = tmp_path_factory.mktemp("bai") / "all_bins.bai"
    parts = [b"BAI\x01", struct.pack("<ii", 1, 37451)]
    for b in range(37450):
        parts.append(struct.pack("<IiQQ", b, 1, b << 16, (b << 16) + 1))
    # pseudo-bin: (ref_beg, ref_end), (n_mapped, n_unmapped)
    parts.append(struct.pack("<IiQQQQ", 37450, 2, 11, 1 << 20, 7, 3))
    parts.append(struct.pack("<i", 0))      # no linear index
    path.write_bytes(b"".join(parts))
    return BAIFile(str(path))


@pytest.fixture
def tiny_bam(test_dir):
    return str(test_dir / "tiny.bam")


# ------------------------------------
# StrandFormatError, MDTagMissingError
# ------------------------------------

def test_StrandFormatError_attributes_and_message():
    e = StrandFormatError("chr1 1 2", "x")
    assert isinstance(e, Exception)
    assert (e.string, e.strand) == ("chr1 1 2", "x")
    assert str(e) == ("'Strand information can not be recognized in this "
                      "line: \"chr1 1 2\",\"x\"'")


def test_MDTagMissingError_attributes_and_message():
    e = MDTagMissingError(b"read1", b"NM:i:0")
    assert isinstance(e, Exception)
    assert (e.name, e.aux) == (b"read1", b"NM:i:0")
    assert str(e) == ("'MD tag is missing! Please use \"samtools calmd\" "
                      "command to add MD tags!\\nName of sequence:read1\\n"
                      "Current auxiliary data section: NM:i:0'")


# ------------------------------------
# BAIFile: constructor, get_metadata_by_refseq
# ------------------------------------

def test_baifile_rejects_non_bai(tmp_path):
    path = tmp_path / "x.bai"
    path.write_bytes(b"chr1\t1\t2\n")
    with pytest.raises(Exception, match=re.escape(
            "Not a BAI file. The first 4 bytes are 'b'chr1''")):
        BAIFile(str(path))


def test_baifile_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        BAIFile(str(tmp_path / "missing.bai"))


def test_metadata_tiny_matches_spec_reader_and_pysam(tiny_bam):
    pysam = pytest.importorskip("pysam")
    bai = BAIFile(tiny_bam + ".bai")
    _, meta = read_bai(tiny_bam + ".bai")
    with pysam.AlignmentFile(tiny_bam, "rb") as fh:
        stats = fh.get_index_statistics()
    for ref_n, s in enumerate(stats):
        got = bai.get_metadata_by_refseq(ref_n)
        assert got == meta[ref_n]
        if got is not None:
            assert (got["n_mapped"], got["n_unmapped"]) == (s.mapped,
                                                            s.unmapped)


def test_metadata_pysam_bam(make_alignments):
    pysam = pytest.importorskip("pysam")
    refs = (("chr1", 100000), ("chr2", 50000), ("chr3", 1000))
    reads = [dict(name="a%d" % i, ref="chr1", pos=100 * i) for i in range(5)]
    reads += [dict(name="u1", ref="chr1", pos=700, flag=4)]     # placed
    reads += [dict(name="b%d" % i, ref="chr2", pos=10 * i) for i in range(3)]
    reads += [dict(name="u2", ref=None, flag=4, cigar="*")]     # unplaced
    path = make_alignments(reads, refs=refs, index=True)
    bai = BAIFile(path + ".bai")
    with pysam.AlignmentFile(path, "rb") as fh:
        stats = {s.contig: (s.mapped, s.unmapped)
                 for s in fh.get_index_statistics()}
    recs = bam_records(path)
    for ref_n, chrom in enumerate(["chr1", "chr2"]):
        meta = bai.get_metadata_by_refseq(ref_n)
        assert (meta["n_mapped"], meta["n_unmapped"]) == stats[chrom]
        on_ref = [r for r in recs if r["ref"] == ref_n]
        # virtual offsets of the first record's start and last record's end
        assert meta["ref_beg"] == on_ref[0]["v_beg"]
        assert meta["ref_end"] == on_ref[-1]["v_end"]
    assert stats["chr1"] == (5, 1)
    assert stats["chr2"] == (3, 0)
    assert bai.get_metadata_by_refseq(2) is None        # no reads on chr3


def test_metadata_synthetic_bai(all_bins_bai):
    assert all_bins_bai.get_metadata_by_refseq(0) == {
        "ref_beg": 11, "ref_end": 1 << 20, "n_mapped": 7, "n_unmapped": 3}


def test_metadata_reference_out_of_range(all_bins_bai):
    with pytest.raises(KeyError):
        all_bins_bai.get_metadata_by_refseq(1)
    with pytest.raises(OverflowError):
        all_bins_bai.get_metadata_by_refseq(-1)


# ------------------------------------
# BAIFile: get_chunks_by_bin, get_chunks_by_list_of_bins
# ------------------------------------

def test_chunks_by_bin_tiny_match_spec_reader(tiny_bam):
    bai = BAIFile(tiny_bam + ".bai")
    bins, _ = read_bai(tiny_bam + ".bai")
    assert any(bins)
    for ref_n, ref_bins in enumerate(bins):
        for b, chunks in ref_bins.items():
            assert bai.get_chunks_by_bin(ref_n, b) == sorted(chunks)


@pytest.mark.parametrize("b", [0, 1, 9, 73, 585, 4681, 37448])
def test_chunks_by_bin_synthetic(all_bins_bai, b):
    assert all_bins_bai.get_chunks_by_bin(0, b) == [(b << 16, (b << 16) + 1)]


@pytest.mark.parametrize("b", [37450, 37451, 60000])
def test_chunks_by_bin_absent_bin_is_empty(all_bins_bai, b):
    # 37450 is the pseudo-bin, removed when the index is loaded
    assert all_bins_bai.get_chunks_by_bin(0, b) == []


def test_chunks_by_bin_reference_without_bins(make_alignments):
    path = make_alignments([dict(name="a", ref="chr1", pos=10)],
                           refs=(("chr1", 1000), ("chr2", 1000)), index=True)
    bai = BAIFile(path + ".bai")
    assert bai.get_chunks_by_bin(1, 4681) == []
    assert bai.get_chunks_by_bin(0, 4681) != []


def test_chunks_by_bin_reference_out_of_range(all_bins_bai):
    with pytest.raises(IndexError):
        all_bins_bai.get_chunks_by_bin(1, 0)
    with pytest.raises(OverflowError):
        all_bins_bai.get_chunks_by_bin(-1, 0)


@pytest.mark.parametrize("bins", [
    [],
    [4681],
    [4682, 4681, 4681],             # duplicates are used once
    [37449, 0, 585, 60000],         # absent bins are ignored
], ids=["empty", "one", "duplicates", "absent"])
def test_chunks_by_list_of_bins_synthetic(all_bins_bai, bins):
    expected = sorted((b << 16, (b << 16) + 1) for b in set(bins)
                      if b < 37450)
    assert all_bins_bai.get_chunks_by_list_of_bins(0, bins) == expected


def test_chunks_by_list_of_bins_tiny(tiny_bam):
    bai = BAIFile(tiny_bam + ".bai")
    bins, _ = read_bai(tiny_bam + ".bai")
    some = sorted(bins[0])
    expected = sorted(c for b in some for c in bins[0][b])
    assert bai.get_chunks_by_list_of_bins(0, some + some[:3]) == expected


# ------------------------------------
# BAIFile: get_chunks_by_region, get_chunks_by_list_of_regions
# (and the C-only get_bins_by_region)
# ------------------------------------

REGIONS = [
    (0, 1), (0, 16384), (0, 16385), (16383, 16385), (16384, 32768),
    (100000, 250000), (1 << 20, (1 << 20) + 5),
    ((1 << 26) - 10, (1 << 26) + 10), (0, (1 << 29) - 1),
]


@pytest.mark.parametrize("beg,end", REGIONS,
                         ids=["%d-%d" % r for r in REGIONS])
def test_chunks_by_region_bin_set(all_bins_bai, beg, end):
    bins = [c[0] >> 16 for c in all_bins_bai.get_chunks_by_region(0, beg, end)]
    # every bin the spec lists for [beg, end) is used ...
    assert set(spec_reg2bins(beg, end)) <= set(bins)
    # ... and exactly the spec's bins for [beg, end]: ``end`` is inclusive
    assert bins == sorted(set(spec_reg2bins(beg, end + 1)))


def test_chunks_by_region_tiny_match_spec(tiny_bam):
    bai = BAIFile(tiny_bam + ".bai")
    bins, _ = read_bai(tiny_bam + ".bai")
    for beg, end in [(0, 100000), (142999, 163000), (557000, 851600),
                     (900000, 1000000), (5000000, 6000000)]:
        expected = sorted(c for b in set(spec_reg2bins(beg, end + 1))
                          for c in bins[0].get(b, []))
        assert bai.get_chunks_by_region(0, beg, end) == expected


def test_chunks_by_region_cover_every_overlapping_read(make_alignments):
    reads = [dict(name="r%04d" % i, ref="chr1", pos=37 * i,
                  flag=16 * (i % 3 == 0)) for i in range(3000)]
    path = make_alignments(reads, refs=(("chr1", 200000),), index=True)
    bai = BAIFile(path + ".bai")
    recs = bam_records(path)
    assert len({r["v_beg"] >> 16 for r in recs}) > 1    # several BGZF blocks
    for beg, end in [(0, 1000), (20000, 21000), (50000, 80000),
                     (110000, 200000)]:
        chunks = bai.get_chunks_by_region(0, beg, end)
        overlapping = [r for r in recs if r["start"] < end and r["end"] > beg]
        assert overlapping
        for r in overlapping:
            assert any(cb <= r["v_beg"] < ce for cb, ce in chunks), r["name"]


@pytest.mark.parametrize("regions", [
    [],
    [(0, 1)],
    [(0, 16384), (16384, 32768)],
    [(100000, 250000), (1 << 20, (1 << 20) + 5), (0, 1)],
], ids=["empty", "one", "adjacent", "three"])
def test_chunks_by_list_of_regions_synthetic(all_bins_bai, regions):
    bins = set()
    for beg, end in regions:
        bins.update(spec_reg2bins(beg, end + 1))
    expected = sorted((b << 16, (b << 16) + 1) for b in bins)
    assert all_bins_bai.get_chunks_by_list_of_regions(0, regions) == expected


def test_chunks_by_list_of_regions_single_equals_region(tiny_bam):
    bai = BAIFile(tiny_bam + ".bai")
    assert (bai.get_chunks_by_list_of_regions(0, [(557000, 851600)]) ==
            bai.get_chunks_by_region(0, 557000, 851600))


# ------------------------------------
# BAIFile: get_coffset_by_region, get_coffsets_by_list_of_regions
# ------------------------------------

def test_coffset_by_region_is_leftmost_chunk_block(tiny_bam):
    bai = BAIFile(tiny_bam + ".bai")
    recs = [r for r in bam_records(tiny_bam) if r["ref"] == 0]
    for beg, end in [(0, 100000), (142999, 163000), (557000, 851600)]:
        chunks = bai.get_chunks_by_region(0, beg, end)
        coffset = bai.get_coffset_by_region(0, beg, end)
        assert coffset == min(c[0] for c in chunks) >> 16
        first = min(r["v_beg"] for r in recs
                    if r["start"] < end and r["end"] > beg)
        assert coffset <= first >> 16


def test_coffset_by_region_without_chunks_is_zero(make_alignments):
    path = make_alignments([dict(name="a", ref="chr1", pos=10)],
                           refs=(("chr1", 1000), ("chr2", 1000)), index=True)
    bai = BAIFile(path + ".bai")
    assert bai.get_coffset_by_region(1, 0, 1000) == 0
    assert bai.get_coffset_by_region(0, 0, 1000) > 0


# ------------------------------------
# BAMaccessor: constructor, close, get_chromosomes, get_rlengths
# ------------------------------------

def test_accessor_chromosomes_and_lengths_tiny(tiny_bam):
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(tiny_bam, "rb") as fh:
        refs = [r.encode() for r in fh.references]
        lengths = dict(zip(refs, fh.lengths))
    acc = BAMaccessor(tiny_bam)
    assert acc.get_chromosomes() == refs
    assert acc.get_rlengths() == lengths
    acc.close()


def test_accessor_chromosomes_include_empty_ones(make_alignments):
    refs = (("chr1", 100000), ("chr2", 1000), ("chrM", 16569))
    path = make_alignments(md_reads([dict(name="a", ref="chr1", pos=10)]),
                           refs=refs, index=True)
    acc = BAMaccessor(path)
    assert acc.get_chromosomes() == [b"chr1", b"chr2", b"chrM"]
    assert acc.get_rlengths() == {b"chr1": 100000, b"chr2": 1000,
                                  b"chrM": 16569}


def test_accessor_requires_bai(make_alignments):
    path = make_alignments(md_reads([dict(name="a", ref="chr1", pos=10)]))
    msg = ("BAI is not available! Please make sure the `%s.bai` file exists "
           "in the same path" % path)
    with pytest.raises(Exception, match="^" + re.escape(msg) + "$"):
        BAMaccessor(path)


@pytest.mark.parametrize("hd", [{"VN": "1.6", "SO": "unsorted"},
                                {"VN": "1.6", "SO": "queryname"},
                                {"VN": "1.6"}],
                         ids=["unsorted", "queryname", "no-SO"])
def test_accessor_requires_coordinate_sorted_header(tmp_path, hd):
    header = {"HD": hd, "SQ": [{"SN": "chr1", "LN": 1000}]}
    path = write_bam_with_header(tmp_path / "x.bam", header,
                                 [("a", 10), ("b", 20)])
    with pytest.raises(Exception,
                       match=r"^BAM should be sorted by coordinates!$"):
        BAMaccessor(path)


def test_accessor_close_then_query_raises(make_alignments):
    path = make_alignments(md_reads([dict(name="a", ref="chr1", pos=10)]),
                           index=True)
    acc = BAMaccessor(path)
    acc.close()
    with pytest.raises(ValueError):
        acc.get_reads_in_region(b"chr1", 0, 100)


# ------------------------------------
# BAMaccessor.get_reads_in_region
# ------------------------------------

@pytest.mark.parametrize("left,right", [
    (100, 900000), (100000, 550000), (142999, 163000), (557000, 851600),
], ids=lambda x: str(x))
@pytest.mark.parametrize("max_dup", [1, 2, 100])
def test_reads_in_region_tiny_match_pysam(tiny_bam, left, right, max_dup):
    acc = BAMaccessor(tiny_bam)
    got = acc.get_reads_in_region(b"chr10", left, right, maxDuplicate=max_dup)
    assert [str(r) for r in got] == expected_reads(tiny_bam, "chr10", left,
                                                   right, max_dup)


def test_reads_in_region_flag_and_mapq_filter(make_alignments):
    cases = [
        # (name, FLAG, MAPQ, kept)
        ("plus", 0, 30, True),
        ("minus", 16, 30, True),
        ("unmapped", 4, 30, False),
        ("secondary", 256, 30, False),
        ("qcfail", 512, 30, False),
        ("supplementary", 2048, 30, False),
        ("duplicate", 1024, 30, True),
        ("mate1", 99, 30, True),
        ("mate2", 147, 30, True),           # both mates are kept here
        ("notproper", 97, 30, False),
        ("mateunmapped", 75, 30, False),
        ("mapq0", 0, 0, False),
        ("mapq1", 0, 1, True),
        ("mapq254", 0, 254, True),
        ("mapq255", 0, 255, False),
    ]
    reads = md_reads([dict(name=n, ref="chr1", pos=100 + 10 * i, flag=f,
                           mapq=q) for i, (n, f, q, _) in enumerate(cases)])
    path = make_alignments(reads, index=True)
    got = BAMaccessor(path).get_reads_in_region(b"chr1", 0, 1000,
                                                 maxDuplicate=100)
    assert names(got) == [n for n, _, _, kept in cases if kept]
    assert [str(r) for r in got] == expected_reads(path, "chr1", 0, 1000, 100)


def test_reads_in_region_read_fields(make_alignments):
    reads = [dict(name="a", ref="chr1", pos=100, flag=16,
                  cigar="3S10M2D20M3S", tags=[("MD", "10^AC20")]),
             dict(name="b", ref="chr1", pos=105,
                  tags=[("NM", 1), ("MD", "17A18")])]
    a, b = BAMaccessor(make_alignments(reads, index=True)
                       ).get_reads_in_region(b"chr1", 0, 1000)
    # rpos = 100 + 10 (M) + 2 (D) + 20 (M); length is the 36 query bases
    assert str(a) == "chr1\t100\t132\ta\t36\t-"
    assert (a["lpos"], a["rpos"], a["strand"], a["chrom"]) == (100, 132, 1,
                                                               b"chr1")
    # CIGAR ops packed as length << 4 | op, with S = 4, M = 0, D = 2
    assert a["cigar"] == (3 << 4 | 4, 10 << 4 | 0, 2 << 4 | 2, 20 << 4 | 0,
                          3 << 4 | 4)
    assert a["MD"] == b"10^AC20"
    assert str(b) == "chr1\t105\t141\tb\t36\t+"
    assert b["MD"] == b"17A18"


def test_reads_in_region_overlap_boundaries(make_alignments):
    reads = md_reads([
        dict(name="ends_at_left", ref="chr1", pos=64),     # [64, 100)
        dict(name="crosses_left", ref="chr1", pos=80),     # [80, 116)
        dict(name="inside", ref="chr1", pos=150),
        dict(name="crosses_right", ref="chr1", pos=199),   # [199, 235)
        dict(name="after", ref="chr1", pos=300)])
    path = make_alignments(reads, index=True)
    got = BAMaccessor(path).get_reads_in_region(b"chr1", 100, 200)
    assert names(got) == ["crosses_left", "inside", "crosses_right"]
    assert [str(r) for r in got] == expected_reads(path, "chr1", 100, 200)


@pytest.mark.parametrize("max_dup,expected", [
    (0, []),
    (1, ["a1", "b", "m", "c"]),
    (2, ["a1", "a2", "b", "m", "c"]),
    (3, ["a1", "a2", "a3", "b", "m", "c"]),
    (100, ["a1", "a2", "a3", "b", "m", "c"]),
])
def test_reads_in_region_max_duplicate(make_alignments, max_dup, expected):
    reads = md_reads([
        dict(name="a1", ref="chr1", pos=100),
        dict(name="a2", ref="chr1", pos=100),
        dict(name="a3", ref="chr1", pos=100),
        dict(name="b", ref="chr1", pos=100, cigar="30M6S"),    # other CIGAR
        dict(name="m", ref="chr1", pos=100, flag=16),          # other strand
        dict(name="c", ref="chr1", pos=120)])
    path = make_alignments(reads, index=True)
    got = BAMaccessor(path).get_reads_in_region(b"chr1", 0, 1000,
                                                 maxDuplicate=max_dup)
    assert names(got) == expected
    assert [str(r) for r in got] == expected_reads(path, "chr1", 0, 1000,
                                                   max_dup)


def test_reads_in_region_default_max_duplicate_is_one(make_alignments):
    reads = md_reads([dict(name="a1", ref="chr1", pos=100),
                      dict(name="a2", ref="chr1", pos=100)])
    got = BAMaccessor(make_alignments(reads, index=True)
                      ).get_reads_in_region(b"chr1", 0, 1000)
    assert names(got) == ["a1"]


def test_reads_in_region_single_read(make_alignments):
    path = make_alignments(md_reads([dict(name="a", ref="chr1", pos=0)]),
                           index=True)
    got = BAMaccessor(path).get_reads_in_region(b"chr1", 0, 1)
    assert [str(r) for r in got] == ["chr1\t0\t36\ta\t36\t+"]


def test_reads_in_region_empty_cases(make_alignments):
    refs = (("chr1", 1000000), ("chr2", 1000))
    path = make_alignments(md_reads([dict(name="a", ref="chr1", pos=100),
                                     dict(name="b", ref="chr1", pos=200)]),
                           refs=refs, index=True)
    acc = BAMaccessor(path)
    assert acc.get_reads_in_region(b"chr2", 0, 1000) == []      # no reads
    assert acc.get_reads_in_region(b"chr1", 500000, 600000) == []   # beyond
    assert acc.get_reads_in_region(b"chr1", 0, 50) == []        # before
    with pytest.raises(ValueError):
        acc.get_reads_in_region(b"chrX", 0, 1000)               # not in header


def big_bam(make_alignments):
    """3060 records on chr1 over several BGZF blocks: a read every 37 bp,
    every third on the - strand, every 97th with MAPQ 0, and a duplicate
    of every 50th."""
    reads = []
    for i in range(3000):
        r = dict(name="r%04d" % i, ref="chr1", pos=37 * i,
                 flag=16 if i % 3 == 0 else 0,
                 mapq=0 if i % 97 == 0 else 30, tags=[("MD", "36")])
        reads.append(r)
        if i % 50 == 0:
            reads.append(dict(r, name="d%04d" % i, mapq=30))
    return make_alignments(reads, refs=(("chr1", 200000),), index=True)


# right ends are not multiples of 37, so no read starts exactly at them
BIG_REGIONS = [(0, 120000), (0, 1000), (20000, 21000), (50000, 80000),
               (55500, 55501), (110000, 200000)]


def test_reads_in_region_across_bgzf_blocks(make_alignments):
    path = big_bam(make_alignments)
    recs = bam_records(path)
    acc = BAMaccessor(path)
    for left, right in BIG_REGIONS:
        expected = expected_reads(path, "chr1", left, right)
        assert expected
        got = acc.get_reads_in_region(b"chr1", left, right)
        assert [str(r) for r in got] == expected
    # the (50000, 80000) region spans more than one BGZF block
    blocks = {r["v_beg"] >> 16 for r in recs
              if r["start"] < 80000 and r["end"] > 50000}
    assert len(blocks) > 1


def test_reads_in_region_repeated_and_reversed_queries(make_alignments):
    path = big_bam(make_alignments)
    expected = {reg: expected_reads(path, "chr1", *reg, max_dup=100)
                for reg in BIG_REGIONS}
    acc = BAMaccessor(path)
    order = BIG_REGIONS + BIG_REGIONS[::-1] + [BIG_REGIONS[2]] * 2
    for reg in order:
        got = acc.get_reads_in_region(b"chr1", *reg, maxDuplicate=100)
        assert [str(r) for r in got] == expected[reg]


def test_reads_in_region_without_md_tag_raises(make_alignments):
    path = make_alignments([dict(name="r1", ref="chr1", pos=100)], index=True)
    with pytest.raises(MDTagMissingError) as info:
        BAMaccessor(path).get_reads_in_region(b"chr1", 0, 1000)
    assert info.value.name == b"r1\x00"     # read name with its NUL byte
    assert info.value.aux == b""
    assert "MD tag is missing!" in str(info.value)


