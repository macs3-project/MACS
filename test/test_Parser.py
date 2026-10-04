#!/usr/bin/env python
# Time-stamp: <2025-04-11 13:46:05 Tao Liu>

import gzip
import logging
import re
import unittest

import numpy as np
import pytest

from MACS3.IO.Parser import (guess_parser,
                             StrandFormatError,
                             GenericParser,
                             BEDParser,
                             BEDPEParser,
                             ELANDResultParser,
                             ELANDMultiParser,
                             ELANDExportParser,
                             SAMParser,
                             BAMParser,
                             BAMPEParser,
                             BowtieParser,
                             FragParser)
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import (PETrackI,
                                         PETrackII)
from MACS3.Utilities.Constants import READ_BUFFER_SIZE


class Test_auto_guess(unittest.TestCase):

    def setUp(self):
        self.bedfile = "test/tiny.bed.gz"
        self.bedpefile = "test/tiny.bedpe.gz"
        self.samfile = "test/tiny.sam.gz"
        self.bamfile = "test/tiny.bam"

    def test_guess_parser_bed(self):
        p = guess_parser(self.bedfile)
        self.assertTrue(p.is_gzipped())
        self.assertTrue(isinstance(p, BEDParser))

    def test_guess_parser_sam(self):
        p = guess_parser(self.samfile)
        self.assertTrue(p.is_gzipped())
        self.assertTrue(isinstance(p, SAMParser))

    def test_guess_parser_bam(self):
        p = guess_parser(self.bamfile)
        self.assertTrue(p.is_gzipped())
        self.assertTrue(isinstance(p, BAMParser))


class Test_parsing(unittest.TestCase):
    def setUp(self):
        self.bedfile = "test/tiny.bed.gz"
        self.bedpefile = "test/tiny.bedpe.gz"
        self.samfile = "test/tiny.sam.gz"
        self.bamfile = "test/tiny.bam"
        self.fragfile = "test/tiny.frag.tsv.gz"
        
    def test_fragment_file(self):
        p = FragParser(self.fragfile)
        petrack = p.build_petrack()
        petrack.finalize()


# ====================================================================
# Comprehensive tests, added below the original ones.
#
# Each test writes its own small input under ``tmp_path`` and checks the
# exact reads stored in the resulting track. Expected coordinates are
# derived by hand from each format's definition:
#
# * BED: 0-based, half-open. The 5' end of a + read is ``start``; the
#   5' end of a - read is ``end`` (the position after its last base).
#   Without a strand column the read is taken as +.
# * BEDPE / fragments: (chrom, left, right[, barcode, count]), 0-based.
# * SAM: POS is 1-based. The 5' end of a + read is POS - 1; of a - read
#   it is POS - 1 plus the reference span of the CIGAR (M, D, N, =, X).
# * BAM: pos is 0-based; same rule as SAM.
# * ELAND result / multi / export: 1-based leftmost position. F reads
#   give pos - 1, R reads give pos - 1 + read length.
# * Bowtie: 0-based leftmost offset; - reads give offset + read length.
#
# SAM/BAM single-end filtering, from the parser comments: drop
# unmapped (0x4), secondary (0x100), QC-fail (0x200) and supplementary
# (0x800) records; for paired records (0x1) keep only proper pairs (0x2)
# whose mate is mapped (no 0x8) and that are not the second mate (no
# 0x80). Duplicates (0x400) and MAPQ are not filtered.
#
# C-only methods are covered through the Python-visible callers:
#   skip_first_commentlines (BED, BEDPE, ELAND*, SAM, Frag, Generic) ->
#       every constructor; the header tests below check its effect;
#   tlen_parse_line (all text parsers) -> tsize, sniff, guess_parser;
#   fw_parse_line (all single-end text parsers) -> build_fwtrack,
#       append_fwtrack;
#   pe_parse_line (BEDPE, Frag) -> build_petrack, append_petrack;
#   bam_fw_binary_parse -> BAMParser.build_fwtrack / append_fwtrack;
#   bampe_pe_binary_parse -> BAMPEParser.build_petrack / append_petrack.
# ====================================================================

# ------------------------------------
# helpers
# ------------------------------------

INT_MAX = 2147483647            # length FWTrack/PETrack give unknown chromosomes
SEQ25 = "ACGTACGTACGTACGTACGTACGTA"     # a 25 bp read for ELAND inputs
SAM_HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"


def write_text(path, text, gz=False):
    """Write ``text`` byte for byte (no newline added); gzip if ``gz``."""
    opener = gzip.open if gz else open
    with opener(path, "wb") as fh:
        fh.write(text.encode())
    return str(path)


def gzip_copy(src, dst):
    """Write a gzipped copy of ``src`` to ``dst``; return ``dst``."""
    with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
        fo.write(fi.read())
    return str(dst)


def fw_locations(track):
    """Finalize an FWTrack; return {chrom: ([+ 5' ends], [- 5' ends])}."""
    track.finalize()
    out = {}
    for chrom in track.get_chr_names():
        plus, minus = track.get_locations_by_chr(chrom)
        out[chrom] = (plus.tolist(), minus.tolist())
    return out


def pe_locations(track):
    """Finalize a PETrackI; return {chrom: [(left, right), ...]} sorted."""
    track.finalize()
    return {chrom: [tuple(x) for x in track.get_locations_by_chr(chrom).tolist()]
            for chrom in track.get_chr_names()}


def frag_records(track):
    """Finalize a PETrackII; return {chrom: sorted [(l, r, count, barcode)]}."""
    track.finalize()
    names = {v: k for k, v in track.barcode_dict.items()}
    out = {}
    for chrom in track.get_chr_names():
        locs = track.get_locations_by_chr(chrom).tolist()
        bcs = track.barcodes[chrom].tolist()
        out[chrom] = sorted((l, r, c, names[b])
                            for (l, r, c), b in zip(locs, bcs))
    return out


def parser_messages(caplog):
    """INFO-or-higher records of the Parser logger, memory prefix removed."""
    return [re.sub(r"^\[\d+ MB\] ", "", r.getMessage())
            for r in caplog.records
            if r.name == "MACS3.IO.Parser" and r.levelno >= logging.INFO]


def bed_rows(n=12, length=36, chrom="chr1"):
    """``n`` BED6 rows, alternating + and -, reads of ``length`` bp."""
    return [(chrom, 1000 + 100 * i, 1000 + 100 * i + length, "r%d" % i, 0,
             "+" if i % 2 == 0 else "-") for i in range(n)]


def bed_rows_expected(n=12, length=36, chrom="chr1"):
    """5' ends of ``bed_rows(n, length)``: starts of +, ends of - reads."""
    return {chrom.encode(): ([1000 + 100 * i for i in range(0, n, 2)],
                             [1000 + 100 * i + length
                              for i in range(1, n, 2)])}


def sam_record(name, flag, rname, pos1, cigar="36M", seqlen=36, mapq=30):
    """One SAM alignment line (POS is 1-based)."""
    return "\t".join([name, str(flag), rname, str(pos1), str(mapq), cigar,
                      "*", "0", "0", "A" * seqlen, "I" * seqlen])


def eland_result_line(name, code, chrom=None, pos=None, strand=None,
                      seq=SEQ25):
    """One ELAND result line: name, sequence, match code, counts of
    exact/1-mismatch/2-mismatch hits, then (if mapped) the chromosome
    file, 1-based position, strand and the base-call fields."""
    fields = [">" + name, seq, code, "1", "0", "0"]
    if chrom is not None:
        fields += [chrom, str(pos), strand, "..", ""]
    return "\t".join(fields)


def eland_multi_line(name, counts, hits=None, seq=SEQ25):
    """One ELAND multi line: name, sequence, ``x:y:z`` counts, hits."""
    return "\t".join([name, seq, counts] + ([hits] if hits else []))


def eland_export_line(chrom, pos, strand, seq=SEQ25):
    """One 22-column ELAND export line; ``pos`` is 1-based, None when
    the read is not mapped (then ``chrom`` is NM or QC)."""
    mapped = pos is not None
    fields = ["HWUSI-EAS100", "1", "2", "3", "1000", "2000", "0", "1",
              seq, "I" * len(seq), chrom, "",
              str(pos) if mapped else "", strand if mapped else "",
              "25" if mapped else "", "40" if mapped else "",
              "", "", "", "", "", "Y"]
    return "\t".join(fields)


def bowtie_line(name, strand, chrom, offset, seq="A" * 36):
    """One Bowtie line: name, strand, reference, 0-based offset,
    sequence, qualities, other-instances count, mismatches (empty)."""
    return "\t".join([name, strand, chrom, str(offset), seq,
                      "I" * len(seq), "0", ""])


def pair_reads(name, start, end, flag1=99, flag2=147, ref="chr1",
               readlen=36):
    """pysam read dicts for a pair spanning [start, end).

    The mate with the reverse bit (0x10) sits at ``end - readlen``, the
    other at ``start``; TLEN is +length for the leftmost mate.
    """
    tlen = end - start

    def one(flag, mate_flag):
        rev = bool(flag & 16)
        pos = end - readlen if rev else start
        mpos = end - readlen if mate_flag & 16 else start
        return dict(name=name, ref=ref, pos=pos, flag=flag,
                    cigar="%dM" % readlen, next_ref=ref, next_pos=mpos,
                    tlen=-tlen if rev else tlen)
    return [one(flag1, flag2), one(flag2, flag1)]


# ------------------------------------
# StrandFormatError
# ------------------------------------

def test_StrandFormatError_attributes_and_message():
    e = StrandFormatError("chr1 1 2", "x")
    assert e.string == "chr1 1 2"
    assert e.strand == "x"
    # __str__ is repr() of the sentence, so the outer quotes are included
    assert str(e) == ("'Strand information can not be recognized in this "
                      "line: \"chr1 1 2\",\"x\"'")


def test_StrandFormatError_raise_and_catch():
    with pytest.raises(StrandFormatError) as info:
        raise StrandFormatError(b"line", b"?")
    assert info.value.string == b"line"
    assert info.value.strand == b"?"


# ------------------------------------
# guess_parser
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_guess_parser_detects_bed(write_bed, caplog, gz):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    path = write_bed(bed_rows(), name="reads.bed.gz" if gz else "reads.bed")
    p = guess_parser(path)
    assert type(p) is BEDParser
    assert p.is_gzipped() is gz
    assert parser_messages(caplog) == (["Detected format is: BED"] +
                                       (["* Input file is gzipped."]
                                        if gz else []))
    # the returned parser is positioned at the first read
    assert fw_locations(p.build_fwtrack()) == bed_rows_expected()


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_guess_parser_detects_sam(make_alignments, tmp_path, caplog, gz):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    # forward reads only: reverse ones hit the bug tested further below
    reads = [dict(name="r%d" % i, ref="chr1", pos=100 * i,
                  flag=147 if i % 4 == 0 else 0) for i in range(12)]
    path = make_alignments(reads, name="reads.sam", fmt="sam")
    if gz:
        path = gzip_copy(path, tmp_path / "reads.sam.gz")
    p = guess_parser(path)
    assert type(p) is SAMParser
    assert p.is_gzipped() is gz
    assert parser_messages(caplog) == (["Detected format is: SAM"] +
                                       (["* Input file is gzipped."]
                                        if gz else []))
    assert fw_locations(p.build_fwtrack()) == {
        b"chr1": ([100 * i for i in range(12) if i % 4], [])}


def test_guess_parser_detects_bam(make_alignments, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    reads = [dict(name="r%d" % i, ref="chr1", pos=100 * i,
                  flag=16 if i % 3 == 0 else 0) for i in range(12)]
    p = guess_parser(make_alignments(reads))
    assert type(p) is BAMParser
    assert p.is_gzipped() is True       # BGZF is gzip-compatible
    assert parser_messages(caplog) == ["Detected format is: BAM",
                                       "* Input file is gzipped."]
    assert fw_locations(p.build_fwtrack()) == {
        b"chr1": ([100 * i for i in range(12) if i % 3],
                  [100 * i + 36 for i in range(0, 12, 3)])}


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_guess_parser_detects_eland_result(tmp_path, caplog, gz):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    lines = [eland_result_line("r%d" % i, "U0", "chr1.fa", 1001 + 100 * i,
                               "F") for i in range(12)]
    path = write_text(tmp_path / ("s_1.txt.gz" if gz else "s_1.txt"),
                      "\n".join(lines) + "\n", gz)
    p = guess_parser(path)
    assert type(p) is ELANDResultParser
    assert p.is_gzipped() is gz
    assert parser_messages(caplog)[0] == "Detected format is: ELAND"
    assert fw_locations(p.build_fwtrack()) == {
        b"chr1": ([1000 + 100 * i for i in range(12)], [])}


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_guess_parser_detects_eland_export(tmp_path, caplog, gz):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    lines = [eland_export_line("chr1", 1001 + 100 * i, "R")
             for i in range(12)]
    path = write_text(tmp_path / ("export.txt.gz" if gz else "export.txt"),
                      "\n".join(lines) + "\n", gz)
    p = guess_parser(path)
    assert type(p) is ELANDExportParser
    assert p.is_gzipped() is gz
    assert parser_messages(caplog)[0] == "Detected format is: ELANDEXPORT"
    assert fw_locations(p.build_fwtrack()) == {
        b"chr1": ([], [1025 + 100 * i for i in range(12)])}


@pytest.mark.parametrize("buffer_size", [1, 3, 100000])
def test_guess_parser_buffer_size_keeps_reads(write_bed, buffer_size):
    p = guess_parser(write_bed(bed_rows(15)), buffer_size=buffer_size)
    track = p.build_fwtrack()
    assert track.buffer_size == buffer_size
    assert fw_locations(track) == bed_rows_expected(15)


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_guess_parser_unrecognised_format_raises(write_bed, gz):
    """Every candidate rejects 5 bp BED reads, so guess_parser raises.

    The starts are 4 modulo 8, so read as a SAM FLAG the second column has
    the 'unmapped' bit and SAMParser skips each line; the fifth column is
    one character, so BowtieParser's tag size is 1; BED's is 5. All are
    outside (10, 10000).
    """
    rows = [("chr1", 1004 + 8 * i, 1009 + 8 * i, "r%d" % i, 0, "+")
            for i in range(12)]
    path = write_bed(rows, name="short.bed.gz" if gz else "short.bed")
    with pytest.raises(Exception, match=r"^Can't detect format!$"):
        guess_parser(path)


def test_guess_parser_one_column_text_raises_index_error(tmp_path):
    """Text that is no supported format, with fewer than 3 tab-separated
    fields per line, raises IndexError rather than "Can't detect format!".

    Pins the current output. BEDParser.sniff, tried second, indexes the
    third column of each line and raises before the remaining candidates
    are tried. The input is not in any supported format and the docstring
    only promises an Exception, which IndexError is.
    """
    path = write_text(tmp_path / "notes.txt", "hello world\nsecond line\n")
    with pytest.raises(IndexError):
        guess_parser(path)


# ------------------------------------
# GenericParser: is_gzipped, close, constructor
# ------------------------------------

TEXT_PARSERS = [GenericParser, BEDParser, BEDPEParser, ELANDResultParser,
                ELANDMultiParser, ELANDExportParser, SAMParser, BowtieParser,
                FragParser]


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
@pytest.mark.parametrize("cls", TEXT_PARSERS, ids=lambda c: c.__name__)
def test_is_gzipped(tmp_path, cls, gz):
    path = write_text(tmp_path / ("x.gz" if gz else "x.txt"),
                      "chr1\t100\t136\tr1\t1\t+\n", gz)
    p = cls(path)
    assert p.is_gzipped() is gz
    p.close()


def test_bam_is_gzipped(make_alignments):
    p = BAMParser(make_alignments([dict(name="r", ref="chr1", pos=1)]))
    assert p.is_gzipped() is True
    p.close()


@pytest.mark.parametrize("cls", TEXT_PARSERS + [BAMParser, BAMPEParser],
                         ids=lambda c: c.__name__)
def test_missing_file_raises(tmp_path, cls):
    with pytest.raises(FileNotFoundError):
        cls(str(tmp_path / "missing.txt"))


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_close_closes_the_stream(write_bed, gz):
    p = BEDParser(write_bed(bed_rows(2),
                            name="r.bed.gz" if gz else "r.bed"))
    p.close()
    with pytest.raises(ValueError):
        p.build_fwtrack()


def test_build_fwtrack_closes_the_stream(write_bed):
    p = BEDParser(write_bed(bed_rows(2)))
    p.build_fwtrack()
    with pytest.raises(ValueError):
        p.build_fwtrack()


def test_generic_parser_parses_no_reads(write_bed):
    # the base class's line parser never returns a valid alignment
    track = GenericParser(write_bed(bed_rows(3))).build_fwtrack()
    assert isinstance(track, FWTrack)
    assert fw_locations(track) == {}


def test_generic_parser_append_fwtrack_returns_track_unchanged(write_bed):
    track = FWTrack()
    track.add_loc(b"chr9", 5, 0)
    out = GenericParser(write_bed(bed_rows(3))).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chr9": ([5], [])}


@pytest.mark.parametrize("method", ["tsize", "sniff"])
def test_generic_parser_tsize_and_sniff_not_implemented(write_bed, method):
    p = GenericParser(write_bed(bed_rows(3)))
    with pytest.raises(NotImplementedError):
        getattr(p, method)()


# ------------------------------------
# GenericParser.tsize and sniff (through BEDParser)
# ------------------------------------

@pytest.mark.parametrize("lengths,expected", [
    ([36] * 10, 36),
    ([36] * 9 + [45], 36),          # 369 / 10 = 36.9, truncated
    ([50, 51], 50),                 # 101 / 2 = 50.5, truncated
    ([36] * 10 + [1000], 36),       # only the first 10 usable lines count
    ([100], 100),                   # a single read
], ids=["equal", "truncated-10", "truncated-2", "first-ten", "single"])
def test_bed_tsize_mean_of_first_ten_lengths(write_bed, lengths, expected):
    rows = [("chr1", 100 * i, 100 * i + n, "r", 0, "+")
            for i, n in enumerate(lengths)]
    assert BEDParser(write_bed(rows)).tsize() == expected


def test_bed_tsize_skips_lines_without_positive_length(write_bed):
    rows = [("chr1", 200, 100, "a", 0, "+"), ("chr1", 300, 300, "b", 0, "+"),
            ("chr1", 1000, 1040, "c", 0, "+"), ("chr1", 2000, 2040, "d", 0, "-")]
    assert BEDParser(write_bed(rows)).tsize() == 40


def test_bed_tsize_without_usable_lines_is_minus_one(write_bed):
    rows = [("chr1", 300, 300, "a", 0, "+")] * 3
    assert BEDParser(write_bed(rows)).tsize() == -1


def test_tsize_rewinds_skips_header_and_is_cached(write_bed):
    rows = ["track name=t", "# comment", ("chr1", 10, 46, "a", 0, "+"),
            ("chr1", 20, 56, "b", 0, "-")]
    p = BEDParser(write_bed(rows))
    assert p.tsize() == 36
    assert fw_locations(p.build_fwtrack()) == {b"chr1": ([10], [56])}
    # cached: still answered after build_fwtrack closed the stream
    assert p.tsize() == 36


@pytest.mark.parametrize("lengths,expected", [
    ([10] * 10, False),
    ([11] * 10, True),
    ([9999] * 10, True),
    ([10000] * 10, False),
    ([11] * 9 + [2], False),        # 101 / 10 = 10.1 -> 10
], ids=["10", "11", "9999", "10000", "10.1"])
def test_bed_sniff_tag_size_bounds(write_bed, lengths, expected):
    rows = [("chr1", 100000 * i, 100000 * i + n, "r", 0, "+")
            for i, n in enumerate(lengths)]
    assert BEDParser(write_bed(rows)).sniff() is expected


def test_sniff_true_keeps_every_read(write_bed):
    p = BEDParser(write_bed(["browser hide all"] + bed_rows(12)))
    assert p.sniff() is True
    assert fw_locations(p.build_fwtrack()) == bed_rows_expected(12)


# ------------------------------------
# BEDParser
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_bed_build_fwtrack(write_bed, gz):
    rows = [("chr1", 100, 136, "r1", 0, "+"),
            ("chr1", 200, 236, "r2", 0, "-"),
            ("chr1", 50, 86, "r3", 0, "+"),
            ("chr2", 10, 46, "r4", 0, "-"),
            ("chr2", 0, 36, "r5", 0, "+"),
            ("chr1", 200, 236, "r6", 0, "-")]      # duplicates are kept
    track = BEDParser(write_bed(rows, name="r.bed.gz" if gz else "r.bed")
                      ).build_fwtrack()
    assert isinstance(track, FWTrack)
    assert fw_locations(track) == {b"chr1": ([50, 100], [236, 236]),
                                   b"chr2": ([0], [46])}
    assert track.total == 6


@pytest.mark.parametrize("ncol", [3, 4, 5])
def test_bed_without_strand_column_is_plus(write_bed, ncol):
    rows = [("chr1", 300, 336, "r1", 0)[:ncol],
            ("chr1", 5, 41, "r2", 0)[:ncol]]
    assert fw_locations(BEDParser(write_bed(rows)).build_fwtrack()) == {
        b"chr1": ([5, 300], [])}


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
@pytest.mark.parametrize("header", [
    ["track name=reads"],
    ["browser position chr1:1-1000"],
    ["# comment"],
    ["track name=a", "browser hide all", "#c1", "#c2"],
], ids=["track", "browser", "hash", "mixed"])
def test_bed_header_lines_skipped(write_bed, header, gz):
    rows = header + [("chr1", 100, 136, "r1", 0, "+"),
                     ("chr1", 200, 236, "r2", 0, "-")]
    path = write_bed(rows, name="r.bed.gz" if gz else "r.bed")
    assert fw_locations(BEDParser(path).build_fwtrack()) == {
        b"chr1": ([100], [236])}


def test_bed_crlf_line_endings(tmp_path):
    path = write_text(tmp_path / "crlf.bed", "chr1\t100\t136\tr1\t0\t+\r\n"
                      "chr1\t200\t236\tr2\t0\t-\r\n")
    assert fw_locations(BEDParser(path).build_fwtrack()) == {
        b"chr1": ([100], [236])}


def test_bed_negative_five_prime_end_and_empty_chromosome_skipped(write_bed):
    rows = [("chr1", -5, 31, "a", 0, "+"),     # 5' end -5 < 0: skipped
            ("chr1", -5, 31, "b", 0, "-"),     # 5' end 31: kept
            ("", 100, 136, "c", 0, "+"),       # no chromosome: skipped
            ("chr1", 40, 76, "d", 0, "+")]
    assert fw_locations(BEDParser(write_bed(rows)).build_fwtrack()) == {
        b"chr1": ([40], [31])}


def test_bed_int32_extreme_positions(write_bed):
    rows = [("chr1", 0, 36, "a", 0, "+"),
            ("chr1", 2147483611, 2147483647, "b", 0, "-"),
            ("chr1", 2147483646, 2147483647, "c", 0, "+")]
    assert fw_locations(BEDParser(write_bed(rows)).build_fwtrack()) == {
        b"chr1": ([0, 2147483646], [2147483647])}


def test_bed_single_read(write_bed):
    track = BEDParser(write_bed([("chrX", 7, 8, "a", 0, "-")])).build_fwtrack()
    assert fw_locations(track) == {b"chrX": ([], [8])}
    assert track.total == 1


@pytest.mark.parametrize("strand", [".", "*", "x"])
def test_bed_unknown_strand_raises(write_bed, strand):
    path = write_bed([("chr1", 100, 136, "r1", 0, strand)])
    with pytest.raises(StrandFormatError) as info:
        BEDParser(path).build_fwtrack()
    assert info.value.strand == strand.encode()
    assert info.value.string == ("chr1\t100\t136\tr1\t0\t" + strand).encode()
    assert "Strand information can not be recognized" in str(info.value)


@pytest.mark.parametrize("buffer_size", [1, 2, 7, 100000])
def test_bed_buffer_size_keeps_reads(write_bed, buffer_size):
    rows = bed_rows(25) + [("chr2", 5 * i, 5 * i + 36, "s", 0, "-")
                           for i in range(9)]
    expected = bed_rows_expected(25)
    expected[b"chr2"] = ([], [5 * i + 36 for i in range(9)])
    track = BEDParser(write_bed(rows), buffer_size=buffer_size).build_fwtrack()
    assert track.buffer_size == buffer_size
    assert fw_locations(track) == expected


def test_bed_track_has_unknown_chromosome_lengths(write_bed):
    rows = [("chr1", 1, 37, "a", 0, "+"), ("chr2", 1, 37, "b", 0, "+")]
    track = BEDParser(write_bed(rows)).build_fwtrack()
    assert track.get_rlengths() == {b"chr1": INT_MAX, b"chr2": INT_MAX}


def test_bed_append_fwtrack(write_bed):
    a = write_bed([("chr1", 100, 136, "a", 0, "+")], name="a.bed")
    b = write_bed([("chr1", 50, 86, "b", 0, "-"), ("chr2", 7, 43, "c", 0, "+")],
                  name="b.bed.gz")
    track = BEDParser(a).build_fwtrack()
    out = BEDParser(b).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chr1": ([100], [86]), b"chr2": ([7], [])}


def test_bed_file_larger_than_read_buffer(tmp_path):
    # lines straddle the READ_BUFFER_SIZE (10 MB) blocks the parser reads
    n = 320000
    starts = np.arange(n, dtype=np.int64) * 3
    text = "".join("chr1\t%d\t%d\tread%07d\t0\t%s\n"
                   % (s, s + 36, i, "+" if i % 2 == 0 else "-")
                   for i, s in enumerate(starts.tolist()))
    path = tmp_path / "big.bed"
    path.write_text(text)
    assert path.stat().st_size > READ_BUFFER_SIZE
    track = BEDParser(str(path)).build_fwtrack()
    track.finalize()
    plus, minus = track.get_locations_by_chr(b"chr1")
    np.testing.assert_array_equal(plus, starts[0::2])
    np.testing.assert_array_equal(minus, starts[1::2] + 36)


# ------------------------------------
# BEDPEParser: build_petrack, append_petrack
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_bedpe_build_petrack(write_bedpe, gz):
    rows = [("chr1", 100, 300), ("chr1", 50, 150), ("chr2", 0, 300),
            ("chr1", 100, 300)]                     # duplicates are kept
    p = BEDPEParser(write_bedpe(rows, name="f.bedpe.gz" if gz else "f.bedpe"))
    track = p.build_petrack()
    assert isinstance(track, PETrackI)
    assert pe_locations(track) == {b"chr1": [(50, 150), (100, 300), (100, 300)],
                                   b"chr2": [(0, 300)]}
    assert p.n == 4
    assert p.d == 200.0             # (200 + 100 + 300 + 200) / 4
    assert track.get_rlengths() == {b"chr1": INT_MAX, b"chr2": INT_MAX}


def test_bedpe_extra_columns_ignored(write_bedpe):
    rows = [("chr1", 100, 300, "chr1", 264, 300, "frag1", 60, "+", "-")]
    assert pe_locations(BEDPEParser(write_bedpe(rows)).build_petrack()) == {
        b"chr1": [(100, 300)]}


@pytest.mark.parametrize("header", [
    ["track name=f"],
    ["browser hide all"],
    ["#chrom\tstart\tend"],
    ["track", "#x", "browser y"],
], ids=["track", "browser", "hash", "mixed"])
def test_bedpe_header_lines_skipped(write_bedpe, header):
    path = write_bedpe(header + [("chr1", 10, 60), ("chr1", 5, 25)])
    assert pe_locations(BEDPEParser(path).build_petrack()) == {
        b"chr1": [(5, 25), (10, 60)]}


def test_bedpe_crlf_line_endings(tmp_path):
    path = write_text(tmp_path / "crlf.bedpe", "chr1\t10\t60\r\nchr1\t5\t25\r\n")
    assert pe_locations(BEDPEParser(path).build_petrack()) == {
        b"chr1": [(5, 25), (10, 60)]}


def test_bedpe_negative_left_and_empty_chromosome_skipped(write_bedpe):
    rows = [("chr1", -1, 100), ("", 10, 100), ("chr1", 0, 10)]
    p = BEDPEParser(write_bedpe(rows))
    assert pe_locations(p.build_petrack()) == {b"chr1": [(0, 10)]}
    assert (p.n, p.d) == (1, 10.0)


def test_bedpe_int32_extreme_positions(write_bedpe):
    p = BEDPEParser(write_bedpe([("chr1", 0, 1),
                                 ("chr1", 2147483000, 2147483647)]))
    assert pe_locations(p.build_petrack()) == {
        b"chr1": [(0, 1), (2147483000, 2147483647)]}
    assert p.d == 324.0             # (1 + 647) / 2


@pytest.mark.parametrize("left,right", [(100, 100), (200, 100)],
                         ids=["equal", "reversed"])
def test_bedpe_right_not_after_left_raises(write_bedpe, left, right):
    path = write_bedpe([("chr1", 1, 50), ("chr1", left, right)])
    line = ("chr1\t%d\t%d" % (left, right)).encode()
    msg = ("Right position must be larger than left position, check your "
           "BED file at line: " + repr(line))
    with pytest.raises(AssertionError, match=re.escape(msg)):
        BEDPEParser(path).build_petrack()


@pytest.mark.parametrize("bad", ["chr1\t100", "chr1"],
                         ids=["two-columns", "one-column"])
def test_bedpe_fewer_than_three_columns_raises(write_bedpe, bad):
    path = write_bedpe([("chr1", 1, 100), bad, ("chr1", 5, 50)])
    msg = "Less than 3 columns found at this line: " + repr(bad.encode())
    with pytest.raises(Exception, match=re.escape(msg)):
        BEDPEParser(path).build_petrack()


@pytest.mark.parametrize("buffer_size", [1, 2, 7, 100000])
def test_bedpe_buffer_size_keeps_fragments(write_bedpe, buffer_size):
    rows = [("chr1", 10 * i, 10 * i + 100 + i) for i in range(20)]
    rows += [("chr2", 5, 50)]
    track = BEDPEParser(write_bedpe(rows), buffer_size=buffer_size
                        ).build_petrack()
    assert track.buffer_size == buffer_size
    assert pe_locations(track) == {
        b"chr1": [(10 * i, 10 * i + 100 + i) for i in range(20)],
        b"chr2": [(5, 50)]}


def test_bedpe_append_petrack(write_bedpe):
    a = write_bedpe([("chr1", 0, 100), ("chr1", 50, 250)], name="a.bedpe")
    b = write_bedpe([("chr2", 10, 310)], name="b.bedpe.gz")
    p1 = BEDPEParser(a)
    track = p1.build_petrack()
    assert (p1.n, p1.d) == (2, 150.0)
    p2 = BEDPEParser(b)
    out = p2.append_petrack(track)
    assert out is track
    # a new parser starts from n = 0, so n and d describe its own file
    assert (p2.n, p2.d) == (1, 300.0)
    assert pe_locations(out) == {b"chr1": [(0, 100), (50, 250)],
                                 b"chr2": [(10, 310)]}
    assert out.get_rlengths() == {b"chr1": INT_MAX, b"chr2": INT_MAX}


def test_bedpe_file_larger_than_read_buffer(tmp_path):
    n = 600000
    lefts = np.arange(n, dtype=np.int64) * 2
    text = "".join("chr1\t%d\t%d\n" % (x, x + 150 + (x % 7))
                   for x in lefts.tolist())
    path = tmp_path / "big.bedpe"
    path.write_text(text)
    assert path.stat().st_size > READ_BUFFER_SIZE
    p = BEDPEParser(str(path))
    track = p.build_petrack()
    track.finalize()
    locs = track.get_locations_by_chr(b"chr1")
    np.testing.assert_array_equal(locs["l"], lefts)
    np.testing.assert_array_equal(locs["r"], lefts + 150 + lefts % 7)
    assert p.n == n


# ------------------------------------
# ELANDResultParser
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_eland_result_build_fwtrack(tmp_path, gz):
    lines = ["# ELAND result", "#second comment",
             eland_result_line("r1", "U0", "chr1.fa", 1001, "F"),
             eland_result_line("r2", "U1", "chr1.fa", 2001, "R"),
             eland_result_line("r3", "U2", "chr2", 501, "F"),
             eland_result_line("r4", "NM"),
             eland_result_line("r5", "QC"),
             eland_result_line("r6", "R0"),
             eland_result_line("r7", "R1", "chr1.fa", 3001, "F"),
             ""]
    path = write_text(tmp_path / ("e.txt.gz" if gz else "e.txt"),
                      "\n".join(lines) + "\n", gz)
    # r1: 1001 - 1; r2: 2001 - 1 + 25; r3: 501 - 1; ".fa" is removed;
    # NM, QC and repeat (R*) codes are skipped, as is the blank line
    assert fw_locations(ELANDResultParser(path).build_fwtrack()) == {
        b"chr1": ([1000], [2025]), b"chr2": ([500], [])}


def test_eland_result_tsize_and_sniff(tmp_path):
    lines = [eland_result_line("r%d" % i, "U0", "chr1.fa", 1001 + i, "F")
             for i in range(3)]
    p = ELANDResultParser(write_text(tmp_path / "e.txt",
                                     "\n".join(lines) + "\n"))
    assert p.tsize() == 25
    assert p.sniff() is True


def test_eland_result_unknown_strand_raises(tmp_path):
    line = eland_result_line("r1", "U0", "chr1.fa", 1001, "X")
    path = write_text(tmp_path / "e.txt", line + "\n")
    with pytest.raises(StrandFormatError) as info:
        ELANDResultParser(path).build_fwtrack()
    assert info.value.strand == b"X"
    assert info.value.string == line.rstrip().encode()


def test_eland_result_append_fwtrack(tmp_path, write_bed):
    track = BEDParser(write_bed([("chr1", 1, 37, "a", 0, "+")])).build_fwtrack()
    path = write_text(tmp_path / "e.txt", eland_result_line(
        "r1", "U0", "chr1.fa", 101, "R") + "\n")
    out = ELANDResultParser(path).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chr1": ([1], [125])}


# ------------------------------------
# ELANDMultiParser
# ------------------------------------

def test_eland_multi_lines_without_hits_give_empty_track(tmp_path):
    lines = ["# header", eland_multi_line("r1", "NM"),
             eland_multi_line("r2", "QC"), ""]
    path = write_text(tmp_path / "m.txt", "\n".join(lines) + "\n")
    assert fw_locations(ELANDMultiParser(path).build_fwtrack()) == {}


def test_eland_multi_tsize_and_sniff(tmp_path):
    lines = [eland_multi_line("r%d" % i, "NM") for i in range(3)]
    p = ELANDMultiParser(write_text(tmp_path / "m.txt",
                                    "\n".join(lines) + "\n"))
    assert p.tsize() == 25
    assert p.sniff() is True


# ------------------------------------
# ELANDExportParser
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_eland_export_build_fwtrack(tmp_path, gz):
    lines = ["#export",
             eland_export_line("chr1", 1001, "F"),
             eland_export_line("chr1", 2001, "R"),
             eland_export_line("NM", None, None),
             eland_export_line("QC", None, None),
             eland_export_line("chr2", 1, "F"),
             ""]
    path = write_text(tmp_path / ("x.txt.gz" if gz else "x.txt"),
                      "\n".join(lines) + "\n", gz)
    # unmapped reads have an empty position column and are skipped
    assert fw_locations(ELANDExportParser(path).build_fwtrack()) == {
        b"chr1": ([1000], [2025]), b"chr2": ([0], [])}


def test_eland_export_tsize_skips_unmapped(tmp_path):
    lines = [eland_export_line("NM", None, None, seq="A" * 50)] * 2
    lines += [eland_export_line("chr1", 11 + i, "F") for i in range(3)]
    p = ELANDExportParser(write_text(tmp_path / "x.txt",
                                     "\n".join(lines) + "\n"))
    assert p.tsize() == 25
    assert p.sniff() is True


def test_eland_export_unknown_strand_raises(tmp_path):
    line = eland_export_line("chr1", 1001, "X")
    path = write_text(tmp_path / "x.txt", line + "\n")
    with pytest.raises(StrandFormatError) as info:
        ELANDExportParser(path).build_fwtrack()
    assert info.value.strand == b"X"
    assert info.value.string == line.encode()


def test_eland_export_append_fwtrack(tmp_path):
    track = FWTrack()
    track.add_loc(b"chr1", 3, 1)
    path = write_text(tmp_path / "x.txt",
                      eland_export_line("chr1", 11, "F") + "\n")
    out = ELANDExportParser(path).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chr1": ([10], [3])}


# ------------------------------------
# SAMParser and BAMParser: shared single-end flag and CIGAR rules
# ------------------------------------

SE_CASES = [
    # (FLAG, CIGAR, expected (strand, 5' end) for a read at 0-based 1000,
    #  or None when the record is filtered out)
    (0, "36M", (0, 1000)),
    (16, "36M", (1, 1036)),
    (16, "5S31M", (1, 1031)),           # soft clips do not cover reference
    (16, "31M5S", (1, 1031)),
    (16, "5H31M", (1, 1031)),           # hard clip
    (16, "10M2D26M", (1, 1038)),        # deletion covers 2 reference bases
    (16, "10M2I24M", (1, 1034)),        # insertion covers none
    (16, "10M100N26M", (1, 1136)),      # skipped region (N)
    (16, "10=1X25=", (1, 1036)),        # sequence match / mismatch
    (0, "5S31M", (0, 1000)),            # + reads start at the first aligned base
    (0, "10M100N26M", (0, 1000)),
    (4, "36M", None),                   # unmapped
    (256, "36M", None),                 # secondary
    (272, "36M", None),                 # secondary, reverse
    (512, "36M", None),                 # QC failure
    (2048, "36M", None),                # supplementary
    (2064, "36M", None),                # supplementary, reverse
    (1024, "36M", (0, 1000)),           # duplicates are kept
    (1040, "36M", (1, 1036)),
    (99, "36M", (0, 1000)),             # first mate of a proper pair
    (83, "36M", (1, 1036)),             # first mate, reverse
    (1123, "36M", (0, 1000)),           # first mate, duplicate
    (147, "36M", None),                 # second mate
    (163, "36M", None),                 # second mate, forward
    (65, "36M", None),                  # paired, not proper
    (97, "36M", None),                  # paired, not proper, mate reverse
    (73, "36M", None),                  # mate unmapped
    (75, "36M", None),                  # proper flag set but mate unmapped
]


def se_params():
    # SAM reverse-strand reads are left out: SAMParser raises TypeError on
    # them in this version.
    params = []
    for fmt in ("sam", "bam"):
        for flag, cigar, expected in SE_CASES:
            if fmt == "sam" and expected is not None and expected[0] == 1:
                continue
            params.append(pytest.param(fmt, flag, cigar, expected,
                                       id="%s-%d-%s" % (fmt, flag, cigar)))
    return params


@pytest.mark.parametrize("fmt,flag,cigar,expected", se_params())
def test_se_flag_and_cigar(make_alignments, fmt, flag, cigar, expected):
    path = make_alignments([dict(name="r1", ref="chr1", pos=1000, flag=flag,
                                 cigar=cigar)],
                           name="reads." + fmt, fmt=fmt)
    parser = SAMParser(path) if fmt == "sam" else BAMParser(path)
    locs = fw_locations(parser.build_fwtrack())
    if expected is None:
        assert locs == {}
    else:
        strand, pos = expected
        assert locs == {b"chr1": ([pos], []) if strand == 0 else ([], [pos])}


@pytest.mark.parametrize("fmt", ["sam", "bam"])
@pytest.mark.parametrize("mapq", [0, 1, 255])
def test_se_mapq_is_not_filtered(make_alignments, fmt, mapq):
    path = make_alignments([dict(name="r1", ref="chr1", pos=10, mapq=mapq),
                            dict(name="r2", ref="chr1", pos=20, flag=1024,
                                 mapq=mapq)],
                           name="reads." + fmt, fmt=fmt)
    parser = SAMParser(path) if fmt == "sam" else BAMParser(path)
    assert fw_locations(parser.build_fwtrack()) == {b"chr1": ([10, 20], [])}


def test_sam_and_bam_give_identical_tracks(make_alignments):
    refs = (("chr1", 10000), ("chr2", 10000))
    # kept reads are forward: SAMParser raises TypeError on reverse reads
    # in this version
    reads = [dict(name="a", ref="chr1", pos=10),
             dict(name="b", ref="chr1", pos=20, cigar="4S30M2S"),
             dict(name="c", ref="chr2", pos=5, flag=99, cigar="20M3D16M"),
             dict(name="d", ref="chr2", pos=7, flag=147),
             dict(name="f", ref="chr2", pos=9, flag=272),
             dict(name="e", ref=None, flag=4, cigar="*")]
    sam = make_alignments(reads, refs=refs, name="r.sam", fmt="sam")
    bam = make_alignments(reads, refs=refs, name="r.bam")
    expected = {b"chr1": ([10, 20], []), b"chr2": ([5], [])}
    assert fw_locations(SAMParser(sam).build_fwtrack()) == expected
    assert fw_locations(BAMParser(bam).build_fwtrack()) == expected


# ------------------------------------
# SAMParser
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_sam_build_fwtrack_from_text(tmp_path, gz):
    text = SAM_HEADER + "@CO\tfree text\n" + "\n".join([
        sam_record("r1", 0, "chr1", 101),
        sam_record("r2", 0, "chr1", 201, "10M5D26M"),
        sam_record("r3", 0, "chr2.fa", 1),
        sam_record("r4", 4, "*", 0, "*"),
        sam_record("r5", 0, "chr1", 301, "3S30M3S"),
        sam_record("r6", 2048, "chr1", 401)]) + "\n"
    path = write_text(tmp_path / ("r.sam.gz" if gz else "r.sam"), text, gz)
    # POS - 1; ".fa" is removed from chr2.fa; r4 (unmapped) and r6
    # (supplementary) are skipped
    assert fw_locations(SAMParser(path).build_fwtrack()) == {
        b"chr1": ([100, 200, 300], []), b"chr2": ([0], [])}


def test_sam_blank_lines_and_crlf(tmp_path):
    text = (SAM_HEADER.replace("\n", "\r\n") +
            sam_record("r1", 0, "chr1", 101) + "\r\n\r\n" +
            sam_record("r2", 0, "chr1", 201) + "\r\n\n")
    path = write_text(tmp_path / "r.sam", text)
    assert fw_locations(SAMParser(path).build_fwtrack()) == {
        b"chr1": ([100, 200], [])}


def test_sam_tsize_skips_filtered_reads(tmp_path):
    lines = [sam_record("u%d" % i, 4, "*", 0, "*", seqlen=50)
             for i in range(3)]
    lines.append(sam_record("s", 256, "chr1", 50, "60M", seqlen=60))
    lines.append(sam_record("m2", 147, "chr1", 60, "70M", seqlen=70))
    lines += [sam_record("r%d" % i, 0, "chr1", 101 + i) for i in range(10)]
    p = SAMParser(write_text(tmp_path / "r.sam",
                             SAM_HEADER + "\n".join(lines) + "\n"))
    assert p.tsize() == 36
    assert p.sniff() is True


def test_sam_too_few_columns_raises_index_error(tmp_path):
    """A record with 4 columns (no CIGAR) cannot be parsed.

    Pins the current output. The SAM spec requires 11 columns but says
    nothing about how a reader fails; the original raises IndexError.
    """
    path = write_text(tmp_path / "r.sam", SAM_HEADER + "r1\t0\tchr1\t101\n")
    with pytest.raises(IndexError):
        SAMParser(path).build_fwtrack()


@pytest.mark.parametrize("buffer_size", [1, 2, 100000])
def test_sam_buffer_size_keeps_reads(make_alignments, buffer_size):
    refs = (("chr1", 100000), ("chr2", 100000))
    reads = [dict(name="r%d" % i, ref="chr1" if i % 2 else "chr2",
                  pos=10 * i) for i in range(15)]
    path = make_alignments(reads, refs=refs, name="r.sam", fmt="sam")
    track = SAMParser(path, buffer_size=buffer_size).build_fwtrack()
    assert track.buffer_size == buffer_size
    assert fw_locations(track) == {
        b"chr1": ([10 * i for i in range(1, 15, 2)], []),
        b"chr2": ([10 * i for i in range(0, 15, 2)], [])}


def test_sam_append_fwtrack(make_alignments, write_bed):
    track = BEDParser(write_bed([("chr1", 5, 41, "a", 0, "+")])).build_fwtrack()
    sam = make_alignments([dict(name="r", ref="chr1", pos=1000),
                           dict(name="s", ref="chr2", pos=7, flag=163)],
                          refs=(("chr1", 10000), ("chr2", 10000)),
                          name="r.sam", fmt="sam")
    out = SAMParser(sam).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chr1": ([5, 1000], [])}


# ------------------------------------
# BAMParser: sniff, tsize, get_references, build_fwtrack, append_fwtrack
# ------------------------------------

def test_bam_get_references_lists_every_header_chromosome(make_alignments):
    refs = (("chr1", 5000), ("chr2", 3000), ("chrM", 16569))
    p = BAMParser(make_alignments([dict(name="r", ref="chr1", pos=10)],
                                  refs=refs))
    assert p.get_references() == ([b"chr1", b"chr2", b"chrM"],
                                  {b"chr1": 5000, b"chr2": 3000,
                                   b"chrM": 16569})


def test_bam_header_chromosomes_without_reads(make_alignments, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    refs = (("chr1", 5000), ("chr2", 3000), ("chr3", 7000))
    reads = [dict(name="a", ref="chr1", pos=10),
             dict(name="b", ref="chr3", pos=20, flag=16),
             dict(name="u", ref=None, flag=4, cigar="*")]    # unplaced
    track = BAMParser(make_alignments(reads, refs=refs)).build_fwtrack()
    assert isinstance(track, FWTrack)
    assert fw_locations(track) == {b"chr1": ([10], []), b"chr3": ([], [56])}
    # lengths are kept only for chromosomes that received reads
    assert track.get_rlengths() == {b"chr1": 5000, b"chr3": 7000}
    assert parser_messages(caplog) == ["2 reads have been read."]


def test_bam_without_alignments_gives_empty_track(make_alignments, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    track = BAMParser(make_alignments([])).build_fwtrack()
    assert fw_locations(track) == {}
    assert parser_messages(caplog) == ["0 reads have been read."]


def test_bam_positions_near_int32_max(make_alignments):
    refs = (("chr1", 2147483647),)
    reads = [dict(name="a", ref="chr1", pos=0),
             dict(name="b", ref="chr1", pos=2147483611, flag=16)]
    track = BAMParser(make_alignments(reads, refs=refs)).build_fwtrack()
    assert fw_locations(track) == {b"chr1": ([0], [2147483647])}
    assert track.get_rlengths() == {b"chr1": 2147483647}


def test_bam_tsize_mean_of_first_ten_records(make_alignments):
    lengths = [36] * 9 + [45] + [100, 100]
    reads = [dict(name="r%d" % i, ref="chr1", pos=1000 * i, cigar="%dM" % n)
             for i, n in enumerate(lengths)]
    p = BAMParser(make_alignments(reads))
    # (9 * 36 + 45) / 10 = 36.9, truncated; records 11 and 12 are not read
    assert p.tsize() == 36
    p.build_fwtrack()
    assert p.tsize() == 36          # cached after the stream was closed


def test_bam_sniff_true_keeps_every_read(make_alignments):
    reads = [dict(name="r%d" % i, ref="chr1", pos=100 * i,
                  flag=16 * (i % 2)) for i in range(12)]
    p = BAMParser(make_alignments(reads))
    assert p.sniff() is True
    assert fw_locations(p.build_fwtrack()) == {
        b"chr1": ([100 * i for i in range(0, 12, 2)],
                  [100 * i + 36 for i in range(1, 12, 2)])}


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_bam_sniff_false_on_text(write_bed, gz):
    p = BAMParser(write_bed(bed_rows(), name="r.bed.gz" if gz else "r.bed"))
    assert p.sniff() is False
    assert p.is_gzipped() is gz
    p.close()


def test_bam_sniff_raises_when_records_have_no_sequence(make_alignments):
    reads = [dict(name="r%d" % i, ref="chr1", pos=100 * i, seq="", qual="")
             for i in range(10)]
    with pytest.raises(Exception,
                       match=r"^File is not of a valid BAM format! 0$"):
        BAMParser(make_alignments(reads)).sniff()


@pytest.mark.parametrize("buffer_size", [1, 2, 7, 100000])
def test_bam_buffer_size_keeps_reads(make_alignments, buffer_size):
    refs = (("chr1", 100000), ("chr2", 100000))
    reads = [dict(name="r%d" % i, ref="chr1" if i < 15 else "chr2",
                  pos=10 * i, flag=16 * (i % 2)) for i in range(20)]
    track = BAMParser(make_alignments(reads, refs=refs),
                      buffer_size=buffer_size).build_fwtrack()
    assert track.buffer_size == buffer_size
    assert fw_locations(track) == {
        b"chr1": ([10 * i for i in range(0, 15, 2)],
                  [10 * i + 36 for i in range(1, 15, 2)]),
        b"chr2": ([10 * i for i in range(16, 20, 2)],
                  [10 * i + 36 for i in range(15, 20, 2)])}


def test_bam_append_fwtrack_to_bed_track(write_bed, make_alignments, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    track = BEDParser(write_bed([("chrB", 5, 41, "a", 0, "+")])).build_fwtrack()
    bam = make_alignments([dict(name="r", ref="chr1", pos=1000, flag=16),
                           dict(name="s", ref="chr1", pos=10),
                           dict(name="t", ref="chr1", pos=20, flag=256)],
                          refs=(("chr1", 5000), ("chr2", 100)))
    out = BAMParser(bam).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chrB": ([5], []),
                                 b"chr1": ([10], [1036])}
    assert out.get_rlengths() == {b"chrB": INT_MAX, b"chr1": 5000}
    assert parser_messages(caplog) == ["2 reads have been read."]


def test_bam_append_takes_lengths_from_appended_header(make_alignments):
    """Appending a BAM whose header lacks chr2 sets chr2's length to
    INT_MAX, although the first BAM gave it as 3000.

    append_fwtrack passes only the appended file's header lengths to
    FWTrack.set_rlengths, whose docstring states that a chromosome missing
    from that mapping is assigned INT_MAX.
    """
    a = make_alignments([dict(name="a", ref="chr1", pos=10),
                         dict(name="b", ref="chr2", pos=10)],
                        refs=(("chr1", 5000), ("chr2", 3000)), name="a.bam")
    b = make_alignments([dict(name="c", ref="chr1", pos=20)],
                        refs=(("chr1", 5000),), name="b.bam")
    track = BAMParser(a).build_fwtrack()
    assert track.get_rlengths() == {b"chr1": 5000, b"chr2": 3000}
    BAMParser(b).append_fwtrack(track)
    assert track.get_rlengths() == {b"chr1": 5000, b"chr2": INT_MAX}


# ------------------------------------
# BAMPEParser: build_petrack, append_petrack
# ------------------------------------

PE_CASES = [
    # (R1 FLAG, R2 FLAG, fragments expected from a pair spanning [1000, 1200))
    (99, 147, [(1000, 1200)]),      # R1 forward, R2 reverse
    (83, 163, [(1000, 1200)]),      # R1 reverse: left end is its mate's start
    (1123, 1171, [(1000, 1200)]),   # duplicates (0x400) are kept
    (97, 145, []),                  # paired but not a proper pair
    (355, 403, []),                 # secondary (0x100)
    (611, 659, []),                 # QC failure (0x200)
    (2147, 2195, []),               # supplementary (0x800)
    (75, 133, []),                  # mate unmapped (0x8); the mate itself (0x4)
]


@pytest.mark.parametrize("flag1,flag2,expected", PE_CASES,
                         ids=["%d-%d" % (a, b) for a, b, _ in PE_CASES])
def test_bampe_flag_filtering(make_alignments, flag1, flag2, expected):
    # a valid pair on chr2 keeps the track non-empty
    reads = (pair_reads("p", 1000, 1200, flag1, flag2) +
             pair_reads("anchor", 10, 310, ref="chr2"))
    p = BAMPEParser(make_alignments(reads, refs=(("chr1", 10000),
                                                 ("chr2", 10000))))
    locs = pe_locations(p.build_petrack())
    assert locs.get(b"chr1", []) == expected
    assert locs[b"chr2"] == [(10, 310)]
    assert p.n == 1 + len(expected)


def test_bampe_build_petrack(make_alignments, caplog):
    caplog.set_level(logging.INFO, logger="MACS3.IO.Parser")
    refs = (("chr1", 5000), ("chr2", 6000), ("chr3", 7000))
    reads = (pair_reads("a", 100, 300) + pair_reads("b", 150, 250, 83, 163) +
             pair_reads("c", 0, 300, ref="chr2"))
    p = BAMPEParser(make_alignments(reads, refs=refs))
    track = p.build_petrack()
    assert isinstance(track, PETrackI)
    assert pe_locations(track) == {b"chr1": [(100, 300), (150, 250)],
                                   b"chr2": [(0, 300)]}
    assert p.n == 3
    assert p.d == 200.0             # (200 + 100 + 300) / 3
    assert track.get_rlengths() == {b"chr1": 5000, b"chr2": 6000}
    assert parser_messages(caplog) == ["3 fragments have been read."]


def test_bampe_sniff_tsize_and_references(make_alignments):
    reads = []
    for i in range(5):
        reads += pair_reads("p%d" % i, 1000 * i, 1000 * i + 200)
    p = BAMPEParser(make_alignments(reads, refs=(("chr1", 9000),
                                                 ("chr2", 10))))
    assert p.sniff() is True
    assert p.tsize() == 36
    assert p.get_references() == ([b"chr1", b"chr2"],
                                  {b"chr1": 9000, b"chr2": 10})


@pytest.mark.parametrize("buffer_size", [1, 2, 7, 100000])
def test_bampe_buffer_size_keeps_fragments(make_alignments, buffer_size):
    reads = []
    for i in range(12):
        reads += pair_reads("p%d" % i, 100 * i, 100 * i + 150 + i)
    track = BAMPEParser(make_alignments(reads), buffer_size=buffer_size
                        ).build_petrack()
    assert track.buffer_size == buffer_size
    assert pe_locations(track) == {
        b"chr1": [(100 * i, 100 * i + 150 + i) for i in range(12)]}


def test_bampe_append_petrack(make_alignments):
    a = make_alignments(pair_reads("a", 100, 300) + pair_reads("b", 400, 500),
                        name="a.bam")
    b = make_alignments(pair_reads("c", 50, 350, ref="chr2"),
                        refs=(("chr1", 100000), ("chr2", 1000)), name="b.bam")
    p1 = BAMPEParser(a)
    track = p1.build_petrack()
    assert (p1.n, p1.d) == (2, 150.0)
    p2 = BAMPEParser(b)
    out = p2.append_petrack(track)
    assert out is track
    # a new parser starts from n = 0, so n and d describe its own file
    assert (p2.n, p2.d) == (1, 300.0)
    assert pe_locations(out) == {b"chr1": [(100, 300), (400, 500)],
                                 b"chr2": [(50, 350)]}
    assert out.get_rlengths() == {b"chr1": 100000, b"chr2": 1000}


# ------------------------------------
# BowtieParser
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_bowtie_build_fwtrack(tmp_path, gz):
    lines = [bowtie_line("r1", "+", "chr1", 100),
             bowtie_line("r2", "-", "chr1", 200),
             bowtie_line("r3", "-", "chr2", 0, "A" * 50),
             "",
             bowtie_line("r4", "+", "chr2", 7)]
    path = write_text(tmp_path / ("b.map.gz" if gz else "b.map"),
                      "\n".join(lines) + "\n", gz)
    # - reads: offset + read length
    assert fw_locations(BowtieParser(path).build_fwtrack()) == {
        b"chr1": ([100], [236]), b"chr2": ([7], [50])}


def test_bowtie_empty_file_gives_empty_track(tmp_path):
    path = write_text(tmp_path / "empty.map", "")
    assert fw_locations(BowtieParser(path).build_fwtrack()) == {}


def test_bowtie_tsize_and_sniff(tmp_path):
    lines = [bowtie_line("r%d" % i, "+", "chr1", 100 * i) for i in range(12)]
    p = BowtieParser(write_text(tmp_path / "b.map", "\n".join(lines) + "\n"))
    assert p.tsize() == 36
    assert p.sniff() is True


def test_bowtie_unknown_strand_raises(tmp_path):
    line = bowtie_line("r1", "*", "chr1", 100)
    path = write_text(tmp_path / "b.map", line + "\n")
    with pytest.raises(StrandFormatError) as info:
        BowtieParser(path).build_fwtrack()
    assert info.value.strand == b"*"
    assert info.value.string == line.rstrip().encode()


def test_bowtie_append_fwtrack(tmp_path):
    track = FWTrack()
    track.add_loc(b"chr1", 1, 0)
    path = write_text(tmp_path / "b.map",
                      bowtie_line("r1", "-", "chr1", 10) + "\n")
    out = BowtieParser(path).append_fwtrack(track)
    assert out is track
    assert fw_locations(out) == {b"chr1": ([1], [46])}


def test_bowtie_comment_lines_skipped(tmp_path):
    lines = ["# bowtie output", bowtie_line("r1", "+", "chr1", 100)]
    path = write_text(tmp_path / "b.map", "\n".join(lines) + "\n")
    assert fw_locations(BowtieParser(path).build_fwtrack()) == {
        b"chr1": ([100], [])}


# ------------------------------------
# FragParser: build_petrack, append_petrack, max_count, barcodes
# ------------------------------------

@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_frag_build_petrack(write_frag, gz):
    rows = ["# id=sample", "# pipeline=cellranger-atac",
            ("chr1", 10, 50, "AAAC-1", 1),
            ("chr1", 100, 160, "AAAG-1", 2),
            ("chr2", 0, 30, "AAAC-1", 5),
            ("chr1", 5, 25, "AAAT-1", 1)]
    p = FragParser(write_frag(rows, name="f.tsv.gz" if gz else "f.tsv"))
    track = p.build_petrack()
    assert isinstance(track, PETrackII)
    assert frag_records(track) == {
        b"chr1": [(5, 25, 1, b"AAAT-1"), (10, 50, 1, b"AAAC-1"),
                  (100, 160, 2, b"AAAG-1")],
        b"chr2": [(0, 30, 5, b"AAAC-1")]}
    assert p.n == 4
    assert p.d == 37.5      # (40 + 60 + 30 + 20) / 4: per line, not per count
    assert track.get_rlengths() == {b"chr1": INT_MAX, b"chr2": INT_MAX}


def test_frag_extra_columns_ignored(write_frag):
    path = write_frag([("chr1", 10, 50, "A", 2, "extra", 9)])
    assert frag_records(FragParser(path).build_petrack()) == {
        b"chr1": [(10, 50, 2, b"A")]}


@pytest.mark.parametrize("max_count,counts", [
    (0, [1, 5, 10]),            # 0 means no cap
    (3, [1, 3, 3]),
    (1, [1, 1, 1]),
    (100, [1, 5, 10]),
])
def test_frag_max_count_caps_counts(write_frag, max_count, counts):
    rows = [("chr1", 10, 50, "A", 1), ("chr1", 20, 60, "B", 5),
            ("chr1", 30, 70, "C", 10)]
    track = FragParser(write_frag(rows)).build_petrack(max_count=max_count)
    assert frag_records(track) == {
        b"chr1": [(10, 50, counts[0], b"A"), (20, 60, counts[1], b"B"),
                  (30, 70, counts[2], b"C")]}


@pytest.mark.parametrize("max_count,expected_count", [(0, 9), (2, 2)])
def test_frag_append_petrack(write_frag, max_count, expected_count):
    a = write_frag([("chr1", 10, 50, "A", 1)], name="a.tsv")
    b = write_frag([("chr1", 0, 100, "B", 9), ("chr3", 5, 15, "A", 1)],
                   name="b.tsv.gz")
    p1 = FragParser(a)
    track = p1.build_petrack()
    p2 = FragParser(b)
    out = p2.append_petrack(track, max_count=max_count)
    assert out is track
    assert frag_records(out) == {
        b"chr1": [(0, 100, expected_count, b"B"), (10, 50, 1, b"A")],
        b"chr3": [(5, 15, 1, b"A")]}
    assert (p1.n, p1.d) == (1, 40.0)
    # a new parser starts from n = 0, so n and d describe its own file
    assert (p2.n, p2.d) == (2, 55.0)        # (100 + 10) / 2
    assert out.get_rlengths() == {b"chr1": INT_MAX, b"chr3": INT_MAX}


def test_frag_barcode_subset(write_frag):
    rows = [("chr1", 10, 50, "A", 1), ("chr1", 20, 60, "B", 2),
            ("chr2", 0, 30, "B", 1), ("chr2", 5, 35, "C", 3)]
    track = FragParser(write_frag(rows)).build_petrack()
    track.finalize()
    sub = track.subset({b"B", b"not-in-file"})
    assert frag_records(sub) == {b"chr1": [(20, 60, 2, b"B")],
                                 b"chr2": [(0, 30, 1, b"B")]}


def test_frag_count_65535_is_kept(write_frag):
    path = write_frag([("chr1", 10, 50, "A", 65535)])
    assert frag_records(FragParser(path).build_petrack()) == {
        b"chr1": [(10, 50, 65535, b"A")]}


def test_frag_negative_left_and_empty_chromosome_skipped(write_frag):
    rows = [("chr1", -3, 50, "A", 1), ("", 10, 50, "A", 1),
            ("chr1", 0, 8, "B", 2)]
    p = FragParser(write_frag(rows))
    assert frag_records(p.build_petrack()) == {b"chr1": [(0, 8, 2, b"B")]}
    assert (p.n, p.d) == (1, 8.0)


@pytest.mark.parametrize("left,right", [(100, 100), (200, 100)],
                         ids=["equal", "reversed"])
def test_frag_right_not_after_left_raises(write_frag, left, right):
    path = write_frag([("chr1", 1, 50, "A", 1), ("chr1", left, right, "B", 1)])
    line = ("chr1\t%d\t%d\tB\t1" % (left, right)).encode()
    msg = ("Right position must be larger than left position, check your "
           "BED file at line: " + repr(line))
    with pytest.raises(AssertionError, match=re.escape(msg)):
        FragParser(path).build_petrack()


@pytest.mark.parametrize("bad", ["chr1\t10\t50\tA", "chr1\t10\t50"],
                         ids=["four-columns", "three-columns"])
def test_frag_fewer_than_five_columns_raises(write_frag, bad):
    path = write_frag([("chr1", 1, 100, "A", 1), bad])
    msg = "Less than 5 columns found at this line: " + repr(bad.encode())
    with pytest.raises(Exception, match=re.escape(msg)):
        FragParser(path).build_petrack()


@pytest.mark.parametrize("buffer_size", [1, 2, 7, 100000])
def test_frag_buffer_size_keeps_fragments(write_frag, buffer_size):
    rows = [("chr1", 10 * i, 10 * i + 40 + i, "BC%d" % (i % 3), 1 + i % 4)
            for i in range(20)]
    track = FragParser(write_frag(rows), buffer_size=buffer_size
                       ).build_petrack()
    assert track.buffer_size == buffer_size
    assert frag_records(track) == {
        b"chr1": [(10 * i, 10 * i + 40 + i, 1 + i % 4, b"BC%d" % (i % 3))
                  for i in range(20)]}


# ------------------------------------
# files without a trailing newline (parsed in a child interpreter only)
# ------------------------------------


