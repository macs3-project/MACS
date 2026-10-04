#!/usr/bin/env python

"""Module Description: Test the ReadAlignment class in
MACS3.Signal.ReadAlignment, which holds one BAM alignment (CIGAR, packed
sequence, base qualities and MD tag) for callvar.

Reads are built the way BAMaccessor builds them: the sequence packed
4 bits per base in BAM order "=ACMGRSVTWYHKDBN", raw Phred qualities,
CIGAR operations encoded as length << 4 | op ("MIDNSHP=X" -> 0..8) and
the MD tag as bytes. Expected reference sequences and per-position bases
are derived by hand from the reference string REF below.

The C-only methods get_n_edits and relative_ref_pos_to_relative_query_pos
are covered through the ``n_edits`` item and get_base_by_ref_pos,
get_bq_by_ref_pos and get_base_bq_by_ref_pos.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import pickle
import re

import pytest

from MACS3.Signal.ReadAlignment import (ReadAlignment)
from MACS3.IO.BAM import (BAMaccessor,
                          MDTagMissingError)

BAM_CODES = "=ACMGRSVTWYHKDBN"
CIGAR_OPS = "MIDNSHP=X"

# reference sequence at genome positions 1000..1033
#        0         1         2         3
#        0123456789012345678901234567890123
REF = "ACGTTGCAAGCTGACCTGATCGGTACATGCAAGT"
OFFSET = 1000


# ------------------------------------
# helpers
# ------------------------------------

def pack_seq(seq):
    """BAM 4-bit packing; an odd-length read ends with a 0 ('=') nibble."""
    codes = [BAM_CODES.index(c) for c in seq]
    if len(codes) % 2:
        codes.append(0)
    return bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))


def cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def encode_cigar(cigar):
    return tuple((n << 4) | CIGAR_OPS.index(op) for n, op in cigar_ops(cigar))


def ref_span(cigar):
    return sum(n for n, op in cigar_ops(cigar) if op in "MDN=X")


def md_tag(cigar, seq, start):
    """MD tag of ``seq`` aligned at REF[start:] (independent check of the
    hand-written MD strings)."""
    md, run, r, q = "", 0, start, 0
    for n, op in cigar_ops(cigar):
        if op in "M=X":
            for _ in range(n):
                if seq[q] == REF[r]:
                    run += 1
                else:
                    md += "%d%s" % (run, REF[r])
                    run = 0
                q += 1
                r += 1
        elif op == "D":
            md += "%d^%s" % (run, REF[r:r + n])
            run = 0
            r += n
        elif op in "IS":
            q += n
        elif op == "N":
            r += n
    return md + str(run)


def make_ra(start, cigar, seq, md, qual=None, strand=0, name="r", chrom=b"chr1"):
    """ReadAlignment at genome position OFFSET + start."""
    if qual is None:
        qual = [30] * len(seq)
    lpos = OFFSET + start
    return ReadAlignment(name.encode(), chrom, lpos, lpos + ref_span(cigar),
                         strand, pack_seq(seq), bytes(qual),
                         encode_cigar(cigar), md.encode())


# (id, start, cigar, read sequence, MD, expected reference, n_edits)
CASES = [
    ("match", 2, "10M", "GTTGCAAGCT", "10", REF[2:12], 0),
    # read base 4 is T where REF[6] is C
    ("mismatch", 2, "10M", "GTTGTAAGCT", "4C5", REF[2:12], 1),
    # REF[6:8] = CA deleted from the read
    ("deletion", 2, "4M2D6M", "GTTGAGCTGA", "4^CA6", REF[2:14], 2),
    # GG inserted after REF[5]
    ("insertion", 2, "4M2I6M", "GTTGGGCAAGCT", "10", REF[2:12], 2),
    # two and one soft-clipped bases, odd read length
    ("softclip", 2, "2S8M1S", "TTGTTGCAAGA", "8", REF[2:10], 3),
    # hard-clipped bases are not in SEQ
    ("hardclip", 2, "3H8M2H", "GTTGCAAG", "8", REF[2:10], 0),
    # 5 skipped reference bases: the MD tag covers aligned bases only, so
    # the reference rebuilt from it is REF[2:5] + REF[10:14]
    ("skip", 2, "3M5N4M", "GTTCTGA", "7", REF[2:5] + REF[10:14], 0),
    # S2 M3(TGC) I1(A) M2(A, T for A) D1(G) M3(CTG) S1
    ("combined", 4, "2S3M1I2M1D3M1S", "CCTGCAATCTGG", "4A0^G3", REF[4:13], 6),
    # mismatches at the first and last aligned base
    ("mismatch_ends", 0, "6M", "TCGTTA", "0A4G0", REF[0:6], 2),
    # deletion followed by a mismatch
    ("del_then_mismatch", 0, "3M2D3M", "ACGGTA", "3^TT1C1", REF[0:8], 3),
    ("single_base", 33, "1M", "T", "1", REF[33], 0),
    # sequence match / mismatch operators
    ("eq_x", 2, "4=1X5=", "GTTGTAAGCT", "4C5", REF[2:12], 1),
]
CASE_IDS = [c[0] for c in CASES]


def case_ra(case_id, **kw):
    _, start, cigar, seq, md, _, _ = CASES[CASE_IDS.index(case_id)]
    return make_ra(start, cigar, seq, md, **kw)


# ------------------------------------
# construction and accessors
# ------------------------------------

@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_hand_written_md_tags(case):
    _, start, cigar, seq, md, _, _ = case
    assert md_tag(cigar, seq, start) == md


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_construction_fields(case):
    _, start, cigar, seq, md, _, n_edits = case
    qual = [20 + i % 20 for i in range(len(seq))]
    ra = make_ra(start, cigar, seq, md, qual=qual, strand=1, name="read1")
    assert ra["readname"] == b"read1"
    assert ra["chrom"] == b"chr1"
    assert ra["lpos"] == OFFSET + start
    assert ra["rpos"] == OFFSET + start + ref_span(cigar)
    assert ra["strand"] == 1
    assert ra["SEQ"] == seq.encode()
    assert ra["QUAL"] == bytes(qual)
    assert ra["l"] == len(seq)
    assert ra["binaryseq"] == pack_seq(seq)
    assert ra["binaryqual"] == bytes(qual)
    assert ra["cigar"] == encode_cigar(cigar)
    assert ra["MD"] == md.encode()
    # n_edits = inserted + soft-clipped bases + letters in MD
    assert ra["n_edits"] == n_edits


def test_n_edits_counts_mismatch_deletion_insertion_softclip():
    # 1 softclip + 2 inserted + 1 mismatch + 3 deleted bases
    seq = "A" + REF[0:3] + "GG" + REF[3:5] + "T" + REF[9:11]
    ra = make_ra(0, "1S3M2I3M3D2M", seq, md_tag("1S3M2I3M3D2M", seq, 0))
    assert ra["n_edits"] == 1 + 2 + 1 + 3


def test_getitem_unknown_key():
    with pytest.raises(KeyError, match="No such key"):
        case_ra("match")["n_edit"]


@pytest.mark.parametrize("strand,s", [(0, "+"), (1, "-")])
def test_str(strand, s):
    ra = case_ra("softclip", strand=strand)
    assert str(ra) == "chr1\t1002\t1010\tr\t11\t%s" % s


def test_seq_and_qual_lengths_must_agree():
    with pytest.raises(AssertionError,
                       match="Lengths of seq and qual are not consistent!"):
        ReadAlignment(b"r", b"chr1", 0, 4, 0, pack_seq("ACGT"), bytes([30] * 3),
                      encode_cigar("4M"), b"4")


def test_odd_length_drops_padding_nibble():
    ra = make_ra(0, "3M", "ACG", "3")
    assert ra["SEQ"] == b"ACG"
    assert len(ra["binaryseq"]) == 2


def test_int32_positions():
    lpos = 2**31 - 11
    ra = ReadAlignment(b"r", b"chr1", lpos, lpos + 10, 0, pack_seq("ACGTACGTAC"),
                       bytes([30] * 10), encode_cigar("10M"), b"10")
    assert ra["rpos"] == 2**31 - 1
    # 2**31 - 2 is the last aligned base, query index 9 of ACGTACGTAC
    assert ra.get_base_by_ref_pos(2**31 - 2) == ord("C")
    assert ra.get_variant_bq_by_ref_pos(2**31 - 2) == (bytearray(b"C"),
                                                       bytearray([30]), 0, True, 9)
    with pytest.raises(OverflowError):
        ReadAlignment(b"r", b"chr1", 2**31, 2**31 + 1, 0, pack_seq("A"),
                      bytes([30]), encode_cigar("1M"), b"1")


def test_state_and_pickle_round_trip():
    ra = case_ra("combined", strand=1)
    clone = pickle.loads(pickle.dumps(ra))
    assert clone.__getstate__() == ra.__getstate__()
    assert str(clone) == str(ra)
    assert clone.get_REFSEQ() == ra.get_REFSEQ()


# ------------------------------------
# get_FASTQ
# ------------------------------------

def test_get_FASTQ_forward():
    ra = case_ra("softclip", qual=list(range(30, 41)), name="q1")
    # Phred + 33: 30..40 -> '?'..'I'
    assert ra.get_FASTQ() == b"@q1\nTTGTTGCAAGA\n+\n?@ABCDEFGHI\n"


def test_get_FASTQ_reverse_strand_is_reverse_complemented():
    ra = case_ra("softclip", qual=list(range(30, 41)), name="q1", strand=1)
    # reverse of TTGTTGCAAGA is AGAACGTTGTT, complement TCTTGCAACAA
    assert ra.get_FASTQ() == b"@q1\nTCTTGCAACAA\n+\nIHGFEDCBA@?\n"


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_get_FASTQ_reverse_matches_python_revcomp(case):
    _, start, cigar, seq, md, _, _ = case
    qual = [2 + 3 * i for i in range(len(seq))]
    ra = make_ra(start, cigar, seq, md, qual=qual, strand=1)
    revcomp = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    qtext = "".join(chr(q + 33) for q in qual[::-1])
    assert ra.get_FASTQ() == ("@r\n%s\n+\n%s\n" % (revcomp, qtext)).encode()


# ------------------------------------
# get_REFSEQ
# ------------------------------------

@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_get_REFSEQ(case):
    _, start, cigar, seq, md, expected, _ = case
    refseq = make_ra(start, cigar, seq, md).get_REFSEQ()
    assert isinstance(refseq, bytearray)
    assert refseq == expected.encode()


@pytest.mark.parametrize("case", [c for c in CASES if c[0] != "skip"],
                         ids=[c for c in CASE_IDS if c != "skip"])
def test_get_REFSEQ_length_is_reference_span(case):
    _, start, cigar, seq, md, _, _ = case
    ra = make_ra(start, cigar, seq, md)
    assert len(ra.get_REFSEQ()) == ra["rpos"] - ra["lpos"]


def test_get_REFSEQ_does_not_change_SEQ():
    ra = case_ra("combined")
    ra.get_REFSEQ()
    assert ra["SEQ"] == b"CCTGCAATCTGG"


def test_get_REFSEQ_rejects_unknown_md_character():
    ra = make_ra(2, "10M", "GTTGCAAGCT", "5a4")
    with pytest.raises(Exception, match="Don't understand this operator in MD: a"):
        ra.get_REFSEQ()


# ------------------------------------
# get_base_by_ref_pos / get_bq_by_ref_pos / get_base_bq_by_ref_pos
# ------------------------------------

# "combined" read at 1004: SEQ CC TGC A AT - CTG G, qualities 20 + index
COMBINED_QUAL = list(range(20, 32))
COMBINED_POS = [
    # (ref_pos, query index or None for the deleted base)
    (1004, 2), (1005, 3), (1006, 4), (1007, 6), (1008, 7), (1009, None),
    (1010, 8), (1011, 9), (1012, 10),
]


@pytest.mark.parametrize("ref_pos,qi", COMBINED_POS)
def test_base_and_bq_by_ref_pos(ref_pos, qi):
    ra = case_ra("combined", qual=COMBINED_QUAL)
    seq = "CCTGCAATCTGG"
    if qi is None:
        assert ra.get_base_by_ref_pos(ref_pos) is None
        assert ra.get_bq_by_ref_pos(ref_pos) is None
        assert ra.get_base_bq_by_ref_pos(ref_pos) is None
    else:
        # the base comes back as its byte value
        assert ra.get_base_by_ref_pos(ref_pos) == ord(seq[qi])
        assert ra.get_bq_by_ref_pos(ref_pos) == COMBINED_QUAL[qi]
        assert ra.get_base_bq_by_ref_pos(ref_pos) == (ord(seq[qi]),
                                                     COMBINED_QUAL[qi])


@pytest.mark.parametrize("case_id,ref_pos,expected", [
    ("match", 1002, "G"),
    ("match", 1011, "T"),
    ("mismatch", 1006, "T"),        # the read base, not the reference
    ("insertion", 1005, "G"),
    ("insertion", 1006, "C"),       # first base after the insertion
    ("softclip", 1002, "G"),        # soft-clipped bases are skipped
    ("hardclip", 1009, "G"),
    ("skip", 1004, "T"),
    ("skip", 1006, None),           # inside the N gap
    ("skip", 1010, "C"),
    ("deletion", 1007, None),
    ("deletion", 1008, "A"),
    ("eq_x", 1006, "T"),
])
def test_get_base_by_ref_pos_cases(case_id, ref_pos, expected):
    result = case_ra(case_id).get_base_by_ref_pos(ref_pos)
    assert result == (None if expected is None else ord(expected))


@pytest.mark.parametrize("ref_pos", [1003, 1013, -1, 2**40])
@pytest.mark.parametrize("method", ["get_base_by_ref_pos", "get_bq_by_ref_pos",
                                    "get_base_bq_by_ref_pos",
                                    "get_variant_bq_by_ref_pos"])
def test_ref_pos_outside_alignment(method, ref_pos):
    ra = case_ra("combined")             # covers [1004, 1013)
    with pytest.raises(AssertionError,
                       match="Given position out of alignment location"):
        getattr(ra, method)(ref_pos)


# ------------------------------------
# get_variant_bq_by_ref_pos
# ------------------------------------

@pytest.mark.parametrize("ref_pos,alleles,bqs,pos", [
    (1004, b"T", [22], 2),
    (1005, b"G", [23], 3),
    (1006, b"CA", [24, 25], 4),     # last base before I1 carries the insertion
    (1007, b"A", [26], 6),
    (1008, b"T", [27], 7),          # SNV: read T, reference A
    (1009, b"*", [93], 8),          # deleted base: '*' with quality 93
    (1010, b"C", [28], 8),
    (1012, b"G", [30], 10),
])
def test_get_variant_bq_by_ref_pos_combined(ref_pos, alleles, bqs, pos):
    ra = case_ra("combined", qual=COMBINED_QUAL, strand=1)
    result = ra.get_variant_bq_by_ref_pos(ref_pos)
    assert result == (bytearray(alleles), bytearray(bqs), 1, False, pos)
    assert isinstance(result[0], bytearray) and isinstance(result[1], bytearray)


@pytest.mark.parametrize("case_id,ref_pos,expected", [
    # tips: the first and last query base of the read
    ("match", 1002, (b"G", [30], 0, True, 0)),
    ("match", 1011, (b"T", [39], 0, True, 9)),
    ("match", 1006, (b"C", [34], 0, False, 4)),
    ("insertion", 1005, (b"GGG", [33, 34, 35], 0, False, 3)),
    ("insertion", 1006, (b"C", [36], 0, False, 6)),
    ("deletion", 1005, (b"G", [33], 0, False, 3)),
    ("deletion", 1006, (b"*", [93], 0, False, 4)),
    ("deletion", 1007, (b"*", [93], 0, False, 4)),
    ("softclip", 1002, (b"G", [32], 0, False, 2)),
    ("softclip", 1009, (b"G", [39], 0, False, 9)),
    ("hardclip", 1002, (b"G", [30], 0, True, 0)),
    ("hardclip", 1009, (b"G", [37], 0, True, 7)),
    ("skip", 1006, (b"*", [93], 0, False, 3)),     # N is treated like D
    ("skip", 1010, (b"C", [33], 0, False, 3)),
    ("single_base", 1033, (b"T", [30], 0, True, 0)),
    ("eq_x", 1006, (b"T", [34], 0, False, 4)),
])
def test_get_variant_bq_by_ref_pos_cases(case_id, ref_pos, expected):
    seq = CASES[CASE_IDS.index(case_id)][3]
    ra = case_ra(case_id, qual=[30 + i for i in range(len(seq))])
    alleles, bqs, strand, tip, pos = expected
    assert ra.get_variant_bq_by_ref_pos(ref_pos) == (bytearray(alleles),
                                                     bytearray(bqs), strand,
                                                     tip, pos)


def test_get_variant_bq_insertion_at_read_end():
    # 5M2I: the inserted bases follow the last aligned base
    ra = make_ra(0, "5M2I", REF[0:5] + "GG", "5", qual=[30, 31, 32, 33, 34, 35, 36])
    assert ra.get_variant_bq_by_ref_pos(1004) == (bytearray(b"TGG"),
                                                  bytearray([34, 35, 36]),
                                                  0, False, 4)


# ------------------------------------
# reads built by BAMaccessor
# ------------------------------------

def write_case_bam(make_alignments, with_md=True, extra_tags=()):
    reads = []
    for i, (cid, start, cigar, seq, md, _, _) in enumerate(CASES):
        qual = "".join(chr(33 + 20 + (j + i) % 20) for j in range(len(seq)))
        tags = ([("MD", md)] if with_md else []) + list(extra_tags)
        reads.append(dict(name=cid, ref="chr1", pos=OFFSET + start,
                          flag=16 if i % 2 else 0, cigar=cigar, seq=seq,
                          qual=qual, tags=tags))
    return make_alignments(reads, refs=(("chr1", 5000),), index=True)


def test_BAMaccessor_builds_the_same_ReadAlignment(make_alignments):
    bam = write_case_bam(make_alignments)
    got = BAMaccessor(bam).get_reads_in_region(b"chr1", 900, 1100,
                                               maxDuplicate=100)
    assert len(got) == len(CASES)
    by_name = {ra["readname"].decode(): ra for ra in got}
    for i, (cid, start, cigar, seq, md, refseq, n_edits) in enumerate(CASES):
        qual = [20 + (j + i) % 20 for j in range(len(seq))]
        expected = make_ra(start, cigar, seq, md, qual=qual,
                           strand=1 if i % 2 else 0, name=cid)
        assert by_name[cid].__getstate__() == expected.__getstate__()
        assert by_name[cid].get_REFSEQ() == refseq.encode()


def test_BAMaccessor_without_MD_raises(make_alignments):
    bam = write_case_bam(make_alignments, with_md=False,
                         extra_tags=[("XX", "abc")])
    with pytest.raises(MDTagMissingError) as excinfo:
        BAMaccessor(bam).get_reads_in_region(b"chr1", 900, 1100)
    message = str(excinfo.value)
    assert "MD tag is missing! Please use \"samtools calmd\" command to add MD tags!" in message
    # the first read in coordinate order is "mismatch_ends" at 1000
    assert excinfo.value.name.rstrip(b"\x00") == b"mismatch_ends"
    assert "Name of sequence:mismatch_ends" in message
    assert excinfo.value.aux == b"XXZabc\x00"
    assert "Current auxiliary data section: XXZabc" in message
