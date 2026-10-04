#!/usr/bin/env python

"""Module Description: Test the UnitigRAs and UnitigCollection classes
in MACS3.Signal.UnitigRACollection, which hold reads re-mapped to
fermi-lite unitigs and the unitig-to-reference alignments used by
callvar after local assembly.

The unitig alignment used throughout is written by hand:

    reference  A C G T A - - C G T T G C A    (2000..2011)
    unitig     A C G T A G G C G T T - C A

so the unitig carries a GG insertion after 2004 and lacks the G at
2009. Reads are exact substrings of the unitig sequence ACGTAGGCGTTCA.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import pickle

import pytest

from MACS3.Signal.ReadAlignment import (ReadAlignment)
from MACS3.Signal.UnitigRACollection import (UnitigRAs,
                                             UnitigCollection)

BAM_CODES = "=ACMGRSVTWYHKDBN"
CHROM = b"chr1"
REF_ALN = b"ACGTA--CGTTGCA"
UTG_ALN = b"ACGTAGGCGTT-CA"
UTG_SEQ = "ACGTAGGCGTTCA"


def pack_seq(seq):
    codes = [BAM_CODES.index(c) for c in seq]
    if len(codes) % 2:
        codes.append(0)
    return bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))


def read(seq, qual, strand=0, name="r"):
    """ReadAlignment with SEQ ``seq`` (only SEQ, l, qualities and strand
    matter once a read is assigned to a unitig)."""
    n = len(seq)
    return ReadAlignment(name.encode(), CHROM, 0, n, strand, pack_seq(seq),
                         bytes(qual), ((n << 4) | 0,), str(n).encode())


# treatment reads R1 = unitig[0:8], R2 = unitig[3:11], R3 = unitig[6:13];
# control read C1 = unitig[4:12]
R1 = (UTG_SEQ[0:8], list(range(31, 39)), 0)
R2 = (UTG_SEQ[3:11], list(range(41, 49)), 1)
R3 = (UTG_SEQ[6:13], list(range(21, 28)), 0)
C1 = (UTG_SEQ[4:12], list(range(51, 59)), 1)


def make_ura(control=True):
    T = [read(*R1, name="R1"), read(*R2, name="R2"), read(*R3, name="R3")]
    C = [read(*C1, name="C1")] if control else []
    return UnitigRAs(CHROM, 2000, 2012, UTG_ALN, REF_ALN, [T, C])


def simple_ura(lpos, aln, reads_T, reads_C=()):
    return UnitigRAs(CHROM, lpos, lpos + len(aln), aln, aln,
                     [list(reads_T), list(reads_C)])


# ------------------------------------
# UnitigRAs
# ------------------------------------

def test_UnitigRAs_items():
    ura = make_ura()
    assert ura["chrom"] == CHROM
    assert (ura["lpos"], ura["rpos"]) == (2000, 2012)
    assert ura["seq"] == UTG_SEQ.encode()
    assert ura["unitig_aln"] == UTG_ALN and ura["reference_aln"] == REF_ALN
    assert ura["unitig_length"] == 13
    assert ura["reference_length"] == 12
    assert ura["aln_length"] == 14
    assert ura["count"] == 4


def test_UnitigRAs_unknown_key():
    with pytest.raises(KeyError, match="Unavailable key"):
        make_ura()["RAlists"]


def test_UnitigRAs_alignment_lengths_must_agree():
    with pytest.raises(AssertionError,
                       match="aln on unitig and reference should be the same length!"):
        UnitigRAs(CHROM, 2000, 2012, UTG_ALN + b"A", REF_ALN, [[], []])


@pytest.mark.parametrize("ref_pos,expected", [
    # (allele, bq_T, bq_C, strand_T, strand_C, tip_T, pos_T, pos_C)
    # 2000: first base, R1 position 0 (a tip)
    (2000, (b"A", [31], [], [0], [], [True], [0], [])),
    # 2004: the GG insertion follows, R1 at 4, R2 at 1, C1 at 0
    (2004, (b"AGG", [35, 42], [51], [0, 1], [1], [False, False], [4, 1], [0])),
    # 2005: first base after the insertion; R1 ends here (a tip)
    (2005, (b"C", [38, 45, 22], [54], [0, 1, 0], [1], [True, False, False],
            [7, 4, 1], [3])),
    (2007, (b"T", [47, 24], [56], [1, 0], [1], [False, False], [6, 3], [5])),
    # 2010: after the deletion, only R3 (position 5) and C1 (7) reach it
    (2010, (b"C", [26], [58], [0], [1], [False], [5], [7])),
    # 2011: last base, R3 position 6 (a tip); C1 ends before it
    (2011, (b"A", [27], [], [0], [], [True], [6], [])),
])
def test_UnitigRAs_get_variant_bq_by_ref_pos(ref_pos, expected):
    result = make_ura().get_variant_bq_by_ref_pos(ref_pos)
    assert isinstance(result[0], bytearray)
    assert result == (bytearray(expected[0]),) + expected[1:]


def test_UnitigRAs_get_variant_bq_deleted_base():
    s, bq_t, bq_c, strand_t, strand_c, tip_t, pos_t, pos_c = \
        make_ura().get_variant_bq_by_ref_pos(2009)
    assert s == bytearray(b"*")
    # deleted bases get quality 93; R3 (position 4) and C1 span the deletion
    assert set(bq_t) == {93} and set(bq_c) == {93}
    assert 4 in pos_t


def test_UnitigRAs_without_reads():
    ura = UnitigRAs(CHROM, 2000, 2012, UTG_ALN, REF_ALN, [[], []])
    assert ura.get_variant_bq_by_ref_pos(2004) == (bytearray(b"AGG"), [], [], [],
                                                   [], [], [], [])
    assert ura["count"] == 0


def test_UnitigRAs_pickle_round_trip():
    ura = make_ura()
    ura2 = pickle.loads(pickle.dumps(ura))
    for key in ("chrom", "lpos", "rpos", "seq", "unitig_aln", "reference_aln",
                "unitig_length", "reference_length", "aln_length", "count"):
        assert ura2[key] == ura[key]
    assert (ura2.get_variant_bq_by_ref_pos(2005)
            == ura.get_variant_bq_by_ref_pos(2005))


# ------------------------------------
# UnitigCollection
# ------------------------------------

def three_uras():
    ura1 = make_ura()
    ura2 = simple_ura(2100, b"TTTTGGGGCC",
                      [read("TTGGGG", [30] * 6, 1, "R4")])
    # overlaps ura1 at 2008..2011
    ura3 = simple_ura(2008, b"TGCA", [read("TGCA", [33] * 4, 1, "R5")])
    return ura1, ura2, ura3


def test_UnitigCollection_sort_and_items():
    ura1, ura2, ura3 = three_uras()
    uc = UnitigCollection(CHROM, {"start": 1990, "end": 2120}, [ura2, ura3, ura1])
    assert uc["chrom"] == CHROM
    assert (uc["left"], uc["right"], uc["length"]) == (1990, 2120, 130)
    assert (uc["URAs_left"], uc["URAs_right"]) == (2000, 2110)
    assert uc["count"] == 3
    assert [u["lpos"] for u in uc["URAs_list"]] == [2000, 2008, 2100]
    uc.sort()
    assert [u["lpos"] for u in uc["URAs_list"]] == [2000, 2008, 2100]


def test_UnitigCollection_unknown_key():
    uc = UnitigCollection(CHROM, {"start": 0, "end": 1}, [make_ura()])
    with pytest.raises(KeyError, match="Unavailable key"):
        uc["peak"]


def test_UnitigCollection_empty_list_raises():
    with pytest.raises(IndexError):
        UnitigCollection(CHROM, {"start": 0, "end": 1}, [])


def test_UnitigCollection_get_PosReadsInfo_ref_pos_insertion():
    ura1, ura2, ura3 = three_uras()
    uc = UnitigCollection(CHROM, {"start": 1990, "end": 2120}, [ura1, ura2, ura3])
    pri = uc.get_PosReadsInfo_ref_pos(2004, b"A")
    s = pri.__getstate__()
    n_T, n_C, bq_T, bq_C, strand = s[6], s[7], s[4], s[5], s[9]
    assert (n_T[b"AGG"], n_C[b"AGG"]) == (2, 1)
    assert (bq_T[b"AGG"], bq_C[b"AGG"]) == ([35, 42], [51])
    assert (strand[0][b"AGG"], strand[1][b"AGG"]) == (1, 1)
    assert pri.raw_read_depth() == 3


def test_UnitigCollection_get_PosReadsInfo_ref_pos_quality_cutoff():
    uc = UnitigCollection(CHROM, {"start": 1990, "end": 2120}, [make_ura()])
    pri = uc.get_PosReadsInfo_ref_pos(2004, b"A", Q=40)
    s = pri.__getstate__()
    # R1 has quality 35 at this base, not > 40
    assert (s[6][b"AGG"], s[4][b"AGG"]) == (1, [42])


def test_UnitigCollection_get_PosReadsInfo_ref_pos_combines_unitigs():
    ura1, ura2, ura3 = three_uras()
    uc = UnitigCollection(CHROM, {"start": 1990, "end": 2120}, [ura3, ura1, ura2])
    pri = uc.get_PosReadsInfo_ref_pos(2010, b"C")
    s = pri.__getstate__()
    # ura1: R3 (quality 26), C1 (58); ura3: R5 at position 2 (33)
    assert sorted(s[4][b"C"]) == [26, 33]
    assert s[5][b"C"] == [58]
    assert (s[6][b"C"], s[7][b"C"]) == (2, 1)


def test_UnitigCollection_get_PosReadsInfo_ref_pos_uncovered():
    ura1, ura2, ura3 = three_uras()
    uc = UnitigCollection(CHROM, {"start": 1990, "end": 2120}, [ura1, ura2, ura3])
    assert uc.get_PosReadsInfo_ref_pos(2050, b"A").raw_read_depth() == 0
    # rpos is exclusive
    assert uc.get_PosReadsInfo_ref_pos(2110, b"A").raw_read_depth() == 0
    pri = uc.get_PosReadsInfo_ref_pos(2103, b"T")
    assert pri.__getstate__()[6][b"T"] == 1


def test_UnitigCollection_pickle_round_trip():
    ura1, ura2, ura3 = three_uras()
    uc = UnitigCollection(CHROM, {"start": 1990, "end": 2120}, [ura1, ura2, ura3])
    uc2 = pickle.loads(pickle.dumps(uc))
    assert [u["lpos"] for u in uc2["URAs_list"]] == [2000, 2008, 2100]
    assert (uc2.get_PosReadsInfo_ref_pos(2010, b"C").__getstate__()
            == uc.get_PosReadsInfo_ref_pos(2010, b"C").__getstate__())
