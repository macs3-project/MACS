#!/usr/bin/env python

"""Module Description: Test the RACollection class in
MACS3.Signal.RACollection, which holds the read alignments overlapping
one peak, rebuilds the peak reference sequence from their MD tags,
collects per-position allele information and runs the fermi-lite local
assembly.

Reads are ReadAlignment objects built as BAMaccessor builds them, from
a deterministic pseudo-random reference REF placed at genome position
G0. Expected reference sequences and allele counts are derived by hand
from REF and the read layout.

The C-only methods fermi_assemble, align_unitig_to_REFSEQ, verify_alns,
remap_RAs_w_unitigs and add_to_unitig_list are covered through
build_unitig_collection; filter_unitig_with_bad_aln is an empty stub
that nothing calls.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import pickle
import re

import pytest

from MACS3.IO.PeakIO import (PeakIO)
from MACS3.Signal.ReadAlignment import (ReadAlignment)
from MACS3.Signal.PosReadsInfo import (PosReadsInfo)
from MACS3.Signal.RACollection import (RACollection)
from MACS3.Signal.UnitigRACollection import (UnitigCollection)

BAM_CODES = "=ACMGRSVTWYHKDBN"
CIGAR_OPS = "MIDNSHP=X"
CHROM = b"chr1"
G0 = 5000


# ------------------------------------
# helpers
# ------------------------------------

def lcg_seq(n, seed=20261003):
    """Deterministic pseudo-random DNA from a linear congruential
    generator (the same on every platform and Python version)."""
    x = seed
    out = []
    for _ in range(n):
        x = (1103515245 * x + 12345) % 2**31
        out.append("ACGT"[(x >> 16) & 3])
    return "".join(out)


REF = lcg_seq(400)


def pack_seq(seq):
    codes = [BAM_CODES.index(c) for c in seq]
    if len(codes) % 2:
        codes.append(0)
    return bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))


def cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def md_tag(cigar, seq, start, ref=REF):
    md, run, r, q = "", 0, start, 0
    for n, op in cigar_ops(cigar):
        if op in "M=X":
            for _ in range(n):
                if seq[q] == ref[r]:
                    run += 1
                else:
                    md += "%d%s" % (run, ref[r])
                    run = 0
                q += 1
                r += 1
        elif op == "D":
            md += "%d^%s" % (run, ref[r:r + n])
            run = 0
            r += n
        elif op in "IS":
            q += n
        elif op == "N":
            r += n
    return md + str(run)


def make_ra(start, cigar, seq, qual=None, strand=0, name="r"):
    """ReadAlignment at G0 + start with the MD tag computed against REF."""
    if qual is None:
        qual = [30] * len(seq)
    span = sum(n for n, op in cigar_ops(cigar) if op in "MDN=X")
    return ReadAlignment(name.encode(), CHROM, G0 + start, G0 + start + span,
                         strand, pack_seq(seq), bytes(qual),
                         tuple((n << 4) | CIGAR_OPS.index(op)
                               for n, op in cigar_ops(cigar)),
                         md_tag(cigar, seq, start).encode())


def other_base(b):
    return "ACGT"[("ACGT".index(b) + 1) % 4]


def make_peak(start, end):
    """A peak as callvar gets it: a PeakContent from PeakIO."""
    peaks = PeakIO()
    peaks.add(CHROM, G0 + start, G0 + end)
    return peaks.get_data_from_chrom(CHROM)[0]


def pri_state(pri):
    """(n_reads_T, n_reads_C, bq_set_T, bq_set_C, n_strand, n_tips)."""
    s = pri.__getstate__()
    return s[6], s[7], s[4], s[5], s[9], s[10]


# four reads around the peak [G0+5, G0+40)
#   r1 [0, 20)   exact
#   r3 [5, 25)   10M2D8M, REF[15:17] deleted
#   r2 [10, 30)  reverse strand, mismatch at REF[15]
#   r4 [30, 45)  exact
ALT15 = other_base(REF[15])


def basic_reads():
    r1 = make_ra(0, "20M", REF[0:20], qual=[30 + i % 10 for i in range(20)],
                 name="r1")
    r2 = make_ra(10, "20M", REF[10:15] + ALT15 + REF[16:30], qual=[40] * 20,
                 strand=1, name="r2")
    r3 = make_ra(5, "10M2D8M", REF[5:15] + REF[17:25], name="r3")
    r4 = make_ra(30, "15M", REF[30:45], name="r4")
    return r1, r2, r3, r4


def basic_collection(control=()):
    r1, r2, r3, r4 = basic_reads()
    return RACollection(CHROM, make_peak(5, 40), [r4, r2, r1, r3], list(control))


def names(collection_fastq):
    """Read names of a FASTQ text (first line of every 4)."""
    return [ln[1:] for ln in bytes(collection_fastq).split(b"\n")[0::4] if ln]


# ------------------------------------
# construction and item access
# ------------------------------------

def test_construction_items():
    c = basic_collection()
    assert c["chrom"] == CHROM
    assert (c["left"], c["right"], c["length"]) == (G0 + 5, G0 + 40, 35)
    assert (c["RAs_left"], c["RAs_right"]) == (G0, G0 + 45)
    assert (c["count"], c["count_T"], c["count_C"]) == (4, 4, 0)


def test_peak_refseq_from_md_tags():
    # r1, r2 (mismatch corrected through MD) and r4 cover [0, 45)
    c = basic_collection()
    assert c["peak_refseq_ext"] == REF[0:45].encode()
    assert c["peak_refseq"] == REF[5:40].encode()


def test_peak_refseq_uncovered_positions_are_N():
    reads = [make_ra(0, "20M", REF[0:20]), make_ra(25, "20M", REF[25:45])]
    c = RACollection(CHROM, make_peak(-10, 50), reads)
    expected = "N" * 10 + REF[0:20] + "N" * 5 + REF[25:45] + "N" * 5
    assert c["peak_refseq_ext"] == expected.encode()
    assert c["peak_refseq"] == expected.encode()


def test_peak_refseq_uses_control_reads():
    reads_T = [make_ra(0, "20M", REF[0:20])]
    reads_C = [make_ra(20, "20M", REF[20:40])]
    c = RACollection(CHROM, make_peak(0, 40), reads_T, reads_C)
    assert c["peak_refseq"] == REF[0:40].encode()
    assert (c["count_T"], c["count_C"]) == (1, 1)


def test_single_read_collection():
    r = make_ra(3, "10M", REF[3:13])
    c = RACollection(CHROM, make_peak(3, 13), [r])
    assert (c["RAs_left"], c["RAs_right"]) == (G0 + 3, G0 + 13)
    assert c["peak_refseq"] == REF[3:13].encode()


def test_empty_treatment_raises():
    with pytest.raises(Exception,
                       match="No reads from ChIP sample to construct RAcollection!"):
        RACollection(CHROM, make_peak(0, 10), [], [make_ra(0, "10M", REF[0:10])])


def test_unknown_key_raises():
    with pytest.raises(KeyError, match="Unavailable key"):
        basic_collection()["peak"]


def test_peak_can_be_a_dict():
    r = make_ra(0, "10M", REF[0:10])
    c = RACollection(CHROM, {"start": G0, "end": G0 + 10}, [r])
    assert c["peak_refseq"] == REF[0:10].encode()


def test_pickle_round_trip():
    c = basic_collection(control=[make_ra(12, "10M", REF[12:22], name="c1")])
    c2 = pickle.loads(pickle.dumps(c))
    for key in ("chrom", "left", "right", "RAs_left", "RAs_right", "length",
                "count_T", "count_C", "peak_refseq", "peak_refseq_ext"):
        assert c2[key] == c[key]
    assert c2.get_FASTQ() == c.get_FASTQ()


# ------------------------------------
# sort / get_FASTQ
# ------------------------------------

def test_sort_by_lpos_and_get_FASTQ_order():
    c = basic_collection(control=[make_ra(12, "10M", REF[12:22], name="c2"),
                                  make_ra(2, "10M", REF[2:12], name="c1")])
    assert names(c.get_FASTQ()) == [b"r1", b"r3", b"r2", b"r4", b"c1", b"c2"]
    c.sort()
    assert names(c.get_FASTQ()) == [b"r1", b"r3", b"r2", b"r4", b"c1", b"c2"]


def test_sort_is_stable_for_equal_lpos():
    reads = [make_ra(0, "10M", REF[0:10], name="b"),
             make_ra(0, "12M", REF[0:12], name="a"),
             make_ra(0, "8M", REF[0:8], name="c")]
    c = RACollection(CHROM, make_peak(0, 12), reads)
    assert names(c.get_FASTQ()) == [b"b", b"a", b"c"]


def test_get_FASTQ_is_concatenation_of_read_FASTQ():
    r1, r2, r3, r4 = basic_reads()
    c = basic_collection()
    fastq = c.get_FASTQ()
    assert isinstance(fastq, bytearray)
    assert fastq == (r1.get_FASTQ() + r3.get_FASTQ() + r2.get_FASTQ()
                     + r4.get_FASTQ())


# ------------------------------------
# remove_outliers / n_edits_sum
# ------------------------------------

def edited_read(start, n_mismatch, name="bad"):
    seq = list(REF[start:start + 10])
    for i in range(n_mismatch):
        seq[2 + 2 * i] = other_base(seq[2 + 2 * i])
    return make_ra(start, "10M", "".join(seq), name=name)


def clean_reads(n, prefix="ok"):
    return [make_ra(i, "10M", REF[i:i + 10], name="%s%d" % (prefix, i))
            for i in range(n)]


def test_remove_outliers_drops_top_5_percent():
    # 21 reads: the 3-edit read is above the 95% quantile (index 19 -> 0)
    c = RACollection(CHROM, make_peak(0, 40), clean_reads(20) + [edited_read(5, 3)])
    c.remove_outliers()
    assert c["count_T"] == 20
    assert b"bad" not in names(c.get_FASTQ())


def test_remove_outliers_includes_control_reads():
    c = RACollection(CHROM, make_peak(0, 40), clean_reads(10),
                     clean_reads(10, "c") + [edited_read(5, 3)])
    c.remove_outliers(percent=5)
    assert (c["count_T"], c["count_C"]) == (10, 10)


def test_remove_outliers_keeps_ties_at_threshold():
    # n_edits [0]*10 + [1]*10 + [2]: threshold index 19 -> 1; only 2 goes
    reads = (clean_reads(10) + [edited_read(10 + i, 1, "one%d" % i) for i in range(10)]
             + [edited_read(5, 2)])
    c = RACollection(CHROM, make_peak(0, 40), reads)
    c.remove_outliers()
    assert c["count_T"] == 20


# ------------------------------------
# get_PosReadsInfo_ref_pos
# ------------------------------------

def test_get_PosReadsInfo_ref_pos_mixed_alleles():
    c = basic_collection()
    pri = c.get_PosReadsInfo_ref_pos(G0 + 15, REF[15].encode())
    assert isinstance(pri, PosReadsInfo)
    nT, nC, bqT, bqC, strand, tips = pri_state(pri)
    ref, alt = REF[15].encode(), ALT15.encode()
    # r1: reference base, quality 30 + 15 % 10; r2: SNV on the reverse
    # strand, quality 40; r3: deleted base, '*' with quality 93; r4 misses
    assert (nT[ref], nT[alt], nT[b"*"]) == (1, 1, 1)
    assert sum(nT.values()) == 3 and sum(nC.values()) == 0
    assert (bqT[ref], bqT[alt], bqT[b"*"]) == ([35], [40], [93])
    assert strand[1][alt] == 1 and strand[0][ref] == 1 and strand[0][b"*"] == 1
    assert sum(tips.values()) == 0


def test_get_PosReadsInfo_ref_pos_quality_cutoff():
    c = basic_collection()
    pri = c.get_PosReadsInfo_ref_pos(G0 + 15, REF[15].encode(), Q=35)
    nT = pri_state(pri)[0]
    assert nT[REF[15].encode()] == 0       # quality 35 is not > 35
    assert nT[ALT15.encode()] == 1 and nT[b"*"] == 1


def test_get_PosReadsInfo_ref_pos_tips_and_control():
    c = basic_collection(control=[make_ra(12, "10M", REF[12:22], name="c1")])
    # G0 + 29 is only covered by r2, as its last base (a tip)
    pri = c.get_PosReadsInfo_ref_pos(G0 + 29, REF[29].encode())
    nT, nC, bqT, bqC, strand, tips = pri_state(pri)
    assert nT[REF[29].encode()] == 1 and tips[REF[29].encode()] == 1
    pri = c.get_PosReadsInfo_ref_pos(G0 + 13, REF[13].encode())
    nT, nC, bqT, bqC, strand, tips = pri_state(pri)
    ref = REF[13].encode()
    assert (nT[ref], nC[ref]) == (3, 1)
    assert bqC[ref] == [30]


def test_get_PosReadsInfo_ref_pos_uncovered():
    c = basic_collection()
    pri = c.get_PosReadsInfo_ref_pos(G0 + 100, REF[100].encode())
    assert pri.raw_read_depth() == 0


# ------------------------------------
# build_unitig_collection (fermi-lite)
# ------------------------------------

def deletion_site(lo):
    """First position >= lo whose base differs from both neighbours, so a
    1-bp deletion there has a single alignment."""
    d = lo
    while REF[d - 1] == REF[d] or REF[d] == REF[d + 1]:
        d += 1
    return d


def deletion_reads(d, step=4, readlen=50, span=300, het=False):
    """Reads tiling [0, span) of the haplotype REF without REF[d]; with
    het=True every other read comes from REF instead."""
    alt = REF[:d] + REF[d + 1:]
    reads = []
    for j, a in enumerate(range(0, span - readlen + 1, step)):
        qual = [30 + (a + i) % 11 for i in range(readlen)]
        name = "r%d" % a
        strand = j % 2
        if het and j % 2:
            reads.append(make_ra(a, "%dM" % readlen, REF[a:a + readlen],
                                 qual=qual, strand=strand, name=name))
        elif a + readlen <= d:
            reads.append(make_ra(a, "%dM" % readlen, alt[a:a + readlen],
                                 qual=qual, strand=strand, name=name))
        elif a >= d:
            reads.append(make_ra(a + 1, "%dM" % readlen, alt[a:a + readlen],
                                 qual=qual, strand=strand, name=name))
        else:
            k = d - a
            reads.append(make_ra(a, "%dM1D%dM" % (k, readlen - k),
                                 alt[a:a + readlen], qual=qual, strand=strand,
                                 name=name))
    return alt, reads


def test_build_unitig_collection_reference_only():
    # 63 error-free reads tile REF[0:298): one unitig equal to that span
    reads = [make_ra(a, "50M", REF[a:a + 50], name="r%d" % a)
             for a in range(0, 251, 4)]
    c = RACollection(CHROM, make_peak(50, 250), reads)
    uc = c.build_unitig_collection(30)
    assert isinstance(uc, UnitigCollection)
    assert uc["count"] == 1
    ura = uc["URAs_list"][0]
    assert (ura["lpos"], ura["rpos"], ura["count"]) == (G0, G0 + 298, 63)
    assert ura["unitig_aln"] == ura["reference_aln"] == REF[0:298].encode()
    pri = uc.get_PosReadsInfo_ref_pos(G0 + 150, REF[150].encode())
    nT = pri.__getstate__()[6]
    # positions G0+150: reads starting at 104..148 step 4 -> 12 reads
    assert nT[REF[150].encode()] == 12
    assert pri.raw_read_depth() == 12


def test_build_unitig_collection_homozygous_deletion():
    d = deletion_site(150)
    alt, reads = deletion_reads(d)
    c = RACollection(CHROM, make_peak(50, 250), reads)
    uc = c.build_unitig_collection(30)
    assert isinstance(uc, UnitigCollection)
    # one unitig: the deletion haplotype over the 298 bases the reads
    # cover, aligned to REF[0:299] with a gap at the deleted base
    assert uc["count"] == 1
    ura = uc["URAs_list"][0]
    assert (ura["lpos"], ura["rpos"], ura["count"]) == (G0, G0 + 299, len(reads))
    assert ura["unitig_aln"] == (alt[:d] + "-" + alt[d:298]).encode()
    assert ura["reference_aln"] == REF[0:299].encode()
    assert ura.get_variant_bq_by_ref_pos(G0 + d)[0] == bytearray(b"*")
    pri = uc.get_PosReadsInfo_ref_pos(G0 + d, REF[d].encode())
    nT = pri.__getstate__()[6]
    assert nT[b"*"] > 0 and nT[REF[d].encode()] == 0


def test_build_unitig_collection_heterozygous_deletion():
    d = deletion_site(150)
    alt, reads = deletion_reads(d, step=3, het=True)
    c = RACollection(CHROM, make_peak(50, 250), reads)
    uc = c.build_unitig_collection(30)
    assert isinstance(uc, UnitigCollection)
    for ura in uc["URAs_list"]:
        assert ura["seq"] in alt.encode() or ura["seq"] in REF.encode()
    # both haplotypes are assembled over the deletion site
    alleles = {bytes(u.get_variant_bq_by_ref_pos(G0 + d)[0])
               for u in uc["URAs_list"] if u["lpos"] <= G0 + d < u["rpos"]}
    assert alleles == {b"*", REF[d].encode()}


def test_build_unitig_collection_no_unitig_returns_0():
    # minimum overlap equal to the read length: no two reads can be joined
    reads = [make_ra(a, "30M", REF[a:a + 30], name="r%d" % a)
             for a in range(0, 31, 10)]
    c = RACollection(CHROM, make_peak(0, 60), reads)
    assert c.build_unitig_collection(30) == 0
