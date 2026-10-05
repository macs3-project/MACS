#!/usr/bin/env python

"""Module Description: Test the Variant and PeakVariants classes in
MACS3.Signal.PeakVariants, which hold the variants called in one peak
and write them as VCF records.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import pickle

import numpy as np
import pytest

from MACS3.Signal.PeakVariants import (Variant,
                                       PeakVariants)

FIELDS = ("ref_allele", "alt_allele", "GQ", "filter", "type",
          "mutation_type", "top1allele", "top2allele", "DPT", "DPC", "DP1T",
          "DP2T", "DP1C", "DP2C", "PLUS1T", "PLUS2T", "MINUS1T", "MINUS2T",
          "deltaBIC", "BIC_homo_major", "BIC_homo_minor", "BIC_heter_noAS",
          "BIC_heter_AS", "AR", "GT", "DP", "PL_00", "PL_01", "PL_11")

# float values are exact in C float, so "%.2f" is unambiguous
DEFAULTS = dict(ref_allele="A", alt_allele="G", GQ=50, filter=".",
                type="heter_noAS", mutation_type="SNV", top1allele="A",
                top2allele="G", DPT=10, DPC=2, DP1T=6, DP2T=4, DP1C=1, DP2C=1,
                PLUS1T=3, PLUS2T=2, MINUS1T=3, MINUS2T=2, deltaBIC=12.5,
                BIC_homo_major=40.25, BIC_homo_minor=60.75,
                BIC_heter_noAS=27.75, BIC_heter_AS=30.0, AR=0.75, GT="0/1",
                DP=12, PL_00=100, PL_01=0, PL_11=50)

DEFAULT_VCF = ("A\tG\t50\t.\tM=heter_noAS;MT=SNV;DPT=10;DPC=2;DP1T=6A;DP2T=4G;"
               "DP1C=1A;DP2C=1G;SB=3,2,3,2;DBIC=12.50;BICHOMOMAJOR=40.25;"
               "BICHOMOMINOR=60.75;BICHETERNOAS=27.75;BICHETERAS=30.00;"
               "AR=0.75\tGT:DP:GQ:PL\t0/1:12:50:100,0,50")

# reference sequence of the peak [100, 116)
START = 100
REFSEQ = b"GATTACAGCTTGCAAC"


def make_variant(**kw):
    d = dict(DEFAULTS)
    d.update(kw)
    return Variant(*[d[f] for f in FIELDS])


def vstate(v):
    return dict(zip(FIELDS, v.__getstate__()))


def snv(ref, alt):
    return make_variant(ref_allele=ref, alt_allele=alt, top1allele=ref,
                        top2allele=alt)


def deletion(ref, top1_is_ref=True):
    """Deletion of the reference base ``ref``: heterozygous when top1 is
    the reference allele, homozygous when the deletion '*' is top1."""
    if top1_is_ref:
        return make_variant(ref_allele=ref, alt_allele="*",
                            mutation_type="Deletion", top1allele=ref,
                            top2allele="*")
    return make_variant(ref_allele=ref, alt_allele="*", type="homo", GT="1/1",
                        mutation_type="Deletion", top1allele="*",
                        top2allele=ref)


def insertion(ref, alt):
    return make_variant(ref_allele=ref, alt_allele=alt,
                        mutation_type="Insertion", top1allele=ref,
                        top2allele=alt)


def peak(**variants):
    pv = PeakVariants("chr1", START, START + len(REFSEQ), REFSEQ)
    for p, v in variants.items():
        pv.add_variant(int(p[1:]), v)
    return pv


# ------------------------------------
# Variant
# ------------------------------------

def test_Variant_state_keeps_every_field():
    v = make_variant()
    assert vstate(v) == DEFAULTS


def test_Variant_toVCF():
    assert make_variant().toVCF() == DEFAULT_VCF


def test_Variant_toVCF_rounds_through_C_float():
    # deltaBIC, BIC values and AR are stored as C floats
    v = make_variant(deltaBIC=23.215, AR=0.715, BIC_heter_AS=1e-3)
    info = dict(x.split("=") for x in v.toVCF().split("\t")[4].split(";"))
    assert info["DBIC"] == "%.2f" % float(np.float32(23.215))
    assert info["AR"] == "%.2f" % float(np.float32(0.715))
    assert info["BICHETERAS"] == "0.00"


def test_Variant_float_GQ_is_truncated():
    v = make_variant(GQ=58.9, PL_00=159.99, PL_11=58.9)
    fields = v.toVCF().split("\t")
    assert fields[2] == "58"
    assert fields[6] == "0/1:12:58:159,0,58"


def test_Variant_toVCF_multi_allelic():
    v = make_variant(ref_allele="A", alt_allele="G,T", top1allele="G",
                     top2allele="T", mutation_type="SNV,SNV", GT="1/2")
    fields = v.toVCF().split("\t")
    assert fields[:2] == ["A", "G,T"]
    assert "MT=SNV,SNV;" in fields[4] and "DP1T=6G;DP2T=4T;" in fields[4]
    assert fields[6].startswith("1/2:")


@pytest.mark.parametrize("value", [0, 2**31 - 1, -(2**31)])
def test_Variant_int32_depths(value):
    v = make_variant(DPT=value)
    assert vstate(v)["DPT"] == value
    assert ";DPT=%d;" % value in v.toVCF()


def test_Variant_depth_overflow():
    with pytest.raises(OverflowError):
        make_variant(DPT=2**31)


@pytest.mark.parametrize("mt,indel,only_del,only_ins", [
    ("SNV", False, False, False),
    ("Deletion", True, True, False),
    ("Insertion", True, False, True),
    ("SNV,Deletion", True, False, False),
    ("Insertion,SNV", True, False, False),
    ("SNV,SNV", False, False, False),
    ("Deletion,Insertion", True, False, False),
])
def test_Variant_indel_flags(mt, indel, only_del, only_ins):
    v = make_variant(mutation_type=mt)
    assert v.is_indel() is indel
    assert v.is_only_del() is only_del
    assert v.is_only_insertion() is only_ins


@pytest.mark.parametrize("AR,top1,ar,expected", [
    (0.85, "A", None, True),       # default cutoff 0.85, >= comparison
    (0.849, "A", None, False),
    (0.99, "G", None, False),      # top1 is not the reference allele
    (0.85, "A", 0.9, False),
    (0.95, "A", 0.9, True),
    (1.0, "A", 1.0, True),
])
def test_Variant_is_refer_biased_01(AR, top1, ar, expected):
    v = make_variant(AR=AR, top1allele=top1)
    result = v.is_refer_biased_01() if ar is None else v.is_refer_biased_01(ar)
    assert result is expected


@pytest.mark.parametrize("top1,top2,is1,is2", [
    ("A", "G", True, False),
    ("G", "A", False, True),
    ("G", "T", False, False),
    ("A", "A", True, True),
])
def test_Variant_top_is_reference(top1, top2, is1, is2):
    v = make_variant(ref_allele="A", top1allele=top1, top2allele=top2)
    assert v.top1isreference() is is1
    assert v.top2isreference() is is2


@pytest.mark.parametrize("key", ["ref_allele", "alt_allele", "top1allele",
                                 "top2allele"])
def test_Variant_getitem_setitem_alleles(key):
    v = make_variant()
    assert v[key] == DEFAULTS[key]
    v[key] = "ACGT"
    assert v[key] == "ACGT"
    assert vstate(v)[key] == "ACGT"


@pytest.mark.parametrize("key", ["GQ", "chrom", ""])
def test_Variant_getitem_setitem_unknown_key(key):
    v = make_variant()
    with pytest.raises(Exception, match="keyname is not accessible"):
        v[key]
    with pytest.raises(Exception, match="keyname is not accessible"):
        v[key] = 1


def test_Variant_rejects_wrong_types():
    with pytest.raises(TypeError):
        make_variant(ref_allele=b"A")
    with pytest.raises(TypeError):
        make_variant(DPT="10")


def test_Variant_pickle_round_trip():
    v = make_variant(deltaBIC=23.215)
    w = pickle.loads(pickle.dumps(v))
    assert w.toVCF() == v.toVCF()
    assert vstate(w) == vstate(v)


# ------------------------------------
# PeakVariants: container methods
# ------------------------------------

def test_PeakVariants_empty():
    pv = peak()
    assert pv.n_variants() == 0
    assert pv.has_indel() is False
    assert pv.has_refer_biased_01() is False
    assert pv.get_refer_biased_01s() == []
    assert pv.toVCF() == ""
    pv.fix_indels()
    assert pv.n_variants() == 0


def test_PeakVariants_add_variant_and_overwrite():
    pv = peak()
    pv.add_variant(105, snv("C", "T"))
    pv.add_variant(102, snv("T", "A"))
    assert pv.n_variants() == 2
    pv.add_variant(105, snv("C", "G"))
    assert pv.n_variants() == 2
    assert pv.toVCF().splitlines()[1].split("\t")[:5] == ["chr1", "106", ".", "C", "G"]


def test_PeakVariants_add_variant_requires_Variant():
    with pytest.raises(TypeError):
        peak().add_variant(105, "not a variant")


def test_PeakVariants_toVCF_sorted_and_one_based():
    pv = peak(p110=snv("T", "C"), p101=snv("A", "G"), p2147483646=snv("C", "A"))
    lines = pv.toVCF().split("\n")
    assert lines[-1] == ""
    assert [ln.split("\t")[1] for ln in lines[:-1]] == ["102", "111", "2147483647"]
    assert lines[0] == "chr1\t102\t.\t" + snv("A", "G").toVCF()


def test_PeakVariants_has_indel():
    assert peak(p101=snv("A", "G"), p103=snv("T", "C")).has_indel() is False
    assert peak(p101=snv("A", "G"), p108=deletion("C")).has_indel() is True
    assert peak(p107=insertion("G", "GA")).has_indel() is True


def test_PeakVariants_refer_biased():
    pv = peak(p110=make_variant(AR=0.9), p101=make_variant(AR=0.86),
              p105=make_variant(AR=0.5), p106=make_variant(AR=0.95,
                                                          top1allele="G"))
    assert pv.has_refer_biased_01() is True
    assert pv.get_refer_biased_01s() == [101, 110]
    assert peak(p105=make_variant(AR=0.5)).has_refer_biased_01() is False


def test_PeakVariants_remove_variant():
    pv = peak(p101=snv("A", "G"), p103=snv("T", "C"))
    pv.remove_variant(101)
    assert pv.n_variants() == 1
    assert pv.toVCF().split("\t")[1] == "104"
    with pytest.raises(AssertionError):
        pv.remove_variant(101)


def test_PeakVariants_replace_variant():
    pv = peak(p101=snv("A", "G"))
    pv.replace_variant(101, snv("A", "T"))
    assert pv.n_variants() == 1
    assert pv.toVCF().split("\t")[4] == "T"
    with pytest.raises(AssertionError):
        pv.replace_variant(102, snv("T", "A"))


def test_PeakVariants_pickle_keeps_variants_and_chrom():
    pv = peak(p101=snv("A", "G"), p110=deletion("T"))
    pv2 = pickle.loads(pickle.dumps(pv))
    assert pv2.n_variants() == 2
    assert pv2.toVCF() == pv.toVCF()


# ------------------------------------
# PeakVariants.fix_indels
# ------------------------------------

def test_fix_indels_leaves_snvs_and_insertions():
    pv = peak(p101=snv("A", "G"), p102=snv("T", "C"), p107=insertion("G", "GA"))
    before = pv.toVCF()
    pv.fix_indels()
    assert pv.toVCF() == before


def test_fix_indels_removes_deletion_after_insertion():
    pv = peak(p101=snv("A", "G"), p107=insertion("G", "GA"), p108=deletion("C"))
    pv.fix_indels()
    assert pv.n_variants() == 1
    assert pv.toVCF().split("\t")[1] == "102"


def test_fix_indels_deletion_at_peak_start_is_kept():
    # p == start: there is no preceding base inside the peak
    pv = peak(p100=deletion("G"))
    before = pv.toVCF()
    pv.fix_indels()
    assert pv.toVCF() == before


def test_fix_indels_deletion_after_snv_is_not_anchored():
    pv = peak(p107=snv("G", "A"), p108=deletion("C"))
    before = pv.toVCF()
    pv.fix_indels()
    assert pv.toVCF() == before


@pytest.mark.parametrize("n_del", [2, 3])
def test_fix_indels_merges_consecutive_heterozygous_deletions(n_del):
    # deletions at 108..: C, T, T; an SNV at 107 keeps them unanchored
    variants = {"p107": snv("G", "A")}
    for i in range(n_del):
        variants["p%d" % (108 + i)] = deletion(chr(REFSEQ[8 + i]))
    pv = peak(**variants)
    pv.fix_indels()
    merged = REFSEQ[8:8 + n_del].decode()
    assert pv.n_variants() == 2
    line = pv.toVCF().splitlines()[1].split("\t")
    assert line[1:5] == ["109", ".", merged, "*"]
    assert ";DP1T=6%s;DP2T=4*;" % merged in line[7]


def test_fix_indels_anchoring_moves_deletion_one_base_left():
    """The anchored record replaces the deletion at p by one at p - 1
    carrying the other fields of the original variant."""
    pv = peak(p108=deletion("C"))
    pv.fix_indels()
    assert pv.n_variants() == 1
    fields = pv.toVCF().rstrip("\n").split("\t")
    assert fields[1] == "108"
    assert fields[5:7] == ["50", "."]
    assert fields[8:] == ["GT:DP:GQ:PL", "0/1:12:50:100,0,50"]
    assert fields[3].endswith("C") and len(fields[3]) == 2
    assert fields[4] == fields[3][0]
