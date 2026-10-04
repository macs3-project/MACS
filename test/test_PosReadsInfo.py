#!/usr/bin/env python

"""Module Description: Test the PosReadsInfo class in
MACS3.Signal.PosReadsInfo, which collects the alleles and base
qualities seen at one reference position and calls the genotype.

Genotype likelihoods are checked against direct numpy/scipy
implementations of the four models (homozygous major, homozygous minor,
heterozygous without and with allele-specific binding); see
test_VariantStat.py for the formulas. The C-only methods SB_score_ChIP
and SB_score_ATAC are not called by any Python-visible method.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import math
import pickle

import numpy as np
import pytest
from scipy.stats import binom

from MACS3.Signal.PosReadsInfo import (PosReadsInfo)
from MACS3.Signal.PeakVariants import (Variant)

LN10 = math.log(10)

STATE_FIELDS = ("ref_pos", "ref_allele", "alt_allele", "filterout",
                "bq_set_T", "bq_set_C", "n_reads_T", "n_reads_C", "n_reads",
                "n_strand", "n_tips", "top1allele", "top2allele",
                "top12alleles_ratio", "lnL_homo_major", "lnL_heter_AS",
                "lnL_heter_noAS", "lnL_homo_minor", "BIC_homo_major",
                "BIC_heter_AS", "BIC_heter_noAS", "BIC_homo_minor",
                "heter_noAS_kc", "heter_noAS_ki", "heter_AS_kc",
                "heter_AS_ki", "heter_AS_alleleratio", "GQ_homo_major",
                "GQ_heter_noAS", "GQ_heter_AS", "GQ_heter_ASsig", "GQ", "GT",
                "type", "hasfermiinfor", "fermiNTs")


# ------------------------------------
# helpers
# ------------------------------------

def st(pri):
    """PosReadsInfo state as a dict keyed by attribute name."""
    state = pri.__getstate__()
    assert len(state) == len(STATE_FIELDS)
    return dict(zip(STATE_FIELDS, state))


def set_state(pri, **changes):
    s = st(pri)
    s.update(changes)
    pri.__setstate__(tuple(s[k] for k in STATE_FIELDS))


def make_pri(ref=b"A", T=(), C=(), pos=100, Q=20):
    """PosReadsInfo with treatment reads T = [(allele, bq, strand, tip)]
    and control reads C = [(allele, bq, strand)]."""
    pri = PosReadsInfo(pos, ref)
    for i, (allele, q, strand, tip) in enumerate(T):
        pri.add_T(i, allele, q, strand, tip, Q=Q)
    for i, (allele, q, strand) in enumerate(C):
        pri.add_C(i, allele, q, strand, Q=Q)
    return pri


def reads(allele, n, q=30, strand=0, tip=False):
    return [(allele, q, strand, tip)] * n


def err_rate(quals):
    return 10.0 ** (-np.asarray(quals, dtype=float) / 10.0)


def ref_homo(t1T, t1C, t2T, t2C):
    """Every top1 base right, every top2 base an error; BIC = -2 lnL."""
    e1 = err_rate(list(t1T) + list(t1C))
    e2 = err_rate(list(t2T) + list(t2C))
    lnL = float(np.log1p(-e1).sum() + np.log(e2).sum())
    return lnL, -2 * lnL


def ref_lnL_k(me, ne, k, p):
    tn = len(me) + len(ne)
    f = k / tn
    em, en = err_rate(me), err_rate(ne)
    return float(binom.logpmf(k, tn, p)
                 + np.log((1 - em) * f + em * (1 - f)).sum()
                 + np.log((1 - en) * (1 - f) + en * f).sum())


def ref_noAS_max(me, ne):
    """Maximum over k with p = 0.5. When one allele has no read, MACS3
    fixes k at the observed count instead of searching (by design)."""
    if len(me) == 0 or len(ne) == 0:
        return ref_lnL_k(me, ne, len(me), 0.5)
    return max(ref_lnL_k(me, ne, k, 0.5) for k in range(len(me) + len(ne) + 1))


def ref_AS_max(me, ne, max_ar):
    """Maximum over k with p = k/tn capped to [1-maxAR, maxAR]; k = tn
    when every read carries top1 (no search, as in MACS3)."""
    tn = len(me) + len(ne)
    mar = float(np.float32(max_ar))
    if len(ne) == 0:
        return ref_lnL_k(me, ne, tn, mar)
    return max(ref_lnL_k(me, ne, k, min(max(k / tn, 1 - mar), mar))
               for k in range(tn + 1))


def ref_models(t1T, t1C, t2T, t2C, max_ar=0.99):
    """(lnL, BIC) of the four genotype models, independent of MACS3."""
    tn_T, tn_C = len(t1T) + len(t2T), len(t1C) + len(t2C)
    homo_major = ref_homo(t1T, t1C, t2T, t2C)
    homo_minor = ref_homo(t2T, t2C, t1T, t1C)
    lnL = ref_noAS_max(t1T, t2T)
    pen = math.log(tn_T)
    lnL_as = ref_AS_max(t1T, t2T, max_ar)
    pen_as = 2 * math.log(tn_T)
    if tn_C:
        lnL_c = ref_noAS_max(t1C, t2C)
        lnL += lnL_c
        lnL_as += lnL_c
        pen += math.log(tn_C)
        pen_as += math.log(tn_C)
    return {"homo_major": homo_major, "homo_minor": homo_minor,
            "heter_noAS": (lnL, -2 * lnL + pen),
            "heter_AS": (lnL_as, -2 * lnL_as + pen_as)}


def phred(lnL):
    return -10.0 * lnL / LN10


def assert_models(s, m):
    for name in ("homo_major", "homo_minor", "heter_noAS", "heter_AS"):
        assert s["lnL_" + name] == pytest.approx(m[name][0], rel=1e-10, abs=1e-12)
        assert s["BIC_" + name] == pytest.approx(m[name][1], rel=1e-10, abs=1e-12)


# ------------------------------------
# construction, filterflag, raw_read_depth
# ------------------------------------

def test_new_PosReadsInfo_defaults():
    pri = PosReadsInfo(12345, b"A")
    s = st(pri)
    assert s["ref_pos"] == 12345
    assert s["ref_allele"] == b"A"
    assert s["alt_allele"] == b"."
    assert s["GT"] == "unsure"
    assert s["GQ"] == 0
    assert pri.filterflag() is False
    assert list(s["n_reads"]) == [b"A", b"C", b"G", b"T", b"N", b"*"]
    assert all(v == 0 for v in s["n_reads"].values())
    assert all(v == [] for v in s["bq_set_T"].values())


def test_new_PosReadsInfo_reference_key_comes_first():
    s = st(PosReadsInfo(0, b"G"))
    assert list(s["n_reads"]) == [b"G", b"A", b"C", b"T", b"N", b"*"]
    s = st(PosReadsInfo(0, b"N"))
    assert list(s["n_reads_T"]) == [b"N", b"A", b"C", b"G", b"T", b"*"]


@pytest.mark.parametrize("pos", [0, 2**31 - 1, 2**40])
def test_PosReadsInfo_position_range(pos):
    assert st(PosReadsInfo(pos, b"C"))["ref_pos"] == pos


@pytest.mark.parametrize("opt,expected", [("all", 5), ("T", 3), ("C", 2)])
def test_raw_read_depth(opt, expected):
    pri = make_pri(T=reads(b"A", 2) + reads(b"G", 1),
                   C=[(b"A", 30, 0), (b"T", 30, 1)])
    assert pri.raw_read_depth(opt=opt) == expected


def test_raw_read_depth_default_and_empty():
    pri = PosReadsInfo(1, b"A")
    assert pri.raw_read_depth() == 0
    pri.add_T(0, b"A", 30, 0, False)
    assert pri.raw_read_depth() == 1


def test_raw_read_depth_bad_option():
    with pytest.raises(Exception, match="opt should be either 'all', 'T' or 'C'."):
        PosReadsInfo(1, b"A").raw_read_depth(opt="X")


# ------------------------------------
# add_T / add_C
# ------------------------------------

def test_add_T_counts_and_strands():
    pri = make_pri(T=[(b"A", 30, 0, False), (b"A", 35, 1, True),
                      (b"G", 40, 1, False)])
    s = st(pri)
    assert s["bq_set_T"][b"A"] == [30, 35]
    assert s["bq_set_T"][b"G"] == [40]
    assert s["n_reads_T"][b"A"] == 2 and s["n_reads"][b"A"] == 2
    assert s["n_strand"][0][b"A"] == 1 and s["n_strand"][1][b"A"] == 1
    assert s["n_strand"][1][b"G"] == 1
    assert s["n_tips"][b"A"] == 1 and s["n_tips"][b"G"] == 0
    assert s["n_reads_C"][b"A"] == 0


@pytest.mark.parametrize("q,Q,counted", [
    (20, 20, False),     # bq must be strictly larger than Q
    (21, 20, True),
    (0, 0, False),
    (1, 0, True),
    (93, 92, True),
    (93, 93, False),
])
def test_add_T_quality_cutoff(q, Q, counted):
    pri = PosReadsInfo(1, b"A")
    pri.add_T(0, b"C", q, 0, False, Q=Q)
    assert pri.raw_read_depth("T") == int(counted)


def test_add_T_default_Q_is_20():
    pri = PosReadsInfo(1, b"A")
    pri.add_T(0, b"C", 20, 0, False)
    pri.add_T(1, b"C", 21, 0, False)
    assert st(pri)["bq_set_T"][b"C"] == [21]


def test_add_T_new_allele_creates_every_key():
    pri = PosReadsInfo(1, b"A")
    pri.add_T(0, b"ATT", 30, 1, True)
    s = st(pri)
    assert s["bq_set_T"][b"ATT"] == [30]
    assert s["bq_set_C"][b"ATT"] == []
    assert s["n_reads_T"][b"ATT"] == 1 and s["n_reads_C"][b"ATT"] == 0
    assert s["n_strand"][0][b"ATT"] == 0 and s["n_strand"][1][b"ATT"] == 1
    assert s["n_tips"][b"ATT"] == 1


def test_add_T_bad_strand():
    with pytest.raises(IndexError):
        PosReadsInfo(1, b"A").add_T(0, b"A", 30, 2, False)


def test_add_C_counts_without_strand_or_tip():
    pri = make_pri(C=[(b"A", 30, 0), (b"A", 31, 1), (b"*", 93, 1)])
    s = st(pri)
    assert s["bq_set_C"][b"A"] == [30, 31]
    assert s["bq_set_C"][b"*"] == [93]
    assert s["n_reads_C"][b"A"] == 2 and s["n_reads"][b"A"] == 2
    assert s["n_reads_T"][b"A"] == 0
    assert s["n_strand"][0][b"A"] == 0 and s["n_strand"][1][b"A"] == 0


def test_add_C_quality_cutoff_and_new_allele():
    pri = PosReadsInfo(1, b"A")
    pri.add_C(0, b"AG", 25, 0, Q=25)
    assert pri.raw_read_depth("C") == 0
    pri.add_C(0, b"AG", 26, 0, Q=25)
    s = st(pri)
    assert s["n_reads_C"][b"AG"] == 1
    assert s["bq_set_T"][b"AG"] == [] and s["n_tips"][b"AG"] == 0


# ------------------------------------
# update_top_alleles
# ------------------------------------

def test_update_top_alleles_heterozygous():
    pri = make_pri(T=reads(b"A", 5) + reads(b"G", 4))
    pri.update_top_alleles()
    s = st(pri)
    assert (s["top1allele"], s["top2allele"]) == (b"A", b"G")
    assert pri.filterflag() is False
    assert s["top12alleles_ratio"] == 1.0


def test_update_top_alleles_tie_follows_key_order():
    # 3 A and 3 T with reference G: sorted() is stable, keys are G,A,C,T,...
    pri = make_pri(ref=b"G", T=reads(b"T", 3) + reads(b"A", 3))
    pri.update_top_alleles()
    s = st(pri)
    assert (s["top1allele"], s["top2allele"]) == (b"A", b"T")


def test_update_top_alleles_reference_only_is_homo_ref():
    pri = make_pri(T=reads(b"A", 5))
    pri.update_top_alleles()
    s = st(pri)
    assert s["top1allele"] == b"A"
    assert s["type"] == "homo_ref"
    assert pri.filterflag() is True


@pytest.mark.parametrize("n_alt,min_count,filtered", [
    (1, 2, True),     # 1 alt read < 2: alt allele dropped -> homo_ref
    (2, 2, False),
    (2, 3, True),
    (3, 3, False),
])
def test_update_top_alleles_altallele_count(n_alt, min_count, filtered):
    pri = make_pri(T=reads(b"A", 8) + reads(b"G", n_alt))
    pri.update_top_alleles(min_altallele_count=min_count)
    s = st(pri)
    assert pri.filterflag() is filtered
    if filtered:
        assert s["type"] == "homo_ref"
        assert s["n_reads_T"][b"G"] == 0 and s["bq_set_T"][b"G"] == []


def test_update_top_alleles_tip_reads_do_not_count():
    # 3 alt reads but all at read ends: 3 - 3 tips < 2
    pri = make_pri(T=reads(b"A", 5) + reads(b"G", 3, tip=True))
    pri.update_top_alleles()
    assert pri.filterflag() is True
    assert st(pri)["n_tips"][b"G"] == 0


@pytest.mark.parametrize("max_ar,filtered", [
    (0.95, True),     # 39 / 41 = 0.951 > 0.95: minor allele dropped
    (0.96, False),
])
def test_update_top_alleles_max_allowed_ar(max_ar, filtered):
    pri = make_pri(T=reads(b"A", 39) + reads(b"G", 2))
    pri.update_top_alleles(max_allowed_ar=max_ar)
    assert pri.filterflag() is filtered


@pytest.mark.parametrize("min_ratio,filtered", [(0.8, True), (0.7, False)])
def test_update_top_alleles_top12_ratio(min_ratio, filtered):
    # top1 A (4) + top2 G (3) = 7 of 10 reads: ratio 0.7
    pri = make_pri(T=reads(b"A", 4) + reads(b"G", 3) + reads(b"T", 3))
    pri.update_top_alleles(min_top12alleles_ratio=min_ratio)
    assert pri.filterflag() is filtered
    assert st(pri)["top12alleles_ratio"] == pytest.approx(0.7, rel=1e-6)


def test_update_top_alleles_homozygous_alt_kept():
    pri = make_pri(T=reads(b"G", 6))
    pri.update_top_alleles()
    s = st(pri)
    assert (s["top1allele"], s["top2allele"]) == (b"G", b"A")
    assert pri.filterflag() is False


def test_update_top_alleles_single_alt_read_is_filtered():
    # top1 G (1 read) is dropped too: 1 - 0 < 2, nothing left
    pri = make_pri(T=reads(b"G", 1))
    pri.update_top_alleles()
    assert pri.filterflag() is True
    assert st(pri)["n_reads_T"][b"G"] == 0


def test_update_top_alleles_control_only_is_filtered():
    pri = make_pri(C=[(b"G", 30, 0)] * 4)
    pri.update_top_alleles()
    assert pri.filterflag() is True


def test_update_top_alleles_insertion_skips_count_filter():
    # multi-base alleles skip the minimum-count rule
    pri = make_pri(T=reads(b"A", 8) + reads(b"AT", 1))
    pri.update_top_alleles()
    s = st(pri)
    assert (s["top1allele"], s["top2allele"]) == (b"A", b"AT")
    assert pri.filterflag() is False


# ------------------------------------
# top12alleles
# ------------------------------------

def test_top12alleles_prints(capsys):
    pri = make_pri(T=[(b"A", 30, 0, False), (b"A", 31, 0, False),
                      (b"G", 40, 1, False), (b"G", 41, 1, False)],
                   C=[(b"G", 25, 0)])
    pri.update_top_alleles()
    pri.top12alleles()
    assert capsys.readouterr().out == (
        "100 b'A'\n"
        "Top1allele b'A' Treatment [30, 31] Control []\n"
        "Top2allele b'G' Treatment [40, 41] Control [25]\n")


def test_top12alleles_before_update_raises():
    with pytest.raises(KeyError):
        make_pri(T=reads(b"A", 2)).top12alleles()


# ------------------------------------
# call_GT
# ------------------------------------

def test_call_GT_homozygous_alt_without_minor_reads():
    T = [(b"G", 30, 0, False)] * 4 + [(b"G", 30, 1, False)] * 2
    pri = make_pri(T=T)
    pri.update_top_alleles()
    pri.call_GT()
    s = st(pri)
    m = ref_models([30] * 6, [], [], [])
    assert_models(s, m)
    # no top2 read at all: 1/1 when min(other BICs) - BIC_homo_major >= 2
    dBIC = (min(m["heter_noAS"][1], m["heter_AS"][1], m["homo_minor"][1])
            - m["homo_major"][1])
    assert dBIC >= 2
    assert (s["type"], s["GT"], s["alt_allele"]) == ("homo", "1/1", b"G")
    PL_11 = phred(m["homo_major"][0])
    PL_00 = max(0, phred(m["homo_minor"][0]) - PL_11)
    PL_01 = max(0, phred(max(m["heter_noAS"][0], m["heter_AS"][0])) - PL_11)
    assert s["GQ"] == pytest.approx(min(PL_00, PL_01), rel=1e-9)
    assert pri.filterflag() is False


def test_call_GT_homozygous_alt_too_shallow_is_filtered():
    # 2 alt reads: deltaBIC < 2
    pri = make_pri(T=reads(b"G", 2))
    pri.update_top_alleles()
    pri.call_GT()
    m = ref_models([30] * 2, [], [], [])
    assert (min(m["heter_noAS"][1], m["heter_AS"][1], m["homo_minor"][1])
            - m["homo_major"][1]) < 2
    assert pri.filterflag() is True
    assert st(pri)["GT"] == "unsure"


def test_call_GT_heter_noAS():
    T = (reads(b"A", 3, strand=0) + reads(b"A", 2, strand=1)
         + reads(b"G", 2, strand=0) + reads(b"G", 3, strand=1))
    pri = make_pri(T=T)
    pri.update_top_alleles()
    pri.call_GT()
    s = st(pri)
    m = ref_models([30] * 5, [], [30] * 5, [])
    assert_models(s, m)
    b = {k: v[1] for k, v in m.items()}
    assert b["heter_noAS"] + 2 <= min(b["homo_major"], b["homo_minor"], b["heter_AS"])
    assert (s["type"], s["GT"], s["alt_allele"]) == ("heter_noAS", "0/1", b"G")
    PL_01 = phred(m["heter_noAS"][0])
    GQ = min(phred(m["homo_minor"][0]) - PL_01, phred(m["homo_major"][0]) - PL_01)
    assert s["GQ"] == pytest.approx(GQ, rel=1e-9)


def test_call_GT_two_alt_alleles_is_1_2():
    pri = make_pri(T=reads(b"G", 5) + reads(b"T", 5))
    pri.update_top_alleles()
    pri.call_GT()
    s = st(pri)
    assert (s["top1allele"], s["top2allele"]) == (b"G", b"T")
    assert (s["type"], s["GT"], s["alt_allele"]) == ("heter_noAS", "1/2", b"G,T")
    assert "MT=SNV,SNV;" in pri.to_vcf()


def test_call_GT_reference_is_top2():
    pri = make_pri(T=reads(b"G", 6) + reads(b"A", 4))
    pri.update_top_alleles()
    pri.call_GT()
    s = st(pri)
    assert (s["top1allele"], s["top2allele"]) == (b"G", b"A")
    assert s["type"].startswith("heter")
    assert (s["GT"], s["alt_allele"]) == ("0/1", b"G")


def test_call_GT_heter_AS():
    pri = make_pri(T=reads(b"A", 30, q=40) + reads(b"G", 8, q=40))
    pri.update_top_alleles()
    pri.call_GT(max_allowed_ar=0.95)
    s = st(pri)
    m = ref_models([40] * 30, [], [40] * 8, [], max_ar=0.95)
    assert_models(s, m)
    b = {k: v[1] for k, v in m.items()}
    assert b["heter_AS"] + 2 <= min(b["homo_major"], b["homo_minor"], b["heter_noAS"])
    assert (s["type"], s["GT"], s["alt_allele"]) == ("heter_AS", "0/1", b"G")
    PL_01 = phred(m["heter_AS"][0])
    GQ = min(phred(m["homo_minor"][0]) - PL_01, phred(m["homo_major"][0]) - PL_01)
    assert s["GQ"] == pytest.approx(GQ, rel=1e-9)


def test_call_GT_heter_unsure():
    pri = make_pri(T=reads(b"A", 7) + reads(b"G", 3))
    pri.update_top_alleles()
    pri.call_GT()
    s = st(pri)
    m = ref_models([30] * 7, [], [30] * 3, [])
    assert_models(s, m)
    b = {k: v[1] for k, v in m.items()}
    # neither heterozygous model beats the other by 2, both beat homo by 2
    assert abs(b["heter_AS"] - b["heter_noAS"]) < 2
    assert b["heter_AS"] + 2 <= min(b["homo_major"], b["homo_minor"])
    assert (s["type"], s["GT"]) == ("heter_unsure", "0/1")
    PL_01 = phred(max(m["heter_noAS"][0], m["heter_AS"][0]))
    GQ = min(phred(m["homo_minor"][0]) - PL_01, phred(m["homo_major"][0]) - PL_01)
    assert s["GQ"] == pytest.approx(GQ, rel=1e-9)


def test_call_GT_homo_ref():
    # 20 reference reads and 2 low-quality alt reads: homo_major wins
    pri = make_pri(T=reads(b"A", 20) + reads(b"G", 2, q=21))
    pri.update_top_alleles()
    assert pri.filterflag() is False
    pri.call_GT()
    s = st(pri)
    m = ref_models([30] * 20, [], [21] * 2, [])
    assert_models(s, m)
    b = {k: v[1] for k, v in m.items()}
    assert b["homo_major"] < min(b["homo_minor"], b["heter_noAS"], b["heter_AS"])
    assert (s["type"], s["GT"]) == ("homo_ref", "0/0")
    assert pri.filterflag() is True


def test_call_GT_uses_control_reads():
    T = reads(b"A", 5) + reads(b"G", 5)
    C = [(b"A", 30, 0)] * 3 + [(b"G", 30, 1)] * 2
    pri = make_pri(T=T, C=C)
    pri.update_top_alleles()
    pri.call_GT()
    assert_models(st(pri), ref_models([30] * 5, [30] * 3, [30] * 5, [30] * 2))


@pytest.mark.parametrize("alt,mt", [(b"*", "Deletion"), (b"AT", "Insertion")])
def test_call_GT_indel_mutation_type(alt, mt):
    q = 93 if alt == b"*" else 30
    pri = make_pri(T=reads(b"A", 5) + reads(alt, 5, q=q))
    pri.update_top_alleles()
    pri.call_GT()
    s = st(pri)
    assert s["alt_allele"] == alt and s["GT"] == "0/1"
    assert ";MT=%s;" % mt in pri.to_vcf()


def test_call_GT_on_filtered_does_nothing():
    pri = make_pri(T=reads(b"A", 5))
    pri.update_top_alleles()          # homo_ref, filtered
    pri.call_GT()
    s = st(pri)
    assert s["GT"] == "unsure" and s["lnL_homo_major"] == 0


# ------------------------------------
# apply_GQ_cutoff / apply_deltaBIC_cutoff
# ------------------------------------

@pytest.mark.parametrize("type_,GQ,filtered", [
    ("homo", 49.9, True),
    ("homo", 50.0, False),
    ("heter_noAS", 99.0, True),
    ("heter_AS", 100.0, False),
    ("heter_unsure", 150.0, False),
])
def test_apply_GQ_cutoff_defaults(type_, GQ, filtered):
    pri = PosReadsInfo(1, b"A")
    set_state(pri, type=type_, GQ=GQ)
    pri.apply_GQ_cutoff()
    assert pri.filterflag() is filtered


def test_apply_GQ_cutoff_custom_and_already_filtered():
    pri = PosReadsInfo(1, b"A")
    set_state(pri, type="homo", GQ=3.0)
    pri.apply_GQ_cutoff(3, 0)
    assert pri.filterflag() is False
    pri.apply_GQ_cutoff(4, 0)
    assert pri.filterflag() is True
    # an already filtered position is left alone, even without a type
    pri2 = PosReadsInfo(1, b"A")
    set_state(pri2, filterout=True)
    pri2.apply_GQ_cutoff()
    assert pri2.filterflag() is True


def test_apply_deltaBIC_cutoff():
    pri = make_pri(T=reads(b"A", 5) + reads(b"G", 5))
    pri.update_top_alleles()
    pri.call_GT()
    m = ref_models([30] * 5, [], [30] * 5, [])
    # deltaBIC of heter_noAS = min(BIC homo major, homo minor) - BIC noAS
    dBIC = min(m["homo_major"][1], m["homo_minor"][1]) - m["heter_noAS"][1]
    pri.apply_deltaBIC_cutoff()            # default 10
    assert pri.filterflag() is False
    pri.apply_deltaBIC_cutoff(dBIC - 0.01)
    assert pri.filterflag() is False
    pri.apply_deltaBIC_cutoff(dBIC + 0.01)
    assert pri.filterflag() is True


@pytest.mark.parametrize("cutoff,filtered", [(0, False), (10, True), (-1, False)])
def test_apply_deltaBIC_cutoff_on_fresh_object(cutoff, filtered):
    # deltaBIC starts at 0
    pri = PosReadsInfo(1, b"A")
    pri.apply_deltaBIC_cutoff(cutoff)
    assert pri.filterflag() is filtered


# ------------------------------------
# to_vcf / toVariant
# ------------------------------------

def expected_noAS_vcf(float_cast=float):
    """Expected to_vcf()/toVCF() text of the 5 A + 5 G example."""
    m = ref_models([30] * 5, [], [30] * 5, [])
    b = {k: float_cast(v[1]) for k, v in m.items()}
    dBIC = float_cast(min(m["homo_major"][1], m["homo_minor"][1])
                      - m["heter_noAS"][1])
    PL_01 = phred(m["heter_noAS"][0])
    PL_00 = phred(m["homo_minor"][0]) - PL_01
    PL_11 = phred(m["homo_major"][0]) - PL_01
    GQ = min(PL_00, PL_11)
    info = ("M=heter_noAS;MT=SNV;DPT=10;DPC=0;DP1T=5A;DP2T=5G;DP1C=0A;DP2C=0G;"
            "SB=3,2,2,3;DBIC=%.2f;BICHOMOMAJOR=%.2f;BICHOMOMINOR=%.2f;"
            "BICHETERNOAS=%.2f;BICHETERAS=%.2f;AR=0.50"
            % (dBIC, b["homo_major"], b["homo_minor"], b["heter_noAS"],
               b["heter_AS"]))
    sample = "0/1:10:%d:%d,0,%d" % (GQ, PL_00, PL_11)
    return "\t".join(("A", "G", "%d" % GQ, ".", info, "GT:DP:GQ:PL", sample))


def called_noAS():
    T = (reads(b"A", 3, strand=0) + reads(b"A", 2, strand=1)
         + reads(b"G", 2, strand=0) + reads(b"G", 3, strand=1))
    pri = make_pri(T=T)
    pri.update_top_alleles()
    pri.call_GT()
    return pri


def test_to_vcf():
    assert called_noAS().to_vcf() == expected_noAS_vcf()


def test_toVariant_returns_Variant_with_same_columns():
    v = called_noAS().toVariant()
    assert isinstance(v, Variant)
    # Variant keeps BIC values in C floats
    assert v.toVCF() == expected_noAS_vcf(lambda x: float(np.float32(x)))
    assert (v["ref_allele"], v["alt_allele"]) == ("A", "G")
    assert (v["top1allele"], v["top2allele"]) == ("A", "G")


def test_to_vcf_with_control_depths():
    T = reads(b"A", 5) + reads(b"G", 5)
    C = [(b"A", 30, 0)] * 3 + [(b"G", 30, 1)] * 2 + [(b"T", 30, 0)]
    pri = make_pri(T=T, C=C)
    pri.update_top_alleles()
    pri.call_GT()
    fields = pri.to_vcf().split("\t")
    info = dict(x.split("=") for x in fields[4].split(";"))
    assert (info["DPT"], info["DPC"]) == ("10", "6")
    assert (info["DP1T"], info["DP2T"]) == ("5A", "5G")
    assert (info["DP1C"], info["DP2C"]) == ("3A", "2G")
    # DP counts every read kept at the position, including the third allele
    assert fields[6].split(":")[1] == "16"


# ------------------------------------
# __getstate__ / __setstate__
# ------------------------------------

def test_state_round_trip():
    pri = called_noAS()
    other = PosReadsInfo(0, b"N")
    other.__setstate__(pri.__getstate__())
    assert other.__getstate__() == pri.__getstate__()


def test_pickle_round_trip():
    pri = called_noAS()
    assert pickle.loads(pickle.dumps(pri)).__getstate__() == pri.__getstate__()
