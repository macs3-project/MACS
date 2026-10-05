#!/usr/bin/env python

"""Module Description: Test functions in MACS3.Signal.VariantStat that
compute genotype log-likelihoods, BIC values and genotype qualities
for callvar.

The expected values come from direct numpy/scipy implementations of the
likelihood formulas: a base with Phred quality q is wrong with
probability e = 10^(-q/10); a read carries the top1 allele with
probability k/tn when k of the tn reads are drawn from the top1 allele;
the number of top1 copies k follows a binomial distribution with
probability 0.5 (no allele-specific binding) or with probability k/tn
capped to [1 - maxAR, maxAR] (allele-specific binding).

The C-only helpers GreedyMaxFunctionAS, GreedyMaxFunctionNoAS and
calculate_ln are covered through CalModel_Heter_noAS and
CalModel_Heter_AS.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import math

import numpy as np
import pytest
from scipy.stats import binom

from MACS3.Signal.VariantStat import (CalModel_Homo,
                                      CalModel_Heter_noAS,
                                      CalModel_Heter_AS,
                                      calculate_GQ,
                                      calculate_GQ_heterASsig)


# ------------------------------------
# reference implementations
# ------------------------------------

def bq(*values):
    """int32 array of base qualities, as PosReadsInfo.call_GT builds it."""
    return np.array(values, dtype="i4")


EMPTY = bq()


def err_rate(quals):
    """Phred quality -> probability that the base call is wrong."""
    return 10.0 ** (-np.asarray(quals, dtype=float) / 10.0)


def ref_homo(t1T, t1C, t2T, t2C):
    """Homozygous model: every top1 base is right, every top2 base is a
    sequencing error. No free parameter, so BIC = -2 lnL."""
    e1 = err_rate(np.concatenate([t1T, t1C]))
    e2 = err_rate(np.concatenate([t2T, t2C]))
    lnL = float(np.log1p(-e1).sum() + np.log(e2).sum())
    return lnL, -2.0 * lnL


def ref_lnL_k(me, ne, k, p):
    """ln L when k of the tn reads come from the top1 allele.

    The count k is binomial(tn, p); a top1 base is seen with probability
    (1-e)*f + e*(1-f) and a top2 base with (1-e)*(1-f) + e*f, f = k/tn.
    """
    tn = len(me) + len(ne)
    f = k / tn
    em = err_rate(me)
    en = err_rate(ne)
    lnL = binom.logpmf(k, tn, p)
    lnL += np.log((1 - em) * f + em * (1 - f)).sum()
    lnL += np.log((1 - en) * (1 - f) + en * f).sum()
    return float(lnL)


def ref_noAS_max(me, ne):
    """Maximum over k of the no-allele-specific likelihood (p = 0.5).

    The log-likelihood is concave in k, so MACS3's hill climb from k = m
    reaches this maximum whenever both alleles have reads.
    """
    tn = len(me) + len(ne)
    return max(ref_lnL_k(me, ne, k, 0.5) for k in range(tn + 1))


def f32(x):
    """Value of a Python float after a round trip through C float, as
    for the ``max_allowed_ar: cython.float`` arguments."""
    return float(np.float32(x))


def ref_AS_max(me, ne, max_ar=0.99):
    """Maximum over k of the allele-specific likelihood, p = k/tn capped
    to [1 - maxAR, maxAR]."""
    tn = len(me) + len(ne)
    mar = f32(max_ar)
    return max(ref_lnL_k(me, ne, k, min(max(k / tn, 1 - mar), mar))
               for k in range(tn + 1))


def ref_noAS_model(t1T, t1C, t2T, t2C):
    """(lnL, BIC) of the heterozygous model without allele-specific
    binding: one free k for treatment and one for control."""
    tn_T = len(t1T) + len(t2T)
    tn_C = len(t1C) + len(t2C)
    lnL = ref_noAS_max(t1T, t2T)
    penalty = math.log(tn_T)
    if tn_C:
        lnL += ref_noAS_max(t1C, t2C)
        penalty += math.log(tn_C)
    return lnL, -2 * lnL + penalty


def ref_AS_model(t1T, t1C, t2T, t2C, max_ar=0.99):
    """(lnL, BIC) of the heterozygous model with allele-specific binding:
    k and the allele ratio are free in treatment (2 log tn_T), control
    follows the no-AS model (log tn_C)."""
    tn_T = len(t1T) + len(t2T)
    tn_C = len(t1C) + len(t2C)
    lnL = ref_AS_max(t1T, t2T, max_ar)
    penalty = 2 * math.log(tn_T)
    if tn_C:
        lnL += ref_noAS_max(t1C, t2C)
        penalty += math.log(tn_C)
    return lnL, -2 * lnL + penalty


def ref_GQ(lnL1, lnL2, lnL3):
    """GQ = -10 log10((L2+L3)/(L1+L2+L3)) with L relative to L1 and each
    of L2, L3 clipped into [1e-110, 1]; truncated to int."""
    L2 = min(max(math.exp(lnL2 - lnL1), 1e-110), 1.0)
    L3 = min(max(math.exp(lnL3 - lnL1), 1e-110), 1.0)
    return int(-10 * math.log10((L2 + L3) / (1 + L2 + L3)))


def ref_GQ_ASsig(lnL1, lnL2):
    """-10 log10(L2/(L1+L2)) with L2/L1 clipped into [1e-110, 1]."""
    L2 = min(max(math.exp(lnL2 - lnL1), 1e-110), 1.0)
    return int(-10 * math.log10(L2 / (1 + L2)))


# ------------------------------------
# CalModel_Homo
# ------------------------------------

HOMO_CASES = [
    # (top1_T, top1_C, top2_T, top2_C)
    (bq(30, 30, 30), EMPTY, EMPTY, EMPTY),
    (bq(30, 35, 40, 21), EMPTY, bq(25, 33), EMPTY),
    (bq(40, 40), bq(30), bq(22), bq(37, 38)),
    (bq(93), EMPTY, bq(93), EMPTY),          # highest Phred in BAM
    (bq(1, 2, 3), bq(4), bq(1), bq(2)),      # lowest non-zero qualities
    (EMPTY, EMPTY, bq(30, 20), EMPTY),       # no top1 base at all
    (bq(*range(21, 61)), bq(*range(30, 40)), bq(*range(21, 31)), EMPTY),
]


@pytest.mark.parametrize("arrays", HOMO_CASES)
def test_CalModel_Homo_matches_reference(arrays):
    lnL, BIC = CalModel_Homo(*arrays)
    exp_lnL, exp_BIC = ref_homo(*arrays)
    assert lnL == pytest.approx(exp_lnL, rel=1e-12, abs=1e-12)
    assert BIC == pytest.approx(exp_BIC, rel=1e-12, abs=1e-12)
    assert BIC == -2 * lnL


def test_CalModel_Homo_by_hand():
    # one top1 base at Q30: ln(1 - 0.001) = -0.0010005003335835335
    lnL, BIC = CalModel_Homo(bq(30), EMPTY, EMPTY, EMPTY)
    assert lnL == pytest.approx(-0.0010005003335835335, rel=1e-12)
    assert BIC == pytest.approx(0.002001000667167067, rel=1e-12)
    # one top2 base at Q20 is an error with probability 0.01:
    # ln(0.01) = -4.605170185988091
    lnL, BIC = CalModel_Homo(EMPTY, EMPTY, bq(20), EMPTY)
    assert lnL == pytest.approx(-4.605170185988091, rel=1e-12)
    assert BIC == pytest.approx(9.210340371976182, rel=1e-12)


def test_CalModel_Homo_empty_is_zero():
    lnL, BIC = CalModel_Homo(EMPTY, EMPTY, EMPTY, EMPTY)
    assert (lnL, BIC) == (0.0, 0.0)


def test_CalModel_Homo_swapping_treatment_and_control_is_symmetric():
    # treatment and control enter the homozygous model the same way
    a = CalModel_Homo(bq(30, 31), bq(40), bq(25), bq(22, 23))
    b = CalModel_Homo(bq(40), bq(30, 31), bq(22, 23), bq(25))
    assert a[0] == pytest.approx(b[0], rel=1e-12)


def test_CalModel_Homo_needs_arrays():
    # the quality vectors must be numpy arrays (their .shape is read)
    with pytest.raises(AttributeError, match="shape"):
        CalModel_Homo([30], EMPTY, EMPTY, EMPTY)


# ------------------------------------
# CalModel_Heter_noAS
# ------------------------------------

HETER_CASES = [
    (bq(30, 30, 30, 30, 30), EMPTY, bq(30, 30, 30, 30, 30), EMPTY),
    (bq(40, 40, 40, 40, 40, 40, 40), EMPTY, bq(40, 40, 40), EMPTY),
    (bq(35, 30, 25, 38, 40, 22), EMPTY, bq(33, 21), EMPTY),
    (bq(30, 30, 30, 30), bq(30, 30), bq(30, 30), bq(30, 30)),
    (bq(40, 40, 40, 40, 40, 40, 40, 40), bq(25), bq(40, 40), bq(25, 26, 27)),
    (bq(*([30] * 30)), EMPTY, bq(*([30] * 10)), EMPTY),
]


@pytest.mark.parametrize("arrays", HETER_CASES)
def test_CalModel_Heter_noAS_matches_reference(arrays):
    lnL, BIC = CalModel_Heter_noAS(*arrays)
    exp_lnL, exp_BIC = ref_noAS_model(*arrays)
    assert lnL == pytest.approx(exp_lnL, rel=1e-10)
    assert BIC == pytest.approx(exp_BIC, rel=1e-10)


def test_CalModel_Heter_noAS_balanced_by_hand():
    # 5 + 5 reads at Q30, best k = 5, f = 0.5: every base term is ln(0.5)
    # lnL = ln C(10,5) + 10 ln 0.5 + 10 ln 0.5 = ln 252 + 20 ln 0.5
    lnL, BIC = CalModel_Heter_noAS(bq(*[30] * 5), EMPTY, bq(*[30] * 5), EMPTY)
    exp = math.log(252) + 20 * math.log(0.5)
    assert lnL == pytest.approx(exp, rel=1e-12)
    assert BIC == pytest.approx(-2 * exp + math.log(10), rel=1e-12)


def test_CalModel_Heter_noAS_single_read():
    # tn = 1: k = 1 wins; lnL = ln 0.5 + ln(1 - 0.001), BIC adds ln 1 = 0
    lnL, BIC = CalModel_Heter_noAS(bq(30), EMPTY, EMPTY, EMPTY)
    exp = math.log(0.5) + math.log1p(-0.001)
    assert lnL == pytest.approx(exp, rel=1e-12)
    assert BIC == pytest.approx(-2 * exp, rel=1e-12)


def test_CalModel_Heter_noAS_all_top1():
    # When one allele has no read, GreedyMaxFunctionNoAS fixes k at the
    # observed count instead of searching over k (by design).
    # m == tn: k = tn, f = 1; lnL = 3 ln 0.5 + 3 ln(1 - 0.001)
    lnL, BIC = CalModel_Heter_noAS(bq(30, 30, 30), EMPTY, EMPTY, EMPTY)
    exp = 3 * math.log(0.5) + 3 * math.log1p(-0.001)
    assert lnL == pytest.approx(exp, rel=1e-12)
    assert BIC == pytest.approx(-2 * exp + math.log(3), rel=1e-12)


def test_CalModel_Heter_noAS_no_top1():
    # m == 0: k = 0, f = 0; each top2 base at Q20 contributes ln(0.99)
    lnL, BIC = CalModel_Heter_noAS(EMPTY, EMPTY, bq(20, 20), EMPTY)
    exp = 2 * math.log(0.5) + 2 * math.log(0.99)
    assert lnL == pytest.approx(exp, rel=1e-12)
    assert BIC == pytest.approx(-2 * exp + math.log(2), rel=1e-12)


def test_CalModel_Heter_noAS_control_adds_penalty():
    t1T, t2T = bq(30, 30, 30), bq(30, 30)
    t1C, t2C = bq(30), bq(30, 30, 30)
    lnL_T, BIC_T = CalModel_Heter_noAS(t1T, EMPTY, t2T, EMPTY)
    lnL_C, BIC_C = CalModel_Heter_noAS(t1C, EMPTY, t2C, EMPTY)
    lnL, BIC = CalModel_Heter_noAS(t1T, t1C, t2T, t2C)
    assert lnL == pytest.approx(lnL_T + lnL_C, rel=1e-12)
    assert BIC == pytest.approx(-2 * (lnL_T + lnL_C) + math.log(5) + math.log(4),
                                rel=1e-12)


@pytest.mark.parametrize("control", [(EMPTY, EMPTY), (bq(30), bq(30))])
def test_CalModel_Heter_noAS_no_treatment_reads_raises(control):
    with pytest.raises(Exception, match="Total number of treatment reads is 0!"):
        CalModel_Heter_noAS(EMPTY, control[0], EMPTY, control[1])


# ------------------------------------
# CalModel_Heter_AS
# ------------------------------------

@pytest.mark.parametrize("max_ar", [0.99, 0.95])
@pytest.mark.parametrize("arrays", HETER_CASES)
def test_CalModel_Heter_AS_matches_reference(arrays, max_ar):
    lnL, BIC = CalModel_Heter_AS(*arrays, max_ar)
    exp_lnL, exp_BIC = ref_AS_model(*arrays, max_ar=max_ar)
    assert lnL == pytest.approx(exp_lnL, rel=1e-10)
    assert BIC == pytest.approx(exp_BIC, rel=1e-10)


def test_CalModel_Heter_AS_default_max_ar_is_0_99():
    arrays = (bq(*[40] * 7), EMPTY, bq(40, 40, 40), EMPTY)
    assert CalModel_Heter_AS(*arrays) == CalModel_Heter_AS(*arrays, 0.99)


@pytest.mark.parametrize("max_ar", [0.99, 0.95, 0.8])
def test_CalModel_Heter_AS_all_top1(max_ar):
    # m == tn: k = tn, r = 1 is capped to maxAR (as a C float):
    # lnL = tn ln(maxAR) + sum ln(1 - e); BIC = -2 lnL + 2 ln tn
    t1 = bq(30, 40, 20)
    lnL, BIC = CalModel_Heter_AS(t1, EMPTY, EMPTY, EMPTY, max_ar)
    exp = 3 * math.log(f32(max_ar)) + float(np.log1p(-err_rate(t1)).sum())
    assert lnL == pytest.approx(exp, rel=1e-12)
    assert BIC == pytest.approx(-2 * exp + 2 * math.log(3), rel=1e-12)


def test_CalModel_Heter_AS_single_read_default():
    # tn = 1: k = 1 with r = 1 capped to 0.99
    lnL, BIC = CalModel_Heter_AS(bq(30), EMPTY, EMPTY, EMPTY)
    exp = math.log(f32(0.99)) + math.log1p(-0.001)
    assert lnL == pytest.approx(exp, rel=1e-12)
    assert BIC == pytest.approx(-2 * exp, rel=1e-12)


def test_CalModel_Heter_AS_control_uses_noAS_model():
    t1T, t2T = bq(*[40] * 6), bq(40, 40)
    t1C, t2C = bq(30, 30), bq(30, 30, 30)
    lnL_T, _ = CalModel_Heter_AS(t1T, EMPTY, t2T, EMPTY, 0.95)
    lnL_C, _ = CalModel_Heter_noAS(t1C, EMPTY, t2C, EMPTY)
    lnL, BIC = CalModel_Heter_AS(t1T, t1C, t2T, t2C, 0.95)
    assert lnL == pytest.approx(lnL_T + lnL_C, rel=1e-12)
    assert BIC == pytest.approx(-2 * (lnL_T + lnL_C) + 2 * math.log(8)
                                + math.log(5), rel=1e-12)


def test_CalModel_Heter_AS_balanced_equals_noAS_likelihood():
    # with 5 + 5 reads the best allele ratio is 0.5, so both heterozygous
    # models reach the same lnL; AS pays one extra ln(tn) in BIC
    arrays = (bq(*[30] * 5), EMPTY, bq(*[30] * 5), EMPTY)
    lnL_AS, BIC_AS = CalModel_Heter_AS(*arrays)
    lnL_no, BIC_no = CalModel_Heter_noAS(*arrays)
    assert lnL_AS == pytest.approx(lnL_no, rel=1e-12)
    assert BIC_AS - BIC_no == pytest.approx(math.log(10), rel=1e-12)


@pytest.mark.parametrize("control", [(EMPTY, EMPTY), (bq(30), bq(30))])
def test_CalModel_Heter_AS_no_treatment_reads_raises(control):
    with pytest.raises(Exception, match="Total number of treatment reads is 0!"):
        CalModel_Heter_AS(EMPTY, control[0], EMPTY, control[1])


# ------------------------------------
# calculate_GQ
# ------------------------------------

@pytest.mark.parametrize("lnLs", [
    (0.0, 0.0, 0.0),            # L2 = L3 = 1: -10 log10(2/3) = 1.76
    (0.0, -10.0, -20.0),        # 43.4
    (-5.0, -7.5, -6.25),
    (0.0, 5.0, -1000.0),        # L2 capped at 1, L3 at 1e-110: 3.01
    (0.0, -1e6, -1e6),          # both capped at 1e-110: 1096.98
    (-123.4, -130.0, -150.0),
])
def test_calculate_GQ(lnLs):
    result = calculate_GQ(*lnLs)
    assert isinstance(result, int)
    assert result == ref_GQ(*lnLs)


def test_calculate_GQ_by_hand():
    assert calculate_GQ(0.0, 0.0, 0.0) == 1
    assert calculate_GQ(0.0, 5.0, -1000.0) == 3
    assert calculate_GQ(0.0, -1e6, -1e6) == 1096


# ------------------------------------
# calculate_GQ_heterASsig
# ------------------------------------

@pytest.mark.parametrize("lnLs", [
    (0.0, 0.0),             # L2 = 1: -10 log10(1/2) = 3.01
    (0.0, 3.0),             # L2 capped at 1
    (0.0, -10.0),           # 43.4
    (-50.0, -52.5),
    (0.0, -250.0),          # L2 = 2.7e-109, still above the 1e-110 floor
])
def test_calculate_GQ_heterASsig(lnLs):
    result = calculate_GQ_heterASsig(*lnLs)
    assert isinstance(result, int)
    assert result == ref_GQ_ASsig(*lnLs)


def test_calculate_GQ_heterASsig_by_hand():
    assert calculate_GQ_heterASsig(0.0, 0.0) == 3
    assert calculate_GQ_heterASsig(0.0, -10.0) == 43


def test_calculate_GQ_heterASsig_at_floor_is_255():
    """Once L2/L1 reaches its 1e-110 floor the score is 255, below the
    scores of weaker evidence just above the floor (1098 at a log
    likelihood difference of 253).

    Pins the current output. The score is int(-4.34294 ln(L2/(1+L2)))
    while L2/(1+L2) > 1e-110 and 255 otherwise, and the clip makes the
    ratio exactly 1e-110 for differences beyond ln(1e110) = 253.3. The
    function has no docstring and nothing in MACS3 calls it, so whether
    255 is meant as a cap for every score is not documented.
    """
    diffs = (10, 100, 250, 253)
    assert [calculate_GQ_heterASsig(0.0, -d) for d in diffs] == \
        [ref_GQ_ASsig(0.0, -d) for d in diffs] == [43, 434, 1085, 1098]
    assert [calculate_GQ_heterASsig(0.0, -d) for d in (254, 1000)] == \
        [255, 255]
