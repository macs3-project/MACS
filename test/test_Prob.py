#!/usr/bin/env python
# Time-stamp: <2025-09-29 15:04:42 Tao Liu>

"""Module Description: Test functions to calculate probabilities.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest

from math import log10
from MACS3.Signal.Prob import (factorial,
                               poisson_cdf,
                               chisq_pvalue_e,
                               chisq_logp_e,
                               binomial_cdf,
                               binomial_cdf_inv,
                               binomial_pdf)

import math

import numpy as np
import pytest
from scipy import stats
from scipy.special import logsumexp

from MACS3.Signal.Prob import (pnorm,
                               pnorm2,
                               poisson_cdf_inv,
                               poisson_cdf_Q_inv,
                               poisson_pdf,
                               binomial_sf,
                               pduplication)

# ------------------------------------
# Main function
# ------------------------------------


class Test_factorial(unittest.TestCase):

    def setUp(self):
        self.n1 = 100
        self.n2 = 10
        self.n3 = 1

    def test_factorial_big_n1(self):
        expect = 9.332622e+157
        result = factorial(self.n1)
        self.assertTrue(abs(result - expect) < 1e-5*result)

    def test_factorial_median_n2(self):
        expect = 3628800
        result = factorial(self.n2)
        self.assertEqual(result, expect)

    def test_factorial_small_n3(self):
        expect = 1
        result = factorial(self.n3)
        self.assertEqual(result, expect)


class Test_poisson_cdf(unittest.TestCase):

    def setUp(self):
        # n, lam
        self.n1 = (80, 100)
        self.n2 = (200, 100)
        self.n3 = (100, 1000)
        self.n4 = (1500, 1000)

    def test_poisson_cdf_n1(self):
        expect = (round(0.9773508, 5), round(0.02264918, 5))
        result = (round(poisson_cdf(self.n1[0], self.n1[1], False), 5),
                  round(poisson_cdf(self.n1[0], self.n1[1], True), 5))
        self.assertEqual(result, expect)

    def test_poisson_cdf_n2(self):
        expect = (round(log10(4.626179e-19), 4),
                  round(log10(1), 4))
        result = (round(log10(poisson_cdf(self.n2[0], self.n2[1], False)), 4),
                  round(log10(poisson_cdf(self.n2[0], self.n2[1], True)), 4))
        self.assertEqual(result, expect)

    def test_poisson_cdf_n3(self):
        expect = (round(log10(1), 2),
                  round(log10(6.042525e-293), 2))
        result = (round(poisson_cdf(self.n3[0], self.n3[1], False, True), 2),
                  round(poisson_cdf(self.n3[0], self.n3[1], True, True), 2))
        self.assertEqual(result, expect)

    def test_poisson_cdf_n4(self):
        expect = (round(log10(2.097225e-49), 4),
                  round(log10(1), 4))
        result = (round(log10(poisson_cdf(self.n4[0], self.n4[1], False)), 4),
                  round(log10(poisson_cdf(self.n4[0], self.n4[1], True)), 4))
        self.assertEqual(result, expect)


class Test_chisq_p_e(unittest.TestCase):
    """Test chisq pvalue calculation -- assuming df is an even number. We
    only implemented even number pchisq for upper tail. Because this
    is the function we need to combine p-values using fisher's method

    """
    def setUp(self):
        # x, k, p(upper), -log p upper, -log10 p upper
        self.c = ((10, 2, 0.006737947, 5, 2.171472),
                  (100, 2, 1.92875e-22, 50, 21.71472),
                  (1000, 22, 1.956374e-197, 452.9382, 196.7085),
                  (10, 4, 0.04042768, 3.208241, 1.393321),
                  (100, 8, 4.269159e-18, 39.99511, 17.36966),
                  (1000, 80, 6.889598e-159, 364.181, 158.1618),
                  (54, 6, 7.377151e-10, 21.02746, 9.132111),
                  (565, 10, 5.518772e-115, 263.0891, 114.2582),
                  (7765, 12, 0, 3845.965, 1670.2814),
                  )

    def test_chisq_p(self):
        expect = [round(x[2], 4) for x in self.c]
        result = [round(chisq_pvalue_e(x[0], x[1]), 4) for x in self.c]
        self.assertEqual(result, expect)

    def test_chisq_logp(self):
        expect = [round(x[3], 4) for x in self.c]
        result = [round(chisq_logp_e(x[0], x[1]), 4) for x in self.c]
        self.assertEqual(result, expect)

    def test_chisq_log10p(self):
        expect = [round(x[4], 4) for x in self.c]
        result = [round(chisq_logp_e(x[0], x[1], log10=True), 4) for x in self.c]
        self.assertEqual(result, expect)


class Test_binomial_cdf(unittest.TestCase):
    def setUp(self):
        # x, a,  b
        self.n1 = (20, 1000, 0.01)
        self.n2 = (200, 1000, 0.01)

    def test_binomial_cdf_n1(self):
        expect = (round(0.001496482, 5), round(0.9985035, 5))
        result = (round(binomial_cdf(self.n1[0], self.n1[1], self.n1[2], False), 5),
                  round(binomial_cdf(self.n1[0], self.n1[1], self.n1[2], True), 5))
        self.assertEqual(result, expect)

    def test_binomial_cdf_n2(self):
        expect = (round(log10(8.928717e-190), 4),
                  round(log10(1), 4))
        result = (round(log10(binomial_cdf(self.n2[0], self.n2[1], self.n2[2], False)), 4),
                  round(log10(binomial_cdf(self.n2[0], self.n2[1], self.n2[2], True)), 4))
        self.assertEqual(result, expect)


class Test_binomial_cdf_inv(unittest.TestCase):

    def setUp(self):
        # x, a, b
        self.n1 = (0.1, 1000, 0.01)
        self.n2 = (0.01, 1000, 0.01)

    def test_binomial_cdf_inv_n1(self):
        expect = 6
        result = binomial_cdf_inv(self.n1[0], self.n1[1], self.n1[2])
        self.assertEqual(result, expect)

    def test_poisson_cdf_inv_n2(self):
        expect = 3
        result = binomial_cdf_inv(self.n2[0], self.n2[1], self.n2[2])
        self.assertEqual(result, expect)


class Test_binomial_pdf(unittest.TestCase):
    def setUp(self):
        # x, a, b
        self.n1 = (20, 1000, 0.01)
        self.n2 = (200, 1000, 0.01)

    def test_binomial_cdf_inv_n1(self):
        expect = round(0.001791878, 5)
        result = round(binomial_pdf(self.n1[0], self.n1[1], self.n1[2]), 5)
        self.assertEqual(result, expect)

    def test_poisson_cdf_inv_n2(self):
        expect = round(log10(2.132196e-188), 4)
        result = binomial_pdf(self.n2[0], self.n2[1], self.n2[2])
        result = round(log10(result), 4)
        self.assertEqual(result, expect)


# ------------------------------------
# Comprehensive tests (added below the original ones)
# ------------------------------------
# Expected values come from scipy.stats (norm, chi2, poisson, binom),
# math and scipy.special.logsumexp. The C-only helpers of Prob are not
# callable from Python; they are covered through these callers:
#   ex20                              -> chisq_pvalue_e
#   logspace_add                      -> chisq_logp_e,
#                                        poisson_cdf(..., log10=True)
#   log10_poisson_cdf_P_large_lambda  -> poisson_cdf(n, lam, True, True)
#   log10_poisson_cdf_Q_large_lambda  -> poisson_cdf(n, lam, False, True)
#   poz, binomial_coef                -> no caller anywhere in MACS3
# The private inline helpers __poisson_cdf, __poisson_cdf_large_lambda,
# __poisson_cdf_Q and __poisson_cdf_Q_large_lambda are reached through
# poisson_cdf(..., log10=False) with lam <= 700 and lam > 700, and
# _binomial_cdf_f/_binomial_cdf_r through binomial_cdf, binomial_sf and
# pduplication.

LN10 = math.log(10)


# ------------------------------------
# pnorm
# ------------------------------------

PNORM_CASES = [
    # x, u, v (v is the variance)
    (0, 0, 1),
    (1, 0, 1),
    (-3, 2, 4),
    (10, 10, 100),
    (5, 3, 2),
    (-100, 100, 10000),
    (30, 0, 1),                 # ~5.8e-197
    (40, 0, 1),                 # exp(-800) underflows to 0
    (0, 0, 2**31 - 1),          # int32 maximum as the variance
]


@pytest.mark.parametrize("x, u, v", PNORM_CASES)
def test_pnorm_is_the_normal_density(x, u, v):
    expected = stats.norm.pdf(x, loc=u, scale=math.sqrt(v))
    result = pnorm(x, u, v)
    assert isinstance(result, float)
    # the variance and x-u pass through float32, hence rel=1e-6
    assert result == pytest.approx(expected, rel=1e-6, abs=1e-300)


def test_pnorm_zero_variance_raises_zerodivision():
    with pytest.raises(ZeroDivisionError, match="float division"):
        pnorm(0, 0, 0)


def test_pnorm_negative_variance_raises_valueerror():
    # the math module's message is worded differently from Python 3.14
    with pytest.raises(ValueError, match="math domain error|expected a (positive|nonnegative) input"):
        pnorm(0, 0, -1)


@pytest.mark.parametrize("args", [(2**31, 0, 1), (0, -2**31 - 1, 1),
                                  (0, 0, 2**31)])
def test_pnorm_rejects_values_outside_int32(args):
    with pytest.raises(OverflowError, match="convert"):
        pnorm(*args)


# ------------------------------------
# pnorm2
# ------------------------------------

PNORM2_CASES = [
    # all values are exactly representable in float32
    (0.0, 0.0, 1.0),
    (1.5, 0.5, 2.0),
    (-2.0, 1.0, 0.25),
    (200.0, 180.0, 400.0),
    (50.0, 200.0, 400.0),
    (0.0, 0.0, 0.0625),
    (100.0, 0.0, 1.0),          # underflows to 0
]


@pytest.mark.parametrize("x, u, v", PNORM2_CASES)
def test_pnorm2_is_the_normal_density(x, u, v):
    expected = stats.norm.pdf(x, loc=u, scale=math.sqrt(v))
    # float32 return value
    assert pnorm2(x, u, v) == pytest.approx(expected, rel=1e-6, abs=1e-300)


def test_pnorm2_returns_float32_rounded_value():
    expected = float(np.float32(stats.norm.pdf(1.0, loc=0.0, scale=1.0)))
    assert pnorm2(1.0, 0.0, 1.0) == expected


def test_pnorm2_negative_variance_exits_with_code_1():
    # math.sqrt raises ValueError, which pnorm2 turns into sys.exit(1)
    with pytest.raises(SystemExit) as excinfo:
        pnorm2(0.0, 0.0, -1.0)
    assert excinfo.value.code == 1


def test_pnorm2_zero_variance_raises_zerodivision():
    with pytest.raises(ZeroDivisionError, match="float division"):
        pnorm2(0.0, 0.0, 0.0)


# ------------------------------------
# factorial
# ------------------------------------

@pytest.mark.parametrize("n", [0, 1, 2, 3, 5, 10, 18, 20, 22])
def test_factorial_exact_while_representable(n):
    # n! for n <= 22 is exactly representable in float64 (its odd part
    # is < 2**53), and so is every partial product
    result = factorial(n)
    assert isinstance(result, float)
    assert result == float(math.factorial(n))


@pytest.mark.parametrize("n", [23, 50, 100, 150, 170])
def test_factorial_large_n_relative_accuracy(n):
    # n successive roundings: relative error <= n * 2**-53 < 1e-13
    assert factorial(n) == pytest.approx(float(math.factorial(n)), rel=1e-13)


@pytest.mark.parametrize("n", [171, 1000])
def test_factorial_overflows_to_inf(n):
    # 171! > DBL_MAX
    assert factorial(n) == math.inf


@pytest.mark.parametrize("n", [-1, 2**32])
def test_factorial_rejects_values_outside_uint32(n):
    with pytest.raises(OverflowError, match="convert"):
        factorial(n)


# ------------------------------------
# chisq_pvalue_e
# ------------------------------------

CHISQ_X = [0.5, 1.0, 5.0, 10.0, 39.0, 41.0, 100.0, 400.0]
CHISQ_DF = [2, 4, 10, 22, 100]


@pytest.mark.parametrize("df", CHISQ_DF)
@pytest.mark.parametrize("x", CHISQ_X)
def test_chisq_pvalue_e_matches_chi2_sf(x, df):
    # For even df = 2k, sf(x) = exp(-a) * sum_{j<k} a**j / j!, a = x/2.
    # ex20 sets each of the k terms to 0 when its exponent is < -20, so
    # the absolute error is at most k * exp(-20).
    expected = stats.chi2.sf(x, df)
    result = chisq_pvalue_e(x, df)
    assert result == pytest.approx(expected, rel=1e-12,
                                   abs=(df // 2) * math.exp(-20))


@pytest.mark.parametrize("x, df", [(0.5, 4), (10.0, 10), (39.0, 22),
                                   (5.0, 100)])
def test_chisq_pvalue_e_exact_when_x_at_most_40(x, df):
    # a = x/2 <= 20, so ex20 is a plain exp and no term is truncated
    assert chisq_pvalue_e(x, df) == pytest.approx(stats.chi2.sf(x, df),
                                                  rel=1e-12, abs=0)


def test_chisq_pvalue_e_df2_cutoff_at_x_40():
    # df=2: sf = exp(-x/2); ex20 returns exp(-20) at a = 20 and 0 above
    assert chisq_pvalue_e(40.0, 2) == math.exp(-20.0)
    assert chisq_pvalue_e(40.000001, 2) == 0.0


@pytest.mark.parametrize("x", [0.0, -1.0, -1e300])
@pytest.mark.parametrize("df", [2, 10])
def test_chisq_pvalue_e_nonpositive_x_is_1(x, df):
    assert chisq_pvalue_e(x, df) == 1.0


@pytest.mark.parametrize("x", [1.0, 7.0, 30.0])
def test_chisq_pvalue_e_odd_df_documented_as_unsupported(x):
    """The docstring states that df must be even and odd df gives a
    wrong result. The loop runs floor((df-1)/2) times, so odd df
    computes exactly what df+1 computes (df=1 computes df=2)."""
    assert chisq_pvalue_e(x, 1) == chisq_pvalue_e(x, 2)
    assert chisq_pvalue_e(x, 3) == chisq_pvalue_e(x, 4)
    assert chisq_pvalue_e(x, 9) == chisq_pvalue_e(x, 10)


def test_chisq_pvalue_e_negative_df_raises_overflow():
    with pytest.raises(OverflowError, match="convert"):
        chisq_pvalue_e(1.0, -2)


# ------------------------------------
# chisq_logp_e
# ------------------------------------

def neg_log_chisq_sf_even(x, df):
    """-ln P(chi2_df > x) for even df, from the Poisson form of the
    chi-square survival function, evaluated in log space."""
    a = 0.5 * x
    k = df // 2
    terms = [j * math.log(a) - math.lgamma(j + 1) for j in range(k)]
    return a - float(logsumexp(terms))


CHISQ_LOGP_X = [0.5, 2.0, 10.0, 39.9, 40.1, 100.0, 1000.0, 7765.0, 1e5]
CHISQ_LOGP_DF = [2, 4, 10, 22, 80, 200]


@pytest.mark.parametrize("df", CHISQ_LOGP_DF)
@pytest.mark.parametrize("x", CHISQ_LOGP_X)
def test_chisq_logp_e_matches_log_space_reference(x, df):
    expected = neg_log_chisq_sf_even(x, df)
    # results near 0 (sf ~ 1) carry an absolute rounding error ~1e-16
    assert chisq_logp_e(x, df) == pytest.approx(expected, rel=1e-10,
                                                abs=1e-12)
    assert chisq_logp_e(x, df, log10=True) == pytest.approx(
        expected / LN10, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("x, df", [(10.0, 4), (100.0, 8), (54.0, 6),
                                   (565.0, 10), (1000.0, 80)])
def test_chisq_logp_e_matches_scipy_logsf(x, df):
    assert chisq_logp_e(x, df) == pytest.approx(-stats.chi2.logsf(x, df),
                                                rel=1e-9)


@pytest.mark.parametrize("x", [0.5, 3.0, 41.0, 1e6])
def test_chisq_logp_e_df2_is_half_x(x):
    # df=2: -ln(exp(-x/2)) = x/2, computed without rounding
    assert chisq_logp_e(x, 2) == 0.5 * x


@pytest.mark.parametrize("x", [0.0, -3.0])
def test_chisq_logp_e_nonpositive_x_is_0(x):
    assert chisq_logp_e(x, 4) == 0.0
    assert chisq_logp_e(x, 4, True) == 0.0


def test_chisq_logp_e_huge_x_where_linear_space_underflows():
    # chi2.sf(7765, 12) underflows to 0, its log does not
    assert chisq_pvalue_e(7765.0, 12) == 0.0
    assert chisq_logp_e(7765.0, 12) == pytest.approx(
        neg_log_chisq_sf_even(7765.0, 12), rel=1e-12)


# ------------------------------------
# poisson_cdf: linear space
# ------------------------------------

def poisson_tail(n, lam, lower):
    """P(X <= n) when lower, else P(X > n), X ~ Poisson(lam)."""
    if lower:
        return stats.poisson.cdf(n, lam)
    return stats.poisson.sf(n, lam)


SMALL_LAMBDAS = [1e-10, 0.5, 3.7, 10.0, 100.0, 700.0]
SMALL_LAMBDA_N = [0, 1, 5, 50, 300, 1000]


@pytest.mark.parametrize("lower", [True, False], ids=["lower", "upper"])
@pytest.mark.parametrize("n", SMALL_LAMBDA_N)
@pytest.mark.parametrize("lam", SMALL_LAMBDAS)
def test_poisson_cdf_small_lambda_matches_scipy(lam, n, lower):
    # lam <= 700: direct summation of the pmf recurrence. Values below
    # ~1e-290 (gradual underflow) are compared on an absolute scale.
    expected = poisson_tail(n, lam, lower)
    assert poisson_cdf(n, lam, lower) == pytest.approx(expected, rel=1e-9,
                                                       abs=1e-290)


LARGE_LAMBDAS = [700.5, 1000.0, 5000.0, 1e5]
LARGE_LAMBDA_Z = [-5, -1, 0, 1, 5]


@pytest.mark.parametrize("lower", [True, False], ids=["lower", "upper"])
@pytest.mark.parametrize("z", LARGE_LAMBDA_Z)
@pytest.mark.parametrize("lam", LARGE_LAMBDAS)
def test_poisson_cdf_large_lambda_matches_scipy(lam, z, lower):
    # lam > 700: the rescaled summation branches; n = lam + z*sqrt(lam)
    n = int(lam + z * math.sqrt(lam))
    expected = poisson_tail(n, lam, lower)
    assert poisson_cdf(n, lam, lower) == pytest.approx(expected, rel=1e-9,
                                                       abs=1e-290)


def test_poisson_cdf_default_is_upper_tail_linear():
    assert poisson_cdf(80, 100.0) == poisson_cdf(80, 100.0, False, False)
    assert poisson_cdf(80, 100.0) == pytest.approx(
        stats.poisson.sf(80, 100.0), rel=1e-12)


@pytest.mark.parametrize("n, lam", [(5000, 100.0), (2000, 3.0),
                                    (100000, 10000.0)])
def test_poisson_cdf_upper_tail_underflows_to_zero(n, lam):
    # the true tail is far below the smallest double (log10 < -350)
    assert log10_poisson_sf_ref(n, lam) < -350
    assert poisson_cdf(n, lam, False) == 0.0


def test_poisson_cdf_tiny_lambda_upper_tail():
    # P(X > 0) = 1 - exp(-lam) = lam for lam = 1e-300
    assert poisson_cdf(0, 1e-300, False) == pytest.approx(1e-300, rel=1e-12,
                                                          abs=0)
    assert poisson_cdf(0, 1e-300, True) == 1.0


# ------------------------------------
# poisson_cdf: log10 space
# ------------------------------------

def log10_poisson_sf_ref(n, lam):
    """log10 P(X > n) by logsumexp over the log pmf."""
    top = max(n, lam)
    hi = int(top + 60 * math.sqrt(top) + 200)
    i = np.arange(n + 1, hi)
    return float(logsumexp(stats.poisson.logpmf(i, lam))) / LN10


LOG10_UPPER_CASES = [
    # n, lam
    (0, 0.5),
    (5, 1.0),
    (10, 3.0),
    (20, 5.0),
    (50, 10.0),
    (100, 30.0),
    (150, 100.0),
    (200, 100.0),
    (1000, 100.0),
    (100, 1000.0),
    (1000, 1000.0),
    (1100, 1000.0),
    (1500, 1000.0),
    (3000, 1000.0),
    (103000, 1e5),
    (200000, 1e5),
]


@pytest.mark.parametrize("n, lam", LOG10_UPPER_CASES)
def test_poisson_cdf_log10_upper_matches_reference(n, lam):
    """The series stops once a term changes the natural-log sum by less
    than 1e-5 and the result is rounded to 5 decimals. The neglected
    tail is about term * lam/(m - lam), which stays below 1e-4 in log10
    for lam <= 1000 and below 1e-3 for lam = 1e5 on these cases."""
    tol = 1e-4 if lam <= 1000 else 1e-3
    expected = log10_poisson_sf_ref(n, lam)
    assert poisson_cdf(n, lam, False, True) == pytest.approx(expected,
                                                             abs=tol)


@pytest.mark.parametrize("n, lam", [(150, 100.0), (1500, 1000.0),
                                    (5000, 100.0), (7, 0.25)])
def test_poisson_cdf_log10_is_rounded_to_5_decimals(n, lam):
    for lower in (False, True):
        result = poisson_cdf(n, lam, lower, True)
        assert result == round(result, 5)


def test_poisson_cdf_log10_upper_finite_where_linear_underflows():
    assert poisson_cdf(5000, 100.0, False) == 0.0
    expected = log10_poisson_sf_ref(5000, 100.0)
    assert expected < -300
    assert poisson_cdf(5000, 100.0, False, True) == pytest.approx(expected,
                                                                  abs=1e-4)


@pytest.mark.parametrize("lam", [0.5, 10.0, 1000.0, 1e5])
def test_poisson_cdf_log10_lower_n0_is_minus_lambda_over_ln10(lam):
    # P(X <= 0) = exp(-lam), so log10 = -lam / ln(10), rounded to 5 places
    assert poisson_cdf(0, lam, True, True) == round(-lam / LN10, 5)


# ------------------------------------
# poisson_cdf: errors
# ------------------------------------

@pytest.mark.parametrize("log10", [False, True])
def test_poisson_cdf_lambda_zero_raises_assertion(log10):
    with pytest.raises(AssertionError,
                       match=r"^Lambda must > 0, however we got 0$"):
        poisson_cdf(1, 0.0, False, log10)


@pytest.mark.parametrize("lam", [-1.0, -2.5, -1e10])
def test_poisson_cdf_negative_lambda_raises_assertion(lam):
    with pytest.raises(AssertionError, match=r"^Lambda must > 0"):
        poisson_cdf(1, lam)


@pytest.mark.parametrize("n", [-1, 2**32])
def test_poisson_cdf_n_outside_uint32_raises_overflow(n):
    with pytest.raises(OverflowError, match="convert"):
        poisson_cdf(n, 1.0)


# ------------------------------------
# poisson_cdf_inv
# ------------------------------------

POISSON_INV_CASES = [
    # q, lam with q > P(X = 0), where the result equals scipy's ppf
    (0.5, 1.0),
    (0.9, 1.0),
    (0.999, 1.0),
    (0.01, 5.0),
    (0.5, 5.0),
    (0.99, 5.0),
    (0.1, 20.0),
    (0.95, 20.0),
    (0.5, 100.0),
    (0.001, 500.0),
    (0.999, 700.0),
]


@pytest.mark.parametrize("q, lam", POISSON_INV_CASES)
def test_poisson_cdf_inv_matches_ppf(q, lam):
    # returns the i >= 1 with CDF(i-1) <= q <= CDF(i)
    assert poisson_cdf_inv(q, lam) == int(stats.poisson.ppf(q, lam))


def test_poisson_cdf_inv_zero_is_zero():
    assert poisson_cdf_inv(0.0, 10.0) == 0


def test_poisson_cdf_inv_caps_at_maximum():
    # ppf(0.999, 500) is ~570; no i <= 100 qualifies
    assert stats.poisson.ppf(0.999, 500.0) > 100
    assert poisson_cdf_inv(0.999, 500.0, 100) == 100
    assert poisson_cdf_inv(0.999, 500.0, maximum=600) == int(
        stats.poisson.ppf(0.999, 500.0))


@pytest.mark.parametrize("q", [-0.1, 1.5, -1e-300])
def test_poisson_cdf_inv_out_of_range_cdf_raises(q):
    with pytest.raises(Exception, match=r"^CDF must >= 0 and <= 1$") as e:
        poisson_cdf_inv(q, 10.0)
    assert type(e.value) is Exception


@pytest.mark.parametrize("lam", [740.0, 1000.0])
def test_poisson_cdf_inv_lambda_740_or_more_raises_assertion(lam):
    with pytest.raises(AssertionError):
        poisson_cdf_inv(0.5, lam)


# ------------------------------------
# poisson_cdf_Q_inv
# ------------------------------------

def test_poisson_cdf_Q_inv_zero_is_zero():
    assert poisson_cdf_Q_inv(0.0, 10.0) == 0


@pytest.mark.parametrize("q", [-0.5, 2.0])
def test_poisson_cdf_Q_inv_out_of_range_cdf_raises(q):
    with pytest.raises(Exception, match=r"^CDF must >= 0 and <= 1$") as e:
        poisson_cdf_Q_inv(q, 10.0)
    assert type(e.value) is Exception


def test_poisson_cdf_Q_inv_lambda_740_raises_assertion():
    with pytest.raises(AssertionError):
        poisson_cdf_Q_inv(0.5, 740.0)


@pytest.mark.parametrize("q, lam", [(0.01, 10.0), (0.2, 50.0)])
def test_poisson_cdf_Q_inv_inverts_lower_tail_like_poisson_cdf_inv(q, lam):
    # Despite the Q in its name (Q marks the upper tail elsewhere in
    # Prob.py), poisson_cdf_Q_inv has the same body and docstring
    # ("cdf : the CDF") as poisson_cdf_inv: it returns the i >= 1 with
    # CDF(i-1) <= q <= CDF(i) of the lower tail, i.e. scipy's ppf.
    assert poisson_cdf_Q_inv(q, lam) == int(stats.poisson.ppf(q, lam))
    assert poisson_cdf_Q_inv(q, lam) == poisson_cdf_inv(q, lam)


# ------------------------------------
# poisson_pdf
# ------------------------------------

POISSON_PDF_CASES = [
    # k, lam
    (0, 0.5),
    (0, 1e-300),
    (1, 1.0),
    (3, 2.5),
    (10, 10.0),
    (50, 20.0),
    (100, 100.0),
    (150, 100.0),
    (170, 50.0),
    (0, 700.0),
]


@pytest.mark.parametrize("k, lam", POISSON_PDF_CASES)
def test_poisson_pdf_matches_scipy(k, lam):
    result = poisson_pdf(k, lam)
    assert isinstance(result, float)
    assert result == pytest.approx(stats.poisson.pmf(k, lam), rel=1e-12, abs=0)


@pytest.mark.parametrize("k, lam", [(0, -1.0), (3, -0.5), (3, 0.0)])
def test_poisson_pdf_nonpositive_lambda_returns_zero(k, lam):
    # documented by the code: a <= 0 returns 0 (correct for k > 0, lam = 0)
    assert poisson_pdf(k, lam) == 0.0


def test_poisson_pdf_lambda_zero_k_zero_returns_zero():
    # lambda = 0 is outside the Poisson parameter domain (poisson_cdf in
    # the same module asserts lam > 0). The explicit `a <= 0` guard
    # returns 0 for every k, including k = 0, where the degenerate limit
    # (scipy's pmf(0, 0)) would be 1.
    assert stats.poisson.pmf(0, 0.0) == 1.0
    assert poisson_pdf(0, 0.0) == 0.0


@pytest.mark.parametrize("k", [-1, 2**32])
def test_poisson_pdf_k_outside_uint32_raises_overflow(k):
    with pytest.raises(OverflowError, match="convert"):
        poisson_pdf(k, 1.0)


# ------------------------------------
# binomial_pdf
# ------------------------------------

BINOMIAL_PDF_CASES = [
    # x, n, p
    (0, 1, 0.5),
    (1, 1, 0.5),
    (3, 10, 0.3),
    (7, 10, 0.3),               # x > n - x branch
    (5, 10, 0.5),
    (0, 100, 0.01),
    (100, 100, 0.99),
    (20, 1000, 0.01),
    (200, 1000, 0.01),          # rescaling below 1e-100
    (500, 1000, 0.5),
    (5000, 10000, 0.5),
    (1, 1000000, 1e-6),
    (999, 1000, 0.999),
    (900, 1000, 0.01),          # underflows to 0
]


@pytest.mark.parametrize("x, n, p", BINOMIAL_PDF_CASES)
def test_binomial_pdf_matches_scipy(x, n, p):
    expected = stats.binom.pmf(x, n, p)
    assert binomial_pdf(x, n, p) == pytest.approx(expected, rel=1e-9,
                                                  abs=1e-300)


@pytest.mark.parametrize("x, n, p, expected", [
    (-1, 10, 0.5, 0.0),         # x < 0
    (11, 10, 0.5, 0.0),         # x > n
    (0, 10, 0.0, 1.0),          # p = 0: all mass at 0
    (3, 10, 0.0, 0.0),
    (10, 10, 1.0, 1.0),         # p = 1: all mass at n
    (9, 10, 1.0, 0.0),
    (0, -5, 0.5, 0.0),          # negative n
])
def test_binomial_pdf_boundaries(x, n, p, expected):
    assert binomial_pdf(x, n, p) == expected


# ------------------------------------
# binomial_cdf
# ------------------------------------

BINOMIAL_CDF_CASES = [
    # x, n, p
    (0, 10, 0.3),
    (3, 10, 0.3),
    (9, 10, 0.3),
    (10, 10, 0.3),
    (5, 1000, 0.01),
    (20, 1000, 0.01),
    (50, 100, 0.5),
    (500, 1000, 0.5),
    (4999, 10000, 0.5),
    (0, 100, 0.99),
    (99, 100, 0.99),
    (200, 1000, 0.01),
]


@pytest.mark.parametrize("lower", [True, False], ids=["lower", "upper"])
@pytest.mark.parametrize("x, n, p", BINOMIAL_CDF_CASES)
def test_binomial_cdf_matches_scipy(x, n, p, lower):
    # lower: P(X <= x); upper: P(X > x)
    if lower:
        expected = stats.binom.cdf(x, n, p)
    else:
        expected = stats.binom.sf(x, n, p)
    assert binomial_cdf(x, n, p, lower) == pytest.approx(expected, rel=1e-9,
                                                         abs=1e-300)


def test_binomial_cdf_default_is_lower_tail():
    assert binomial_cdf(3, 10, 0.3) == binomial_cdf(3, 10, 0.3, True)


@pytest.mark.parametrize("x, n, p, lower_expected", [
    (-1, 10, 0.5, 0.0),         # x < 0
    (11, 10, 0.5, 1.0),         # x > n
    (3, 10, 0.0, 1.0),          # p = 0
    (3, 10, 1.0, 0.0),          # p = 1, x < n
])
def test_binomial_cdf_boundaries(x, n, p, lower_expected):
    assert binomial_cdf(x, n, p, True) == lower_expected
    assert binomial_cdf(x, n, p, False) == 1.0 - lower_expected


# ------------------------------------
# binomial_sf
# ------------------------------------

@pytest.mark.parametrize("x, n, p", BINOMIAL_CDF_CASES)
def test_binomial_sf_lower_is_one_minus_cdf(x, n, p):
    # lower=True: 1 - P(X <= x) = P(X > x); computed by subtraction,
    # so the comparison is on an absolute scale
    assert binomial_sf(x, n, p) == pytest.approx(stats.binom.sf(x, n, p),
                                                 abs=1e-12)
    assert binomial_sf(x, n, p, True) == 1.0 - binomial_cdf(x, n, p, True)


@pytest.mark.parametrize("x, n, p", BINOMIAL_CDF_CASES)
def test_binomial_sf_upper_flag_returns_lower_cdf(x, n, p):
    # lower=False: 1 - P(X > x) = P(X <= x)
    assert binomial_sf(x, n, p, False) == pytest.approx(
        stats.binom.cdf(x, n, p), abs=1e-12)
    assert binomial_sf(x, n, p, False) == 1.0 - binomial_cdf(x, n, p, False)


# ------------------------------------
# pduplication
# ------------------------------------

def pduplication_ref(pmf, n_obs):
    """mean over the pmf of P(Binomial(n_obs, p) > 2), p rounded to
    float32 as pduplication does."""
    ps = [float(np.float32(p)) for p in pmf]
    return float(np.mean([stats.binom.sf(2, n_obs, p) for p in ps]))


@pytest.mark.parametrize("pmf, n_obs", [
    ([0.1, 0.2, 0.3, 0.4], 10),
    ([0.25, 0.25, 0.25, 0.25], 4),
    ([0.5], 3),
    ([0.001, 0.002, 0.997], 1000),
    ([1e-4] * 5, 100),
])
def test_pduplication_is_mean_binomial_sf_at_2(pmf, n_obs):
    result = pduplication(np.array(pmf, dtype="f8"), n_obs)
    # float32 accumulation and return value
    assert result == pytest.approx(pduplication_ref(pmf, n_obs), rel=1e-5,
                                   abs=1e-7)


@pytest.mark.parametrize("n_obs", [0, 1, 2])
def test_pduplication_fewer_than_3_fragments_is_zero(n_obs):
    # P(X > 2) = 0 when there are at most 2 trials
    result = pduplication(np.array([0.3, 0.7]), n_obs)
    assert result == pytest.approx(0.0, abs=1e-6)


def test_pduplication_accepts_float32_and_int_arrays():
    a = pduplication(np.array([0.5, 0.5], dtype="f4"), 10)
    b = pduplication(np.array([0.5, 0.5], dtype="f8"), 10)
    assert a == b
    assert pduplication(np.array([0, 0], dtype="i8"), 10) == 0.0


def test_pduplication_empty_pmf_raises_zerodivision():
    with pytest.raises(ZeroDivisionError, match="float division"):
        pduplication(np.array([], dtype="f8"), 10)


def test_pduplication_rejects_list():
    with pytest.raises(TypeError, match="incorrect type"):
        pduplication([0.5, 0.5], 10)


# ------------------------------------
# binomial_cdf_inv
# ------------------------------------

@pytest.mark.parametrize("q, n, p", [
    (0.1, 1000, 0.01),
    (0.01, 1000, 0.01),
    (0.5, 10, 0.3),
    (0.99, 100, 0.5),
    (0.01, 50, 0.9),
    (0.999, 20, 0.05),
    (0.3, 1, 0.5),
])
def test_binomial_cdf_inv_matches_ppf(q, n, p):
    # for q strictly between CDF values, smallest x with CDF(x) > q
    assert binomial_cdf_inv(q, n, p) == int(stats.binom.ppf(q, n, p))


def test_binomial_cdf_inv_uses_strict_inequality():
    # n=2, p=0.5: CDF = 0.25, 0.75, 1.0 (exact in binary). The function
    # returns the smallest x with q < CDF(x).
    assert binomial_cdf_inv(0.2, 2, 0.5) == 0
    assert binomial_cdf_inv(0.25, 2, 0.5) == 1
    assert binomial_cdf_inv(0.75, 2, 0.5) == 2


@pytest.mark.parametrize("q, n, p, expected", [
    (0.0, 10, 0.5, 0),          # 0 < CDF(0)
    (1.0, 10, 0.5, 10),         # sum reaches exactly 1.0, falls through
    (0.5, 0, 0.5, 0),           # n = 0
])
def test_binomial_cdf_inv_boundaries(q, n, p, expected):
    assert binomial_cdf_inv(q, n, p) == expected


@pytest.mark.parametrize("q", [-0.01, 1.01])
def test_binomial_cdf_inv_out_of_range_raises(q):
    with pytest.raises(Exception, match=r"^CDF must >= 0 or <= 1$") as e:
        binomial_cdf_inv(q, 10, 0.5)
    assert type(e.value) is Exception


@pytest.mark.parametrize("func, args", [
    (binomial_pdf, (0, 2**63, 0.5)),
    (binomial_pdf, (-2**63 - 1, 10, 0.5)),
    (binomial_cdf, (2**63, 10, 0.5)),
    (binomial_sf, (0, 2**63, 0.5)),
    (binomial_cdf_inv, (0.5, 2**63, 0.5)),
])
def test_binomial_functions_reject_values_outside_int64(func, args):
    with pytest.raises(OverflowError, match="convert"):
        func(*args)
