# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-04 15:09:10 Tao Liu>

"""Module Description: statistics functions to calculate p-values

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from math import fabs
from math import log1p          # as py_log1p
from math import sqrt
import sys
# ------------------------------------
# Other modules
# ------------------------------------
# import numpy as np
import cython
from cython.cimports import numpy as cnp
from cython.cimports.numpy import (uint32_t, uint64_t, int32_t,
                                   int64_t, float32_t, float64_t)
from cython.cimports.cpython import bool

# ------------------------------------
# C lib
# ------------------------------------
from cython.cimports.libc.math import exp, log, M_LN10  # ,fabs,log1p

# ------------------------------------
# constants
# ------------------------------------
LSTEP: int32_t = 200
EXPTHRES: float64_t = exp(LSTEP)
EXPSTEP: float64_t = exp(-1*LSTEP)
bigx: float64_t = 20

# ------------------------------------
# Normal distribution functions
# ------------------------------------
# x is the value, u is the mean, v is the variance


@cython.ccall
def pnorm(x: int32_t, u: int32_t, v: int32_t):
    """The probability of X=x when X=Norm(u,v)

    v is int32_t
    """
    i: float32_t = cython.cast(float32_t, v)
    m: float32_t = cython.cast(float32_t, (x-u))
    n: float32_t = cython.cast(float32_t, v)
    # 6.283185307179586 is 2 * pi
    return 1.0/sqrt(6.283185307179586 * i) * exp(-m**2 / (2.0 * n))


@cython.ccall
def pnorm2(x: float32_t, u: float32_t, v: float32_t) -> float32_t:
    """The probability of X=x when X=Norm(u,v)

    v is float32_t
    """
    ret: float32_t
    try:
        # 6.283185307179586 is 2 * pi
        ret = 1.0/sqrt(6.283185307179586 * v) * exp(-(x-u)**2 / (2.0 * v))
    except ValueError:
        sys.exit(1)
    return ret

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.ccall
def factorial(n: uint32_t) -> float64_t:
    """Calculate N!.

    """
    fact: float64_t = 1
    i: uint64_t
    if n < 0:
        return 0
    for i in range(2, n+1):
        fact = fact * i
    return fact


@cython.cfunc
def poz(z: float64_t) -> float64_t:
    """ probability of normal z value
    """
    y: float64_t
    x: float64_t
    w: float64_t
    Z_MAX: float64_t = 6.0      # Maximum meaningful z value

    if z == 0.0:
        x = 0.0
    else:
        y = 0.5 * fabs(z)
        if y >= (Z_MAX * 0.5):
            x = 1.0
        elif y < 1.0:
            w = y * y
            x = ((((((((0.000124818987 * w
                        - 0.001075204047) * w + 0.005198775019) * w
                      - 0.019198292004) * w + 0.059054035642) * w
                    - 0.151968751364) * w + 0.319152932694) * w
                  - 0.531923007300) * w + 0.797884560593) * y * 2.0
        else:
            y -= 2.0
            x = (((((((((((((-0.000045255659 * y
                             + 0.000152529290) * y - 0.000019538132) * y
                           - 0.000676904986) * y + 0.001390604284) * y
                         - 0.000794620820) * y - 0.002034254874) * y
                       + 0.006549791214) * y - 0.010557625006) * y
                     + 0.011630447319) * y - 0.009279453341) * y
                   + 0.005353579108) * y - 0.002141268741) * y
                 + 0.000535310849) * y + 0.999936657524
    if z > 0.0:
        return ((x + 1.0) * 0.5)
    else:
        return ((1.0 - x) * 0.5)


@cython.cfunc
def ex20(x: float64_t) -> float64_t:
    """Wrapper on exp function. It will return 0 if x is smaller than -20.
    """
    if x < -20.0:
        return 0.0
    else:
        return exp(x)


# ------------------------------------
# chi-square distribution functions
# ------------------------------------


@cython.ccall
def chisq_pvalue_e(x: float64_t, df: uint32_t) -> float64_t:
    """Chisq CDF function for upper tail and even degree of freedom.
    Good for p-value calculation and designed for combining pvalues.

    Note: we assume df to be an even number larger than 1! Do not
    violate this assumption and there is no checking.

    df has to be even number! if df is odd, result will be wrong!

    """
    a: float64_t
    y: float64_t
    s: float64_t
    e: float64_t
    c: float64_t
    z: float64_t

    if x <= 0.0:
        return 1.0

    a = 0.5 * x
    # even = ((2*(df/2)) == df)
    y = ex20(-a)
    s = y
    if df > 2:
        x = 0.5 * (df - 1.0)
        z = 1.0
        if a > bigx:            # approximation for big number
            e = 0.0
            c = log(a)
            while z <= x:
                e = log(z) + e
                s += ex20(c*z-a-e)
                z += 1.0
            return s
        else:
            e = 1.0
            c = 0.0
            while z <= x:
                e = e * (a / z)
                c = c + e
                z += 1.0
            return c * y + s
    else:
        return s


@cython.ccall
def chisq_logp_e(x: float64_t, df: uint32_t, log10: bool = False) -> float64_t:
    """Chisq CDF function in log space for upper tail and even degree
    of freedom Good for p-value calculation and designed for combining
    pvalues.

    Note: we assume df to be an even number larger than 1! Do not
    violate this assumption and there is no checking.

    Return value is -logp. If log10 is set as True, return -log10p
    instead.

    """
    a: float64_t
    y: float64_t
    s: float64_t
    e: float64_t
    c: float64_t
    z: float64_t

    if x <= 0.0:
        return 0.0

    a = 0.5 * x
    y = exp(-a)             # y is for small number calculation
    # initialize s
    s = -a
    if df > 2:
        x = 0.5 * (df - 1.0)    # control number of iterations
        z = 1.0             # variable for iteration
        if a > bigx:            # approximation for big number
            e = 0.0         # constant
            c = log(a)         # constant
            while z <= x:       # iterations
                e += log(z)
                s = logspace_add(s, c*z-a-e)
                z += 1.0
        else:                   # for small number
            e = 1.0             # not a constant
            c = 0.0             # not a constant
            while z <= x:
                e = e * (a / z)
                c = c + e
                z += 1.0
            s = log(y+c*y)      # logspace_add( s, log(c) ) - a
    # return
    if log10:
        return -s/log(10)
    else:
        return -s


# ------------------------------------
# Poisson distribution functions
# ------------------------------------ 


@cython.ccall
def poisson_cdf(n: uint32_t, lam: float64_t,
                lower: bool = False, log10: bool = False) -> float64_t:
    """Poisson CDF evaluater.

    This is a more stable CDF function. It can tolerate large lambda
    value. While the lambda is larger than 700, the function will be a
    little slower.

    Parameters:
    n     : your observation
    lam   : lambda of poisson distribution
    lower : if lower is False, calculate the upper tail CDF, otherwise, to
            calculate lower tail; Default is False.
    log10 : if log10 is True, calculation will be in log space. Default is
            False.
    """
    assert lam > 0.0, "Lambda must > 0, however we got %d" % lam

    if log10:
        if lower:
            # lower tail
            return log10_poisson_cdf_P_large_lambda(n, lam)
        else:
            # upper tail
            return log10_poisson_cdf_Q_large_lambda(n, lam)

    if lower:
        if lam > 700:           # may be problematic
            return __poisson_cdf_large_lambda(n, lam)
        else:
            return __poisson_cdf(n, lam)
    else:
        # upper tail
        if lam > 700: # may be problematic
            return __poisson_cdf_Q_large_lambda(n, lam)
        else:
            return __poisson_cdf_Q(n, lam)


@cython.cfunc
@cython.inline
def __poisson_cdf(k: uint32_t, a: float64_t) -> float64_t:
    """Poisson CDF For small lambda. If a > 745, this will return
    incorrect result.

    Parameters:
    k	: observation
    a	: lambda
    """
    nextcdf: float64_t
    cdf: float64_t
    i: uint32_t
    lastcdf: float64_t

    if k < 0:
        return 0.0                        # special cases

    nextcdf = exp(-1 * a)
    cdf = nextcdf

    for i in range(1, k + 1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf = cdf + nextcdf
    if cdf > 1.0:
        return 1.0
    else:
        return cdf


@cython.cfunc
@cython.inline
def __poisson_cdf_large_lambda(k: uint32_t, a: float64_t) -> float64_t:
    """Slower poisson cdf for large lambda ( > 700 )

    Parameters:
    k	: observation
    a	: lambda
    """
    num_parts: int32_t
    lastexp: float64_t
    nextcdf: float64_t
    cdf: float64_t = 0
    i: uint32_t
    lastcdf: float64_t

    assert a > 700
    if k < 0:
        return 0.0                        # special cases

    num_parts = cython.cast(int32_t, (a / LSTEP))
    lastexp = exp(-1 * (a % LSTEP))
    nextcdf = EXPSTEP

    num_parts -= 1

    for i in range(1, k + 1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        if nextcdf > EXPTHRES or cdf > EXPTHRES:
            if num_parts >= 1:
                cdf *= EXPSTEP
                nextcdf *= EXPSTEP
                num_parts -= 1
            else:
                cdf *= lastexp
                lastexp = 1

    for i in range(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf


@cython.cfunc
@cython.inline
def __poisson_cdf_Q(k: uint32_t, a: float64_t) -> float64_t:
    """internal Poisson CDF evaluater for upper tail with small
    lambda.

    Parameters:
    k	: observation
    a	: lambda
    """
    i: uint32_t
    nextcdf: float64_t
    lastcdf: float64_t
    cdf: float64_t = 0.0

    if k < 0:
        return 1.0                        # special cases
    nextcdf = exp(-1 * a)

    for i in range(1, k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i

    i = k+1
    while nextcdf > 0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i += 1
    return cdf


@cython.cfunc
@cython.inline
def __poisson_cdf_Q_large_lambda(k: uint32_t, a: float64_t) -> float64_t:
    """Slower internal Poisson CDF evaluater for upper tail with large
    lambda.

    Parameters:
    k	: observation
    a	: lambda
    """
    num_parts: uint32_t = cython.cast(int32_t, (a/LSTEP))
    lastexp: float64_t = exp(-1 * (a % LSTEP))
    nextcdf: float64_t = EXPSTEP
    i: uint32_t
    lastcdf: float64_t
    cdf: float64_t = 0.0

    assert a > 700
    if k < 0:
        return 1.0                        # special cases

    num_parts -= 1

    for i in range(1, k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        if nextcdf > EXPTHRES:
            if num_parts >= 1:
                nextcdf *= EXPSTEP
                num_parts -= 1
            else:
                # simply raise an error
                raise Exception("Unexpected error")
                # cdf *= lastexp
                # lastexp = 1

    i = k+1
    while nextcdf > 0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i += 1
        if nextcdf > EXPTHRES or cdf > EXPTHRES:
            if num_parts >= 1:
                cdf *= EXPSTEP
                nextcdf *= EXPSTEP
                num_parts -= 1
            else:
                cdf *= lastexp
                lastexp = 1

    for i in range(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf


@cython.cfunc
@cython.inline
def log10_poisson_cdf_P_large_lambda(k: uint32_t, lbd: float64_t) -> float64_t:
    """Slower Poisson CDF evaluater for lower tail which allow
    calculation in log space. Better for the pvalue < 10^-310.

    Parameters:
    k	: observation
    lbd	: lambda

    ret = -lambda + \\ln( \\sum_{i=k+1}^{\\inf} {lambda^i/i!} = -lambda + \\ln(sum{ exp{ln(F)} }), where F=lambda^m/m!
    \\ln{F(m)} = m*ln{lambda} - \\sum_{x=1}^{m}\\ln(x)
    Calculate \\ln( sum{exp{N} ) by logspace_add function

    Return the log10(pvalue)
    """
    residue: float64_t = 0
    logx: float64_t = 0
    ln_lbd: float64_t = log(lbd)
    m: int32_t = k
    sum_ln_m: float64_t = 0
    i: int32_t = 0

    for i in range(1, m+1):
        sum_ln_m += log(i)
    logx = m*ln_lbd - sum_ln_m
    residue = logx

    while m > 1:
        m -= 1
        logy = logx-ln_lbd+log(m)
        pre_residue = residue
        residue = logspace_add(pre_residue, logy)
        if fabs(pre_residue-residue) < 1e-10:
            break
        logx = logy

    return round((residue-lbd)/M_LN10, 5)


@cython.cfunc
@cython.inline
def log10_poisson_cdf_Q_large_lambda(k: uint32_t, lbd: float64_t) -> float64_t:
    """Slower Poisson CDF evaluater for upper tail which allow
    calculation in log space. Better for the pvalue < 10^-310.

    Parameters:
    k	: observation
    lbd	: lambda

    ret = -lambda + \\ln( \\sum_{i=k+1}^{\\inf} {lambda^i/i!} = -lambda + \\ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
    \\ln{F(m)} = m*ln{lambda} - \\sum_{x=1}^{m}\\ln(x)
    Calculate \\ln( sum{exp{N} ) by logspace_add function

    Return the log10(pvalue)
    """
    residue: float64_t = 0
    logx: float64_t = 0
    ln_lbd: float64_t = log(lbd)
    m: int32_t = k+1
    sum_ln_m: float64_t = 0
    i: int32_t = 0

    for i in range(1, m+1):
        sum_ln_m += log(i)
    logx = m*ln_lbd - sum_ln_m
    residue = logx

    while True:
        m += 1
        logy = logx+ln_lbd-log(m)
        pre_residue = residue
        residue = logspace_add(pre_residue, logy)
        if fabs(pre_residue-residue) < 1e-5:
            break
        logx = logy

    return round((residue-lbd)/log(10), 5)


@cython.cfunc
@cython.inline
def logspace_add(logx: float64_t, logy: float64_t) -> float64_t:
    """addition in log space.

    Given two log values, such as logx and logy, return
    log(exp(logx)+exp(logy)).

    """
    if logx > logy:
        return logx + log1p( exp ( logy - logx ) )
    else:
        return logy + log1p( exp ( logx - logy ) )


@cython.ccall
def poisson_cdf_inv ( cdf: float64_t, lam: float64_t, maximum: int32_t = 1000 ) -> int32_t:
    """inverse poisson distribution.

    cdf : the CDF
    lam : the lambda of poisson distribution

    note: maxmimum return value is 1000
    and lambda must be smaller than 740.
    """
    sum2: float64_t = 0
    newval: float64_t = exp( -1*lam )
    i: int32_t
    sumold: float64_t 
    lastval: float64_t 

    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception ("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0

    sum2 = newval

    for i in range(1,maximum+1):
        sumold = sum2
        lastval = newval
        newval = lastval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i

    return maximum


@cython.ccall
def poisson_cdf_Q_inv(cdf: float64_t, lam: float64_t, maximum: int32_t = 1000) -> int32_t:
    """inverse poisson distribution.

    cdf : the CDF
    lam : the lambda of poisson distribution

    note: maxmimum return value is 1000
    and lambda must be smaller than 740.
    """
    sum2: float64_t = 0
    newval: float64_t = exp(-1 * lam)
    i: int32_t
    sumold: float64_t
    lastval: float64_t

    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0
    sum2 = newval

    for i in range(1, maximum+1):
        sumold = sum2
        lastval = newval
        newval = lastval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i

    return maximum


@cython.ccall
def poisson_pdf(k: uint32_t, a: float64_t) -> float64_t:
    """Poisson PDF.

    PDF(K,A) is the probability that the number of events observed in
    a unit time period will be K, given the expected number of events
    in a unit time as A.
    """
    if a <= 0:
        return 0
    return exp(-a) * pow(a, k, None) / factorial(k)


# ------------------------------------
# binomial distribution functions
# ------------------------------------


@cython.cfunc
def binomial_coef(n: int64_t, k: int64_t) -> float64_t:
    """BINOMIAL_COEF computes the Binomial coefficient C(N,K)

    n,k are integers.
    """
    mn: int64_t = min(k, n-k)
    mx: int64_t
    cnk: float64_t
    i: int64_t
    if mn < 0:
        return 0
    elif mn == 0:
        return 1
    else:
        mx = max(k, n-k)
        cnk = cython.cast(float32_t, (mx+1))
        for i in range(2, mn+1):
            cnk = cnk * (mx+i) / i
    return cnk


@cython.ccall
def binomial_cdf(x: int64_t, a: int64_t, b: float64_t, lower: bool = True) -> float64_t:
    """ BINOMIAL_CDF compute the binomial CDF.

    CDF(x)(A,B) is the probability of at most X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return _binomial_cdf_f(x, a, b)
    else:
        return _binomial_cdf_r(x, a, b)


@cython.ccall
def binomial_sf(x: int64_t, a: int64_t, b: float64_t, lower: bool = True) -> float64_t:
    """ BINOMIAL_SF compute the binomial survival function (1-CDF)

    SF(x)(A,B) is the probability of more than X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return 1.0 - _binomial_cdf_f(x, a, b)
    else:
        return 1.0 - _binomial_cdf_r(x, a, b)


@cython.ccall
def pduplication(pmf: cnp.ndarray, N_obs: int32_t) -> float32_t:
    """return the probability of a duplicate fragment given a pmf
    and a number of observed fragments N_obs
    """
    n: int32_t = pmf.shape[0]
    p: float32_t
    sf: float32_t = 0.0

    for p in pmf:
        sf += binomial_sf(2, N_obs, p)
    return sf / cython.cast(float32_t, n)


@cython.cfunc
def _binomial_cdf_r(x: int64_t, a: int64_t, b: float64_t) -> float64_t:
    """ Binomial CDF for upper tail.

    """
    argmax: int64_t = cython.cast(int32_t, (a*b))
    seedpdf: float64_t
    cdf: float64_t
    pdf: float64_t
    i: int64_t

    if x < 0:
        return 1
    elif a < x:
        return 0
    elif b == 0:
        return 0
    elif b == 1:
        return 1

    if x < argmax:
        seedpdf = binomial_pdf(argmax, a, b)
        pdf = seedpdf
        cdf = pdf
        for i in range(argmax-1, x, -1):
            pdf /= (a-i)*b/(1-b)/(i+1)
            if pdf == 0.0:
                break
            cdf += pdf

        pdf = seedpdf
        i = argmax
        while True:
            pdf *= (a-i)*b/(1-b)/(i+1)
            if pdf == 0.0:
                break
            cdf += pdf
            i += 1
        cdf = min(1, cdf)
        # cdf = float("%.10e" %cdf)
        return cdf
    else:
        pdf = binomial_pdf(x+1, a, b)
        cdf = pdf
        i = x+1
        while True:
            pdf *= (a-i)*b/(1-b)/(i+1)
            if pdf == 0.0:
                break
            cdf += pdf
            i += 1
        cdf = min(1, cdf)
        # cdf = float("%.10e" %cdf)
        return cdf


@cython.cfunc
def _binomial_cdf_f(x: int64_t, a: int64_t, b: float64_t) -> float64_t:
    """ Binomial CDF for lower tail.

    """
    argmax: int64_t = cython.cast(int32_t, (a*b))
    seedpdf: float64_t
    cdf: float64_t
    pdf: float64_t
    i: int64_t

    if x < 0:
        return 0
    elif a < x:
        return 1
    elif b == 0:
        return 1
    elif b == 1:
        return 0

    if x > argmax:
        seedpdf = binomial_pdf(argmax, a, b)
        pdf = seedpdf
        cdf = pdf
        for i in range(argmax-1, -1, -1):
            pdf /= (a-i)*b/(1-b)/(i+1)
            if pdf == 0.0:
                break
            cdf += pdf

        pdf = seedpdf
        for i in range(argmax, x):
            pdf *= (a-i)*b/(1-b)/(i+1)
            if pdf == 0.0:
                break
            cdf += pdf
        cdf = min(1, cdf)
        # cdf = float("%.10e" %cdf)
        return cdf
    else:
        pdf = binomial_pdf(x, a, b)
        cdf = pdf
        for i in range(x-1, -1, -1):
            pdf /= (a-i)*b/(1-b)/(i+1)
            if pdf == 0.0:
                break
            cdf += pdf
        cdf = min(1, cdf)
        # cdf = float("%.10e" %cdf)
        return cdf


@cython.ccall
def binomial_cdf_inv(cdf: float64_t, a: int64_t, b: float64_t) -> int64_t:
    """BINOMIAL_CDF_INV inverts the binomial CDF. For lower tail only!

    """
    cdf2: float64_t = 0.0
    x: int64_t

    if cdf < 0 or cdf > 1:
        raise Exception("CDF must >= 0 or <= 1")

    for x in range(0, a+1):
        pdf = binomial_pdf(x, a, b)
        cdf2 = cdf2 + pdf
        if cdf < cdf2:
            return x
    return a


@cython.ccall
def binomial_pdf(x: int64_t, a: int64_t, b: float64_t) -> float64_t:
    """binomial PDF by H. Gene Shin

    """
    p: float64_t
    mn: int64_t
    mx: int64_t
    pdf: float64_t
    t: int64_t
    q: int64_t

    if a < 1:
        return 0.0
    elif x < 0 or a < x:
        return 0.0
    elif b == 0:
        if x == 0:
            return 1.0
        else:
            return 0.0
    elif b == 1:
        if x == a:
            return 1.0
        else:
            return 0.0

    if x > a-x:
        p = 1-b
        mn = a-x
        mx = x
    else:
        p = b
        mn = x
        mx = a-x
    pdf = 1
    t = 0
    for q in range(1, mn+1):
        pdf *= (a-q+1)*p/(mn-q+1)
        if pdf < 1e-100:
            while pdf < 1e-3:
                pdf /= 1-p
                t -= 1
        if pdf > 1e+100:
            while pdf > 1e+3 and t < mx:
                pdf *= 1-p
                t += 1

    for i in range(mx-t):
        pdf *= 1-p

    # pdf=float("%.10e" % pdf)
    return pdf

# cdef normal_01_cdf (float64_t x ):
#     """NORMAL_01_CDF evaluates the Normal 01 CDF.
#     """
#     cdef float64_t a1 = 0.398942280444
#     cdef float64_t a2 = 0.399903438504
#     cdef float64_t a3 = 5.75885480458
#     cdef float64_t a4 = 29.8213557808
#     cdef float64_t a5 = 2.62433121679
#     cdef float64_t a6 = 48.6959930692
#     cdef float64_t a7 = 5.92885724438
#     cdef float64_t b0 = 0.398942280385
#     cdef float64_t b1 = 3.8052E-08
#     cdef float64_t b2 = 1.00000615302
#     cdef float64_t b3 = 3.98064794E-04
#     cdef float64_t b4 = 1.98615381364
#     cdef float64_t b5 = 0.151679116635
#     cdef float64_t b6 = 5.29330324926
#     cdef float64_t b7 = 4.8385912808
#     cdef float64_t b8 = 15.1508972451
#     cdef float64_t b9 = 0.742380924027
#     cdef float64_t b10 = 30.789933034
#     cdef float64_t b11 = 3.99019417011
#     cdef float64_t cdf

#     if abs ( x ) <= 1.28:
#         y = 0.5 * x * x
#         q = 0.5 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) )
#     elif abs ( x ) <= 12.7:
#         y = 0.5 * x * x

#         q = exp ( - y ) * b0 / ( abs ( x ) - b1 \
#                                 + b2 / ( abs ( x ) + b3 \
#                                 + b4 / ( abs ( x ) - b5 \
#                                 + b6 / ( abs ( x ) + b7 \
#                                 - b8 / ( abs ( x ) + b9 \
#                                 + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
#     else :
#         q = 0.0
#     #
#     #  Take account of negative X.
#     #
#     if x < 0.0:
#         cdf = q
#     else:
#         cdf = 1.0 - q

#     return cdf

# def normal_cdf_inv(p, mu = None, sigma2 = None, lower=True):

#     upper = not lower
#     if p < 0 or p > 1:
#         raise Exception("Illegal argument %f for qnorm(p)." % p)

#     split = 0.42
#     a0 = 2.50662823884
#     a1 = -18.61500062529
#     a2 = 41.39119773534
#     a3 = -25.44106049637
#     b1 = -8.47351093090
#     b2 = 23.08336743743
#     b3 = -21.06224101826
#     b4 = 3.13082909833
#     c0 = -2.78718931138
#     c1 = -2.29796479134
#     c2 = 4.85014127135
#     c3 = 2.32121276858
#     d1 = 3.54388924762
#     d2 = 1.63706781897
#     q = p - 0.5

#     r = 0.0
#     ppnd = 0.0

#     if abs(q) <= split:
#         r = q * q
#         ppnd = q * (((a3 * r + a2) * r + a1) * r + a0) / ((((b4 * r + b3) * r + b2) * r + b1) * r + 1)
#     else:
#         r = p
#         if q > 0:
#             r = 1 - p

#         if r > 0:
#             r = math.sqrt(- math.log(r))
#             ppnd = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + 1)

#             if q < 0:
#                 ppnd = - ppnd
#         else:
#             ppnd = 0

#     if upper:
#         ppnd = - ppnd

#     if mu != None and sigma2 != None:
#         return ppnd * math.sqrt(sigma2) + mu
#     else:
#         return ppnd

# def normal_cdf (z, mu = 0.0, sigma2 = 1.0, lower=True):

#     upper = not lower

#     z = (z - mu) / math.sqrt(sigma2)

#     ltone = 7.0
#     utzero = 18.66
#     con = 1.28
#     a1 = 0.398942280444
#     a2 = 0.399903438504
#     a3 = 5.75885480458
#     a4 = 29.8213557808
#     a5 = 2.62433121679
#     a6 = 48.6959930692
#     a7 = 5.92885724438
#     b1 = 0.398942280385
#     b2 = 3.8052e-8
#     b3 = 1.00000615302
#     b4 = 3.98064794e-4
#     b5 = 1.986153813664
#     b6 = 0.151679116635
#     b7 = 5.29330324926
#     b8 = 4.8385912808
#     b9 = 15.1508972451
#     b10 = 0.742380924027
#     b11 = 30.789933034
#     b12 = 3.99019417011

#     y = 0.0
#     alnorm = 0.0

#     if z < 0:
#         upper = not upper
#         z = - z

#     if z <= ltone or upper and z <= utzero:
#         y = 0.5 * z * z
#         if z > con:
#             alnorm = b1 * math.exp(- y) / (z - b2 + b3 / (z + b4 + b5 / (z - b6 + b7 / (z + b8 - b9 / (z + b10 + b11 / (z + b12))))))
#         else:
#             alnorm = 0.5 - z * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))))
#     else:
#         alnorm = 0.0
#     if not upper:
#         alnorm = 1.0 - alnorm
#     return alnorm
