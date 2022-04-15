# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-04-14 17:47:31 Tao Liu>

"""Module Description: statistics functions to calculate p-values

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from math import fabs
from math import log1p #as py_log1p
from math import sqrt
import sys
# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
from cpython cimport bool

# ------------------------------------
# C lib
# ------------------------------------
from libc.math cimport exp,log,log10, M_LN10 #,fabs,log1p

# ------------------------------------
# constants
# ------------------------------------
cdef int32_t LSTEP = 200
cdef float64_t EXPTHRES = exp(LSTEP)
cdef float64_t EXPSTEP  = exp(-1*LSTEP)
cdef float64_t bigx = 20

# ------------------------------------
# Normal distribution functions
# ------------------------------------
# x is the value, u is the mean, v is the variance
cpdef pnorm(int32_t x, int32_t u, int32_t v):
    """The probability of X=x when X=Norm(u,v)
    """
    return 1.0/sqrt(2.0 * 3.141592653589793 * <float32_t>(v)) * exp(-<float32_t>(x-u)**2 / (2.0 * <float32_t>(v)))

cpdef float32_t pnorm2(float32_t x, float32_t u, float32_t v):
    """The probability of X=x when X=Norm(u,v)
    """
    cdef:
        float32_t ret

    try:
        ret =1.0/sqrt(2.0 * 3.141592653589793 * v) * exp(-(x-u)**2 / (2.0 * v))
    except ValueError:
        sys.exit(1)
    return ret

# ------------------------------------
# Misc functions
# ------------------------------------

cpdef factorial ( uint32_t n ):
    """Calculate N!.

    """
    cdef float64_t fact = 1
    cdef uint64_t i
    if n < 0:
        return 0
    for i in range( 2,n+1 ):
        fact = fact * i
    return fact

cdef float64_t poz(float64_t z):
    """ probability of normal z value
    """
    cdef:
        float64_t y, x, w
        float64_t Z_MAX = 6.0 # Maximum meaningful z value
    if z == 0.0:
        x = 0.0;
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

cdef float64_t ex20 ( float64_t x ):
    """Wrapper on exp function. It will return 0 if x is smaller than -20.
    """
    if x < -20.0:
        return 0.0
    else:
        return exp(x)

cpdef float64_t chisq_pvalue_e ( float64_t x, uint32_t df ):
    """Chisq CDF function for upper tail and even degree of freedom.
    Good for p-value calculation and designed for combining pvalues.

    Note: we assume df to be an even number larger than 1! Do not
    violate this assumption and there is no checking.

    df has to be even number! if df is odd, result will be wrong!

    """
    cdef:
        float64_t  a, y, s
        float64_t  e, c, z

    if x <= 0.0:
        return 1.0

    a = 0.5 * x
    even = ((2*(df/2)) == df)
    y = ex20(-a)
    s = y
    if df > 2:
        x = 0.5 * (df - 1.0)
        z = 1.0
        if a > bigx:            # approximation for big number
            e = 0.0
            c = log (a)
            while z <= x:
                e = log (z) + e
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

cpdef float64_t chisq_logp_e ( float64_t x, uint32_t df, bool log10 = False ):
    """Chisq CDF function in log space for upper tail and even degree of freedom
    Good for p-value calculation and designed for combining pvalues.

    Note: we assume df to be an even number larger than 1! Do not
    violate this assumption and there is no checking.

    Return value is -logp. If log10 is set as True, return -log10p
    instead.

    """
    cdef:
        float64_t a, y
        float64_t s                # s is the return value
        float64_t e, c, z

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
            c = log (a)         # constant
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
            s = log(y+c*y) #logspace_add( s, log(c) ) - a
    # return
    if log10:
        return -s/log(10)
    else:
        return -s

cpdef float64_t poisson_cdf ( uint32_t n, float64_t lam, bool lower=False, bool log10=False ):
    """Poisson CDF evaluater.

    This is a more stable CDF function. It can tolerate large lambda
    value. While the lambda is larger than 700, the function will be a
    little slower.

    Parameters:
    n     : your observation
    lam   : lambda of poisson distribution
    lower : if lower is False, calculate the upper tail CDF, otherwise, to calculate lower tail; Default is False.
    log10 : if log10 is True, calculation will be in log space. Default is False.
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
        if lam > 700: # may be problematic
            return __poisson_cdf_large_lambda (n, lam)
        else:
            return __poisson_cdf(n,lam)
    else:
        # upper tail
        if lam > 700: # may be problematic
            return __poisson_cdf_Q_large_lambda (n, lam)
        else:
            return __poisson_cdf_Q(n,lam)

cdef inline float64_t __poisson_cdf ( uint32_t k, float64_t a ):
    """Poisson CDF For small lambda. If a > 745, this will return
    incorrect result.

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef:
        float64_t nextcdf
        float64_t cdf
        uint32_t i
        float64_t lastcdf

    if k < 0:
        return 0.0                        # special cases

    nextcdf = exp( -1 * a )
    cdf = nextcdf

    for i in range( 1, k + 1 ):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf = cdf + nextcdf
    if cdf > 1.0:
        return 1.0
    else:
        return cdf

cdef inline float64_t __poisson_cdf_large_lambda ( uint32_t k, float64_t a ):
    """Slower poisson cdf for large lambda ( > 700 )

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef:
        int32_t num_parts
        float64_t lastexp
        float64_t nextcdf
        float64_t cdf
        uint32_t i
        float64_t lastcdf

    assert a > 700
    if k < 0:
        return 0.0                        # special cases

    num_parts = <int32_t>( a / LSTEP )
    lastexp = exp( -1 * ( a % LSTEP ) )
    nextcdf = EXPSTEP

    num_parts -= 1

    for i in range( 1 , k + 1 ):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf = cdf + nextcdf
        if nextcdf > EXPTHRES or cdf > EXPTHRES:
           if num_parts>=1:
               cdf *= EXPSTEP
               nextcdf *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1

    for i in range( num_parts ):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

cdef inline float64_t __poisson_cdf_Q ( uint32_t k, float64_t a ):
    """internal Poisson CDF evaluater for upper tail with small
    lambda.

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef uint32_t i

    if k < 0:
        return 1.0                        # special cases
    cdef float64_t nextcdf
    nextcdf = exp( -1 * a )
    cdef float64_t lastcdf

    for i in range(1,k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i

    cdef float64_t cdf = 0.0
    i = k+1
    while nextcdf >0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i+=1
    return cdf

cdef inline float64_t __poisson_cdf_Q_large_lambda ( uint32_t k, float64_t a ):
    """Slower internal Poisson CDF evaluater for upper tail with large
    lambda.

    Parameters:
    k	: observation
    a	: lambda
    """
    assert a > 700
    if k < 0: return 1.0                        # special cases
    cdef uint32_t num_parts = <int32_t>(a/LSTEP)
    cdef float64_t lastexp = exp(-1 * (a % LSTEP) )
    cdef float64_t nextcdf = EXPSTEP
    cdef uint32_t i
    cdef float64_t lastcdf

    num_parts -= 1

    for i in range(1,k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        if nextcdf > EXPTHRES:
           if num_parts>=1:
               nextcdf *= EXPSTEP
               num_parts -= 1
           else:
               # simply raise an error
               raise Exception("Unexpected error")
               #cdf *= lastexp
               #lastexp = 1
    cdef float64_t cdf = 0.0
    i = k+1
    while nextcdf >0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i+=1
        if nextcdf > EXPTHRES or cdf > EXPTHRES:
           if num_parts>=1:
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

cdef inline float64_t log10_poisson_cdf_P_large_lambda ( uint32_t k, float64_t lbd ):
    """Slower Poisson CDF evaluater for lower tail which allow
    calculation in log space. Better for the pvalue < 10^-310.

    Parameters:
    k	: observation
    lbd	: lambda

    ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
    \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
    Calculate \ln( sum{exp{N} ) by logspace_add function

    Return the log10(pvalue)
    """
    cdef float64_t residue = 0
    cdef float64_t logx = 0
    cdef float64_t ln_lbd = log(lbd)

    # first residue
    cdef int32_t m = k
    cdef float64_t sum_ln_m = 0
    cdef int32_t i = 0
    for i in range(1,m+1):
        sum_ln_m += log(i)
    logx = m*ln_lbd - sum_ln_m
    residue = logx

    while m > 1:
        m -= 1
        logy = logx-ln_lbd+log(m)
        pre_residue = residue
        residue = logspace_add(pre_residue,logy)
        if fabs(pre_residue-residue) < 1e-10:
            break
        logx = logy

    return round((residue-lbd)/M_LN10,5)

cdef inline float64_t log10_poisson_cdf_Q_large_lambda ( uint32_t k, float64_t lbd ):
    """Slower Poisson CDF evaluater for upper tail which allow
    calculation in log space. Better for the pvalue < 10^-310.

    Parameters:
    k	: observation
    lbd	: lambda

    ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
    \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
    Calculate \ln( sum{exp{N} ) by logspace_add function

    Return the log10(pvalue)
    """
    cdef float64_t residue = 0
    cdef float64_t logx = 0
    cdef float64_t ln_lbd = log(lbd)

    # first residue
    cdef int32_t m = k+1
    cdef float64_t sum_ln_m = 0
    cdef int32_t i = 0
    for i in range(1,m+1):
        sum_ln_m += log(i)
    logx = m*ln_lbd - sum_ln_m
    residue = logx

    while True:
        m += 1
        logy = logx+ln_lbd-log(m)
        pre_residue = residue
        residue = logspace_add(pre_residue,logy)
        if fabs(pre_residue-residue) < 1e-5:
            break
        logx = logy

    return round((residue-lbd)/log(10),5)

cdef inline float64_t logspace_add ( float64_t logx, float64_t logy ):
    """addition in log space.

    Given two log values, such as logx and logy, return
    log(exp(logx)+exp(logy)).

    """
    if logx > logy:
        return logx + log1p( exp ( logy - logx ) )
    else:
        return logy + log1p( exp ( logx - logy ) )

cpdef poisson_cdf_inv ( float64_t cdf, float64_t lam, int32_t maximum=1000 ):
    """inverse poisson distribution.

    cdf : the CDF
    lam : the lambda of poisson distribution

    note: maxmimum return value is 1000
    and lambda must be smaller than 740.
    """
    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception ("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0
    cdef float64_t sum2 = 0
    cdef float64_t newval = exp( -1*lam )
    sum2 = newval

    cdef int32_t i
    cdef float64_t sumold
    cdef float64_t lastval

    for i in range(1,maximum+1):
        sumold = sum2
        lastval = newval
        newval = lastval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i

    return maximum

cpdef poisson_cdf_Q_inv ( float64_t cdf, float64_t lam, int32_t maximum=1000 ):
    """inverse poisson distribution.

    cdf : the CDF
    lam : the lambda of poisson distribution

    note: maxmimum return value is 1000
    and lambda must be smaller than 740.
    """
    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception ("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0
    cdef float64_t sum2 = 0
    cdef float64_t newval = exp( -1 * lam )
    sum2 = newval

    cdef int32_t i
    cdef float64_t lastval
    cdef float64_t sumold

    for i in range(1,maximum+1):
        sumold = sum2
        lastval = newval
        newval = lastval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i

    return maximum

cpdef poisson_pdf ( uint32_t k, float64_t a ):
    """Poisson PDF.

    PDF(K,A) is the probability that the number of events observed in
    a unit time period will be K, given the expected number of events
    in a unit time as A.
    """
    if a <= 0:
        return 0
    return exp(-a) * pow (a, k, None) / factorial (k)


cdef binomial_coef ( int64_t n, int64_t k ):
    """BINOMIAL_COEF computes the Binomial coefficient C(N,K)

    n,k are integers.
    """
    cdef int64_t mn = min (k, n-k)
    cdef int64_t mx
    cdef float64_t cnk
    cdef int64_t i
    if mn < 0:
        return 0
    elif mn == 0:
        return 1
    else:
        mx = max(k,n-k)
        cnk = <float32_t>(mx+1)
        for i in range(2,mn+1):
            cnk = cnk * (mx+i) / i
    return cnk

cpdef binomial_cdf ( int64_t x, int64_t a, float64_t b, bool lower=True ):
    """ BINOMIAL_CDF compute the binomial CDF.

    CDF(x)(A,B) is the probability of at most X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return _binomial_cdf_f (x,a,b)
    else:
        return _binomial_cdf_r (x,a,b)

cpdef binomial_sf ( int64_t x, int64_t a, float64_t b, bool lower=True ):
    """ BINOMIAL_SF compute the binomial survival function (1-CDF)

    SF(x)(A,B) is the probability of more than X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return 1.0 - _binomial_cdf_f (x,a,b)
    else:
        return 1.0 - _binomial_cdf_r (x,a,b)

cpdef pduplication (np.ndarray[np.float64_t] pmf, int32_t N_obs):
    """return the probability of a duplicate fragment given a pmf
    and a number of observed fragments N_obs
    """
    cdef:
        n = pmf.shape[0]
        float32_t p, sf = 0.0
    for p in pmf:
        sf += binomial_sf(2, N_obs, p)
    return sf / <float32_t>n

cdef _binomial_cdf_r ( int64_t x, int64_t a, float64_t b ):
    """ Binomial CDF for upper tail.

    """
    cdef int64_t argmax=<int32_t>(a*b)
    cdef float64_t seedpdf
    cdef float64_t cdf
    cdef float64_t pdf
    cdef int64_t i

    if x < 0:
        return 1
    elif a < x:
        return 0
    elif b == 0:
        return 0
    elif b == 1:
        return 1

    if x<argmax:
        seedpdf=binomial_pdf(argmax,a,b)
        pdf=seedpdf
        cdf = pdf
        for i in range(argmax-1,x,-1):
            pdf/=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf += pdf

        pdf = seedpdf
        i = argmax
        while True:
            pdf*=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf+=pdf
            i+=1
        cdf = min(1,cdf)
        #cdf = float("%.10e" %cdf)
        return cdf
    else:
        pdf=binomial_pdf(x+1,a,b)
        cdf = pdf
        i = x+1
        while True:
            pdf*=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf += pdf
            i+=1
        cdf = min(1,cdf)
        #cdf = float("%.10e" %cdf)
        return cdf


cdef _binomial_cdf_f ( int64_t x, int64_t a, float64_t b ):
    """ Binomial CDF for lower tail.

    """
    cdef int64_t argmax=<int32_t>(a*b)
    cdef float64_t seedpdf
    cdef float64_t cdf
    cdef float64_t pdf
    cdef int64_t i

    if x < 0:
        return 0
    elif a < x:
        return 1
    elif b == 0:
        return 1
    elif b == 1:
        return 0

    if x>argmax:
        seedpdf=binomial_pdf(argmax,a,b)
        pdf=seedpdf
        cdf = pdf
        for i in range(argmax-1,-1,-1):
            pdf/=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf += pdf

        pdf = seedpdf
        for i in range(argmax,x):
            pdf*=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf+=pdf
        cdf=min(1,cdf)
        #cdf = float("%.10e" %cdf)
        return cdf
    else:
        pdf=binomial_pdf(x,a,b)
        cdf = pdf
        for i in range(x-1,-1,-1):
            pdf/=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf += pdf
        cdf=min(1,cdf)
        #cdf = float("%.10e" %cdf)
        return cdf

cpdef binomial_cdf_inv ( float64_t cdf, int64_t a, float64_t b ):
    """BINOMIAL_CDF_INV inverts the binomial CDF. For lower tail only!

    """
    if cdf < 0 or cdf >1:
        raise Exception("CDF must >= 0 or <= 1")

    cdef float64_t cdf2 = 0
    cdef int64_t x

    for x in range(0,a+1):
        pdf = binomial_pdf (x,a,b)
        cdf2 = cdf2 + pdf
        if cdf < cdf2:
            return x
    return a

cpdef binomial_pdf( int64_t x, int64_t a, float64_t b ):
    """binomial PDF by H. Gene Shin

    """

    if a<1:
        return 0
    elif x<0 or a<x:
        return 0
    elif b==0:
        if x==0:
            return 1
        else:
            return 0
    elif b==1:
        if x==a:
            return 1
        else:
            return 0

    cdef float64_t p
    cdef int64_t mn
    cdef int64_t mx
    cdef float64_t pdf
    cdef int64_t t
    cdef int64_t q

    if x>a-x:
        p=1-b
        mn=a-x
        mx=x
    else:
        p=b
        mn=x
        mx=a-x
    pdf=1
    t = 0
    for q in range(1,mn+1):
        pdf*=(a-q+1)*p/(mn-q+1)
        if pdf < 1e-100:
            while pdf < 1e-3:
                pdf /= 1-p
                t-=1
        if pdf > 1e+100:
            while pdf > 1e+3 and t<mx:
                pdf *= 1-p
                t+=1

    for i in range(mx-t):
        pdf *= 1-p

    #pdf=float("%.10e" % pdf)
    return pdf

# cdef normal_01_cdf ( float64_t x ):
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
