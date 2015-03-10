# Time-stamp: <2013-09-15 21:55:17 Tao Liu>

"""Module Description

Copyright (c) 2008,2009,2010,2011 Hyunjin Shin, Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Hyunjin Gene Shin, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
from libc.math cimport exp,log,log10, M_LN10 #,fabs,log1p
from math import fabs
from math import log1p #as py_log1p
from math import sqrt

import numpy as np
cimport numpy as np

from cpython cimport bool
# ------------------------------------
# constants
# ------------------------------------
cdef int LSTEP = 200
cdef double EXPTHRES = exp(LSTEP)
cdef double EXPSTEP  = exp(-1*LSTEP)

# ------------------------------------
# Normal distribution functions
# ------------------------------------
# x is the value, u is the mean, v is the variance
cpdef pnorm(int x, int u, int v):
    """The probability of X=x when X=Norm(u,v)
    """
    return 1.0/sqrt(2.0 * 3.141592653589793 * <float>v) * exp(-<float>(x-u)**2 / (2.0 * <float>v))

# ------------------------------------
# Misc functions
# ------------------------------------

cpdef factorial ( unsigned int n ):
    """Calculate N!.
    
    """
    cdef double fact = 1
    cdef unsigned long i
    if n < 0:
        return 0
    for i in xrange( 2,n+1 ):
        fact = fact * i
    return fact

cpdef double poisson_cdf ( unsigned int n, double lam, bool lower=False, bool log10=False ):
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
        if lam > 700:
            return __poisson_cdf_large_lambda (n, lam)
        else:
            return __poisson_cdf(n,lam)
    else:
        # upper tail
        if lam > 700:
            return __poisson_cdf_Q_large_lambda (n, lam)
        else:
            return __poisson_cdf_Q(n,lam)

cdef inline double __poisson_cdf ( unsigned int k, double a ):
    """Poisson CDF For small lambda. If a > 745, this will return
    incorrect result.

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef:
        double nextcdf
        double cdf
        unsigned int i
        double lastcdf

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
    
cdef inline double __poisson_cdf_large_lambda ( unsigned int k, double a ):
    """Slower poisson cdf for large lambda ( > 700 )

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef:
        int num_parts
        double lastexp
        double nextcdf
        double cdf
        unsigned int i
        double lastcdf
    
    assert a > 700
    if k < 0:
        return 0.0                        # special cases

    num_parts = int( a / LSTEP )
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

cdef inline double __poisson_cdf_Q ( unsigned int k, double a ):
    """internal Poisson CDF evaluater for upper tail with small
    lambda.

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef unsigned int i

    if k < 0:
        return 1.0                        # special cases
    cdef double nextcdf
    nextcdf = exp( -1 * a )
    cdef double lastcdf

    for i in xrange(1,k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i

    cdef double cdf = 0.0
    i = k+1
    while nextcdf >0.0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i+=1
    return cdf

cdef inline double __poisson_cdf_Q_large_lambda ( unsigned int k, double a ):
    """Slower internal Poisson CDF evaluater for upper tail with large
    lambda.

    Parameters:
    k	: observation
    a	: lambda    
    """
    assert a > 700
    if k < 0:
        return 1.0                        # special cases
    cdef unsigned int num_parts = int(a/LSTEP)
    cdef double lastexp = exp(-1 * (a % LSTEP) )
    cdef double nextcdf = EXPSTEP
    cdef unsigned int i
    cdef double lastcdf

    num_parts -= 1

    for i in xrange(1,k+1):
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
    cdef double cdf = 0.0
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

    for i in xrange(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

cdef inline double log10_poisson_cdf_P_large_lambda ( unsigned int k, double lbd ):
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
    cdef double residue = 0
    cdef double logx = 0
    cdef double ln_lbd = log(lbd)

    # first residue
    cdef int m = k
    cdef double sum_ln_m = 0
    cdef int i = 0
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

cdef inline double log10_poisson_cdf_Q_large_lambda ( unsigned int k, double lbd ):
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
    cdef double residue = 0
    cdef double logx = 0
    cdef double ln_lbd = log(lbd)

    # first residue
    cdef int m = k+1
    cdef double sum_ln_m = 0
    cdef int i = 0
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

cdef inline double logspace_add ( double logx, double logy ):
    return max (logx, logy) + log1p (exp (-fabs (logx - logy)))

cpdef poisson_cdf_inv ( double cdf, double lam, int maximum=1000 ):
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
    cdef double sum2 = 0
    cdef double newval = exp( -1*lam )
    sum2 = newval

    cdef int i
    cdef double sumold
    cdef double lastval

    for i in xrange(1,maximum+1):
        sumold = sum2
        lastval = newval
        newval = lastval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    
    return maximum

cpdef poisson_cdf_Q_inv ( double cdf, double lam, int maximum=1000 ):
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
    cdef double sum2 = 0
    cdef double newval = exp( -1 * lam )
    sum2 = newval

    cdef int i
    cdef double lastval
    cdef double sumold

    for i in xrange(1,maximum+1):
        sumold = sum2
        lastval = newval
        newval = lastval * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    
    return maximum

cpdef poisson_pdf ( unsigned int k, double a ):
    """Poisson PDF.

    PDF(K,A) is the probability that the number of events observed in
    a unit time period will be K, given the expected number of events
    in a unit time as A.
    """
    if a <= 0:
        return 0
    return exp(-a) * pow (a, k, None) / factorial (k)
    

cdef binomial_coef ( long n, long k ):
    """BINOMIAL_COEF computes the Binomial coefficient C(N,K)

    n,k are integers.
    """
    cdef long mn = min (k, n-k)
    cdef long mx
    cdef double cnk
    cdef long i
    if mn < 0:
        return 0
    elif mn == 0:
        return 1
    else:
        mx = max(k,n-k)
        cnk = float(mx+1)
        for i in xrange(2,mn+1):
            cnk = cnk * (mx+i) / i
    return cnk

cpdef binomial_cdf ( long x, long a, double b, bool lower=True ):
    """ BINOMIAL_CDF compute the binomial CDF.

    CDF(x)(A,B) is the probability of at most X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return _binomial_cdf_f (x,a,b)
    else:
        return _binomial_cdf_r (x,a,b)

cpdef binomial_sf ( long x, long a, double b, bool lower=True ):
    """ BINOMIAL_SF compute the binomial survival function (1-CDF)

    SF(x)(A,B) is the probability of more than X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return 1.0 - _binomial_cdf_f (x,a,b)
    else:
        return 1.0 - _binomial_cdf_r (x,a,b)

cpdef pduplication (np.ndarray[np.float64_t] pmf, int N_obs):
    """return the probability of a duplicate fragment given a pmf
    and a number of observed fragments N_obs
    """
    cdef:
        n = pmf.shape[0]
        float p, sf = 0.0
    for p in pmf:
        sf += binomial_sf(2, N_obs, p)
    return sf / <float>n

cdef _binomial_cdf_r ( long x, long a, double b ):
    """ Binomial CDF for upper tail.

    """
    if x < 0:
        return 1
    elif a < x:
        return 0
    elif b == 0:
        return 0
    elif b == 1:
        return 1

    cdef long argmax=int(a*b)
    cdef double seedpdf
    cdef double cdf
    cdef double pdf
    cdef long i
    
    if x<argmax:
        seedpdf=binomial_pdf(argmax,a,b)
        pdf=seedpdf
        cdf = pdf
        for i in xrange(argmax-1,x,-1):
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
        cdf=min(1,cdf)
        cdf = float("%.10e" %cdf)
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
        cdf=min(1,cdf)
        cdf = float("%.10e" %cdf)
        return cdf


cdef _binomial_cdf_f ( long x, long a, double b ):
    """ Binomial CDF for lower tail.

    """    
    if x < 0:
        return 0
    elif a < x:
        return 1
    elif b == 0:
        return 1
    elif b == 1:
        return 0

    cdef long argmax=int(a*b)
    cdef double seedpdf
    cdef double cdf
    cdef double pdf
    cdef long i

    if x>argmax:
        seedpdf=binomial_pdf(argmax,a,b)
        pdf=seedpdf
        cdf = pdf
        for i in xrange(argmax-1,-1,-1):
            pdf/=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf += pdf
            
        pdf = seedpdf
        for i in xrange(argmax,x):
            pdf*=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf+=pdf
        cdf=min(1,cdf)
        cdf = float("%.10e" %cdf)
        return cdf
    else:
        pdf=binomial_pdf(x,a,b)
        cdf = pdf
        for i in xrange(x-1,-1,-1):
            pdf/=(a-i)*b/(1-b)/(i+1)
            if pdf==0.0: break
            cdf += pdf
        cdf=min(1,cdf)
        cdf = float("%.10e" %cdf)
        return cdf

cpdef binomial_cdf_inv ( double cdf, long a, double b ):
    """BINOMIAL_CDF_INV inverts the binomial CDF. For lower tail only!

    """
    if cdf < 0 or cdf >1:
        raise Exception("CDF must >= 0 or <= 1")

    cdef double cdf2 = 0
    cdef long x

    for x in xrange(0,a+1):
        pdf = binomial_pdf (x,a,b)
        cdf2 = cdf2 + pdf
        if cdf < cdf2:
            return x
    return a
    
cpdef binomial_pdf( long x, long a, double b ):
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

    cdef double p
    cdef long mn
    cdef long mx
    cdef double pdf
    cdef long t
    cdef long q

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
    for q in xrange(1,mn+1):
        pdf*=(a-q+1)*p/(mn-q+1)
        if pdf < 1e-100:
            while pdf < 1e-3:
                pdf /= 1-p
                t-=1
        if pdf > 1e+100:
            while pdf > 1e+3 and t<mx:
                pdf *= 1-p
                t+=1

    for i in xrange(mx-t):
        pdf *= 1-p
        
    pdf=float("%.10e" % pdf)
    return pdf

# cdef normal_01_cdf ( double x ):
#     """NORMAL_01_CDF evaluates the Normal 01 CDF.
#     """
#     cdef double a1 = 0.398942280444
#     cdef double a2 = 0.399903438504
#     cdef double a3 = 5.75885480458
#     cdef double a4 = 29.8213557808
#     cdef double a5 = 2.62433121679
#     cdef double a6 = 48.6959930692
#     cdef double a7 = 5.92885724438
#     cdef double b0 = 0.398942280385
#     cdef double b1 = 3.8052E-08
#     cdef double b2 = 1.00000615302
#     cdef double b3 = 3.98064794E-04
#     cdef double b4 = 1.98615381364
#     cdef double b5 = 0.151679116635
#     cdef double b6 = 5.29330324926
#     cdef double b7 = 4.8385912808
#     cdef double b8 = 15.1508972451
#     cdef double b9 = 0.742380924027
#     cdef double b10 = 30.789933034
#     cdef double b11 = 3.99019417011
#     cdef double cdf

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
