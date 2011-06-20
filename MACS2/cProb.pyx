# Time-stamp: <2011-06-20 17:34:31 Tao Liu>

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
from libc.math cimport exp,log,log10 #,fabs,log1p
from math import fabs
from math import log1p #as py_log1p

from cpython cimport bool
# ------------------------------------
# constants
# ------------------------------------
cdef int LSTEP = 200
cdef double EXPTHRES = exp(LSTEP)
cdef double EXPSTEP  = exp(-1*LSTEP)
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

cpdef poisson_cdf ( unsigned int n, double lam, bool lower=False, bool log10=False ):
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

cdef __poisson_cdf ( unsigned int k, double a ):
    """Poisson CDF For small lambda. If a > 745, this will return
    incorrect result.

    Parameters:
    k	: observation
    a	: lambda
    """
    if k < 0:
        return 0                        # special cases
    cdef double nextcdf = exp( -1 * a )
    cdef double cdf = nextcdf
    cdef unsigned int i
    cdef double lastcdf
    for i in xrange(1,k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf = cdf + nextcdf
    if cdf > 1:
        return 1
    else:
        return cdf
    
cdef __poisson_cdf_large_lambda ( unsigned int k, double a ):
    """Slower poisson cdf for large lambda ( > 700 )

    Parameters:
    k	: observation
    a	: lambda
    """
    assert a > 700
    if k < 0:
        return 0                        # special cases
    cdef int num_parts = int(a/LSTEP)
    cdef double lastexp = exp(-1 * (a % LSTEP) )
    cdef double nextcdf = EXPSTEP

    num_parts -= 1
    cdef double cdf = nextcdf
    cdef unsigned int i
    cdef double lastcdf
    for i in xrange(1,k+1):
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

    for i in xrange(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

cdef __poisson_cdf_Q ( unsigned int k, double a ):
    """internal Poisson CDF evaluater for upper tail with small
    lambda.

    Parameters:
    k	: observation
    a	: lambda
    """
    cdef unsigned int i

    if k < 0:
        return 1                        # special cases
    cdef double nextcdf
    nextcdf = exp( -1 * a )
    cdef double lastcdf

    for i in xrange(1,k+1):
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i

    cdef double cdf = 0
    i = k+1
    while nextcdf >0:
        lastcdf = nextcdf
        nextcdf = lastcdf * a / i
        cdf += nextcdf
        i+=1
    return cdf

cdef __poisson_cdf_Q_large_lambda ( unsigned int k, double a ):
    """Slower internal Poisson CDF evaluater for upper tail with large
    lambda.

    Parameters:
    k	: observation
    a	: lambda    
    """
    assert a > 700
    if k < 0:
        return 1                        # special cases
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
    cdef double cdf = 0
    i = k+1
    while nextcdf >0:
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

cdef double log10_poisson_cdf_P_large_lambda ( unsigned int k, double lbd ):
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

    return round((residue-lbd)/log(10),2)

cdef double log10_poisson_cdf_Q_large_lambda ( unsigned int k, double lbd ):
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

    return round((residue-lbd)/log(10),2)

cdef double logspace_add ( double logx, double logy ):
    return max (logx, logy) + log1p (exp (-fabs (logx - logy)));

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

