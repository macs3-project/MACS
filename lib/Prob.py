# Time-stamp: <2011-02-14 15:49:14 Tao Liu>

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
from math import exp
from math import log
# ------------------------------------
# constants
# ------------------------------------
LSTEP = 200
EXPTHRES = exp(LSTEP)
EXPSTEP  = exp(-LSTEP)
# ------------------------------------
# Misc functions
# ------------------------------------

import math

def normal_01_cdf ( x ):
    """NORMAL_01_CDF evaluates the Normal 01 CDF.
    """
    a1 = 0.398942280444
    a2 = 0.399903438504
    a3 = 5.75885480458
    a4 = 29.8213557808
    a5 = 2.62433121679
    a6 = 48.6959930692
    a7 = 5.92885724438
    b0 = 0.398942280385
    b1 = 3.8052E-08
    b2 = 1.00000615302
    b3 = 3.98064794E-04
    b4 = 1.98615381364
    b5 = 0.151679116635
    b6 = 5.29330324926
    b7 = 4.8385912808
    b8 = 15.1508972451
    b9 = 0.742380924027
    b10 = 30.789933034
    b11 = 3.99019417011

    if abs ( x ) <= 1.28:
        y = 0.5 * x * x
        q = 0.5 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) )
    elif abs ( x ) <= 12.7:
        y = 0.5 * x * x

        q = exp ( - y ) * b0 / ( abs ( x ) - b1 \
                                + b2 / ( abs ( x ) + b3 \
                                + b4 / ( abs ( x ) - b5 \
                                + b6 / ( abs ( x ) + b7 \
                                - b8 / ( abs ( x ) + b9 \
                                + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
    else :
        q = 0.0
    #
    #  Take account of negative X.
    #
    if x < 0.0:
        cdf = q
    else:
        cdf = 1.0 - q

    return cdf

def normal_cdf_inv(p, mu = None, sigma2 = None, lower=True):

    upper = not lower
    if p < 0 or p > 1:
        raise Exception("Illegal argument %f for qnorm(p)." % p)
    
    split = 0.42
    a0 = 2.50662823884
    a1 = -18.61500062529
    a2 = 41.39119773534
    a3 = -25.44106049637
    b1 = -8.47351093090
    b2 = 23.08336743743
    b3 = -21.06224101826
    b4 = 3.13082909833
    c0 = -2.78718931138
    c1 = -2.29796479134
    c2 = 4.85014127135
    c3 = 2.32121276858
    d1 = 3.54388924762
    d2 = 1.63706781897
    q = p - 0.5
    
    r = 0.0
    ppnd = 0.0
    
    if abs(q) <= split:
        r = q * q
        ppnd = q * (((a3 * r + a2) * r + a1) * r + a0) / ((((b4 * r + b3) * r + b2) * r + b1) * r + 1)
    else:
        r = p
        if q > 0:
            r = 1 - p
        
        if r > 0:
            r = math.sqrt(- math.log(r))
            ppnd = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + 1)
            
            if q < 0:
                ppnd = - ppnd
        else:
            ppnd = 0
            
    if upper:
        ppnd = - ppnd
    
    if mu != None and sigma2 != None:
        return ppnd * math.sqrt(sigma2) + mu
    else:
        return ppnd

def normal_cdf (z, mu = 0.0, sigma2 = 1.0, lower=True):
    
    upper = not lower

    z = (z - mu) / math.sqrt(sigma2)
    
    ltone = 7.0
    utzero = 18.66
    con = 1.28
    a1 = 0.398942280444
    a2 = 0.399903438504
    a3 = 5.75885480458
    a4 = 29.8213557808
    a5 = 2.62433121679
    a6 = 48.6959930692
    a7 = 5.92885724438
    b1 = 0.398942280385
    b2 = 3.8052e-8
    b3 = 1.00000615302
    b4 = 3.98064794e-4
    b5 = 1.986153813664
    b6 = 0.151679116635
    b7 = 5.29330324926
    b8 = 4.8385912808
    b9 = 15.1508972451
    b10 = 0.742380924027
    b11 = 30.789933034
    b12 = 3.99019417011

    y = 0.0
    alnorm = 0.0
    
    if z < 0:
        upper = not upper
        z = - z
    
    if z <= ltone or upper and z <= utzero:
        y = 0.5 * z * z
        if z > con:
            alnorm = b1 * math.exp(- y) / (z - b2 + b3 / (z + b4 + b5 / (z - b6 + b7 / (z + b8 - b9 / (z + b10 + b11 / (z + b12))))))
        else:
            alnorm = 0.5 - z * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))))
    else:
        alnorm = 0.0
    if not upper:
        alnorm = 1.0 - alnorm
    return alnorm

def poisson_cdf (n, lam,lower=True):
    """Poisson CDF evaluater.

    This is a more stable CDF function. It can tolerate large lambda
    value. While the lambda is larger than 700, the function will be a
    little slower.

    Parameters:
    n     : your observation
    lam   : lambda of poisson distribution
    lower : if lower is False, calculate the upper tail CDF
    """
    k = int(n)
    if lam <= 0.0:
        raise Exception("Lambda must > 0")

    if lower:
        if lam > 700:
            return __poisson_cdf_large_lambda (k, lam)
        else:
            return __poisson_cdf(k,lam)
    else:
        if lam > 700:
            return __poisson_cdf_Q_large_lambda (k, lam)
        else:
            return __poisson_cdf_Q(k,lam)

def __poisson_cdf (k,a):
    """Poisson CDF For small lambda. If a > 745, this will return
    incorrect result.

    """
    if k < 0:
        return 0                        # special cases
    next = exp( -a )
    cdf = next
    for i in xrange(1,k+1):
        last = next
        next = last * a / i
        cdf = cdf + next
    if cdf > 1:
        return 1
    else:
        return cdf
    
def __poisson_cdf_large_lambda ( k,a ):
    """Slower poisson cdf for large lambda.
    
    """
    if k < 0:
        return 0                        # special cases
    num_parts = int(a/LSTEP)
    last_part = a % LSTEP
    lastexp = exp(-last_part)
    next = EXPSTEP
    num_parts -= 1
    cdf = next
    for i in xrange(1,k+1):
        last = next
        next = last * a / i
        cdf = cdf + next
        if next > EXPTHRES or cdf > EXPTHRES:
           if num_parts>=1:
               cdf *= EXPSTEP
               next *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1

    for i in xrange(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

def __poisson_cdf_Q (k,a):
    """internal Poisson CDF evaluater for upper tail with small
    lambda.

    """
    if k < 0:
        return 1                        # special cases
    next = exp( -a )

    for i in xrange(1,k+1):
        last = next
        next = last * a / i

    cdf = 0
    i = k+1
    while next >0:
        last = next
        next = last * a / i
        cdf += next
        i+=1
    return cdf

def __poisson_cdf_Q_large_lambda (k,a):
    """Slower internal Poisson CDF evaluater for upper tail with large
    lambda.
    
    """
    if k < 0:
        return 1                        # special cases
    num_parts = int(a/LSTEP)
    last_part = a % LSTEP
    lastexp = exp(-last_part)
    next = EXPSTEP
    num_parts -= 1

    for i in xrange(1,k+1):
        last = next
        next = last * a / i
        if next > EXPTHRES:
           if num_parts>=1:
               next *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1
    cdf = 0
    i = k+1
    while next >0:
        last = next
        next = last * a / i
        cdf += next
        i+=1
        if next > EXPTHRES or cdf > EXPTHRES:
           if num_parts>=1:
               cdf *= EXPSTEP
               next *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1

    for i in xrange(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

def poisson_cdf_inv ( cdf, lam, maximum=1000):
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
    sum2 = 0
    newval = exp( -lam )
    sum2 = newval
#     if cdf <= sum2:
#         return i

    for i in xrange(1,maximum+1):
        sumold = sum2
#         if i == 0:
#             newval = exp( -a )
#             if newval==0:
#                 newval = 4.9406564584124654e-324
#             sum2 = newval
#         else:
        last = newval
        newval = last * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    
    return maximum

def poisson_cdf_Q_inv ( cdf, lam, maximum=1000):
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
    sum2 = 0
    newval = exp( -lam )
    sum2 = newval

    for i in xrange(1,maximum+1):
        sumold = sum2
        last = newval
        newval = last * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    
    return maximum

def poisson_pdf ( k, a ):
    """Poisson PDF.

    PDF(K,A) is the probability that the number of events observed in
    a unit time period will be K, given the expected number of events
    in a unit time.
    """
    if a <= 0:
        return 0
    return exp(-a) * pow (a, k) / factorial (k)
    

def binomial_coef (n,k):
    """BINOMIAL_COEF computes the Binomial coefficient C(N,K)

    n,k are integers.
    """
    mn = min (k, n-k)
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

def binomial_cdf (x, a, b, lower=True):
    """ BINOMIAL_CDF compute the binomial CDF.

    CDF(x)(A,B) is the probability of at most X successes in A trials,
    given that the probability of success on a single trial is B.
    """
    if lower:
        return _binomial_cdf_f (x,a,b)
    else:
        return _binomial_cdf_r (x,a,b)

def _binomial_cdf_r (x,a,b):
    if x < 0:
        return 1
    elif a < x:
        return 0
    elif b == 0:
        return 0
    elif b == 1:
        return 1
    else:
        argmax=int(a*b)
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
            #for i in xrange(argmax,x):
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
            #for i in xrange(x-1,-1,-1):
                pdf*=(a-i)*b/(1-b)/(i+1)
                if pdf==0.0: break
                cdf += pdf
                i+=1
            cdf=min(1,cdf)
            cdf = float("%.10e" %cdf)
            return cdf



def _binomial_cdf_f (x,a,b):
    if x < 0:
        return 0
    elif a < x:
        return 1
    elif b == 0:
        return 1
    elif b == 1:
        return 0
    else:
        argmax=int(a*b)
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


def binomial_cdf_inv ( cdf, a, b ):
    """BINOMIAL_CDF_INV inverts the binomial CDF.

    """
    if cdf < 0 or cdf >1:
        raise Exception("CDF must >= 0 or <= 1")
    cdf2 = 0
    a = int(a)
    for x in xrange(0,a+1):
        pdf = binomial_pdf (x,a,b)
        cdf2 = cdf2 + pdf
        if cdf < cdf2:
            return x
    return x

    
def binomial_pdf( x, a, b):
    """binomial PDF by H. Gene Shin
    
    """
    
    if a<1:
        return 0
    elif x<0 or a<x:
        return 0
    elif b==0:
        if x==0: return 1
        else: return 0
    elif b==1:
        if x==a: return 1
        else: return 0
    else:
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
            #print "pdf:", pdf
        for i in xrange(mx-t):
            pdf *= 1-p
        #print pdf,1-p,mx-t
        #pdf*=(1-p)**(mx-t)
        
        pdf=float("%.10e" %pdf)
        return pdf

def facotrial (n):
    """N!.
    
    """
    if n < 0:
        return 0
    fact = 1
    for i in xrange(2,n+1):
        fact = fact * i
    return fact


# ------------------------------------
# Main function
# ------------------------------------
def test():

#    print binomial_cdf(4553,12988,3.58261e-01)
#    print binomial_cdf(6048,12988,4.84496e-01)
#    print binomial_pdf(2000,5000,0.358)
#    print binomial_pdf(35,1000,0.0358)
#    print binomial_pdf(2,10,0.002)
#    print binomial_cdf(6520,12988,5.96064e-01)
#    print binomial_cdf(15,20,5.96064e-01)
    print ( 1-binomial_cdf(11679,12988,1.77023e-01) )
                       
if __name__ == '__main__':
    test()
    
