# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-12-04 18:41:28 Tao Liu>

"""Module for SAPPER BAMParser class

Copyright (c) 2017 Tao Liu <tliu4@buffalo.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
from cpython cimport bool

cimport cython

import numpy as np
cimport numpy as np

ctypedef np.float32_t float32_t
ctypedef np.int32_t int32_t

#from libc.math cimport log10, log, exp, M_LN10 #,fabs,log1p
#from libc.math cimport M_LN10
from math import log1p, exp, log

LN10 = 2.3025850929940458
LN10_tenth = 0.23025850929940458

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef tuple CalModel_Homo( np.ndarray[int32_t, ndim=1] top1_bq_T, np.ndarray[int32_t, ndim=1] top1_bq_C, np.ndarray[int32_t, ndim=1] top2_bq_T, np.ndarray[int32_t, ndim=1] top2_bq_C):
    """Return (lnL, BIC).

    """
    cdef:
        int i
        double lnL, BIC

    lnL=0
    # Phred score is Phred = -10log_{10} E, where E is the error rate.
    # to get the 1-E: 1-E = 1-exp( Phred/-10*M_LN10 ) = 1-exp( Phred * -LOG10_E_tenth )
    for i in range( top1_bq_T.shape[0] ):
        lnL += log1p( -exp(-top1_bq_T[ i ]*LN10_tenth) )
    for i in range( top1_bq_C.shape[0] ):
        lnL += log1p( -exp(-top1_bq_C[ i ]*LN10_tenth) )
        
    for i in range( top2_bq_T.shape[0] ):
        lnL += log( exp(-top2_bq_T[ i ]*LN10_tenth) )
    for i in range( top2_bq_C.shape[0] ):
        lnL += log( exp(-top2_bq_C[ i ]*LN10_tenth) )        

    BIC = -2*lnL                # no free variable, no penalty
    return (lnL, BIC)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef tuple CalModel_Heter_noAS( np.ndarray[int32_t, ndim=1] top1_bq_T,np.ndarray[int32_t, ndim=1] top1_bq_C,np.ndarray[int32_t, ndim=1] top2_bq_T,np.ndarray[int32_t, ndim=1] top2_bq_C ):
    """Return (lnL, BIC)

    k_T
    k_C
    """
    cdef:
        int k_T, k_C
        double lnL, BIC
        int i
        int tn_T, tn_C, tn  # total observed NTs
        double lnL_T, lnL_C # log likelihood for treatment and control

    lnL = 0
    BIC = 0
    #for k_T
    # total oberseved treatment reads from top1 and top2 NTs
    tn_T = top1_bq_T.shape[0] + top2_bq_T.shape[0]
    
    if tn_T == 0:
        raise Exception("Total number of treatment reads is 0!")
    else:
        ( lnL_T, k_T ) = GreedyMaxFunctionNoAS( top1_bq_T.shape[0], top2_bq_T.shape[0], tn_T, top1_bq_T, top2_bq_T )
        lnL += lnL_T
        BIC += -2*lnL_T

    #for k_C
    tn_C = top1_bq_C.shape[0] + top2_bq_C.shape[0]

    if tn_C == 0:
        pass
    else:
        ( lnL_C, k_C )  = GreedyMaxFunctionNoAS( top1_bq_C.shape[0], top2_bq_C.shape[0], tn_C, top1_bq_C, top2_bq_C )
        lnL += lnL_C
        BIC += -2*lnL_C

    tn = tn_C + tn_T

    # we penalize big model depending on the number of reads/samples
    if tn_T == 0:
        BIC += log( tn_C )
    elif tn_C == 0:
        BIC += log( tn_T )
    else:
        BIC += log( tn_T ) + log( tn_C )

    return ( lnL, BIC )


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef tuple CalModel_Heter_AS( np.ndarray[int32_t, ndim=1] top1_bq_T, np.ndarray[int32_t, ndim=1] top1_bq_C, np.ndarray[int32_t, ndim=1] top2_bq_T, np.ndarray[int32_t, ndim=1] top2_bq_C, float max_allowed_ar = 0.99 ):
    """Return (lnL, BIC)

    kc
    ki
    AS_alleleratio
    """
    cdef:
        int k_T, k_C
        double lnL, BIC
        int i
        int tn_T, tn_C, tn  # total observed NTs
        double lnL_T, lnL_C # log likelihood for treatment and control
        double AS_alleleratio   # allele ratio

    lnL = 0
    BIC = 0

    #assert top2_bq_T.shape[0] + top2_bq_C.shape[0] > 0, "Total number of top2 nt should not be zero while using this function: CalModel_Heter_AS!"

    # Treatment
    tn_T = top1_bq_T.shape[0] + top2_bq_T.shape[0]

    if tn_T == 0:
        raise Exception("Total number of treatment reads is 0!")
    else:
        ( lnL_T, k_T, AS_alleleratio ) = GreedyMaxFunctionAS( top1_bq_T.shape[0], top2_bq_T.shape[0], tn_T, top1_bq_T, top2_bq_T, max_allowed_ar)
        #print ">>>",lnL_T, k_T, AS_alleleratio
        lnL += lnL_T
        BIC += -2*lnL_T

    # control
    tn_C = top1_bq_C.shape[0] + top2_bq_C.shape[0]

    if tn_C == 0:
        pass
    else:
        # We assume control will not have allele preference
        ( lnL_C, k_C ) = GreedyMaxFunctionNoAS ( top1_bq_C.shape[0], top2_bq_C.shape[0], tn_C, top1_bq_C, top2_bq_C)
        lnL += lnL_C
        BIC += -2*lnL_C

    tn = tn_C + tn_T

    # we penalize big model depending on the number of reads/samples
    if tn_T == 0:
        BIC += log( tn_C )
    elif tn_C == 0:
        BIC += 2 * log( tn_T )
    else:
        BIC += 2 * log( tn_T ) + log( tn_C )
    
    return (lnL, BIC)


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple GreedyMaxFunctionAS( int m, int n, int tn, np.ndarray[int32_t, ndim=1] me, np.ndarray[int32_t, ndim=1] ne, float max_allowed_ar = 0.99 ):
    """Return lnL, k and alleleratio in tuple.

    Note: I only translate Liqing's C++ code into pyx here. Haven't done any review.
    """
    cdef:
        double dnew, dold, rold, rnew
        int kold, knew
        bool btemp
        int k0
        double dl, dr, d0, d1l, d1r

    assert m+n == tn
    btemp = False
    if tn == 1:                 # only 1 read; I don't expect this to be run...
        dl=calculate_ln(m,n,tn,me,ne,0,0);
        dr=calculate_ln(m,n,tn,me,ne,1,1);

        if dl>dr:
            k = 0
            return ( dl, 0, 0 )
        else:
            k = 1
            return ( dr, 1, 1 )
    elif m == 0:                          #no top1 nt
        return ( calculate_ln( m, n, tn, me, ne, 0, m, max_allowed_ar ), m, 1-max_allowed_ar )
        #k0 = m + 1
    elif m == tn:                         #all reads are top1
        return ( calculate_ln( m, n, tn, me, ne, 1, m, max_allowed_ar ), m, max_allowed_ar )
    else:
        k0 = m

    d0  = calculate_ln( m, n, tn, me, ne, float(k0)/tn, k0, max_allowed_ar )
    d1l = calculate_ln( m, n, tn, me, ne, float(k0-1)/tn, k0-1, max_allowed_ar )
    d1r = calculate_ln( m, n, tn, me, ne, float(k0+1)/tn, k0+1, max_allowed_ar )

    if d0 > d1l-1e-8 and d0 > d1r-1e-8:
        k = k0
        ar = float(k0)/tn
        return ( d0, k, ar )
    elif d1l > d0:
        dold = d1l
        kold = k0-1
        rold = float(k0-1)/tn
        while kold > 1: #disable: when kold=1 still run, than knew=0 is the final run
            knew = kold - 1
            rnew = float(knew)/tn

            dnew = calculate_ln( m,n,tn,me,ne,rnew,knew, max_allowed_ar )

            if(dnew-1e-8 < dold) :
                btemp=True
                break
            kold=knew
            dold=dnew
            rold=rnew

        if btemp: #maximum L value is in [1,m-1];
            k = kold 
            ar= rold
            return ( dold, k, ar )
        else: #L(k=0) is the max for [0,m-1]
            k = kold
            ar = rold
            return ( dold, k, ar )

    elif d1r > d0:
        dold = d1r
        kold = k0 + 1
        rold = float(k0 + 1)/tn
        while kold < tn - 1: #//disable: when kold=tn-1 still run, than knew=tn is the final run
            knew = kold + 1
            
            rnew = float(knew)/tn
            
            dnew = calculate_ln( m,n,tn,me,ne,rnew,knew, max_allowed_ar )
            
            if dnew - 1e-8 < dold: 
                btemp = True
                break
            kold = knew
            dold = dnew
            rold = rnew

        if btemp: #maximum L value is in [m+1,tn-1]
            k = kold
            ar= rold
            return ( dold, k, ar )
        else: #L(k=tn) is the max for [m+1,tn]
            k = kold
            ar = rold
            return ( dold, k, ar )
    else:
        raise Exception("error in GreedyMaxFunctionAS")


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef tuple GreedyMaxFunctionNoAS (int m, int n, int tn, np.ndarray[int32_t, ndim=1] me, np.ndarray[int32_t, ndim=1] ne ):
    """Return lnL, and k in tuple.

    Note: I only translate Liqing's C++ code into pyx here. Haven't done any review.
    """
    cdef:
        double dnew, dold
        int kold, knew
        bool btemp
        int k0
        double bg_r, dl, dr, d0, d1l, d1r

    btemp = False
    bg_r = 0.5

    if tn == 1:
        dl = calculate_ln( m, n, tn, me, ne, bg_r, 0)
        dr= calculate_ln( m, n, tn, me, ne, bg_r, 1)
        if dl > dr:
            k = 0
            return ( dl, 0 )
        else:
            k = 1
            return ( dr, 1 )
    elif m == 0:                          #no top1 nt
        return ( calculate_ln( m, n, tn, me, ne, bg_r, m ), m )
        #k0 = m + 1
    elif m == tn:                         #all reads are top1
        return ( calculate_ln( m, n, tn, me, ne, bg_r, m ), m )
    #elif m == 0:
    #    k0 = m + 1
    #elif m == tn:
    #    k0 = m - 1
    else:
        k0 = m

    d0  = calculate_ln( m, n, tn, me, ne, bg_r, k0)
    d1l = calculate_ln( m, n, tn, me, ne, bg_r, k0 - 1)
    d1r = calculate_ln( m, n, tn, me, ne, bg_r, k0 + 1)

    if d0 > d1l - 1e-8 and d0 > d1r - 1e-8:
        k = k0
        return ( d0, k )
    elif d1l > d0:
        dold = d1l
        kold=k0 - 1
        while kold >= 1:   #//when kold=1 still run, than knew=0 is the final run
            knew = kold - 1
            dnew = calculate_ln( m, n, tn, me, ne, bg_r, knew )
            if dnew - 1e-8 < dold:
                btemp = True
                break
            kold=knew
            dold=dnew
                        
        if btemp: #//maximum L value is in [1,m-1];
            k = kold
            return ( dold, k )
        else: #//L(k=0) is the max for [0,m-1]
            k = kold
            return ( dold, k )
    elif d1r > d0:
        dold = d1r
        kold = k0 + 1
        while kold <= tn - 1: #//when kold=tn-1 still run, than knew=tn is the final run
            knew = kold + 1
            dnew = calculate_ln( m, n, tn, me, ne, bg_r, knew )
            if dnew - 1e-8 < dold:
                btemp = True
                break
            kold = knew
            dold = dnew

        if btemp: #//maximum L value is in [m+1,tn-1]
            k = kold
            return ( dold, k )
        else: #//L(k=tn) is the max for [m+1,tn]
            k = kold
            return ( dold, k )
    else:
        raise Exception("error in GreedyMaxFunctionNoAS")

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef calculate_ln( int m, int n, int tn, np.ndarray[int32_t, ndim=1] me, np.ndarray[int32_t, ndim=1] ne, double r, int k, float max_allowed_r = 0.99):
    """Calculate log likelihood given quality of top1 and top2, the ratio r and the observed k.

    """
    cdef:
        int i
        double lnL
        double e

    lnL = 0

    if r > max_allowed_r or r < 1 - max_allowed_r: # r is extremely high or
        #print "here1"
        lnL += k*log( max_allowed_r ) + (tn-k)*log( 1- max_allowed_r) #-10000
    else:
        lnL += k*log( r ) + (tn-k)*log(1-r)

    if k == 0 or k == tn:       # it's entirely biased toward 1 allele
        #print "here2"
        pass
        #lnL += k*log( max_allowed_r ) #-10000
        #lnL += -10000
    elif k <= tn/2:
        for i in range( k ):
            lnL += log(float(tn-i)/(k-i))
    else:
        for i in range( tn-k ):
            lnL += log(float(tn-i)/(tn-k-i))

    for i in range( m ):
        e = exp( - me[ i ] * LN10_tenth )
        lnL += log((1-e)*(float(k)/tn) + e*(1-float(k)/tn))
    
    for i in range( n ):
        e = exp( - ne[ i ] * LN10_tenth )
        lnL += log((1-e)*(1-float(k)/tn) + e*(float(k)/tn))

    #print r,k,lnL
    return lnL

cpdef int calculate_GQ ( double lnL1, double lnL2, double lnL3):
    """GQ1 = -10*log_{10}((L2+L3)/(L1+L2+L3))
       
       
    """
    cdef:
        double L1, L2, L3, sum, tmp
        int GQ_score
        
    #L1 = exp(lnL1-lnL1)
    L1 = 1
    L2 = exp(lnL2-lnL1)
    L3 = exp(lnL3-lnL1)

    #if L1 > 1:
    #    L1 = 1

    if L2 > 1:
        L2 = 1
    if L3 > 1:
        L3 = 1
    #if(L1<1e-110) L1=1e-110;
    if L2 < 1e-110:
        L2=1e-110
    if L3 < 1e-110:
        L3 = 1e-110

    sum = L1 + L2 + L3
    tmp = ( L2 + L3 )/sum
    if tmp > 1e-110:
        GQ_score = (int)(-4.34294*log(tmp))
    else:
        GQ_score = 255

    return GQ_score;

cpdef int calculate_GQ_heterASsig( double lnL1, double lnL2):
    """
    """
    cdef:
        double L1, L2, sum, tmp
        int ASsig_score

    #L1=exp(2.7182818,lnL1-lnL1)
    L1 = 1
    L2 = exp( lnL2 - lnL1 )

    #if L1 > 1:
    #    L1 = 1
    if L2 > 1:
        L2 = 1
    #if L1 < 1e-110:
    #    L1 = 1e-110
    if L2 < 1e-110:
        L2 = 1e-110

    sum = L1 + L2
    tmp = L2/sum
    if tmp > 1e-110:
        ASsig_score = (int)(-4.34294*log(tmp))
    else:
        ASsig_score = 255

    return ASsig_score

