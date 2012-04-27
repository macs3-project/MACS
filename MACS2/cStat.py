# Time-stamp: <2012-02-29 15:09:27 Tao Liu>

"""Module Description

Copyright (c) 2012 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

from array import array as pyarray
from MACS2.Constants import *
from random import gammavariate as rgamma
from random import seed as rseed
from math import log
import pymc
from pymc import deterministic
# ------------------------------------
# constants
# ------------------------------------
import numpy.random as numpyrand

LOG2E = log(2.718281828459045,2)        # for converting natural log to log2

gfold_dict = {}                         # temporarily save all precomputed gfold

# ------------------------------------
# Misc functions
# ------------------------------------

# BUGFIX FOR PYMC ARGUMENT CHANGE
from inspect import getargspec
PROGRESS_BAR_ENABLED = 'progress_bar' in getargspec(pymc.MCMC.sample)[0]

def MCMCPoissonPosteriorRatio (sample_number, burn, count1, count2):
    """MCMC method to calculate ratio distribution of two Posterior Poisson distributions.

    sample_number: number of sampling. It must be greater than burn, however there is no check.
    burn: number of samples being burned.
    count1: observed counts of condition 1
    count2: observed counts of condition 2

    return: list of log2-ratios
    """
    lam1 = pymc.Uniform('U1',0,10000)   # prior of lambda is uniform distribution
    lam2 = pymc.Uniform('U2',0,10000)   # prior of lambda is uniform distribution    
    poi1 = pymc.Poisson('P1',lam1,value=count1,observed=True) # Poisson with observed value count1
    poi2 = pymc.Poisson('P2',lam2,value=count2,observed=True) # Poisson with observed value count2
    @deterministic
    def ratio (l1=lam1,l2=lam2):
        return log(l1,2) - log(l2,2)
    mcmcmodel  = pymc.MCMC([ratio,lam1,poi1,lam2,poi2])
    mcmcmodel.use_step_method(pymc.AdaptiveMetropolis,[ratio,lam1,lam2,poi1,poi2], delay=20000)
    if PROGRESS_BAR_ENABLED:
        mcmcmodel.sample(iter=sample_number, progress_bar=False, burn=burn)    
    else:
        mcmcmodel.sample(iter=sample_number, burn=burn)    
    return ratio.trace()


def MLEPoissonPosteriorRatio (sample_number, burn, count1, count2):
    """MLE method to calculate ratio distribution of two Posterior Poisson distributions.

    MLE of Posterior Poisson is Gamma(k+1,1) if there is only one observation k.

    sample_number: number of sampling. It must be greater than burn, however there is no check.
    burn: number of samples being burned.
    count1: observed counts of condition 1
    count2: observed counts of condition 2

    return: list of log2-ratios
    """
    rseed(1)
    ratios = pyarray('f',[])
    ra = ratios.append
    for i in xrange(sample_number):
        x1 = rgamma(count1+1,1)
        x2 = rgamma(count2+1,1)
        ra( log(x1,2) - log(x2,2) )
    return ratios[int(burn):]

def get_gfold ( v1, v2, precompiled_get=None, cutoff=0.01, sample_number=10000, burn=500, offset=0, mcmc=False):    
    # try cached gfold in this module first
    if gfold_dict.has_key((v1,v2)):
        return gfold_dict[(v1,v2)]

    # calculate ratio+offset

    # first, get the value from precompiled table
    try:
        V = precompiled_get( v1, v2 )
        if v1 > v2:
            # X >= 0
            ret = max(0,V+offset)
        elif v1 < v2:
            # X < 0
            ret = min(0,V+offset)
        else:
            ret = 0.0
        
    except IndexError:
        if mcmc:
            numpyrand.seed([10])
            P_X = MCMCPoissonPosteriorRatio(sample_number,burn,v1,v2)
            i = int( (sample_number-burn) * cutoff)
        else:
            P_X = MLEPoissonPosteriorRatio(sample_number,0,v1,v2)
            i = int(sample_number * cutoff)            

        P_X = map(lambda x:x+offset,sorted(P_X))
        P_X_mean = float(sum(P_X))/len(P_X)
        
        if P_X_mean >= 0:
            # X >= 0
            ret = max(0,P_X[i])
        elif P_X_mean < 0:
            # X < 0
            ret = min(0,P_X[-1*i])
            
        #print v1,v2,P_X_mean,'-',offset,ret,i,P_X[i],P_X[-i]
    

    gfold_dict[(v1,v2)] = ret
    return ret

#def convert_gfold ( v, cutoff = 0.01, precompiled_gfold=None, mcmc=False ):
def convert_gfold ( v, precompiled_gfold, sample_number=15000, burn=5000, offset=0, cutoff=0.01, mcmc=False):
    """Take (name, count1, count2), try to extract precompiled gfold
    from precompiled_gfold.get; if failed, calculate the gfold using
    MCMC if mcmc is True, or simple MLE solution if mcmc is False.
    """
    ret = []
    retadd = ret.append
    get_func = precompiled_gfold.get
    for i in xrange(len(v[0])):
        rid= v[0][i]
        v1 = int(v[1][i])
        v2 = int(v[2][i])
        # calculate gfold from precompiled table, MCMC or MLE
        gf = get_gfold(v1,v2,precompiled_get=get_func,cutoff=cutoff,sample_number=sample_number,burn=burn,offset=offset,mcmc=mcmc)
        retadd([rid,gf])
    return ret

# ------------------------------------
# Classes
# ------------------------------------
