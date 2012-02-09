# Time-stamp: <2012-02-09 17:08:36 Tao Liu>

"""Module Description

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

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

#from collections import Counter
from array import array as pyarray
from MACS2.Constants import *
from random import gammavariate as rgamma
from random import seed as rseed
from math import log
import pymc
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
# def histogram ( vl, breaks=None, minv=None, maxv=None, binsize=None):
#     """Return histogram statistics.

#     Parameters:

#     vl: 2D numpy.array as [ [value, length], [value, length], ...]
    
#     breaks: if breaks is not None and a valid integar, split [min,max]
#     of values in vl into number of equal sized bins. Otherwise, no
#     binning is involved.

#     Return Value:
#     Counter object


#     when breaks is not None, key values in Counter is the start points
#     of each bin.
    
#     """
#     assert breaks == None or isinstance(breaks,int)
    
#     ret = Counter()

#     if breaks == None and binsize == None:
#         for (v,l) in vl:
#             ret[v] += int(l)
#     else:
#         if maxv == None:
#             maxv = vl[:,0].max()
#         if minv == None:
#             minv = vl[:,0].min()
#         if binsize == None:
#             binsize = (maxv-minv)/breaks
#         for (v,l) in vl:
#             k = (v - minv)//binsize*binsize + minv
#             #print k
#             ret[ k ] += int(l)

#     return ret

# def histogram2D ( md ):
#     """Return histogram statistics.

#     Parameters:

#     vl: 2D numpy.array as [ [value, length], [value, length], ...]
    
#     breaks: if breaks is not None and a valid integar, split [min,max]
#     of values in vl into number of equal sized bins. Otherwise, no
#     binning is involved.

#     Return Value:
#     Counter object


#     when breaks is not None, key values in Counter is the start points
#     of each bin.
    
#     """
#     ret = Counter()

#     for (m, d, l) in md:
#         ret[ (m,d) ] += int(l)

#     return ret

def MCMCGammaSamplingRatio (sample_number, alpha1, alpha2, beta1, beta2):
    gamma1 = pymc.Gamma('G1',alpha1,beta1)
    gamma2 = pymc.Gamma('G2',alpha2,beta2)
    logratio = pymc.log(gamma1)-pymc.log(gamma2)
    model  = pymc.MCMC([gamma1,gamma2])
    model.seed()
    #model.sample(iter=int(sample_number*1.5), burn=int(sample_number/2), progress_bar=False)
    model.sample(iter=sample_number, progress_bar=False)    
    x1 = gamma1.trace()
    x2 = gamma2.trace()
    return map(lambda x,y:log(x,2)-log(y,2), x1, x2)

def SimpleGammaSamplingRatio ( sample_number, alpha1, alpha2, beta1, beta2):
    ret = pyarray(FBYTE4,[])

    for i in xrange(sample_number):
        y1 = rgamma(alpha1,beta1)
        y2 = rgamma(alpha2,beta2)
        ret.append(log(y1,2)-log(y2,2))

    return ret

gfold_dict = {}
rseed(10)

def get_gfold ( v1, v2, cutoff, mcmc=False ):
    sample_number = 1000

    if gfold_dict.has_key((v1,v2)):
        return gfold_dict[(v1,v2)]

    if mcmc:
        P_X = MCMCGammaSamplingRatio(sample_number,v1+1,v2+1,1,1)
    else:
        P_X = SimpleGammaSamplingRatio(sample_number,v1+1,v2+1,1,1)

    P_X = sorted(P_X)


    i = int(sample_number * cutoff)
    
    if v1 > v2:
        # X >= 0
        ret = max(0,P_X[i])
    else:
        # X < 0
        ret = min(0,P_X[-1*i])        
    gfold_dict[(v1,v2)] = ret
    return ret

def convert_gfold ( v, cutoff = 0.01, mcmc=False ):
    ret = []
    for i in xrange(len(v[0])):
        rid= v[0][i]
        v1 = v[1][i]
        v2 = v[2][i]
        gf = get_gfold (v1,v2,cutoff, mcmc=mcmc)
        ret.append((rid,gf))
    return ret

# ------------------------------------
# Classes
# ------------------------------------
