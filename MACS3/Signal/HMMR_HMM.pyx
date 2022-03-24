# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-02-23 17:37:36 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
from math import sqrt
import logging
debug   = logging.debug
info    = logging.info
# ------------------------------------
# Other modules
# ------------------------------------

import numpy as np
cimport numpy as np
from cpython cimport bool
from hmmlearn import hmm
# from hmmlearn cimport hmm
# from sklearn.cluster import KMeans 
# from sklearn.cluster cimport KMeans

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.Prob import pnorm2

# ------------------------------------
# Misc functions
# ------------------------------------

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
cdef inline float get_weighted_density( int x, float m, float v, w ):
    """Description:
    
    parameters:
      1. x: the observed value
      2. m: the mean of gaussian
      3. v: the variance of the gaussian
      4. w: the weight
    return value:
    """
    return w * pnorm2( float(x), m, v )

# ------------------------------------
# Classes
# ------------------------------------
# ------------------------------------
# public functions
# ------------------------------------

# cpdef initial_state_kmeans( list training_data, k = 3 ):
    # trial for kmeans:
    # kmeans = KMeans(n_clusters = k, random_state = 0).fit(training_data)
    # labels = kmeans.labels_.tolist()
    # ls = []
    # for i in range (k):
    #     ls.append(labels.count(i))
    # initial_state = [x / sum(ls) for x in ls]
    # initial_means = kmeans.cluster_centers_
    # initial_cov = np.cov(initial_means)
    # return initial_state, initial_means, initial_cov

cpdef hmm_training( list training_data, k = 3 ):
    hmm_model = hmm.GaussianHMM( n_components = k )
    hmm_model.fit( training_data )
    return hmm_model

cpdef hmm_predict( list signals, hmm_model, binsize = 10 ):
    predictions = hmm_model.predict( signals )
    return predictions


