# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-04-15 13:43:07 Tao Liu>

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


cpdef hmm_training( list training_data, list training_data_lengths, int n_states = 3, int random_seed = 12345 ):
    # training data should be in array like format: X = np.array([[.5, .3, .1, .1], [.6, .4, 0, 0]])
    # if we do not want init_prob to be updated through learning, set params = 'tmc' and init_prob = initial_state otherwise it will be overwritten
    # according to base documentation, if init_prob not stated, it is set to be equally likely for any state (1/ # of components)
    # if we have other known parameters, we should set these (ie: means_weights, covariance_type etc.)
    hmm_model = hmm.GaussianHMM( n_components=3, covariance_type = 'full', random_state = random_seed )
    #hmm_model = hmm.GMMHMM( n_components = n_states, covariance_type = 'full' )
    hmm_model.fit( training_data, training_data_lengths )
    return hmm_model

cpdef hmm_predict( list signals, list lens, hmm_model ):
    predictions = hmm_model.predict_proba( signals, lens )
    return predictions
