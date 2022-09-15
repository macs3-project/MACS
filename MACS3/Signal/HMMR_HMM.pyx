# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-09-15 14:29:10 Tao Liu>

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
import json
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
    rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(random_seed)))
    hmm_model = hmm.GaussianHMM( n_components=3, covariance_type = 'full', random_state = rs )
    #hmm_model = hmm.GMMHMM( n_components = n_states, covariance_type = 'full' )
    hmm_model = hmm_model.fit( training_data, training_data_lengths )
    hmm_model.transmat_ = np.around(hmm_model.transmat_, decimals = 6)
    hmm_model.means_ = np.around(hmm_model.means_, decimals = 6)
    hmm_model.covars_ = np.around(hmm_model.covars_, decimals = 6)
    return hmm_model

cpdef hmm_predict( list signals, list lens, hmm_model ):
    predictions = hmm_model.predict_proba( signals, lens )
    #print( len(predictions), len(signals) )
    #print( sum( lens ) )
    #print( predictions, signals )
    return predictions

cpdef hmm_model_init( model_file ):
    f = open(model_file, 'r')
    model_txt = f.read()
    model_txt = model_txt.replace('  ', ' ').replace('[ ', '[').replace(' ', ',').replace('\n\n\n', ' $ ').replace('\n', '')
    a,b,c,d,e,f,g,h,i = model_txt.split(" $ ")[0:9]
    startprob = np.array(json.loads(a))
    transmat = np.array(json.loads(b))
    means = np.array(json.loads(c))
    covars = np.array(json.loads(d))
    n_features = int(e)
    i_open_region = int(f)
    i_background_region = int(g)
    i_nucleosomal_region = int(h)
    binsize = int(i)

    hmm_model = hmm.GaussianHMM( n_components=3, covariance_type='full' ) #change 3 to variable
    hmm_model.startprob_ = startprob
    hmm_model.transmat_ = transmat
    hmm_model.means_ = means
    hmm_model.covars_ = covars
    hmm_model.n_features = n_features
    return hmm_model, i_open_region, i_background_region, i_nucleosomal_region, binsize
