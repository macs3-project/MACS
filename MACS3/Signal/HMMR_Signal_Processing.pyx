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

cpdef list generate_weight_mapping( list fraglen_list, list means, list stddevs ):
    """
    return: list of dict, with key as fraglen, value as [ w_s, w_m, w_d, w_t ]
    """
    cdef:
        list ret_mapping
        list variances
        int l
        float m_s, m_m, m_d, m_t
        float v_s, v_m, v_d, v_t
        float p_s, p_m, p_d, p_t
        float w_s, w_m, w_d, w_t
        float s
        int i, j
    assert len(means) == 4
    assert len(stddevs) == 4
    [m_s, m_m, m_d, m_t] = means
    [v_s, v_m, v_d, v_t] = [ x**2 for x in stddevs ]
    ret_mapping = [ {}, {}, {}, {} ]
    for i in range( len(fraglen_list) ):
        l = fraglen_list[ i ]
        p_s = pnorm2( float(l), m_s, v_s )
        p_m = pnorm2( float(l), m_m, v_m )
        p_d = pnorm2( float(l), m_d, v_d )
        p_t = pnorm2( float(l), m_t, v_t )
        s = p_s + p_m + p_d + p_t
        w_s = p_s / s
        w_m = p_m / s
        w_d = p_d / s
        w_t = p_t / s
        ret_mapping[ 0 ][ l ] = w_s
        ret_mapping[ 1 ][ l ] = w_m
        ret_mapping[ 2 ][ l ] = w_d
        ret_mapping[ 3 ][ l ] = w_t
    return ret_mapping

cpdef dict generate_digested_signals( object petrack, list weight_mapping ):
    cdef:
        dict ret_digested_signals
    ret_digested_signals = {}
    return ret_digested_signals

cpdef dict extract_signals_from_training_regions( dict signals, object peaks, binsize = 10 ):
    cdef:
        dict ret_training_data
    ret_training_data = {}
    return ret_training_data


