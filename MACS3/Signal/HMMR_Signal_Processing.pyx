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

cpdef dict generate_weight_mapping( list fraglen_list, list means, list stddevs ):
    cdef:
        dict ret_mapping
    ret_mapping = {}
    return ret_mapping

cpdef dict generate_digested_signals( object petrack, dict weight_mapping ):
    cdef:
        dict ret_digested_signals
    ret_digested_signals = {}
    return ret_digested_signals

cpdef dict extract_signals_from_training_regions( dict signals, object peaks, binsize = 10 ):
    cdef:
        dict ret_training_data
    ret_training_data = {}
    return ret_training_data


