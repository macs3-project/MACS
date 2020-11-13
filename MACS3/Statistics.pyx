# cython: language_level=3
# cython: profile=True
# Time-stamp: <2019-10-30 17:27:50 taoliu>

"""Module Description: Statistics function wrappers.

NOTE: This file is no longer used in any other MACS2 codes.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

from libc.math cimport log10, log
from MACS2.Prob import poisson_cdf

cdef class P_Score_Upper_Tail:
    """Unit to calculate -log10 poisson_CDF of upper tail and cache
    results in a hashtable.
    """
    cdef:
        dict pscore_dict
        
    def __init__ ( self ):
        self.pscore_dict = dict()
            
    cpdef float get_pscore ( self, int x, float l ):
        cdef:
            float val

        if ( x, l ) in self.pscore_dict:
            return self.pscore_dict [ (x, l ) ]
        else:
            # calculate and cache
            val = -1 * poisson_cdf ( x, l, False, True )
            self.pscore_dict[ ( x, l ) ] = val
            return val

cdef class LogLR_Asym:
    """Unit to calculate asymmetric log likelihood, and cache
    results in a hashtable.
    """
    cdef:
        dict logLR_dict

    def __init__ ( self ):
        self.logLR_dict = dict()
            
    cpdef float get_logLR_asym ( self, float x, float y ):
        cdef:
            float val

        if ( x, y ) in self.logLR_dict:
            return self.logLR_dict[ ( x, y ) ]
        else:
            # calculate and cache
            if x > y:
                val = (x*(log10(x)-log10(y))+y-x)
            elif x < y:
                val = (x*(-log10(x)+log10(y))-y+x)
            else:
                val = 0
            self.logLR_dict[ ( x, y ) ] = val
            return val


cdef class LogLR_Sym:
    """Unit to calculate symmetric log likelihood, and cache
    results in a hashtable.

    "symmetric" means f(x,y) = -f(y,x)
    """
    cdef:
        dict logLR_dict
        
    def __init__ ( self ):
        self.logLR_dict = dict()
            
    cpdef float get_logLR_sym ( self, float x, float y ):
        cdef:
            float val

        if ( x, y ) in self.logLR_dict:
            return self.logLR_dict[ ( x, y ) ]
        else:
            # calculate and cache
            if x > y:
                val = (x*(log10(x)-log10(y))+y-x)
            elif y > x:
                val = (y*(log10(x)-log10(y))+y-x)
            else:
                val = 0

            self.logLR_dict[ ( x, y ) ] = val
            return val


cdef class LogLR_Diff:
    """Unit to calculate log likelihood for differential calling, and cache
    results in a hashtable.

    here f(x,y) = f(y,x) and f() is always positive.
    """
    cdef:
        dict logLR_dict
        
    def __init__ ( self ):
        self.logLR_dict = dict()
            
    cpdef float get_logLR_diff ( self, float x, float y ):
        cdef:
            float val

        if ( x, y ) in self.logLR_dict:
            return self.logLR_dict[ ( x, y ) ]
        else:
            # calculate and cache
            if y > x: y, x = x, y
            if x == y:
                val = 0
            else:
                val = (x*(log10(x)-log10(y))+y-x)

            self.logLR_dict[ ( x, y ) ] = val
            return val
