# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-03-10 18:39:52 Tao Liu>

"""Module Description: New pileup algorithm based on p-v array
idea. It's modified from original algorithm in MACS2 proposed by Jie
Wang. Now we allow different weights for different genomic ranges, and
simplify the approach.

For a genomic range r_i covering genomic postions from s_i to e_i to
be piled up, we assign the start position a value w_i, and the end
-w_i, so we will have (s_i, w_i) and (e_i, -w_i). Then all N ranges
will be made into an array (2D) of position and weights as:

PV = [ (s_0, w_0), (e_0, -w_0), (s_1, w_1), (e_1, -w_1), ... (s_i, w_i), (e_i, -w_i), ..., (s_N, w_N), (e_N, -w_N) ]

Then the array PV will be sorted by the first dimension, aka the position, no matter the position is from start or end positions from ranges.

PV_sorted = [ (p_0, v_0), (p_1, v_1), ... , (p_i, v_i), ..., (p_{2N}, v_{2N}) ]

The pileup algorithm to produce a bedGraph style pileup (another p-v
array as in Pileup.pyx) can be simply described as:

set the initial pileup z as 0 or a given value, and a start position s
as 0, and an end position e as not-defined.

for i from 0 to 2N in PV_sorted:
    z = z + v_i
    e = p_i
    save the pileup from position s to e is z --  in bedGraph style is to only save (e, z)
    s = e

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).

"""

# ------------------------------------
# python modules
# ------------------------------------

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Signal.BedGraph import bedGraphTrackI

# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
cimport numpy as np
from numpy cimport int32_t, float32_t, int64_t, float64_t, uint32_t, uint64_t
from cpython cimport bool
from cpython cimport PyObject

# ------------------------------------
# C lib
# ------------------------------------
from libc.stdlib cimport malloc, free, qsort

# ------------------------------------
# utility internal functions
# ------------------------------------

cdef float32_t __mapping_function_always_1 ( int32_t L, int32_t R ):
    return 1

cdef void clean_up_ndarray ( np.ndarray x ):
    """ Clean up numpy array in two steps
    """
    cdef:
        long i
    i = x.shape[0] // 2
    x.resize( 100000 if i > 100000 else i, refcheck=False)
    x.resize( 0, refcheck=False)
    return

# ------------------------------------
# public python functions
# ------------------------------------
cpdef np.ndarray pileup_from_LR_hmmratac ( np.ndarray LR_array, mapping_dict ):
    cdef:
        uint64_t l_LR, l_PV, i
        int32_t L, R
        np.ndarray PV, pileup

    l_LR = LR_array.shape[ 0 ]
    l_PV = 2 * l_LR
    PV = np.zeros( shape=l_PV, dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )
    for i in range( l_LR ):
        ( L, R ) = LR_array[ i ]
        PV[ i*2 ] = ( L, mapping_dict[ R - L ] )
        PV[ i*2 + 1 ] = ( R, -1 * mapping_dict[ R - L ] )
    PV.sort( order = 'p' )
    pileup = pileup_PV( PV )
    clean_up_ndarray( PV )
    return pileup    

cpdef np.ndarray pileup_from_LR ( np.ndarray LR_array, mapping_func = __mapping_function_always_1 ):
    """This function will pile up the ndarray containing left and
    right positions, which is typically from PETrackI object. It's
    useful when generating the pileup of a single chromosome is
    needed.

    """
    cdef:
        np.ndarray PV_array, pileup

    PV_array = make_PV_from_LR( LR_array, mapping_func = mapping_func )
    pileup = pileup_PV( PV_array )
    clean_up_ndarray( PV_array )
    return pileup

cpdef np.ndarray pileup_from_PN ( np.ndarray P_array, np.ndarray N_array, int extsize ):
    """This function will pile up the ndarray containing plus and
    minus positions of all reads, which is typically from FWTrackI
    object. It's useful when generating the pileup of a single
    chromosome is needed.

    """
    cdef:
        np.ndarray PV_array, pileup

    PV_array = make_PV_from_PN( P_array, N_array, extsize )
    pileup = pileup_PV( PV_array )
    clean_up_ndarray( PV_array )
    return pileup

# ------------------------------------
# C functions
# ------------------------------------

cdef np.ndarray make_PV_from_LR ( np.ndarray LR_array, mapping_func = __mapping_function_always_1 ):
    """Make sorted PV array from a LR array for certain chromosome in a
    PETrackI object. The V/weight will be assigned as 
    mapping_func( L, R ) or simply 1 if mapping_func is the default.

    LR array is an np.ndarray as with dtype
    [('l','int32'),('r','int32')] with length of N

    PV array is an np.ndarray with
    dtype=[('p','uint32'),('v','float32')] with length of 2N """
    cdef:
        uint64_t l_LR, l_PV, i
        int32_t L, R
        np.ndarray PV

    l_LR = LR_array.shape[ 0 ]
    l_PV = 2 * l_LR
    PV = np.zeros( shape=l_PV, dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )
    for i in range( l_LR ):
        ( L, R ) = LR_array[ i ]
        PV[ i*2 ] = ( L, mapping_func( L, R ) )
        PV[ i*2 + 1 ] = ( R, -1 * mapping_func( L, R ) )
    PV.sort( order = 'p' )
    return PV

cdef np.ndarray make_PV_from_PN ( np.ndarray P_array, np.ndarray N_array, int extsize ):
    """Make sorted PV array from two arrays for certain chromosome in
    a FWTrack object. P_array is for the 5' end positions in plus
    strand, and N_array is for minus strand. We don't support weight
    in this case since all positions should be extended with a fixed
    'extsize'.

    P_array or N_array is an np.ndarray with dtype='int32'

    PV array is an np.ndarray with
    dtype=[('p','uint32'),('v','float32')] with length of 2N """
    cdef:
        uint64_t l_PN, l_PV, i, t
        int32_t L, R
        np.ndarray PV

    l_PN = P_array.shape[ 0 ]
    assert l_PN == N_array.shape[ 0 ]
    l_PV = 4 * l_PN
    PV = np.zeros( shape=l_PV, dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )
    for i in range( l_PN ):
        L = P_array[ i ]
        R = L + extsize
        PV[ i*2 ] = ( L, 1 )
        PV[ i*2 + 1 ] = ( R, -1 )
    for i in range( l_PN ):
        R = N_array[ i ]
        L = R - extsize
        PV[ (l_PN + i)*2 ] = ( L, 1 )
        PV[ (l_PN + i)*2 + 1 ] = ( R, -1 )
    PV.sort( order = 'p' )
    return PV

cdef np.ndarray pileup_PV ( np.ndarray PV_array ):
    """The pileup algorithm to produce a bedGraph style pileup (another
    p-v array as in Pileup.pyx) can be simply described as:
    
    set the initial pileup z as 0 or a given value, and a start
    position s as 0, and an end position e as not-defined.
    
    for i from 0 to 2N in PV_sorted:
        z = z + v_i
        e = p_i
        save the pileup from position s to e is z --  in bedGraph style is to only save (e, z)
        s = e
    """
    cdef:
        float32_t z, v, pre_z
        uint64_t s, e, i, c
        np.ndarray pileup_PV    # this is in bedGraph style as in Pileup.pyx, p is the end of a region, and v is the pileup value
    z = 0
    pre_z = -10000
    s = 0
    pileup_PV = np.zeros( shape=PV_array.shape[0], dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )
    c = 0
    for i in range( PV_array.shape[0] ):
        e = PV_array[i]['p']
        v = PV_array[i]['v']
        if e != s:              # make sure only to record the final value for the same position
            if z == pre_z:      # merge the p-v pair with the previous pair if the same v is found
                pileup_PV[ c-1 ][ 'p' ] = e
            else:
                pileup_PV[ c ] = ( e, z )
                c += 1
                pre_z = z
        z +=  v
        s = e
    pileup_PV.resize( c, refcheck=False )
    #assert z == 0
    return pileup_PV
