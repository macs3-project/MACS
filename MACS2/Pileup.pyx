# Time-stamp: <2015-04-20 14:26:53 Tao Liu>

"""Module Description: For pileup functions.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

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
from libc.stdlib cimport malloc, free, qsort

from MACS2.IO.FixWidthTrack import FWTrack
from MACS2.IO.PairedEndTrack import PETrackI
from MACS2.IO.BedGraph import bedGraphTrackI
from MACS2.Constants import *

import numpy as np
cimport numpy as np

from numpy cimport int32_t
ctypedef np.float32_t float32_t

from cpython cimport bool
from cpython cimport PyObject

from cPosValCalculation cimport single_end_pileup as c_single_end_pileup
from cPosValCalculation cimport write_pv_array_to_bedGraph as c_write_pv_array_to_bedGraph
from cPosValCalculation cimport PosVal
from cPosValCalculation cimport quick_pileup_simple

from cython.parallel import *

from time import time as ttime

# ------------------------------------
# functions
# ------------------------------------

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline long long_max(long a, long b): return a if a >= b else b
cdef inline float float_max(float a, float b): return a if a >= b else b

# This function uses pure C code for pileup
cpdef pileup_and_write( trackI,
                        output_filename,
                        int d,
                        float scale_factor,
                        float baseline_value = 0.0,
                        bint directional = True, 
                        bint halfextension = True ):

    cdef:
        long five_shift, three_shift, l, i
        list chroms
        int n_chroms
        str chrom
        np.ndarray[np.int32_t, ndim=1] plus_tags, minus_tags
        dict chrlengths = trackI.get_rlengths ()
        int * plus_tags_pos
        int * minus_tags_pos
        long rlength
        long fl
        bytes py_bytes
        char * chrom_char
        PosVal ** _data
        long * l_data

    
    # This block should be reused to determine the actual shift values
    if directional:
        # only extend to 3' side
        if halfextension:
            five_shift = d/-4  # five shift is used to move cursor towards 5' direction to find the start of fragment
            three_shift = d*3/4 # three shift is used to move cursor towards 3' direction to find the end of fragment
        else:
            five_shift = 0
            three_shift = d
    else:
        # both sides
        if halfextension:
            five_shift = d/4
            three_shift = five_shift
        else:
            five_shift = d/2
            three_shift = d - five_shift
    # end of the block

    chroms = chrlengths.keys()
    n_chroms = len( chroms )
    _data = < PosVal ** > malloc( n_chroms * sizeof( PosVal * ) )
    l_data = < long * > malloc( n_chroms * sizeof( long ) )
  
    for i in range( n_chroms ):
        chrom = chroms[ i ]
        (plus_tags, minus_tags) = trackI.get_locations_by_chr(chrom)
        rlength = <long> chrlengths[ chrom ]
        plus_tags_pos = <int *> plus_tags.data
        minus_tags_pos = <int *> minus_tags.data
        _data[ i ] = c_single_end_pileup( plus_tags_pos, plus_tags.shape[0], minus_tags_pos, minus_tags.shape[0], five_shift, three_shift, 0, rlength, scale_factor, baseline_value, &l_data[ i ] )

    for i in range( n_chroms ):
        chrom = chroms[ i ]
        py_bytes = chrom.encode()
        chrom_char = py_bytes
        c_write_pv_array_to_bedGraph( _data[ i ], l_data[ i ], chrom_char, output_filename, 1 )

    free( l_data )
    free( _data )

# Unified pileup function #
cpdef unified_pileup_bdg(track,
                         ds,
                         scale_factors,
                         float baseline_value = 0.0,
                         bint directional = True, 
                         bint halfextension = True):
    """This function will call corresponding function for FWTrack
    or PETrackI to pileup fragments.

    It will call pileup_w_multiple_d* functions for control input, then
    calculate maximum values; or call normal pileup functions for
    treatment data.

    """

    chrs = track.get_chr_names()
    if type(ds) is list:
        # multiple pileup (e.g. control with 1k, 10k and d extension)
        if isinstance(track, FWTrack):
            return pileup_w_multiple_d_bdg(track, ds, scale_factors,
                                           baseline_value,
                                           directional, halfextension)
        elif isinstance(track, PETrackI):
            return pileup_w_multiple_d_bdg_pe(track, ds, scale_factors,
                                              baseline_value)
        else:
            raise ValueError("track must be of type FWTrack or PETrackI")
    else:
        # single extension (e.g. treatment data)
        if isinstance(track, FWTrack):
            return pileup_bdg_se(track, ds, scale_factors,
                                 baseline_value,
                                 directional, halfextension)
        elif isinstance(track, PETrackI):
            if ds is None:      # do not extend pair end library
                return pileup_bdg_pe(track, scale_factors, baseline_value)
            else:               # still extend pair end library centered at middle point of fragment.
                return pileup_bdg_pe_w_ext(track, ds, scale_factors,
                                           baseline_value)
        else:
            raise ValueError("track must be of type FWTrack or PETrackI")

## Fixed-width functions for single end library##
cdef pileup_bdg_se(object trackI, int d,
                    float scale_factor = 1.0,
                    float baseline_value = 0.0,
                    bool directional = True,
                    bool halfextension = True ):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards 3' side with size of d if directional is Ture,
    or both sides with d/2 if directional is False.

    A tag is a single genomic location.

    trackI  : A FWTrack object with raw plus and minus 5' end positions
    d       : tag will be extended to this value to 3' direction, unless directional is False.
    baseline_value : a value to be filled for missing values.
    directional: if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
    halfextension: only make a fragment of d/2 size centered at fragment center

    Return a bedGraphTrackI object.
    """
    cdef:
        long five_shift, three_shift, l
        int rlength
        str chrom
        Ends ends
        np.ndarray[np.int32_t, ndim=1] plus_tags, minus_tags
        dict chrlengths = trackI.get_rlengths ()

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.


    if directional:
        # only extend to 3' side
        if halfextension:
            five_shift = d/-4  # five shift is used to move cursor towards 5' direction to find the start of fragment
            three_shift = d*3/4 # three shift is used to move cursor towards 3' direction to find the end of fragment
        else:
            five_shift = 0
            three_shift = d
    else:
        # both sides
        if halfextension:
            five_shift = d/4
            three_shift = five_shift
        else:
            five_shift = d/2
            three_shift = d - five_shift

    for chrom in sorted(chrlengths.keys()):
        rlength = chrlengths[chrom]
        (plus_tags, minus_tags) = trackI.get_locations_by_chr(chrom)

        ends = start_and_end_poss( plus_tags, minus_tags, five_shift, three_shift , rlength)

        ret.add_a_chromosome( chrom, quick_pileup ( ends.startposs, ends.endposs, scale_factor, baseline_value ) )

        # free mem
        ends.startposs.resize(100000, refcheck=False)
        ends.startposs.resize(0, refcheck=False)
        ends.endposs.resize(100000, refcheck=False)
        ends.endposs.resize(0, refcheck=False)                

    return ret

cdef pileup_w_multiple_d_bdg(object trackI, list d_s, list scale_factor_s = [],
                              float baseline_value = 0.0,
                              bool directional = True,
                              bool halfextension = True):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards 3' side with size of d if directional is Ture,
    or both sides with d/2 if directional is False.

    A tag is a single genomic location.

    trackI  : A FWTrack object with raw plus and minus 5' end positions
    d       : tag will be extended to this value to 3' direction, unless directional is False.
    baseline_value : a value to be filled for missing values.
    directional: if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
    halfextension: only make a fragment of d/2 size centered at fragment center

    Return a bedGraphTrackI object.
    """
    cdef:
        long d, five_shift, three_shift, l, i  
        float scale_factor
        dict chrlengths = trackI.get_rlengths()
        Ends ends

    assert len(d_s) == len(scale_factor_s), "Arguments d_s and scale_factor_s should have the same length!"

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    chrs = trackI.get_chr_names()       

    five_shift_s = []
    three_shift_s = []

    for d in d_s:
        if directional:
            # only extend to 3' side
            if halfextension:
                five_shift_s.append(d/-4)  # five shift is used to move cursor towards 5' direction to find the start of fragment
                three_shift_s.append(d*3/4) # three shift is used to move cursor towards 3' direction to find the end of fragment
            else:
                five_shift_s.append(0)
                three_shift_s.append(d)
        else:
            # both sides
            if halfextension:
                five_shift_s.append(d/4)
                three_shift_s.append(d/4)
            else:
                five_shift_s.append(d/2)
                three_shift_s.append(d - d/2)

    for chrom in sorted(chrlengths.keys()):
        rlength = chrlengths[chrom]
        (plus_tags,minus_tags) = trackI.get_locations_by_chr(chrom)

        prev_pileup = None

        for i in range(len(d_s)):
            
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]

            ends = start_and_end_poss( plus_tags, minus_tags, five_shift, three_shift, rlength )

            tmp_pileup = quick_pileup ( ends.startposs, ends.endposs, scale_factor, baseline_value )

            # free mem
            ends.startposs.resize(100000, refcheck=False)
            ends.startposs.resize(0, refcheck=False)
            ends.endposs.resize(100000, refcheck=False)
            ends.endposs.resize(0, refcheck=False)                            
            
            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        ret.add_a_chromosome( chrom, prev_pileup )

    return ret

## Paired-end functions ##

# baseline_value needs to be float not int, otherwise we cause error in 
# poisson CDF
cdef pileup_bdg_pe(object trackI, float scale_factor, float baseline_value):
    """Pileup fragments into bedGraphTrackI object.

    trackI  : A PETrackI object with genomic locations
    baseline_value : a value to be filled for missing values.
    scale_factor : value to be multiplied at each position

    Return a bedGraphTrackI object.
    """
    cdef:
        int rlength
        str chrom
        np.ndarray[np.int32_t, ndim=2] locs
        dict chrlengths = trackI.get_rlengths ()
        
    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.
    for chrom in sorted(chrlengths.keys()):
        rlength = chrlengths[chrom]
        locs = trackI.get_locations_by_chr(chrom)
        ret.add_a_chromosome(chrom, quick_pileup(locs[:,0], locs[:,1],
                                                 scale_factor, 
                                                 baseline_value))
    return ret

cdef pileup_bdg_pe_w_ext (object trackI, int d, float scale_factor = 1.0,
                           float baseline_value = 0.0):
    """Pileup fragments into bedGraphTrackI object with extension. Fragment will
    be extended both directions from midpoint by distance d/2, or the original
    width will be used if d = 0

    trackI  : A PETrackI object with genomic locations
    d       : fragments will be extended by 1/2 this value in both directions
              from their midpoint
    baseline_value : a value to be filled for missing values.
    scale_factor : value to be multiplied at each position

    Return a bedGraphTrackI object.
    """
    cdef:
        int five_shift, three_shift
        int rlength
        str chrom
        np.ndarray[np.int32_t, ndim=2] locs
        np.ndarray[np.int32_t, ndim=1] start_poss, end_poss
        dict chrlengths = trackI.get_rlengths ()
        
    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    five_shift = d/2
    three_shift = d - five_shift

    for chrom in sorted(chrlengths.keys()):
        rlength = chrlengths[chrom]
        locs = trackI.get_locations_by_chr(chrom)
        midpoints = locs[:,0] + (locs[:,1] - locs[:,0]) / 2

        # fix negative coordinations
        start_poss = midpoints - five_shift
        start_poss = fix_coordinates(start_poss, rlength)
        end_poss = midpoints + three_shift
        end_poss = fix_coordinates( end_poss, rlength) 

        pileup = quick_pileup ( start_poss, end_poss, scale_factor,
                                baseline_value )
        ret.add_a_chromosome( chrom, pileup )

        # free mem
        start_poss.resize(100000, refcheck=False)
        start_poss.resize(0, refcheck=False)
        end_poss.resize(100000, refcheck=False)
        end_poss.resize(0, refcheck=False)                

    return ret

cdef pileup_w_multiple_d_bdg_pe ( object trackI, list d_s = [],  
                                   list scale_factor_s = [],
                                   float baseline_value = 0):
    """Pileup fragments into bedGraphTrackI object with extension. Fragment will
    be extended by d / 2 in both directions from midpoint

    trackI  : A PETrackI object
    d_s       : tag will be extended by 1/2 this value in both directions by
                each d and multiplied by corresponding scale_factor
    baseline_value : a value to be filled for missing values.
    scale_factor_s: same length as d_s, scale factor for each d
    scale_factor_0: scale factor for original fragments

    Return a bedGraphTrackI object.
    """
    cdef:
        long d, five_shift, three_shift, l, i
        float scale_factor
        dict chrlengths = trackI.get_rlengths ()

    assert len(d_s) == len(scale_factor_s), "Arguments d_s and scale_factor_s should have the same length!"

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    chrs = trackI.get_chr_names()       

    five_shift_s = [d / 2 for d in d_s[1:]]
    three_shift_s = [d - d / 2 for d in d_s[1:]]

    for chrom in sorted(chrlengths.keys()):
        rlength = chrlengths[chrom]
        locs = trackI.get_locations_by_chr(chrom)
        midpoints = locs[:,0] + (locs[:,1] - locs[:,0]) / 2

        prev_pileup = quick_pileup(locs[:,0], locs[:,1],
                                   scale_factor_s[0], baseline_value)

        for i in range(len(five_shift_s)):
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i + 1]

            # fix negative coordinations
            start_poss = midpoints - five_shift
            start_poss = fix_coordinates(start_poss, rlength)
            end_poss = midpoints + three_shift
            end_poss = fix_coordinates( end_poss, rlength) 

            tmp_pileup = quick_pileup(start_poss, end_poss, scale_factor,
                                      baseline_value)

            # free mem
            start_poss.resize(100000, refcheck=False)
            start_poss.resize(0, refcheck=False)
            end_poss.resize(100000, refcheck=False)
            end_poss.resize(0, refcheck=False)                            
            
            prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            
        ret.add_a_chromosome( chrom, prev_pileup )

    return ret

cdef class Ends:
    cdef:
        np.ndarray startposs
        np.ndarray endposs 

cdef start_and_end_poss ( np.ndarray plus_tags, np.ndarray minus_tags,
                          long five_shift, long three_shift, int rlength):
    cdef:
        long i
        long lp = plus_tags.shape[0]
        long lm = minus_tags.shape[0]
        long l = lp + lm

    start_poss = np.concatenate( ( plus_tags-five_shift, minus_tags-three_shift ) )
    end_poss   = np.concatenate( ( plus_tags+three_shift, minus_tags+five_shift ) )    

    # sort
    start_poss.sort()
    end_poss.sort()
    
    # fix negative coordinations and those extends over end of chromosomes
    ends = Ends()
    ends.startposs = fix_coordinates(start_poss, rlength)
    ends.endposs = fix_coordinates(end_poss, rlength)

    return ends

cdef np.ndarray[np.int32_t, ndim=1] fix_coordinates(np.ndarray[np.int32_t, ndim=1] poss, int rlength):
    cdef:
        long i
        int32_t * ptr
    
    ptr = <int32_t *> poss.data

    for i in range( poss.shape[0] ):
        if ptr[i] < 0:
            ptr[i] = 0
        else:
            break
        
    for i in range( poss.shape[0]-1, -1, -1 ):
        if ptr[i] > rlength:
            ptr[i] = rlength
        else:
            break
    
    return poss

cdef int * fix_coordinates_2 ( int * poss, int l_of_poss, int rlength) nogil:
    cdef long i
    
    for i in range( l_of_poss ):
        if poss[i] < 0:
            poss[i] = 0
        else:
            break
        
    for i in range( l_of_poss-1, -1, -1 ):
        if poss[i] > rlength:
            poss[i] = rlength
        else:
            break
    
    return poss

# general pileup function
cpdef se_all_in_one_pileup ( np.ndarray[np.int32_t, ndim=1] plus_tags, np.ndarray[np.int32_t, ndim=1] minus_tags, long five_shift, long three_shift, int rlength, float scale_factor, float baseline_value ):
    """Return pileup given 5' end of fragment at plus or minus strand
    separately, and given shift at both direction to recover a
    fragment. This function is for single end sequencing library
    only. Please directly use 'quick_pileup' function for Pair-end
    library.
    
    It contains a super-fast and simple algorithm proposed by Jie
    Wang. It will take sorted start positions and end positions, then
    compute pileup values.

    It will return a pileup result in similar structure as
    bedGraph. There are two python arrays:
    
    [end positions, values] or '[p,v] array' in other description for
    functions within MACS2.

    Two arrays have the same length and can be matched by index. End
    position at index x (p[x]) record continuous value of v[x] from
    p[x-1] to p[x].

    """
    cdef:
        long i_s, i_e, i, I
        int a, b, p, pre_p, pileup
        list tmp
        long lp = plus_tags.shape[0]
        long lm = minus_tags.shape[0]
        long l = lp + lm
        np.ndarray[np.int32_t, ndim=1] start_poss, end_poss, ret_p
        np.ndarray[np.float32_t, ndim=1] ret_v
        # pointers are used for numpy arrays
        int32_t * start_poss_ptr
        int32_t * end_poss_ptr
        int32_t * ret_p_ptr     # pointer for position array 
        float32_t * ret_v_ptr     # pointer for value array
        

    start_poss = np.concatenate( ( plus_tags-five_shift, minus_tags-three_shift ) )
    end_poss   = np.concatenate( ( plus_tags+three_shift, minus_tags+five_shift ) )    

    # sort
    start_poss.sort()
    end_poss.sort()
    
    # fix negative coordinations and those extends over end of chromosomes
    start_poss = fix_coordinates(start_poss, rlength)
    end_poss = fix_coordinates(end_poss, rlength)

    lx = start_poss.shape[0]

    start_poss_ptr = <int32_t *> start_poss.data
    end_poss_ptr = <int32_t *> end_poss.data

    ret_p = np.zeros( 2 * lx, dtype="int32" )
    ret_v = np.zeros( 2 * lx, dtype="float32" )

    ret_p_ptr = <int32_t *> ret_p.data
    ret_v_ptr = <float32_t *> ret_v.data

    tmp = [ret_p, ret_v] # for (endpos,value)
    # #tmppadd = tmp[0].append
    # #tmpvadd = tmp[1].append
    i_s = 0                         # index of start_poss
    i_e = 0                         # index of end_poss
    I = 0

    pileup = 0
    if start_poss.shape[0] == 0: return tmp
    pre_p = min(start_poss_ptr[0],end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = float_max(0,baseline_value) 
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1
        
    pre_v = pileup

    assert start_poss.shape[0] == end_poss.shape[0]
    lx = start_poss.shape[0]

    while i_s < lx and i_e < lx:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = float_max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup += 1
            i_s += 1
            start_poss_ptr += 1
        elif start_poss_ptr[0] > end_poss_ptr[0]:
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = float_max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            i_e += 1
            end_poss_ptr += 1
        else:
            i_s += 1
            i_e += 1
            start_poss_ptr += 1
            end_poss_ptr += 1

    if i_e < lx:
        # add rest of end positions
        for i in range(i_e, lx):
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = float_max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    # clean mem
    start_poss.resize(100000, refcheck=False)
    start_poss.resize(0, refcheck=False)
    end_poss.resize(100000, refcheck=False)
    end_poss.resize(0, refcheck=False)   

    # resize
    ret_p.resize( I, refcheck=False )
    ret_v.resize( I, refcheck=False )

    return tmp

cdef int compare(const void * a, const void * b) nogil:
    if a - b > 0: return 1
    return 0

cpdef quick_pileup ( np.ndarray[np.int32_t, ndim=1] start_poss, np.ndarray[np.int32_t, ndim=1] end_poss, float scale_factor, float baseline_value ):
    """Return pileup given plus strand and minus strand positions of fragments.
    
    A super-fast and simple algorithm proposed by Jie Wang. It will
    take sorted start positions and end positions, then compute pileup
    values. 

    It will return a pileup result in similar structure as
    bedGraph. There are two python arrays:
    
    [end positions, values] or [p,v]

    Two arrays have the same length and can be matched by index. End
    position at index x (p[x]) record continuous value of v[x] from
    p[x-1] to p[x].

    """
    cdef:
        long i_s, i_e, i, I
        int a, b, p, pre_p, pileup
        list tmp
        long ls = start_poss.shape[0]
        long le = end_poss.shape[0]
        long l = ls + le
        np.ndarray[np.int32_t, ndim=1] ret_p
        np.ndarray[np.float32_t, ndim=1] ret_v
        # pointers are used for numpy arrays
        int32_t * start_poss_ptr
        int32_t * end_poss_ptr
        int32_t * ret_p_ptr     # pointer for position array 
        float32_t * ret_v_ptr     # pointer for value array
        #int max_pileup = 0

    start_poss_ptr = <int32_t *> start_poss.data
    end_poss_ptr = <int32_t *> end_poss.data

    ret_p = np.zeros( l, dtype="int32" )
    ret_v = np.zeros( l, dtype="float32" )

    ret_p_ptr = <int32_t *> ret_p.data
    ret_v_ptr = <float32_t *> ret_v.data

    tmp = [ret_p, ret_v] # for (endpos,value)

    i_s = 0                         # index of plus_tags
    i_e = 0                         # index of minus_tags
    I = 0

    pileup = 0
    if ls == 0: return tmp
    pre_p = min(start_poss_ptr[0], end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = float_max(0,baseline_value) 
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1
        
    pre_v = pileup
    
    while i_s < ls and i_e < le:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = float_max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup += 1
            #if pileup > max_pileup:
            #    max_pileup = pileup
            i_s += 1
            start_poss_ptr += 1
        elif start_poss_ptr[0] > end_poss_ptr[0]:
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = float_max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            i_e += 1
            end_poss_ptr += 1
        else:
            i_s += 1
            i_e += 1
            start_poss_ptr += 1
            end_poss_ptr += 1

    if i_e < le:
        # add rest of end positions
        for i in range(i_e, le):
            p = end_poss_ptr[0]
            #for p in minus_tags[i_e:]:
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = float_max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    ret_p.resize( I, refcheck=False )
    ret_v.resize( I, refcheck=False )

    return tmp

# general function to calculate maximum between two arrays.

cpdef list max_over_two_pv_array ( list tmparray1, list tmparray2 ):
    """Merge two position-value arrays. For intersection regions, take
    the maximum value within region.

    tmparray1 and tmparray2 are [p,v] type lists, same as the output
    from quick_pileup function. 'p' and 'v' are numpy arrays of int32
    and float32.
    """
    
    cdef:
        int pre_p, p1, p2
        float v1, v2
        np.ndarray[np.int32_t, ndim=1] a1_pos, a2_pos, ret_pos
        np.ndarray[np.float32_t, ndim=1] a1_v, a2_v, ret_v
        int32_t * a1_pos_ptr
        int32_t * a2_pos_ptr
        int32_t * ret_pos_ptr
        float32_t * a1_v_ptr
        float32_t * a2_v_ptr
        float32_t * ret_v_ptr
        long l1, l2, l, i1, i2, I

    [ a1_pos, a1_v ] = tmparray1
    [ a2_pos, a2_v ] = tmparray2
    ret_pos = np.zeros( a1_pos.shape[0] + a2_pos.shape[0], dtype="int32" )
    ret_v   = np.zeros( a1_pos.shape[0] + a2_pos.shape[0], dtype="float32" )

    a1_pos_ptr = <int32_t *> a1_pos.data
    a1_v_ptr = <float32_t *> a1_v.data
    a2_pos_ptr = <int32_t *> a2_pos.data
    a2_v_ptr = <float32_t *> a2_v.data
    ret_pos_ptr = <int32_t *> ret_pos.data
    ret_v_ptr = <float32_t *> ret_v.data

    l1 = a1_pos.shape[0]
    l2 = a2_pos.shape[0]

    i1 = 0
    i2 = 0
    I = 0

    pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret

    while i1 < l1 and i2 < l2:
        if a1_pos_ptr[0] < a2_pos_ptr[0]:
            # clip a region from pre_p to p1, then set pre_p as p1.
            ret_pos_ptr[0] = a1_pos_ptr[0]
            ret_v_ptr[0] =  float_max( a1_v_ptr[0], a2_v_ptr[0] )
            ret_pos_ptr += 1
            ret_v_ptr += 1
            I += 1
            pre_p = a1_pos_ptr[0]
            # call for the next p1 and v1
            a1_pos_ptr += 1
            a1_v_ptr += 1
            i1 += 1
        elif a1_pos_ptr[0] > a2_pos_ptr[0]:
            # clip a region from pre_p to p2, then set pre_p as p2.
            ret_pos_ptr[0] = a2_pos_ptr[0]
            ret_v_ptr[0] =  float_max( a1_v_ptr[0], a2_v_ptr[0] )
            ret_pos_ptr += 1
            ret_v_ptr += 1
            I += 1
            pre_p = a2_pos_ptr[0]
            # call for the next p1 and v1
            a2_pos_ptr += 1
            a2_v_ptr += 1
            i2 += 1
        else:
            # from pre_p to p1 or p2, then set pre_p as p1 or p2.
            ret_pos_ptr[0] = a1_pos_ptr[0]
            ret_v_ptr[0] =  float_max( a1_v_ptr[0], a2_v_ptr[0] )
            ret_pos_ptr += 1
            ret_v_ptr += 1
            I += 1
            pre_p = a1_pos_ptr[0]
            # call for the next p1, v1, p2, v2.
            a1_pos_ptr += 1
            a1_v_ptr += 1
            i1 += 1
            a2_pos_ptr += 1
            a2_v_ptr += 1
            i2 += 1

    ret_pos.resize( I, refcheck=False )
    ret_v.resize( I, refcheck=False )
    return [ret_pos, ret_v]


