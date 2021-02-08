# cython: language_level=3
# cython: profile=True
# Time-stamp: <2021-02-05 16:31:53 Tao Liu>

"""Module Description: For pileup functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from time import time as ttime

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.cPosValCalculation cimport single_end_pileup as c_single_end_pileup
from MACS3.Signal.cPosValCalculation cimport write_pv_array_to_bedGraph as c_write_pv_array_to_bedGraph
from MACS3.Signal.cPosValCalculation cimport PosVal
from MACS3.Signal.cPosValCalculation cimport quick_pileup as c_quick_pileup

# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
cimport numpy as np
from numpy cimport int32_t, float32_t
from cpython cimport bool
from cpython cimport PyObject

# ------------------------------------
# C lib
# ------------------------------------
from libc.stdlib cimport malloc, free, qsort

# ------------------------------------
# utility internal functions
# ------------------------------------
cdef inline float mean( float a, float b ):
    return ( a + b ) / 2

cdef void clean_up_ndarray ( np.ndarray x ):
    """ Clean up numpy array in two steps
    """
    cdef:
        long i
    i = x.shape[0] // 2
    x.resize( 100000 if i > 100000 else i, refcheck=False)
    x.resize( 0, refcheck=False)
    return

cdef np.ndarray[np.int32_t, ndim=1] fix_coordinates(np.ndarray[np.int32_t, ndim=1] poss, int rlength):
    """Fix the coordinates.
    """
    cdef:
        long i
        int32_t * ptr

    ptr = <int32_t *> poss.data

    # fix those negative coordinates
    for i in range( poss.shape[0] ):
        if ptr[i] < 0:
            ptr[i] = 0
        else:
            break

    # fix those over-boundary coordinates
    for i in range( poss.shape[0]-1, -1, -1 ):
        if ptr[i] > rlength:
            ptr[i] = rlength
        else:
            break
    return poss

# ------------------------------------
# functions
# ------------------------------------

# ------------------------------------
# functions for pileup_cmd only
# ------------------------------------

# This function uses pure C code for pileup
cpdef void pileup_and_write_se( trackI,
                        bytes output_filename,
                        int d,
                        float scale_factor,
                        float baseline_value = 0.0,
                        bint directional = True,
                        bint halfextension = True ):
    """ Pileup a FWTrackI object and write the pileup result into a
    bedGraph file.

    This function calls pure C functions from cPosValCalculation.

    This function is currently only used `macs3 pileup` cmd.
    """
    cdef:
        long five_shift, three_shift, l, i
        list chroms
        int n_chroms
        bytes chrom
        np.ndarray[np.int32_t, ndim=1] plus_tags, minus_tags
        dict chrlengths = trackI.get_rlengths ()
        int * plus_tags_pos
        int * minus_tags_pos
        long rlength
        long fl
        bytes py_bytes
        char * chrom_char
        PosVal * _data
        long l_data

    # This block should be reused to determine the actual shift values
    if directional:
        # only extend to 3' side
        if halfextension:
            five_shift = d//-4  # five shift is used to move cursor towards 5' direction to find the start of fragment
            three_shift = d*3//4 # three shift is used to move cursor towards 3' direction to find the end of fragment
        else:
            five_shift = 0
            three_shift = d
    else:
        # both sides
        if halfextension:
            five_shift = d//4
            three_shift = five_shift
        else:
            five_shift = d//2
            three_shift = d - five_shift
    # end of the block

    chroms = list(chrlengths.keys())
    n_chroms = len( chroms )

    fh = open(output_filename, "w")
    fh.write("")
    fh.close()

    for i in range( n_chroms ):
        chrom = chroms[ i ]
        (plus_tags, minus_tags) = trackI.get_locations_by_chr(chrom)
        rlength = <long> chrlengths[ chrom ]
        plus_tags_pos = <int *> plus_tags.data
        minus_tags_pos = <int *> minus_tags.data

        _data = c_single_end_pileup( plus_tags_pos, plus_tags.shape[0], minus_tags_pos, minus_tags.shape[0], five_shift, three_shift, 0, rlength, scale_factor, baseline_value, &l_data )

        # write
        py_bytes = chrom
        chrom_char = py_bytes
        c_write_pv_array_to_bedGraph( _data, l_data, chrom_char, output_filename, 1 )

        # clean
        free( _data )
    return

# function to pileup BAMPE/BEDPE stored in PETrackI object and write to a BEDGraph file
# this function uses c function
cpdef pileup_and_write_pe( petrackI,
                           bytes output_filename,
                           float scale_factor = 1,
                           float baseline_value = 0.0):
    """ Pileup a PETrackI object and write the pileup result into a
    bedGraph file.

    This function calls pure C functions from cPosValCalculation.

    This function is currently only used `macs3 pileup` cmd.
    """
    cdef:
        dict chrlengths = petrackI.get_rlengths ()
        list chroms
        int n_chroms
        int i
        bytes chrom

        np.ndarray locs
        np.ndarray[np.int32_t, ndim=1] locs0
        np.ndarray[np.int32_t, ndim=1] locs1

        int * start_pos
        int * end_pos
        long fl
        bytes py_bytes
        char * chrom_char
        PosVal * _data
        long l_data

    chroms = list(chrlengths.keys())
    n_chroms = len( chroms )

    fh = open(output_filename, "w")
    fh.write("")
    fh.close()

    for i in range( n_chroms ):
        chrom = chroms[ i ]
        locs = petrackI.get_locations_by_chr(chrom)

        locs0 = np.sort(locs['l'])
        locs1 = np.sort(locs['r'])
        start_pos = <int *> locs0.data
        end_pos   = <int *> locs1.data

        _data = c_quick_pileup ( start_pos, end_pos, locs0.shape[0], scale_factor, baseline_value, &l_data )

        # write
        py_bytes = chrom
        chrom_char = py_bytes
        c_write_pv_array_to_bedGraph( _data, l_data, chrom_char, output_filename, 1 )

        # clean
        free( _data )        
    return

# ------------------------------------
# functions for other codes
# ------------------------------------

# general pileup function implemented in cython
cpdef list se_all_in_one_pileup ( np.ndarray[np.int32_t, ndim=1] plus_tags, np.ndarray[np.int32_t, ndim=1] minus_tags, long five_shift, long three_shift, int rlength, float scale_factor, float baseline_value ):
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
    functions within MACS3.

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

    tmp = [ret_p, ret_v]        # for (endpos,value)
    i_s = 0                     # index of start_poss
    i_e = 0                     # index of end_poss
    I = 0

    pileup = 0
    if start_poss.shape[0] == 0: return tmp
    pre_p = min(start_poss_ptr[0],end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = max(0,baseline_value)
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
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
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
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
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
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    # clean mem
    clean_up_ndarray( start_poss )
    clean_up_ndarray( end_poss )

    # resize
    ret_p.resize( I, refcheck=False )
    ret_v.resize( I, refcheck=False )

    return tmp

# quick pileup implemented in cython
cpdef list quick_pileup ( np.ndarray[np.int32_t, ndim=1] start_poss, np.ndarray[np.int32_t, ndim=1] end_poss, float scale_factor, float baseline_value ):
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
        ret_v_ptr[0] = max(0,baseline_value)
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1

    pre_v = pileup

    while i_s < ls and i_e < le:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
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
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
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
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    ret_p.resize( I, refcheck=False )
    ret_v.resize( I, refcheck=False )

    return tmp

# quick pileup implemented in cython
cpdef list naive_quick_pileup ( np.ndarray[np.int32_t, ndim=1] sorted_poss, int extension):
    """Simple pileup, every tag will be extended left and right with length `extension`.

    Note: Assumption is that `poss` has to be pre-sorted! There is no
    check on whether it's sorted.
    """
    cdef:
        long i_s, i_e, i, I
        int a, b, p, pre_p, pileup
        long l = sorted_poss.shape[0]
        np.ndarray[np.int32_t, ndim=1] start_poss
        np.ndarray[np.int32_t, ndim=1] end_poss
        np.ndarray[np.int32_t, ndim=1] ret_p
        np.ndarray[np.float32_t, ndim=1] ret_v
        # pointers are used for numpy arrays
        int32_t * start_poss_ptr
        int32_t * end_poss_ptr
        int32_t * ret_p_ptr     # pointer for position array
        float32_t * ret_v_ptr     # pointer for value array

    start_poss = sorted_poss - extension
    start_poss[ start_poss < 0 ] = 0
    end_poss = sorted_poss + extension

    start_poss_ptr = <int32_t *> start_poss.data
    end_poss_ptr = <int32_t *> end_poss.data

    
    ret_p = np.zeros( 2*l, dtype="int32" )
    ret_v = np.zeros( 2*l, dtype="float32" )

    ret_p_ptr = <int32_t *> ret_p.data
    ret_v_ptr = <float32_t *> ret_v.data

    i_s = 0                         # index of plus_tags
    i_e = 0                         # index of minus_tags
    I = 0
    
    pileup = 0
    if l == 0: return ret
    pre_p = min(start_poss_ptr[0], end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = 0
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1

    pre_v = pileup

    while i_s < l and i_e < l:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = pileup
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
                ret_v_ptr[0] = pileup
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

    # add rest of end positions
    if i_e < l:
        for i in range(i_e, l):
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = pileup
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    ret_p.resize( I, refcheck=False )
    ret_v.resize( I, refcheck=False )

    return [ ret_p, ret_v ]

# general function to compare two pv arrays in cython.
cpdef list over_two_pv_array ( list pv_array1, list pv_array2, func="max" ):
    """Merge two position-value arrays. For intersection regions, take
    the maximum value within region.

    pv_array1 and pv_array2 are [p,v] type lists, same as the output
    from quick_pileup function. 'p' and 'v' are numpy arrays of int32
    and float32.

    available operations are 'max', 'min', and 'mean'
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

    if func == "max":
        f = max
    elif func == "min":
        f = min
    elif func == "mean":
        f = mean
    else:
        raise Exception("Invalid function")

    [ a1_pos, a1_v ] = pv_array1
    [ a2_pos, a2_v ] = pv_array2
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
        ret_v_ptr[0] =  f( a1_v_ptr[0], a2_v_ptr[0] )
        I += 1
        if a1_pos_ptr[0] < a2_pos_ptr[0]:
            # clip a region from pre_p to p1, then set pre_p as p1.
            ret_pos_ptr[0] = a1_pos_ptr[0]
            ret_pos_ptr += 1
            ret_v_ptr += 1
            pre_p = a1_pos_ptr[0]
            # call for the next p1 and v1
            a1_pos_ptr += 1
            a1_v_ptr += 1
            i1 += 1
        elif a1_pos_ptr[0] > a2_pos_ptr[0]:
            # clip a region from pre_p to p2, then set pre_p as p2.
            ret_pos_ptr[0] = a2_pos_ptr[0]
            ret_pos_ptr += 1
            ret_v_ptr += 1
            pre_p = a2_pos_ptr[0]
            # call for the next p1 and v1
            a2_pos_ptr += 1
            a2_v_ptr += 1
            i2 += 1
        else:
            # from pre_p to p1 or p2, then set pre_p as p1 or p2.
            ret_pos_ptr[0] = a1_pos_ptr[0]
            ret_pos_ptr += 1
            ret_v_ptr += 1
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

cpdef naive_call_peaks ( list pv_array, float min_v, float max_v = 1e30, int max_gap = 50, int min_length = 200 ):
    cdef:
        int peak_length, pre_p, p, i, summit, tstart, tend
        long x                  # index used for searching the first peak
        double v, summit_value, tvalue
        bytes chrom
        set chrs
        list peak_content       # (pre_p, p, v) for each region in the peak
        list ret_peaks = []     # returned peak summit and height

    peak_content = []
    ( ps,vs ) = pv_array
    psn = iter(ps).__next__         # assign the next function to a viable to speed up
    vsn = iter(vs).__next__
    x = 0
    pre_p = 0                   # remember previous position
    while True:
        # find the first region above min_v
        try:                    # try to read the first data range for this chrom
            p = psn()
            v = vsn()
        except:
            break
        x += 1                  # index for the next point
        if v > min_v:
            peak_content = [ ( pre_p , p, v ), ]
            pre_p = p
            break               # found the first range above min_v
        else:
            pre_p = p

    for i in range( x, len( ps ) ):
        # continue scan the rest regions
        p = psn()
        v = vsn()
        if v <= min_v: # not be detected as 'peak'
            pre_p = p
            continue
        # for points above min_v
        # if the gap is allowed
        # gap = pre_p - peak_content[-1][1] or the dist between pre_p and the last p
        if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
            peak_content.append( ( pre_p, p, v ) )
        else:
            # when the gap is not allowed, close this peak IF length is larger than min_length
            if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:
                __close_peak( peak_content, ret_peaks, max_v, min_length )
            # reset and start a new peak
            peak_content = [ ( pre_p, p, v ), ]
        pre_p = p

    # save the last peak
    if peak_content:
        if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:        
            __close_peak( peak_content, ret_peaks, max_v, min_length )
    return ret_peaks

cdef void __close_peak( peak_content, peaks, float max_v, int min_length ):
    """Internal function to find the summit and height

    If the height is larger than max_v, skip
    """
    tsummit = []
    summit = 0
    summit_value = 0
    for (tstart,tend,tvalue) in peak_content:
        if not summit_value or summit_value < tvalue:
            tsummit = [int((tend+tstart)/2),]
            summit_value = tvalue
        elif summit_value == tvalue:
            tsummit.append( int((tend+tstart)/2) )
    summit = tsummit[int((len(tsummit)+1)/2)-1 ]
    if summit_value < max_v:
        peaks.append( (summit, summit_value) )
    return
