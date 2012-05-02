 # cython: profile=True
# Time-stamp: <2012-05-01 18:34:03 Tao Liu>

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
import sys
from array import array

from MACS2.IO.cBedGraph import bedGraphTrackI
from MACS2.Constants import *

from libc.stdint cimport int32_t

import numpy as np
cimport numpy as np

from cpython cimport bool

# ------------------------------------
# functions
# ------------------------------------

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline long long_max(long a, long b): return a if a >= b else b
cdef inline float float_max(float a, float b): return a if a >= b else b

## Fixed-width functions ##
def pileup_bdg (trackI, int d, float baseline_value = 0, bool directional = True, bool halfextension = True, float scale_factor = 1, rlengths = None):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards 3' side with size of d if directional is Ture,
    or both sides with d/2 if directional is False.

    A tag is a single genomic location.

    trackI  : A FWTrackIII object with raw plus and minus 5' end positions
    d       : tag will be extended to this value to 3' direction, unless directional is False.
    baseline_value : a value to be filled for missing values.
    directional: if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
    halfextension: only make a fragment of d/2 size centered at fragment center

    Return a bedGraphTrackI object.
    """
    cdef long five_shift, three_shift, l, i, j, i_s, i_e, p, pre_p
    cdef int32_t rlength = 2147483647
    #cdef int * start_poss
    #cdef int * end_poss    

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    chrs = trackI.get_chr_names()       

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

    for chrom in chrs:
        if rlengths is not None: rlength = rlengths[chrom] 
        (plus_tags,minus_tags) = trackI.get_locations_by_chr(chrom)

        l = len(plus_tags)+len(minus_tags)        

        ( start_poss, end_poss ) = start_and_end_poss( plus_tags, minus_tags, five_shift, three_shift , rlength)

        ret.add_a_chromosome( chrom, pileup_a_chromosome ( start_poss, end_poss, l, scale_factor, baseline_value ) )

        # free mem
        start_poss.resize(100000, refcheck=False)
        start_poss.resize(0, refcheck=False)
        end_poss.resize(100000, refcheck=False)
        end_poss.resize(0, refcheck=False)                

    return ret

def pileup_w_multiple_d_bdg ( trackI, d_s, float baseline_value = 0, bool directional = True, bool halfextension = True, scale_factor_s = [] , rlengths = None):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards 3' side with size of d if directional is Ture,
    or both sides with d/2 if directional is False.

    A tag is a single genomic location.

    trackI  : A FWTrackIII object with raw plus and minus 5' end positions
    d       : tag will be extended to this value to 3' direction, unless directional is False.
    baseline_value : a value to be filled for missing values.
    directional: if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
    halfextension: only make a fragment of d/2 size centered at fragment center

    Return a bedGraphTrackI object.
    """
    cdef long d, five_shift, three_shift, l, i, j, i_s, i_e, p, pre_p
    cdef int32_t rlength = 2147483647
    cdef float scale_factor
    #cdef int * start_poss
    #cdef int * end_poss

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

    for chrom in chrs:
        if rlengths is not None: rlength = rlengths[chrom]
        (plus_tags,minus_tags) = trackI.get_locations_by_chr(chrom)

        l = len(plus_tags)+len(minus_tags)

        prev_pileup = None

        for i in range(len(d_s)):
            
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]

            #start_poss = build_start_poss( plus_tags, minus_tags, five_shift, three_shift, l )
            #end_poss = build_end_poss( plus_tags, minus_tags, five_shift, three_shift, l )
            (start_poss, end_poss) = start_and_end_poss( plus_tags, minus_tags, five_shift, three_shift, rlength )

            #print chrom, baseline_value
            tmp_pileup = pileup_a_chromosome ( start_poss, end_poss, l, scale_factor, baseline_value )

            # free mem
            #del(start_poss)
            #del(end_poss)
            start_poss.resize(100000, refcheck=False)
            start_poss.resize(0, refcheck=False)
            end_poss.resize(100000, refcheck=False)
            end_poss.resize(0, refcheck=False)                            
            
            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        ret.add_a_chromosome( chrom, prev_pileup )

    return ret

## Paired-end functions ##

# baseline_value needs to be float not int, otherwise we cause error in 
# poisson CDF
def pileup_frag_bdg (trackI, float baseline_value = 0, float scale_factor = 1.0, rlengths = None):
    """Pileup fragments into bedGraphTrackI object.

    trackI  : A PETrackI object with genomic locations
    baseline_value : a value to be filled for missing values.
    scale_factor : value to be multiplied at each position

    Return a bedGraphTrackI object.
    """
    cdef long five_shift, three_shift, l, i, j, i_s, i_e, p, pre_p
    cdef int32_t rlength = 2147483647
    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.
    chrs = trackI.get_chr_names()
    for chrom in chrs:
        starts = trackI.get_locations_by_chr(chrom)[:,0]
        ends = trackI.get_locations_by_chr(chrom)[:,1]
        l = len(starts)
        assert len(ends) == l, "tracks do not have the same length"
        pileup = pileup_a_chromosome(starts, ends, l, scale_factor, 
                                     baseline_value)
        ret.add_a_chromosome(chrom, pileup)
    return ret

def pileup_and_ext_frag_bdg (trackI, int d, float baseline_value = 0, float scale_factor = 1, rlengths = None):
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
    cdef long five_shift, three_shift, l, i, j, i_s, i_e, p, pre_p
    cdef int32_t rlength = 2147483647
    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.
    chrs = trackI.get_chr_names()       

    five_shift = d/2
    three_shift = d - five_shift

    for chrom in chrs:
        if rlengths is not None: rlength = rlengths[chrom]
        starts = trackI.get_locations_by_chr(chrom)[:,0]
        ends = trackI.get_locations_by_chr(chrom)[:,1]
        l = len(starts)
        assert len(ends) == l, "tracks do not have the same length"        
        midpoints = (starts + ends) / 2

        start_poss = midpoints - five_shift
        end_poss   = midpoints + three_shift    
    
        # fix negative coordinations
        start_poss = fix_coordinates(start_poss, rlength)
        end_poss = fix_coordinates(end_poss, rlength) 

        pileup = pileup_a_chromosome ( start_poss, end_poss, l, scale_factor,
                                       baseline_value )
        ret.add_a_chromosome( chrom, pileup )

        # free mem
        start_poss.resize(100000, refcheck=False)
        start_poss.resize(0, refcheck=False)
        end_poss.resize(100000, refcheck=False)
        end_poss.resize(0, refcheck=False)                

    return ret

def pileup_frag_w_multiple_d_bdg ( trackI, list d_s = [], float baseline_value = 0, 
                                   list scale_factor_s = [], float scale_factor_0 = 1.0 ,
                                   rlengths = None):
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
    cdef long d, five_shift, three_shift, l, i, j, i_s, i_e, p, pre_p
    cdef float scale_factor
    cdef int32_t rlength = 2147483647
    #cdef int * start_poss
    #cdef int * end_poss

    assert len(d_s) == len(scale_factor_s), "Arguments d_s and scale_factor_s should have the same length!"

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    chrs = trackI.get_chr_names()       

    five_shift_s = [d / 2 for d in d_s]
    three_shift_s = [d - d / 2 for d in d_s]

    for chrom in chrs:
        if rlengths is not None: rlength = rlengths[chrom]
        starts = trackI.get_locations_by_chr(chrom)[:,0]
        ends = trackI.get_locations_by_chr(chrom)[:,1]
        l = len(starts)
        assert len(ends) == l, "tracks do not have the same length"        
        midpoints = (starts + ends) / 2

        prev_pileup = pileup_a_chromosome ( starts, ends, l, scale_factor_0, baseline_value )

        for i in range(len(d_s)):
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]

            start_poss = midpoints - five_shift
            end_poss   = midpoints + three_shift    
        
            # fix negative coordinations
            start_poss = fix_coordinates(start_poss, rlength)
            end_poss = fix_coordinates(end_poss, rlength)

            #print chrom, baseline_value
            tmp_pileup = pileup_a_chromosome ( start_poss, end_poss, l, scale_factor, baseline_value )

            # free mem
            #del(start_poss)
            #del(end_poss)
            start_poss.resize(100000, refcheck=False)
            start_poss.resize(0, refcheck=False)
            end_poss.resize(100000, refcheck=False)
            end_poss.resize(0, refcheck=False)                            
            
            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        ret.add_a_chromosome( chrom, prev_pileup )

    return ret


cdef start_and_end_poss ( np.ndarray plus_tags, np.ndarray minus_tags, long five_shift, long three_shift, int32_t rlength = 2147483647):
    cdef long i
    cdef long lp = plus_tags.shape[0]
    cdef long lm = minus_tags.shape[0]
    cdef long l = lp + lm

    start_poss = np.concatenate( ( plus_tags-five_shift, minus_tags-three_shift ) )
    end_poss   = np.concatenate( ( plus_tags+three_shift, minus_tags+five_shift ) )    

    # sort
    start_poss.sort()
    end_poss.sort()
    
    # fix negative coordinations and over end of chromosomes
    start_poss = fix_coordinates(start_poss, rlength)
    end_poss = fix_coordinates(end_poss, rlength)

    return (start_poss, end_poss)

cdef fix_coordinates(np.ndarray poss, int32_t rlength = 2147483647):
    cdef long i
    
    for i in range( poss.shape[0] ):
        if poss[i] < 0:
            poss[i] = 0
        else:
            break
        
    i = -1
    while True:
        if poss[i] > rlength:
            poss[i] = rlength
            i -= 1
        else:
            break
    return poss

cdef pileup_a_chromosome ( np.ndarray start_poss, np.ndarray end_poss, long l, float scale_factor = 1, float baseline_value = 0 ):
    """Return pileup of one chromosome.

    """
    cdef long i_s, i_e, p, pileup, pre_p, i
    cdef int a, b
    
    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append
    i_s = 0                         # index of start_poss
    i_e = 0                         # index of end_poss

    pileup = 0
    if start_poss.shape[0] == 0: return tmp
    pre_p = min(start_poss[0],end_poss[0])
    #print pre_p
    if pre_p != 0:
        # the first chunk of 0
        tmppadd( pre_p )
        tmpvadd( float_max(0,baseline_value) )
        
    pre_v = pileup
    
    while i_s < l and i_e < l:
        a = start_poss[i_s]
        b = end_poss[i_e]
        if a < b:
            p = a
            if p != pre_p:
                tmppadd( p )
                tmpvadd( float_max(pileup * scale_factor, baseline_value) )
                pre_p = p
            pileup += 1
            i_s += 1
        elif a > b:
            p = b
            if p != pre_p:
                tmppadd( p )
                tmpvadd( float_max(pileup * scale_factor, baseline_value) ) 
                pre_p = p
            pileup -= 1
            i_e += 1
        else:
            i_s += 1
            i_e += 1
    if i_e < l:
        # add rest of end positions
        for i in range(i_e, l):
            p = end_poss[i]
            #for p in end_poss[i_e:]:
            if p != pre_p:
                tmppadd( p )
                tmpvadd( float_max(pileup * scale_factor, baseline_value) )
                #ret.add_loc(chrom,pre_p,p,pileup)
                pre_p = p
            pileup -= 1
    if i_s < l:
        # add rest of start positions ( I don't think this will happen )
        raise Exception("start positions can't be the only things left!")

    return tmp

cdef max_over_two_pv_array ( tmparray1, tmparray2 ):
    cdef int pre_p, p1, p2
    cdef double v1, v2

    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append    

    (p1s,v1s) = tmparray1
    p1n = iter(p1s).next         # assign the next function to a viable to speed up
    v1n = iter(v1s).next

    (p2s,v2s) = tmparray2
    p2n = iter(p2s).next         # assign the next function to a viable to speed up
    v2n = iter(v2s).next

    pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
    try:
        p1 = p1n()
        v1 = v1n()

        p2 = p2n()
        v2 = v2n()

        while True:
            if p1 < p2:
                # clip a region from pre_p to p1, then set pre_p as p1.
                tmppadd( p1 )
                tmpvadd( float_max(v1,v2) )
                pre_p = p1
                # call for the next p1 and v1
                p1 = p1n()
                v1 = v1n()
            elif p2 < p1:
                # clip a region from pre_p to p2, then set pre_p as p2.
                tmppadd( p2 )
                tmpvadd( float_max(v1,v2) )
                pre_p = p2
                # call for the next p2 and v2
                p2 = p2n()
                v2 = v2n()
            elif p1 == p2:
                # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                tmppadd( p1 )
                tmpvadd( float_max(v1,v2) )
                pre_p = p1
                # call for the next p1, v1, p2, v2.
                p1 = p1n()
                v1 = v1n()
                p2 = p2n()
                v2 = v2n()
    except StopIteration:
        # meet the end of either bedGraphTrackI, simply exit
        pass
    return tmp


