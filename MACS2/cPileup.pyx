# Time-stamp: <2013-04-09 15:32:25 Tao Liu>

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

from MACS2.IO.cFixWidthTrack import FWTrackIII
from MACS2.IO.cPairedEndTrack import PETrackI
from MACS2.IO.cBedGraph import bedGraphTrackI
from MACS2.Constants import *

import numpy as np
cimport numpy as np

from numpy cimport int32_t

from cpython cimport bool

# ------------------------------------
# functions
# ------------------------------------

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline long long_max(long a, long b): return a if a >= b else b
cdef inline float float_max(float a, float b): return a if a >= b else b

# Unified pileup function #
cpdef unified_pileup_bdg(track,
                         ds,
                         scale_factors,
                         float baseline_value = 0.0,
                         bint directional = True, 
                         bint halfextension = True):
    """This function will call corresponding function for FWTrackIII
    or PETrackI to pileup fragments.

    It will call pileup_w_multiple_d* functions for control input, then
    calculate maximum values; or call normal pileup functions for
    treatment data.

    """

    chrs = track.get_chr_names()
    if type(ds) is list:
        # multiple pileup (e.g. control with 1k, 10k and d extension)
        if isinstance(track, FWTrackIII):
            return pileup_w_multiple_d_bdg(track, ds, scale_factors,
                                           baseline_value,
                                           directional, halfextension)
        elif isinstance(track, PETrackI):
            return pileup_w_multiple_d_bdg_pe(track, ds, scale_factors,
                                              baseline_value)
        else:
            raise ValueError("track must be of type FWTrackIII or PETrackI")
    else:
        # single extension (e.g. treatment data)
        if isinstance(track, FWTrackIII):
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
            raise ValueError("track must be of type FWTrackIII or PETrackI")

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

    trackI  : A FWTrackIII object with raw plus and minus 5' end positions
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

    trackI  : A FWTrackIII object with raw plus and minus 5' end positions
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

            #print chrom, baseline_value
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
    cdef long i
    
    for i in range( poss.shape[0] ):
        if poss[i] < 0:
            poss[i] = 0
        else:
            break
        
    for i in range( poss.shape[0]-1, -1, -1 ):
        if poss[i] > rlength:
            poss[i] = rlength
        else:
            break
    
    return poss

# general pileup function
cpdef se_all_in_one_pileup ( np.ndarray plus_tags, np.ndarray minus_tags, long five_shift, long three_shift, int rlength, float scale_factor, float baseline_value ):
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
        long i_s, i_e, i
        int a, b, p, pre_p, pileup
        list tmp
        long lp = plus_tags.shape[0]
        long lm = minus_tags.shape[0]
        long l = lp + lm
        np.ndarray start_poss, end_poss

    start_poss = np.concatenate( ( plus_tags-five_shift, minus_tags-three_shift ) )
    end_poss   = np.concatenate( ( plus_tags+three_shift, minus_tags+five_shift ) )    

    # sort
    start_poss.sort()
    end_poss.sort()

    # fix negative coordinations and those extends over end of chromosomes
    start_poss = fix_coordinates(start_poss, rlength)
    end_poss = fix_coordinates(end_poss, rlength)

    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append
    i_s = 0                         # index of start_poss
    i_e = 0                         # index of end_poss

    pileup = 0
    if start_poss.shape[0] == 0: return tmp
    pre_p = min(start_poss[0],end_poss[0])

    if pre_p != 0:
        # the first chunk of 0
        tmppadd( pre_p )
        tmpvadd( float_max(0,baseline_value) )
        
    pre_v = pileup

    a = start_poss[i_s]
    b = end_poss[i_e]

    try:
        while 1:
            if a < b:
                p = a
                if p != pre_p:
                    tmppadd( p )
                    tmpvadd( float_max(pileup * scale_factor, baseline_value) )
                    pre_p = p
                pileup += 1
                i_s += 1
                a = start_poss[i_s]            
            elif a > b:
                p = b
                if p != pre_p:
                    tmppadd( p )
                    tmpvadd( float_max(pileup * scale_factor, baseline_value) ) 
                    pre_p = p
                pileup -= 1
                i_e += 1
                b = end_poss[i_e]
            else:
                i_s += 1
                i_e += 1
                a = start_poss[i_s]
                b = end_poss[i_e]
    except IndexError:
        pass

    if i_e < l:
        # add rest of end positions
        for i in range(i_e, l):
            p = end_poss[i]
            if p != pre_p:
                tmppadd( p )
                tmpvadd( float_max(pileup * scale_factor, baseline_value) )
                pre_p = p
            pileup -= 1

    # clean mem
    start_poss.resize(100000, refcheck=False)
    start_poss.resize(0, refcheck=False)
    end_poss.resize(100000, refcheck=False)
    end_poss.resize(0, refcheck=False)                            
    
    return tmp


cpdef se_all_in_one_pileup2 ( np.ndarray plus_tags, np.ndarray minus_tags, long five_shift, long three_shift, int rlength, float scale_factor, float baseline_value ):
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
        long i_s, i_e, i, j
        int a, b, p, pre_p, pileup
        list tmp
        long lp = plus_tags.shape[0]
        long lm = minus_tags.shape[0]
        np.ndarray start_poss, end_poss, alldata, parray, varray

    alldata = np.zeros( 2*(lp+lm), dtype = [("pos", "int32"),("v","int8")] )
    alldata["pos"][:lp] = plus_tags-five_shift
    alldata["pos"][lp:lp+lm] = minus_tags-three_shift
    alldata["v"][:lp+lm] = 1

    alldata["pos"][lp+lm:lp+lp+lm] = plus_tags+three_shift
    alldata["pos"][lp+lp+lm:] = minus_tags+five_shift
    alldata["v"][lp+lm:] = -1

    # sort
    alldata.sort( order = "pos" )

    #print alldata

    # fix coordinates
    for i in range( alldata.shape[0] ):
        if alldata[i][0] < 0:
            alldata[i][0] = 0
        else:
            break
        
    for i in range( alldata.shape[0]-1, -1, -1 ):
        if alldata[i][0] > rlength:
            alldata[i][0] = rlength
        else:
            break
    
    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append
    i_s = 0                         # index of start_poss
    i_e = 0                         # index of end_poss

    parray = alldata["pos"]
    varray = alldata["v"]

    pileup = 0
    if alldata.shape[0] == 0: return tmp

    pre_p = parray[0]

    if pre_p != 0:
        # the first chunk of 0
        tmppadd( pre_p )
        tmpvadd( float_max(0,baseline_value) )

    for i in range( alldata.shape[0] ):
        p = parray[ i ]
        if p != pre_p:
            tmppadd( p )
            tmpvadd( float_max(pileup * scale_factor, baseline_value) )
            pre_p = p
        pileup += varray[ i ]

    alldata.resize(100000, refcheck=False)
    alldata.resize(0, refcheck=False)
    
    return tmp


cpdef quick_pileup ( np.ndarray plus_tags, np.ndarray minus_tags, float scale_factor, float baseline_value ):
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
        long i_s, i_e, i
        int a, b, p, pre_p, pileup
        list tmp
        long lp = plus_tags.shape[0]
        long lm = minus_tags.shape[0]
        long l = lp + lm

    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append
    i_s = 0                         # index of plus_tags
    i_e = 0                         # index of minus_tags

    pileup = 0
    if plus_tags.shape[0] == 0: return tmp
    pre_p = min(plus_tags[0],minus_tags[0])
    #print pre_p
    if pre_p != 0:
        # the first chunk of 0
        tmppadd( pre_p )
        tmpvadd( float_max(0,baseline_value) )
        
    pre_v = pileup
    
    while i_s < lp and i_e < lm:
        a = plus_tags[i_s]
        b = minus_tags[i_e]
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
    if i_e < lm:
        # add rest of end positions
        for i in range(i_e, lm):
            p = minus_tags[i]
            #for p in minus_tags[i_e:]:
            if p != pre_p:
                tmppadd( p )
                tmpvadd( float_max(pileup * scale_factor, baseline_value) )
                #ret.add_loc(chrom,pre_p,p,pileup)
                pre_p = p
            pileup -= 1
    #    if i_s < l:
    #        # add rest of start positions ( I don't think this will happen )
    #        raise Exception("start positions can't be the only things left!")
    return tmp

# general function to calculate maximum between two arrays.

cpdef max_over_two_pv_array ( tmparray1, tmparray2 ):
    """Merge two position-value arrays. For intersection regions, take
    the maximum value within region.

    tmparray1 and tmparray2 are [p,v] type lists, same as the output
    from quick_pileup function.
    """
    
    cdef:
        int pre_p, p1, p2
        float v1, v2
        list tmp
        long l1, l2, l

    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append    

    (p1s,v1s) = tmparray1
    p1n = iter(p1s).next         # assign the next function to a viable to speed up
    v1n = iter(v1s).next
    l1  = len(p1s)

    (p2s,v2s) = tmparray2
    p2n = iter(p2s).next         # assign the next function to a viable to speed up
    v2n = iter(v2s).next
    l2  = len(p2s)
    
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


