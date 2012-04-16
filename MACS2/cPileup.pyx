# Time-stamp: <2012-04-12 11:06:08 Tao Liu>

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

from cpython cimport bool

# ------------------------------------
# constants
# ------------------------------------
# to determine the byte size

def pileup_bdg (trackI, int d, int baseline_value = 0, bool directional = True, bool halfextension = True, float scale_factor = 1):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards 3' side with size of d if directional is Ture,
    or both sides with d/2 if directional is False.

    A tag is a single genomic location.

    trackI  : A FWTrackII object with raw plus and minus 5' end positions
    d       : tag will be extended to this value to 3' direction, unless directional is False.
    baseline_value : a value to be filled for missing values.
    directional: if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
    halfextension: only make a fragment of d/2 size centered at fragment center

    Return a bedGraphTrackI object.
    """
    cdef long five_shift, three_shift, l, i, j, i_s, i_e, p, pileup, pre_p

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
        (plus_tags,minus_tags) = trackI.get_locations_by_chr(chrom)

        l = len(plus_tags)+len(minus_tags)        
        
        ( start_poss, end_poss ) = start_and_end_poss( plus_tags, minus_tags, five_shift, three_shift )

        ret.add_a_chromosome( chrom, pileup_a_chromosome ( start_poss, end_poss, l, scale_factor ) )

        # free mem?
        start_poss = None
        end_poss = None        

    return ret

def pileup_w_multiple_d_bdg ( trackI, d_s, int baseline_value = 0, bool directional = True, bool halfextension = True, scale_factor_s = [] ):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards 3' side with size of d if directional is Ture,
    or both sides with d/2 if directional is False.

    A tag is a single genomic location.

    trackI  : A FWTrackII object with raw plus and minus 5' end positions
    d       : tag will be extended to this value to 3' direction, unless directional is False.
    baseline_value : a value to be filled for missing values.
    directional: if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
    halfextension: only make a fragment of d/2 size centered at fragment center

    Return a bedGraphTrackI object.
    """
    cdef long d, five_shift, three_shift, l, i, j, i_s, i_e, p, pileup, pre_p
    cdef float scale_factor

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
        (plus_tags,minus_tags) = trackI.get_locations_by_chr(chrom)

        l = len(plus_tags)+len(minus_tags)

        prev_pileup = None

        for i in range(len(d_s)):
            
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]

            ( start_poss, end_poss ) = start_and_end_poss( plus_tags, minus_tags, five_shift, three_shift )

            tmp_pileup = pileup_a_chromosome ( start_poss, end_poss, l, scale_factor )
            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        ret.add_a_chromosome( chrom, prev_pileup )
        start_poss = None
        end_poss = None        

    return ret

cdef start_and_end_poss ( plus_tags, minus_tags, long five_shift, long three_shift ):
    cdef long i
    
    start_poss = array(BYTE4,[])    # store all start positions
    end_poss   = array(BYTE4,[])    # store all end positions

    # for plus tags
    for i in xrange(len(plus_tags)):
        # shift to get start positions. To 5' side.
        start_poss.append(plus_tags[i]-five_shift) 
        # shift to get end positions by extending to d. To 3' side.
        end_poss.append(plus_tags[i]+three_shift)

    # for minus tags
    for i in xrange(len(minus_tags)):
        # shift to get start positions by extending to d. To 3' side.
        start_poss.append(minus_tags[i]-three_shift)
        # shift to get end positions. To 5' side.
        end_poss.append(minus_tags[i]+five_shift)
            
    # sort
    start_poss = sorted(start_poss)
    end_poss = sorted(end_poss)

    return ( start_poss, end_poss )

cdef pileup_a_chromosome ( start_poss, end_poss, long l, float scale_factor = 1 ):
    """Return pileup of one chromosome.

    """
    cdef long i_s, i_e, p, pileup, pre_p
    
    tmp = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
    tmppadd = tmp[0].append
    tmpvadd = tmp[1].append    
    i_s = 0                         # index of start_poss
    i_e = 0                         # index of end_poss

    pileup = 0
    pre_p = min(start_poss[0],end_poss[0])
    if pre_p != 0:
        # the first chunk of 0
        tmppadd( pre_p )
        tmpvadd( 0 )
        
    pre_v = pileup
    
    while i_s < l and i_e < l:
        if start_poss[i_s] < end_poss[i_e]:
            p = start_poss[i_s]
            if p != pre_p:
                tmppadd( p )
                tmpvadd( pileup * scale_factor )
                #ret.add_loc(chrom,pre_p,p,pileup)
                pre_p = p
            pileup += 1
            i_s += 1
        elif start_poss[i_s] > end_poss[i_e]:
            p = end_poss[i_e]
            if p != pre_p:
                tmppadd( p )
                tmpvadd( pileup * scale_factor ) 
                #ret.add_loc(chrom,pre_p,p,pileup)
                pre_p = p
            pileup -= 1
            i_e += 1
        else:
            i_s += 1
            i_e += 1
    if i_e < l:
        # add rest of end positions
        for p in end_poss[i_e:]:
            if p != pre_p:
                tmppadd( p )
                tmpvadd( pileup * scale_factor )
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
                tmpvadd( max(v1,v2) )
                pre_p = p1
                # call for the next p1 and v1
                p1 = p1n()
                v1 = v1n()
            elif p2 < p1:
                # clip a region from pre_p to p2, then set pre_p as p2.
                tmppadd( p2 )
                tmpvadd( max(v1,v2) )
                pre_p = p2
                # call for the next p2 and v2
                p2 = p2n()
                v2 = v2n()
            elif p1 == p2:
                # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                tmppadd( p1 )
                tmpvadd( max(v1,v2) )
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
