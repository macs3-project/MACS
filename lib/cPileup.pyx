# Time-stamp: <2011-06-15 22:43:05 Tao Liu>

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
from array import array

from MACS2.IO.cFeatIO import bedGraphTrackI
from time import time
from MACS2.Constants import *
# ------------------------------------
# constants
# ------------------------------------
# to determine the byte size

def pileup_bdg (trackI, d, baseline_value = 0, directional=True, halfextension=True):
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
    #step = 10000000 + 2*d               # step to cache data points.

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    chrs = trackI.get_chr_names()       

    if directional:
        # only extend to 3' side
        if halfextension:
            five_shift = int(d*-0.25)  # five shift is used to move cursor towards 5' direction to find the start of fragment
            three_shift = int(d*0.75) # three shift is used to move cursor towards 3' direction to find the end of fragment
        else:
            five_shift = 0
            three_shift = d
    else:
        # both sides
        if halfextension:
            five_shift = int(d*0.25)
            three_shift = five_shift
        else:
            five_shift = int(d/2)
            three_shift = d - five_shift            

    for chrom in chrs:
        (plus_tags,minus_tags) = trackI.get_locations_by_chr(chrom)
        
        l = len(plus_tags)+len(minus_tags)

        start_poss = array(BYTE4,[])    # store all start positions
        end_poss   = array(BYTE4,[])    # store all end positions

        # for plus tags
        for i in xrange(len(plus_tags)):
            # shift to get start positions. To 5' side.
            #start_poss.append(max(0,plus_tags[i]-five_shift)) # prevent coordinates < 0
            start_poss.append(plus_tags[i]-five_shift) 
            # shift to get end positions by extending to d. To 3' side.
            #end_poss.append(max(0,plus_tags[i]+three_shift))
            end_poss.append(plus_tags[i]+three_shift)

        # for minus tags
        for i in xrange(len(minus_tags)):
            # shift to get start positions by extending to d. To 3' side.
            #start_poss.append(max(0,minus_tags[i]-three_shift)) # prevent coordinates < 0
            start_poss.append(minus_tags[i]-three_shift)
            # shift to get end positions. To 5' side.
            #end_poss.append(max(0,minus_tags[i]+five_shift))
            end_poss.append(minus_tags[i]+five_shift)
            
        # sort
        start_poss = sorted(start_poss)
        end_poss = sorted(end_poss)

        # Pileup by go through start positions and end positions,
        # while seeing start position, pileup ++
        # while seeing end position, pileup --
        #
        i_s = 0                         # index of start_poss
        i_e = 0                         # index of end_poss

        pileup = 0
        pre_p = min(start_poss[0],end_poss[0])
        while i_s < l and i_e < l:
            if start_poss[i_s] < end_poss[i_e]:
                p = start_poss[i_s]
                if p != pre_p:
                    ret.add_loc(chrom,pre_p,p,pileup)
                    pre_p = p
                pileup += 1
                i_s += 1
            elif start_poss[i_s] > end_poss[i_e]:
                p = end_poss[i_e]
                if p != pre_p:
                    ret.add_loc(chrom,pre_p,p,pileup)
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
                    ret.add_loc(chrom,pre_p,p,pileup)
                    pre_p = p
                pileup -= 1
        if i_s < l:
            # add rest of start positions ( I don't think this will happen )
            raise Exception("start positions can't be the only things left!")
            for p in start_poss[i_s:]:
                if p != pre_p:
                    ret.add_loc(chrom,pre_p,p,pileup)
                    pre_p = p
                pileup += 1

    return ret
