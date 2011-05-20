# Time-stamp: <2011-05-19 23:20:16 Tao Liu>

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

def pileup_bdg (trackI, d, baseline_value = 0):
    """Pileup tags into bedGraphTrackI object with extension. Tag will
    be extended towards both sides with 1/2 of d.

    A tag is a single genomic location.

    trackI  : A FWTrackII object. For example, the shifted tags from PeakDetect object.
    d       : tag will be extended to 1/2 of this value to each side, i.e. a 'd' size fragment.
    baseline_value : a value to be filled for missing values.

    Return a bedGraphTrackI object.
    """
    #step = 10000000 + 2*d               # step to cache data points.

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    cdef int half_d1 = int(d//2)
    cdef int half_d2 = d - int(d//2)

    chrs = trackI.get_chr_names()       
    for chrom in chrs:
        tags = trackI.get_locations_by_chr(chrom)[0]
        l = len(tags)

        start_poss = array(BYTE4,[])
        end_poss   = array(BYTE4,[])

        for i in xrange(len(tags)):
            # shift to get start positions
            start_poss.append(max(0,tags[i]-half_d1)) # prevent coordinates < 0
            # shift to get end positions
            end_poss.append(max(0,tags[i]+half_d2))

        # combine start and end positions
        i_s = 0                         # index of start_poss
        i_e = 0                         # index of end_poss

        pileup = 0
        pre_p = 0
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
