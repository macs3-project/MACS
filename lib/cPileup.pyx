# Time-stamp: <2011-05-18 14:20:32 Tao Liu>

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
    step = 10000000 + 2*d               # step to cache data points.

    ret = bedGraphTrackI(baseline_value=baseline_value) # bedGraphTrackI object to be returned.

    chrs = trackI.get_chr_names()       
    for chrom in chrs:
        tags = trackI.get_locations_by_chr(chrom)[0]
        l = len(tags)
        window_counts = array(BYTE4,[0]*step)
        startp = -1*d
        endp   = startp+step
        index_tag = 0
		
        while index_tag<l:
            s = tags[index_tag]-d/2     # start of tag
            e = s+d                     # end of tag
            
            if e < endp:
                # project tag to window_counts line
                ps = s-startp # projection start
                pe = ps+d     # projection end
                for i in xrange(ps,pe):
                    window_counts[i] += 1
                index_tag += 1
            else:
                # write it to zbdg file then reset parameters
                # keep this tag for next window
                prev = window_counts[d]
                left = startp+d
                right = left+1
                for i in xrange(d+1,step-d):
                    if window_counts[i] == prev:
                        # same value, extend
                        right += 1
                    else:
                        # diff value, close
                        if prev != 0:
                            ret.add_loc(chrom,left,right,prev)
                        prev = window_counts[i]
                        left = right
                        right = left + 1
                # last bin
                if prev != 0:                
                    ret.add_loc(chrom,left,right,prev)
                    
                # reset
                window_counts_next = array(BYTE4,[0]*step)
                # copy d values from the tail of previous window to next window
                for n,i in enumerate(xrange(step-2*d,step)): # debug
                    window_counts_next[n] = window_counts[i]
                window_counts = window_counts_next
                startp = endp - 2*d
                endp = startp+step
        # last window
        prev = window_counts[d]
        left = startp+d
        right = left+1
        for i in xrange(d+1,step-d):
            if window_counts[i] == prev:
                # same value, exrend
                right += 1
            else:
                # diff value, close
                if prev != 0:                
                    ret.add_loc(chrom,left,right,prev)
                prev = window_counts[i]
                left = right
                right = left + 1
        # last bin
        if prev != 0:        
            ret.add_loc(chrom,left,right,prev)
            
    return ret
