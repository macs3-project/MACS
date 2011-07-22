# Time-stamp: <2011-07-22 08:08:19 Tao Liu>

"""Module for Composite Score Track IO classes.

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import numpy as np
from numpy import int64,int32,float32

from libc.math cimport sqrt,log10

from MACS2.Constants import *
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO

# ------------------------------------
# constants
# ------------------------------------
__version__ = "scoreTrackI $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "scoreTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class compositeScoreTrackI:
    """Class for composite scoreGraph type data. Modified from
    bedGraphTrackI. The only difference is that we store pvalue score,
    qvalue score and foldchange together.

    In bedGraph, data are represented as continuous non-overlapping
    regions in the whole genome. I keep this assumption in all the
    functions. If data has overlaps, some functions will definitely
    give incorrect results.

    1. Continuous: the next region should be after the previous one
    unless they are on different chromosomes;
    
    2. Non-overlapping: the next region should never have overlaps
    with preceding region.

    The way to memorize bedGraph data is to remember the transition
    points together with values of their preceding regions. The last
    data point may exceed chromosome end, unless a chromosome
    dictionary is given. Remember the coordinations in bedGraph and
    this class is 0-indexed and right-open.
    
    """
    def __init__ (self):
        """Different with bedGraphTrackI, missing values are simply
        replaced with 0.
        
        """
        self.data = {}
        self.pointer = {}

    def add_chromosome ( self, chrom, chrom_max_len ):
        if not self.data.has_key(chrom):
            self.data[chrom] = np.zeros(chrom_max_len,dtype=[('pos','int64'),
                                                             ('sc1i1','float32'),
                                                             ('sc2i2','float32'),
                                                             ('sc1c2','float32'),
                                                             ('sc2c1','float32')])
            self.pointer[chrom] = 0

    def add (self,chromosome,endpos,sc1i1,sc2i2,sc1c2,sc2c1):
        """Add a chr-endpos-score-score-score-score block into data
        dictionary.

        """
        c = self.data[chromosome]
        i = self.pointer[chromosome]
        # get the preceding region
        c[i] = (endpos,sc1i1,sc2i2,sc1c2,sc2c1)
        self.pointer[chromosome] += 1

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if self.data.has_key(chromosome):
            return self.data[chromosome]
        else:
            return None

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.data.keys())
        return l

    def call_consistent (self, cutoff=50, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        Consistent peaks are those met all the following criteria

        1. sc1i1 >= cutoff
        2. sc2i2 >= cutoff
        3. sc1c2 <= cutoff
        4. sc2c1 <= cutoff
        """
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        #condition1_unique_peaks = PeakIO()                      # dictionary to save peaks
        #condition2_unique_peaks = PeakIO()                      # dictionary to save peaks        
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc1i1   = chrom_d[ 'sc1i1' ]
            chrom_sc2i2  = chrom_d[ 'sc2i2' ]
            chrom_sc1c2 = chrom_d[ 'sc1c2' ]
            chrom_sc2c1  = chrom_d[ 'sc2c1' ]

            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                vc1i1 = chrom_sc1i1[ x ]
                vc2i2 = chrom_sc2i2[ x ]
                vc1c2 = chrom_sc1c2[ x ]
                vc2c1 = chrom_sc2c1[ x ]                
                x += 1                  # index for the next point
                if vc1i1 >= cutoff and vc2i2 >= cutoff and vc1c2 <= cutoff and vc2c1 <= cutoff:
                    peak_content = [ ( pre_p, p, 0, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                v = chrom_score[ i ]
                if vc1i1 < cutoff or vc2i2 < cutoff or vc1c2 > cutoff or vc2c1 > cutoff:
                    pre_p = p
                    continue
                # for points met all criteria
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, 0, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = 0,
                                   peak_score  = 0,
                                   pileup      = 0,
                                   pscore      = 0,
                                   fold_change = 0,
                                   qscore      = 0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, 0, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = 0,
                           peak_score  = 0,
                           pileup      = 0,
                           pscore      = 0,
                           fold_change = 0,
                           qscore      = 0,
                           )
            
        return peaks

    def call_condition1_unique (self, cutoff=50, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        Condition 1 unique peaks are those met all the following criteria

        1. sc1i1 >= cutoff
        2. sc1c2 >= cutoff
        """
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc1i1   = chrom_d[ 'sc1i1' ]
            chrom_sc1c2 = chrom_d[ 'sc1c2' ]

            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                vc1i1 = chrom_sc1i1[ x ]
                vc1c2 = chrom_sc1c2[ x ]
                x += 1                  # index for the next point
                if vc1i1 >= cutoff and vc1c2 >= cutoff:
                    peak_content = [ ( pre_p, p, 0, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                v = chrom_score[ i ]
                if vc1i1 < cutoff or vc1c2 < cutoff:
                    pre_p = p
                    continue
                # for points met all criteria
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, 0, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = 0,
                                   peak_score  = 0,
                                   pileup      = 0,
                                   pscore      = 0,
                                   fold_change = 0,
                                   qscore      = 0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, 0, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = 0,
                           peak_score  = 0,
                           pileup      = 0,
                           pscore      = 0,
                           fold_change = 0,
                           qscore      = 0,
                           )
            
        return peaks

    def call_condition2_unique (self, cutoff=50, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        Condition 2 unique peaks are those met all the following criteria

        1. sc2i2 >= cutoff
        2. sc2c1 >= cutoff
        """
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc2i2   = chrom_d[ 'sc2i2' ]
            chrom_sc2c1 = chrom_d[ 'sc2c1' ]

            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                vc2i2 = chrom_sc2i2[ x ]
                vc2c1 = chrom_sc2c1[ x ]
                x += 1                  # index for the next point
                if vc2i2 >= cutoff and vc2c1 >= cutoff:
                    peak_content = [ ( pre_p, p, 0, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                v = chrom_score[ i ]
                if vc2i2 < cutoff or vc2c1 < cutoff:
                    pre_p = p
                    continue
                # for points met all criteria
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, 0, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = 0,
                                   peak_score  = 0,
                                   pileup      = 0,
                                   pscore      = 0,
                                   fold_change = 0,
                                   qscore      = 0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, 0, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = 0,
                           peak_score  = 0,
                           pileup      = 0,
                           pscore      = 0,
                           fold_change = 0,
                           qscore      = 0,
                           )
        return peaks

    def total ( self ):
        """Return the number of regions in this object.

        """
        t = 0
        for chrom in self.data.keys():
            t += self.pointer[chrom]
        return t


def make_compositeScoreTrack (bdgTrack1, bdgTrack2, bdgTrack3, bdgTrack4 ):
    """A modified overlie function for MACS DIFF.
    
    """
    assert isinstance(bdgTrack1,bedGraphTrackI), "bdgTrack1 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack3,bedGraphTrackI), "bdgTrack3 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack4,bedGraphTrackI), "bdgTrack4 is not a bedGraphTrackI object"    
    
    ret = compositeScoreTrackI()
    retadd = ret.add

    chr1 = set(bdgTrack1.get_chr_names())
    chr2 = set(bdgTrack2.get_chr_names())
    chr3 = set(bdgTrack1.get_chr_names())
    chr4 = set(bdgTrack2.get_chr_names())    
    
    common_chr = chr1.intersection(chr2).intersection(chr3).intersection(chr4)
    for chrom in common_chr:
            
        (p1s,v1s) = bdgTrack1.get_data_by_chr(chrom) # arrays for position and values
        p1n = iter(p1s).next         # assign the next function to a viable to speed up
        v1n = iter(v1s).next

        (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
        p2n = iter(p2s).next         # assign the next function to a viable to speed up
        v2n = iter(v2s).next

        (p3s,v3s) = bdgTrack3.get_data_by_chr(chrom) # arrays for position and values
        p1n = iter(p1s).next         # assign the next function to a viable to speed up
        v1n = iter(v1s).next

        (p4s,v4s) = bdgTrack4.get_data_by_chr(chrom) # arrays for position and values
        p2n = iter(p4s).next         # assign the next function to a viable to speed up
        v2n = iter(v4s).next

        chrom_max_len = len(p1s)+len(p2s)+len(p3s)+len(p4s) # this is the maximum number of locations needed to be recorded in scoreTrackI for this chromosome.
            
        ret.add_chromosome(chrom,chrom_max_len)

        pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
        try:
            p1 = p1n()
            v1 = v1n()

            p2 = p2n()
            v2 = v2n()

            p3 = p3n()
            v3 = v3n()

            p4 = p4n()
            v4 = v4n()
            
            while True:
                min_p = min( p1, p2, p3, p4 )
                retadd( chrom, p1, v1, v2, v3, v4 )
                pre_p = min_p

                if p1 == min_p:
                    p1 = p1n()
                    v1 = v1n()
                if p2 == min_p:
                    p2 = p2n()
                    v2 = v2n()
                if p3 == min_p:
                    p3 = p3n()
                    v3 = v3n()
                if p4 == min_p:
                    p4 = p4n()
                    v4 = v4n()                                        
        except StopIteration:
            # meet the end of either bedGraphTrackI, simply exit
            pass
        
    return ret
