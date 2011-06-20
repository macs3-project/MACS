# Time-stamp: <2011-06-19 17:09:12 Tao Liu>

"""Module for Feature IO classes.

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
import logging

from array import array

import numpy as np

from libc.math cimport sqrt

from MACS2.Constants import *
from MACS2.cProb import poisson_cdf
from MACS2.IO.cScoreTrack import scoreTrackI

# ------------------------------------
# constants
# ------------------------------------
__version__ = "BedGraph $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "bedGraphTrackI class"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
class bedGraphTrackI:
    """Class for bedGraph type data.

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
    def __init__ (self, baseline_value=0):
        """
        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with baseline_value.
        
        """
        self.__data = {}
        self.maxvalue =-10000
        self.minvalue = 10000
        self.baseline_value = baseline_value

    def add_loc (self,chromosome,startpos,endpos,value):
        """Add a chr-start-end-value block into __data dictionary.

        """
        # basic assumption, end pos should > start pos
        assert endpos > startpos, "endpos %d can't be smaller than start pos %d" % (endpos,startpos)

        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0
        
        if not self.__data.has_key(chromosome):
            self.__data[chromosome] = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
            c = self.__data[chromosome]
            if startpos:
                # start pos is not 0, then add two blocks, the first
                # with "baseline_value"; the second with "value"
                c[0].append(startpos)
                c[1].append(self.baseline_value)
            c[0].append(endpos)
            c[1].append(value)
        else:
            c = self.__data[chromosome]            
            # get the preceding region
            pre_pos = c[0][-1]
            pre_v   = c[1][-1]
            # to check 1. continuity; 2. non-overlapping
            assert pre_pos < endpos , "bedGraph regions are not continuous."
            assert pre_pos <= startpos , "bedGraph regions have overlappings."
            
            if startpos != pre_pos:
                # there is a gap, so fill it with baseline_value
                c[0].append(startpos)
                c[1].append(self.baseline_value)
                # then add this region
                c[0].append(endpos)
                c[1].append(value)
            else:
                # if this region is next to the previous one.
                if pre_v == value:
                    # if value is the same, simply extend it.
                    c[0][-1] = endpos
                else:
                    # otherwise, add a new region
                    c[0].append(endpos)
                    c[1].append(value)

        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if self.__data.has_key(chromosome):
            return self.__data[chromosome]
        else:
            return None

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.__data.keys())
        return l

    def write_bedGraph (self, fhd, name, description):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        shift will be used to shift the coordinates. default: 0
        """
        #fhd.write("track type=bedGraph name=\"%s\" description=\"%s\"\n" % (name,description))
        chrs = self.get_chr_names()
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).next
            vnext = iter(v).next
            pre = 0
            for i in xrange(len(p)):
                pos = pnext()
                value = vnext()
                #if value != self.baseline_value:
                # never write baseline_value
                fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,pos,value))
                pre = pos

    def reset_baseline (self, baseline_value):
        """Reset baseline value to baseline_value.

        So any region between self.baseline_value and baseline_value
        will be set to baseline_value.
        
        """
        self.baseline_value = baseline_value
        self.filter_score(cutoff=baseline_value)
        self.merge_regions()

    def merge_regions (self):
        """Merge nearby regions with the same value.
        
        """
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).next
            vnext = iter(v).next

            # new arrays
            new_pos = array(BYTE4,[pnext(),])
            new_value = array(FBYTE4,[vnext(),])

            newpa = new_pos.append
            newva = new_value.append
            
            new_pre_pos = new_pos[0]
            new_pre_value = new_value[0]

            for i in xrange(1,len(p)):
                pos = pnext()
                value = vnext()
                if value == new_pre_value:
                    new_pos[-1] = pos
                else:
                    # add new region
                    newpa(pos)
                    newva(value)
                    new_pre_pos = pos
                    new_pre_value = value
            self.__data[chrom] = [new_pos,new_value]
        return True
                
    def filter_score (self, cutoff=0):
        """Filter using a score cutoff. Any region lower than score
        cutoff will be set to self.baseline_value.

        Self will be modified.
        """
        ret = bedGraphTrackI()
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).next
            vnext = iter(v).next

            # new arrays
            new_pos = array(BYTE4,[])
            new_value = array(FBYTE4,[])
            new_pre_pos = 0
            new_pre_value = None

            for i in xrange(len(p)):
                pos = pnext()
                value = vnext()

                if value < cutoff:
                    # this region will be set to baseline_value
                    if new_pre_value == self.baseline_value:
                        # if preceding region is at baseline, extend it
                        new_pos[-1] = pos
                    else:
                        # else add a new baseline region
                        new_pos.append(pos)
                        new_value.append(self.baseline_value)
                else:
                    # put it into new arrays
                    new_pos.append(pos)
                    new_value.append(value)
                new_pre_pos = new_pos[-1]
                new_pre_value = new_value[-1]
            self.__data[chrom]=[new_pos,new_value]
        return True

    def summary (self):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).
        
        """
        n_v = 0
        sum_v = 0
        max_v = -100000
        min_v = 100000
        for (p,v) in self.__data.values():
            # for each chromosome
            pre_p = 0
            for i in range(len(p)):
                # for each region
                l = p[i]-pre_p
                sum_v += v[i]*l
                n_v += l
                pre_p = p[i]
            max_v = max(max(v),max_v)
            min_v = min(min(v),min_v)
        mean_v = float(sum_v)/n_v
        variance = 0.0
        for (p,v) in self.__data.values():
            for i in range(len(p)):
                # for each region
                tmp = v[i]-mean_v
                l = p[i]-pre_p
                variance += tmp*tmp*l
                pre_p = p[i]

        variance /= float(n_v-1)
        std_v = sqrt(variance)
        return (sum_v, max_v, min_v, mean_v, std_v)

    def call_peaks (self, cutoff=1, up_limit=1e310, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 1.
        up_limit: the highest acceptable value. Default 10^{310}
          * so only allow peak with value >=cutoff and <=up_limit
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        """
        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            (ps,vs) = self.get_data_by_chr(chrom) # arrays for position and values
            psn = iter(ps).next         # assign the next function to a viable to speed up
            vsn = iter(vs).next
            x = 0
            pre_p = 0                   # remember previous position
            while True:
                # find the first region above cutoff
                try:                    # try to read the first data range for this chrom
                    p = psn()
                    v = vsn()
                except:
                    break
                x += 1                  # index for the next point
                if v >= cutoff and v <= up_limit:
                    peak_content = [(pre_p,p,v),]
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in range(x,len(ps)):
                # continue scan the rest regions
                p = psn()
                v = vsn()
                if v < cutoff or v > up_limit: # not be detected as 'peak'
                    pre_p = p
                    continue
                # for points above cutoff
                # if the gap is allowed
                if pre_p - peak_content[-1][1] <= max_gap:
                    peak_content.append((pre_p,p,v))
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[-1][1]-peak_content[0][0]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        tsummit = []
                        summit = None
                        summit_value = None
                        for (tstart,tend,tvalue) in peak_content:
                            if not summit_value or summit_value < tvalue:
                                tsummit = [int((tend+tstart)/2),]
                                summit_value = tvalue
                            elif summit_value == tvalue:
                                tsummit.append( int((tend+tstart)/2) )
                        summit = tsummit[int((len(tsummit)+1)/2)-1 ]
                        peaks.add(chrom,peak_content[0][0],peak_content[-1][1],
                                  summit=summit,peak_height=summit_value)
                    # start a new peak
                    peak_content = [(pre_p,p,v),]
                pre_p = p
                
            # save the last peak
            if peak_length >= min_length: # if the peak is too small, reject it
                summit = None
                summit_value = None
                for (tstart,tend,tvalue) in peak_content:
                    if not summit_value or summit_value < tvalue:
                        summit = int((tend+tstart)/2)
                        summit_value = tvalue
                peaks.add(chrom,peak_content[0][0],peak_content[-1][1],
                          summit=summit-peak_content[0][0],peak_height=summit_value)
            
        return peaks

    def total (self):
        """Return the number of regions in this object.

        """
        t = 0
        for (p,s) in self.__data.values():
            t += len(p)
        return t


    def overlie (self, bdgTrack2, func=max ):
        """Calculate two bedGraphTrackI objects by letting self
        overlying bdgTrack2, with user-defined functions.

        Transition positions from both bedGraphTrackI objects will be
        considered and combined. For example:

           #1 bedGraph (self)   |      #2 bedGraph
        -----------------------------------------------
        chr1  0    100  0       | chr1    0    150  1
        chr1  100  200  3       | chr1    150  250  2
        chr1  200  300  4       | chr1    250  300  4

        these two bedGraphs will be combined to have five transition
        points: 100, 150, 200, 250, and 300. So in order to calculate
        two bedGraphs, I pair values within the following regions
        like:

        chr   s   e     (#1,#2)   applied_func_max
        -----------------------------------------------
        chr1  0   100   (0,1)     1
        chr1  100 150   (3,1)     3
        chr1  150 200   (3,2)     3
        chr1  200 250   (4,2)     4
        chr1  250 300   (4,4)     4

        Then the given 'func' will be applied on each 2-tuple as func(#1,#2)


        Return value is a bedGraphTrackI object.
        """
        assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = bedGraphTrackI()
        retadd = ret.add_loc
        
        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for chrom in common_chr:
            (p1s,v1s) = self.get_data_by_chr(chrom) # arrays for position and values
            p1n = iter(p1s).next         # assign the next function to a viable to speed up
            v1n = iter(v1s).next

            (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
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
                        retadd(chrom,pre_p,p1,func(v1,v2))
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        retadd(chrom,pre_p,p2,func(v1,v2))
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        retadd(chrom,pre_p,p1,func(v1,v2))
                        pre_p = p1
                        # call for the next p1, v1, p2, v2.
                        p1 = p1n()
                        v1 = v1n()
                        p2 = p2n()
                        v2 = v2n()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass
        
        ret.merge_regions()
        return ret
                       
    def apply_func ( self, func ):
        """Apply function 'func' to every value in this bedGraphTrackI object.

        *Two adjacent regions with same value after applying func will
        not be merged.
        """
        t = 0
        for (p,s) in self.__data.values():
            for i in xrange(len(s)):
                s[i] = func(s[i])
        self.maxvalue = func(self.maxvalue)
        self.minvalue = func(self.minvalue)
        return True

    def make_scoreTrack_for_macs (self, bdgTrack2 ):
        """A modified overlie function for MACS v2.

        Return value (-1000*log10pvalue as integer, a trick) is a bedGraphTrackI object.
        """
        assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = scoreTrackI()
        retadd = ret.add
        
        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for chrom in common_chr:
            
            (p1s,v1s) = self.get_data_by_chr(chrom) # arrays for position and values
            p1n = iter(p1s).next         # assign the next function to a viable to speed up
            v1n = iter(v1s).next

            (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
            p2n = iter(p2s).next         # assign the next function to a viable to speed up
            v2n = iter(v2s).next

            chrom_max_len = len(p1s)+len(p2s) # this is the maximum number of locations needed to be recorded in scoreTrackI for this chromosome.
            
            ret.add_chromosome(chrom,chrom_max_len)

            pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
            try:
                p1 = p1n()
                v1 = v1n()

                p2 = p2n()
                v2 = v2n()

                while True:
                    if p1 < p2:
                        # clip a region from pre_p to p1, then set pre_p as p1.
                        retadd( chrom, p1, v1, v2 )
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        retadd( chrom, p2, v1, v2 )
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        retadd( chrom, p1, v1, v2 )
                        pre_p = p1
                        # call for the next p1, v1, p2, v2.
                        p1 = p1n()
                        v1 = v1n()
                        p2 = p2n()
                        v2 = v2n()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass
        
        #ret.merge_regions()
        return ret
