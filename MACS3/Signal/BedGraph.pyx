# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-11-24 17:14:17 Tao Liu>

"""Module for BedGraph data class.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import logging
from array import array

# ------------------------------------
# MACS3 modules
# ------------------------------------

from MACS3.Utilities.Constants import *
from MACS3.Signal.ScoreTrack import scoreTrackII,CombinedTwoTrack
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO
from MACS3.Signal.Prob import chisq_logp_e

# ------------------------------------
# Other modules
# ------------------------------------

import numpy as np
from cpython cimport bool

# ------------------------------------
# C lib
# ------------------------------------

from libc.math cimport sqrt, log, log1p, exp, log10

# ------------------------------------
# constants
# ------------------------------------
__version__ = "BedGraph $Revision$"
__author__ = "Tao Liu <vladimir.liu@gmail.com>"
__doc__ = "bedGraphTrackI class"

# ------------------------------------
# Misc functions
# ------------------------------------
LOG10_E = 0.43429448190325176

cdef inline float fisher_method_combining_two_log10pvalues ( float p1, float p2 ):
    return ( p1 + p2 ) - log1p( ( p1 + p2 ) / LOG10_E ) * LOG10_E

cdef inline float mean ( float p1, float p2 ):
    return ( p1 + p2 ) / 2

# ------------------------------------
# Classes
# ------------------------------------
cdef class bedGraphTrackI:
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
    cdef:
        dict __data
        double maxvalue
        double minvalue
        double baseline_value

    def __init__ (self, double baseline_value=0):
        """
        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with baseline_value.

        """
        self.__data = {}
        self.maxvalue = -10000000 # initial maximum value is tiny since I want safe_add_loc to update it
        self.minvalue = 10000000  # initial minimum value is large since I want safe_add_loc to update it
        self.baseline_value = baseline_value

    def add_a_chromosome ( self, bytes chrom, d ):
        """Unsafe method. Only to be used by cPileup.pyx.
        """
        self.__data[chrom] = d

    cpdef add_loc ( self, bytes chromosome, int startpos, int endpos, float value ):
        """Add a chr-start-end-value block into __data dictionary.

        Difference between safe_add_loc: no check, but faster. Save
        time while being called purely within MACS, so that regions
        are continuous without gaps.

        """
        cdef float pre_v
        # basic assumption, end pos should > start pos

        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0

        if chromosome not in self.__data:
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
            pre_v   = c[1][-1]

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


    def destroy ( self ):
        """ destroy content, free memory.
        """
        cdef:
            set chrs
            bytes chromosome

        chrs = self.get_chr_names()
        for chromosome in chrs:
            if chromosome in self.__data:
                self.__data[chromosome] = [None, None]
                self.__data.pop(chromosome)
        return True

    def safe_add_loc ( self, bytes chromosome, int startpos, int endpos, double value):
        """Add a chr-start-end-value block into __data dictionary.

        """
        # basic assumption, end pos should > start pos
        assert endpos > startpos, "endpos %d can't be smaller than start pos %d" % (endpos,startpos)

        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0

        if chromosome not in self.__data:
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

    def get_data_by_chr (self, bytes chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if chromosome in self.__data:
            return self.__data[chromosome]
        else:
            return None

    cpdef set get_chr_names (self):
        """Return all the chromosome names stored.

        """
        return set(sorted(self.__data.keys()))

    def write_bedGraph (self, fhd, str name, str description, bool trackline=True):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        shift will be used to shift the coordinates. default: 0
        """
        cdef:
            int pre, pos, i
            double value
            bytes chrom
            set chrs

        if trackline:
            trackcontents = (name.replace("\"", "\\\""), description.replace("\"", "\\\""))
            fhd.write("track type=bedGraph name=\"%s\" description=\"%s\" visibility=2 alwaysZero=on\n" % trackcontents)
        chrs = self.get_chr_names()
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).__next__
            vnext = iter(v).__next__
            pre = 0

            for i in range(len(p)):
                pos = pnext()
                value = vnext()
                #if value != self.baseline_value:
                # never write baseline_value
                fhd.write("%s\t%d\t%d\t%.5f\n" % (chrom.decode(),pre,pos,value))
                pre = pos

    def reset_baseline (self, double baseline_value):
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
        cdef:
            int new_pre_pos, pos, i
            double new_pre_value, value
            bytes chrom
            set chrs

        chrs = self.get_chr_names()
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).__next__
            vnext = iter(v).__next__

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

    def filter_score (self, double cutoff=0):
        """Filter using a score cutoff. Any region lower than score
        cutoff will be set to self.baseline_value.

        Self will be modified.
        """
        cdef:
            int new_pre_pos, pos, i
            double new_pre_value, value
            bytes chrom
            set chrs

        chrs = self.get_chr_names()
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).__next__
            vnext = iter(v).__next__

            # new arrays
            new_pos = array(BYTE4,[])
            new_value = array(FBYTE4,[])
            new_pre_pos = 0
            new_pre_value = 0

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
        cdef:
            long n_v
            double sum_v, max_v, min_v, mean_v, variance, tmp, std_v
            int pre_p, l, i

        pre_p = 0
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
        mean_v = sum_v/n_v
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
        return (sum_v, n_v, max_v, min_v, mean_v, std_v)

    def call_peaks (self, double cutoff=1, double up_limit=1e310, int min_length=200, int max_gap=50,
                    bool call_summits=False):
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
        cdef:
            int peak_length, x, pre_p, p, i, summit, tstart, tend
            double v, summit_value, tvalue
            bytes chrom
            set chrs

        #if call_summits: close_peak = self.__close_peak2
        #else: close_peak = self.__close_peak
        close_peak = self.__close_peak
        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            peak_content = None
            peak_length = 0
            (ps,vs) = self.get_data_by_chr(chrom) # arrays for position and values
            psn = iter(ps).__next__         # assign the next function to a viable to speed up
            vsn = iter(vs).__next__
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
                    close_peak(peak_content, peaks, min_length, chrom) #, smoothlen=max_gap / 2 )
                    # start a new peak
                    peak_content = [(pre_p,p,v),]
                pre_p = p

            # save the last peak
            if not peak_content:
                continue
            close_peak(peak_content, peaks, min_length, chrom) #, smoothlen=max_gap / 2 )
        return peaks

    def __close_peak( self, peak_content, peaks, int min_length, bytes chrom ):

        peak_length = peak_content[-1][1]-peak_content[0][0]
        if peak_length >= min_length: # if the peak is too small, reject it
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
            peaks.add( chrom,
                       peak_content[0][0],
                       peak_content[-1][1],
                       summit      = summit,
                       peak_score  = summit_value,
                       pileup      = 0,
                       pscore      = 0,
                       fold_change = 0,
                       qscore      = 0
                       )
            return True

    def call_broadpeaks (self, double lvl1_cutoff=500, double lvl2_cutoff=100, int min_length=200,
                         int lvl1_max_gap=50, int lvl2_max_gap=400):
        """This function try to find enriched regions within which,
        scores are continuously higher than a given cutoff for level
        1, and link them using the gap above level 2 cutoff with a
        maximum length of lvl2_max_gap.

        lvl1_cutoff:  cutoff of value at enriched regions, default 500.
        lvl2_cutoff:  cutoff of value at linkage regions, default 100.
        min_length :  minimum peak length, default 200.
        lvl1_max_gap   :  maximum gap to merge nearby enriched peaks, default 50.
        lvl2_max_gap   :  maximum length of linkage regions, default 400.
        colname: can be 'sample','control','-100logp','-100logq'. Cutoff will be applied to the specified column.

        Return both general PeakIO object for highly enriched regions
        and gapped broad regions in BroadPeakIO.
        """
        cdef bytes chrom
        cdef int i, j
        cdef set chrs
        #cdef int tmp_n

        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."
        lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff, min_length=min_length, max_gap=lvl1_max_gap)
        lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff, min_length=min_length, max_gap=lvl2_max_gap)
        chrs = lvl1_peaks.get_chr_names()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in chrs:
            #tmp_n = 0
            lvl1peakschrom = lvl1_peaks.get_data_from_chrom(chrom)
            lvl2peakschrom = lvl2_peaks.get_data_from_chrom(chrom)
            lvl1peakschrom_next = iter(lvl1peakschrom).__next__
            tmppeakset = []             # to temporarily store lvl1 region inside a lvl2 region
            # our assumption is lvl1 regions should be included in lvl2 regions
            try:
                lvl1 = lvl1peakschrom_next()
                for i in range( len(lvl2peakschrom) ):
                    # for each lvl2 peak, find all lvl1 peaks inside
                    lvl2 = lvl2peakschrom[i]
                    while True:
                        if lvl2["start"] <= lvl1["start"]  and lvl1["end"] <= lvl2["end"]:
                            tmppeakset.append(lvl1)
                            lvl1 = lvl1peakschrom_next()
                        else:
                            self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
                            #tmp_n += 1
                            tmppeakset = []
                            break
            except StopIteration:
                self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
                #tmp_n += 1
                tmppeakset = []
                for j in range( i+1, len(lvl2peakschrom) ):
                    self.__add_broadpeak ( broadpeaks, chrom, lvl2peakschrom[j], tmppeakset)
                    #tmp_n += 1

            #print len(lvl1peakschrom), len(lvl2peakschrom), tmp_n

        return broadpeaks

    def __add_broadpeak (self, bpeaks, chrom, lvl2peak, lvl1peakset):
        """Internal function to create broad peak.

        """
        cdef:
            long start, end, blockNum
            bytes blockSizes, blockStarts, thickStart, thickEnd

        start      = lvl2peak["start"]
        end        = lvl2peak["end"]

        # the following code will add those broad/lvl2 peaks with no strong/lvl1 peaks inside
        if not lvl1peakset:
            # try:
            # will complement by adding 1bps start and end to this region
            # may change in the future if gappedPeak format was improved.
            bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=(b"%d" % start), thickEnd=(b"%d" % end),
                       blockNum = 2, blockSizes = b"1,1", blockStarts = (b"0,%d" % (end-start-1)), pileup = lvl2peak["pileup"],
                       pscore = lvl2peak["pscore"], fold_change = lvl2peak["fc"],
                       qscore = lvl2peak["qscore"] )
            return bpeaks

        thickStart = b"%d" % lvl1peakset[0]["start"]
        thickEnd   = b"%d" % lvl1peakset[-1]["end"]
        blockNum   = len(lvl1peakset)
        blockSizes = b",".join( [b"%d" % x["length"] for x in lvl1peakset] )
        blockStarts = b",".join( [b"%d" % (x["start"]-start) for x in lvl1peakset] )

        if int(thickStart) != start:
            # add 1bp left block
            thickStart = b"%d" % start
            blockNum += 1
            blockSizes = b"1,"+blockSizes
            blockStarts = b"0,"+blockStarts
        if int(thickEnd) != end:
            # add 1bp right block
            thickEnd = b"%d" % end
            blockNum += 1
            blockSizes = blockSizes+b",1"
            blockStarts = blockStarts + b"," + (b"%d" % (end-start-1))

        bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=thickStart, thickEnd=thickEnd,
                   blockNum = blockNum, blockSizes = blockSizes, blockStarts = blockStarts,  pileup = lvl2peak["pileup"],
                   pscore = lvl2peak["pscore"], fold_change = lvl2peak["fc"],
                   qscore = lvl2peak["qscore"] )
        return bpeaks


    def total (self):
        """Return the number of regions in this object.

        """
        cdef long t
        t = 0
        for (p,s) in self.__data.values():
            t += len(p)
        return t

    def set_single_value (self, double new_value):
        """Change all the values in bedGraph to the same new_value,
        return a new bedGraphTrackI.

        """
        cdef:
            bytes chrom
            int max_p

        ret = bedGraphTrackI()
        chroms = set(self.get_chr_names())
        for chrom in chroms:
            (p1,v1) = self.get_data_by_chr(chrom) # arrays for position and values
            # maximum p
            max_p = max(p1)
            # add a region from 0 to max_p
            ret.add_loc(chrom,0,max_p,new_value)
        return ret

    def overlie (self, bdgTracks, func="max" ):
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

        Supported 'func' are "max", "mean" and "fisher".

        Return value is a bedGraphTrackI object.
        """
        cdef:
            int pre_p, p1, p2
            double v1, v2
            bytes chrom

        nr_tracks = len(bdgTracks) + 1  # +1 for self
        assert nr_tracks >= 2, "Specify at least two replicates"
        for i, bdgTrack in enumerate(bdgTracks):
            assert isinstance(bdgTrack, bedGraphTrackI), "bdgTrack{} is not a bedGraphTrackI object".format(i + 1)

        if func == "max":
            f = max
        elif func == "mean":
            #f = np_mean
            def f(*args):
                return sum(args)/nr_tracks
        elif func == "fisher":
            # combine -log10pvalues
            def f(*args):
                # chisq statistics = sum(-log10p)/log10(e)*2, chisq df = 2*number_of_reps
                return chisq_logp_e(2*sum(args)/LOG10_E, 2*nr_tracks, log10=True)
        else:
            raise Exception("Invalid function")

        ret = bedGraphTrackI()
        retadd = ret.add_loc

        common_chr = set(self.get_chr_names())
        for track in bdgTracks:
            common_chr = common_chr.intersection(set(track.get_chr_names()))

        for chrom in common_chr:
            datas = [self.get_data_by_chr(chrom)]
            datas.extend([bdgTracks[i].get_data_by_chr(chrom) for i in range(len(bdgTracks))])

            ps, vs, pn, vn = [], [], [], []
            for data in datas:
                ps.append(data[0])
                pn.append(iter(ps[-1]).__next__)
                vs.append(data[1])
                vn.append(iter(vs[-1]).__next__)

            pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            try:
                ps_cur = [pn[i]() for i in range(len(pn))]
                vs_cur = [vn[i]() for i in range(len(pn))]

                while True:
                    # get the lowest position
                    lowest_p = min(ps_cur)

                    # at least one lowest position, could be multiple
                    locations = [i for i in range(len(ps_cur)) if ps_cur[i] == lowest_p]

                    # add the data until the interval
                    ret.add_loc(chrom, pre_p, ps_cur[locations[0]], f(*vs_cur))

                    pre_p = ps_cur[locations[0]]
                    for index in locations:
                        ps_cur[index] = pn[index]()
                        vs_cur[index] = vn[index]()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass
        return ret

    def apply_func ( self, func ):
        """Apply function 'func' to every value in this bedGraphTrackI object.

        *Two adjacent regions with same value after applying func will
        not be merged.
        """
        cdef int i

        for (p,s) in self.__data.values():
            for i in xrange(len(s)):
                s[i] = func(s[i])
        self.maxvalue = func(self.maxvalue)
        self.minvalue = func(self.minvalue)
        return True

    def p2q ( self ):
        """Convert pvalue scores to qvalue scores.

        *Assume scores in this bedGraph are pvalue scores! Not work
         for other type of scores.
        """
        cdef:
            bytes chrom
            object pos_array, pscore_array
            dict pvalue_stat = {}
            dict pqtable = {}
            long n, pre_p, this_p, length, j, pre_l, l, i
            double this_v, pre_v, v, q, pre_q, this_t, this_c
            long N, k, this_l
            double f
            long nhcal = 0
            long npcal = 0
            list unique_values
            double t0, t1, t

        # calculate frequencies of each p-score
        for chrom in self.get_chr_names():
            pre_p = 0

            [pos_array, pscore_array] = self.__data[ chrom ]

            pn = iter(pos_array).__next__
            vn = iter(pscore_array).__next__

            for i in range( len( pos_array ) ):
                this_p = pn()
                this_v = vn()
                this_l = this_p - pre_p
                if this_v in pvalue_stat:
                    pvalue_stat[ this_v ] += this_l
                else:
                    pvalue_stat[ this_v ] = this_l
                pre_p = this_p

            nhcal += len( pos_array )

        nhval = 0

        N = sum(pvalue_stat.values()) # total length
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value

        # calculate qscore for each pscore
        pqtable = {}
        unique_values = sorted(pvalue_stat.keys(), reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = pvalue_stat[v]
            q = v + (log10(k) + f)
            q = max(0,min(pre_q,q))           # make q-score monotonic
            pqtable[ v ] = q
            pre_v = v
            pre_q = q
            k+=l
            nhcal += 1

        # convert pscore to qscore
        for chrom in self.get_chr_names():
            [pos_array, pscore_array] = self.__data[ chrom ]

            for i in range( len( pos_array ) ):
                pscore_array[ i ] = pqtable[ pscore_array[ i ] ]

        self.merge_regions()
        return


    def extract_value ( self, bdgTrack2 ):
        """It's like overlie function. THe overlapped regions between
        bdgTrack2 and self, will be recorded. The values from self in
        the overlapped regions will be outputed in a single array for
        follow statistics.

        """
        cdef:
            int pre_p, p1, p2, i
            double v1, v2
            bytes chrom

        assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = [[],array(FBYTE4,[]),array(BYTE4,[])] # 1: region in
                                                    # bdgTrack2; 2:
                                                    # value; 3: length
                                                    # with the value
        radd = ret[0].append
        vadd = ret[1].append
        ladd = ret[2].append

        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for i in range( len( common_chr ) ):
            chrom = common_chr.pop()
            (p1s,v1s) = self.get_data_by_chr(chrom) # arrays for position and values
            p1n = iter(p1s).__next__         # assign the next function to a viable to speed up
            v1n = iter(v1s).__next__

            (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
            p2n = iter(p2s).__next__         # assign the next function to a viable to speed up
            v2n = iter(v2s).__next__

            pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret

            try:
                p1 = p1n()
                v1 = v1n()

                p2 = p2n()
                v2 = v2n()

                while True:
                    if p1 < p2:
                        # clip a region from pre_p to p1, then set pre_p as p1.
                        if v2>0:
                            radd(chrom+"."+str(pre_p)+"."+str(p1))
                            vadd(v1)
                            ladd(p1-pre_p)
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        if v2>0:
                            radd(chrom+"."+str(pre_p)+"."+str(p2))
                            vadd(v1)
                            ladd(p2-pre_p)
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        if v2>0:
                            radd(chrom+"."+str(pre_p)+"."+str(p1))
                            vadd(v1)
                            ladd(p1-pre_p)
                        pre_p = p1
                        # call for the next p1, v1, p2, v2.
                        p1 = p1n()
                        v1 = v1n()
                        p2 = p2n()
                        v2 = v2n()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass

        return ret

    def make_scoreTrackII_for_macs (self, bdgTrack2, float depth1 = 1.0, float depth2 = 1.0 ):
        """A modified overlie function for MACS v2.

        effective_depth_in_million: sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        Return value is a bedGraphTrackI object.
        """
        cdef:
            int pre_p, p1, p2
            double v1, v2
            bytes chrom

        assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = scoreTrackII( treat_depth = depth1, ctrl_depth = depth2 )
        retadd = ret.add

        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for chrom in common_chr:

            (p1s,v1s) = self.get_data_by_chr(chrom) # arrays for position and values
            p1n = iter(p1s).__next__         # assign the next function to a viable to speed up
            v1n = iter(v1s).__next__

            (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
            p2n = iter(p2s).__next__         # assign the next function to a viable to speed up
            v2n = iter(v2s).__next__

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

        ret.finalize()
        #ret.merge_regions()
        return ret

    cpdef str cutoff_analysis ( self, int max_gap, int min_length, int steps = 100 ):
        cdef:
            set chrs
            list tmplist, peak_content
            bytes  chrom
            str ret
            float cutoff
            long total_l, total_p, i, n, ts, te, lastp, tl, peak_length
            dict cutoff_npeaks, cutoff_lpeaks
            float s, midvalue

        chrs = self.get_chr_names()

        midvalue = self.minvalue/2 + self.maxvalue/2
        s = float(self.minvalue - midvalue)/steps

        tmplist = list( np.arange( midvalue, self.minvalue - s, s ) )

        cutoff_npeaks = {}
        cutoff_lpeaks = {}

        for chrom in chrs:
            ( pos_array, score_array ) = self.__data[ chrom ]
            pos_array = np.array( self.__data[ chrom ][ 0 ] )
            score_array = np.array( self.__data[ chrom ][ 1 ] )

            for n in range( len( tmplist ) ):
                cutoff = round( tmplist[ n ], 3 )
                total_l = 0           # total length of peaks
                total_p = 0           # total number of peaks

                # get the regions with scores above cutoffs
                above_cutoff = np.nonzero( score_array > cutoff )[0]# this is not an optimized method. It would be better to store score array in a 2-D ndarray?
                above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
                above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

                if above_cutoff_endpos.size == 0:
                    continue

                # first bit of region above cutoff
                acs_next = iter(above_cutoff_startpos).__next__
                ace_next = iter(above_cutoff_endpos).__next__

                ts = acs_next()
                te = ace_next()
                peak_content = [( ts, te ), ]
                lastp = te

                for i in range( 1, above_cutoff_startpos.size ):
                    ts = acs_next()
                    te = ace_next()
                    tl = ts - lastp
                    if tl <= max_gap:
                        peak_content.append( ( ts, te ) )
                    else:
                        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                        if peak_length >= min_length: # if the peak is too small, reject it
                            total_l +=  peak_length
                            total_p += 1
                        peak_content = [ ( ts, te ), ]
                    lastp = te

                if peak_content:
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        total_l +=  peak_length
                        total_p += 1
                cutoff_lpeaks[ cutoff ] = cutoff_lpeaks.get( cutoff, 0 ) + total_l
                cutoff_npeaks[ cutoff ] = cutoff_npeaks.get( cutoff, 0 ) + total_p

        # write pvalue and total length of predicted peaks
        ret = "pscore\tnpeaks\tlpeaks\tavelpeak\n"
        for cutoff in sorted(cutoff_lpeaks.keys(), reverse=True):
            if cutoff_npeaks[ cutoff ] > 0:
                ret += "%.2f\t%d\t%d\t%.2f\n" % ( cutoff, cutoff_npeaks[ cutoff ], cutoff_lpeaks[ cutoff ], cutoff_lpeaks[ cutoff ]/cutoff_npeaks[ cutoff ] )
        return ret

def scoreTracktoBedGraph (scoretrack, str colname):
    """Produce a bedGraphTrackI object with certain column as scores.

    colname: can be 'sample','control','-100logp','-100logq'

    """
    cdef:
        int pre, i
        bytes chrom
        set chrs

    bdgtrack = bedGraphTrackI( baseline_value = 0 )
    if colname not in ['sample','control','-100logp','-100logq']:
        raise Exception("%s not supported!" % colname)
    if colname in ['-100logp', '-100logq']:
        flag100 = True              # for pvalue or qvalue, divide them by 100 while writing to bedGraph file
    else:
        flag100 = False
    chrs = scoretrack.get_chr_names()
    for chrom in chrs:
        d = scoretrack.data[chrom]
        l = scoretrack.pointer[chrom]
        pre = 0
        pos   = d['pos']
        if flag100:
            value = d[colname]/100.0
        else:
            value = d[colname]
        for i in xrange( l ):
            bdgtrack.add_loc( chrom, pre, pos[i] ,value[i] )
            pre = pos[i]

    return bdgtrack

class bedRegionTrackI (bedGraphTrackI):
    """A similar class to bedGraphTrackI, but is designed to save
    traditional 3-fields BED format data.

    """
    def __init__ (self):
        self.__data = {}
        self.maxvalue = 1
        self.minvalue = 0
        self.baseline_value = 0

    def safe_add_loc (self, bytes chromosome, int startpos, int endpos):
        """Add a chr-start-end-value block into __data dictionary.

        """
        cdef:
            int pre_pos
            double pre_v

        # basic assumption, end pos should > start pos
        assert endpos > startpos, "endpos %d can't be smaller than start pos %d" % (endpos,startpos)

        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0

        if chromosome not in self.__data:
            self.__data[chromosome] = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
            c = self.__data[chromosome]
            if startpos:
                # start pos is not 0, then add two blocks, the first
                # with "baseline_value"; the second with "value"
                c[0].append(startpos)
                c[1].append(self.baseline_value)
            c[0].append(endpos)
            c[1].append(1)
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
                c[1].append(1)
            else:
                # if this region is next to the previous one.
                if pre_v == 1:
                    # if value is the same, simply extend it.
                    c[0][-1] = endpos
                else:
                    # otherwise, add a new region
                    c[0].append(endpos)
                    c[1].append(1)


