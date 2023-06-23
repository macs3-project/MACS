# cython: language_level=3
# cython: profile=True
# Time-stamp: <2023-06-08 10:52:59 Tao Liu>

"""Module for BedGraph data class.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
#from array import array
from cpython cimport array
from array import array as pyarray
from math import prod
# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.ScoreTrack import ScoreTrackII
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO, RegionIO
from MACS3.Signal.Prob import chisq_logp_e

# ------------------------------------
# Other modules
# ------------------------------------

from cpython cimport bool
import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t

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

cdef inline mean_func( x ):
    return sum( x )/len( x )

cdef inline fisher_func( x ):
    # combine -log10pvalues
    return chisq_logp_e( 2*sum (x )/LOG10_E, 2*len( x ), log10=True )

cdef inline subtract_func( x ):
    # subtraction of two items list
    return x[1] - x[0]

cdef inline divide_func( x ):
    # division of two items list
    return x[1] / x[2]

cdef inline product_func( x ):
    # production of a list of values
    # only python 3.8 or above
    return prod( x )
    
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
        float32_t maxvalue
        float32_t minvalue
        float32_t baseline_value

    def __init__ (self, float32_t baseline_value=0 ):
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

    cpdef add_loc ( self, bytes chromosome, int32_t startpos, int32_t endpos, float32_t value):
        """Add a chr-start-end-value block into __data dictionary.

        Difference between safe_add_loc: no check, but faster. Save
        time while being called purely within MACS, so that regions
        are continuous without gaps. Note: I forgot that I removed
        safe_add_loc :P

        """
        cdef float32_t pre_v
        # basic assumption, end pos should > start pos

        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0

        if chromosome not in self.__data:
            self.__data[chromosome] = [ pyarray('i',[]), pyarray('f',[]) ]
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

    cpdef add_loc_wo_merge ( self, bytes chromosome, int32_t startpos, int32_t endpos, float32_t value):
        """Add a chr-start-end-value block into __data dictionary.

        Difference between safe_add_loc: no check, but faster. Save
        time while being called purely within MACS, so that regions
        are continuous without gaps. Note: I forgot that I removed
        safe_add_loc :P

        This one won't merge nearby ranges with the same value
        """
        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0

        if value < self.baseline_value:
            value = self.baseline_value
            
        if chromosome not in self.__data:
            self.__data[chromosome] = [ pyarray('i',[]), pyarray('f',[]) ]
            c = self.__data[chromosome]
            if startpos:
                # start pos is not 0, then add two blocks, the first
                # with "baseline_value"; the second with "value"
                c[0].append(startpos)
                c[1].append(self.baseline_value)
        c = self.__data[chromosome]
        c[0].append(endpos)
        c[1].append(value)
        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value        

    cpdef add_chrom_data( self, bytes chromosome, object p, object v ):
        """Add a pv data to a chromosome. Replace the previous data.

        p: a pyarray object 'i' for positions
        v: a pyarray object 'f' for values

        Note: no checks for error, use with caution
        """
        cdef:
            float32_t maxv, minv

        self.__data[ chromosome ] = [ p, v ]
        maxv = max( v )
        minv = min( v )
        if maxv > self.maxvalue:
            self.maxvalue = maxv
        if minv < self.minvalue:
            self.minvalue = minv
        return

    cpdef add_chrom_data_hmmr_PV( self, bytes chromosome, object pv ):
        """Add a pv data to a chromosome. Replace the previous data.

        This is a kinda silly function to waste time and convert a PV
        array (2-d named numpy array) into two python arrays for this
        BedGraph class. May have better function later.

        Note: no checks for error, use with caution
        """
        cdef:
            float32_t maxv, minv
            int32_t i

        self.__data[ chromosome ] = [ pyarray('i', pv['p']), pyarray('f',pv['v']) ]
        minv = pv['v'].min()
        maxv = pv['p'].max()
        if maxv > self.maxvalue:
            self.maxvalue = maxv
        if minv < self.minvalue:
            self.minvalue = minv
        return
    
    cpdef bool destroy ( self ):
        """ destroy content, free memory.
        """
        cdef:
            set chrs
            bytes chrom

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            if chrom in self.__data:
                self.__data[chrom] = [None, None]
                self.__data.pop(chrom)
        return True

    cpdef list get_data_by_chr (self, bytes chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if chromosome in self.__data:
            return self.__data[chromosome]
        else:
            return []

    cpdef set get_chr_names (self):
        """Return all the chromosome names stored.

        """
        return set(sorted(self.__data.keys()))

    cpdef void write_bedGraph (self, fhd, str name, str description, bool trackline=True):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        shift will be used to shift the coordinates. default: 0
        """
        cdef:
            int32_t pre, pos, i
            float32_t value
            bytes chrom
            set chrs

        if trackline:
            trackcontents = (name.replace("\"", "\\\""), description.replace("\"", "\\\""))
            fhd.write("track type=bedGraph name=\"%s\" description=\"%s\" visibility=2 alwaysZero=on\n" % trackcontents)
        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
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
        return

    cpdef void reset_baseline (self, float32_t baseline_value):
        """Reset baseline value to baseline_value.

        So any region between self.baseline_value and baseline_value
        will be set to baseline_value.

        """
        self.baseline_value = baseline_value
        self.filter_score(cutoff=baseline_value)
        self.merge_regions()
        return

    cdef merge_regions (self):
        """Merge nearby regions with the same value.

        """
        cdef:
            int32_t new_pre_pos, pos, i
            float32_t new_pre_value, value
            bytes chrom
            set chrs

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            (p,v) = self.__data[chrom]
            pnext = iter(p).__next__
            vnext = iter(v).__next__

            # new arrays
            new_pos = pyarray('L',[pnext(),])
            new_value = pyarray('f',[vnext(),])

            newpa = new_pos.append
            newva = new_value.append

            new_pre_pos = new_pos[0]
            new_pre_value = new_value[0]

            for i in range(1,len(p)):
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

    cpdef bool filter_score (self, float32_t cutoff=0):
        """Filter using a score cutoff. Any region lower than score
        cutoff will be set to self.baseline_value.

        Self will be modified.
        """
        cdef:
            int32_t new_pre_pos, pos, i
            float32_t new_pre_value, value
            bytes chrom
            set chrs

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            (p,v) = self.__data[chrom]
            pnext = iter(p).__next__
            vnext = iter(v).__next__

            # new arrays
            new_pos = pyarray('L',[])
            new_value = pyarray('f',[])
            new_pre_pos = 0
            new_pre_value = 0

            for i in range(len(p)):
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

    cpdef tuple summary (self):
        """Calculate the sum, total_length, max, min, mean, and std. 

        Return a tuple for (sum, total_length, max, min, mean, std).
        """
        cdef:
            int64_tn_v
            float32_t sum_v, max_v, min_v, mean_v, variance, tmp, std_v
            int32_t pre_p, l, i

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

    cpdef object call_peaks (self, float32_t cutoff=1,
                             int32_t min_length=200, int32_t max_gap=50,
                             bool call_summits=False):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 1.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.

        Removed option:

        up_limit: the highest acceptable value. Default 10^{310}
          * so only allow peak with value >=cutoff and <=up_limit

        This does not work. The region above upper limit may still be
        included as `gap` .

        """
        cdef:
            int32_t peak_length, x, pre_p, p, i, summit, tstart, tend
            float32_t v, summit_value, tvalue
            bytes chrom
            set chrs
            object peaks

        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in sorted(chrs):
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
                if v >= cutoff:
                    peak_content = [(pre_p,p,v),]
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in range(x,len(ps)):
                # continue scan the rest regions
                p = psn()
                v = vsn()
                if v < cutoff: # not be detected as 'peak'
                    pre_p = p
                    continue
                # for points above cutoff
                # if the gap is allowed
                if pre_p - peak_content[-1][1] <= max_gap:
                    peak_content.append((pre_p,p,v))
                else:
                    # when the gap is not allowed, close this peak
                    self.__close_peak(peak_content, peaks, min_length, chrom) #, smoothlen=max_gap / 2 )
                    # start a new peak
                    peak_content = [(pre_p,p,v),]
                pre_p = p

            # save the last peak
            if not peak_content:
                continue
            self.__close_peak(peak_content, peaks, min_length, chrom) #, smoothlen=max_gap / 2 )
        return peaks

    cdef bool __close_peak( self, list peak_content, object peaks, int32_t min_length, bytes chrom ):
        cdef:
            list tsummit        # list for temporary summits
            int32_t peak_length, summit, tstart, tend
            float32_t summit_value, tvalue
            
        peak_length = peak_content[-1][1]-peak_content[0][0]
        if peak_length >= min_length: # if the peak is too small, reject it
            tsummit = []
            summit = 0
            summit_value = 0
            for (tstart,tend,tvalue) in peak_content:
                if not summit_value or summit_value < tvalue:
                    tsummit = [<int32_t>((tend+tstart)/2),]
                    summit_value = tvalue
                elif summit_value == tvalue:
                    tsummit.append( <int32_t>((tend+tstart)/2) )
            summit = tsummit[<int32_t>((len(tsummit)+1)/2)-1 ]
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

    cpdef object call_broadpeaks (self, float32_t lvl1_cutoff=500, float32_t lvl2_cutoff=100,
                                  int32_t min_length=200, int32_t lvl1_max_gap=50, int32_t lvl2_max_gap=400):
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
        cdef:
            bytes chrom
            int32_t i, j
            set chrs
            object lvl1, lvl2   # PeakContent class object
            list temppeakset, lvl1peakschrom, lvl2peakschrom
            
        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."
        lvl1_peaks = self.call_peaks( cutoff=lvl1_cutoff, min_length=min_length, max_gap=lvl1_max_gap, call_summits=False )
        lvl2_peaks = self.call_peaks( cutoff=lvl2_cutoff, min_length=min_length, max_gap=lvl2_max_gap, call_summits=False )
        chrs = lvl1_peaks.get_chr_names()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in sorted(chrs):
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
                            tmppeakset = []
                            break
            except StopIteration:
                self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
                tmppeakset = []
                for j in range( i+1, len(lvl2peakschrom) ):
                    self.__add_broadpeak ( broadpeaks, chrom, lvl2peakschrom[j], tmppeakset)
        return broadpeaks

    cdef object __add_broadpeak (self, object bpeaks, bytes chrom, object lvl2peak, list lvl1peakset):
        """Internal function to create broad peak.

        """
        cdef:
            int32_t start, end, blockNum
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


    cpdef int32_t total (self):
        """Return the number of regions in this object.

        """
        cdef:
            int32_t t
        t = 0
        for ( p, v ) in self.__data.values():
            t += len(p)
        return t

    cpdef object set_single_value (self, float32_t new_value):
        """Change all the values in bedGraph to the same new_value,
        return a new bedGraphTrackI.

        """
        cdef:
            bytes chrom
            int32_t max_p
            object ret

        ret = bedGraphTrackI()
        chroms = set(self.get_chr_names())
        for chrom in sorted(chroms):
            (p1,v1) = self.get_data_by_chr(chrom) # arrays for position and values
            # maximum p
            max_p = max(p1)
            # add a region from 0 to max_p
            ret.add_loc(chrom,0,max_p,new_value)
        return ret

    cpdef object overlie (self, object bdgTracks, str func="max" ):
        """Calculate two or more bedGraphTrackI objects by letting self
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

        Supported 'func' are "sum", "subtract" (only for two bdg
        objects), "product", "divide" (only for two bdg objects),
        "max", "mean" and "fisher".

        Return value is a new bedGraphTrackI object.

        Option: bdgTracks can be a list of bedGraphTrackI objects
        """
        cdef:
            int32_t pre_p, p1, p2
            float32_t v1, v2
            bytes chrom

        nr_tracks = len(bdgTracks) + 1  # +1 for self
        assert nr_tracks >= 2, "Specify at least one more bdg objects."
        for i, bdgTrack in enumerate(bdgTracks):
            assert isinstance(bdgTrack, bedGraphTrackI), "bdgTrack{} is not a bedGraphTrackI object".format(i + 1)

        if func == "max":
            f = max
        elif func == "mean":
            f = mean_func
        elif func == "fisher":
            f = fisher_func
        elif func == "sum":
            f = sum
        elif func == "product":
            f = product_func
        elif func == "subtract":
            if nr_tracks == 2:
                f = subtract_func
            else:
                raise Exception(f"Only one more bdg object is allowed, but provided {nr_tracks-1}")
        elif func == "divide":
            if nr_tracks == 2:
                f = divide_func
            else:
                raise Exception(f"Only one more bdg object is allowed, but provided {nr_tracks-1}")
        else:
            raise Exception("Invalid function {func}! Choose from 'sum', 'subtract' (only for two bdg objects), 'product', 'divide' (only for two bdg objects), 'max', 'mean' and 'fisher'. ")

        ret = bedGraphTrackI()
        retadd = ret.add_loc

        common_chr = set(self.get_chr_names())
        for track in bdgTracks:
            common_chr = common_chr.intersection(set(track.get_chr_names()))

        for chrom in sorted(common_chr):
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
                    ret.add_loc(chrom, pre_p, ps_cur[locations[0]], f(vs_cur))

                    pre_p = ps_cur[locations[0]]
                    for index in locations:
                        ps_cur[index] = pn[index]()
                        vs_cur[index] = vn[index]()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass
        return ret

    cpdef bool apply_func ( self, func ):
        """Apply function 'func' to every value in this bedGraphTrackI object.

        *Two adjacent regions with same value after applying func will
        not be merged.
        """
        cdef int32_t i

        for (p,s) in self.__data.values():
            for i in range(len(s)):
                s[i] = func(s[i])
        self.maxvalue = func(self.maxvalue)
        self.minvalue = func(self.minvalue)
        return True

    cpdef p2q ( self ):
        """Convert pvalue scores to qvalue scores.

        *Assume scores in this bedGraph are pvalue scores! Not work
         for other type of scores.
        """
        cdef:
            bytes chrom
            object pos_array, pscore_array
            dict pvalue_stat = {}
            dict pqtable = {}
            int64_t n, pre_p, this_p, length, j, pre_l, l, i
            float32_t this_v, pre_v, v, q, pre_q, this_t, this_c
            int64_t N, k, this_l
            float32_t f
            int64_t nhcal = 0
            int64_t npcal = 0
            list unique_values
            float32_t t0, t1, t

        # calculate frequencies of each p-score
        for chrom in sorted(self.get_chr_names()):
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
        for chrom in sorted(self.get_chr_names()):
            [pos_array, pscore_array] = self.__data[ chrom ]

            for i in range( len( pos_array ) ):
                pscore_array[ i ] = pqtable[ pscore_array[ i ] ]

        self.merge_regions()
        return


    cpdef object extract_value ( self, object bdgTrack2 ):
        """Extract values from regions defined in bedGraphTrackI class object
        `regions`.

        """
        cdef:
            int32_t pre_p, p1, p2, i
            float32_t v1, v2
            bytes chrom
            object ret

        assert isinstance(bdgTrack2,bedGraphTrackI), "not a bedGraphTrackI object"

        ret = [ [], pyarray('f',[]), pyarray('L',[]) ] # 1: region in bdgTrack2; 2: value; 3: length with the value
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
                            radd(str(chrom)+"."+str(pre_p)+"."+str(p1))
                            vadd(v1)
                            ladd(p1-pre_p)
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        if v2>0:
                            radd(str(chrom)+"."+str(pre_p)+"."+str(p2))
                            vadd(v1)
                            ladd(p2-pre_p)
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        if v2>0:
                            radd(str(chrom)+"."+str(pre_p)+"."+str(p1))
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

    cpdef object extract_value_hmmr ( self, object bdgTrack2 ):
        """Extract values from regions defined in bedGraphTrackI class object
        `regions`.

        I will try to tweak this function to output only the values of
        bdgTrack1 (self) in the regions in bdgTrack2

        This is specifically for HMMRATAC. bdgTrack2 should be a
        bedgraph object containing the bins with value set to
        'mark_bin' -- the bins in the same region will have the same
        value.
        """
        cdef:
             int32_t pre_p, p1, p2, i
             float32_t v1, v2
             bytes chrom
             list ret

        assert isinstance(bdgTrack2,bedGraphTrackI), "not a bedGraphTrackI object"

        ret = [ [], pyarray('f',[]), pyarray('i',[]) ] # 0: bin location (chrom, position); 1: value; 2: number of bins in this region
        padd = ret[0].append
        vadd = ret[1].append
        ladd = ret[2].append

        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = sorted(list(chr1.intersection(chr2)))
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
                        # in this case, we don't output any
                        #if v2>0:
                        #    radd(str(chrom)+"."+str(pre_p)+"."+str(p1))
                        #    vadd(v1)
                        #    ladd(p1-pre_p)
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        if v2 != 0: #0 means it's a gap region, we should have value > 1
                            padd( (chrom, p2) )
                            vadd(v1)
                            ladd(int(v2))
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        if v2 != 0: #0 means it's a gap region, we should have 1 or -1
                            padd( (chrom, p2) )
                            vadd(v1)
                            ladd(int(v2))
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

    cpdef make_ScoreTrackII_for_macs (self, object bdgTrack2, float32_t depth1 = 1.0, float32_t depth2 = 1.0 ):
        """A modified overlie function for MACS v2.

        effective_depth_in_million: sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        Return value is a ScoreTrackII object.
        """
        cdef:
            int32_t pre_p, p1, p2
            float32_t v1, v2
            bytes chrom
            object ret

        assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = ScoreTrackII( treat_depth = depth1, ctrl_depth = depth2 )
        retadd = ret.add

        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for chrom in sorted(common_chr):

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

    cpdef str cutoff_analysis ( self, int32_t max_gap, int32_t min_length, int32_t steps = 100 ):
        """
        Cutoff analysis function for bedGraphTrackI object.
    
        This function will try all possible cutoff values on the score column to call peaks. Then 
        will give a report of a number of metrics (number of peaks, total length of peaks, average
        length of peak) at varying score cutoffs. For each score cutoff, the function finds the 
        positions where the score exceeds the cutoff, then groups those positions into "peaks" 
        based on the maximum allowed gap (max_gap) between consecutive positions. If a peak's length
        exceeds the minimum length (min_length), the peak is counted.

        Parameters
        ----------
        max_gap : int32_t
        Maximum allowed gap between consecutive positions above cutoff
        
        min_length : int32_t
        Minimum length of peak
        
        steps: int32_t

        Returns
        -------
        Cutoff analysis report in 'str'
        
        """
        cdef:
            set chrs
            list peak_content, ret_list, cutoff_list, cutoff_npeaks, cutoff_lpeaks
            bytes  chrom
            str ret
            float32_t cutoff
            int64_t total_l, total_p, i, n, ts, te, lastp, tl, peak_length
            #dict cutoff_npeaks, cutoff_lpeaks
            float32_t s, midvalue

        chrs = self.get_chr_names()

        #midvalue = self.minvalue/2 + self.maxvalue/2
        #s = float(self.minvalue - midvalue)/steps
        minv = self.minvalue
        maxv = self.maxvalue

        s = float(maxv - minv)/steps

        # a list of possible cutoff values from minv to maxv with step of s
        cutoff_list = [round(value, 3) for value in np.arange(minv, maxv, s)]

        cutoff_npeaks = [0] * len( cutoff_list )
        cutoff_lpeaks = [0] * len( cutoff_list )

        for chrom in sorted(chrs):
            ( pos_array, score_array ) = self.__data[ chrom ]
            pos_array = np.array( self.__data[ chrom ][ 0 ] )
            score_array = np.array( self.__data[ chrom ][ 1 ] )

            for n in range( len( cutoff_list ) ):
                cutoff = cutoff_list[ n ]
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
                cutoff_lpeaks[ n ] += total_l
                cutoff_npeaks[ n ] += total_p
                
        # prepare the returnning text
        ret_list = ["score\tnpeaks\tlpeaks\tavelpeak\n"]
        for n in range( len( cutoff_list )-1, -1, -1 ):
            cutoff = cutoff_list[ n ]
            if cutoff_npeaks[ n ] > 0:
                ret_list.append("%.2f\t%d\t%d\t%.2f\n" % ( cutoff, cutoff_npeaks[ n ], \
                                                          cutoff_lpeaks[ n ], \
                                                          cutoff_lpeaks[ n ]/cutoff_npeaks[ n ] ))
        ret = ''.join(ret_list)
        return ret


cdef np.ndarray calculate_elbows( np.ndarray values, float32_t threshold=0.01):
    # although this function is supposed to find elbow pts for cutoff analysis, 
    # however, in reality, it barely works...
    cdef: 
        np.ndarray deltas, slopes, delta_slopes, elbows
        np.float32_t avg_delta_slope
        
    # Calculate the difference between each point and the first point
    deltas = values - values[0]
    
    # Calculate the slope between each point and the last point
    slopes = deltas / (values[-1] - values[0])
    
    # Calculate the change in slope
    delta_slopes = np.diff(slopes)
    
    # Calculate the average change in slope
    avg_delta_slope = np.mean(delta_slopes)
    
    # Find all points where the change in slope is significantly larger than the average
    elbows = np.where(delta_slopes > avg_delta_slope + threshold)[0]
    
    return elbows