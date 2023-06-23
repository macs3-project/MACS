# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-09-15 17:17:37 Tao Liu>

"""Module for FWTrack classes.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import io
from copy import copy
from collections import Counter

# ------------------------------------
# MACS3 modules
# ------------------------------------

from MACS3.Utilities.Constants import *
from MACS3.Signal.SignalProcessing import *
from MACS3.IO.PeakIO import PeakIO
from MACS3.Signal.Pileup import se_all_in_one_pileup, over_two_pv_array

# ------------------------------------
# Other modules
# ------------------------------------
from cpython cimport bool
cimport cython
import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t

# ------------------------------------
# constants
# ------------------------------------
__version__ = "FixWidthTrack $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "FWTrack class"

cdef INT_MAX = <int32_t>((<uint32_t>(-1))>>1)

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class FWTrack:
    """Fixed Width Locations Track class  along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    cdef:
        dict __locations
        dict __pointer
        dict __buf_size
        bool __sorted
        bool __destroyed
        dict rlengths
        public int64_t buffer_size
        public int64_t total
        public object annotation
        public object dups
        public int32_t fw
        public int64_t length

    def __init__ (self, int32_t fw=0, char * anno="", int64_t buffer_size = 100000 ):
        """fw is the fixed-width for all locations.

        """
        self.fw = fw
        self.__locations = {}    # location pairs: two strands
        self.__pointer = {}      # location pairs
        self.__buf_size = {}     # location pairs
        self.__sorted = False
        self.total = 0           # total tags
        self.annotation = anno   # need to be figured out
        self.rlengths = {}       # lengths of reference sequences, e.g. each chromosome in a genome
        self.buffer_size = buffer_size
        self.length = 0
        self.__destroyed = False

    cpdef void destroy ( self ):
        """Destroy this object and release mem.
        """
        cdef:
            set chrs
            bytes chromosome

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.__locations:
                self.__locations[chromosome][0].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome][0].resize( 0, refcheck=False )
                self.__locations[chromosome][1].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome][1].resize( 0, refcheck=False )
                self.__locations[chromosome] = [None, None]
                self.__locations.pop(chromosome)
        self.__destroyed = True
        return

    cpdef void add_loc ( self, bytes chromosome, int32_t fiveendpos, int32_t strand ):
        """Add a location to the list according to the sequence name.

        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for minus strand
        strand     -- 0: plus, 1: minus
        """
        cdef:
            int32_t i
            int32_t b
            np.ndarray arr

        if chromosome not in self.__locations:
            self.__buf_size[chromosome] = [ self.buffer_size, self.buffer_size ]
            self.__locations[chromosome] = [ np.zeros(self.buffer_size, dtype='int32'), np.zeros(self.buffer_size, dtype='int32') ] # [plus,minus strand]
            self.__pointer[chromosome] = [ 0, 0 ]
            self.__locations[chromosome][strand][0] = fiveendpos
            self.__pointer[chromosome][strand] = 1
        else:
            i = self.__pointer[chromosome][strand]
            b = self.__buf_size[chromosome][strand]
            arr = self.__locations[chromosome][strand]
            if b == i:
                b += self.buffer_size
                arr.resize( b, refcheck = False )
                self.__buf_size[chromosome][strand] = b
            arr[i]= fiveendpos
            self.__pointer[chromosome][strand] += 1
        return

    cpdef void finalize ( self ):
        """ Resize np arrays for 5' positions and sort them in place

        Note: If this function is called, it's impossible to append more files to this FWTrack object. So remember to call it after all the files are read!
        """

        cdef:
            int32_t i
            bytes c
            set chrnames

        self.total = 0

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.__locations[c][0].resize( self.__pointer[c][0], refcheck=False )
            self.__locations[c][0].sort()
            self.__locations[c][1].resize( self.__pointer[c][1], refcheck=False )
            self.__locations[c][1].sort()
            self.total += self.__locations[c][0].size + self.__locations[c][1].size

        self.__sorted = True
        self.length = self.fw * self.total
        return

    cpdef bint set_rlengths ( self, dict rlengths ):
        """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this fwtrack object will be updated.

        If chromosome in this fwtrack is not covered by given
        rlengths, and it has no associated length, it will be set as
        maximum integer.

        """
        cdef:
            set valid_chroms, missed_chroms, extra_chroms
            bytes chrom

        valid_chroms = set(self.__locations.keys()).intersection(rlengths.keys())
        for chrom in sorted(valid_chroms):
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.__locations.keys()).difference(rlengths.keys())
        for chrom in sorted(missed_chroms):
            self.rlengths[chrom] = INT_MAX
        return True

    cpdef dict get_rlengths ( self ):
        """Get reference chromosome lengths dictionary.

        If self.rlength is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.__locations.keys()])
        return self.rlengths

    cpdef get_locations_by_chr ( self, bytes chromosome ):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if chromosome in self.__locations:
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    cpdef set get_chr_names ( self ):
        """Return all the chromosome names stored in this track object.
        """
        return set(sorted(self.__locations.keys()))

    cpdef void sort ( self ):
        """Naive sorting for locations.

        """
        cdef:
            int32_t i
            bytes c
            set chrnames

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.__locations[c][0].sort()
            self.__locations[c][1].sort()

        self.__sorted = True
        return

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef uint64_t filter_dup ( self, int32_t maxnum = -1):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.

        Note, this function will *throw out* duplicates
        permenantly. If you want to keep them, use separate_dups
        instead.
        """
        cdef:
            int32_t p, m, n, current_loc
            # index for old array, and index for new one
            uint64_t i_old, i_new, size, new_size
            bytes k
            np.ndarray[int32_t, ndim=1] plus, new_plus, minus, new_minus
            set chrnames

        if maxnum < 0: return self.total         # do nothing

        if not self.__sorted:
            self.sort()

        self.total = 0
        self.length = 0

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            # + strand
            i_new = 0
            plus = self.__locations[k][0]
            size = plus.shape[0]
            if len(plus) <= 1:
                new_plus = plus         # do nothing
            else:
                new_plus = np.zeros( self.__pointer[k][0] + 1,dtype='int32' )
                new_plus[ i_new ] = plus[ i_new ] # first item
                i_new += 1
                n = 1                # the number of tags in the current location
                current_loc = plus[0]
                for i_old in range( 1, size ):
                    p = plus[ i_old ]
                    if p == current_loc:
                        n += 1
                    else:
                        current_loc = p
                        n = 1
                    if n <= maxnum:
                        new_plus[ i_new ] = p
                        i_new += 1
                new_plus.resize( i_new, refcheck=False )
                self.total +=  i_new
                self.__pointer[k][0] = i_new
                # free memory?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                plus.resize( self.buffer_size, refcheck=False )
                plus.resize( 0, refcheck=False )
                # hope there would be no mem leak...

            # - strand
            i_new = 0
            minus = self.__locations[k][1]
            size = minus.shape[0]
            if len(minus) <= 1:
                new_minus = minus         # do nothing
            else:
                new_minus = np.zeros( self.__pointer[k][1] + 1,dtype='int32' )
                new_minus[ i_new ] = minus[ i_new ] # first item
                i_new += 1
                n = 1                # the number of tags in the current location
                current_loc = minus[0]
                for i_old in range( 1, size ):
                    p = minus[ i_old ]
                    if p == current_loc:
                        n += 1
                    else:
                        current_loc = p
                        n = 1
                    if n <= maxnum:
                        new_minus[ i_new ] = p
                        i_new += 1
                new_minus.resize( i_new, refcheck=False )
                self.total +=  i_new
                self.__pointer[k][1] = i_new
                # free memory ?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                minus.resize( self.buffer_size, refcheck=False )
                minus.resize( 0, refcheck=False )
                # hope there would be no mem leak...

            self.__locations[k]=[new_plus,new_minus]

        self.length = self.fw * self.total
        return self.total

    cpdef void sample_percent (self, float32_t percent, int32_t seed = -1 ):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            int32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            bytes k
            set chrnames

        self.total = 0
        self.length = 0

        chrnames = self.get_chr_names()

        if seed >= 0:
            np.random.seed(seed)

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            num = <int32_t>round(self.__locations[k][0].shape[0] * percent, 5 )
            np.random.shuffle( self.__locations[k][0] )
            self.__locations[k][0].resize( num, refcheck=False )
            self.__locations[k][0].sort()
            self.__pointer[k][0] = self.__locations[k][0].shape[0]

            num = <int32_t>round(self.__locations[k][1].shape[0] * percent, 5 )
            np.random.shuffle( self.__locations[k][1] )
            self.__locations[k][1].resize( num, refcheck=False )
            self.__locations[k][1].sort()
            self.__pointer[k][1] = self.__locations[k][1].shape[0]

            self.total += self.__pointer[k][0] + self.__pointer[k][1]

        self.length = self.fw * self.total
        return

    cpdef void sample_num (self, uint64_t samplesize, int32_t seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            float32_t percent

        percent = <float32_t>(samplesize)/self.total
        self.sample_percent ( percent, seed )
        return

    cpdef void print_to_bed (self, fhd=None):
        """Output FWTrack to BED format files. If fhd is given,
        write to a file, otherwise, output to standard output.

        """
        cdef:
            int32_t i, i_chrom, p
            bytes k
            set chrnames

        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd,io.IOBase)
        assert self.fw > 0, "FWTrack object .fw should be set larger than 0!"

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            plus = self.__locations[k][0]

            for i in range(plus.shape[0]):
                p = plus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k.decode(),p,p+self.fw,"+") )

            minus = self.__locations[k][1]

            for i in range(minus.shape[0]):
                p = minus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k.decode(),p-self.fw,p,"-") )
        return

    cpdef tuple extract_region_tags ( self, bytes chromosome, int32_t startpos, int32_t endpos ):
        cdef:
            int32_t i, pos
            np.ndarray[int32_t, ndim=1] rt_plus, rt_minus
            list temp
            set chrnames

        if not self.__sorted: self.sort()

        chrnames = self.get_chr_names()
        assert chromosome in chrnames, "chromosome %s can't be found in the FWTrack object." % chromosome

        (plus, minus) = self.__locations[chromosome]

        temp = []
        for i in range(plus.shape[0]):
            pos = plus[i]
            if pos < startpos:
                continue
            elif pos > endpos:
                break
            else:
                temp.append(pos)
        rt_plus = np.array(temp)

        temp = []
        for i in range(minus.shape[0]):
            pos = minus[i]
            if pos < startpos:
                continue
            elif pos > endpos:
                break
            else:
                temp.append(pos)
        rt_minus = np.array(temp)
        return (rt_plus, rt_minus)

    cpdef list compute_region_tags_from_peaks ( self, peaks, func, int32_t window_size = 100, float32_t cutoff = 5 ):
        """Extract tags in peak, then apply func on extracted tags.

        peaks: redefined regions to extract raw tags in PeakIO type: check cPeakIO.pyx.

        func:  a function to compute *something* from tags found in a predefined region

        window_size: this will be passed to func.

        cutoff: this will be passed to func.

        func needs the fixed number of parameters, so it's not flexible. Here is an example:

        wtd_find_summit(chrom, plus, minus, peak_start, peak_end, name , window_size, cutoff):

        """

        cdef:
            int32_t m, i, j, pre_i, pre_j, pos, startpos, endpos
            np.ndarray[int32_t, ndim=1] plus, minus, rt_plus, rt_minus
            bytes chrom, name
            list temp, retval
            set pchrnames, chrnames

        pchrnames = peaks.get_chr_names()
        retval = []

        # this object should be sorted
        if not self.__sorted: self.sort()
        # PeakIO object should be sorted
        peaks.sort()

        chrnames = self.get_chr_names()

        for chrom in sorted(pchrnames):
            assert chrom in chrnames, "chromosome %s can't be found in the FWTrack object." % chrom
            (plus, minus) = self.__locations[chrom]
            cpeaks = peaks.get_data_from_chrom(chrom)
            prev_i = 0
            prev_j = 0
            for m in range(len(cpeaks)):
                startpos = cpeaks[m]["start"] - window_size
                endpos   = cpeaks[m]["end"] + window_size
                name     = cpeaks[m]["name"]

                temp = []
                for i in range(prev_i,plus.shape[0]):
                    pos = plus[i]
                    if pos < startpos:
                        continue
                    elif pos > endpos:
                        prev_i = i
                        break
                    else:
                        temp.append(pos)
                rt_plus = np.array(temp, dtype="int32")

                temp = []
                for j in range(prev_j,minus.shape[0]):
                    pos = minus[j]
                    if pos < startpos:
                        continue
                    elif pos > endpos:
                        prev_j = j
                        break
                    else:
                        temp.append(pos)
                rt_minus = np.array(temp, dtype="int32")

                retval.append( func(chrom, rt_plus, rt_minus, startpos, endpos, name = name, window_size = window_size, cutoff = cutoff) )
                # rewind window_size
                for i in range(prev_i, 0, -1):
                    if plus[prev_i] - plus[i] >= window_size:
                        break
                prev_i = i

                for j in range(prev_j, 0, -1):
                    if minus[prev_j] - minus[j] >= window_size:
                        break
                prev_j = j
                # end of a loop

        return retval

    cpdef list pileup_a_chromosome ( self, bytes chrom, list ds, list scale_factor_s, float32_t baseline_value = 0.0, bint directional = True, int32_t end_shift = 0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.

        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        directional    : if False, the strand or direction of tag will be ignored, so that extension will be both sides with d/2.
        end_shift      : move cutting ends towards 5->3 direction if value is positive, or towards 3->5 direction if negative. Default is 0 -- no shift at all.


        p and v are numpy.ndarray objects.
        """
        cdef:
            int64_t d
            int64_t five_shift, three_shift  # adjustment to 5' end and 3' end positions to make a fragment
            dict chrlengths = self.get_rlengths ()
            int64_t rlength = chrlengths[chrom]
            object ends
            list five_shift_s = []
            list three_shift_s = []
            list tmp_pileup, prev_pileup

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        # adjust extension length according to 'directional' and 'halfextension' setting.
        for d in ds:
            if directional:
                # only extend to 3' side
                five_shift_s.append(  - end_shift )
                three_shift_s.append( end_shift + d)
            else:
                # both sides
                five_shift_s.append( d//2 - end_shift )
                three_shift_s.append( end_shift + d - d//2)

        prev_pileup = None

        for i in range(len(ds)):
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]
            tmp_pileup = se_all_in_one_pileup ( self.__locations[chrom][0], self.__locations[chrom][1], five_shift, three_shift, rlength, scale_factor, baseline_value )

            if prev_pileup:
                prev_pileup = over_two_pv_array ( prev_pileup, tmp_pileup, func="max" )
            else:
                prev_pileup = tmp_pileup

        return prev_pileup

cdef inline int32_t left_sum ( data, int32_t pos, int32_t width ):
    """
    """
    return sum([data[x] for x in data if x <= pos and x >= pos - width])

cdef inline int32_t right_sum ( data, int32_t pos, int32_t width ):
    """
    """
    return sum([data[x] for x in data if x >= pos and x <= pos + width])

cdef inline int32_t left_forward ( data, int32_t pos, int32_t window_size ):
    return data.get(pos,0) - data.get(pos-window_size, 0)

cdef inline int32_t right_forward ( data, int32_t pos, int32_t window_size ):
    return data.get(pos + window_size, 0) - data.get(pos, 0)

