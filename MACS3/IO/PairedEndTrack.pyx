# cython: language_level=3
# cython: profile=True
# Time-stamp: <2019-11-04 13:23:29 taoliu>

"""Module for filter duplicate tags from paired-end data

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

from MACS3.Constants import *
import io
import sys
from logging import debug, info

import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
#ctypedef np.float64_t float64_t
#ctypedef np.float32_t float32_t
#ctypedef np.int64_t int64_t
#ctypedef np.int32_t int32_t
#ctypedef np.uint64_t uint64_t
#ctypedef np.uint32_t uint32_t

from copy import copy

from cpython cimport bool
cimport cython
#from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

from MACS3.Constants import *
from MACS3.Pileup import quick_pileup, max_over_two_pv_array, se_all_in_one_pileup

cdef INT_MAX = <int32_t>((<uint32_t>(-1))>>1)

cdef packed struct peLoc:
    int32_t l
    int32_t r

# Let numpy enforce PE-ness using ndarray, gives bonus speedup when sorting
# PE data doesn't have strandedness
cdef class PETrackI:
    """Paired End Locations Track class I along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    cdef:
        public dict __locations
        public dict __pointer
        public bool __sorted
        public uint64_t total
        public dict __dup_locations
        public dict __dup_pointer
        public bool __dup_sorted
        public uint64_t dup_total
        public object annotation
        public dict rlengths
        public object dups
        public int64_t buffer_size
        public int64_t length
        public float32_t average_template_length
        bool   __destroyed

    def __init__ (self, char * anno="", int64_t buffer_size = 100000 ):
        """fw is the fixed-width for all locations.

        """
        self.__locations = {}    # location pairs
        self.__pointer = {}      # location pairs
        self.__dup_locations = {}    # location pairs
        self.__dup_pointer = {}      # location pairs
        self.__sorted = False
        self.__dup_sorted = False
        self.total = 0           # total fragments
        self.dup_total = 0           # total fragments
        self.annotation = anno   # need to be figured out
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.average_template_length = 0.0

    cpdef void add_loc ( self, bytes chromosome, int32_t start, int32_t end):
        """Add a location to the list according to the sequence name.

        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        """
        cdef:
            int32_t i

        if chromosome not in self.__locations:
            self.__locations[chromosome] = np.zeros(shape=self.buffer_size, dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            self.__locations[chromosome][ 0 ] = ( start, end )
            self.__pointer[chromosome] = 1
        else:
            i = self.__pointer[chromosome]
            if i % self.buffer_size == 0:
                self.__expand__ ( self.__locations[chromosome] )
            self.__locations[chromosome][ i ] = ( start, end )
            self.__pointer[chromosome] = i + 1
        self.length += end - start
        return

    cpdef void destroy ( self ):
        """Destroy this object and release mem.
        """
        cdef:
            set chrs
            bytes chromosome

        chrs = self.get_chr_names()
        for chromosome in chrs:
            if chromosome in self.__locations:
                self.__locations[chromosome].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome].resize( 0, refcheck=False )
                self.__locations[chromosome] = None
                self.__locations.pop(chromosome)
            if chromosome in self.__dup_locations:
                self.__dup_locations[chromosome].resize( self.buffer_size, refcheck=False )
                self.__dup_locations[chromosome].resize( 0, refcheck=False )
                self.__dup_locations[chromosome] = None
                self.__dup_locations.pop(chromosome)
        self.__destroyed = True

        return

    cpdef void __expand__ ( self, np.ndarray arr ):
        arr.resize((arr.shape[0] + self.buffer_size), refcheck = False )
        return

    cpdef bint set_rlengths ( self, dict rlengths ):
        """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this fwtrack object will be updated.

        If chromosome in this fwtrack is not covered by given
        rlengths, and it has no associated length, it will be set as
        maximum integer.
        """
        cdef:
            set valid_chroms, missed_chroms
            bytes chrom

        valid_chroms = set(self.__locations.keys()).intersection(rlengths.keys())
        for chrom in valid_chroms:
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.__locations.keys()).difference(rlengths.keys())
        for chrom in missed_chroms:
            self.rlengths[chrom] = INT_MAX
        return True

    cpdef dict get_rlengths ( self ):
        """Get reference chromosome lengths dictionary.

        If self.rlengths is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.__locations.keys()])
        return self.rlengths

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
            self.__locations[c].resize((self.__pointer[c]), refcheck=False)
            self.__locations[c].sort( order=['l', 'r'] )
            self.total += self.__locations[c].shape[0]

        self.__sorted = True
        self.average_template_length = <float32_t>( self.length ) / self.total
        return

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

    # cpdef length ( self ):
    #     """Total sequenced length = sum(end-start) over chroms

    #     TL: efficient?
    #     """
    #     cdef:
    #         int64_t l = 0
    #         np.ndarray[np.int32_t, ndim=2] v
    #     for v in self.__locations.values():
    #         l += (v[:,1] - v[:,0]).sum()
    #     return l

    cpdef void sort ( self ):
        """Naive sorting for locations.

        """
        cdef:
            uint32_t i
            bytes c
            set chrnames

        chrnames = self.get_chr_names()

        for c in chrnames:
            #print "before", self.__locations[c][0:100]
            self.__locations[c].sort( order=['l', 'r'] ) # sort by the leftmost location
            #print "before", self.__locations[c][0:100]
        self.__sorted = True
        return

#    def centered_fake_fragments(track, int32_t d):
#        """Return a copy of the PETrackI object so that its locations are all
#        the same fixed distance d about the center of each fragment
#        """
#        new = PETrackI()
#        new.__pointer = deepcopy(self.__pointer)
#        new.rlengths = deepcopy(self.rlengths)
#        new.total = self.total
#        new.__sorted = self.__sorted
#        new.__locations = {}
#
#        five_shift = d / 2
#        three_shift = d - five_shift
#        for k, v in self.__locations.items():
#            midpoints = (v[:,0] + v[:,1]) / 2
#            new.__locations[k] = np.vstack((midpoints - five_shift,
#                                             midpoints + three_shift))
#        return new

    def pmf( self ):
        """return a 1-D numpy array of the probabilities of observing each
        fragment size, indices are the bin number (first bin is 0)
        """
        cdef:
           np.ndarray[np.int64_t, ndim=1] counts = np.zeros(self.buffer_size, dtype='int64')
           np.ndarray[np.float64_t, ndim=1] pmf
           np.ndarray[np.int64_t, ndim=1] bins
           np.ndarray[np.int32_t, ndim=1] sizes, locs
           bytes c
           int32_t i, bins_len
           int32_t max_bins = 0
           set chrnames

        chrnames = self.get_chr_names()

        for c in chrnames:
            locs = self.__locations[c]
            sizes = locs['r'] - locs['l'] # +1 ?? irrelevant for use
            bins = np.bincount(sizes).astype('int64')
            bins_len = bins.shape[0]
            if bins_len > max_bins: max_bins = bins_len
            if counts.shape[0] < max_bins:
                    counts.resize(max_bins, refcheck=False)
            counts += bins

        counts.resize(max_bins, refcheck=False)
        pmf = counts.astype('float64') / counts.astype('float64').sum()
        return pmf

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef void separate_dups ( self , int32_t maxnum = 1 ):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        cdef:
            int32_t i_chrom, n, start, end, size
#            np.ndarray[np.int32_t, ndim=1] loc #= np.zeros([1,2], np.int32)
#            np.ndarray[np.int32_t, ndim=1] current_loc #= np.zeros([1,2], np.int32)
            int32_t loc_start, loc_end, current_loc_start, current_loc_end
            np.ndarray locs, new_locs, dup_locs
            uint64_t i_old, i_new, i_dup, new_size, dup_size
            set chrnames
            bytes k

        if maxnum < 0: return # condition to return if not filtering

        chrnames = self.get_chr_names()

        if not self.__sorted: self.sort()

        self.__dup_pointer = copy(self.__pointer)
        self.dup_total = 0
        self.total = 0
        self.length = 0
        self.average_template_length = 0.0

        for k in chrnames: # for each chromosome
            locs = self.__locations[k]
            size = locs.shape[0]
            if size == 1:
                # do nothing and continue
                self.total +=  size
                self.__pointer[k] = size
                self.length += locs[0][1] - locs[0][0]
                self.__dup_locations[k] = []
                continue
            # save duplicate raeds to dup_locations[k] and update locations[k]
            new_locs = np.zeros(self.__pointer[k] + 1, dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            dup_locs = np.zeros(self.__pointer[k] + 1, dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            current_loc_start = locs[0][0] # same as locs[0]['l']
            current_loc_end = locs[0][1]# same as locs[0]['r']
            new_locs[0][0] = current_loc_start
            new_locs[0][1] = current_loc_end
            self.length += current_loc_end - current_loc_start
            i_new = 1           # index of new_locs
            i_dup = 0           # index of dup_locs
            n = 1
            for i_old in range(1, size):
                loc_start = locs[i_old][0]
                loc_end = locs[i_old][1]
                if (loc_start == current_loc_start) and (loc_end == current_loc_end) :
                    # both ends are the same
                    n += 1
                else:
                    # not the same, update current_loc and reset n
                    current_loc_start = loc_start
                    current_loc_end = loc_end
                    n = 1
                if n > maxnum:
                    # more than maxnum duplicates, put the duplicates in dup_locs
                    dup_locs[i_dup][0] = loc_start
                    dup_locs[i_dup][1] = loc_end
                    i_dup += 1
                else:
                    # otherwise, put the fragment in new_locs
                    new_locs[i_new][0] = loc_start
                    new_locs[i_new][1] = loc_end
                    self.length += loc_end - loc_start
                    i_new += 1
            new_locs.resize( i_new , refcheck = False)
            dup_locs.resize( i_dup , refcheck = False)
            self.total += i_new
            self.dup_total += i_dup
            self.__pointer[k] = i_new
            self.__dup_pointer[k] = i_dup
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            locs.resize( self.buffer_size, refcheck=False )
            locs.resize( 0, refcheck=False )
            # hope there would be no mem leak...
            self.__locations[k] = new_locs
            self.__dup_locations[k] = dup_locs
        self.average_template_length = self.length / self.total
        return

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef void filter_dup ( self, int32_t maxnum=-1):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        cdef:
            int32_t i_chrom, n, start, end
            int32_t loc_start, loc_end, current_loc_start, current_loc_end
            uint64_t i_old, i_new, size, new_size
            bytes k
            np.ndarray locs, new_locs
            set chrnames
            #peLoc current_range, this_range

        if maxnum < 0: return # condition to return if not filtering

        if not self.__sorted: self.sort()

        self.total = 0
        self.length = 0
        self.average_template_length = 0.0

        chrnames = self.get_chr_names()

        for k in chrnames: # for each chromosome
            locs = self.__locations[k]
            size = locs.shape[0]
            if size == 1:
                # do nothing and continue
                self.total += size
                self.__pointer[k] = size
                self.length += locs[0][1] - locs[0][0]
                continue
            # discard duplicate reads and make a new __locations[k]
            new_locs = np.zeros( self.__pointer[k] + 1, dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            # get the first loc
            #current_loc_start = locs[0][0]
            #current_loc_end = locs[0][1]
            ( current_loc_start, current_loc_end ) = locs[0]
            #new_locs[0][0] = current_loc_start
            #new_locs[0][1] = current_loc_end
            new_locs[0] = ( current_loc_start, current_loc_end )
            self.length += current_loc_end - current_loc_start
            i_new = 1           # index of new_locs
            n = 1               # the number of tags in the current genomic location
            for i_old in range(1, size):
                #loc_start = locs[i_old][0]
                #loc_end = locs[i_old][1]
                ( loc_start, loc_end ) = locs[i_old]
                if loc_start == current_loc_start and loc_end == current_loc_end:
                    # both ends are the same, add 1 to duplicate number n
                    n += 1
                else:
                    # not the same, update currnet_loc, reset n
                    current_loc_start = loc_start
                    current_loc_end = loc_end
                    n = 1
                if n <= maxnum:
                    # smaller than maxnum, then add to new_locs,
                    # otherwise, discard
                    #new_locs[i_new][0] = loc_start
                    #new_locs[i_new][1] = loc_end
                    new_locs[i_new] = (loc_start, loc_end)
                    self.length += loc_end - loc_start
                    i_new += 1
            new_locs.resize( i_new, refcheck = False )
            new_size = new_locs.shape[0]
            self.__pointer[k] = new_size
            self.total += new_size
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            locs.resize( self.buffer_size, refcheck=False )
            locs.resize( 0, refcheck=False )
            # hope there would be no mem leak...
            self.__locations[k] = new_locs
        self.average_template_length = self.length / self.total
        return

    cpdef void sample_percent (self, float32_t percent, int32_t seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            uint32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            bytes k
            set chrnames

        self.total = 0
        self.length = 0
        self.average_template_length = 0.0

        chrnames = self.get_chr_names()

        if seed >= 0:
            np.random.seed(seed)

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            num = <uint32_t>round(self.__locations[k].shape[0] * percent, 5 )
            np.random.shuffle( self.__locations[k] )
            self.__locations[k].resize( num, refcheck = False )
            self.__locations[k].sort( order = ['l', 'r'] ) # sort by leftmost positions
            self.__pointer[k] = self.__locations[k].shape[0]
            self.length += ( self.__locations[k]['r'] - self.__locations[k]['l'] ).sum()
            self.total += self.__pointer[k]
        self.average_template_length = <float32_t>( self.length )/ self.total
        return

    cpdef void sample_num (self, uint64_t samplesize, int32_t seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef float32_t percent

        percent = <float32_t>(samplesize)/self.total
        self.sample_percent ( percent, seed )
        return

    cpdef void print_to_bed (self, fhd=None):
        """Output to BEDPE format files. If fhd is given, write to a
        file, otherwise, output to standard output.

        """
        cdef:
            int32_t i, i_chrom, s, e
            bytes k
            set chrnames


        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, io.IOBase)

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            locs = self.__locations[k]

            for i in range(locs.shape[0]):
                s, e = locs[ i ]
                fhd.write("%s\t%d\t%d\n" % (k.decode(), s, e))

        return

    cpdef list pileup_a_chromosome ( self, bytes chrom, list scale_factor_s, float32_t baseline_value = 0.0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.

        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float32_t scale_factor

        #if not self.__sorted:
        #    self.sort()

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            scale_factor = scale_factor_s[i]

            tmp_pileup = quick_pileup ( np.sort(self.__locations[chrom]['l']), np.sort(self.__locations[chrom]['r']), scale_factor, baseline_value ) # Can't directly pass partial nparray there since that will mess up with pointer calculation.

            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        return prev_pileup

    cpdef list pileup_a_chromosome_c ( self, bytes chrom, list ds, list scale_factor_s, float32_t baseline_value = 0.0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.

        This function is for control track. Basically, here is a
        simplified function from FixWidthTrack. We pretend the PE is
        SE data and left read is on plus strand and right read is on
        minus strand.

        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float32_t scale_factor
            int64_t d, five_shift, three_shift
            int64_t rlength = self.get_rlengths()[chrom]

        if not self.__sorted: self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            five_shift = d//2
            three_shift= d//2

            tmp_pileup = se_all_in_one_pileup ( self.__locations[chrom]['l'], self.__locations[chrom]['r'], five_shift, three_shift, rlength, scale_factor, baseline_value )

            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        return prev_pileup


