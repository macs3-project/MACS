# Time-stamp: <2015-04-20 14:22:09 Tao Liu>

"""Module for filter duplicate tags from paired-end data

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
from MACS2.Constants import *
from logging import debug, info
import numpy as np
cimport numpy as np
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
from copy import copy
from cpython cimport bool
cimport cython

from MACS2.Pileup import quick_pileup, max_over_two_pv_array, se_all_in_one_pileup
import sys

cdef INT_MAX = <int>((<unsigned int>-1)>>1)


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
        public unsigned long total
        public dict __dup_locations
        public dict __dup_pointer
        public bool __dup_sorted
        public unsigned long dup_total
        public object annotation
        public dict rlengths
        public object dups
        public long buffer_size
        public long length
        public float average_template_length
        bool   __destroyed
    
    def __init__ (self, char * anno="", long buffer_size = 100000 ):
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
        
    cpdef add_loc ( self, str chromosome, int start, int end):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        """
        cdef:
            long i

        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = np.zeros(shape=self.buffer_size, dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            self.__locations[chromosome][ 0 ] = ( start, end )
            self.__pointer[chromosome] = 1
        else:
            i = self.__pointer[chromosome]
            if i % self.buffer_size == 0:
                self.__expand__ ( self.__locations[chromosome] )
            self.__locations[chromosome][ i ] = ( start, end )
            self.__pointer[chromosome] += 1
        self.length += end - start

    cpdef destroy ( self ):
        """Destroy this object and release mem.
        """
        cdef:
            set chrs
            str chromosome
            
        chrs = set(self.get_chr_names())
        for chromosome in chrs:
            if self.__locations.has_key(chromosome):
                self.__locations[chromosome].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome].resize( 0, refcheck=False )
                self.__locations[chromosome] = None
                self.__locations.pop(chromosome)
            if self.__dup_locations.has_key(chromosome):
                self.__dup_locations[chromosome].resize( self.buffer_size, refcheck=False )
                self.__dup_locations[chromosome].resize( 0, refcheck=False )
                self.__dup_locations[chromosome] = None
                self.__dup_locations.pop(chromosome)
        self.__destroyed = True

        return True

    cpdef __expand__ ( self, np.ndarray arr ):
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
            str chrom

        valid_chroms = set(self.__locations.keys()).intersection(rlengths.keys())
        for chrom in valid_chroms:
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.__locations.keys()).difference(rlengths.keys())
        for chrom in missed_chroms:
            self.rlength[chrom] = INT_MAX
        return True

    cpdef dict get_rlengths ( self ):
        """Get reference chromosome lengths dictionary.

        If self.rlength is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.__locations.keys()])
        return self.rlengths
    

    def finalize ( self ):
        """ Resize np arrays for 5' positions and sort them in place """
        
        cdef int32_t i
        cdef str c
        
        self.total = 0

        chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            self.__locations[c].resize((self.__pointer[c]), refcheck=False)
            self.__locations[c].sort(order='l')
            self.total += self.__locations[c].shape[0]

        self.__sorted = True
        self.average_template_length = float( self.length ) / self.total
        return

    def get_locations_by_chr ( self, str chromosome ):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if self.__locations.has_key(chromosome):
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    cpdef list get_chr_names ( self ):
        """Return all the chromosome names stored in this track object.
        """
        l = self.__locations.keys()
        l.sort()
        return l

    # cpdef length ( self ):
    #     """Total sequenced length = sum(end-start) over chroms 

    #     TL: efficient?
    #     """
    #     cdef:
    #         long l = 0
    #         np.ndarray[np.int32_t, ndim=2] v
    #     for v in self.__locations.values():
    #         l += (v[:,1] - v[:,0]).sum() 
    #     return l

    cpdef sort ( self ):
        """Naive sorting for locations.
        
        """
        cdef uint32_t i
        cdef str c
        cdef list chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            #print "before", self.__locations[c][0:100]
            self.__locations[c].sort(order='l') # sort by the leftmost location
            #print "before", self.__locations[c][0:100]
        self.__sorted = True

#    def centered_fake_fragments(track, int d):
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
           str c
           int i, bins_len
           int max_bins = 0
           list chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
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
    def separate_dups ( self , int maxint = 1 ):
        """Filter the duplicated reads.
    
        Run it right after you add all data into this object.
        """
        cdef:
            int i_chrom, n, start, end, size
#            np.ndarray[np.int32_t, ndim=1] loc #= np.zeros([1,2], np.int32)
#            np.ndarray[np.int32_t, ndim=1] current_loc #= np.zeros([1,2], np.int32)
            int loc_start, loc_end, current_loc_start, current_loc_end
            np.ndarray locs, new_locs, dup_locs
            unsigned long i_old, i_new, i_dup, new_size, dup_size
            list chrnames = self.get_chr_names()
            str k
                        
        if not self.__sorted: self.sort()
        
        self.__dup_pointer = copy(self.__pointer)
        self.dup_total = 0
        self.total = 0
        self.length = 0
        self.average_template_length = 0.0

        for i_chrom in range(len(chrnames)): # for each chromosome
            k = chrnames [ i_chrom ]
#            dups.__locations[k] = self.__locations[k].copy()
            i_new = 0
            i_dup = 0
            locs = self.__locations[k]
            size = locs.shape[0]
            if size < 1:
                new_locs = locs
            else:
                new_locs = np.zeros(self.__pointer[k], dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
                dup_locs = np.zeros(self.__pointer[k], dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
                i_new += 1
                n = 1
            
                current_loc_start = locs[0][0] # same as locs[0]['l']
                current_loc_end = locs[0][1]# same as locs[0]['r']
                new_locs[i_new][0] = current_loc_start
                new_locs[i_new][1] = current_loc_end
                self.length += current_loc_end - current_loc_start
                for i_old in range(1, size):
                    loc_start = locs[i_old][0]
                    loc_end = locs[i_old][1]
                    all_same = ((loc_start == current_loc_start) and
                                (loc_end == current_loc_end)) 
                    if all_same:
                        n += 1
                    else:
                        current_loc_start = loc_start
                        current_loc_end = loc_end
                        n = 1
                    if n > maxint:
                        dup_locs[i_dup][0] = loc_start
                        dup_locs[i_dup][1] = loc_end
                        i_dup += 1
                    else:
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
                # unnecessary
#                new_size = new_locs.shape[0]
#                dup_size = dup_locs.shape[0]
#                self.__pointer[k] = new_size
#                dups.__pointer[k] = dup_size
           # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            locs.resize( self.buffer_size, refcheck=False )
            locs.resize( 0, refcheck=False )
            # hope there would be no mem leak...
    
            self.__locations[k] = new_locs
            self.__dup_locations[k] = dup_locs
        self.average_template_length = float( self.length ) / self.total
        return
    
    @cython.boundscheck(False) # do not check that np indices are valid
    def filter_dup ( self, int maxnum=-1):
        """Filter the duplicated reads.
    
        Run it right after you add all data into this object.
        """
        cdef:
            int i_chrom, n, start, end
#            np.ndarray[np.int32_t, ndim=1] loc #= np.zeros([1,2], np.int32)
#            np.ndarray[np.int32_t, ndim=1] current_loc #= np.zeros([1,2], np.int32)
            int loc_start, loc_end, current_loc_start, current_loc_end
            unsigned long i_old, i_new, size, new_size
            str k
            np.ndarray locs, new_locs
                
        if maxnum < 0: return # condition to return if not filtering
        
        if not self.__sorted: self.sort()
        
        self.total = 0
        self.length = 0
        self.average_template_length = 0.0
        
        chrnames = self.get_chr_names()
        
        for i_chrom in range(len(chrnames)): # for each chromosome
            k = chrnames [ i_chrom ]
            i_new = 0
            locs = self.__locations[k]
            size = locs.shape[0]
            if size < 1:
                new_locs = locs
            else:
                new_locs = np.zeros( self.__pointer[k], dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
                i_new += 1
                n = 1                # the number of tags in the current location
            
                current_loc_start = locs[0][0]
                current_loc_end = locs[0][1]
                new_locs[i_new][0] = current_loc_start
                new_locs[i_new][1] = current_loc_end
                self.length += current_loc_end - current_loc_start
                for i_old in range(1, size):
                    loc_start = locs[i_old][0]
                    loc_end = locs[i_old][1]
                    all_same = ((loc_start == current_loc_start) and
                                (loc_end == current_loc_end))
                    if all_same:
                        n += 1
                        if n <= maxnum:
                            new_locs[i_new][0] = loc_start
                            new_locs[i_new][1] = loc_end
                            self.length += loc_end - loc_start                            
                            i_new += 1
                    else:
                        current_loc_start = loc_start
                        current_loc_end = loc_end
                        new_locs[i_new][0] = loc_start
                        new_locs[i_new][1] = loc_end
                        self.length += loc_end - loc_start                        
                        i_new += 1
                        n = 1
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
        self.average_template_length = float( self.length ) / self.total
        return

    def sample_percent (self, float percent, int seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            uint32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            str key
        
        self.total = 0
        self.length = 0
        self.average_template_length = 0.0
        
        chrnames = self.get_chr_names()

        if seed >= 0:
            np.random.seed(seed)
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            key = chrnames[ i_chrom ]
        
            num = <uint32_t>round(self.__locations[key].shape[0] * percent, 5 )
            np.random.shuffle( self.__locations[key] )
            self.__locations[key].resize( num, refcheck = False )
            self.__locations[key].sort( order = 'l' ) # sort by leftmost positions
            self.__pointer[key] = self.__locations[key].shape[0]
            self.length += ( self.__locations[key]['r'] - self.__locations[key]['l'] ).sum()
            self.total += self.__pointer[key]
        self.average_template_length = float( self.length )/ self.total
        return

    def sample_num (self, uint64_t samplesize, int seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef float percent

        percent = float(samplesize)/self.total
        self.sample_percent ( percent, seed )
        return

    def print_to_bed (self, fhd=None):
        """Output to BED format files. If fhd is given, write to a
        file, otherwise, output to standard output.
        
        """
        cdef int32_t i, i_chrom, s, e
        cdef str k
        
        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, file)
        assert self.fw > 0, "width should be set larger than 0!"

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            k = chrnames[ i_chrom ]

            locs = self.__locations[k]

            for i in range(locs.shape[0]):
                s, e = locs[ i ]
                fhd.write("%s\t%d\t%d\t.\t.\t.\n" % (k, s, e))

        return
    
    cpdef pileup_a_chromosome ( self, str chrom, list scale_factor_s, float baseline_value = 0.0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.
        
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float scale_factor
            
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

    cpdef pileup_a_chromosome_c ( self, str chrom, list ds, list scale_factor_s, float baseline_value = 0.0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.

        This function is for control track. Basically, here is a simplified function from FixWidthTrack.
        
        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float scale_factor
            long d, five_shift, three_shift
            long rlength = self.get_rlengths()[chrom]

        if not self.__sorted: self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            five_shift = d/2
            three_shift= d/2

            tmp_pileup = se_all_in_one_pileup ( self.__locations[chrom]['l'], self.__locations[chrom]['r'], five_shift, three_shift, rlength, scale_factor, baseline_value )

            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup

        return prev_pileup


