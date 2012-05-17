# cython: profile=True
# Time-stamp: <2012-04-29 18:04:37 Tao Liu>

"""Module for filter duplicate tags from paired-end data

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
from MACS2.Constants import *
from logging import debug
import numpy as np
cimport numpy as np
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
from copy import copy
from cpython cimport bool

import sys

# Let numpy enforce PE-ness using ndarray, gives bonus speedup when sorting
# PE data doesn't have strandedness
class PETrackI:
    """Paired End Locations Track class I along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    def __init__ (self, char * anno=""):
        """fw is the fixed-width for all locations.
        
        """
        self.__locations = {}    # location pairs
        self.__pointer = {}      # location pairs
        self.__sorted = False
        self.total = 0           # total tags
        self.annotation = anno   # need to be figured out
        self.rlengths = None


    def add_loc ( self, str chromosome, np.ndarray[int32_t, ndim=2] loc):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        """
        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = np.zeros(shape=(BUFFER_SIZE, 2),
                                                    dtype='int32')
            self.__pointer[chromosome] = 0
        try:
            self.__locations[chromosome][self.__pointer[chromosome],:] = loc 
            self.__pointer[chromosome] += 1
        except IndexError:
            self.__expand__ ( self.__locations[chromosome] )
            self.__locations[chromosome][self.__pointer[chromosome],:] = loc
            self.__pointer[chromosome] += 1

    def __expand__ ( self, np.ndarray arr ):
        arr.resize((arr.shape[0] + BUFFER_SIZE, 2), refcheck = False )
        return

    def finalize ( self ):
        """ Resize np arrays for 5' positions and sort them in place """
        
        cdef int32_t i
        cdef str c
        
        self.total += 0 #??

        chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            self.__locations[c].resize((self.__pointer[c], 2), refcheck=False)
            self.__locations[c].sort(0)
            self.total += self.__locations[c].shape[0]

        self.__sorted = True
        return

    def get_locations_by_chr ( self, str chromosome ):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if self.__locations.has_key(chromosome):
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names ( self ):
        """Return all the chromosome names stored in this track object.
        """
        l = self.__locations.keys()
        l.sort()
        return l

    def length ( self ):
        """Total sequenced length = sum(end-start+1) over chroms   
        """
        cdef uint64_t l
        for v in self.__locations.values():
            l += (v[:,1] - v[:,0] + 1).sum() 
        return l

    def sort ( self ):
        """Naive sorting for locations.
        
        """
        cdef uint32_t i
        cdef str c

        chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            self.__locations[c].sort(0)

        self.__sorted = True


    def pmf( self ):
        """return a 1-D numpy array of the probabilities of observing each
        fragment size, indices are the bin number (first bin is 0)
        """
        cdef:
           np.ndarray[np.int64_t, ndim=1] counts = np.zeros(BUFFER_SIZE, dtype='int64')
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
            sizes = locs[:, 1] - locs[:, 0] # +1 ?? irrelevant for use
            bins = np.bincount(sizes).astype('int64')
            bins_len = bins.shape[0]
            if bins_len > max_bins: max_bins = bins_len
            if counts.shape[0] < max_bins:
                    counts.resize(max_bins, refcheck=False)
            counts += bins
        
        counts.resize(max_bins, refcheck=False)
        pmf = counts.astype('float64') / counts.astype('float64').sum()
        return 

    def get_dups(self ):
        """return a track of only the duplicated reads
        """
        cdef:
            int32_t i_chrom, n, start, end
            np.ndarray loc = np.zeros([1,2], np.int32)
            np.ndarray current_loc = np.zeros([1,2], np.int32)
            uint64_t i_old, i_new
            str k
                
        if not self.__sorted: self.sort()
        
        selfcopy = copy(self)
        selfcopy.__locations = {}
        selfcopy.__pointer = {}
            
        selfcopy.total = 0
        chrnames = self.get_chr_names()
        
        for i_chrom in range(len(chrnames)): # for each chromosome
            k = chrnames [ i_chrom ]
            selfcopy.__locations[k] = self.__locations[k].copy()
            selfcopy.__pointer[k] = self.__locations[k].copy()
            i_new = 0
            locs = selfcopy.__locations[k]
            size = locs.shape[0]
            if size < 1:
                new_locs = locs
            else:
                new_locs = np.zeros((selfcopy.__pointer[k], 2), dtype='int32')
                new_locs[i_new, :] = locs[i_new, :]
            
                current_loc = locs[0,:]
                for i_old in range(1, size):
                    loc = locs[i_old, :]
                    if (loc == current_loc).all():
                        new_locs[i_new, :] = loc
                        i_new += 1
                    else:
                        current_loc = loc
                new_locs.resize( (i_new, 2) )
                new_size = new_locs.shape[0]
                selfcopy.total += new_size
           # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            locs.resize( 100000, refcheck=False )
            locs.resize( 0, refcheck=False )
            # hope there would be no mem leak...
    
            selfcopy.__locations[k] = new_locs
        return selfcopy

    def filter_dup ( self, int maxnum=-1, bool keep_original=False):
        """Filter the duplicated reads.
    
        Run it right after you add all data into this object.
        """
        cdef:
            int32_t i_chrom, n, start, end
            np.ndarray loc = np.zeros([1,2], np.int32)
            np.ndarray current_loc = np.zeros([1,2], np.int32)
            uint64_t i_old, i_new
            str k
                
        if maxnum < 0: return # condition to return if not filtering
        
        if not self.__sorted: self.sort()
        
        if keep_original:
            selfcopy = copy(self)
            selfcopy.__locations = {}
            selfcopy.__pointer = {}
        else:
            selfcopy = self
            
        selfcopy.total = 0
        chrnames = self.get_chr_names()
        
        for i_chrom in range(len(chrnames)): # for each chromosome
            k = chrnames [ i_chrom ]
            if keep_original:
                selfcopy.__locations[k] = self.__locations[k].copy()
                selfcopy.__pointer[k] = self.__locations[k].copy()
            i_new = 0
            locs = selfcopy.__locations[k]
            size = locs.shape[0]
            if size < 1:
                new_locs = locs
            else:
                new_locs = np.zeros((selfcopy.__pointer[k], 2), dtype='int32')
                new_locs[i_new, :] = locs[i_new, :]
                n = 1                # the number of tags in the current location
            
                current_loc = locs[0,:]
                for i_old in range(1, size):
                    loc = locs[i_old, :]
                    if (loc == current_loc).all():
                        n += 1
                        if n <= maxnum:
                            new_locs[i_new, :] = loc
                            i_new += 1
                        else:
                            start, end = loc
                            debug("Duplicate fragments found at %s:%d-%d" % (k, start, end) )
                    else:
                        current_loc = loc
                        new_locs[i_new, :] = loc
                        i_new += 1
                        n = 1
                new_locs.resize( (i_new, 2) )
                new_size = new_locs.shape[0]
                selfcopy.total += new_size
           # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            locs.resize( 100000, refcheck=False )
            locs.resize( 0, refcheck=False )
            # hope there would be no mem leak...
    
            selfcopy.__locations[k] = new_locs
        return selfcopy

    def sample_percent (self, float percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef uint32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
        cdef str key
        
        self.total = 0

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            key = chrnames[ i_chrom ]
        
            num = <uint32_t>round(self.__locations[key].shape[0] * percent, 2 )
            np.random.shuffle( self.__locations[key] )
            self.__locations[key].resize( (num, 2) )
            self.__locations[key].sort(0)
            self.__pointer[key] = self.__locations[key].shape[0]
            
            self.total += self.__pointer[key]
        return

    def sample_num (self, uint64_t samplesize):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef float percent

        percent = float(samplesize)/self.total
        self.sample_percent ( percent )
        return

    def print_to_bed (self, fhd=None):
        """Output FWTrackIII to BED format files. If fhd is given,
        write to a file, otherwise, output to standard output.
        
        """
        cdef int32_t i, i_chrom, s, e
        cdef str k
        
        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, file)
        assert self.fw > 0, "FWTrackIII object .fw should be set larger than 0!"

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            k = chrnames[ i_chrom ]

            locs = self.__locations[k]

            for i in range(locs.shape[0]):
                s, e = locs[i, :]
                fhd.write("%s\t%d\t%d\t.\t.\t.\n" % (k, s, e))

        return
    