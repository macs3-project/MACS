# cython: profile=True
# Time-stamp: <2012-08-01 18:08:12 Tao Liu>

"""Module for FWTrack classes.

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

#from array import array
#from random import sample as random_sample
import sys
from copy import copy
from MACS2.Constants import *
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
from cpython cimport bool
cimport cython

import numpy as np
cimport numpy as np

# ------------------------------------
# constants
# ------------------------------------
__version__ = "FixWidthTrack $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "FWTrackII class"

cdef INT_MAX = <int>((<unsigned int>-1)>>1)

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class FWTrackIII:
    """Fixed Width Locations Track class III along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    cdef:
        dict __locations
        dict __pointer
        bool __sorted
        dict __dup_locations
        dict __dup_pointer
        bool __dup_sorted
        dict rlengths
        public long total
        public unsigned long dup_total
        public object annotation

        public object dups
        public int fw
    
    def __init__ (self, int32_t fw=0, char * anno=""):
        """fw is the fixed-width for all locations.
        
        """
        self.fw = fw
        self.__locations = {}    # location pairs
        self.__pointer = {}      # location pairs
        self.__dup_locations = {}    # location pairs
        self.__dup_pointer = {}      # location pairs
        self.__sorted = False
        self.__dup_sorted = False
        self.total = 0           # total tags
        self.dup_total = 0           # total tags
        self.annotation = anno   # need to be figured out
        self.rlengths = {}


    cpdef add_loc ( self, str chromosome, int32_t fiveendpos, int strand ):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = [ np.zeros(BUFFER_SIZE, dtype='int32'), np.zeros(BUFFER_SIZE, dtype='int32') ] # [plus,minus strand]
            self.__pointer[chromosome] = [ 0, 0 ]
        try:
            self.__locations[chromosome][strand][self.__pointer[chromosome][strand]] = fiveendpos
            self.__pointer[chromosome][strand] += 1
        except IndexError:
            self.__expand__ ( self.__locations[chromosome][strand] )
            self.__locations[chromosome][strand][self.__pointer[chromosome][strand]] = fiveendpos
            self.__pointer[chromosome][strand] += 1
        #print chromosome

    cpdef __expand__ ( self, np.ndarray arr ):
        arr.resize( arr.size + BUFFER_SIZE, refcheck = False )
        return

    cpdef finalize ( self ):
        """ Resize np arrays for 5' positions and sort them in place """
        
        cdef:
            int32_t i
            str c
        
        self.total = 0

        chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            self.__locations[c][0].resize( self.__pointer[c][0], refcheck=False )
            self.__locations[c][0].sort()
            self.__locations[c][1].resize( self.__pointer[c][1], refcheck=False )
            self.__locations[c][1].sort()
            self.total += self.__locations[c][0].size + self.__locations[c][1].size

        self.__sorted = True
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

    cpdef get_locations_by_chr ( self, str chromosome ):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if self.__locations.has_key(chromosome):
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    cpdef get_chr_names ( self ):
        """Return all the chromosome names stored in this track object.
        """
        l = self.__locations.keys()
        l.sort()
        return l

    cpdef length ( self ):
        """Total sequenced length = total number of tags * width of tag		
        """
        return self.total*self.fw

    cpdef sort ( self ):
        """Naive sorting for locations.
        
        """
        cdef:
            int32_t i
            str c

        chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            self.__locations[c][0].sort()
            self.__locations[c][1].sort()

        self.__sorted = True

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef separate_dups( self, maxint = 1 ):
        """Separate the duplicated reads into a different track
        stored at self.dup
        """
        cdef:
            int p, m, n, current_loc, i_chrom
            unsigned long i_old, i_new          # index for old array, and index for new one
            unsigned long i_dup, size, new_size, dup_size
            str k
            np.ndarray plus, new_plus, dup_plus, minus, new_minus, dup_minus

        if not self.__sorted:
            self.sort()

        self.__dup_pointer = copy(self.__pointer)
        self.dup_total = 0
        self.total = 0

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            k = chrnames[ i_chrom ]
#            dups.__locations[k] = self.__locations[k].copy()
            # + strand
            i_new = 0
            i_dup = 0
            plus = self.__locations[k][0]
            size = plus.shape[0]
            if len(plus) < 1:
                new_plus = plus         # do nothing
            else:
                new_plus = np.zeros( self.__pointer[k][0],dtype='int32' )
                dup_plus = np.zeros( self.__pointer[k][0],dtype='int32' )
                new_plus[ i_new ] = plus[ i_new ] # first item
                i_new += 1
                current_loc = plus[0]
                for i_old in range( 1, size ):
                    p = plus[ i_old ]
                    if p == current_loc:
                        n += 1
                    else:
                        current_loc = p
                        n = 1
                    if n > maxint:
                        dup_plus [ i_dup ] = p
                        i_dup += 1
                    else:
                        new_plus[ i_new ] = p
                        i_new += 1           
                new_plus.resize( i_new )
                dup_plus.resize( i_dup )
                self.total += i_new
                self.dup_total += i_dup
                self.__pointer[k][0] = i_new
                self.__dup_pointer[k][0] = i_dup
                # unnecessary shape calls
#                self.total +=  new_plus.shape[0]
#                dups.total += dup_plus.shape[0]
#                self.__pointer[k][0] = new_plus.shape[0]
#                dups.__pointer[k][0] = dup_plus.shape[0]
                # free memory?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                plus.resize( 100000, refcheck=False )
                plus.resize( 0, refcheck=False )
                # hope there would be no mem leak...

            # - strand
            i_new = 0
            i_dup = 0
            minus = self.__locations[k][1]
            size = minus.shape[0]
            if len(minus) < 1:
                new_minus = minus         # do nothing
            else:
                new_minus = np.zeros( self.__pointer[k][1],dtype='int32' )
                dup_minus = np.zeros( self.__pointer[k][1],dtype='int32' )
                new_minus[ i_new ] = minus[ i_new ] # first item
                i_new += 1
                current_loc = minus[0]
                n = 1
                for i_old in range( 1, size ):
                    p = minus[ i_old ]
                    if p == current_loc:
                        n += 1
                    else:
                        current_loc = p
                        n = 1
                    if n > maxint:
                        dup_minus [ i_dup ] = p
                        i_dup += 1
                    else:
                        new_minus[ i_new ] = p
                        i_new += 1                        
                new_minus.resize( i_new , refcheck = False) 
                dup_minus.resize( i_dup , refcheck = False) 
                # shape calls unnecessary                      
                self.total +=  i_new
                self.dup_total +=  i_dup
                self.__pointer[k][1] = i_new                
                self.__dup_pointer[k][1] = i_dup                
#                self.total +=  new_minus.shape[0]
#                dups.total +=  dup_minus.shape[0]
#                self.__pointer[k][1] = new_minus.shape[0]                
#                dups.__pointer[k][1] = dup_minus.shape[0]                
                # free memory ?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                minus.resize( 100000, refcheck=False )
                minus.resize( 0, refcheck=False )
                # hope there would be no mem leak...                
            
            self.__locations[k]=[new_plus, new_minus]
            self.__dup_locations[k]=[dup_plus, dup_minus]
        return

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef filter_dup ( self, int32_t maxnum = -1):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        cdef:
            int p, m, n, current_loc, i_chrom
            # index for old array, and index for new one
            unsigned long i_old, i_new, size, new_size 
            str k
            np.ndarray plus, new_plus, dup_plus, minus, new_minus, dup_minus

        if maxnum < 0: return           # do nothing

        if not self.__sorted:
            self.sort()

        self.total = 0

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            k = chrnames[ i_chrom ]
            # + strand
            i_new = 0
            plus = self.__locations[k][0]
            size = plus.shape[0]
            if len(plus) < 1:
                new_plus = plus         # do nothing
            else:
                new_plus = np.zeros( self.__pointer[k][0],dtype='int32' )
                new_plus[ i_new ] = plus[ i_new ] # first item
                i_new += 1
                n = 1                # the number of tags in the current location
                current_loc = plus[0]
                for i_old in range( 1, size ):
                    p = plus[ i_old ]
                    if p == current_loc:
                        n += 1
                        if n <= maxnum:
                            new_plus[ i_new ] = p
                            i_new += 1
                        else:
                            logging.debug("Duplicate reads found at %s:%d at + strand" % (k,p) )
                    else:
                        current_loc = p
                        new_plus[ i_new ] = p
                        i_new += 1                        
                        n = 1
                new_plus.resize( i_new )
                self.total +=  i_new
                self.__pointer[k][0] = i_new
#                self.total +=  new_plus.shape[0]
#                self.__pointer[k][0] = new_plus.shape[0]
                # free memory?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                plus.resize( 100000, refcheck=False )
                plus.resize( 0, refcheck=False )
                # hope there would be no mem leak...

            # - strand
            i_new = 0
            minus = self.__locations[k][1]
            size = minus.shape[0]
            if len(minus) < 1:
                new_minus = minus         # do nothing
            else:
                new_minus = np.zeros( self.__pointer[k][1],dtype='int32' )
                new_minus[ i_new ] = minus[ i_new ] # first item
                i_new += 1
                n = 1                # the number of tags in the current location
                current_loc = minus[0]
                for i_old in range( 1, size ):
                    p = minus[ i_old ]
                    if p == current_loc:
                        n += 1
                        if n <= maxnum:
                            new_minus[ i_new ] = p
                            i_new += 1
                        else:
                            logging.debug("Duplicate reads found at %s:%d at + strand" % (k,p) )
                    else:
                        current_loc = p
                        new_minus[ i_new ] = p
                        i_new += 1                        
                        n = 1
                new_minus.resize( i_new )
                self.total +=  i_new
                self.__pointer[k][1] = i_new                
#                self.total +=  new_minus.shape[0]
#                self.__pointer[k][1] = new_minus.shape[0]                
                # free memory ?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                minus.resize( 100000, refcheck=False )
                minus.resize( 0, refcheck=False )
                # hope there would be no mem leak...                
            
            self.__locations[k]=[new_plus,new_minus]
        return self

    cpdef sample_percent (self, float percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            int32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            str key
        
        self.total = 0

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            key = chrnames[ i_chrom ]
        
            num = <int32_t>round(self.__locations[key][0].shape[0] * percent, 5 )
            np.random.shuffle( self.__locations[key][0] )
            self.__locations[key][0].resize( num )
            self.__locations[key][0].sort()
            self.__pointer[key][0] = self.__locations[key][0].shape[0]

            num = <int32_t>round(self.__locations[key][1].shape[0] * percent, 5 )
            np.random.shuffle( self.__locations[key][1] )
            self.__locations[key][1].resize( num )
            self.__locations[key][1].sort()
            self.__pointer[key][1] = self.__locations[key][1].shape[0]            
            
            self.total += self.__pointer[key][0] + self.__pointer[key][1]
        return

    cpdef sample_num (self, uint64_t samplesize):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            float percent

        percent = float(samplesize)/self.total
        self.sample_percent ( percent )
        return

    cpdef print_to_bed (self, fhd=None):
        """Output FWTrackIII to BED format files. If fhd is given,
        write to a file, otherwise, output to standard output.
        
        """
        cdef:
            int32_t i, i_chrom, p
            str k
        
        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, file)
        assert self.fw > 0, "FWTrackIII object .fw should be set larger than 0!"

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            k = chrnames[ i_chrom ]

            plus = self.__locations[k][0]

            for i in range(plus.shape[0]):
                p = plus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,p,p+self.fw,"+") )

            minus = self.__locations[k][1]
            
            for i in range(minus.shape[0]):
                p = minus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,p-self.fw,p,"-") )
        return
    
