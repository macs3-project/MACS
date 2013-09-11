# Time-stamp: <2013-09-11 17:42:07 Tao Liu>

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
from collections import Counter

from MACS2.Constants import *
from MACS2.cSignal import *
from MACS2.IO.cPeakIO import PeakIO
from MACS2.cPileup import se_all_in_one_pileup, max_over_two_pv_array

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
        bool __destroyed
        bool __dup_separated
        dict rlengths
        public long buffer_size
        public long total
        public unsigned long dup_total
        public object annotation
        public object dups
        public int fw
    
    def __init__ (self, int32_t fw=0, char * anno="", long buffer_size = 100000 ):
        """fw is the fixed-width for all locations.
        
        """
        self.fw = fw
        self.__locations = {}    # location pairs
        self.__pointer = {}      # location pairs
        self.__dup_locations = {}    # location pairs
        self.__dup_pointer = {}      # location pairs
        self.__sorted = False
        self.__dup_sorted = False
        self.__dup_separated = False
        self.total = 0           # total tags
        self.dup_total = 0           # total tags
        self.annotation = anno   # need to be figured out
        self.rlengths = {}       # lengths of reference sequences, e.g. each chromosome in a genome
        self.buffer_size = buffer_size
        self.__destroyed = False

    cpdef destroy ( self ):
        """Destroy this object and release mem.
        """
        cdef:
            set chrs
            str chromosome
            
        chrs = set(self.get_chr_names())
        for chromosome in chrs:
            if self.__locations.has_key(chromosome):
                self.__locations[chromosome][0].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome][0].resize( 0, refcheck=False )
                self.__locations[chromosome][1].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome][1].resize( 0, refcheck=False )
                self.__locations[chromosome] = [None, None]
                self.__locations.pop(chromosome)
            if self.__dup_locations.has_key(chromosome):
                self.__dup_locations[chromosome][0].resize( self.buffer_size, refcheck=False )
                self.__dup_locations[chromosome][0].resize( 0, refcheck=False )
                self.__dup_locations[chromosome][1].resize( self.buffer_size, refcheck=False )
                self.__dup_locations[chromosome][1].resize( 0, refcheck=False )
                self.__dup_locations[chromosome] = [None, None]
                self.__dup_locations.pop(chromosome)
        self.__destroyed = True
        return True


    cpdef add_loc ( self, str chromosome, int32_t fiveendpos, int strand ):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = [ np.zeros(self.buffer_size, dtype='int32'), np.zeros(self.buffer_size, dtype='int32') ] # [plus,minus strand]
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
        arr.resize( arr.size + self.buffer_size, refcheck = False )
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

    cpdef list get_chr_names ( self ):
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

    @cython.boundscheck(False)
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
            dup_plus = np.zeros( self.__pointer[k][0],dtype='int32' )
            if len(plus) < 1:
                new_plus = plus         # do nothing
            else:
                new_plus = np.zeros( self.__pointer[k][0],dtype='int32' )
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
                plus.resize( self.buffer_size, refcheck=False )
                plus.resize( 0, refcheck=False )
                # hope there would be no mem leak...

            # - strand
            i_new = 0
            i_dup = 0
            minus = self.__locations[k][1]
            size = minus.shape[0]
            dup_minus = np.zeros( self.__pointer[k][1],dtype='int32' )
            if len(minus) < 1:
                new_minus = minus         # do nothing
            else:
                new_minus = np.zeros( self.__pointer[k][1],dtype='int32' )
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
                minus.resize( self.buffer_size, refcheck=False )
                minus.resize( 0, refcheck=False )
                # hope there would be no mem leak...                
            
            self.__locations[k]=[new_plus, new_minus]
            self.__dup_locations[k]=[dup_plus, dup_minus]

        self.__dup_separated = True
        
        return

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef addback_dups( self ):
        """Add back the duplicate reads stored in self.__dup_locations to self.__locations
        """
        cdef:
            int p, m, n, current_loc, i_chrom
            unsigned long i_old, i_new          # index for old array, and index for new one
            unsigned long i_dup, size, new_size, dup_size
            str k
            np.ndarray plus, new_plus, dup_plus, minus, new_minus, dup_minus

        if not self.__sorted:
            self.sort()

        assert self.__dup_separated == True, "need to run separate_dups first."
        self.total = 0

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            k = chrnames[ i_chrom ]
            plus = self.__locations[k][0]
            dup_plus = self.__dup_locations[k][0]
            minus = self.__locations[k][1]
            dup_minus = self.__dup_locations[k][1]
            
            # concatenate
            new_plus = np.concatenate((plus, dup_plus))
            new_minus= np.concatenate((minus, dup_minus))

            # clean old data
            plus.resize( self.buffer_size, refcheck=False )
            plus.resize( 0, refcheck=False )
            dup_plus.resize( self.buffer_size, refcheck=False )
            dup_plus.resize( 0, refcheck=False )
            minus.resize( self.buffer_size, refcheck=False )
            minus.resize( 0, refcheck=False )
            dup_minus.resize( self.buffer_size, refcheck=False )
            dup_minus.resize( 0, refcheck=False )            

            # sort then assign
            new_plus.sort()
            new_minus.sort()
            self.__locations[k][0] = new_plus
            self.__locations[k][1] = new_minus
            self.__dup_locations[k][0] = None
            self.__dup_locations[k][1] = None            

            self.__pointer[k][0] = plus.shape[0]
            self.__pointer[k][1] = minus.shape[0]
            self.__dup_pointer[k][0] = 0
            self.__dup_pointer[k][1] = 0            
            self.total +=  plus.shape[0] + minus.shape[0]

        self.dup_total =  0
        self.__dup_separated = False
        return

    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef filter_dup ( self, int32_t maxnum = -1):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.

        Note, this function will *throw out* duplicates
        permenantly. If you want to keep them, use separate_dups
        instead.
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
                plus.resize( self.buffer_size, refcheck=False )
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
                minus.resize( self.buffer_size, refcheck=False )
                minus.resize( 0, refcheck=False )
                # hope there would be no mem leak...                
            
            self.__locations[k]=[new_plus,new_minus]
        return self

    cpdef sample_percent (self, float percent, int seed = -1 ):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            int32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            str key
        
        self.total = 0

        chrnames = self.get_chr_names()
        
        if seed >= 0:
            np.random.seed(seed)

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

    cpdef sample_num (self, uint64_t samplesize, int seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef:
            float percent

        percent = float(samplesize)/self.total
        self.sample_percent ( percent, seed )
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
    
    cpdef tuple extract_region_tags ( self, str chromosome, int32_t startpos, int32_t endpos ):
        cdef:
            int32_t i, pos
            np.ndarray rt_plus, rt_minus
            list temp

        if not self.__sorted: self.sort()
        
        chrnames = self.get_chr_names()
        assert chromosome in chrnames, "chromosome %s can't be found in the FWTrackIII object." % chromosome
        
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

    cpdef compute_region_tags_from_peaks ( self, peaks, func, int window_size = 100, float cutoff = 5 ):
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
            np.ndarray plus, minus, rt_plus, rt_minus
            str chrom, name
            list temp, retval, pchrnames

        pchrnames = peaks.get_chr_names()
        retval = []

        # this object should be sorted
        if not self.__sorted: self.sort()
        # PeakIO object should be sorted
        peaks.sort()
        
        chrnames = self.get_chr_names()

        for chrom in pchrnames:
            assert chrom in chrnames, "chromosome %s can't be found in the FWTrackIII object." % chrom
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
                rt_plus = np.array(temp)
                
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
                rt_minus = np.array(temp)

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

    cpdef refine_peak_from_tags_distribution ( self, peaks, int window_size = 100, float cutoff = 5 ):
        """Extract tags in peak, then apply func on extracted tags.
        
        peaks: redefined regions to extract raw tags in PeakIO type: check cPeakIO.pyx.

        window_size: this will be passed to func.

        cutoff: this will be passed to func.

        func needs the fixed number of parameters, so it's not flexible. Here is an example:

        wtd_find_summit(chrom, plus, minus, peak_start, peak_end, name , window_size, cutoff):

        """
        
        cdef:
            int32_t m, i, j, pre_i, pre_j, pos, startpos, endpos #, n_peaks
            np.ndarray plus, minus, rt_plus, rt_minus
            str chrom #, peak_name
            list temp, retval, pchrnames, cpeaks
            np.ndarray adjusted_summits, passflags

        pchrnames = sorted(peaks.get_chr_names())
        retval = []

        # this object should be sorted
        if not self.__sorted: self.sort()
        # PeakIO object should be sorted
        peaks.sort()
        
        chrnames = self.get_chr_names()

        #n_peaks = 1
        ret_peaks = PeakIO()
        
        for chrom in pchrnames:
            assert chrom in chrnames, "chromosome %s can't be found in the FWTrackIII object. %s" % (chrom, str(chrnames))
            (plus, minus) = self.__locations[chrom]
            cpeaks = peaks.get_data_from_chrom(chrom)
            #ret_peaks.peaks[chrom] = []
            #npeaks = ret_peaks.peaks[chrom]
            
            prev_i = 0
            prev_j = 0
            for m in range(len(cpeaks)):
                thispeak = cpeaks[m]
                startpos = thispeak["start"] - window_size
                endpos   = thispeak["end"] + window_size
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
                rt_plus = np.array(temp)
                
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
                rt_minus = np.array(temp)

                #peak_name = name + "_" + str(n_peaks)
                (adjusted_summits, passflags) = wtd_find_summit(chrom, rt_plus, rt_minus, startpos, endpos, window_size, cutoff)
                # those local maxima above cutoff will be defined as good summits
                for i in range(len(adjusted_summits)):
                    adjusted_summit = adjusted_summits[i]
                    passflag = passflags[i]
                    if passflag:
                        tmppeak = copy(thispeak)
                        tmppeak["summit"] = adjusted_summit
                        ret_peaks.add_PeakContent(chrom, tmppeak)
                    
                #thispeak["summit"] = adjusted_summit
                #if passflag:
                #    thispeak["name"] = "passed"
                #else:
                #    thispeak["name"] = "failed"
                #retval.append( wtd_find_summit(chrom, rt_plus, rt_minus, startpos, endpos, peak_name, window_size, cutoff) )
                #n_peaks += 1
                
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
        return ret_peaks

    cpdef pileup_a_chromosome ( self, str chrom, list ds, list scale_factor_s, float baseline_value = 0.0, bint directional = True, bint halfextension = True ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.
        
        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        directional    : if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
        halfextension  : only make a fragment of d/2 size centered at fragment center        
        """
        cdef:
            long d
            long five_shift, three_shift  # adjustment to 5' end and 3' end positions to make a fragment
            dict chrlengths = self.get_rlengths ()
            long rlength = chrlengths[chrom]
            object ends
            list five_shift_s = []
            list three_shift_s = []
            list tmp_pileup, prev_pileup

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        # adjust extension length according to 'directional' and 'halfextension' setting.
        for d in ds:
            if directional:
                # only extend to 3' side
                if halfextension:
                    five_shift_s.append(d/-4)  # five shift is used to move cursor towards 5' direction to find the start of fragment
                    three_shift_s.append(d*3/4) # three shift is used to move cursor towards 3' direction to find the end of fragment
                else:
                    five_shift_s.append(0)
                    three_shift_s.append(d)
            else:
                # both sides
                if halfextension:
                    five_shift_s.append(d/4)
                    three_shift_s.append(d/4)
                else:
                    five_shift_s.append(d/2)
                    three_shift_s.append(d - d/2)

        prev_pileup = None

        for i in range(len(ds)):
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]
            
            tmp_pileup = se_all_in_one_pileup ( self.__locations[chrom][0], self.__locations[chrom][1], five_shift, three_shift, rlength, scale_factor, baseline_value )

            if prev_pileup:
                prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup )
            else:
                prev_pileup = tmp_pileup
        return prev_pileup

cdef inline int32_t left_sum ( data, int pos, int width ):
    """
    """
    return sum([data[x] for x in data if x <= pos and x >= pos - width])

cdef inline int32_t right_sum ( data, int pos, int width ):
    """
    """
    return sum([data[x] for x in data if x >= pos and x <= pos + width])

cdef inline int32_t left_forward ( data, int pos, int window_size ):
    return data.get(pos,0) - data.get(pos-window_size, 0)

cdef inline int32_t right_forward ( data, int pos, int window_size ):
    return data.get(pos + window_size, 0) - data.get(pos, 0)

cdef wtd_find_summit(chrom, np.ndarray plus, np.ndarray minus, int32_t search_start, int32_t search_end, int32_t window_size, float cutoff):
    """internal function to be called by refine_peak_from_tags_distribution()

    """
    cdef:
        int32_t i, j, watson_left, watson_right, crick_left, crick_right, wtd_max_pos
        float wtd_max_val
        np.ndarray wtd_list, wtd_other_max_pos, wtd_other_max_val
        
    watson, crick = (Counter(plus), Counter(minus))
    watson_left = left_sum(watson, search_start, window_size)
    crick_left = left_sum(crick, search_start, window_size)
    watson_right = right_sum(watson, search_start, window_size)
    crick_right = right_sum(crick, search_start, window_size)

    wtd_list = np.zeros( search_end - search_start + 1, dtype="float32")
    i = 0
    for j in range(search_start, search_end+1):
        wtd_list[i] = max((2 * (watson_left * crick_right)**0.5 - watson_right - crick_left),0) # minimum score is 0
        watson_left += left_forward(watson, j, window_size)
        watson_right += right_forward(watson, j, window_size)
        crick_left += left_forward(crick, j, window_size)
        crick_right += right_forward(crick, j, window_size)
        i += 1

    #wtd_max_val = max(wtd_list)
    #wtd_max_pos = wtd_list.index(wtd_max_val) + search_start

    # smooth
    #wtd_list = smooth(wtd_list, window="flat") # window size is by default 11.
    #wtd_max_pos = np.where(wtd_list==max(wtd_list))[0][0]
    #wtd_max_val = wtd_list[wtd_max_pos]
    #wtd_max_pos += search_start
    # search for other local maxima
    #wtd_other_max_pos = np.arange(len(wtd_list))[np.r_[False, wtd_list[1:] > wtd_list[:-1]] & np.r_[wtd_list[:-1] > wtd_list[1:], False]]
    #wtd_other_max_val = wtd_list[wtd_other_max_pos]
    #wtd_other_max_pos = wtd_other_max_pos + search_start

    wtd_other_max_pos = maxima(wtd_list, window_size = window_size)
    wtd_other_max_pos = enforce_peakyness( wtd_list, wtd_other_max_pos )
    wtd_other_max_val = wtd_list[wtd_other_max_pos]
    wtd_other_max_pos = wtd_other_max_pos + search_start

    #return (chrom, wtd_max_pos, wtd_max_pos+1, wtd_max_val)

    return (wtd_other_max_pos, wtd_other_max_val > cutoff)

    #if wtd_max_val > cutoff:
    #    return (wtd_max_pos, True)
    #    #return (chrom, wtd_max_pos, wtd_max_pos+1, name+"_R" , wtd_max_val) # 'R'efined
    #else:
    #    return (wtd_max_pos, False)
    #    #return (chrom, wtd_max_pos, wtd_max_pos+1, name+"_F" , wtd_max_val) # 'F'ailed
