# cython: profile=True
# Time-stamp: <2012-05-02 01:07:31 Tao Liu>

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
    cdef public dict __locations
    cdef public dict __pointer
    cdef public bool __sorted
    cdef public long total
    cdef public dict __dup_locations
    cdef public dict __dup_pointer
    cdef public bool __dup_sorted
    cdef public unsigned long dup_total
    cdef public object annotation
    cdef public object rlengths
    cdef public object dups
    cdef public int fw
    
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
        self.rlengths = None


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

    cpdef __expand__ ( self, np.ndarray arr ):
        arr.resize( arr.size + BUFFER_SIZE, refcheck = False )
        return

    def finalize ( self ):
        """ Resize np arrays for 5' positions and sort them in place """
        
        cdef int32_t i
        cdef str c
        
        self.total+=0 # ??

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
        """Total sequenced length = total number of tags * width of tag		
        """
        return self.total*self.fw

    def sort ( self ):
        """Naive sorting for locations.
        
        """
        cdef int32_t i
        cdef str c

        chrnames = self.get_chr_names()

        for i in range(len(chrnames)):
            c = chrnames[i]
            self.__locations[c][0].sort()
            self.__locations[c][1].sort()

        self.__sorted = True

    @cython.boundscheck(False) # do not check that np indices are valid
    def separate_dups( self, maxint = 1 ):
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
    def filter_dup ( self, int32_t maxnum = -1):
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
        return

    def sample_percent (self, float percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        cdef int32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
        cdef str key
        
        self.total = 0

        chrnames = self.get_chr_names()
        
        for i_chrom in range( len(chrnames) ):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            
            key = chrnames[ i_chrom ]
        
            num = <int32_t>round(self.__locations[key][0].shape[0] * percent, 2 )
            np.random.shuffle( self.__locations[key][0] )
            self.__locations[key][0].resize( num )
            self.__locations[key][0].sort()
            self.__pointer[key][0] = self.__locations[key][0].shape[0]

            num = <int32_t>round(self.__locations[key][1].shape[0] * percent, 2 )
            np.random.shuffle( self.__locations[key][1] )
            self.__locations[key][1].resize( num )
            self.__locations[key][1].sort()
            self.__pointer[key][1] = self.__locations[key][1].shape[0]            
            
            self.total += self.__pointer[key][0] + self.__pointer[key][1]
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
        cdef int32_t i, i_chrom, p
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

            plus = self.__locations[k][0]

            for i in range(plus.shape[0]):
                p = plus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,p,p+self.fw,"+") )

            minus = self.__locations[k][1]
            
            for i in range(minus.shape[0]):
                p = minus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,p-self.fw,p,"-") )
        return
    
#class FWTrackII:
#    """Fixed Width Locations Track class II along the whole genome
#    (commonly with the same annotation type), which are stored in a
#    dict.
#
#    Locations are stored and organized by sequence names (chr names) in a
#    dict. They can be sorted by calling self.sort() function.
#    """
#    def __init__ (self, int fw=0, char * anno=""):
#        """fw is the fixed-width for all locations.
#        
#        """
#        self.fw = fw
#        self.__locations = {}           # locations
#        self.__sorted = False
#        self.total = 0                  # total tags
#        self.annotation = anno          # need to be figured out
#
#    def add_loc (self, str chromosome, int fiveendpos, int strand):
#        """Add a location to the list according to the sequence name.
#        
#        chromosome -- mostly the chromosome name
#        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
#        strand     -- 0: plus, 1: minus
#        """
#        if not self.__locations.has_key(chromosome):
#            self.__locations[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
#            #self.__locations[chromosome] = [ plus , minus] # for (+strand, -strand)
#        self.__locations[chromosome][strand].append(fiveendpos)
#        self.total+=1
#
#    def finalize( self ):
#        return
#
#    def get_locations_by_chr (self, str chromosome):
#        """Return a tuple of two lists of locations for certain chromosome.
#
#        """
#        if self.__locations.has_key(chromosome):
#            return self.__locations[chromosome]
#        else:
#            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))
#
#    def get_chr_names ( self ):
#        """Return all the chromosome names stored in this track object.
#        """
#        l = self.__locations.keys()
#        l.sort()
#        return l
#
#    def length ( self ):
#        """Total sequenced length = total number of tags * width of tag		
#        """
#        return self.total*self.fw
#
#    def sort ( self ):
#        """Naive sorting for locations.
#        
#        """
#        cdef str k
#
#        for k in self.__locations.keys():
#            (tmparrayplus,tmparrayminus) = self.get_locations_by_chr(k)
#            self.__locations[k][0] = sorted(tmparrayplus)
#            if len(tmparrayplus) < 1:
#                logging.warning("NO records for chromosome %s, plus strand!" % (k))
#            self.__locations[k][1] = sorted(tmparrayminus)
#            if len(tmparrayminus) < 1:            
#                logging.warning("NO records for chromosome %s, minus strand!" % (k))
#        self.__sorted = True
#
#    def filter_dup ( self, maxnum ):
#        """Filter the duplicated reads.
#
#        Run it right after you add all data into this object.
#        """
#        cdef int p, m, n, current_loc
#        cdef str k
#
#        if maxnum < 0: return # condition to return if not filtering
#        if not self.__sorted:
#            self.sort()
#        self.total = 0
#        for k in self.__locations.keys(): # for each chromosome
#            # + strand
#            plus = self.__locations[k][0]
#            if len(plus) <1:
#                new_plus = []
#            else:
#                new_plus = array(BYTE4,[plus[0]])
#                pappend = new_plus.append
#                n = 1                # the number of tags in the current location
#                current_loc = plus[0]
#                for p in plus[1:]:
#                    if p == current_loc:
#                        n += 1
#                        if n <= maxnum:
#                            pappend(p)
#                        else:
#                            logging.debug("Duplicate reads found at %s:%d at + strand" % (k,p) )
#                    else:
#                        current_loc = p
#                        pappend(p)
#                        n = 1
#                self.total +=  len(new_plus)
#
#            # - strand
#            minus = self.__locations[k][1]
#            if len(minus) <1:
#                new_minus = []
#            else:
#                new_minus = array(BYTE4,[minus[0]])
#                mappend = new_minus.append
#                n = 1                # the number of tags in the current location
#                current_loc = minus[0]
#                for p in minus[1:]:
#                    if p == current_loc:
#                        n += 1
#                        if n <= maxnum:
#                            mappend(p)
#                        else:
#                            logging.debug("Duplicate reads found at %s:%d at - strand" % (k,p) )                            
#                    else:
#                        current_loc = p
#                        mappend(p)
#                        n = 1
#                self.total +=  len(new_minus)
#            self.__locations[k]=[new_plus,new_minus]
#
#    def merge_plus_minus_locations_naive ( self ):
#        """Merge plus and minus strand locations
#        
#        """
#        cdef str chrom
#        
#        for chrom in self.__locations.keys():
#            #(plus_tags,minus_tags) = self.__locations[chrom]
#            self.__locations[chrom][0].extend(self.__locations[chrom][1])
#            self.__locations[chrom][0] = sorted(self.__locations[chrom][0])
#            self.__locations[chrom][1] = []
#
#    def merge_plus_minus_locations ( self ):
#        """Merge plus and minus strand locations.
#
#        Tao: Amazingly, this function for merging two sorted lists is
#        slower than merge_plus_minus_locations_naive which only
#        concatenate the two lists then sort it again! I am so discouraged!
#        """
#        cdef long ip, im, lenp, lenm
#        cdef str chrom
#        
#        if not self.__sorted:
#            self.sort()
#        for chrom in self.__locations.keys():
#            (plus_tags,minus_tags) = self.__locations[chrom]
#            new_plus_tags = array(BYTE4,[])
#            ip = 0
#            im = 0
#            lenp = len(plus_tags)
#            lenm = len(minus_tags)
#            while ip < lenp and im < lenm:
#                if plus_tags[ip] < minus_tags[im]:
#                    new_plus_tags.append(plus_tags[ip])
#                    ip += 1
#                else:
#                    new_plus_tags.append(minus_tags[im])
#                    im += 1
#            if im < lenm:
#                # add rest of minus tags
#                new_plus_tags.extend(minus_tags[im:])
#            if ip < lenp:
#                # add rest of plus tags
#                new_plus_tags.extend(plus_tags[ip:])
#                    
#            self.__locations[chrom] = [new_plus_tags,[]]
#            self.total += len(new_plus_tags)
#
#		
#    def sample_percent (self, double percent):
#        """Sample the tags for a given percentage.
#
#        Warning: the current object is changed!
#        """
#        cdef long num
#        cdef str key
#        
#        self.total = 0
#        for key in self.__locations.keys():
#            num = long( len( self.__locations[key][0] ) * percent )
#            self.__locations[key][0] = array( BYTE4,
#                                              sorted( random_sample( self.__locations[key][0], num ) ) )
#            num = long( len( self.__locations[key][1] ) * percent )
#            self.__locations[key][1] = array( BYTE4,
#                                              sorted( random_sample( self.__locations[key][1], num ) ) )
#            self.total += len( self.__locations[key][0] ) + len( self.__locations[key][1] )
#
#    def sample_num (self, long samplesize):
#        """Sample the tags for a given percentage.
#
#        Warning: the current object is changed!
#        """
#        cdef double percent
#        cdef long num
#        cdef str key
#
#        percent = double(samplesize)/self.total
#        self.total = 0
#        for key in self.__locations.keys():
#            num = long( len( self.__locations[key][0] ) * percent )
#            self.__locations[key][0] = array( BYTE4,
#                                              sorted( random_sample( self.__locations[key][0], num ) ) )
#            num = long( len( self.__locations[key][1] ) * percent )
#            self.__locations[key][1] = array( BYTE4,
#                                              sorted( random_sample( self.__locations[key][1], num ) ) )
#            self.total += len( self.__locations[key][0] ) + len( self.__locations[key][1] )
#            
#    def __str__ ( self ):
#        return self.__to_wiggle()
#        
#    def __to_wiggle ( self ):
#        """Use a lot of memory!
#        
#        """
#        cdef int i
#        cdef str k, t
#        
#        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
#        for k in self.__locations.keys():
#            if self.__locations[k][0]:
#                t += "variableStep chrom=%s span=%d strand=0\n" % (k,self.fw)
#                for i in self.__locations[k][0]:
#                    t += "%d\t1\n" % i
#            if self.__locations[k][1]:
#                t += "variableStep chrom=%s span=%d strand=1\n" % (k,self.fw)
#                for i in self.__locations[k][1]:
#                    t += "%d\t1\n" % i
#        return t
#
#    def print_to_bed (self, fhd=None):
#        """Output FWTrackII to BED format files. If fhd is given,
#        write to a file, otherwise, output to standard output.
#        
#        """
#        cdef int i
#        cdef str k
#        
#        if not fhd:
#            fhd = sys.stdout
#        assert isinstance(fhd, file)
#        assert self.fw > 0, "FWTrackII object .fw should be set larger than 0!"
#        for k in self.__locations.keys():
#            if self.__locations[k][0]:
#                for i in self.__locations[k][0]:
#                    fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,i,int(i+self.fw),"+") )
#            if self.__locations[k][1]:
#                for i in self.__locations[k][1]:
#                    fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,int(i-self.fw),i,"-") )
#        return
#
#
#    #def __del__ ( self ):
#    #    for k in self.__locations.keys():
#    #        self.__locations[k][0] = None
#    #        self.__locations[k][1] = None
#    
#            
