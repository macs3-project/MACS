# Time-stamp: <2011-09-06 16:57:34 Tao Liu>

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

from array import array
from random import sample as random_sample
from sys import stdout as sysstdout
from MACS2.Constants import *

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

class FWTrackII:
    """Fixed Width Locations Track class II along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    def __init__ (self,fw=0,anno=""):
        """fw is the fixed-width for all locations.
        
        """
        self.fw = fw
        self.__locations = {}           # locations
        self.__sorted = False
        self.total = 0                  # total tags
        self.annotation = anno   # need to be figured out

    def add_loc (self, char * chromosome, long fiveendpos, int strand):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
            #self.__locations[chromosome] = [ plus , minus] # for (+strand, -strand)
        self.__locations[chromosome][strand].append(fiveendpos)
        self.total+=1

    def get_locations_by_chr (self, chromosome):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if self.__locations.has_key(chromosome):
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored in this track object.
        """
        l = self.__locations.keys()
        l.sort()
        return l

    def length (self):
        """Total sequenced length = total number of tags * width of tag		
        """
        return self.total*self.fw

    def sort (self):
        """Naive sorting for locations.
        
        """
        for k in self.__locations.keys():
            (tmparrayplus,tmparrayminus) = self.get_locations_by_chr(k)
            self.__locations[k][0] = sorted(tmparrayplus)
            if len(tmparrayplus) < 1:
                logging.warning("NO records for chromosome %s, plus strand!" % (k))
            self.__locations[k][1] = sorted(tmparrayminus)
            if len(tmparrayminus) < 1:            
                logging.warning("NO records for chromosome %s, minus strand!" % (k))
        self.__sorted = True

    def filter_dup (self,maxnum):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        if not self.__sorted:
            self.sort()
        self.total = 0
        for k in self.__locations.keys(): # for each chromosome
            # + strand
            plus = self.__locations[k][0]
            if len(plus) <1:
                new_plus = []
            else:
                new_plus = array(BYTE4,[plus[0]])
                pappend = new_plus.append
                n = 1                # the number of tags in the current location
                current_loc = plus[0]
                for p in plus[1:]:
                    if p == current_loc:
                        n += 1
                        if n <= maxnum:
                            pappend(p)
                        else:
                            logging.debug("Duplicate reads found at %s:%d at + strand" % (k,p) )
                    else:
                        current_loc = p
                        pappend(p)
                        n = 1
                self.total +=  len(new_plus)

            # - strand
            minus = self.__locations[k][1]
            if len(minus) <1:
                new_minus = []
            else:
                new_minus = array(BYTE4,[minus[0]])
                mappend = new_minus.append
                n = 1                # the number of tags in the current location
                current_loc = minus[0]
                for p in minus[1:]:
                    if p == current_loc:
                        n += 1
                        if n <= maxnum:
                            mappend(p)
                        else:
                            logging.debug("Duplicate reads found at %s:%d at - strand" % (k,p) )                            
                    else:
                        current_loc = p
                        mappend(p)
                        n = 1
                self.total +=  len(new_minus)
            self.__locations[k]=[new_plus,new_minus]

    def merge_plus_minus_locations_naive (self):
        """Merge plus and minus strand locations
        
        """
        for chrom in self.__locations.keys():
            #(plus_tags,minus_tags) = self.__locations[chrom]
            self.__locations[chrom][0].extend(self.__locations[chrom][1])
            self.__locations[chrom][0] = sorted(self.__locations[chrom][0])
            self.__locations[chrom][1] = []

    def merge_plus_minus_locations (self):
        """Merge plus and minus strand locations.

        Tao: Amazingly, this function for merging two sorted lists is
        slower than merge_plus_minus_locations_naive which only
        concatenate the two lists then sort it again! I am so discouraged!
        """
        if not self.__sorted:
            self.sort()
        for chrom in self.__locations.keys():
            (plus_tags,minus_tags) = self.__locations[chrom]
            new_plus_tags = array(BYTE4,[])
            ip = 0
            im = 0
            lenp = len(plus_tags)
            lenm = len(minus_tags)
            while ip < lenp and im < lenm:
                if plus_tags[ip] < minus_tags[im]:
                    new_plus_tags.append(plus_tags[ip])
                    ip += 1
                else:
                    new_plus_tags.append(minus_tags[im])
                    im += 1
            if im < lenm:
                # add rest of minus tags
                new_plus_tags.extend(minus_tags[im:])
            if ip < lenp:
                # add rest of plus tags
                new_plus_tags.extend(plus_tags[ip:])
                    
            self.__locations[chrom] = [new_plus_tags,[]]
            self.total += len(new_plus_tags)

		
    def sample (self, percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        self.total = 0
        for key in self.__locations.keys():
            num = int(len(self.__locations[key][0])*percent)
            self.__locations[key][0]=array(BYTE4,sorted(random_sample(self.__locations[key][0],num)))
            num = int(len(self.__locations[key][1])*percent)
            self.__locations[key][1]=array(BYTE4,sorted(random_sample(self.__locations[key][1],num)))
            self.total += len(self.__locations[key][0]) + len(self.__locations[key][1])
            
    def __str__ (self):
        return self.__to_wiggle()
        
    def __to_wiggle (self):
        """Use a lot of memory!
        
        """
        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
        for k in self.__locations.keys():
            if self.__locations[k][0]:
                t += "variableStep chrom=%s span=%d strand=0\n" % (k,self.fw)
                for i in self.__locations[k][0]:
                    t += "%d\t1\n" % i
            if self.__locations[k][1]:
                t += "variableStep chrom=%s span=%d strand=1\n" % (k,self.fw)
                for i in self.__locations[k][1]:
                    t += "%d\t1\n" % i
        return t

    def print_to_bed (self, fhd=None):
        """Output FWTrackII to BED format files. If fhd is given,
        write to a file, otherwise, output to standard output.
        
        """
        if not fhd:
            fhd = sysstdout
        assert isinstance(fhd, file)
        assert self.fw > 0, "FWTrackII object .fw should be set larger than 0!"
        for k in self.__locations.keys():
            if self.__locations[k][0]:
                for i in self.__locations[k][0]:
                    fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,i,int(i+self.fw),"+") )
            if self.__locations[k][1]:
                for i in self.__locations[k][1]:
                    fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k,int(i-self.fw),i,"-") )
        return

