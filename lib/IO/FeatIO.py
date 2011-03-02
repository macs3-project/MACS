# Time-stamp: <2011-03-02 17:28:15 Tao Liu>

"""Module for Feature IO classes.

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
import re
import sys
import logging
import struct
from array import array
from random import sample as random_sample
from operator import itemgetter
from math import sqrt

from MACS14.Constants import *

# ------------------------------------
# constants
# ------------------------------------
__version__ = "FeatIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "PeakIO, FWTrackI, TrackI, and WigTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------


# ------------------------------------
# Classes
# ------------------------------------

class PeakIO:
    """IO for peak information.

    """
    def __init__ (self):
        self.peaks = {}
    
    def add (self, chromosome, start, end, summit=None, 
             peak_height=None, number_tags=None, 
             pvalue=None, fold_enrichment=None, fdr=None):
        """items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append((start,end,end-start,summit,
                                       peak_height,number_tags,
                                       pvalue,fold_enrichment,fdr))

    def filter_pvalue (self, pvalue_cut ):
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fc_low and p[7]<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
    
    def ave_fold_enrichment (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        t = 0
        for chrom in chrs:
            x += len(peaks[chrom])
            for p in peaks[chrom]:
                t+=p[7]
        return t/x

    def max_fold_enrichment (self):
        """Return the maximum fc.
        
        """
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            if peaks[chrom]:
                m = max([i[7] for i in peaks[chrom]])
                if m>x:
                    x=m
        return x
        
    
    def tobed (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (chrom,peak[0],peak[1])
        return text

    def towig (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (peak[0],peak[1])
        return text
        
    def init_from_dict (self, data):
        """Initialize the data from a dictionary. Improper assignment
        will damage the data structure.
        
        """
        self.peaks = {}
        chrs = data.keys()
        chrs.sort()
        for chrom in chrs:
            self.peaks[chrom]=[]
            a = self.peaks[chrom].append
            for i in data[chrom]:
                a(i)

    def merge_overlap ( self ):
        """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            n_append = new_peaks[chrom].append
            prev_peak = None
            peaks_chr = peaks[chrom]
            for i in xrange(len(peaks_chr)):
                if not prev_peak:
                    prev_peak = peaks_chr[i]
                    continue
                else:
                    if peaks_chr[i][0] <= prev_peak[1]:
                        s_new_peak = prev_peak[0]
                        e_new_peak = peaks_chr[i][1]
                        l_new_peak = e_new_peak-s_new_peak
                        if peaks_chr[i][4] > prev_peak[4]:
                            summit_new_peak = peaks_chr[i][3]
                            h_new_peak = peaks_chr[i][4]
                        else:
                            summit_new_peak = prev_peak[3]
                            h_new_peak = prev_peak[4]
                        prev_peak = (s_new_peak,e_new_peak,l_new_peak,summit_new_peak,h_new_peak,peaks_chr[i][5]+prev_peak[5])
                    else:
                        n_append(prev_peak)
                        prev_peak = peaks_chr[i]
            if prev_peak:
                n_append(prev_peak)
        del peaks
        self.peaks = new_peaks
        return True

    def overlap_with_other_peaks (self, peaks2, cover=0):
        """Peaks2 is a PeakIO object or dictionary with can be
        initialzed as a PeakIO. check __init__ for PeakIO for detail.

        return how many peaks are intersected by peaks2 by percentage
        coverage on peaks2(if 50%, cover = 0.5).
        """
        peaks1 = self.peaks
        if isinstance(peaks2,PeakIO):
            peaks2 = peaks2.peaks
        total_num = 0
        chrs1 = peaks1.keys()
        chrs2 = peaks2.keys()
        for k in chrs1:
            if not chrs2.count(k):
                continue
            rl1_k = iter(peaks1[k])
            rl2_k = iter(peaks2[k])
            tmp_n = False
            try:
                r1 = rl1_k.next()
                r2 = rl2_k.next()
                while (True):
                    if r2[0] < r1[1] and r1[0] < r2[1]:
                        a = sorted([r1[0],r1[1],r2[0],r2[1]])
                        if float(a[2]-a[1]+1)/r2[2] > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                    if r1[1] < r2[1]:
                        r1 = rl1_k.next()
                        tmp_n = False
                    else:
                        r2 = rl2_k.next()
            except StopIteration:
                continue
        return total_num
        

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

    def add_loc (self, chromosome, fiveendpos, strand):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
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

class FWTrackI:
    """Fixed Width Ranges along the whole genome (commonly with the
    same annotation type), which are stored in a dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.

    Example:
       >>> tabfile = TabFile('tag.bed',format='bed',mode='r')
       >>> track = FWTrackI()
       >>> for (chromosome,rg) in tabfile:
       ...    track.add_location(chromosome,rg)
       >>> track.get_locations_by_chr['chr1'] # all locations in chr1 
    """
    def __init__ (self,fw=0,anno=""):
        """fw is the fixed-width for all locations
        """
        self.fw = fw
        self.__locations = {}           # locations
        self.__counts = {}              # number of tags at the same location
        self.__well_merged = False
        self.total = 0					# total tags
        self.total_unique = 0		# total unique tags
        self.annotation = anno   # need to be figured out

    def add_loc (self, chromosome, fiveendpos, strand):
        """Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__locations.has_key(chromosome):
            self.__locations[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
            self.__counts[chromosome] = [array(UBYTE2,[]),array(UBYTE2,[])] # for (+,-)
        self.__locations[chromosome][strand].append(fiveendpos)
        self.__counts[chromosome][strand].append(1)
        self.total+=1

    def get_locations_by_chr (self, chromosome):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if self.__locations.has_key(chromosome):
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_counts_by_chr (self, chromosome):
        """Return a tuple of two lists of counts for certain chromosome.

        """
        if self.__counts.has_key(chromosome):
            return self.__counts[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_loc_counts_by_chr (self, chromosome):
        """Return a tuple of two lists of (loc,count) for certain chromosome.

        """
        if self.__counts.has_key(chromosome):
            return (zip(self.__locations[chromosome][0],self.__counts[chromosome][0]),
                    zip(self.__locations[chromosome][1],self.__counts[chromosome][1]))
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
            g0 = itemgetter(0)
            g1 = itemgetter(1)
            (tmparrayplus,tmparrayminus) = self.get_loc_counts_by_chr(k)
            sortedtmparrayplus = sorted(tmparrayplus,key=g0)
            self.__locations[k][0] = [g0(x) for x in sortedtmparrayplus]
            self.__counts[k][0] = [g1(x) for x in sortedtmparrayplus]
            sortedtmparrayminus = sorted(tmparrayminus,key=g0)
            self.__locations[k][1] = [g0(x) for x in sortedtmparrayminus]
            self.__counts[k][1] = [g1(x) for x in sortedtmparrayminus]

    def merge_overlap (self):
        """merge the SAME locations. Record the duplicate number in self.__counts{}

        Run it right after you add all data into this object.
        
        *Note: different with the merge_overlap() in TrackI class,
        which merges the overlapped locations.
        """
        self.total = 0
        self.total_unique = 0
        for k in self.__locations.keys(): # for each chromosome
            # + strand
            plus = sorted(self.__locations[k][0])
            if len(plus) <1:
                logging.warning("NO records for chromosome %s, plus strand!" % (k))
                new_plus = []
                new_plus_c = []
            else:
                (new_plus,new_plus_c) = (array(BYTE4,[plus[0]]),array(UBYTE2,[1]))
            
                pappend = new_plus.append
                pcappend = new_plus_c.append
                n = 0                # the position in new list
                for p in plus[1:]:
                    if p == new_plus[n]:
                        try:
                            new_plus_c[n]+=1
                        except OverflowError:
                            logging.warning("> 65535 + strand tags mapped to position %d on chromosome %s!" % (p,k))
                            new_plus_c[n]=65535

                    else:
                        pappend(p)
                        pcappend(1)
                        n += 1
                self.total_unique +=  len(new_plus)
                self.total += sum(new_plus_c)
            # - strand
            minus = sorted(self.__locations[k][1])
            if len(minus) <1:
                logging.warning("NO records for chromosome %s, minus strand!" % (k))
                new_minus = []
                new_minus_c = []
            else:
                (new_minus,new_minus_c) = (array(BYTE4,[minus[0]]),array(UBYTE2,[1]))
            
                mappend = new_minus.append
                mcappend = new_minus_c.append
                n = 0                # the position in new list
                for p in minus[1:]:
                    if p == new_minus[n]:
                        try:
                            new_minus_c[n]+=1
                        except OverflowError:
                            logging.warning("> 65535 - strand tags mapped to position %d on chromosome %s!" % (p,k))
                            new_minus_c[n]=65535
                    else:
                        mappend(p)
                        mcappend(1)
                        n += 1
                self.total_unique +=  len(new_minus)
                self.total += sum(new_minus_c)

            self.__locations[k]=[new_plus,new_minus]
            self.__counts[k]=[new_plus_c,new_minus_c]
            self.__well_merged = True
		
    def merge_plus_minus_locations_w_duplicates (self,maxnum=None):
        """Merge minus strand locations to plus strand. The duplicates
        on a single strand is defined in 'maxnum'. If maxnum is None,
        then keep all the duplicates.
                
        Side effect: Reset the counts. self.total_unique is set to
        None. This object is changed.
        """
        self.total_unique = None
        self.total = 0
        for chrom in self.__locations.keys():
            (plus_tags,minus_tags) = self.__locations[chrom]
            (plus_counts,minus_counts) = self.__counts[chrom]

            new_plus_tags = array(BYTE4,[])
            new_plus_counts = array(UBYTE2,[])

            ip = 0
            im = 0
            lenp = len(plus_tags)
            lenm = len(minus_tags)
            while ip < lenp and im < lenm:
                if plus_tags[ip] < minus_tags[im]:
                    for i in xrange(plus_counts[ip]):
                        if maxnum and i+1>maxnum:
                            break
                        new_plus_tags.append(plus_tags[ip])
                        new_plus_counts.append(1)
                    ip += 1
                else:
                    for i in xrange(minus_counts[im]):
                        if maxnum and i+1>maxnum:
                            break
                        new_plus_tags.append(minus_tags[im])
                        new_plus_counts.append(1)                        
                    im += 1
            for im2 in xrange(im,lenm):
                # add rest of minus tags
                for i in xrange(minus_counts[im2]):
                    if maxnum and i+1>maxnum:
                        break
                    new_plus_tags.append(minus_tags[im2])
                    new_plus_counts.append(1)
                    
            for ip2 in xrange(ip,lenp):
                # add rest of minus tags
                for i in xrange(plus_counts[ip2]):
                    if maxnum and i+1>maxnum:
                        break
                    new_plus_tags.append(plus_tags[ip2])
                    new_plus_counts.append(1)
                    
            self.__locations[chrom] = [new_plus_tags,[]]
            self.__counts[chrom] = [new_plus_counts,[]]
            self.total += len(new_plus_tags)

    def sample (self, percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        
        Side effect: self.total_unique is set to None, and counts are unset.
        """
        self.total = 0
        self.total_unique = None
        for key in self.__locations.keys():
            num = int(len(self.__locations[key][0])*percent)
            self.__locations[key][0]=array(BYTE4,sorted(random_sample(self.__locations[key][0],num)))
            num = int(len(self.__locations[key][1])*percent)
            self.__locations[key][1]=array(BYTE4,sorted(random_sample(self.__locations[key][1],num)))
            self.total += len(self.__locations[key][0]) + len(self.__locations[key][1])
            self.__counts[key] = [[],[]]
            
    def __str__ (self):
        return self.__to_wiggle()
        
    def __to_wiggle (self):
        """Use a lot of memory!
        
        """
        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
        for k in self.__locations.keys():
            (tmparrayplus,tmparrayminus) = self.get_loc_counts_by_chr(k)

            if self.__locations[k][0]:
                t += "variableStep chrom=%s span=%d strand=0\n" % (k,self.fw)
                for (i,j) in tmparrayplus:
                    t += "%d\t%d\n" % (i,j)
            if self.__locations[k][1]:
                t += "variableStep chrom=%s span=%d strand=1\n" % (k,self.fw)
                for (i,j) in tmparrayminus:
                    t += "%d\t%d\n" % (i,j)
        return t


class WigTrackI:
    """Designed only for wig files generated by MACS/pMA2C/MAT(future
    version). The limitation is 'span' parameter for every track must
    be the same.
    
    """
    def __init__ (self):
        self.__data = {}
        self.span=0
        self.maxvalue =-10000
        self.minvalue = 10000

    def add_loc (self,chromosome,pos,value):
        if not self.__data.has_key(chromosome):
            self.__data[chromosome] = [array(BYTE4,[]),array(FBYTE4,[])] # for (pos,value)
        self.__data[chromosome][0].append(pos)
        self.__data[chromosome][1].append(value)
        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value
        
    def sort (self):
        """Naive sorting for tags. After sorting, counts are massed
        up.

        """
        for k in self.__data.keys():
            (p,v) = self.__data[k]
            pv = zip(p,v)
            pv = sorted(pv)
            self.__data[k] = [array(BYTE4,[]),array(FBYTE4,[])]
            pappend = self.__data[k][0].append
            vappend = self.__data[k][1].append
            for (tp,tv) in pv:
                pappend(tp)
                vappend(tv)

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([pos],[value])
        """
        if self.__data.has_key(chromosome):
            return self.__data[chromosome]
        else:
            return None
            #raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.__data.keys())
        return l

    def write_wig (self, fhd, name, shift=0):
        """Write all data to fhd in Wiggle Format.

        shift will be used to shift the coordinates. default: 0
        """
        fhd.write("track type=wiggle_0 name=\"%s\" description=\"%s\"\n" % (name,name))
        chrs = self.get_chr_names()
        for chrom in chrs:
            fhd.write("variableStep chrom=%s span=%d\n" % (chrom,self.span))
            (p,s) = self.__data[chrom]
            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                fhd.write("%d\t%.4f\n" % (pos+shift,score))                

    def filter_score (self, cutoff=0):
        """Filter using a score cutoff. Return a new WigTrackI.
        
        """
        ret = WigTrackI()
        ret.span = self.span
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(FBYTE4,[])]
            (np,ns) = ret.__data[chrom]
            npa = np.append
            nsa = ns.append

            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                if score > cutoff:
                    npa(pos)
                    nsa(score)
        return ret

    def filter_score_below (self, cutoff=0):
        """Keep points below a score cutoff. Return a new WigTrackI.
        
        """
        ret = WigTrackI()
        ret.span = self.span
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(FBYTE4,[])]
            (np,ns) = ret.__data[chrom]
            npa = np.append
            nsa = ns.append

            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                if score < cutoff:
                    npa(pos)
                    nsa(score)
        return ret

    def write_gff (self, fhd, shift=0, source=".", feature="."):
        """Write all data to fhd in GFF format.

        shift will be used to shift the coordinates. default: 0
        """
        assert isinstance(fhd,file)
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                fhd.write(
                    "\t".join( (chrom,source,feature,
                                str(pi-shift),str(pi-shift+self.span-1),
                                str(si),'+','.'
                                ) )+"\n"
                    )

    def write_bed (self, fhd):
        """Write all data to fhd in BED format.
        
        """
        pass

    def remove_redundant (self):
        """Remove redundant position, keep the highest score.
        
        """
        chrs = set(self.__data.keys())
        ndata = {}
        for chrom in chrs:
            ndata[chrom] = [array(BYTE4,[]),array(FBYTE4,[])] # for (pos,value)
            nd_p_append = ndata[chrom][0].append
            nd_s_append = ndata[chrom][1].append
            (p,s) = self.__data[chrom]
            prev_p = None
            prev_s = None
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if not prev_p:
                    prev_p = pi
                    prev_s = si
                else:
                    if pi == prev_p:
                        if si>prev_s:
                            prev_s = si
                    else:
                       nd_p_append (prev_p)
                       nd_s_append (prev_s)
            nd_p_append (prev_p)
            nd_s_append (prev_s)
        del self.__data
        self.__data = ndata

    def find_peaks (self, bw=None):
        """A naive peak finding algorithm to find all local maximum
        points.

        bw : Bandwidth. Two local maximums will compete if their
        distance is less than bw. The higher one will be kept, or if
        they are equal, the last point will be kept. Default is 4*span.

        Return type:
        
        find_peak will return a new WigTrackI object with position
        and score for each local maximum(peaks).
        """
        if not bw:
            bw = 4*self.span
        bw = max(1,bw)
        ret = WigTrackI()
        ret.span = 1
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            prev_peak_p = -10000
            prev_peak_s = -10000
            prev_peak_summits = []
            prev_p = -10000
            prev_s = -10000
            increase = False
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if si > prev_s:
                    # increase
                    increase = True
                elif si < prev_s:
                    # decrease
                    if increase:
                        # prev_p is a summit
                        #print prev_p,prev_s,pi,si
                        if prev_p-prev_peak_p < bw:
                            # two summits are near
                            if prev_s > prev_peak_s:
                                # new summit is high
                                prev_peak_s = prev_s
                                prev_peak_p = prev_p
                        else:
                            # it's a new summit
                            #print chrom,prev_peak_p,prev_peak_s
                            try:
                                ret.add_loc(chrom,prev_peak_p,prev_peak_s)
                            except OverflowError:
                                pass
                            prev_peak_s = prev_s
                            prev_peak_p = prev_p
                            
                    increase = False
                prev_p = pi
                prev_s = si
            try:
                ret.add_loc(chrom,prev_peak_p,prev_peak_s)            
            except OverflowError:
                pass
        return ret

    def find_valleys (self, bw=None):
        """A naive peak finding algorithm to find all local minimum
        points.

        bw : Bandwidth. Two local maximums will compete if their
        distance is less than bw. The higher one will be kept, or if
        they are equal, the last point will be kept. Default is 4*span.

        Return type:
        
        find_peak will return a new WigTrackI object with position
        and score for each local maximum(peaks).
        """
        if not bw:
            bw = 4*self.span
        bw = max(1,bw)
        ret = WigTrackI()
        ret.span = 1
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            prev_peak_p = -10000
            prev_peak_s = -10000
            prev_peak_summits = []
            prev_p = -10000
            prev_s = -10000
            decrease = False
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if si < prev_s:
                    # decrease
                    decrease = True
                elif si > prev_s:
                    # increase
                    if decrease:
                        # prev_p is a valley
                        #print prev_p,prev_s,pi,si
                        if prev_p-prev_peak_p < bw:
                            # two summits are near
                            if prev_s < prev_peak_s:
                                # new summit is lower
                                prev_peak_s = prev_s
                                prev_peak_p = prev_p
                        else:
                            # it's a new valley
                            #print chrom,prev_peak_p,prev_peak_s
                            try:
                                ret.add_loc(chrom,prev_peak_p,prev_peak_s)
                            except OverflowError:
                                pass
                            prev_peak_s = prev_s
                            prev_peak_p = prev_p
                            
                    decrease = False
                prev_p = pi
                prev_s = si
            try:
                ret.add_loc(chrom,prev_peak_p,prev_peak_s)            
            except OverflowError:
                pass
        return ret

    def summary (self):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).
        
        """
        n_v = 0
        sum_v = 0
        max_v = -100000
        min_v = 100000
        for (p,v) in self.__data.values():
            sum_v += sum(v)
            n_v += len(v)
            max_v = max(max(v),max_v)
            min_v = min(min(v),min_v)
        mean_v = float(sum_v)/n_v
        variance = 0.0
        for (p,v) in self.__data.values():
            for vv in v:
                tmp = vv-mean_v
                variance += tmp*tmp
        variance /= float(n_v-1)
        std_v = sqrt(variance)
        return (sum_v, max_v, min_v, mean_v, std_v)

    def null_model_summary (self, sample=10):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).

        This is for NULL model which is a symetric normal distribution
        abased on sample of the whole data set.
        """
        # sample data step
        data_step = int(100/sample)

        null_list = array(FBYTE4,[])
        na = null_list.append
        for (p,v) in self.__data.values():
            i = 0
            for vv in v:
                i+=1
                if i==data_step:
                    na(vv)
                    i=0
                
        
        sum_v = sum(null_list)
        mean_v = sum_v/float(len(null_list))

        null_list_len = len(null_list)
        null_list = sorted(null_list)
        median_index1 = (null_list_len - 1) / 2
        median_index2 = null_list_len / 2
        median_v = (null_list[median_index1]+null_list[median_index2])/2.0

        # make right part of nullList

        for i in xrange(null_list_len/2):
            null_list[null_list_len-i-1] = 2* median_v - null_list[i]
        
        std_v = std(null_list)

        return (sum_v,max(null_list),min(null_list),median_v,std_v)

    def normalize (self,null=False,sample_percent=10):
        """Normalize values centered at 0 and variance as 1.

        If null is True, it will use the null list to calculate mean and std.
        When null is True, sample_percent will be passed to null_model to sample the data.
        """
        if null:
            (sum_v,max_v,min_v,mean_v,std_v) = self.null_model_summary(sample=sample_percent)
        else:
            (sum_v,max_v,min_v,mean_v,std_v) = self.summary()
        for (p,v) in self.__data.values():
            for i in range(len(v)):
                v[i] = float(v[i]-mean_v)/std_v
        return (sum_v, max_v, min_v, mean_v, std_v)
                

    def call_peaks (self, cutoff=1, up_limit=1e310, min_length=200, max_gap=50):
        """This function try to find some region within which, scores
        are continuously higher than a cutoff.

        cutoff:  cutoff of value, default 1
        min_length :  minimum peak length, default 200
        gap   :  maximum gap to merge nearby peaks
        """
        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            (ps,ss) = self.get_data_by_chr(chrom)
            psn = iter(ps).next         # assign the next function to a virable to speed up
            ssn = iter(ss).next
            x = 0
            while True:
                # find the first point above cutoff
                try:
                    p = psn()
                    s = ssn()
                except:
                    break
                x += 1                  # index for the next point
                if s >= cutoff and s<=up_limit:
                    peak_content = [(p,s),]
                    break               # found the first point above cutoff

            for i in range(x,len(ps)):
                p = psn()
                s = ssn()
                if s < cutoff or s > up_limit:
                    continue
                # for points above cutoff
                if p - peak_content[-1][0] <= max_gap:
                    peak_content.append((p,s))
                else:
                    # a new peak
                    peak_length = peak_content[-1][0]-peak_content[0][0]+self.span
                    if peak_length >= min_length:
                        summit = None
                        summit_score = None
                        for (m,n) in peak_content:
                            if not summit_score or summit_score < n:
                                summit = m
                                summit_score = n
                        peaks.add(chrom,peak_content[0][0],peak_content[-1][0]+self.span,
                                  summit=summit-peak_content[0][0],peak_height=summit_score)
                        #print chrom,peak_content[0][0],peak_content[-1][0]+self.span,peak_length
                    peak_content = [(p,s),]
        return peaks

    def total (self):
        t = 0
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            t += len(p)
        return t


class bedGraphTrackI:
    """Class for bedGraph type data.

    In bedGraph, data are represented as non-overlapping regions in
    the whole genome. I keep this assumption in all the functions. If
    data has overlaps, some functions will definitely give incorrect result.
    
    """
    def __init__ (self):
        self.__data = {}
        self.maxvalue =-10000
        self.minvalue = 10000

    def add_loc (self,chromosome,startpos,endpos,value):
        if not self.__data.has_key(chromosome):
            self.__data[chromosome] = [array(BYTE4,[]),array(BYTE4,[]),array(FBYTE4,[])] # for (startpos,endpos,value)
        self.__data[chromosome][0].append(startpos)
        self.__data[chromosome][1].append(endpos)
        self.__data[chromosome][2].append(value)        
        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value
        
    def sort (self):
        """Naive sorting for tags. Sort the data by the start
        position.

        """
        for k in self.__data.keys():
            (ps,pe,v) = self.__data[k]
            pv = zip(ps,pe,v)
            pv = sorted(pv)
            self.__data[k] = [array(BYTE4,[]),array(BYTE4,[]),array(FBYTE4,[])]
            psappend = self.__data[k][0].append
            peappend = self.__data[k][1].append            
            vappend = self.__data[k][2].append
            for (tps,tpe,tv) in pv:
                psappend(tps)
                peappend(tpe)                
                vappend(tv)

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([start pos],[end pos],[value])
        """
        if self.__data.has_key(chromosome):
            return self.__data[chromosome]
        else:
            return None

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.__data.keys())
        return l

    def write_bedGraph (self, fhd, name):
        """Write all data to fhd in Wiggle Format.

        shift will be used to shift the coordinates. default: 0
        """
        fhd.write("track type=bedGraph name=\"%s\" description=\"%s\"\n" % (name,name))
        chrs = self.get_chr_names()
        for chrom in chrs:
            (ps,pe,v) = self.__data[chrom]
            psnext = iter(ps).next
            penext = iter(pe).next            
            vnext = iter(v).next            
            for i in xrange(len(ps)):
                startpos = psnext()
                endpos = penext()
                score = snext()
                fhd.write("%d\t%d\t%.4f\n" % (startpos,endpos,score))                

    def filter_score (self, cutoff=0):
        """Filter using a score cutoff. Return a new bedGraph.
        
        """
        ret = bedGraphTrackI()
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (ps,pe,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(BYTE4,[]),array(FBYTE4,[])]
            (nps,npe,ns) = ret.__data[chrom] # new pos start, new pos end, new score
            npsa = nps.append
            npea = npe.append
            nsa = ns.append

            psnext = iter(ps).next
            penext = iter(pe).next            
            snext = iter(s).next            
            for i in xrange(len(ps)):
                startpos = psnext()
                endpos = penext()                
                score = snext()
                if score > cutoff:
                    npsa(startpos)
                    npea(endpos)
                    nsa(score)
        return ret

    def filter_score_below (self, cutoff=0):
        """Keep points below a score cutoff. Return a new bedGraphTrackI.
        
        """
        ret = bedGraphTrackI()
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (ps,pe,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(BYTE4,[]),array(FBYTE4,[])]
            (nps,npe,ns) = ret.__data[chrom] # new pos start, new pos end, new score
            npsa = nps.append
            npea = npe.append
            nsa = ns.append

            psnext = iter(ps).next
            penext = iter(pe).next            
            snext = iter(s).next            
            for i in xrange(len(ps)):
                startpos = psnext()
                endpos = penext()                
                score = snext()
                if score < cutoff:
                    npsa(startpos)
                    npea(endpos)
                    nsa(score)
        return ret

    def summary (self):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).
        
        """
        n_v = 0
        sum_v = 0
        max_v = -100000
        min_v = 100000
        for (s,e,v) in self.__data.values():
            # for each chromosome
            for i in range(len(s)):
                # for each region
                l = e[i]-s[i]
                sum_v += v[i]*l
                n_v += l
            max_v = max(max(v),max_v)
            min_v = min(min(v),min_v)
        mean_v = float(sum_v)/n_v
        variance = 0.0
        for (s,e,v) in self.__data.values():
            for vv in v:
                tmp = vv-mean_v
                variance += tmp*tmp
        variance /= float(n_v-1)
        std_v = sqrt(variance)
        return (sum_v, max_v, min_v, mean_v, std_v)

    def null_model_summary (self, sample_step=100):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).

        This is for NULL model which is a symetric normal distribution
        based on sample% of the whole data set.
        """
        # sample data step
        data_step = sample_step

        null_list = array(FBYTE4,[])
        na = null_list.append
        for (s,e,v) in self.__data.values():
            i = 0
            for vv in v:
                i+=1
                if i==data_step:
                    na(vv)
                    i=0
        
        sum_v = sum(null_list)
        mean_v = sum_v/float(len(null_list))

        null_list_len = len(null_list)
        null_list = sorted(null_list)
        median_index1 = (null_list_len - 1) / 2
        median_index2 = null_list_len / 2
        median_v = (null_list[median_index1]+null_list[median_index2])/2.0

        # make right part of nullList

        for i in xrange(null_list_len/2):
            null_list[null_list_len-i-1] = 2* median_v - null_list[i]
        
        std_v = std(null_list)

        return (sum_v,max(null_list),min(null_list),median_v,std_v)

    def normalize (self,null=False,sample_step=10):
        """Normalize values centered at 0 and variance as 1.

        If null is True, it will use the null list to calculate mean and std.
        When null is True, sample_percent will be passed to null_model to sample the data.
        """
        if null:
            (sum_v,max_v,min_v,mean_v,std_v) = self.null_model_summary(sample_step=sample_step)
        else:
            (sum_v,max_v,min_v,mean_v,std_v) = self.summary()
        for (ps,pe,v) in self.__data.values():
            for i in range(len(v)):
                v[i] = float(v[i]-mean_v)/std_v
        return (sum_v, max_v, min_v, mean_v, std_v)
                

    def call_peaks (self, cutoff=1, up_limit=1e310, min_length=200, max_gap=50):
        """This function try to find some region within which, scores
        are continuously higher than a cutoff.

        cutoff:  cutoff of value, default 1
        min_length :  minimum peak length, default 200
        gap   :  maximum gap to merge nearby peaks
        """
        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            (pss,pse,ss) = self.get_data_by_chr(chrom)
            pssn = iter(pss).next         # assign the next function to a virable to speed up
            pesn = iter(pes).next            
            ssn = iter(ss).next
            x = 0
            while True:
                # find the first region above cutoff
                try:                    # try to read the first data range for this chrom
                    ps = pssn()
                    pe = pesn()                    
                    s = ssn()
                except:
                    break
                x += 1                  # index for the next point
                if s >= cutoff and s<=up_limit:
                    peak_content = [(ps,pe,s),]
                    break               # found the first range above cutoff

            for i in range(x,len(ps)):
                ps = pssn()
                pe = pesn()
                s = ssn()                
                if s < cutoff or s > up_limit:
                    continue
                # for points above cutoff
                # if the gap is allowed
                if ps - peak_content[-1][1] <= max_gap:
                    peak_content.append((ps,pe,s))
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[-1][1]-peak_content[0][0]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        summit = None
                        summit_score = None
                        for (tstart,tend,tscore) in peak_content:
                            if not summit_score or summit_score < tscore:
                                summit = int((tend+tstart)/2)
                                summit_score = tscore
                        peaks.add(chrom,peak_content[0][0],peak_content[-1][1],
                                  summit=summit-peak_content[0][0],score=summit_score)
                    # start a new peak
                    peak_content = [(ps,pe,s),]
            # save the last peak
            if peak_length >= min_length: # if the peak is too small, reject it
                summit = None
                summit_score = None
                for (tstart,tend,tscore) in peak_content:
                    if not summit_score or summit_score < tscore:
                        summit = int((tend+tstart)/2)
                        summit_score = tscore
                    peaks.add(chrom,peak_content[0][0],peak_content[-1][1],
                              summit=summit-peak_content[0][0],score=summit_score)
            
        return peaks

    def total (self):
        t = 0
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (ps,pe,s) = self.__data[chrom]
            t += len(ps)
        return t
