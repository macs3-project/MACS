# Time-stamp: <2011-05-20 00:18:01 Tao Liu>

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

from array import array
from random import sample as random_sample
from operator import itemgetter

from libc.math cimport sqrt

from MACS2.Constants import *

# ------------------------------------
# constants
# ------------------------------------
__version__ = "FeatIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "PeakIO, FWTrackII, TrackI, and WigTrackI classes"

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
    
    def add (self, char * chromosome, long start, long end, long summit = 0, 
             double peak_height=0, int number_tags=0, 
             double pvalue=0, double fold_enrichment=0, double fdr=0):
        """items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append((start,end,end-start,summit,
                                       peak_height,number_tags,
                                       pvalue,fold_enrichment,fdr))

    def filter_pvalue (self, double pvalue_cut ):
        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
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
        """Print out peaks in BED5 format.

        Five columns are chromosome, peak start, peak end, peak name, and peak height.

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                text+= "%s\t%d\t%d\tpeak_%d\t%.2f\n" % (chrom,peak[0],peak[1],n_peak,peak[4])
        return text

    def to_summits_bed (self):
        """Print out peak summits in BED5 format.

        Five columns are chromosome, summit start, summit end, peak name, and peak height.

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                summit_p = peak[0]+peak[3]
                text+= "%s\t%d\t%d\tpeak_%d\t%.2f\n" % (chrom,summit_p,summit_p+1,n_peak,peak[4])
        return text

    def write_to_bed (self, fhd, name_prefix="peak_", score_column=4):
        """Write peaks in BED5 format in a file handler. Score (5th
        column) is decided by score_column setting. Check the
        following list. Name column ( 4th column) is made by putting
        name_prefix together with an ascending number.

        Five columns are chromosome, peak start, peak end, peak name, and peak score.

        items in peak object:
        0. peak start
        1. peak end
        2. peak length
        3. peak summit
        4. peak height
        5. number of tags in peak region
        6. peak pvalue
        7. peak fold_enrichment
        8. fdr
        """
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                fhd.write( "%s\t%d\t%d\t%s%d\t%.2f\n" % (chrom,peak[0],peak[1],name_prefix,n_peak,peak[score_column]) )


    def write_to_summit_bed (self, fhd, name_prefix="peak_", score_column=4):
        """Write peak summits in BED5 format in a file handler. Score
        (5th column) is decided by score_column setting. Check the
        following list. Name column ( 4th column) is made by putting
        name_prefix together with an ascending number.

        Five columns are chromosome, summit start, summit end, peak name, and peak score.

        items in peak object:
        0. peak start
        1. peak end
        2. peak length
        3. peak summit
        4. peak height
        5. number of tags in peak region
        6. peak pvalue
        7. peak fold_enrichment
        8. fdr
        """
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                summit_p = peak[0]+peak[3]
                fhd.write( "%s\t%d\t%d\t%s%d\t%.2f\n" % (chrom,summit_p,summit_p+1,name_prefix,n_peak,peak[score_column]) )

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
        #del peaks
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
    def __init__ (self, baseline_value=0):
        """
        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with baseline_value.
        
        """
        self.__data = {}
        self.maxvalue =-10000
        self.minvalue = 10000
        self.baseline_value = baseline_value

    def add_loc (self,chromosome,startpos,endpos,value):
        """Add a chr-start-end-value block into __data dictionary.

        """
        # basic assumption, end pos should > start pos
        assert endpos > startpos, "endpos %d can't be smaller than start pos %d" % (endpos,startpos)
        
        if not self.__data.has_key(chromosome):
            self.__data[chromosome] = [array(BYTE4,[]),array(FBYTE4,[])] # for (endpos,value)
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
            pre_pos = c[0][-1]
            pre_v   = c[1][-1]
            # to check 1. continuity; 2. non-overlapping
            assert pre_pos < endpos , "bedGraph regions are not continuous."
            assert pre_pos <= startpos , "bedGraph regions have overlappings."
            
            if startpos != pre_pos:
                # there is a gap, so fill it with baseline_value
                c[0].append(startpos)
                c[1].append(self.baseline_value)
                # then add this region
                c[0].append(endpos)
                c[1].append(value)
            else:
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

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
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

    def write_bedGraph (self, fhd, name, description):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        shift will be used to shift the coordinates. default: 0
        """
        #fhd.write("track type=bedGraph name=\"%s\" description=\"%s\"\n" % (name,description))
        chrs = self.get_chr_names()
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).next
            vnext = iter(v).next
            pre = 0
            for i in xrange(len(p)):
                pos = pnext()
                value = vnext()
                #if value != self.baseline_value:
                # never write baseline_value
                fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,pos,value))
                pre = pos

    def reset_baseline (self, baseline_value):
        """Reset baseline value to baseline_value.

        So any region between self.baseline_value and baseline_value
        will be set to baseline_value.
        
        """
        self.baseline_value = baseline_value
        self.filter_score(cutoff=baseline_value)
        self.merge_regions()

    def merge_regions (self):
        """Merge nearby regions with the same value.
        
        """
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).next
            vnext = iter(v).next

            # new arrays
            new_pos = array(BYTE4,[pnext(),])
            new_value = array(FBYTE4,[vnext(),])

            newpa = new_pos.append
            newva = new_value.append
            
            new_pre_pos = new_pos[0]
            new_pre_value = new_value[0]

            for i in xrange(1,len(p)):
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
                
    def filter_score (self, cutoff=0):
        """Filter using a score cutoff. Any region lower than score
        cutoff will be set to self.baseline_value.

        Self will be modified.
        """
        ret = bedGraphTrackI()
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,v) = self.__data[chrom]
            pnext = iter(p).next
            vnext = iter(v).next

            # new arrays
            new_pos = array(BYTE4,[])
            new_value = array(FBYTE4,[])
            new_pre_pos = 0
            new_pre_value = None

            for i in xrange(len(p)):
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

    def summary (self):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).
        
        """
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
        mean_v = float(sum_v)/n_v
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
        return (sum_v, max_v, min_v, mean_v, std_v)

    def call_peaks (self, cutoff=1, up_limit=1e310, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 1.
        up_limit: the highest acceptable value. Default 10^{310}
          * so only allow peak with value >=cutoff and <=up_limit
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        """
        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            (ps,vs) = self.get_data_by_chr(chrom) # arrays for position and values
            psn = iter(ps).next         # assign the next function to a viable to speed up
            vsn = iter(vs).next
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
                if v >= cutoff and v <= up_limit:
                    peak_content = [(pre_p,p,v),]
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in range(x,len(ps)):
                # continue scan the rest regions
                p = psn()
                v = vsn()
                if v < cutoff or v > up_limit: # not be detected as 'peak'
                    pre_p = p
                    continue
                # for points above cutoff
                # if the gap is allowed
                if pre_p - peak_content[-1][1] <= max_gap:
                    peak_content.append((pre_p,p,v))
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[-1][1]-peak_content[0][0]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        tsummit = []
                        summit = None
                        summit_value = None
                        for (tstart,tend,tvalue) in peak_content:
                            if not summit_value or summit_value < tvalue:
                                tsummit = [int((tend+tstart)/2),]
                                summit_value = tvalue
                            elif summit_value == tvalue:
                                tsummit.append( int((tend+tstart)/2) )
                        summit = tsummit[int((len(tsummit)+1)/2)-1 ]
                        peaks.add(chrom,peak_content[0][0],peak_content[-1][1],
                                  summit=summit-peak_content[0][0],peak_height=summit_value) # summit is the relative position to peak start
                    # start a new peak
                    peak_content = [(pre_p,p,v),]
                pre_p = p
                
            # save the last peak
            if peak_length >= min_length: # if the peak is too small, reject it
                summit = None
                summit_value = None
                for (tstart,tend,tvalue) in peak_content:
                    if not summit_value or summit_value < tvalue:
                        summit = int((tend+tstart)/2)
                        summit_value = tvalue
                peaks.add(chrom,peak_content[0][0],peak_content[-1][1],
                          summit=summit-peak_content[0][0],peak_height=summit_value)
            
        return peaks

    def total (self):
        """Return the number of regions in this object.

        """
        t = 0
        for (p,s) in self.__data.values():
            t += len(p)
        return t


    def overlie (self, bdgTrack2, func=max ):
        """Calculate two bedGraphTrackI objects by letting self
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


        Return value is a bedGraphTrackI object.
        """
        assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = bedGraphTrackI()
        retadd = ret.add_loc
        
        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for chrom in common_chr:
            (p1s,v1s) = self.get_data_by_chr(chrom) # arrays for position and values
            p1n = iter(p1s).next         # assign the next function to a viable to speed up
            v1n = iter(v1s).next

            (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
            p2n = iter(p2s).next         # assign the next function to a viable to speed up
            v2n = iter(v2s).next

            pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
            try:
                p1 = p1n()
                v1 = v1n()

                p2 = p2n()
                v2 = v2n()

                while True:
                    if p1 < p2:
                        # clip a region from pre_p to p1, then set pre_p as p1.
                        retadd(chrom,pre_p,p1,func(v1,v2))
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        retadd(chrom,pre_p,p2,func(v1,v2))
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        retadd(chrom,pre_p,p1,func(v1,v2))
                        pre_p = p1
                        # call for the next p1, v1, p2, v2.
                        p1 = p1n()
                        v1 = v1n()
                        p2 = p2n()
                        v2 = v2n()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass
        
        ret.merge_regions()
        return ret
                       
    def apply_func ( self, func ):
        """Apply function 'func' to every value in this bedGraphTrackI object.

        *Two adjacent regions with same value after applying func will
        not be merged.
        """
        t = 0
        for (p,s) in self.__data.values():
            for i in xrange(len(s)):
                s[i] = func(s[i])
        self.maxvalue = func(self.maxvalue)
        self.minvalue = func(self.minvalue)
        return True
        
