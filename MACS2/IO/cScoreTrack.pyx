# cython: profile=True
# Time-stamp: <2012-04-13 17:41:47 Tao Liu>

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
import numpy as np
cimport numpy as np
from np import int64,int32,float32

from libc.math cimport log10,log

from MACS2.Constants import *
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO, BroadPeakIO

from MACS2.hashtable import Int64HashTable

import logging

#from time import time as ttime

#from MACS2.IO.cBedGraph import bedGraphTrackI

# ------------------------------------
# constants
# ------------------------------------
__version__ = "scoreTrackI $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "scoreTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b


pscore_dict = {}
LOG10_E = 0.43429448190325176


cdef get_pscore ( int observed, double expectation ):
    """Get p-value score from Poisson test. First check existing
    table, if failed, call poisson_cdf function, then store the result
    in table.
    
    """
    cdef int score
    
    key_value = ( observed, expectation )
    if pscore_dict.has_key(key_value):
        return pscore_dict[key_value]
    else:
        score = int(-100*poisson_cdf(observed,expectation,False,True))
        pscore_dict[(observed,expectation)] = score
        return score

logLR_dict = {}

cdef logLR ( double x, double y ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Then store the values in integar form =
    100*logLR. Set minus sign for depletion.
    
    """
    cdef int s
    #key_value = ( x, y )
    #if logLR_dict.has_key(key_value):
    #    return logLR_dict[key_value]
    #else:
    if x > y:
        s = int( (x*(log(x+1)-log(y+1))+y-x)*LOG10_E*100 )
    elif x < y:
        s = int( (-1*x*(log(x+1)-log(y+1))-y+x)*LOG10_E*100 )
    else:
        s = 0
    #logLR_dict[key_value] = s
    return s

# ------------------------------------
# Classes
# ------------------------------------

class scoreTrackI:
    """Class for scoreGraph type data. Modified from
    bedGraphTrackI. The only difference is that we store pvalue score,
    qvalue score and foldchange together.

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
    def __init__ (self):
        """Different with bedGraphTrackI, missing values are simply
        replaced with 0.
        
        """
        self.data = {}
        self.pointer = {}

    def add_chromosome ( self, str chrom, int chrom_max_len ):
        if not self.data.has_key(chrom):
            self.data[chrom] = np.zeros(chrom_max_len,dtype=[('pos','int32'),
                                                             ('sample','float32'),
                                                             ('control','float32'),
                                                             ('-100logp','int32'),
                                                             ('-100logq','int32'),
                                                             ('100logLR','int32'),])
            self.pointer[chrom] = 0

    def add (self, str chromosome, int endpos, int sample, double control):
        """Add a chr-endpos-sample-control block into data
        dictionary. At the mean time, calculate pvalues.

        """
        cdef int i
        c = self.data[chromosome]
        i = self.pointer[chromosome]
        # get the preceding region
        #c[i] = (endpos,sample,control,int(-100*poisson_cdf(sample,control,False,True)),0)
        a = c[i]
        a['pos'] = endpos
        a['sample'] = sample
        a['control'] = control
        a['-100logp'] = get_pscore(sample,control)
        a['100logLR'] = logLR(sample,control)
        #c[i] = (endpos,sample,control,get_pscore(sample,control),0,logLR(sample,control))
        self.pointer[chromosome] += 1

    def finalize (self):
        cdef str chrom
        cdef int l

        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.pointer[chrom]
            d.resize(l,refcheck=False)

    def get_data_by_chr (self, str chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if self.data.has_key(chromosome):
            return self.data[chromosome]
        else:
            return None

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.data.keys())
        return l

    def write_bedGraph (self, fhd, str name, str description, str colname):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        colname: can be 'sample','control','-100logp','-100logq'

        """
        cdef str chrom
        cdef int l, pre, i, p 
        cdef double pre_v, v
        
        if colname not in ['sample','control','-100logp','-100logq','100logLR']:
            raise Exception("%s not supported!" % colname)
        if colname in ['-100logp', '-100logq','100logLR']:
            flag100 = True              # for pvalue or qvalue, divide them by 100 while writing to bedGraph file
        else:
            flag100 = False
        chrs = self.get_chr_names()
        for chrom in chrs:
            d = self.data[chrom]
            l = self.pointer[chrom]
            pre = 0
            pos   = d['pos']
            if flag100:
                value = d[colname]/100.0
            else:
                value = d[colname]
            pre_v = value[0]
            for i in range( 1, l ):
                v = value[i]
                p = pos[i-1]
                if pre_v != v: 
                    fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,p,pre_v))
                    pre_v = v
                    pre = p
            p = pos[-1]
            # last one
            fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,p,pre_v))            

        return True

    def make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-100log10pvalue:(-100log10qvalue,rank,basepairs)} relationships.
        """
        cdef long n, pre_p, this_p, length, j, pre_l, l, this_v, pre_v, v
        cdef long N, k, q, pre_q
        cdef double f
        cdef str chrom
        
        #logging.info("####test#### start make_pq")
        n = self.total()
        #value_list = np.empty( n, dtype = [('v', '<f4'), ('l', '<i4')])
        value_dict = {}
        #i = 0                           # index for value_list
        for chrom in self.data.keys():
            # for each chromosome
            pre_p  = 0
            pos    = iter(self.data[chrom][ 'pos' ]).next
            value  = iter(self.data[chrom][ '-100logp' ]).next
            length = self.pointer[chrom]
            j = 0
            while j<length:
                this_p = pos()
                this_v = value()
                #value_list[i] = (this_v,this_p-pre_p)
                #i += 1
                if value_dict.has_key(this_v):
                    value_dict[this_v] += this_p-pre_p
                else:
                    value_dict[this_v] = this_p-pre_p
                j += 1
                pre_p = this_p

        N = sum(value_dict.values())
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        pvalue2qvalue = {pre_v:[0,k,0]}              # pvalue:[qvalue,rank,bp_with_this_pvalue]
        #logging.info("####test#### start matching pvalue to qvalue")
        for v in sorted(value_dict.keys(),reverse=True):
            l = value_dict[v]
            q = v+int((log10(k)+f)*100) # we save integars here.
            q = max(0,min(pre_q,q))           # make q-score monotonic
            pvalue2qvalue[v] = [q, k, 0]
            pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
            pre_v = v
            pre_q = q
            k+=l
        pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
        #logging.info("####test#### finish building pqtable")        
        # pop the first -1e100 one
        pvalue2qvalue.pop(-2147483647)

        return pvalue2qvalue

    def assign_qvalue ( self , dict pvalue2qvalue ):
        """Assign -100log10qvalue to every point.

        pvalue2qvalue: a dictionary of -100log10pvalue:-100log10qvalue
        """
        cdef long i,l,j,p
        cdef str chrom
        chroms = self.data.keys()

        # convert pvalue2qvalue to a simple dict
        s_p2q = Int64HashTable()
        g = pvalue2qvalue.get
        for i in pvalue2qvalue.keys():
            s_p2q.set_item(i,g(i)[0])

        g = s_p2q.get_item
        
        for j in range( len(chroms) ):
            chrom = chroms[j]
            pvalue = self.data[chrom]['-100logp']
            qvalue = self.data[chrom]['-100logq']
            l = self.pointer[chrom]
            for i in range( l ):
                qvalue[i] = g(pvalue[i])
        return True

    def call_peaks (self, int cutoff=500, int min_length=200, int max_gap=50, str colname='-100logp'):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 1.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        colname: can be 'sample','control','-100logp','-100logq'. Cutoff will be applied to the specified column.
        ptrack:  an optional track for pileup heights. If it's not None, use it to find summits. Otherwise, use self/scoreTrack.
        """
        cdef int i
        cdef str chrom
        
        assert (colname in [ 'sample', 'control', '-100logp', '-100logq', '100logLR' ]), "%s not supported!" % colname

        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        for chrom in chrs:
            peak_content = []           # to store points above cutoff

            above_cutoff = np.nonzero( self.data[chrom][colname] >= cutoff )[0] # indices where score is above cutoff
            above_cutoff_v = self.data[chrom][colname][above_cutoff] # scores where score is above cutoff

            above_cutoff_endpos = self.data[chrom]['pos'][above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = self.data[chrom]['pos'][above_cutoff-1] # start positions of regions where score is above cutoff
            above_cutoff_sv= self.data[chrom]['sample'][above_cutoff] # sample pileup height where score is above cutoff

            if above_cutoff_v.size == 0:
                continue

            if above_cutoff[0] == 0:
                # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
                above_cutoff_startpos[0] = 0

            # first bit of region above cutoff
            peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], above_cutoff_v[0], above_cutoff_sv[0], above_cutoff[0]) )
            for i in range( 1,above_cutoff_startpos.size ):
                if above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                    # append
                    peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]) )
                else:
                    # close
                    self.__close_peak(peak_content, peaks, min_length, chrom, colname )
                    peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]
            
            # save the last peak
            if not peak_content:
                continue
            else:
                self.__close_peak(peak_content, peaks, min_length, chrom, colname )

        return peaks

    def __close_peak (self, peak_content, peaks, int min_length, str chrom, str colname):
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        """
        cdef int summit_pos, tmpstart, tmpend, tmpindex, middle_summit, summit_index, i
        cdef double summit_value, tmpvalue, tmpsummitvalue

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tmpsummit = []
            summit_pos   = 0
            summit_value = summit_pos
            for i in range(len(peak_content)):
                (tmpstart,tmpend,tmpvalue,tmpsummitvalue, tmpindex) = peak_content[i]
                #for (tmpstart,tmpend,tmpvalue,tmpsummitvalue, tmpindex) in peak_content:
                if not summit_value or summit_value < tmpsummitvalue:
                    tmpsummit = [ ( tmpend+tmpstart )/2, ]
                    tmpsummit_index = [ tmpindex, ]
                    summit_value = tmpsummitvalue
                elif summit_value == tmpsummitvalue:
                    # remember continuous summit values
                    tmpsummit.append( int( (tmpend+tmpstart)/2 ) )
                    tmpsummit_index.append( tmpindex )
            middle_summit = int( ( len(tmpsummit)+1 )/2 )-1 # the middle of all highest points in peak region is defined as summit
            summit_pos    = tmpsummit[ middle_summit ]
            summit_index  = tmpsummit_index[ middle_summit ]
            peaks.add( chrom,
                       peak_content[0][0],
                       peak_content[-1][1],
                       summit      = summit_pos,
                       peak_score  = self.data[chrom][colname][ summit_index ],
                       pileup      = self.data[chrom]['sample'][ summit_index ], # should be the same as summit_value
                       pscore      = self.data[chrom]['-100logp'][ summit_index ]/100.0,
                       fold_change = self.data[chrom]['sample'][ summit_index ]/self.data[chrom]['control'][ summit_index ],
                       qscore      = self.data[chrom]['-100logq'][ summit_index ]/100.0,
                       )
            # start a new peak
            return True

    def call_broadpeaks (self, int lvl1_cutoff=500, int lvl2_cutoff=100,
                         int min_length=200, int lvl1_max_gap=50, int lvl2_max_gap=400,
                         str colname='-100logq'):
        """This function try to find enriched regions within which,
        scores are continuously higher than a given cutoff for level
        1, and link them using the gap above level 2 cutoff with a
        maximum length of lvl2_max_gap.

        lvl1_cutoff:  cutoff of value at enriched regions, default 500.
        lvl2_cutoff:  cutoff of value at linkage regions, default 100.        
        min_length :  minimum peak length, default 200.
        lvl1_max_gap   :  maximum gap to merge nearby enriched peaks, default 50.
        lvl2_max_gap   :  maximum length of linkage regions, default 400.        
        colname: can be 'sample','control','-100logp','-100logq'. Cutoff will be applied to the specified column.

        Return both general PeakIO object for highly enriched regions
        and gapped broad regions in BroadPeakIO.
        """
        cdef int i
        cdef str chrom
        
        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."        
        lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff, min_length=min_length, max_gap=lvl1_max_gap, colname=colname)
        lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff, min_length=min_length, max_gap=lvl2_max_gap, colname=colname)
        chrs = lvl1_peaks.peaks.keys()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in chrs:
            lvl1peakschrom = lvl1_peaks.peaks[chrom]
            lvl2peakschrom = lvl2_peaks.peaks[chrom]
            lvl1peakschrom_next = iter(lvl1peakschrom).next
            tmppeakset = []             # to temporarily store lvl1 region inside a lvl2 region
            # our assumption is lvl1 regions should be included in lvl2 regions
            try:
                lvl1 = lvl1peakschrom_next()
            except StopIteration:
                break
            for i in range( len(lvl2peakschrom) ):
                # for each lvl2 peak, find all lvl1 peaks inside
                lvl2 = lvl2peakschrom[i]
                try:
                    while True:
                        if lvl2["start"] <= lvl1["start"]  and lvl1["end"] <= lvl2["end"]:
                            tmppeakset.append(lvl1)
                        else:
                            if tmppeakset:
                                self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
                            tmppeakset = []
                            break
                        lvl1 = lvl1peakschrom_next()
                except StopIteration:
                    if tmppeakset:
                        self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)                    
                    break
        
        return lvl1_peaks, broadpeaks

    # def call_broadpeaks2 (self, int lvl1_cutoff=500, int lvl2_cutoff=100, int min_length=200,
    #                       int lvl1_max_gap=50, int lvl2_max_gap=400, str colname='-100logq'):
    #     """This function try to find enriched regions within which,
    #     scores are continuously higher than a given cutoff for level
    #     1, and link them using the gap above level 2 cutoff with a
    #     maximum length of lvl2_max_gap.

    #     lvl1_cutoff:  cutoff of value at enriched regions, default 500.
    #     lvl2_cutoff:  cutoff of value at linkage regions, default 100.        
    #     min_length :  minimum peak length, default 200.
    #     lvl1_max_gap   :  maximum gap to merge nearby enriched peaks, default 50.
    #     lvl2_max_gap   :  maximum length of linkage regions, default 400.        
    #     colname: can be 'sample','control','-100logp','-100logq'. Cutoff will be applied to the specified column.

    #     Return both general PeakIO object for highly enriched regions
    #     and gapped broad regions in BroadPeakIO.
    #     """

    #     assert (colname in [ 'sample', 'control', '-100logp', '-100logq', '100logLR' ]), "%s not supported!" % colname

    #     chrs  = self.get_chr_names()
    #     bpeaks = BroadPeakIO()                      # dictionary to save broad peaks

    #     lvl2_cutoff = int(lvl2_cutoff)  # for broad regions
        
    #     for chrom in chrs:
    #         chrom_pointer = self.pointer[chrom]
    #         lvl2_peak_content = []           # to store points above cutoff

    #         above_cutoff = np.nonzero( self.data[chrom][colname] >= lvl2_cutoff )[0] # indices where score is above cutoff
    #         above_cutoff_v = self.data[chrom][colname][above_cutoff] # scores where score is above cutoff

    #         above_cutoff_endpos = self.data[chrom]['pos'][above_cutoff] # end positions of regions where score is above cutoff
    #         above_cutoff_startpos = self.data[chrom]['pos'][above_cutoff-1] # start positions of regions where score is above cutoff
    #         above_cutoff_sv= self.data[chrom]['sample'][above_cutoff] # sample pileup height where score is above cutoff

    #         if above_cutoff_v.size == 0:
    #             continue

    #         if above_cutoff[0] == 0:
    #             # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
    #             above_cutoff_startpos[0] = 0

    #         # first bit of region above cutoff
    #         lvl2_peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], above_cutoff_v[0], above_cutoff_sv[0], above_cutoff[0]) )
    #         for i in xrange(1,above_cutoff_startpos.size):
    #             if above_cutoff_startpos[i] - peak_content[-1][1] <= lvl2_max_gap:
    #                 # append
    #                 lvl2_peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]) )
    #             else:
    #                 # close
    #                 self.__close_broad_peak(lvl2_peak_content, bpeaks, lvl1_max_gap, lvl1_cutoff, colname, min_length, chrom )
    #                 lvl2_peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]
            
    #         # save the last peak
    #         if not peak_content:
    #             continue
    #         else:
    #             self.__close_broad_peak(lvl2_peak_content, bpeaks, lvl1_max_gap, lvl1_cutoff, colname, min_length, chrom )                
    #     return bpeaks

    # def __close_broad_peak( self, lvl2_peak_content, bpeaks, lvl1_max_gap, lvl1_cutoff, colname, min_length, chrom ):
    #     """ Finalize broad peak. Use 

    #     """
    #     pass


    def __add_broadpeak (self, bpeaks, str chrom, lvl2peak, lvl1peakset):
        """Internal function to create broad peak.
        
        """
        cdef int blockNum, thickStart, thickEnd, start, end
        cdef str blockSizes, blockStarts
        
        start      = lvl2peak["start"]
        end        = lvl2peak["end"]
        thickStart = lvl1peakset[0]["start"]
        thickEnd   = lvl1peakset[-1]["end"]
        blockNum   = int(len(lvl1peakset))
        blockSizes = ",".join( map(lambda x:str(x["length"]),lvl1peakset) )
        blockStarts = ",".join( map(lambda x:str(x["start"]-start),lvl1peakset) )
        if lvl2peak["start"] != thickStart:
            # add 1bp mark for the start of lvl2 peak
            blockNum += 1
            blockSizes = "1,"+blockSizes
            blockStarts = "0,"+blockStarts
        if lvl2peak["end"] != thickEnd:
            # add 1bp mark for the end of lvl2 peak            
            blockNum += 1
            blockSizes = blockSizes+",1"
            blockStarts = blockStarts+","+str(end-start-1)
        
        bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=thickStart, thickEnd=thickEnd,
                   blockNum = blockNum, blockSizes = blockSizes, blockStarts = blockStarts)
        return bpeaks

    def total ( self ):
        """Return the number of regions in this object.

        """
        cdef long t
        cdef str chrom
        
        t = 0
        for chrom in self.data.keys():
            t += self.pointer[chrom]
        return t


class CombinedTwoTrack:
    """ For differential peak calling.
    
    """
    def __init__ (self):
        """Different with bedGraphTrackI, missing values are simply
        replaced with 0.
        
        """
        self.data = {}
        self.pointer = {}

    def add_chromosome ( self, str chrom, int chrom_max_len ):
        if not self.data.has_key(chrom):
            self.data[chrom] = np.zeros(chrom_max_len,dtype=[('pos','int32'),
                                                             ('V1','float32'), # value for the first track
                                                             ('V2','float32'), # value for the second track
                                                             ])
            self.pointer[chrom] = 0

    def add ( self, str chromosome, int endpos, double V1, double V2 ):
        """Add a chr-endpos-sample-control block into data
        dictionary. At the mean time, calculate pvalues.

        """
        cdef list c
        cdef int i
        c = self.data[chromosome]
        i = self.pointer[chromosome]
        # get the preceding region
        c[i] = (endpos,V1,V2)
        self.pointer[chromosome] += 1

    def finalize ( self ):
        cdef str chrom
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.pointer[chrom]
            d.resize(l,refcheck=False)

    def get_data_by_chr (self, str chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if self.data.has_key(chromosome):
            return self.data[chromosome]
        else:
            return None

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        cdef set l
        l = set(self.data.keys())
        return l

    def write_bedGraph (self, fhd, str name, str description, str colname):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        colname: can be 'sample','control','-100logp','-100logq'

        """
        cdef str chrom
        cdef set chrs
        cdef int pre, i, l
        
        if colname not in ['V1','V2']:
            raise Exception("%s not supported!" % colname)
        chrs = self.get_chr_names()
        for chrom in chrs:
            d = self.data[chrom]
            l = self.pointer[chrom]
            pre = 0
            pos   = d['pos']
            value = d[colname]
            for i in range( l ):
                fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,pos[i],value[i]))
                pre = pos[i]

        return True

    def total ( self ):
        """Return the number of regions in this object.

        """
        cdef long t
        cdef str chrom
        t = 0
        for chrom in self.data.keys():
            t += self.pointer[chrom]
        return t

    def extract_value ( self, bdgTrack2 ):
        """It's like overlie function. THe overlapped regions between
        bdgTrack2 and self, will be recorded. The values from self in
        the overlapped regions will be outputed in a single array for
        follow statistics.

        """
        cdef set chr1, chr2, common_chr
        cdef str chrom
        cdef int pre_p, p1, p2
        cdef double v11, v21, v2
        
        #assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"

        ret = [[],array(FBYTE4,[]),array(FBYTE4,[]),array(BYTE4,[])] # region,V1,V1,length
        radd = ret[0].append
        v1add = ret[1].append
        v2add = ret[2].append        
        ladd = ret[3].append
        
        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        common_chr = chr1.intersection(chr2)
        for chrom in common_chr:
            chrom_data = self.get_data_by_chr(chrom) # arrays for position and values
            p1n = chrom_data['pos'].flat.next
            v11n = chrom_data['V1'].flat.next
            v21n = chrom_data['V2'].flat.next

            (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
            p2n = iter(p2s).next         # assign the next function to a viable to speed up
            v2n = iter(v2s).next

            pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
            try:
                p1 = p1n()
                v11 = v11n()
                v21 = v21n()                

                p2 = p2n()
                v2 = v2n()

                while True:
                    if p1 < p2:
                        # clip a region from pre_p to p1, then set pre_p as p1.
                        if v2>0:
                            radd(chrom+"."+str(pre_p)+"."+str(p1))
                            v1add(v11)
                            v2add(v21)                            
                            ladd(p1-pre_p)                        
                        pre_p = p1
                        # call for the next p1 and v1
                        p1 = p1n()
                        v11 = v11n()
                        v21 = v21n()
                    elif p2 < p1:
                        # clip a region from pre_p to p2, then set pre_p as p2.
                        if v2>0:
                            radd(chrom+"."+str(pre_p)+"."+str(p2))
                            v1add(v11)
                            v2add(v21)                            
                            ladd(p2-pre_p)                        
                        pre_p = p2
                        # call for the next p2 and v2
                        p2 = p2n()
                        v2 = v2n()
                    elif p1 == p2:
                        # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                        if v2>0:
                            radd(chrom+"."+str(pre_p)+"."+str(p1))
                            v1add(v11)
                            v2add(v21)                            
                            ladd(p1-pre_p)                        
                        pre_p = p1
                        # call for the next p1, v1, p2, v2.
                        p1 = p1n()
                        v11 = v11n()
                        v21 = v21n()
                        p2 = p2n()
                        v2 = v2n()
            except StopIteration:
                # meet the end of either bedGraphTrackI, simply exit
                pass

        # convert to np.array
        #ret = np.array([ret[0],ret[1],ret[2]]).transpose()
        #ret = ret[ret[0,0,:].argsort()]
        return ret


    def extract_average (self, bdgTrack2):
        cdef int i, l
        cdef str chrom, start, end
        
        (rarray,v1array,v2array,larray)  = self.extract_value(bdgTrack2)
        ret = [[],array(FBYTE4,[]),array(FBYTE4,[])] # region,V1,V1
        radd = ret[0].append
        v1add = ret[1].append
        v2add = ret[2].append
        cur_region = [None,None,None,None,None]      # chrom, start, end, s1, s2
        for i in range(len(rarray)):
            (chrom,start,end) = rarray[i].split('.')
            if chrom == cur_region[0] and start == cur_region[2]:
                cur_region[2] =  end
                cur_region[3] += v1array[i]*larray[i]
                cur_region[4] += v2array[i]*larray[i]
            else:
                if cur_region[0]:
                    l = int(cur_region[2])-int(cur_region[1])
                    radd(cur_region[0]+"."+str(cur_region[1])+"."+str(cur_region[2]))
                    v1add(cur_region[3]/float(l))
                    v2add(cur_region[4]/float(l))                    
                cur_region = [chrom, start, end, v1array[i]*larray[i], v2array[i]*larray[i]]

        radd(cur_region[0]+"."+str(cur_region[1])+"."+str(cur_region[2]))
        v1add(cur_region[3]/float(l))
        v2add(cur_region[4]/float(l))
        return ret
