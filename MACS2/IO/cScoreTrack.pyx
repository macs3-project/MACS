# Time-stamp: <2012-01-06 15:10:29 Tao Liu>

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
from numpy import int64,int32,float32

from libc.math cimport sqrt,log10

from MACS2.Constants import *
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO, BroadPeakIO
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

pscore_dict = {}

def get_pscore ( observed, expectation ):
    key_value = (observed, expectation)
    if pscore_dict.has_key(key_value):
        return pscore_dict[key_value]
    else:
        score = int(-100*poisson_cdf(observed,expectation,False,True))
        pscore_dict[(observed,expectation)] = score
        return score

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

    def add_chromosome ( self, chrom, chrom_max_len ):
        if not self.data.has_key(chrom):
            self.data[chrom] = np.zeros(chrom_max_len,dtype=[('pos','int32'),
                                                             ('sample','float32'),
                                                             ('control','float32'),
                                                             ('-100logp','int32'),
                                                             ('-100logq','int32')])
            self.pointer[chrom] = 0

    def add (self,chromosome,endpos,sample,control):
        """Add a chr-endpos-sample-control block into data
        dictionary. At the mean time, calculate pvalues.

        """
        c = self.data[chromosome]
        i = self.pointer[chromosome]
        # get the preceding region
        #c[i] = (endpos,sample,control,int(-100*poisson_cdf(sample,control,False,True)),0)
        c[i] = (endpos,sample,control,get_pscore(sample,control),0)
        self.pointer[chromosome] += 1

    def get_data_by_chr (self, chromosome):
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

    def write_bedGraph (self, fhd, name, description, colname):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        colname: can be 'sample','control','-100logp','-100logq'

        """
        if colname not in ['sample','control','-100logp','-100logq']:
            raise Exception("%s not supported!" % colname)
        if colname in ['-100logp', '-100logq']:
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
            for i in xrange( l ):
                fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,pos[i],value[i]))
                pre = pos[i]

        return True

    def __calculate_fold_change ( self, chrom, index ):
        """From 'sample' and 'control' columns, calculate foldchanges.

        chrom: chromosome name
        index: index in data[chrom]
        
        """
        return self.data[chrom]['sample'][index]/self.data[chrom]['control'][index]

    def make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-100log10pvalue:(-100log10qvalue,rank)} relationships.
        """
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
        pre_v = -1e100
        pre_l = 0
        pre_q = 1e100                       # save the previous q-value
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
        pvalue2qvalue.pop(-1e100)

        return pvalue2qvalue

    def assign_qvalue ( self , pvalue2qvalue ):
        """Assign -100log10qvalue to every point.

        pvalue2qvalue: a dictionary of -100log10pvalue:-100log10qvalue
        """
        t = 0
        for chrom in self.data.keys():
            pvalue = self.data[chrom]['-100logp']
            qvalue = self.data[chrom]['-100logq']
            for i in xrange( self.pointer[chrom] ):
                qvalue[i] = pvalue2qvalue[pvalue[i]][0]
        return True

    def call_peaks (self, cutoff=500, min_length=200, max_gap=50, colname='-100logp'):
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
        assert (colname in [ 'sample', 'control', '-100logp', '-100logq' ]), "%s not supported!" % colname

        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        #tloop = 0

        cutoff = int(cutoff)
        
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            peak_content = []           # to store points above cutoff

            #t0 = ttime()
            above_cutoff = np.nonzero( self.data[chrom][colname] >= cutoff )[0] # indices where score is above cutoff
            above_cutoff_flag = self.data[chrom][colname] >= cutoff
            above_cutoff_v = self.data[chrom][colname][above_cutoff] # scores where score is above cutoff

            above_cutoff_endpos = self.data[chrom]['pos'][above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = self.data[chrom]['pos'][above_cutoff_flag[1:]] # start positions of regions where score is above cutoff
            above_cutoff_sv= self.data[chrom]['sample'][above_cutoff] # sample pileup height where score is above cutoff

            if above_cutoff_v.size == 0:
                continue

            if above_cutoff[0] == 0:
                # first element > cutoff, insert first point in the chromosome
                np.insert(above_cutoff_startpos,0,0)

            # first bit of region above cutoff
            peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], above_cutoff_v[0], above_cutoff_sv[0], above_cutoff[0]) )
            for i in xrange(1,above_cutoff_startpos.size):
                if above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                    # append
                    peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]) )
                else:
                    # close
                    self.__close_peak(peak_content, peaks, min_length, chrom, colname )
                    peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]
            
            #tloop += ttime() - t0
            # save the last peak
            if not peak_content:
                continue
            else:
                self.__close_peak(peak_content, peaks, min_length, chrom, colname )

        #print "loop: %.2f" % tloop
        return peaks

    def __close_peak (self, peak_content, peaks, min_length, chrom, colname):
        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tmpsummit = []
            summit_pos   = None
            summit_value = None
            for (tmpstart,tmpend,tmpvalue,tmpsummitvalue, tmpindex) in peak_content:
                if not summit_value or summit_value < tmpsummitvalue:
                    tmpsummit = [ int(( tmpend+tmpstart )/2), ]
                    tmpsummit_index = [ tmpindex, ]
                    summit_value = tmpsummitvalue
                elif summit_value == tmpsummitvalue:
                    # remember continuous summit values
                    tmpsummit.append( int( (tmpend+tmpstart)/2 ) )
                    tmpsummit_index.append( tmpindex )
            middle_summit = int( ( len(tmpsummit)+1 )/2 )-1 # the middle of all highest points in peak region is defined as summit
            summit_pos    = tmpsummit[ middle_summit ]
            summit_index  = tmpsummit_index[ middle_summit ]
            # char * chromosome, long start, long end, long summit = 0, 
            # double peak_height=0, int pileup=0, 
            # double pvalue=0, double fold_change=0, double qvalue=0
            peaks.add( chrom,
                       peak_content[0][0],
                       peak_content[-1][1],
                       summit      = summit_pos,
                       peak_score  = self.data[chrom][colname][ summit_index ],
                       pileup      = self.data[chrom]['sample'][ summit_index ], # should be the same as summit_value
                       pscore      = self.data[chrom]['-100logp'][ summit_index ]/100.0,
                       fold_change = self.data[chrom]['sample'][ summit_index ]/self.data[chrom]['control'][ summit_index ],
                       qscore      = self.data[chrom]['-100logq'][ summit_index ]/100.0,
                       #peak_score  = chrom_score [ summit_index ],
                       #pileup      = chrom_sample[ summit_index ], # should be the same as summit_value
                       #pscore      = chrom_pvalue[ summit_index ]/100.0,
                       #fold_change = chrom_sample[ summit_index ]/chrom_control[ summit_index ],
                       #qscore      = chrom_qvalue[ summit_index ]/100.0,
                       )
            # start a new peak
            return True

    def call_broadpeaks (self, lvl1_cutoff=500, lvl2_cutoff=100, min_length=200, lvl1_max_gap=50, lvl2_max_gap=400, colname='-100logq'):
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
            for lvl2 in lvl2peakschrom:
                # for each lvl2 peak, find all lvl1 peaks inside
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
        
    def __add_broadpeak (self, bpeaks, chrom, lvl2peak, lvl1peakset):
        """Internal function to create broad peak.
        
        """
        start      = lvl2peak["start"]
        end        = lvl2peak["end"]
        thickStart = lvl1peakset[0]["start"]
        thickEnd   = lvl1peakset[-1]["end"]
        blockNum   = len(lvl1peakset)
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
        t = 0
        for chrom in self.data.keys():
            t += self.pointer[chrom]
        return t
