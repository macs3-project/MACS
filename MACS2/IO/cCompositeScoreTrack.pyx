# Time-stamp: <2012-08-01 18:09:29 Tao Liu>

"""Module for Composite Score Track IO classes.

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
from MACS2.IO.cPeakIO import PeakIO
from MACS2.IO.cBedGraph import bedGraphTrackI

# ------------------------------------
# constants
# ------------------------------------
__version__ = "scoreTrackI $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "scoreTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class compositeScoreTrackII:
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
                                                             ('sample1','float32'), 
                                                             ('sample2','float32'),
                                                             ('control1','float32'),
                                                             ('control2','float32'),                                                             
                                                             ('-100logq1','int32'),
                                                             ('-100logq2','int32')])
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
                fhd.write("%s\t%d\t%d\t%.5f\n" % (chrom,pre,pos[i],value[i]))
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
            above_cutoff_v = self.data[chrom][colname][above_cutoff] # scores where score is above cutoff

            above_cutoff_endpos = self.data[chrom]['pos'][above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = self.data[chrom]['pos'][above_cutoff-1] # start positions of regions where score is above cutoff
            above_cutoff_sv= self.data[chrom]['sample'][above_cutoff] # sample pileup height where score is above cutoff

            if above_cutoff_v.size == 0:
                continue

            if above_cutoff[0] == 0:
                above_cutoff_startpos[0] = 0

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


class compositeScoreTrackI:
    """Class for composite scoreGraph type data for two
    conditions. Modified from bedGraphTrackI. The only difference is
    that we store pvalue score, qvalue score and foldchange together.

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
                                                             ('sc1i1','float32'),
                                                             ('sc2i2','float32'),
                                                             ('sc1c2','float32'),
                                                             ('sc2c1','float32')])
            self.pointer[chrom] = 0

    def add (self,chromosome,endpos,sc1i1,sc2i2,sc1c2,sc2c1):
        """Add a chr-endpos-score-score-score-score block into data
        dictionary.

        """
        c = self.data[chromosome]
        i = self.pointer[chromosome]
        # get the preceding region
        c[i] = (endpos,sc1i1,sc2i2,sc1c2,sc2c1)
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

    def call_consistent (self, cutoff=50, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        Consistent peaks are those met all the following criteria

        1. sc1i1 >= cutoff
        2. sc2i2 >= cutoff
        3. sc1c2 <= cutoff
        4. sc2c1 <= cutoff
        """
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        #condition1_unique_peaks = PeakIO()                      # dictionary to save peaks
        #condition2_unique_peaks = PeakIO()                      # dictionary to save peaks        
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc1i1   = chrom_d[ 'sc1i1' ]
            chrom_sc2i2  = chrom_d[ 'sc2i2' ]
            chrom_sc1c2 = chrom_d[ 'sc1c2' ]
            chrom_sc2c1  = chrom_d[ 'sc2c1' ]

            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                vc1i1 = chrom_sc1i1[ x ]
                vc2i2 = chrom_sc2i2[ x ]
                vc1c2 = chrom_sc1c2[ x ]
                vc2c1 = chrom_sc2c1[ x ]                
                x += 1                  # index for the next point
                if vc1i1 >= cutoff and vc2i2 >= cutoff and vc1c2 <= cutoff and vc2c1 <= cutoff:
                    peak_content = [ ( pre_p, p, 0, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                vc1i1 = chrom_sc1i1[ i ]
                vc2i2 = chrom_sc2i2[ i ]
                vc1c2 = chrom_sc1c2[ i ]
                vc2c1 = chrom_sc2c1[ i ]
                if vc1i1 < cutoff or vc2i2 < cutoff or vc1c2 > cutoff or vc2c1 > cutoff:
                    pre_p = p
                    continue

                # for points met all criteria
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, 0, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = 0,
                                   peak_score  = 0,
                                   pileup      = 0,
                                   pscore      = 0,
                                   fold_change = 0,
                                   qscore      = 0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, 0, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = 0,
                           peak_score  = 0,
                           pileup      = 0,
                           pscore      = 0,
                           fold_change = 0,
                           qscore      = 0,
                           )
            
        return peaks

    def call_condition1_unique (self, cutoff=50, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        Condition 1 unique peaks are those met all the following criteria

        1. sc1i1 >= cutoff
        2. sc1c2 >= cutoff
        """
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc1i1   = chrom_d[ 'sc1i1' ]
            chrom_sc1c2 = chrom_d[ 'sc1c2' ]

            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                vc1i1 = chrom_sc1i1[ x ]
                vc1c2 = chrom_sc1c2[ x ]
                x += 1                  # index for the next point
                if vc1i1 >= cutoff and vc1c2 >= cutoff:
                    peak_content = [ ( pre_p, p, 0, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                vc1i1 = chrom_sc1i1[ i ]
                vc1c2 = chrom_sc1c2[ i ]
                if vc1i1 < cutoff or vc1c2 < cutoff:
                    pre_p = p
                    continue
                # for points met all criteria
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, 0, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = 0,
                                   peak_score  = 0,
                                   pileup      = 0,
                                   pscore      = 0,
                                   fold_change = 0,
                                   qscore      = 0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, 0, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = 0,
                           peak_score  = 0,
                           pileup      = 0,
                           pscore      = 0,
                           fold_change = 0,
                           qscore      = 0,
                           )
            
        return peaks

    def call_condition2_unique (self, cutoff=50, min_length=200, max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        Condition 2 unique peaks are those met all the following criteria

        1. sc2i2 >= cutoff
        2. sc2c1 >= cutoff
        """
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc2i2   = chrom_d[ 'sc2i2' ]
            chrom_sc2c1 = chrom_d[ 'sc2c1' ]
            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                vc2i2 = chrom_sc2i2[ x ]
                vc2c1 = chrom_sc2c1[ x ]
                x += 1                  # index for the next point
                if vc2i2 >= cutoff and vc2c1 >= cutoff:
                    peak_content = [ ( pre_p, p, 0, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                vc2i2 = chrom_sc2i2[ i ]
                vc2c1 = chrom_sc2c1[ i ]
                if vc2i2 < cutoff or vc2c1 < cutoff:
                    pre_p = p
                    continue
                # for points met all criteria
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, 0, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = 0,
                                   peak_score  = 0,
                                   pileup      = 0,
                                   pscore      = 0,
                                   fold_change = 0,
                                   qscore      = 0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, 0, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = 0,
                           peak_score  = 0,
                           pileup      = 0,
                           pscore      = 0,
                           fold_change = 0,
                           qscore      = 0,
                           )
        return peaks


    def call_diff_regions (self, cutoff=50, min_length=200, max_gap=50):
        """A function to call differential regions and common regions together.

        Return: common_regions, condition1_unique, and condition2_unique regions
        
        """
        chrs  = self.get_chr_names()
        consistent_peaks = PeakIO()                      # dictionary to save peaks
        condition1_peaks = PeakIO()                      # dictionary to save peaks
        condition2_peaks = PeakIO()                      # dictionary to save peaks        
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_sc1i1   = chrom_d[ 'sc1i1' ]
            chrom_sc2i2  = chrom_d[ 'sc2i2' ]
            chrom_sc1c2 = chrom_d[ 'sc1c2' ]
            chrom_sc2c1  = chrom_d[ 'sc2c1' ]

            x     = 0                   # index in compositeScoreTrackI
            pre_p = 0                   # remember previous position
            consistent_peak_content = None         # to store points above cutoff
            condition1_peak_content = None         # to store points above cutoff
            condition2_peak_content = None         # to store points above cutoff            
            
            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                vc1i1 = chrom_sc1i1[ i ]
                vc2i2 = chrom_sc2i2[ i ]
                vc1c2 = chrom_sc1c2[ i ]
                vc2c1 = chrom_sc2c1[ i ]

                if vc1i1 >= cutoff and vc2i2 >= cutoff and vc1c2 <= cutoff and vc2c1 <= cutoff:
                    # for points met all criteria
                    # if the gap is allowed
                    if not consistent_peak_content:
                        consistent_peak_content = [ ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ), ] # for consistent region, summit is decided by condition 1
                    if pre_p - consistent_peak_content[ -1 ][ 1 ] <= max_gap:
                        consistent_peak_content.append( ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ) )
                    else:
                        # when the gap is not allowed, close this peak
                        # this is common region.
                        
                        peak_length = consistent_peak_content[ -1 ][ 1 ] - consistent_peak_content[ 0 ][ 0 ]
                        if peak_length >= min_length: # if the peak is too small, reject it

                            # the absolute of peak score or the diff
                            # score is the maximum of max[vc1c2,] and
                            # max[vc2c1,]. If max[vc1c2,] is bigger,
                            # the sign for peak score is
                            # '+'. Otherwise, the sign is '-'

                            m_c1c2 = max([x[4] for x in consistent_peak_content ])
                            m_c2c1 = max([x[5] for x in consistent_peak_content ])                            
                            if m_c1c2 >= m_c2c1:
                                diff_score = m_c1c2
                            else:
                                diff_score = -1* m_c2c1

                            consistent_peaks.add( chrom,
                                                  consistent_peak_content[0][0],
                                                  consistent_peak_content[-1][1],
                                                  summit      = 0,
                                                  peak_score  = diff_score,
                                                  pileup      = 0,
                                                  pscore      = 0,
                                                  fold_change = 0,
                                                  qscore      = 0,
                                                  )
                        # start a new peak
                        consistent_peak_content = [ ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ), ]
                elif vc1i1 >= cutoff and vc1c2 >= cutoff:
                    if not condition1_peak_content:
                        condition1_peak_content = [ ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ),  ]
                    if pre_p - condition1_peak_content[ -1 ][ 1 ] <= max_gap:
                        condition1_peak_content.append( ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i )  )
                    else:
                        # when the gap is not allowed, close this peak
                        # this is condition1 unique region
                        peak_length = condition1_peak_content[ -1 ][ 1 ] - condition1_peak_content[ 0 ][ 0 ]
                        if peak_length >= min_length: # if the peak is too small, reject it
                            # the absolute of peak score or the diff
                            # score is the maximum of max[vc1c2,] and
                            # max[vc2c1,]. If max[vc1c2,] is bigger,
                            # the sign for peak score is
                            # '+'. Otherwise, the sign is '-'

                            diff_score = max([x[4] for x in condition1_peak_content ])
                            #m_c2c1 = max([x[5] in condition2_peak_content ])                            
                            #if m_c1c2 >= m_c2c1:
                            ##    diff_score = m_c1c2
                            #else:
                            #    diff_score = -1* m_c2c1

                            condition1_peaks.add( chrom,
                                                  condition1_peak_content[0][0],
                                                  condition1_peak_content[-1][1],
                                                  summit      = 0,
                                                  peak_score  = diff_score,
                                                  pileup      = 0,
                                                  pscore      = 0,
                                                  fold_change = 0,
                                                  qscore      = 0,
                                                  )
                        # start a new peak
                        condition1_peak_content = [ ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ), ]
                elif vc2i2 >= cutoff and vc2c1 >= cutoff:
                    if not condition2_peak_content:
                        condition2_peak_content = [ ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ), ]
                    if pre_p - condition2_peak_content[ -1 ][ 1 ] <= max_gap:
                        condition2_peak_content.append( ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ) )
                    else:
                        # when the gap is not allowed, close this peak
                        # condition 2 unique peaks
                        peak_length = condition2_peak_content[ -1 ][ 1 ] - condition2_peak_content[ 0 ][ 0 ]
                        if peak_length >= min_length: # if the peak is too small, reject it
                            diff_score = -1 * max([x[5] for x in condition2_peak_content ])
                            condition2_peaks.add( chrom,
                                                  condition2_peak_content[0][0],
                                                  condition2_peak_content[-1][1],
                                                  summit      = 0,
                                                  peak_score  = diff_score,
                                                  pileup      = 0,
                                                  pscore      = 0,
                                                  fold_change = 0,
                                                  qscore      = 0,
                                                  )
                        # start a new peak
                        condition2_peak_content = [ ( pre_p, p, vc1i1, vc2i2, vc1c2, vc2c1, i ), ]
                pre_p = p
            
            # save the last regions
            if consistent_peak_content:
                peak_length = consistent_peak_content[ -1 ][ 1 ] - consistent_peak_content[ 0 ][ 0 ]
                if peak_length >= min_length: # if the peak is too small, reject it
                    m_c1c2 = max([x[4] for x in consistent_peak_content ])
                    m_c2c1 = max([x[5] for x in consistent_peak_content ])                            
                    if m_c1c2 >= m_c2c1:
                        diff_score = m_c1c2
                    else:
                        diff_score = -1* m_c2c1                    
                    consistent_peaks.add( chrom,
                                          consistent_peak_content[0][0],
                                          consistent_peak_content[-1][1],
                                          summit      = 0,
                                          peak_score  = diff_score,
                                          pileup      = 0,
                                          pscore      = 0,
                                          fold_change = 0,
                                          qscore      = 0,
                                          )
            elif condition1_peak_content:
                peak_length = condition1_peak_content[ -1 ][ 1 ] - condition1_peak_content[ 0 ][ 0 ]
                if peak_length >= min_length: # if the peak is too small, reject it
                    diff_score = max([x[4] for x in condition1_peak_content ])
                    condition1_peaks.add( chrom,
                                          condition1_peak_content[0][0],
                                          condition1_peak_content[-1][1],
                                          summit      = 0,
                                          peak_score  = diff_score,
                                          pileup      = 0,
                                          pscore      = 0,
                                          fold_change = 0,
                                          qscore      = 0,
                                          )
            elif condition2_peak_content:
                peak_length = condition2_peak_content[ -1 ][ 1 ] - condition2_peak_content[ 0 ][ 0 ]
                if peak_length >= min_length: # if the peak is too small, reject it
                    diff_score = -1 * max([x[5] for x in condition2_peak_content ])                    
                    condition2_peaks.add( chrom,
                                          condition2_peak_content[0][0],
                                          condition2_peak_content[-1][1],
                                          summit      = 0,
                                          peak_score  = 0,
                                          pileup      = 0,
                                          pscore      = 0,
                                          fold_change = 0,
                                          qscore      = 0,
                                          )

        return ( consistent_peaks, condition1_peaks, condition2_peaks )

    def total ( self ):
        """Return the number of regions in this object.

        """
        t = 0
        for chrom in self.data.keys():
            t += self.pointer[chrom]
        return t

    #def dump ( self ):
    #
    #

def make_compositeScoreTrack (bdgTrack1, bdgTrack2, bdgTrack3, bdgTrack4 ):
    """A modified overlie function for MACS DIFF.
    
    """
    assert isinstance(bdgTrack1,bedGraphTrackI), "bdgTrack1 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack3,bedGraphTrackI), "bdgTrack3 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack4,bedGraphTrackI), "bdgTrack4 is not a bedGraphTrackI object"    
    
    ret = compositeScoreTrackI()
    retadd = ret.add

    chr1 = set(bdgTrack1.get_chr_names())
    chr2 = set(bdgTrack2.get_chr_names())
    chr3 = set(bdgTrack3.get_chr_names())
    chr4 = set(bdgTrack4.get_chr_names())    
    
    common_chr = chr1.intersection(chr2).intersection(chr3).intersection(chr4)
    for chrom in common_chr:
            
        (p1s,v1s) = bdgTrack1.get_data_by_chr(chrom) # arrays for position and values
        p1n = iter(p1s).next         # assign the next function to a viable to speed up
        v1n = iter(v1s).next

        (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
        p2n = iter(p2s).next         # assign the next function to a viable to speed up
        v2n = iter(v2s).next

        (p3s,v3s) = bdgTrack3.get_data_by_chr(chrom) # arrays for position and values
        p3n = iter(p3s).next         # assign the next function to a viable to speed up
        v3n = iter(v3s).next

        (p4s,v4s) = bdgTrack4.get_data_by_chr(chrom) # arrays for position and values
        p4n = iter(p4s).next         # assign the next function to a viable to speed up
        v4n = iter(v4s).next

        chrom_max_len = len(p1s)+len(p2s)+len(p3s)+len(p4s) # this is the maximum number of locations needed to be recorded in scoreTrackI for this chromosome.
            
        ret.add_chromosome(chrom,chrom_max_len)

        pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
        try:
            p1 = p1n()
            v1 = v1n()

            p2 = p2n()
            v2 = v2n()

            p3 = p3n()
            v3 = v3n()

            p4 = p4n()
            v4 = v4n()
            
            while True:
                min_p = min( p1, p2, p3, p4 )
                retadd( chrom, min_p, v1, v2, v3, v4 )
                pre_p = min_p

                if p1 == min_p:
                    p1 = p1n()
                    v1 = v1n()
                if p2 == min_p:
                    p2 = p2n()
                    v2 = v2n()
                if p3 == min_p:
                    p3 = p3n()
                    v3 = v3n()
                if p4 == min_p:
                    p4 = p4n()
                    v4 = v4n()                                        
        except StopIteration:
            # meet the end of either bedGraphTrackI, simply exit
            pass
        
    return ret



def make_compositeScoreTrack2 (bdgTrack1, bdgTrack2, bdgTrack3, bdgTrack4 ):
    """A modified overlie function for MACS DIFF.
    
    """
    assert isinstance(bdgTrack1,bedGraphTrackI), "bdgTrack1 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack2,bedGraphTrackI), "bdgTrack2 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack3,bedGraphTrackI), "bdgTrack3 is not a bedGraphTrackI object"
    assert isinstance(bdgTrack4,bedGraphTrackI), "bdgTrack4 is not a bedGraphTrackI object"    
    
    ret = compositeScoreTrackII()
    retadd = ret.add

    score_btrack1 = bdgTrack1.make_scoreTrack_for_macs( bdgTrack3 )
    pqtable1 = score_btrack1.make_pq_table()
    score_btrack1.assign_qvalue(pqtable1)
    
    score_btrack2 = bdgTrack2.make_scoreTrack_for_macs( bdgTrack4 )
    pqtable2 = score_btrack2.make_pq_table()
    score_btrack2.assign_qvalue(pqtable2)    

    
    
    common_chr = chr1.intersection(chr2).intersection(chr3).intersection(chr4)
    for chrom in common_chr:
            
        (p1s,v1s) = bdgTrack1.get_data_by_chr(chrom) # arrays for position and values
        p1n = iter(p1s).next         # assign the next function to a viable to speed up
        v1n = iter(v1s).next

        (p2s,v2s) = bdgTrack2.get_data_by_chr(chrom) # arrays for position and values
        p2n = iter(p2s).next         # assign the next function to a viable to speed up
        v2n = iter(v2s).next

        (p3s,v3s) = bdgTrack3.get_data_by_chr(chrom) # arrays for position and values
        p3n = iter(p3s).next         # assign the next function to a viable to speed up
        v3n = iter(v3s).next

        (p4s,v4s) = bdgTrack4.get_data_by_chr(chrom) # arrays for position and values
        p4n = iter(p4s).next         # assign the next function to a viable to speed up
        v4n = iter(v4s).next

        chrom_max_len = len(p1s)+len(p2s)+len(p3s)+len(p4s) # this is the maximum number of locations needed to be recorded in scoreTrackI for this chromosome.
            
        ret.add_chromosome(chrom,chrom_max_len)

        pre_p = 0                   # remember the previous position in the new bedGraphTrackI object ret
            
        try:
            p1 = p1n()
            v1 = v1n()

            p2 = p2n()
            v2 = v2n()

            p3 = p3n()
            v3 = v3n()

            p4 = p4n()
            v4 = v4n()
            
            while True:
                min_p = min( p1, p2, p3, p4 )
                retadd( chrom, min_p, v1, v2, v3, v4 )
                pre_p = min_p

                if p1 == min_p:
                    p1 = p1n()
                    v1 = v1n()
                if p2 == min_p:
                    p2 = p2n()
                    v2 = v2n()
                if p3 == min_p:
                    p3 = p3n()
                    v3 = v3n()
                if p4 == min_p:
                    p4 = p4n()
                    v4 = v4n()                                        
        except StopIteration:
            # meet the end of either bedGraphTrackI, simply exit
            pass
        
    return ret

