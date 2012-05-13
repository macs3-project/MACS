# cython: profile=True
# Time-stamp: <2012-05-10 19:05:19 Tao Liu>

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

from array import array as pyarray

from cpython cimport bool
#from scipy.signal import fftconvolve
from MACS2.cSignal import maxima, enforce_valleys
#np_convolve = np.convolve

from libc.math cimport log10,log

from MACS2.Constants import *
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO, BroadPeakIO

from MACS2.hashtable import Int64HashTable

import logging

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

LOG10_E = 0.43429448190325176
pscore_khashtable = Int64HashTable()

cdef int get_pscore ( int observed, double expectation ):
    """Get p-value score from Poisson test. First check existing
    table, if failed, call poisson_cdf function, then store the result
    in table.
    
    """
    cdef:
        int score
        long key_value
    
    #key_value = ( observed, expectation )
    key_value = hash( (observed, expectation ) )
    try:
        return pscore_khashtable.get_item(key_value)
    except KeyError:
        score = int(-100*poisson_cdf(observed,expectation,False,True))
        pscore_khashtable.set_item(key_value, score)
        return score
    #if pscore_dict.has_key(key_value):
    #    return pscore_dict[key_value]
    #else:
    #    score = int(-100*poisson_cdf(observed,expectation,False,True))
    #    pscore_dict[(observed,expectation)] = score
    #return score

logLR_khashtable = Int64HashTable()

cdef int logLR ( double x, double y ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Then store the values in integer form =
    100*logLR. Set minus sign for depletion.
    
    """
    cdef:
        int s
        long key_value
    
    #key_value = ( x, y )
    key_value = hash( (x, y ) )
    try:
        return logLR_khashtable.get_item( key_value )
    except KeyError:
        if x > y:
            s = int( (x*(log(x+1)-log(y+1))+y-x)*LOG10_E*100 )
        elif x < y:
            s = int( (-1*x*(log(x+1)-log(y+1))-y+x)*LOG10_E*100 )
        else:
            s = 0
        logLR_khashtable.set_item(key_value, s)
        return s

cdef int get_foldenrichment ( float x, float y ):
    """ return fold enrichment with +1 pseudocount.
    """
    return int( (x+1)/(y+1) )

cdef int get_substraction ( float x, float y):
    """ return substraction.
    """
    return int( x - y )

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
    def __init__ (self, float effective_depth_in_million = 1.0):
        """Different with bedGraphTrackI, missing values are simply
        replaced with 0.

        effective_depth_in_million: sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.
        """
        self.data = {}
        self.pointer = {}
        self.trackline = False
        self.effective_depth_in_million = effective_depth_in_million
        
    def enable_trackline(self):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True

    def add_chromosome ( self, str chrom, int chrom_max_len ):
        if not self.data.has_key(chrom):
            self.data[chrom] = { 'pos': np.zeros(chrom_max_len, dtype="int32"),
                                 'sample': np.zeros(chrom_max_len, dtype="float32"),
                                 'control': np.zeros(chrom_max_len, dtype="float32"),
                                 '-100logp': np.zeros(chrom_max_len, dtype="int32"),
                                 '-100logq': np.zeros(chrom_max_len, dtype="int32"),
                                 '100logLR': np.zeros(chrom_max_len, dtype="int32") }
            self.pointer[chrom] = 0

    def add (self, str chromosome, int endpos, float sample, float control):
        """Add a chr-endpos-sample-control block into data
        dictionary. At the mean time, calculate pvalues and log
        likelihood.

        """
        cdef int i
        #print chromosome, endpos, sample, control
        i = self.pointer[chromosome]
        c = self.data[chromosome]
        c['pos'][i] = endpos
        c['sample'][i] = sample
        c['control'][i] = control
        c['-100logp'][i] = get_pscore(int(sample),control)
        c['100logLR'][i] = logLR(sample,control)        
        self.pointer[chromosome] += 1

    def finalize (self):
        cdef:
            str chrom, k
            int l

        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.pointer[chrom]
            for k in d.keys():
                d[k].resize(l,refcheck=False)

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

    def write_bedGraph ( self, fhd, str name, str description, str colname, bool do_SPMR = False ):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 'sample','control','-100logp','-100logq', '100logLR'

        do_SPMR: only effective when writing sample/control tracks. When True, save SPMR instead.
        
        """
        cdef:
            str chrom
            int l, pre, i, p 
            float pre_v, v, scale_factor
            
        if self.trackline:
            # this line is REQUIRED by the wiggle format for UCSC browser
            fhd.write( "track type=bedGraph name=\"%s\" description=\"%s\"\n" % ( name,description ) )
        
        if colname not in [ 'sample', 'control', '-100logp', '-100logq', '100logLR' ]:
            raise Exception( "%s not supported!" % colname )

        if colname in [ '-100logp', '-100logq', '100logLR' ]:
            scale_factor = 0.01              # for pvalue or qvalue, divide them by 100 while writing to bedGraph file
        elif colname in [ 'sample', 'control' ]:
            if do_SPMR:
                logging.info( "MACS will save SPMR for fragment pileup using effective depth of %.2f million" % self.effective_depth_in_million )
                scale_factor = 1.0/self.effective_depth_in_million
            else:
                scale_factor = 1
        
        chrs = self.get_chr_names()
        for chrom in chrs:
            d = self.data[ chrom ]
            l = self.pointer[ chrom ]
            pre = 0
            pos   = d[ 'pos' ]
            value = d[ colname ] * scale_factor
            if value.shape[ 0 ] == 0: continue # skip if there's no data
            pre_v = value[ 0 ]
            for i in range( 1, l ):
                v = value[ i ]
                p = pos[ i-1 ]
                if pre_v != v: 
                    fhd.write( "%s\t%d\t%d\t%.2f\n" % ( chrom, pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            fhd.write( "%s\t%d\t%d\t%.2f\n" % ( chrom, pre, p, pre_v ) )
            
        return True

    def make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-100log10pvalue:(-100log10qvalue,rank,basepairs)} relationships.
        """
        cdef:
            long n, pre_p, this_p, length, j, pre_l, l, this_v, pre_v, v
            long N, k, q, pre_q
            double f
            str chrom
        
        #logging.info("####test#### start make_pq")
        n = self.total()
        value_dict = Int64HashTable()
        unique_values = pyarray(BYTE4,[])
        # this is a table of how many positions each p value occurs at
        for chrom in self.data.keys():
            # for each chromosome
            pre_p  = 0
            pos    = self.data[chrom][ 'pos' ]
            value  = self.data[chrom][ '-100logp' ]
            length = self.pointer[chrom]
            for j in xrange(length):
                this_p = pos[j]
                this_v = value[j]
                assert this_v == this_v, "NaN at %d" % pos
                if value_dict.has_key(this_v):
                    value_dict.set_item(this_v, value_dict.get_item(this_v) + this_p - pre_p)
                else:
                    value_dict.set_item(this_v, this_p - pre_p)
                    unique_values.append(this_v)
                pre_p = this_p

        N = 0
        for i in xrange(len(unique_values)):
            N += value_dict.get_item(unique_values[i])
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        pvalue2qvalue = {pre_v:[0,k,0]}              # pvalue:[qvalue,rank,bp_with_this_pvalue]
        #pvalue2qvalue = np.zeros( (len(unique_values)+1,4), dtype='int64' )
        #pvalue2qvalue[0] = (pre_v, 0, k, 0)
        #logging.info("####test#### start matching pvalue to qvalue")
        unique_values = sorted(unique_values,reverse=True)
        for i in xrange(len(unique_values)):
            v = unique_values[i]
            l = value_dict.get_item(v)
            q = v + int((log10(k) + f) * 100) # we save integers here.
            q = max(0,min(pre_q,q))           # make q-score monotonic
            #pvalue2qvalue[i+1] = (v, q, k, 0)
            #pvalue2qvalue[i][3] = k - pvalue2qvalue[i][2]
            pvalue2qvalue[v] = [q, k, 0]
            pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
            pre_v = v
            pre_q = q
            k+=l
        #pvalue2qvalue[i+1][3] = k - pvalue2qvalue[i][2]
        pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
        #logging.info("####test#### finish building pqtable")        
        # pop the first -1e100 one
        pvalue2qvalue.pop(-2147483647)
        #pvalue2qvalue = pvalue2qvalue[1:]

        return pvalue2qvalue

    def assign_qvalue ( self , dict pvalue2qvalue ):
        """Assign -100log10qvalue to every point.

        pvalue2qvalue: a dictionary of -100log10pvalue:-100log10qvalue
        """
        cdef:
            long i,l,j,p
            str chrom
            
        chroms = self.data.keys()

        # convert pvalue2qvalue to a simple dict
        s_p2q = Int64HashTable()
        #g = pvalue2qvalue.get
        for i in pvalue2qvalue.keys():
        #for i in range(pvalue2qvalue.shape[0]):
            s_p2q.set_item(i,pvalue2qvalue[i][0])
            #s_p2q.set_item(pvalue2qvalue[i][0],pvalue2qvalue[i][1])

        g = s_p2q.get_item
        
        for j in range( len(chroms) ):
            chrom = chroms[j]
            pvalue = self.data[chrom]['-100logp']
            qvalue = self.data[chrom]['-100logq']
            l = self.pointer[chrom]
            for i in range( l ):
                qvalue[i] = g(pvalue[i])
        return True

    def call_peaks (self, int cutoff=500, int min_length=200, int max_gap=50, str colname='-100logp',
                    bool call_summits=False, int smoothlen=51):
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
        cdef:
            int i
            str chrom
        
        assert (colname in [ 'sample', 'control', '-100logp', '-100logq', '100logLR' ]), "%s not supported!" % colname

        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        if call_summits: close_peak = self.__close_peak2
        else: close_peak = self.__close_peak
        
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
                    close_peak(peak_content, peaks, min_length, chrom, colname, smoothlen=smoothlen )
                    peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]
            
            # save the last peak
            if not peak_content:
                continue
            else:
                close_peak(peak_content, peaks, min_length, chrom,
                           colname, smoothlen=max_gap/2 )

        return peaks
       
    def __close_peak2 (self, peak_content, peaks, int min_length, str chrom, str colname, int smoothlen=51,
                       float min_valley = 0.8):
        cdef:
            int summit_pos, tstart, tend, tmpindex, summit_index, summit_offset
            int start, end, i, j
            double summit_value, tvalue, tsummitvalue
#            np.ndarray[np.float32_t, ndim=1] w
            np.ndarray[np.float32_t, ndim=1] peakdata
            np.ndarray[np.int32_t, ndim=1] peakindices, valid_summit_indices, summit_offsets, valid_summit_offsets
            
        # this is where the summits are called, need to fix this
        end = peak_content[ -1 ][ 1 ]
        start = peak_content[ 0 ][ 0 ]
        peak_length = end - start
        if end - start < min_length: return # if the region is too small, reject it

        peakdata = np.zeros(end - start, dtype='float32')
        peakindices = np.zeros(end - start, dtype='int32')
        for (tstart,tend,tvalue,tsummitvalue, tmpindex) in peak_content:
            i = tstart - start
            j = tend - start
            peakdata[i:j] = self.data[chrom]['sample'][tmpindex]
            peakindices[i:j] = tmpindex
        # apply smoothing window of smoothlen
#        w = np.ones(smoothlen, dtype='float32') / smoothlen
#        if smoothlen > 0:
#            smoothdata = np_convolve(w, peakdata, mode='same')
#        else:
#            smoothdata = peakdata.copy()
        summit_offsets = maxima(peakdata, smoothlen)
        # ***failsafe if no summits so far*** #
        if summit_offsets.shape[0] == 0:
            summit_offsets = np.asarray([peak_length / 2], dtype='int32')
        valid_summit_offsets = enforce_valleys(peakdata, summit_offsets, min_valley = min_valley)
        valid_summit_indices = peakindices[valid_summit_offsets]
        # this case shouldn't occur anymore because we've disallowed plateaus
        # purge offsets that have the same summit_index
#        unique_offsets = []
#        for index in np.unique(summit_indices):
#            those_index_indices = np.where(summit_indices == index)[0]
#            those_offsets = summit_offsets[those_index_indices]
#            unique_offsets.append(int(those_offsets.mean()))
           
        ## DISABLE PEAKSPLITTER BEHAVIOR FOR NOW ##
        # I think requiring spatial consistency is better than requiring valleys
        
        # also require a valley of at least 0.6 * taller peak
        # in every adjacent two peaks or discard the lesser one
        # this behavior is like PeakSplitter
#        better_offsets = []
#        previous_offset = unique_offsets.pop()
#        while True:
#            if len(unique_offsets) == 0:
#                better_offsets.append(previous_offset)
#                break
#            else:
#                this_offset = unique_offsets.pop()
#                this_h, prev_h = peakdata[[this_offset, previous_offset]]
#                if this_h > prev_h:
#                    prev_is_taller = False
#                    min_valley = 0.6 * this_h
#                else:
#                    prev_is_taller = True
#                    min_valley = 0.6 * prev_h
#                s = slice(this_offset, previous_offset)
#                valley = np.where(peakdata[s] < min_valley)[0]
#                if len(valley) > 0: better_offsets.append(previous_offset)
#                else:
#                    if prev_is_taller: continue # discard this peak
#                    # else: discard previous peak by ignoring it
#                previous_offset = this_offset
#        better_offsets.reverse()
#        better_indices = peakindices[better_offsets]
#        assert len(better_offsets) > 0, "Lost peak summit(s) near %s %d" % (chrom, start) 
#        for summit_offset, summit_index in zip(better_offsets, better_indices):
        for summit_offset, summit_index in zip(valid_summit_offsets, valid_summit_indices):
            peaks.add( chrom,
                       start,
                       end,
                       summit      = start + summit_offset,
                       peak_score  = self.data[chrom][colname][ summit_index ],
                       pileup      = self.data[chrom]['sample'][ summit_index ], # should be the same as summit_value
                       pscore      = self.data[chrom]['-100logp'][ summit_index ]/100.0,
                       fold_change = self.data[chrom]['sample'][ summit_index ]/self.data[chrom]['control'][ summit_index ],
                       qscore      = self.data[chrom]['-100logq'][ summit_index ]/100.0,
                       )
        # start a new peak
        return True

    def __close_peak (self, peak_content, peaks, int min_length, str chrom, str colname, int smoothlen=0):
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        """
        cdef:
            int summit_pos, tstart, tend, tmpindex, summit_index, i, midindex
            double summit_value, tvalue, tsummitvalue

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tsummit = []
            summit_pos   = 0
            summit_value = summit_pos
            for i in xrange(len(peak_content)):
                (tstart,tend,tvalue,tsummitvalue, tindex) = peak_content[i]
                #for (tstart,tend,tvalue,tsummitvalue, tindex) in peak_content:
                if not summit_value or summit_value < tsummitvalue:
                    tsummit = [(tend + tstart) / 2, ]
                    tsummit_index = [ tindex, ]
                    summit_value = tsummitvalue
                elif summit_value == tsummitvalue:
                    # remember continuous summit values
                    tsummit.append(int((tend + tstart) / 2))
                    tsummit_index.append( tindex )
            # the middle of all highest points in peak region is defined as summit
            midindex = int((len(tsummit) + 1) / 2) - 1
            summit_pos    = tsummit[ midindex ]
            summit_index  = tsummit_index[ midindex ]
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
        cdef:
            int i
            str chrom
        
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
        cdef:
            int blockNum, thickStart, thickEnd, start, end
            str blockSizes, blockStarts
        
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
        cdef:
            long t
            str chrom
        
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
        cdef:
            list c
            int i
            
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
        cdef:
            str chrom
            set chrs
            int pre, i, l
        
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
        cdef:
            long t
            str chrom
            
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
        cdef:
            set chr1, chr2, common_chr
            str chrom
            int pre_p, p1, p2
            double v11, v21, v2
        
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
        cdef:
            int i, l
            str chrom, start, end
        
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

    def extract_sum (self, bdgTrack2):
        """Get sum values in each region defined in bdgTrack2.
        
        """
        cdef:
            int i
            str chrom, start, end
        
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
                    radd(cur_region[0]+"."+str(cur_region[1])+"."+str(cur_region[2]))
                    v1add(cur_region[3])
                    v2add(cur_region[4])                    
                cur_region = [chrom, start, end, v1array[i]*larray[i], v2array[i]*larray[i]]

        radd(cur_region[0]+"."+str(cur_region[1])+"."+str(cur_region[2]))
        v1add(cur_region[3])
        v2add(cur_region[4])
        return ret

#cdef find_maxima(np.ndarray[np.float32_t, ndim=1] smoothdata, int min_length):
#    cdef:
#        float peaksplit_p = 0.6
#        int pos, this_peak_pos
#        int n_subpeak = 0
#        bool going_up = True
#        bool direction = (smoothdata[1] > smoothdata[0])
#        bool new_direction
#        int fw_i, rev_i  # how long we're going the right/wrong direction
#        int plateau_i = 0
#        int last_peak_pos = -1
#        float last_peak_value, this_peak_value
#        int last_peak_i = 0
#        int half_l = min_length / 2
#        int quarter_l = min_length / 4
#        np.ndarray[np.int32_t, ndim=1] ret = np.zeros(smoothdata.shape[0], 'int32')
#        int rsize = ret.shape[0]
#        
#    for pos in range(1, smoothdata.shape[0] - 1):
#        # find new_direction
#        if smoothdata[pos + 1] > smoothdata[pos]:
#            new_direction = True
#        elif smoothdata[pos + 1] < smoothdata[pos]:
#            new_direction = False
#        else: # unlikely for floats
#            # we're on a plateau
#            plateau_i += 1
#            if plateau_i > half_l: # don't allow very long plateaus
#                fw_i = 0
#                rev_i = 0
#                last_peak_pos = -1
#            continue
#       
#        # increment fw/rev counter
#        if direction == going_up:
#            fw_i += 1
#        else:
#            rev_i += 1
#       
#        # decide if we're still going up or down 
#        if rev_i > fw_i:
#            last_peak_pos = -1
#            fw_i = rev_i
#            rev_i = 0
##            print 'reversed at %d' % pos
#            going_up = not going_up
#                
#        if (not direction) and new_direction: # decreasing -> increasing
#            direction = new_direction
#            # minimum (decreasing to increasing)
#            if not going_up and not last_peak_pos == -1:               
#                if fw_i >= quarter_l:                # right side long enough                
#                    # try calling a peak
#                    if last_peak_i >= quarter_l:       # left side long enough
#                        # save this peak
#                        ret[n_subpeak] = last_peak_pos
#                        last_peak_pos = -1
#                        n_subpeak += 1
#                        fw_i = 0
#                        rev_i = 0
#                        going_up = True
#                    else:                             # left side too short
#                        last_peak_pos = -1            # discard this peak
#        elif direction and (not new_direction):   # increasing -> decreasing
#            direction = new_direction
#            # maxima plateaus need to be centered
#            this_peak_pos = pos - plateau_i / 2
#            plateau_i = 0
#            # maximum (increasing to decreasing)
#            if not going_up:
#                continue
#            this_peak_value = smoothdata[this_peak_pos]
#            if last_peak_pos == -1 or this_peak_value > last_peak_value:
#                if last_peak_pos == -1: last_peak_i = fw_i
#                else: last_peak_i += fw_i
#                last_peak_pos = this_peak_pos
#                last_peak_value = this_peak_value
#                fw_i = 0
#                rev_i = 0
#                going_up = False
#            else: # this is not a peak at all
#                continue
#        plateau_i = 0
#            
#    if not going_up and fw_i >= quarter_l and not last_peak_pos == -1:
#            ret[n_subpeak] = last_peak_pos
#            n_subpeak += 1
#    ret.resize(n_subpeak, refcheck = False)
#    return ret

cdef class scoreTrackII:
    """Class for scoreGraph type data. Modified from scoreTrackI. The
    difference is that we store a single score data, not
    p/q/loglikelihood altogether. Score (the 4th) column is calculated
    through calling change_method() function. Individual scores for
    filling PeakIO object are calculated by prepare_peakIO_scores()
    function.

    I want to save mem and simplify calculation in this new class.

    """
    cdef:
        dict data                       # dictionary for data of each chromosome
        dict data_stderr                # dictionary for data stderr for each chromosome
        dict datalength                 # length of data array of each chromosome
        bool trackline                  # whether trackline should be saved in bedGraph
        bool stderr_on                  # whether to calculate stderr
        double treat_edm                 # seq depth in million of treatment
        double ctrl_edm                  # seq depth in million of control
        char scoring_method              # method for calculating scores.
        char normalization_method        # scale to control? scale to treatment? both scale to 1million reads?

    
    def __init__ (self, float treat_depth, float ctrl_depth, bool stderr_on = False ):
        """Initialize.

        effective_depth_in_million: sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        
        """
        self.data = {}           # for each chromosome, there is a l*4
                                 # matrix. First column: end position
                                 # of a region; Second: treatment
                                 # pileup * 100; third: control pileup
                                 # * 100; forth: score * 100 ( can be
                                 # p/q-value/likelihood
                                 # ratio/fold-enrichment/substraction
                                 # depending on -c setting)
        self.stderr_on = stderr_on
        self.datalength = {}
        self.trackline = False
        self.treat_edm = treat_depth
        self.ctrl_edm = ctrl_depth
        #scoring_method:  p: -log10 pvalue;
        #                 q: -log10 qvalue;
        #                 l: log10 likelihood ratio ( minus for depletion )
        #                 f: log10 fold enrichment
        #                 d: substraction
        #                 m: fragment pileup per million reads
        #                 N: not set
        self.scoring_method = ord("N")

        #normalization_method: T: scale to depth of treatment;
        #                      C: scale to depth of control;
        #                      M: scale to depth of 1 million;
        #                      N: not set/ raw pileup
        self.normalization_method = ord("N")
        
    cpdef enable_trackline(self):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True

    cpdef add_chromosome ( self, str chrom, int chrom_max_len ):
        """
        chrom: chromosome name
        chrom_max_len: maximum number of data points in this chromosome
        
        """
        if not self.data.has_key(chrom):
            self.data[chrom] = np.zeros( ( chrom_max_len, 4 ), dtype="int32" ) # remember col #2-4 is actual value * 100, I use integer here.
            self.datalength[chrom] = 0

    cpdef add (self, str chromosome, int endpos, int sample, int control):
        """Add a chr-endpos-sample-control block into data
        dictionary.

        *Warning* Need to add regions continuously.
        """
        cdef int i
        i = self.datalength[chromosome]
        c = self.data[chromosome]
        c[ i, 0 ] = endpos
        c[ i, 1 ] = sample * 100
        c[ i, 2 ] = control * 100
        self.datalength[chromosome] += 1

    cpdef finalize ( self ):
        """
        Adjust array size of each chromosome.

        """
        cdef:
            str chrom, k
            int l

        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for k in d.keys():
                d[k].resize( (l,4), refcheck = False )
        return

    cpdef sort ( self, int column = 1 ):
        """ Sort data for each chromosome, by certain column.

        column: 1: position, 2: sample, 3: control, 4: score

        Default: sort by positions.
        """
        for chrom in self.data.keys():
            d = self.data[chrom]
            d.view('int32,int32,int32,int32').sort(axis=0,order=column-1)
        return

    cpdef get_data_by_chr (self, str chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if self.data.has_key(chromosome):
            return self.data[chromosome]
        else:
            return None

    cpdef get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.data.keys())
        return l

    cpdef change_normalization_method ( self, char normalization_method ):
        """Change/set normalization method. However, I do not
        recommend change this back and forward, since some precision
        issue will happen -- I only keep two digits.
        
        normalization_method: T: scale to depth of treatment;
                             C: scale to depth of control;
                             M: scale to depth of 1 million;
                             N: not set/ raw pileup        
        """
        if normalization_method == 'T':
            if self.normalization_method == 'T': # do nothing
                pass
            elif self.normalization_method == 'C':
                self.normalize( self.treat_edm/self.ctrl_edm, self.treat_edm/self.ctrl_edm )
            elif  self.normalization_method == 'M':
                self.normalize( self.treat_edm, self.treat_edm )
            elif self.normalization_method == 'N':
                self.normalize( 1, self.treat_edm/self.ctrl_edm )
            else:
                raise NotImplemented
            self.normalization_method = 'T'
        elif normalization_method == 'C':
            if self.normalization_method == 'T':
                self.normalize( self.ctrl_edm/self.treat_edm, self.ctrl_edm/self.treat_edm )
            elif self.normalization_method == 'C': # do nothing
                pass
            elif  self.normalization_method == 'M':
                self.normalize( self.ctrl_edm, self.ctrl_edm )
            elif self.normalization_method == 'N':
                self.normalize( self.ctrl_edm/self.treat_edm, 1 )
            else:
                raise NotImplemented
            self.normalization_method = 'C'                
        elif normalization_method == 'M':
            if self.normalization_method == 'T':
                self.normalize( 1/self.treat_edm, 1/self.treat_edm )
            elif self.normalization_method == 'C':
                self.normalize( 1/self.ctrl_edm, 1/self.ctrl_edm )
            elif  self.normalization_method == 'M': # do nothing
                pass
            elif self.normalization_method == 'N':
                self.normalize( 1/self.treat_edm, 1/self.ctrl_edm )
            else:
                raise NotImplemented
            self.normalization_method = 'M'                
        elif normalization_method == 'N':
            if self.normalization_method == 'T':
                self.normalize( self.treat_edm, self.treat_edm )
            elif self.normalization_method == 'C':
                self.normalize( self.ctrl_edm, self.ctrl_edm )
            elif  self.normalization_method == 'M':
                self.normalize( self.treat_edm, self.ctrl_edm )
            elif self.normalization_method == 'N': # do nothing
                pass
            else:
                raise NotImplemented
            self.normalization_method = 'N'            

    cdef normalize ( self, double treat_scale, double control_scale ):
        cdef:
            np.ndarray d
            long l, i
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                d[ i, 1 ] *= treat_scale
                d[ i, 2 ] *= control_scale                
        return

    cpdef change_score_method (self, char scoring_method):
        """
        scoring_method:  p: -log10 pvalue;
                         q: -log10 qvalue;
                         l: log10 likelihood ratio ( minus for depletion )
                         f: log10 fold enrichment
                         d: substraction
                         m: fragment pileup per million reads
        """
        if scoring_method == 'p':
            self.compute_pvalue()
        elif scoring_method == 'q':
            #if not already calculated p, compute pvalue first
            if self.scoring_method != 'p':
                self.compute_pvalue()
            self.compute_qvalue()
        elif scoring_method == 'l':
            self.compute_likelihood()
        elif scoring_method == 'f':
            self.compute_foldenrichment()
        elif scoring_method == 'd':
            self.compute_substraction()
        elif scoring_method == 'm':
            self.compute_SPMR()
        else:
            raise NotImplemented
            
    cdef compute_pvalue ( self ):
        """Compute -log_{10}(pvalue)
        """
        cdef:
            np.ndarray d
            long l, i
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                d[ i, 3 ] =  get_pscore( d[ i, 1] / 100, d[ i, 2] / 100.0 )
        self.scoring_method = 'p'
        return 

    cdef compute_qvalue ( self ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            dict pqtable
            long i,l,j,p
            str chrom
            
        # pvalue should be computed first!
        assert self.scoring_method == 'p'
        # make pqtable
        pqtable = self.make_pq_table()
        
        # convert p to q

        # convert pvalue2qvalue to a simple dict
        s_p2q = Int64HashTable()
        #g = pvalue2qvalue.get
        for i in pqtable.keys():
        #for i in range(pvalue2qvalue.shape[0]):
            s_p2q.set_item(i,pqtable[i][0])

        g = s_p2q.get_item
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                d[ i, 3 ] =  g( d[ i, 3 ])
        
        self.scoring_method = 'q'
        return

    cdef dict make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-100log10pvalue:(-100log10qvalue,rank,basepairs)} relationships.
        """
        cdef:
            long n, pre_p, this_p, length, j, pre_l, l, this_v, pre_v, v
            long N, k, q, pre_q
            double f
            str chrom
            np.ndarray d_chrom
            dict pvalue2qvalue

        assert self.scoring_method == 'p'
        
        n = self.total()
        value_dict = Int64HashTable()
        unique_values = pyarray(BYTE4,[])
        # this is a table of how many positions each p value occurs at
        for chrom in self.data.keys():
            # for each chromosome
            pre_p  = 0
            d_chrom = self.data[chrom]
            length = self.datalength[chrom]
            for j in xrange(length):
                this_p = d_chrom[ j, 0 ]
                this_v = d_chrom[ j, 3 ]
                assert this_v == this_v, "NaN at %d" % pos
                if value_dict.has_key(this_v):
                    value_dict.set_item(this_v, value_dict.get_item(this_v) + this_p - pre_p)
                else:
                    value_dict.set_item(this_v, this_p - pre_p)
                    unique_values.append(this_v)
                pre_p = this_p

        N = 0
        for i in xrange(len(unique_values)):
            N += value_dict.get_item(unique_values[i])
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        pvalue2qvalue = {pre_v:[0,k,0]}              # pvalue:[qvalue,rank,bp_with_this_pvalue]
        unique_values = sorted(unique_values,reverse=True)
        for i in xrange(len(unique_values)):
            v = unique_values[i]
            l = value_dict.get_item(v)
            q = v + int((log10(k) + f) * 100) # we save integers here.
            q = max(0,min(pre_q,q))           # make q-score monotonic
            pvalue2qvalue[v] = [q, k, 0]
            pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
            pre_v = v
            pre_q = q
            k+=l
        pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
        # pop the first -1e100 one
        pvalue2qvalue.pop(-2147483647)

        return pvalue2qvalue

    cdef compute_likelihood ( self ):
        cdef:
            np.ndarray d
            long l, i
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                d[ i, 3 ] =  logLR( d[ i, 1]/100.0, d[ i, 2]/100.0 )
        self.scoring_method = 'l'
        return 

    cdef compute_foldenrichment ( self ):
        cdef:
            np.ndarray d
            long l, i
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                #d[ i, 3] = get_foldenrichment ( 100.0 * d[ i, 1], d[ i, 2] ):
                # add pseudo count 1 = 100 in the #2 and #3 column
                d[ i, 3 ] =  int ( 100.0 * (d[ i, 1] + 100) / (d[ i, 2] + 100)  )
        self.scoring_method = 'f'
        return

    cdef compute_substraction ( self ):
        cdef:
            np.ndarray d
            long l, i
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                d[ i, 3 ] =  int ( d[ i, 1] - d[ i, 2] )
        self.scoring_method = 'd'
        return

    cdef compute_SPMR ( self ):
        cdef:
            np.ndarray d
            long l, i
            float scale
        if self.normalization_method == 'T' or self.normalization_method == 'N':
            scale = self.treat_edm
        elif self.normalization_method == 'C':
            scale = self.ctrl_edm
        elif self.normalization_method == 'M':
            scale = 1
        
        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                d[ i, 3 ] =  d[ i, 1] / scale # two digit precision may not be enough...
        self.scoring_method = 'm'
        return

    cpdef write_bedGraph ( self, fhd, str name, str description, short column = 3 ):
        """Write all data to fhd in bedGraph Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 1: chip, 2: control, 3: score

        """
        cdef:
            str chrom
            int l, pre, i, p 
            float pre_v, v
            float scale = 100.0
            np.ndarray d, pos, value

        assert column in range( 1, 4 ), "column should be between 1, 2 or 3."
        
        write = fhd.write

        if self.trackline:
            # this line is REQUIRED by the wiggle format for UCSC browser
            write( "track type=bedGraph name=\"%s\" description=\"%s\"\n" % ( name, description ) )
        
        chrs = self.get_chr_names()
        for chrom in chrs:
            d = self.data[ chrom ]
            pos = d[ :, 0 ]
            value = d[ :, column ]
            l = self.datalength[ chrom ]
            pre = 0
            if d.shape[ 0 ] == 0: continue # skip if there's no data
            pre_v = value[ 0 ] / scale
            for i in range( 1, l ):
                v = value[ i ] / scale
                p = pos[ i-1 ]
                if pre_v != v: 
                    write( "%s\t%d\t%d\t%.2f\n" % ( chrom, pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            write( "%s\t%d\t%d\t%.2f\n" % ( chrom, pre, p, pre_v ) )
            
        return True

    cpdef call_peaks (self, int cutoff=500, int min_length=200, int max_gap=50, bool call_summits=False):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 500. Note, the values stored in this class are actual value * 100, so 500 means 5.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        ptrack:  an optional track for pileup heights. If it's not None, use it to find summits. Otherwise, use self/scoreTrack.
        """
        cdef:
            int i
            str chrom
            np.ndarray pos, sample, control, value, above_cutoff, above_cutoff_v, above_cutoff_endpos, above_cutoff_startpos, above_cutoff_sv
            list peak_content
        
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        #if call_summits: close_peak = self.__close_peak2
        #else:
        # temporarily not use subpeak method
        #close_peak = self.__close_peak
        
        for chrom in chrs:
            peak_content = []           # to store points above cutoff

            pos = self.data[chrom][ :, 0 ]
            sample = self.data[chrom][ :, 1 ]
            control = self.data[chrom][ :, 2 ]            
            value = self.data[chrom][ :, 3 ]

            above_cutoff = np.nonzero( value >= cutoff )[0] # indices where score is above cutoff
            above_cutoff_v = value[above_cutoff] # scores where score is above cutoff

            above_cutoff_endpos = pos[above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = pos[above_cutoff-1] # start positions of regions where score is above cutoff
            above_cutoff_sv= sample[above_cutoff] # sample pileup height where score is above cutoff

            if above_cutoff_v.size == 0:
                # nothing above cutoff
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
                    self.__close_peak(peak_content, peaks, min_length, chrom, max_gap/2 )
                    peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]
            
            # save the last peak
            if not peak_content:
                continue
            else:
                self.__close_peak(peak_content, peaks, min_length, chrom, max_gap/2 )

        return peaks
       
    # def __close_peak2 (self, peak_content, peaks, min_length, chrom, colname, smoothlen=50):
    #     # this is where the summits are called, need to fix this
    #     end, start = peak_content[ -1 ][ 1 ], peak_content[ 0 ][ 0 ]
    #     if end - start < min_length: return # if the peak is too small, reject it
    #     #for (start,end,value,summitvalue,index) in peak_content:
    #     peakdata = np.zeros(end - start, dtype='float32')
    #     peakindices = np.zeros(end - start, dtype='int32')
    #     for (tmpstart,tmpend,tmpvalue,tmpsummitvalue, tmpindex) in peak_content:
    #         i, j = tmpstart-start, tmpend-start
    #         peakdata[i:j] = self.data[chrom]['sample'][tmpindex]
    #         peakindices[i:j] = tmpindex
    #     # apply smoothing window of tsize / 2
    #     w = np.ones(smoothlen, dtype='float32')
    #     smoothdata = np_convolve(w/w.sum(), peakdata, mode='same')
    #     # find maxima and minima
    #     local_extrema = np.where(np.diff(np.sign(np.diff(smoothdata))))[0]+1
    #     # get only maxima by requiring it be greater than the mean
    #     # might be better to take another derivative instead
    #     plateau_offsets = np.intersect1d(local_extrema,
    #                                      np.where(peakdata>peakdata.mean())[0])
    #     # sometimes peak summits are plateaus, so check for adjacent coordinates
    #     # and take the middle ones if needed
    #     if len(plateau_offsets)==0:
    #     #####################################################################
    #     # ***failsafe if no summits so far***                               #
    #         summit_offset_groups = [[(end - start) / 2]]                    #
    #     ##################################################################### 
    #     elif len(plateau_offsets) == 1:
    #         summit_offset_groups = [[plateau_offsets[0]]]
    #     else:
    #         previous_offset = plateau_offsets[0]
    #         summit_offset_groups = [[previous_offset]]
    #         for offset in plateau_offsets:
    #             if offset == previous_offset + 1:
    #                 summit_offset_groups[-1].append(offset)
    #             else:
    #                 summit_offset_groups.append([offset])
    #     summit_offsets = []
    #     for offset_group in summit_offset_groups:
    #         summit_offsets.append(offset_group[len(offset_group) / 2])
    #     summit_indices = peakindices[summit_offsets]
    #     # also purge offsets that have the same summit_index
    #     unique_offsets = []
    #     summit_offsets = np.fromiter(summit_offsets, dtype='int32')
    #     for index in np.unique(summit_indices):
    #         those_index_indices = np.where(summit_indices == index)[0]
    #         those_offsets = summit_offsets[those_index_indices]
    #         unique_offsets.append(int(those_offsets.mean()))
    #     # also require a valley of at least 0.6 * taller peak
    #     # in every adjacent two peaks or discard the lesser one
    #     # this behavior is like PeakSplitter
    #     better_offsets = []
    #     previous_offset = unique_offsets.pop()
    #     while True:
    #         if len(unique_offsets) == 0:
    #             better_offsets.append(previous_offset)
    #             break
    #         else:
    #             this_offset = unique_offsets.pop()
    #             this_h, prev_h = peakdata[[this_offset, previous_offset]]
    #             if this_h > prev_h:
    #                 prev_is_taller = False
    #                 min_valley = 0.6 * this_h
    #             else:
    #                 prev_is_taller = True
    #                 min_valley = 0.6 * prev_h
    #             s = slice(this_offset, previous_offset)
    #             valley = np.where(peakdata[s] < min_valley)[0]
    #             if len(valley) > 0: better_offsets.append(previous_offset)
    #             else:
    #                 if prev_is_taller: continue # discard this peak
    #                 # else: discard previous peak by ignoring it
    #             previous_offset = this_offset
    #     better_offsets.reverse()
    #     better_indices = peakindices[better_offsets]
    #     assert len(better_offsets) > 0, "Lost peak summit(s) near %s %d" % (chrom, start) 
    #     for summit_offset, summit_index in zip(better_offsets, better_indices):
    #         peaks.add( chrom,
    #                    start,
    #                    end,
    #                    summit      = start + summit_offset,
    #                    peak_score  = self.data[chrom][colname][ summit_index ],
    #                    pileup      = self.data[chrom]['sample'][ summit_index ], # should be the same as summit_value
    #                    pscore      = self.data[chrom]['-100logp'][ summit_index ]/100.0,
    #                    fold_change = self.data[chrom]['sample'][ summit_index ]/self.data[chrom]['control'][ summit_index ],
    #                    qscore      = self.data[chrom]['-100logq'][ summit_index ]/100.0,
    #                    )
    #     # start a new peak
    #     return True

    cdef __close_peak (self, list peak_content, peaks, int min_length, str chrom, int smoothlen=0):
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peaks: a PeakIO object

        """
        cdef:
            int summit_pos, tstart, tend, tmpindex, summit_index, i, midindex
            double summit_value, tvalue, tsummitvalue

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tsummit = []
            summit_pos   = 0
            summit_value = summit_pos
            for i in range(len(peak_content)):
                (tstart,tend,tvalue,tsummitvalue, tindex) = peak_content[i]
                #for (tstart,tend,tvalue,tsummitvalue, tindex) in peak_content:
                if not summit_value or summit_value < tsummitvalue:
                    tsummit = [(tend + tstart) / 2, ]
                    tsummit_index = [ tindex, ]
                    summit_value = tsummitvalue
                elif summit_value == tsummitvalue:
                    # remember continuous summit values
                    tsummit.append(int((tend + tstart) / 2))
                    tsummit_index.append( tindex )
            # the middle of all highest points in peak region is defined as summit
            midindex = int((len(tsummit) + 1) / 2) - 1
            summit_pos    = tsummit[ midindex ]
            summit_index  = tsummit_index[ midindex ]
            if self.scoring_method == 'q':
                qscore = self.data[chrom][ summit_index, 3 ] / 100.0
            else:
                # if q value is not computed, use -1
                qscore = -1

            peaks.add( chrom,
                       peak_content[0][0],
                       peak_content[-1][1],
                       summit      = summit_pos,
                       peak_score  = self.data[chrom][ summit_index, 3 ] / 100.0,
                       pileup      = self.data[chrom][ summit_index, 1 ] / 100.0, # should be the same as summit_value
                       pscore      = get_pscore(self.data[chrom][ summit_index, 1 ] / 100, self.data[chrom][ summit_index, 2 ] / 100.0 ) /100.0,
                       fold_change = float ( self.data[chrom][ summit_index, 1 ] + 100 ) / ( self.data[chrom][ summit_index, 2 ] + 100 ),
                       qscore      = qscore,
                       )
            # start a new peak
            return

    # def call_broadpeaks (self, int lvl1_cutoff=500, int lvl2_cutoff=100,
    #                      int min_length=200, int lvl1_max_gap=50, int lvl2_max_gap=400,
    #                      str colname='-100logq'):
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
    #     cdef:
    #         int i
    #         str chrom
        
    #     assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
    #     assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."        
    #     lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff, min_length=min_length, max_gap=lvl1_max_gap, colname=colname)
    #     lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff, min_length=min_length, max_gap=lvl2_max_gap, colname=colname)
    #     chrs = lvl1_peaks.peaks.keys()
    #     broadpeaks = BroadPeakIO()
    #     # use lvl2_peaks as linking regions between lvl1_peaks
    #     for chrom in chrs:
    #         lvl1peakschrom = lvl1_peaks.peaks[chrom]
    #         lvl2peakschrom = lvl2_peaks.peaks[chrom]
    #         lvl1peakschrom_next = iter(lvl1peakschrom).next
    #         tmppeakset = []             # to temporarily store lvl1 region inside a lvl2 region
    #         # our assumption is lvl1 regions should be included in lvl2 regions
    #         try:
    #             lvl1 = lvl1peakschrom_next()
    #         except StopIteration:
    #             break
    #         for i in range( len(lvl2peakschrom) ):
    #             # for each lvl2 peak, find all lvl1 peaks inside
    #             lvl2 = lvl2peakschrom[i]
    #             try:
    #                 while True:
    #                     if lvl2["start"] <= lvl1["start"]  and lvl1["end"] <= lvl2["end"]:
    #                         tmppeakset.append(lvl1)
    #                     else:
    #                         if tmppeakset:
    #                             self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
    #                         tmppeakset = []
    #                         break
    #                     lvl1 = lvl1peakschrom_next()
    #             except StopIteration:
    #                 if tmppeakset:
    #                     self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)                    
    #                 break
        
    #     return lvl1_peaks, broadpeaks


    # def __add_broadpeak (self, bpeaks, str chrom, lvl2peak, lvl1peakset):
    #     """Internal function to create broad peak.
        
    #     """
    #     cdef:
    #         int blockNum, thickStart, thickEnd, start, end
    #         str blockSizes, blockStarts
        
    #     start      = lvl2peak["start"]
    #     end        = lvl2peak["end"]
    #     thickStart = lvl1peakset[0]["start"]
    #     thickEnd   = lvl1peakset[-1]["end"]
    #     blockNum   = int(len(lvl1peakset))
    #     blockSizes = ",".join( map(lambda x:str(x["length"]),lvl1peakset) )
    #     blockStarts = ",".join( map(lambda x:str(x["start"]-start),lvl1peakset) )
    #     if lvl2peak["start"] != thickStart:
    #         # add 1bp mark for the start of lvl2 peak
    #         blockNum += 1
    #         blockSizes = "1,"+blockSizes
    #         blockStarts = "0,"+blockStarts
    #     if lvl2peak["end"] != thickEnd:
    #         # add 1bp mark for the end of lvl2 peak            
    #         blockNum += 1
    #         blockSizes = blockSizes+",1"
    #         blockStarts = blockStarts+","+str(end-start-1)
        
    #     bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=thickStart, thickEnd=thickEnd,
    #                blockNum = blockNum, blockSizes = blockSizes, blockStarts = blockStarts)
    #     return bpeaks

    cdef long total ( self ):
        """Return the number of regions in this object.

        """
        cdef:
            long t
            str chrom
        
        t = 0
        for chrom in self.data.keys():
            t += self.datalength[chrom]
        return t

