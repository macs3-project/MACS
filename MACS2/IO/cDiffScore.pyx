# Time-stamp: <2013-09-11 22:52:23 Tao Liu>

"""Module for Feature IO classes.

Copyright (c) 2013 Tao Liu <vladimir.liu@gmail.com>, Ben Schiller <>

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

from cpython cimport bool

from scipy.stats import chi2

from operator import itemgetter

from cython.parallel import parallel, prange

cimport cython

from libc.math cimport log10, log, floor, ceil

from MACS2.hashtable import Int64HashTable, Float64HashTable
from MACS2.cProb cimport poisson_cdf


import logging

# ------------------------------------
# constants
# ------------------------------------
__version__ = "DiffScoreTrackI $Revision$"
__author__ = "Tao Liu <vladimir.liu@gmail.com>"
__doc__ = "DiffScoreTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------

def do_nothing(*args, **kwargs):
    pass

LOG10_E = 0.43429448190325176

pscore_khashtable = Int64HashTable()

cdef inline double get_pscore ( int observed, double expectation ):
    """Get p-value score from Poisson test. First check existing
    table, if failed, call poisson_cdf function, then store the result
    in table.
    
    """
    cdef:
        double score
        long key_value
    
    #key_value = ( observed, expectation )
    key_value = hash( (observed, expectation ) )
    try:
        return pscore_khashtable.get_item(key_value)
    except KeyError:
        score = -1*poisson_cdf(observed,expectation,False,True)
        pscore_khashtable.set_item(key_value, score)
        return score

cdef double get_interpolated_pscore ( double observed, double expectation ):
    cdef:
        double pscore
        double observed_floor, observed_ceil, step
        double floor_pscore, ceil_pscore
    observed_floor = floor(observed)
    observed_ceil = ceil(observed)
    step = observed - observed_floor
    floor_pscore = get_pscore( int(observed_floor), expectation)
    ceil_pscore = get_pscore( int(observed_ceil), expectation)
    pscore = (ceil_pscore - floor_pscore) * (observed - observed_floor) + \
             floor_pscore
    return pscore

diff_logLR_khashtable = Int64HashTable()

cdef inline double logLR_4diff (double x, double y):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ).
    
    * Always positive, used for evaluating difference only.*
    
    """
    cdef:
        double s
        long key_value
    
    key_value = hash( (x, y) )
    try:
        return diff_logLR_khashtable.get_item( key_value )
    except KeyError:
        if y > x: y, x = x, y
        if x==y: s = 0
        
        else: s = (x*(log(x)-log(y))+y-x)*LOG10_E
        diff_logLR_khashtable.set_item(key_value, s)
        return s

# ------------------------------------
# Classes
# ------------------------------------

cdef class DiffScoreTrackI:
    cdef:
        dict pos, t1, c1, t2, c2, tvsc1, tvsc2, t1vs2, tlogLR
        dict where_peaks, diff_peaks, diff_qvalues
        dict which_peaks1, which_peaks2
        object p1io, p2io
        dict datalength                 # length of data array of each chromosome
        float cond1_depth               # seq depth in million of treatment
        float cond2_depth               # seq depth in million of control
        int pseudocount                 # the pseudocount used to calcuate LLR
        float cutoff
        dict pvalue_stat1, pvalue_stat2
        str track_scoring_method
        str diff_scoring_method
    
    def __init__ (self, t1bdg, c1bdg, t2bdg, c2bdg, float cond1_depth = 1.0, float cond2_depth = 1.0, int pseudocount = 1 ):
        self.pos = {}
        self.t1 = {}
        self.c1 = {}
        self.tvsc1 = {}
        self.t2 = {}
        self.c2 = {}
        self.tvsc2 = {}
        self.t1vs2 = {}
        self.diff_qvalues = {}
        self.tlogLR = {}
        self.where_peaks = {}
        self.which_peaks1 = {}
        self.which_peaks2 = {}
        self.diff_peaks = {}
        self.datalength = {}
        self.cond1_depth = cond1_depth
        self.cond2_depth = cond2_depth
        self.pseudocount = pseudocount
        self.pvalue_stat1 = {}
        self.pvalue_stat2 = {}
        t1chrs = t1bdg.get_chr_names()
        c1chrs = c1bdg.get_chr_names()
        t2chrs = t2bdg.get_chr_names()
        c2chrs = c2bdg.get_chr_names()      
        common_chrs = reduce(lambda x,y:x.intersection(y), (t1chrs,c1chrs,t2chrs,c2chrs))

        for chrname in common_chrs:
            (cond1_treat_ps, cond1_treat_vs) = t1bdg.get_data_by_chr(chrname)
            (cond1_control_ps, cond1_control_vs) = c1bdg.get_data_by_chr(chrname)
            (cond2_treat_ps, cond2_treat_vs) = t2bdg.get_data_by_chr(chrname)
            (cond2_control_ps, cond2_control_vs) = c2bdg.get_data_by_chr(chrname)
            chrom_max_len = len(cond1_treat_ps) + len(cond1_control_ps) +\
                            len(cond2_treat_ps) + len(cond2_control_ps)
            self.add_chromosome( chrname, chrom_max_len )
            self.build_chromosome( chrname,
                                   cond1_treat_ps, cond1_control_ps,
                                   cond2_treat_ps, cond2_control_ps,
                                   cond1_treat_vs, cond1_control_vs,
                                   cond2_treat_vs, cond2_control_vs )
            
    cpdef set_pseudocount( self, int pseudocount ):
        self.pseudocount = pseudocount

    cdef build_chromosome( self, chrname,
                           cond1_treat_ps, cond1_control_ps,
                           cond2_treat_ps, cond2_control_ps,
                           cond1_treat_vs, cond1_control_vs,
                           cond2_treat_vs, cond2_control_vs ):
                                           

        c1tpn = iter(cond1_treat_ps).next
        c1cpn = iter(cond1_control_ps).next
        c2tpn = iter(cond2_treat_ps).next
        c2cpn = iter(cond2_control_ps).next
        c1tvn = iter(cond1_treat_vs).next
        c1cvn = iter(cond1_control_vs).next
        c2tvn = iter(cond2_treat_vs).next
        c2cvn = iter(cond2_control_vs).next

        pre_p = 0

        try:
            c1tp = c1tpn()
            c1tv = c1tvn()
            
            c1cp = c1cpn()
            c1cv = c1cvn()

            c2tp = c2tpn()
            c2tv = c2tvn()
            
            c2cp = c2cpn()
            c2cv = c2cvn()            

            while True:
                minp = min(c1tp, c1cp, c2tp, c2cp)
                self.add( chrname, pre_p, c1tv, c1cv, c2tv, c2cv )
                pre_p = minp
                if c1tp == minp:
                    c1tp = c1tpn()
                    c1tv = c1tvn()
                if c1cp == minp:
                    c1cp = c1cpn()
                    c1cv = c1cvn()
                if c2tp == minp:
                    c2tp = c2tpn()
                    c2tv = c2tvn()
                if c2cp == minp:
                    c2cp = c2cpn()
                    c2cv = c2cvn()                    
        except StopIteration:
            # meet the end of either bedGraphTrackI, simply exit
            pass
    
    # REQUIRES Cython >= 0.16
    @cython.boundscheck(False)
    cpdef rebuild_chromosomes( self ):
        cdef:
#            np.ndarray[np.int32_t] pos, pos_copy
#            np.ndarray[np.float32_t] t1, c1, t2, c2
#            np.ndarray[np.float32_t] t1_copy, c1_copy, t2_copy, c2_copy
            int[:] pos, pos_copy
            float[:] t1, c1, t2, c2
            float[:] t1_copy, c1_copy, t2_copy, c2_copy
            int i, j, i_new, n, k
            str chrom
            np.ndarray[np.int32_t] peaks
            int datalength
        for chrom in self.pos.keys():
            try:
                peaks1 = self.p1io.get_data_from_chrom(chrom)
            except KeyError: peaks1 = []
            try:
                peaks2 = self.p2io.get_data_from_chrom(chrom)
            except KeyError: peaks2 = []
            peaks = np.unique(np.array(map(itemgetter('start'), peaks1) +
                             map(itemgetter('end'), peaks1) +
                             map(itemgetter('start'), peaks2) +
                             map(itemgetter('end'), peaks2),
                             dtype='int32'))
            n = peaks.size
            datalength = self.datalength[chrom]
            # resize the originals
            self.pos[chrom].resize(datalength + n, refcheck=False)
            self.t1[chrom].resize(datalength + n, refcheck=False)
            self.c1[chrom].resize(datalength + n, refcheck=False)
            self.t2[chrom].resize(datalength + n, refcheck=False)
            self.c2[chrom].resize(datalength + n, refcheck=False)
            pos = self.pos[chrom]
            t1 = self.t1[chrom]
            c1 = self.c1[chrom]
            t2 = self.t2[chrom]
            c2 = self.c2[chrom]
            pos_copy = np.ndarray(datalength, 'int32')
            t1_copy = np.ndarray(datalength, 'float32')
            c1_copy = np.ndarray(datalength, 'float32')
            t2_copy = np.ndarray(datalength, 'float32')
            c2_copy = np.ndarray(datalength, 'float32')
            # GIL not needed here, provides speedup
            with cython.boundscheck(False):
                with nogil, parallel():
                    for i in prange(datalength):
#            for i in range(datalength):
                        pos_copy[i] = pos[i]
                        t1_copy[i] = t1[i]
                        c1_copy[i] = c1[i]
                        t2_copy[i] = t2[i]
                        c2_copy[i] = c2[i]
            
                j = 0
                i_new = 0
                for i in range(datalength):
                    while True:
                        if j == n: break
                        if peaks[j] < pos_copy[i]:
                            pos[i_new] = peaks[j]
                            j += 1
                            t1[i_new] = t1_copy[i]
                            t2[i_new] = t2_copy[i]
                            c1[i_new] = c1_copy[i]
                            c2[i_new] = c2_copy[i]
                            i_new += 1
                        elif peaks[j] == pos[i]:
                            j += 1
                            break
                        else: # pos[i] < peaks[j], keep increasing pos
                            break
                    # this would just copy the array without the above loop
                    if j == n: break
                    pos[i_new] = pos_copy[i]
                    t1[i_new] = t1_copy[i]
                    t2[i_new] = t2_copy[i]
                    c1[i_new] = c1_copy[i]
                    c2[i_new] = c2_copy[i]
                    i_new += 1
#            pos_copy.resize(0, refcheck=False)
#            t1_copy.resize(0, refcheck=False)
#            t2_copy.resize(0, refcheck=False)
#            c1_copy.resize(0, refcheck=False)
#            c1_copy.resize(0, refcheck=False)
            # j should never be greater than n - 1 if we used MACS2
            self.pos[chrom].resize(i_new, refcheck=False)
            self.t1[chrom].resize(i_new, refcheck=False)
            self.c1[chrom].resize(i_new, refcheck=False)
            self.t2[chrom].resize(i_new, refcheck=False)
            self.c2[chrom].resize(i_new, refcheck=False)
            self.datalength[chrom] = i_new # resize everything last
        self.finalize()
        
    cdef add_chromosome ( self, str chrom, int chrom_max_len ):
        """
        chrom: chromosome name
        chrom_max_len: maximum number of data points in this chromosome
        
        """
        if not self.pos.has_key(chrom):
            self.pos[chrom] = np.zeros( chrom_max_len, dtype="int32" ) # pos
            self.t1[chrom] = np.zeros( chrom_max_len, dtype="float32" ) # t1
            self.c1[chrom] = np.zeros( chrom_max_len, dtype="float32" ) # c1
            self.tvsc1[chrom] = np.zeros( chrom_max_len, dtype="float32" ) # c1
            self.t2[chrom] = np.zeros( chrom_max_len, dtype="float32" ) # t2
            self.c2[chrom] = np.zeros( chrom_max_len, dtype="float32" ) # c2
            self.tvsc2[chrom] = np.zeros( chrom_max_len, dtype="float32" ) # c1
            self.t1vs2[chrom] = np.zeros(chrom_max_len, dtype="float32" ) # c1
            self.tlogLR[chrom] = np.ndarray( chrom_max_len, dtype="float32")
            self.datalength[chrom] = 0

    cdef add (self, str chromosome, int endpos, float t1, float c1, float t2, float c2):
        """Add a chr-endpos-sample-control block into data
        dictionary.

        chromosome: chromosome name in string
        endpos    : end position of each interval in integer
        chip      : ChIP pileup value of each interval in float
        control   : Control pileup value of each interval in float

        *Warning* Need to add regions continuously.
        """
        cdef int i
        i = self.datalength[chromosome]
        self.pos[chromosome][ i ] = endpos
        self.t1[chromosome][ i ] = t1 * self.cond1_depth
        self.c1[chromosome][ i ] = c1 * self.cond1_depth
        self.t2[chromosome][ i ] = t2 * self.cond2_depth
        self.c2[chromosome][ i ] = c2 * self.cond2_depth
        self.datalength[chromosome] += 1

    cpdef finalize ( self ):
        """
        Adjust array size of each chromosome.

        """
        cdef:
            str chrom, k
            int l

        for chrom in self.pos.keys():
            l = self.datalength[chrom]
            self.pos[chrom].resize( l, refcheck = False )
            self.t1[chrom].resize( l, refcheck = False )
            self.c1[chrom].resize( l, refcheck = False )
            self.tvsc1[chrom].resize( l, refcheck = False )
            self.t2[chrom].resize( l, refcheck = False )
            self.c2[chrom].resize( l, refcheck = False )
            self.tvsc2[chrom].resize( l, refcheck = False )
            self.t1vs2[chrom].resize( l, refcheck = False )
            self.tlogLR[chrom].resize( l, refcheck = False )
        return
    
    cpdef set_track_score_method (self, str scoring_method):
        """
        scoring_method:  p: -log10 pvalue;
                         q: -log10 qvalue;
                         l: log10 likelihood ratio ( minus for depletion )
                         f: log10 fold enrichment
                         F: linear fold enrichment
                         d: subtraction
                         m: fragment pileup per million reads
        """
        if scoring_method == 'p':
            self.compute_treatcontrol_pvalues()
        elif scoring_method == 'q':
            #if not already calculated p, compute pvalue first
            if self.track_scoring_method != 'p':
                self.compute_treatcontrol_pvalues()
            self.compute_treatcontrol_qvalues()
        else:
            raise NotImplemented
        
    cdef compute_treatcontrol_pvalues ( self ):
        """Compute -log_{10}(pvalue)
        """
        cdef:
            np.ndarray[np.float32_t] p, c, v
            np.ndarray[np.int32_t] pos
            long l, i, prev_pos
            str chrom
        
        pseudocount = self.pseudocount
        for chrom in self.pos.keys():
            prev_pos = 0
            pos = self.pos[chrom]
            p = self.t1[chrom]
            c = self.c1[chrom]
            v = self.tvsc1[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                if c[i] > p[i]:
                    v[i] = 1
                else:
                    v[ i ] =  get_pscore(p[i] + pseudocount, c[i])
                try:
                    self.pvalue_stat1[v[ i ]] += pos[ i ] - prev_pos
                except:
                    self.pvalue_stat1[v[ i ]] = pos[ i ] - prev_pos
                prev_pos = pos[ i ]
                    
        for chrom in self.pos.keys():
            prev_pos = 0
            pos = self.pos[chrom]
            p = self.t2[chrom]
            c = self.c2[chrom]
            v = self.tvsc2[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                if c[i] > p[i]:
                    v[i] = 1
                else:
                    v[ i ] =  get_interpolated_pscore(p[i], c[i])
                try:
                    self.pvalue_stat2[v[ i ]] += pos[ i ] - prev_pos
                except:
                    self.pvalue_stat2[v[ i ]] = pos[ i ] - prev_pos
                prev_pos = pos[ i ]
        
        self.track_scoring_method = 'p'
        return 

    cdef compute_treatcontrol_qvalues ( self ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            dict pqtable1, pqtable2
            long i,l,j
            double k
            str chrom
            np.ndarray p, c, v
            
        # pvalue should be computed first!
        assert self.track_scoring_method == 'p'
        # make pqtable
        (pqtable1, pqtable2) = self.make_treatcontrol_pq_tables()
        
        # convert p to q

        # convert pvalue2qvalue to a simple dict based on khash
        # khash has big advantage while checking keys for millions of times.
        s_p2q1 = Float64HashTable()
        for k in pqtable1.keys():
            s_p2q1.set_item(k,pqtable1[k])

        s_p2q2 = Float64HashTable()
        for k in pqtable2.keys():
            s_p2q2.set_item(k,pqtable2[k])

        g1 = s_p2q1.get_item
        g2 = s_p2q2.get_item
        
        for chrom in self.pos.keys():
            v1 = self.tvsc1[chrom]
            v2 = self.tvsc2[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                v1[ i ] =  g1( v1[ i ])
                v2[ i ] =  g2( v2[ i ])
        
        self.track_scoring_method = 'q'
        return
    
    cpdef break_on_peaks(self, p1io, p2io):
        """Introduce breaks at peak regions as needed
        """
    
    cpdef call_peaks (self, float cutoff=2.0, int min_length=200, int max_gap=50):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 2. For -log10pvalue, it means 10^-5.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        """
        cdef:
            int i, first_i, length
            int first_start, this_start, this_end, last_end
            str chrom
            np.ndarray pos, sample, control, value, above_cutoff, above_cutoff_v, above_cutoff_endpos, above_cutoff_startpos, above_cutoff_sv
            np.ndarray in_peaks
        
        self.cutoff = cutoff
        for chrom in self.pos.keys():
            in_peaks = np.zeros(self.datalength[chrom], dtype=np.bool)

            pos = self.pos[chrom]
            value = self.tvsc1[chrom]
            above_cutoff = np.nonzero( value >= cutoff )[0] # indices where score is above cutoff
            above_cutoff_endpos = pos[above_cutoff] # end positions of regions where score is above cutoff

#            print "Regions > cutoff: %d" % above_cutoff.size
            if above_cutoff.size > 1:
                # Do zero manually
                first_i = 0
                this_start = pos[above_cutoff[0] - 1]
                this_end = above_cutoff_endpos[0]
                if this_start > this_end:
                    this_start = 0
                first_start = this_start
                last_end = this_end
                for i in range(1, above_cutoff.size):
                    this_start = above_cutoff_endpos[i - 1]
                    this_end = above_cutoff_endpos[i]
                    if first_i == -1:
                        first_i = i
                        first_start = this_start
                    elif (this_end - last_end) > max_gap:
                        if (last_end - first_start) >= min_length:
                            in_peaks[above_cutoff[first_i]:above_cutoff[i - 1]] = True
#                        else:
#                            print "Rejected", pos[above_cutoff[first_i]-1], pos[above_cutoff[i - 1]]
                        first_i = -1
                    last_end = this_end
            
                if not first_i == -1:
                    if last_end - first_start >= min_length:
                        in_peaks[above_cutoff[first_i]:above_cutoff[i]] = True
            
            value = self.tvsc2[chrom]
            above_cutoff = np.nonzero( value >= cutoff )[0] # indices where score is above cutoff
            above_cutoff_endpos = pos[above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = pos[above_cutoff-1] # start positions of regions where score is above cutoff

            if above_cutoff.size > 1:
                # Do zero manually
                first_i = 0
                this_start = pos[above_cutoff[0] - 1]
                this_end = above_cutoff_endpos[0]
                if this_start > this_end:
                    this_start = 0
                first_start = this_start
                last_end = this_end
                for i in range(1, above_cutoff.size):
                    this_start = above_cutoff_endpos[i - 1]
                    this_end = above_cutoff_endpos[i]
                    if first_i == -1:
                        first_i = i
                        first_start = this_start
                    elif (this_end - last_end) > max_gap:
                        if (last_end - first_start) >= min_length:
                            in_peaks[above_cutoff[first_i]:above_cutoff[i - 1]] = True
                        first_i = -1
                    last_end = this_end
            
                if not first_i == -1:
                    if last_end - first_start >= min_length:
                        in_peaks[above_cutoff[first_i]:above_cutoff[i]] = True

            self.where_peaks[chrom] = np.where(in_peaks)[0].astype('int32')
            print "Total peakage in bp", in_peaks.sum()
        
        return

    cpdef store_peaks(self, p1io, p2io):
        self.p1io = p1io
        self.p2io = p2io
        
    cpdef annotate_peaks(self):
        cdef:
            str chrom
            int i, i_max, ii, j, j_max
            np.ndarray[np.int32_t] pos, where_peaks, which_peaks1, which_peaks2
            np.ndarray[np.int32_t] p1starts, p1ends, p2starts, p2ends
        for chrom in self.pos.keys():
            pos = self.pos[chrom]
            i_max = self.datalength[chrom]
            
            which_peaks1 = -np.ones(self.datalength[chrom], dtype='int32')
            which_peaks2 = -np.ones(self.datalength[chrom], dtype='int32')
            try:
                data = self.p1io.get_data_from_chrom(chrom)
            except KeyError: data = []
            p1starts = np.array(map(itemgetter("start"), data), 'int32')
            p1ends = np.array(map(itemgetter("end"), data), 'int32')
            j = 0
            j_max = p1starts.size
            for i in range(i_max - 1):
                if j == j_max: break
                if pos[i] == p1starts[j]:
                # then the i + 1 fragment starts with the correct value
                    # find the end
                    end = p1ends[j]
                    for ii in range(i + 1, i_max):
                        if pos[ii] == end:
                            break
                        # assert pos[ii] < end, "something went wrong"
                    # ii end = true end, but ii is still part of peak
                    which_peaks1[(i + 1):(ii + 1)] = j
                    j += 1
                    # skip additional subpeaks
                    while True:
                        if j == j_max: break
                        if p1ends[j] == p1ends[j - 1]: j += 1
                        else: break
                    
            try:
                data = self.p2io.get_data_from_chrom(chrom)
            except KeyError: data = []
            p2starts = np.array(map(itemgetter("start"), data), 'int32')
            p2ends = np.array(map(itemgetter("end"), data), 'int32')
            j = 0
            j_max = p2starts.size
            for i in range(i_max - 1):
                if j == j_max: break
                if pos[i] == p2starts[j]:
                # then the i + 1 fragment starts with the correct value
                    # find the end
                    end = p2ends[j]
                    for ii in range(i + 1, i_max - 1):
                        if pos[ii] == end:
                            break
                        # assert pos[ii] < end, "something went wrong"
                    # ii end = true end, but ii is still part of peak
                    which_peaks2[(i + 1):(ii + 1)] = j
                    j += 1
                    # skip additional subpeaks
                    while True:
                        if j == j_max: break
                        if p2ends[j] == p2ends[j - 1]: j += 1
                        else: break
            
            where_peaks = np.where(np.logical_or(which_peaks1 >= 0,
                                                 which_peaks2 >= 0))[0].astype('int32')
            self.where_peaks[chrom] = where_peaks
            self.which_peaks1[chrom] = which_peaks1[where_peaks]
            self.which_peaks2[chrom] = which_peaks2[where_peaks]
            # note that we skipped peaks which have the same start end
            # we'll find them again later using groupby
        

    cpdef tuple make_treatcontrol_pq_tables ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-log10pvalue:(-log10qvalue,rank,basepairs)} relationships.
        """
        cdef:
            long pre_l, l, i, N, k
            double this_v, pre_v, v, q, pre_q
            double f
            str chrom
            np.ndarray v_chrom, pos_chrom
            dict pvalue2qvalue1, pvalue2qvalue2
            dict value_dict
            list unique_values

        assert self.track_scoring_method == 'p'

        value_dict = self.pvalue_stat1
        N = sum(value_dict.values())
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        pvalue2qvalue1 = {}#Float64HashTable()
        unique_values = sorted(value_dict.keys(), reverse=True) #sorted(unique_values,reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = value_dict[v]
            q = v + (log10(k) + f)
            q = max(0,min(pre_q,q))           # make q-score monotonic
            pvalue2qvalue1[ v ] = q
            pre_v = v
            pre_q = q
            k+=l
        
        value_dict = self.pvalue_stat2
        N = sum(value_dict.values())
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        pvalue2qvalue2 = {}#Float64HashTable()
        unique_values = sorted(value_dict.keys(), reverse=True) #sorted(unique_values,reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = value_dict[v]
            q = v + (log10(k) + f)
            q = max(0,min(pre_q,q))           # make q-score monotonic
            pvalue2qvalue2[ v ] = q
            pre_v = v
            pre_q = q
            k+=l
            
        return (pvalue2qvalue1, pvalue2qvalue2)

    cpdef compute_diff_pvalues( self ):
        """Compute -log_{10}(pvalue)
        """
        cdef:
            str chrom
            int i, chrom_length
            int pseudocount
            np.ndarray[np.float32_t] t1, t2, v
        rv = chi2(1) # a chi-squared distribution with one degree of freedom
        sf = rv.sf
        log10 = np.log10
        pseudocount = self.pseudocount
        for chrom in self.pos.keys():
            t1 = self.t1[chrom]
            c1 = self.c1[chrom]
            t2 = self.t2[chrom]
            c2 = self.c2[chrom]
            v = self.tlogLR[chrom]
            with cython.boundscheck(False):
                for i in range(self.pos[chrom].size): 
                    v[i] = logLR_4diff( t1[i] + pseudocount,
                                        t2[i] + pseudocount )
                self.t1vs2[chrom] = -log10(sf(2 * v / LOG10_E)).astype('float32')
    
    cpdef compute_diff_qvalues ( self ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            dict value_dict, pqtable
            list unique_values
            long i,l,j,prev_i
            np.ndarray p, c, v, data
            long pre_l
            double pre_v, unique_v, q, pre_q
            long N, k
            double f, pvalue
        
        # make pqtable for category
        value_dict = {}
        for chrom in self.t1vs2.keys():
            if self.t1vs2[chrom].size == 0: continue

            pos = self.pos[chrom]
            stat = self.t1vs2[chrom]
            where_peaks = self.where_peaks[chrom]
            prev_i = -1
            for j in range(where_peaks.size):
                i = where_peaks[j]
                try:
#                    if (prev_i + 1) == i:
#                        if pos[i]-pos[prev_i] > 1000: print pos[i], pos[i]-pos[prev_i]
#                    else:
#                        if pos[i]-pos[i-1] > 1000: print pos[i], pos[i]-pos[i-1], "i-"
                    value_dict[stat[i]] += pos[i] - pos[i - 1]
                except IndexError:
                    if not value_dict.has_key(stat[i]):
                        value_dict[stat[i]] = 0
                except KeyError:
                    value_dict[stat[i]] = pos[i] - pos[i - 1]
                prev_i = i

        N = sum(value_dict.values())
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        pqtable = {}#Float64HashTable()
        unique_values = sorted(value_dict.keys(), reverse=True) #sorted(unique_values,reverse=True)
        for i in range(len(unique_values)):
            unique_v = unique_values[i]
            l = value_dict[unique_v]
            q = unique_v + (log10(k) + f)
            q = max(0,min(pre_q,q))           # make q-score monotonic
            pqtable[ unique_v ] = q
            pre_v = unique_v
            pre_q = q
            k+=l
        
        # convert p to q

        # convert pvalue2qvalue to a simple dict based on khash
        # khash has big advantage while checking keys for millions of times.
        s_p2q = Float64HashTable()
        for pvalue in pqtable.keys():
            s_p2q.set_item(pvalue,pqtable[pvalue])

        g = s_p2q.get_item
        
        for chrom in self.t1vs2.keys():
            v = self.t1vs2[chrom][self.where_peaks[chrom]]
            qvalues = v.copy()
            for i in range(v.size):
                qvalues[ i ] =  g(v[ i ])
            self.diff_qvalues[chrom] = qvalues
        
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
        l = set(self.pos.keys())
        return l

    cpdef write_bedGraph ( self, fhd, str name, str description, short column = 3):
        """Write all data to fhd in bedGraph Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 1: chip, 2: control, 3: score

        """
        cdef:
            str chrom
            int l, pre, i, p 
            float pre_v, v
            np.ndarray pos, value

        assert column in range( 1, 4 ), "column should be between 1, 2 or 3."
        
        write = fhd.write

        if self.trackline:
            # this line is REQUIRED by the wiggle format for UCSC browser
            write( "track type=bedGraph name=\"%s\" description=\"%s\"\n" % ( name, description ) )
        
        chrs = self.get_chr_names()
        for chrom in chrs:
            pos = self.data[ chrom ][ 0 ]
            value = self.data[ chrom ][ column ]
            l = self.datalength[ chrom ]
            pre = 0
            if pos.shape[ 0 ] == 0: continue # skip if there's no data
            pre_v = value[ 0 ]
            for i in range( 1, l ):
                v = value[ i ]
                p = pos[ i-1 ]
                if abs(pre_v - v)>=1e-6: 
                    write( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            write( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre, p, pre_v ) )
            
        return True

#    cpdef write_matrix ( self, fhd, str name, str description ):
#        """Write all data to fhd into five columns Format:
#
#        col1: chr_start_end
#        col2: t1 vs c1
#        col3: t2 vs c2
#        col4: t1 vs t2
#        col5: t2 vs t1
#
#        fhd: a filehandler to save the matrix.
#
#        """
#        cdef:
#            str chrom
#            int l, pre, i, p 
#            float v1, v2, v3, v4
#            np.ndarray pos, value
#
#        write = fhd.write
#
#        chrs = self.get_chr_names()
#        for chrom in chrs:
#            pos = self.data[ chrom ][ 0 ]
#            value = self.data[ chrom ][ column ]
#            l = self.datalength[ chrom ]
#            pre = 0
#            if pos.shape[ 0 ] == 0: continue # skip if there's no data
#            for i in range( 0, l ):
#                v1 = self.data[ i ][ 1 ]
#                v2 = self.data[ i ][ 2 ]
#                v3 = self.data[ i ][ 3 ]
#                v4 = self.data[ i ][ 4 ]                
#                p = pos[ i ]
#                write( "%s:%d_%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % ( chrom, pre, p, v1, v2, v3, v4 ) )
#                pre = p
#            
#        return True

    @cython.boundscheck(False)
    cpdef call_diff_peaks (self, float cutoff=2.0, int min_length=50,
                           str score_method='q'):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 3. For log10 LR, it means 1000 or -1000.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        ptrack:  an optional track for pileup heights. If it's not None, use it to find summits. Otherwise, use self/scoreTrack.
        """
        cdef:
            int i, first_i, length
            str chrom
            np.ndarray[np.int32_t] pos, above_cutoff, qpos
            np.ndarray[np.float32_t] t1vst2, diff_qvalues
            np.ndarray[np.int32_t, ndim=2] diff_peaks
            int max_gap
            int n_peaks, bigN
        max_gap = min_length / 4
        
        chrs  = self.get_chr_names()

        self.cutoff = cutoff

        # qvalue conversion for each category to improve statistical power
        self.compute_diff_pvalues()
        self.compute_diff_qvalues()
        self.diff_scoring_method = score_method
        if not (score_method == 'q' or score_method == 'p'):
            raise NotImplementedError

        bigN = 1
        for chrom in sorted(chrs):
            n_diff_peaks = 0
            pos = self.pos[chrom]
            qpos = self.where_peaks[chrom]
            t1vst2 = self.t1vs2[chrom]

            if score_method == 'q':
                diff_qvalues = self.diff_qvalues[chrom]
            elif score_method == 'p':
                diff_qvalues = self.t1vs2[chrom][self.where_peaks[chrom]]
            else:
                raise NotImplementedError

            above_cutoff = qpos[np.nonzero(diff_qvalues >= cutoff)[0]]
            # we're going to recalculate this a few times, hopefully it's fast

            # peaks are stored as start_i, end_i (0-based, genomic half-open)
            diff_peaks = np.ndarray((above_cutoff.size, 3), dtype='int32')
            if above_cutoff.size > 1:
                # Do zero manually
                first_i = 0
                this_start = pos[above_cutoff[0] - 1]
                this_end = pos[above_cutoff[0]]
                if this_start > this_end:
                    this_start = 0
                first_start = this_start
                last_end = this_end
#                print "%d (%d), %d (%d)" %(this_end, i, first_start, first_i)
                for i in range(1, above_cutoff.size):
                    this_start = pos[above_cutoff[i] - 1]
                    this_end = pos[above_cutoff[i]]
                    if first_i == -1:
                        first_i = i
                        first_start = this_start
                    elif (this_end - last_end) > max_gap:
                        if (last_end - first_start) >= min_length:
                            diff_peaks[n_diff_peaks,0] = first_i
                            diff_peaks[n_diff_peaks,1] = i - 1
                            diff_peaks[n_diff_peaks,2] = bigN + n_diff_peaks
                            n_diff_peaks += 1
#                        else:
#                            print "Rejectedd", pos[above_cutoff[first_i]-1], pos[above_cutoff[i - 1]]
                        first_i = -1
                    last_end = this_end
            
                if not first_i == -1:
                    if last_end - first_start >= min_length:
                        diff_peaks[n_diff_peaks,0] = first_i
                        diff_peaks[n_diff_peaks,1] = i
                        diff_peaks[n_diff_peaks,2] = bigN + n_diff_peaks
                        n_diff_peaks += 1

            bigN += n_diff_peaks
            diff_peaks.resize((n_diff_peaks, 3), refcheck=False)
#            print n_diff_peaks, "diff peaks"
            self.diff_peaks[chrom] = diff_peaks
        return
    
    cpdef print_some_peaks(self, int max = 10):
        """for testing only
        """
        cdef:
            int start, end, i
            str chrom
        for chrom in sorted(self.diff_peaks.keys()):
            i = 0
            for end in self.where_peaks[chrom]:
                i += 1
                if i == max: break
                print '%s\t%d (%d)' % (chrom, self.pos[chrom][end], end)
#            i = self.pos[chrom].searchsorted(49551000)
#            j = self.pos[chrom].searchsorted(49553000)
#            print i, j
#            for x in range(i,j):
#                print self.pos[chrom][x], self.t1[chrom][x], self.c1[chrom][x], self.tvsc1[chrom][x], \
#                self.t2[chrom][x], self.c2[chrom][x], self.tvsc2[chrom][x]
        return
    
    cpdef print_diff_peaks(self):
        """ for testing only
        """
        cdef:
            int start, end
            str chrom
        for chrom in sorted(self.diff_peaks.keys()):
            qpos = self.where_peaks[chrom]
            for i, j in self.diff_peaks[chrom]:
                above_cutoff = qpos[np.nonzero(self.diff_qvalues[chrom] >= self.cutoff)[0]]
                print '%s\t%d\t%d' % (chrom, 
                                      self.pos[chrom][above_cutoff[i]-1],
                                      self.pos[chrom][above_cutoff[j]] )
        return

    def write_peaks(self, xls=None, bed=None, name_prefix="%s_peak_", name="MACS",
                    description='%s', trackline=True):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        if self.has_peakio: 
            return self.write_peaks2(xls=xls, bed=bed, name_prefix=name_prefix, 
                              name=name, description=description,
                              trackline=trackline)
        if xls is not None:
            xlswrite = xls.write
        else:
            xlswrite = do_nothing
        if bed is not None:
            bedwrite = bed.write
        else:
            bedwrite = do_nothing
        xlswrite("# values are maxmimum in region\n")
        xlswrite("# Depth multiplier used: %f (treat/control values are after multiplication)\n") % self.cond1_depth
        xlswrite("# log10_fold_change is positive if t1 > t2\n")
        tc_method = self.track_scoring_method
        xlswrite("\t".join(("chr", "start", "end", "length",
                         "log2.fold.change", "-log10.diff.pvalue",
                         "-log10.diff.qvalue",
                         "diff.log10LR", "diff.peakname",
                         "treat1", "control1", "log2.fold.enrichment1",
                         "-log10.%svalue1" % tc_method,
                         "treat2", "control2", "log2.fold.enrichment2",
                         "-log10.%svalue2" % tc_method))+"\n")
        
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        trackcontents = (name.replace("\"", "\\\""), desc.replace("\"", "\\\""))
        if trackline:
            try: bedwrite('track name="%s (peaks)" description="%s" visibility=1\n' % trackcontents)
            except: bedwrite('track name=MACS description=Unknown\n') 

        log10 = np.log10
        log2 = np.log2
        median = np.median
        for chrom in sorted(self.diff_peaks.keys()):
            pos = self.pos[chrom]
            t1 = self.t1[chrom]
            t2 = self.t2[chrom]
            c1 = self.c1[chrom]
            c2 = self.c2[chrom]
            tvsc1 = self.tvsc1[chrom]
            tvsc2 = self.tvsc2[chrom]
            diff_pvalues = self.t1vs2[chrom]
            diff_qvalues = self.diff_qvalues[chrom]
            diff_logLR = self.tlogLR[chrom]
            qpos = self.where_peaks[chrom]
#            above_cutoff = np.where(diff_qvalues >= self.cutoff)[0]
            if self.diff_scoring_method == 'q':
                above_cutoff = np.where(self.diff_qvalues[chrom] >= 
                                        self.cutoff)[0].astype('int32')
            elif self.diff_scoring_method == 'p':
                above_cutoff = np.where(self.t1vs2[chrom][self.where_peaks[chrom]] >= self.cutoff)[0].astype('int32')
            for first_i, last_i, n_peak in self.diff_peaks[chrom]:
                start_i = above_cutoff[first_i]
                end_i = above_cutoff[last_i]
                pos_start = qpos[start_i] - 1
                pos_end = qpos[end_i]
                
                start = pos[pos_start]
                end = pos[pos_end]
                if start > end:
                    start = 0
                    pos_start = 0
                t1s = t1[pos_start:(pos_end+1)]
                c1s = c1[pos_start:(pos_end+1)]
                t2s = t2[pos_start:(pos_end+1)]
                c2s = c2[pos_start:(pos_end+1)]
                fold_changes = (t1s+self.pseudocount) / (t2s+self.pseudocount)
                if log2(median(fold_changes)) > 0:
                    log2_fold_change = log2(fold_changes).max()
                else:
                    log2_fold_change = log2(fold_changes).min()
                this_dpvalue = diff_pvalues[pos_start:(pos_end+1)].max()
                this_dqvalue = diff_qvalues[first_i:(last_i+1)].max()
                this_dlogLR = diff_logLR[pos_start:(pos_end+1)].max()
                peakname = "%s%d" % (peakprefix, n_peak)
                max_t1 = t1s.max()
                max_c1 = c1s.max()
                if max_t1 > max_c1: log2_fe1 = log2(t1s / c1s).max()
                else: log2_fe1 = log2(t1s / c1s).min()
                max_t2 = t2s.max()
                max_c2 = c2s.max()
                if max_t1 > max_c1: log2_fe2 = log2(t2s / c2s).max()
                else: log2_fe2 = log2(t2s / c2s).min()
                tc_value1 = tvsc1[pos_start:(pos_end+1)].max()
                tc_value2 = tvsc2[pos_start:(pos_end+1)].max()
                #chr,start,end,length, log10fold_change, diff.pvalue, diff.qvalue,
                #diff.logLR, name,
                #treat1, control1, fold_enrichment1, -log10(p/qvalue1)
                #treat2, control2, fold_enrichment2, -log10(p/qvalue2)
                xlswrite("%s\t%d\t%d\t%d" % (chrom, start+1, end, end - start))
                xlswrite("\t%.5f" % log2_fold_change)
                xlswrite("\t%.5f" % this_dpvalue)
                xlswrite("\t%.5f" % this_dqvalue)
                xlswrite("\t%.5f" % this_dlogLR)
                xlswrite("\t%s" % peakname)
                xlswrite("\t%.5f" % max_t1)
                xlswrite("\t%.5f" % max_c1)
                xlswrite("\t%.5f" % log2_fe1)
                xlswrite("\t%.5f" % tc_value1)
                xlswrite("\t%.5f" % max_t2)
                xlswrite("\t%.5f" % max_c2)
                xlswrite("\t%.5f" % log2_fe2)
                xlswrite("\t%.5f" % tc_value2)
                xlswrite("\n")
                bedwrite("%s\t%d\t%d\t%s\t%.5f\n" %
                         (chrom, start, end, peakname,
                          this_dqvalue))
        return
    
    
    def write_peaks2(self, xls=None, bed=None, name_prefix="%s_peak_", name="MACS",
                    description='%s', trackline=True):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        cdef:
            list peaks1, peaks2
            object peak1, peak2
            np.ndarray[np.int32_t] peaks1_selection, peaks2_selection
            np.ndarray[np.int32_t] i_peaks1, i_peaks2
            np.ndarray[np.int32_t] pos, which_peaks1, which_peaks2
            np.ndarray[np.float32_t] t1, c1, t2, c2
            np.ndarray[np.float32_t] diff_pvalues, diff_qvalues, diff_logLR
            np.ndarray[np.int32_t] qpos, above_cutoff
            int start_i, end_i, first_i, last_i
            int i1, i2
            int qpos_i
            int peak1_pos_i, peak2_pos_i
            int peak_i, peak_ii, peak1_i, peak2_i
            int npeaks1, npeaks2
            int peak1_selection_i, peak2_selection_i
            int n_peak
            bool peak1_present, peak2_present
            float this_dpvalue, this_dqvalue, this_dlogLR
            float pseudocount
        assert self.has_peakio(), "No information on peaks"
        logging.captureWarnings(True) 
        pseudocount = float(self.pseudocount)
        if xls is not None:
            xlswrite = xls.write
        else:
            xlswrite = do_nothing
        if bed is not None:
            bedwrite = bed.write
        else:
            bedwrite = do_nothing
        xlswrite("# summit is defined as greatest summit from greatest sample in region \n")
        xlswrite("# values are reported for the summit, except for \n")
        xlswrite("# fold_enrichment, sample/control qvalues, peaknames, which are copied from the peak info \n")
        xlswrite("# Depth multiplier used: %f (treat/control values are after multiplication)\n" % self.cond1_depth)
        xlswrite("# differential values are reported at the taller sample peak\n")
        xlswrite("# log2_fold_change is positive if t1 > t2\n")
        xlswrite("\t".join(("chr", "start", "end", "length", "summit",
                         "log2.fold.change","log2.fold.change.w.psuedocounts",
                         "-log10.diff.pvalue",
                         "-log10.diff.qvalue",
                         "diff.log10LR", "diff.peakname",
                         "treat1", "control1", "log2.fold.enrichment1",
                         "-log10.qvalue1", "peakname1", "summit1",
                         "treat2", "control2", "log2.fold.enrichment2",
                         "-log10.qvalue2", "peakname2", "summit2",
                         ))+"\n")
        
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        trackcontents = (name.replace("\"", "\\\""), desc.replace("\"", "\\\""))
        if trackline:
            try: bedwrite('track name="%s (peaks)" description="%s" visibility=1\n' % trackcontents)
            except: bedwrite('track name=MACS description=Unknown\n') 

        log10 = np.log10
        log2 = np.log2
        median = np.median
        for chrom in sorted(self.diff_peaks.keys()):
            peak1_pos_i = 0
            peak2_pos_i = 0
            qpos_i = 0
            pos = self.pos[chrom]
            qpos = self.where_peaks[chrom]
            which_peaks1 = self.which_peaks1[chrom]
            which_peaks2 = self.which_peaks2[chrom]
            try: peaks1 = self.p1io.get_data_from_chrom(chrom)
            except KeyError: peaks1 = []
            npeaks1 = len(peaks1)
            try: peaks2 = self.p2io.get_data_from_chrom(chrom)
            except KeyError: peaks2 = []
            npeaks2 = len(peaks2)
            t1 = self.t1[chrom]
            t2 = self.t2[chrom]
            c1 = self.c1[chrom]
            c2 = self.c2[chrom]
            diff_pvalues = self.t1vs2[chrom]
            diff_qvalues = self.diff_qvalues[chrom]
            diff_logLR = self.tlogLR[chrom]
            if self.diff_scoring_method == 'q':
                above_cutoff = np.where(self.diff_qvalues[chrom] >= 
                                        self.cutoff)[0].astype('int32')
            elif self.diff_scoring_method == 'p':
                above_cutoff = np.where(self.t1vs2[chrom][self.where_peaks[chrom]] >= self.cutoff)[0].astype('int32')
            # use some extra memory so we don't have to reallocate every time
            peaks1_selection = np.ndarray(npeaks1, 'int32')
            peaks2_selection = np.ndarray(npeaks2, 'int32')
            for first_i, last_i, n_peak in self.diff_peaks[chrom]:
                # this is going to be one entry
                start_i = above_cutoff[first_i]
                end_i = above_cutoff[last_i]
                pos_start = qpos[start_i] - 1
                pos_end = qpos[end_i]
                i_peaks1 = np.unique(which_peaks1[start_i:(end_i + 1)]).astype('int32')
                try:
                    if i_peaks1[0] == -1: i_peaks1 = i_peaks1[1:]
                except IndexError: pass
                i_peaks2 = np.unique(which_peaks2[start_i:(end_i + 1)]).astype('int32')
                try:
                    if i_peaks2[0] == -1: i_peaks2 = i_peaks2[1:]
                except IndexError: pass

                # get additional subpeaks as needed
                # i1 is the number of subpeaks in peaks1
                i1 = _get_all_subpeaks(peaks1, i_peaks1, peaks1_selection)
#                i1 = peaks1_selection
                i2 = _get_all_subpeaks(peaks2, i_peaks2, peaks2_selection)
#                i2 = peaks2_selection.size
#                i2 = 0
#                for peak_i in i_peaks2:
#                    start = peaks2[peak_i]["start"]
##                    these_peaks2.append( peaks2[peak_i] )
#                    peaks2_selection[i2] = peak_i
#                    i2 += 1
#                    for peak_ii in range(peak_i + 1, npeaks2):
#                        if start == peaks2[peak_ii]["start"]:
##                            these_peaks1.append( peaks2[peak_ii] )
#                            peaks2_selection[i2] = peak_i
#                            i2 += 1
#                        else: break
                # find the best overlapping subpeak
                if i1 > 0:
                    peak1_present = True
                    peak1_selection_i = 0
                    for j in range(1, i1):
                        if peaks1[peaks1_selection[j]]["pileup"] > \
                           peaks1[peaks1_selection[peak1_selection_i]]["pileup"]:
                            peak1_selection_i = j
                    peak1_i = peaks1_selection[peak1_selection_i]
                    peak1 = peaks1[peak1_i]
                    peak1_pos_i = peak1_pos_i + pos[peak1_pos_i:].searchsorted(peak1["summit"])
                else: peak1_present = False
                
                if i2 > 0:
                    peak2_present = True
                    peak2_selection_i = 0
                    for j in range(1, i2):
                        if peaks2[peaks2_selection[j]]["pileup"] > \
                           peaks2[peaks2_selection[peak2_selection_i]]["pileup"]:
                            peak2_selection_i = j
                    peak2_i = peaks2_selection[peak2_selection_i]
                    peak2 = peaks2[peak2_i]
                    peak2_pos_i = peak2_pos_i + pos[peak2_pos_i:].searchsorted(peak2["summit"])
                else: peak2_present = False
                
                if not peak1_present:
                    peak_pos_i = peak2_pos_i 
                elif not peak2_present:
                    peak_pos_i = peak1_pos_i
                else: # peak in both samples
                    if t1[peak1_pos_i] > t2[peak2_pos_i]:
                        peak_pos_i = peak1_pos_i
                    else:
                        peak_pos_i = peak2_pos_i
                
                # back to differential peak region
                start = pos[pos_start]
                end = pos[pos_end]
                if start > end:
                    start = 0
                    pos_start = 0
                #log2_fold_change = log2(t1[peak_pos_i]) - log2(t2[peak_pos_i]) # it will fail if some value is zero.
                
                log2_fold_change_w_pc = log2(t1[peak_pos_i] + pseudocount) - \
                                        log2(t2[peak_pos_i] + pseudocount)
                this_dpvalue = diff_pvalues[peak_pos_i]
                qpos_i = qpos[qpos_i:].searchsorted(peak_pos_i)
                this_dqvalue = diff_qvalues[qpos_i]
#                this_dqvalue = diff_qvalues[qpos.searchsorted(peak_pos_i)]
                this_dlogLR = diff_logLR[peak_pos_i]
#                fold_changes = t1s / t2s
#                if log2(median(fold_changes)) > 0:
#                    log2_fold_change = log2(fold_changes).max()
#                else:
#                    log2_fold_change = log2(fold_changes).min()
#                this_dpvalue = diff_pvalues[pos_start:(pos_end+1)].max()
#                this_dqvalue = diff_qvalues[first_i:(last_i+1)].max()
#                this_dlogLR = diff_logLR[pos_start:(pos_end+1)].max()
                peakname = "%s%d" % (peakprefix, n_peak)
                
#                max_t1 = t1s.max()
#                max_c1 = c1s.max()
#                if max_t1 > max_c1: log2_fe1 = log2(t1s / c1s).max()
#                else: log2_fe1 = log2(t1s / c1s).min()
#                max_t2 = t2s.max()
#                max_c2 = c2s.max()
#                if max_t1 > max_c1: log2_fe2 = log2(t2s / c2s).max()
#                else: log2_fe2 = log2(t2s / c2s).min()
#                tc_value1 = tvsc1[pos_start:(pos_end+1)].max()
#                tc_value2 = tvsc2[pos_start:(pos_end+1)].max()
                #chr,start,end,length, log10fold_change, diff.pvalue, diff.qvalue,
                #diff.logLR, name,
                #treat1, control1, fold_enrichment1, -log10(p/qvalue1)
                #treat2, control2, fold_enrichment2, -log10(p/qvalue2)
                xlswrite("%s\t%d\t%d\t%d" % (chrom, start+1, end, end - start))
                xlswrite("\t%d" % pos[peak_pos_i])
                #xlswrite("\t%.5f" % log2_fold_change)
                xlswrite("\t%.5f" % log2_fold_change_w_pc)
                xlswrite("\t%.5f" % this_dpvalue)
                xlswrite("\t%.5f" % this_dqvalue)
                xlswrite("\t%.5f" % this_dlogLR)
                xlswrite("\t%s" % peakname)
                xlswrite("\t%.5f" % t1[peak_pos_i])
                xlswrite("\t%.5f" % c1[peak_pos_i])
                if peak1_present:
                    xlswrite("\t%.5f" % log2(peak1["fc"]))
                    xlswrite("\t%.5f" % peak1["qscore"])
                    xlswrite("\t%s" % peak1["name"])
                    xlswrite("\t%d" % peak1["summit"])
                else: xlswrite("\tNA\tNA\tNA\tNA")
                xlswrite("\t%.5f" % t2[peak_pos_i])
                xlswrite("\t%.5f" % c2[peak_pos_i])
                if peak2_present:
                    xlswrite("\t%.5f" % log2(peak2["fc"]))
                    xlswrite("\t%.5f" % peak2["qscore"])
                    xlswrite("\t%s" % peak2["name"])
                    xlswrite("\t%d" % peak2["summit"])
                else: xlswrite("\tNA\tNA\tNA\tNA")
                xlswrite("\n")
                bedwrite("%s\t%d\t%d\t%s\t%.5f\n" %
                         (chrom, start, end, peakname,
                          this_dqvalue))
        logging.captureWarnings(False)
        return

    def has_peakio(self):
        """see whether peaks1 or peaks2 are stored"""
        return (len(self.which_peaks1.keys()) + len(self.which_peaks2.keys())) > 0
    

    def write_peaks_by_summit(self, xls1, xls2,
                                name_prefix="%s_peak_", name="MACS"):
        """write a file with information about each summit and differential
        occupancy (separately for condition 1 and condition 2)
        """
        assert self.has_peakio(), "No information on peaks"
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        xls1.write("# Depth multiplier used: %f (treat/control values are after multiplication)\n" % self.cond1_depth)
        xls1.write("# the peak with the closest summit from the other sample is reported if a peak overlaps the summit\n")
        xls1.write("# log2_fold_change is positive if t1 > t2\n")
        xls2.write("# Depth multiplier used: %f (treat/control values are after multiplication)\n" % self.cond2_depth)
        xls2.write("# the peak with the closest summit from the other sample is reported if a peak overlaps the summit\n")
        xls2.write("# log2_fold_change is positive if t2 > t1\n")
        xls1.write("\t".join(("chr", "start", "end", "length", "summit",
                         "log2.fold.change","log2.fold.change.w.psuedocounts",
                         "-log10.diff.pvalue",
                         "-log10.diff.qvalue",
                         "diff.log10LR", "diff.peakname",
                         "treat1", "control1", "log2.fold.enrichment1",
                         "-log10.qvalue1", "peakname1", "summit1",
                         "treat2", "control2", "log2.fold.enrichment2",
                         "-log10.qvalue2", "peakname2", "summit2",
                         ))+"\n")
        xls2.write("\t".join(("chr", "start", "end", "length", "summit",
                         "log2.fold.change","log2.fold.change.w.psuedocounts",
                         "-log10.diff.pvalue",
                         "-log10.diff.qvalue",
                         "diff.log10LR", "diff.peakname",
                         "treat2", "control2", "log2.fold.enrichment2",
                         "-log10.qvalue2", "peakname2", "summit2",
                         "treat1", "control1", "log2.fold.enrichment1",
                         "-log10.qvalue1", "peakname1", "summit1",
                         ))+"\n")
        self._write_peaks_by_summit(xls1, self.p1io, self.p2io,
                                    self.which_peaks2, False,
                                    self.t1, self.c1,
                                    self.t2, self.c2, peakprefix)
        self._write_peaks_by_summit(xls2, self.p2io, self.p1io,
                                    self.which_peaks1, True,
                                    self.t2, self.c2,
                                    self.t1, self.c1, peakprefix)
        
    cdef _write_peaks_by_summit(self, xls, p1io, p2io,
                                dict which_peaks2,
                                bool flip_fc,
                                dict t1_by_chrom,
                                dict c1_by_chrom,
                                dict t2_by_chrom,
                                dict c2_by_chrom,
                                str peakprefix):
        cdef:
            list peaks1, peaks2
            str chrom
            int i, j, d, w, d_max, w_max, p, n_peaks, n
            int datalength
            int peak_i, peak_j
            int summit, start, end, length
            float this_dpvalue, this_dqvalue, this_dlogLR, score_value
            float pseudocount
            float log2_fold_change, log2_fold_change_w_pc
            np.ndarray[np.float32_t] t1, t2, c1, c2
            np.ndarray[np.int32_t] pos, which_peaks, summits1
            np.ndarray[np.float32_t] diff_pvalues, diff_qvalues, diff_logLR
            np.ndarray[np.int32_t, ndim=2] diff_peaks
        log2 = np.log2
        logging.captureWarnings(True) 
        pseudocount = float(self.pseudocount)
        for chrom in sorted(p1io.get_chr_names()):
            write = xls.write
            try: datalength = self.datalength[chrom]
            except KeyError: datalength = 0
            if datalength < 2:
                # we are missing data on this chromosome, just write original peaks
                peaks1 = p1io.get_data_from_chrom(chrom)
                for peak in peaks1:
                    start = peak["start"] + 1
                    end = peak["end"]
                    length = peak["length"]
                    write("%s\t%d\t%d\t%d" % (chrom,start,end,length))
                    write("\t%d" % peak["summit"] + 1) # summit position
                    write("\tNA\tNA\tNA\tNA\tNA\tNA")
                    # this part spits out stuff for sample 1
                    write("\t%.5f" % peak["pileup"])
                    write("\tNA")
                    write("\t%.5f" % log2(peak["fc"]))
                    write("\t%.5f" % peak["qscore"])
                    write("\t%s" % peak["name"])
                    write("\t%d" % peak["summit"] + 1)
                    # this stuff for peak 2
                    write("\tNA\tNA\tNA\tNA\tNA\tNA")
                    write("\n")
                continue
            peaks1 = p1io.get_data_from_chrom(chrom)
            try: peaks2 = p2io.get_data_from_chrom(chrom)
            except KeyError: peaks2 = []
            summits1 = np.array(map(itemgetter('summit'), peaks1),
                                dtype='int32')
            # note, summits should be unique and should be in order
            pos = self.pos[chrom]
            where_peaks = self.where_peaks[chrom]
            w = 0
            w_max = where_peaks.size
            which_peaks = which_peaks2[chrom]
            n_peaks = len(peaks2)
            t1 = t1_by_chrom[chrom]
            t2 = t2_by_chrom[chrom]
            c1 = c1_by_chrom[chrom]
            c2 = c2_by_chrom[chrom]
            diff_peaks = self.diff_peaks[chrom].copy()
            d_max = diff_peaks.shape[0]
            diff_pvalues = self.t1vs2[chrom]
            diff_qvalues = self.diff_qvalues[chrom]
            diff_logLR = self.tlogLR[chrom]
            if self.diff_scoring_method == 'q':
                above_cutoff = np.where(self.diff_qvalues[chrom] >= 
                                        self.cutoff)[0].astype('int32')
            elif self.diff_scoring_method == 'p':
                above_cutoff = np.where(self.t1vs2[chrom][where_peaks] >= self.cutoff)[0].astype('int32')
            # change to coordinates for diff peak now
            for d in range(d_max):
                start = diff_peaks[d, 0]
                end = diff_peaks[d, 1]
                diff_peaks[d, 0] = pos[where_peaks[above_cutoff[start]] - 1]
                diff_peaks[d, 1] = pos[where_peaks[above_cutoff[end]]]
            if d_max > 0:
                if diff_peaks[0, 0] > diff_peaks[0, 1]:
                    diff_peaks[0, 0] = 0
            d = 0
            i = 0
            for j in range(summits1.size):
                if i == datalength: break #shouldn't happen
                summit = summits1[j]
                while True:
                    if i == datalength: break # shouldn't happen
                    elif summit >= pos[i]: i += 1
                    else:
                        # write record
                        peak = peaks1[j]
                        # check if this is a differentially occupied peak
                        while True:
                            if w == w_max - 1: break
                            elif summit < pos[where_peaks[w]]:
                                break
                            else:
                                w += 1
#                        print w, where_peaks[w], summit
                        start = peak["start"] + 1
                        end = peak["end"]
                        length = peak["length"]
                        write("%s\t%d\t%d\t%d" % (chrom,start,end,length))
                        write("\t%d" % (summit + 1)) # summit position
                        # only if summit >= diff_peaks[d, 0] and
                        # summit < diff_peaks[d, 1]
                        this_dpvalue = diff_pvalues[i]
                        this_dqvalue = diff_qvalues[w]
                        if self.diff_scoring_method == 'p':
                            score_value = this_dpvalue
                        elif self.diff_scoring_method == 'q':
                            score_value = this_dqvalue       
                        #log2_fold_change = log2(t1[i]) - log2(t2[i])
                        log2_fold_change_w_pc = log2(t1[i] + pseudocount) - \
                                                log2(t2[i] + pseudocount)
                        this_dlogLR = diff_logLR[i]
                        if flip_fc:
                            #write("\t%.5f" % -log2_fold_change)
                            write("\t%.5f" % -log2_fold_change_w_pc)
                        else:
                            #write("\t%.5f" % log2_fold_change)
                            write("\t%.5f" % log2_fold_change_w_pc)
                        write("\t%.5f" % this_dpvalue)
                        write("\t%.5f" % this_dqvalue)
                        write("\t%.5f" % this_dlogLR)
#                        print i, w, w_max, d, d_max
                        # this part to figure out which differential peak
                        while True:
                            if d == d_max - 1 or d_max == 0: break
                            elif summit < diff_peaks[d, 1]:
                                break
                            else:
                                d += 1
#                        print d, d_max, diff_peaks[d,:], summit
                        if score_value >= self.cutoff:
                            if d_max == 0:
                                write("\tNA")
                            elif summit >= diff_peaks[d, 0]:
                                diffpeakname = "%s%d" % (peakprefix,
                                                         diff_peaks[d, 2])
                                write("\t%s" % diffpeakname)
                            else:
                                write("\tNA")
                        else:
                            write("\tNA")
                        # this part spits out stuff for sample 1
                        write("\t%.5f" % t1[i])
                        write("\t%.5f" % c1[i])
                        write("\t%.5f" % log2(peak["fc"]))
                        write("\t%.5f" % peak["qscore"])
                        write("\t%s" % peak["name"])
                        write("\t%d" % (summit + 1))
                        # this stuff for peak 2
                        write("\t%.5f" % t2[i])
                        write("\t%.5f" % c2[i])
                        if pos[where_peaks[w]] > end:
                            peak_i = -1
                        else:
                            peak_i = which_peaks[w]
                        if peak_i == -1:
                            write("\tNA\tNA\tNA\tNA")
                        else:
                            # find the closest summit2 to this region
                            peak2_summit = peaks2[peak_i]["summit"]
                            peak_j = peak_i + 1
                            for peak_j in range(peak_i + 1, n_peaks):
                                if peaks2[peak_j]["end"] != end: break
                            peak2_summits = np.array(map(itemgetter("summit"),
                                                         peaks2[peak_i:peak_j]),
                                                     dtype='int32')
                            peak_i += np.abs(peak2_summits - summit).argmin()
                            peak2 = peaks2[peak_i]
                            write("\t%.5f" % log2(peak2["fc"]))
                            write("\t%.5f" % peak2["qscore"])
                            write("\t%s" % peak2["name"])
                            write("\t%d" % (peak2["summit"] + 1))
                                
                        write("\n")
                        break
                

    def write_bedgraphs(self, logLR=None, pvalue=None, logFC=None,
                          str name="MACS",
                          str description='%s', bool trackline=True):
        """Write logLR and diff pvalue data to in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        shift will be used to shift the coordinates. default: 0
        """
        cdef:
            int i
            str chrom
            np.ndarray[np.int32_t] pos 
            np.ndarray[np.float32_t] value
        
        isfinite = np.isfinite
        if trackline:
            trackcontents = (name.replace("\"", "\\\""), description.replace("\"", "\\\""))
            if logLR is not None:
                logLR.write("track type=bedGraph name=\"%s\" description=\"log10-likelihood ratio %s\" visibility=2 alwaysZero=on\n" % trackcontents)
            if pvalue is not None:
                pvalue.write("track type=bedGraph name=\"%s\" description=\"-log10(pvalue)%s\" visibility=2 alwaysZero=on\n" % trackcontents)
            if logFC is not None:
                logFC.write("track type=bedGraph name=\"%s\" description=\"log10(sample1/sample2) %s\" visibility=2 alwaysZero=on\n" % trackcontents)
                
        if logLR is not None:
            for chrom in sorted(self.tlogLR.keys()):
                pos = self.pos[chrom]
                value = self.tlogLR[chrom]
                if pos.size > 0:
                    if isfinite(value[0]):
                        logLR.write("%s\t%d\t%d\t%.5f\n" % (chrom, 0, pos[0],
                                                      value[0]))
                    for i in range(1, pos.size):
                        if isfinite(value[i]):
                            logLR.write("%s\t%d\t%d\t%.5f\n" % (chrom, pos[i-1], pos[i],
                                                          value[i]))
                    
        if pvalue is not None:
            for chrom in sorted(self.t1vs2.keys()):
                pos = self.pos[chrom]
                value = self.t1vs2[chrom]
                if pos.size > 0:
                    if isfinite(value[0]):
                        pvalue.write("%s\t%d\t%d\t%.5f\n" % (chrom, 0, pos[0],
                                                      value[0]))
                    for i in range(1, pos.size):
                        if isfinite(value[i]):
                            pvalue.write("%s\t%d\t%d\t%.5f\n" % (chrom, pos[i-1], pos[i],
                                                                 value[i]))
                    
        if logFC is not None:
            for chrom in sorted(self.pos.keys()):
                pos = self.pos[chrom]
                value = np.log2(self.t1[chrom] / self.t2[chrom])
                if pos.size > 0:
                    if isfinite(value[0]):
                        logFC.write("%s\t%d\t%d\t%.5f\n" % (chrom, 0, pos[0],
                                                      value[0]))
                    for i in range(1, pos.size):
                        if isfinite(value[i]):
                            logFC.write("%s\t%d\t%d\t%.5f\n" % (chrom, pos[i-1], pos[i],
                                                          value[i]))
                
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

cdef inline int _get_all_subpeaks(list peaks,
                                  np.ndarray[np.int32_t] i_peaks,
                                  np.ndarray[np.int32_t] peaks_selection):
    """return the number of subpeaks, modify peaks_selection in place"""
    cdef:
        int i
        int peak_i, peak_ii
        int start, start2
    i = 0
    for peak_i in i_peaks:
        start = peaks[peak_i]["start"]
        #these_peaks1.append( peaks1[peak_i] )
        peaks_selection[i] = peak_i
        i += 1
        peak_ii = peak_i + 1
        while True:
            try: start2 = peaks[peak_ii]["start"]
            except IndexError: break
            if start == start2:
                peaks_selection[i] = peak_ii
                i += 1
                peak_ii += 1
            else: break
    return i
