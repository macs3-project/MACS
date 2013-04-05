# cython: profile=True
# Time-stamp: <2013-04-05 15:40:48 Tao Liu>

"""Module for Calculate Scores.

Copyright (c) 2013 Tao Liu <taoliu@jimmy.harvard.edu>

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
#from copy import copy

from operator import itemgetter

from cpython cimport bool
#from scipy.signal import fftconvolve
from MACS2.cSignal import maxima, enforce_valleys, enforce_peakyness
#np_convolve = np.convolve

# Experimental
#from scipy.stats import chi2

#from cython.parallel import parallel, prange
cimport cython

from libc.math cimport log10,log, floor, ceil

from MACS2.Constants import BYTE4, FBYTE4, array
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO, BroadPeakIO, parse_peakname
from MACS2.IO.cFixWidthTrack import FWTrackIII
from MACS2.IO.cPairedEndTrack import PETrackI

from MACS2.hashtable import Int64HashTable, Float64HashTable

import logging

# ------------------------------------
# constants
# ------------------------------------
__version__ = "scoreCalculate $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "scoreTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
def do_nothing(*args, **kwargs):
    pass

LOG10_E = 0.43429448190325176
pscore_khashtable = Int64HashTable()

cdef inline np.ndarray apply_multiple_cutoffs ( list multiple_score_arrays, list multiple_cutoffs ):
    cdef:
        int i
        np.ndarray ret

    ret = multiple_score_arrays[0] > multiple_cutoffs[0]
    
    for i in range(1,len(multiple_score_arrays)):
        ret += multiple_score_arrays[i] > multiple_cutoffs[i]

    return ret

cdef inline list get_from_multiple_scores ( list multiple_score_arrays, int index ):
    cdef:
        list ret = []
        int i

    for i in range(len(multiple_score_arrays)):
        ret.append(multiple_score_arrays[i][index])
    return ret


cdef inline double get_pscore ( int observed, double expectation ):
    """Get p-value score from Poisson test. First check existing
    table, if failed, call poisson_cdf function, then store the result
    in table.
    
    """
    cdef:
        double score
        long key_value
    
    key_value = hash( (observed, expectation ) )
    try:
        return pscore_khashtable.get_item(key_value)
    except KeyError:
        score = -1*poisson_cdf(observed,expectation,False,True)
        pscore_khashtable.set_item(key_value, score)
        return score

logLR_khashtable = Int64HashTable()

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
    
cdef inline double oldlogLR ( double x, double y ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Then store the values in integer form =
    100*logLR. Set minus sign for depletion.
    
    """
    cdef:
        double s
        long key_value
    
    key_value = hash( (x, y ) )
    try:
        return logLR_khashtable.get_item( key_value )
    except KeyError:
        if x > y:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
        elif x < y:
            s = (x*(-log(x)+log(y))-y+x)*LOG10_E
        else:
            s = 0
        logLR_khashtable.set_item(key_value, s)
        return s

cdef inline double logLR_adjusted (int x, int y):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Then store the values in integer form =
    100*logLR. Set minus sign for depletion.
    
    """
    cdef:
        double s
        long key_value
    
    key_value = hash( (x, y) )
    try:
        return logLR_khashtable.get_item( key_value )
    except KeyError:
        if y > x: y, x = x, y
        if x==y: s = 0
        
        else: s = (x*(log(x)-log(y))+y-x)*LOG10_E
        logLR_khashtable.set_item(key_value, s)
        return s
    
cdef inline double logLR ( double x, double y ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Then store the values in integer form =
    100*logLR. Set minus sign for depletion.
    
    """
    cdef:
        double s
        long key_value
    
    key_value = hash( (x, y ) )
    try:
        return logLR_khashtable.get_item( key_value )
    except KeyError:
        if y > x: y, x = x, y
        if x==y: s = 0
        
        else: s = (x*(log(x)-log(y))+y-x)*LOG10_E
        logLR_khashtable.set_item(key_value, s)
        return s
    
cdef inline double get_logFE ( float x, float y ):
    """ return 100* log10 fold enrichment with +1 pseudocount.
    """
    return log10( x/y )

cdef inline float get_subtraction ( float x, float y):
    """ return subtraction.
    """
    return x - y

# ------------------------------------
# Classes
# ------------------------------------
cdef class ScoreCalculator:
    """A unit to calculate scores. 

    It will use hashtable to accelerate re-calculation of scores. 

    """
    cdef:
        str chrom                        # name of chromosome
        char scoring_method              # method for calculating scores.
        dict pqtable                     # remember pvalue->qvalue convertion
        bool pvalue_all_done             # whether the pvalue of whole genome is all calculated. If yes, it's OK to calculate q-value.
        bool trackline                   # whether trackline should be saved in bedGraph
        double treat_edm                 # seq depth in million of treatment
        double ctrl_edm                  # seq depth in million of control
        char normalization_method        # scale to control? scale to treatment? both scale to 1million reads?
        float pseudocount                # the pseudocount used to calcuate logLR, FE or logFE
        float cutoff                     # cutoff to call peaks
        list chr_pos_treat_ctrl              # tmp [position, treat_pileup, ctrl_pileup] for a given chromosome
        int  d
        list ctrl_d_s
        float treat_scale_factor
        list ctrl_scale_factor_s
        bool halfextension
        float lambda_bg
        bool shiftcontrol
        list chromosomes
        object treat
        object ctrl
        bool save_bedGraph
        str bedGraph_filename_prefix # will add _pileup.bdg for treatment and _lambda.bdg for control
        object bedGraph_treat
        object bedGraph_ctrl
        bool no_lambda_flag
        bool PE_mode

    def __init__ (self, treat, ctrl, float treat_depth, float ctrl_depth, 
                  int d = 200, list ctrl_d_s = [200, 1000, 10000], 
                  float treat_scale_factor = 1.0, list ctrl_scale_factor_s = [1.0, 0.2, 0.02], 
                  bool stderr_on = False, 
                  float pseudocount = 1.0, 
                  bool halfextension = False, 
                  float lambda_bg = 0, 
                  bool shiftcontrol = False,
                  bool save_bedGraph = False,
                  str  bedGraph_filename_prefix = ""):
        """Initialize.

        A calculator is unique to each comparison of treat and
        control. Treat_depth and ctrl_depth should not be changed
        during calculation.

        treat and ctrl are two FWTrackIII or PETrackI objects.

        treat_depth and ctrl_depth are effective depth in million:
                                    sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        d, sregion, lregion: d is the fragment size, sregion is the
                             small region size, lregion is the large
                             region size
                                    
        pseudocount: a pseudocount used to calculate logLR, FE or
                     logFE. Please note this value will not be changed
                     with normalization method. So if you really want
                     to set pseudocount 1 per million reads, set it
                     after you normalize treat and control by million
                     reads by `change_normalizetion_method(ord('M'))`.

        """
        cdef:
            set chr1, chr2

        if isinstance(treat, FWTrackIII):
            self.PE_mode = False
        elif isinstance(ctrl, PETrackI):
            self.PE_mode = True
        else:
            raise Exception("Should be FWTrackIII or PETrackI object!")
        
        self.treat = treat
        if ctrl:
            self.ctrl = ctrl
        else:                   # while there is no control
            self.ctrl = treat
        self.trackline = False
        self.treat_edm = treat_depth
        self.ctrl_edm = ctrl_depth
        self.d = d              # note, self.d doesn't make sense in PE mode
        self.ctrl_d_s = ctrl_d_s# note, self.d doesn't make sense in PE mode
        self.treat_scale_factor = treat_scale_factor
        self.ctrl_scale_factor_s= ctrl_scale_factor_s
        self.halfextension = halfextension
        self.shiftcontrol = shiftcontrol
        self.lambda_bg = lambda_bg
        self.pqtable = {}
        self.save_bedGraph = save_bedGraph
        self.bedGraph_filename_prefix =  bedGraph_filename_prefix

        if not self.ctrl_d_s or not self.ctrl_scale_factor_s:
            self.no_lambda_flag = True
        else:
            self.no_lambda_flag = False

        #scoring_method:  p: -log10 pvalue;
        #                 q: -log10 qvalue;
        #                 l: log10 likelihood ratio ( minus for depletion )
        #                 f: log10 fold enrichment
        #                 F: linear fold enrichment        
        #                 d: subtraction
        #                 m: fragment pileup per million reads
        #                 N: not set
        self.scoring_method = ord("N")

        #normalization_method: T: scale to depth of treatment;
        #                      C: scale to depth of control;
        #                      M: scale to depth of 1 million;
        #                      N: not set/ raw pileup
        self.normalization_method = ord("N")

        self.pseudocount = pseudocount

        chr1 = set(self.treat.get_chr_names())
        chr2 = set(self.ctrl.get_chr_names())
        self.chromosomes = list(chr1.intersection(chr2))

    cpdef set_pseudocount( self, float pseudocount ):
        self.pseudocount = pseudocount
        
    cpdef enable_trackline( self ):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True

    cdef __pileup_treat_ctrl_a_chromosome ( self, str chrom ):
        """After this function is called, self.chr_pos_treat_ctrl will
        be reset and assigned to the pileup values of the given
        chromosome.
        
        """
        cdef:
            list treat_pv, ctrl_pv
            long i

        assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom

        # reset or clean existing self.chr_pos_treat_ctrl
        if self.chr_pos_treat_ctrl:     # not a beautiful way to clean
            self.chr_pos_treat_ctrl[0].resize(10000,refcheck=False)
            self.chr_pos_treat_ctrl[1].resize(10000,refcheck=False)
            self.chr_pos_treat_ctrl[2].resize(10000,refcheck=False)
            self.chr_pos_treat_ctrl[0].resize(0,refcheck=False)
            self.chr_pos_treat_ctrl[1].resize(0,refcheck=False)
            self.chr_pos_treat_ctrl[2].resize(0,refcheck=False)            
        if self.PE_mode:
            treat_pv = self.treat.pileup_a_chromosome ( chrom, [self.treat_scale_factor,], baseline_value = 0.0 )
        else:
            treat_pv = self.treat.pileup_a_chromosome( chrom, [self.d,], [self.treat_scale_factor,], baseline_value = 0.0,
                                                       directional = True, halfextension = self.halfextension )
        if not self.no_lambda_flag:
            if self.PE_mode:
                ctrl_pv = self.ctrl.pileup_a_chromosome( chrom, self.ctrl_scale_factor_s, baseline_value = self.lambda_bg )
            else:
                ctrl_pv = self.ctrl.pileup_a_chromosome( chrom, self.ctrl_d_s, self.ctrl_scale_factor_s,
                                                         baseline_value = self.lambda_bg,
                                                         directional = self.shiftcontrol,
                                                         halfextension = self.halfextension )
        else:
            ctrl_pv = [treat_pv[0][-1:], pyarray(FBYTE4,[self.lambda_bg,])] # set a global lambda

        self.chr_pos_treat_ctrl = self.__chrom_pair_treat_ctrl( treat_pv, ctrl_pv)
        # clean treat_pv and ctrl_pv
        treat_pv = []
        ctrl_pv  = []  
        return

    cdef list __chrom_pair_treat_ctrl ( self, treat_pv, ctrl_pv ):
        """*private* Pair treat and ctrl pileup for each region.

        return [p, t, c] list
        """
        cdef:
            list ret
            long pre_p, p1, p2, index_ret
            float v1, v2

        p1n = iter(treat_pv[0]).next         # assign the next function to a viable to speed up
        v1n = iter(treat_pv[1]).next

        p2n = iter(ctrl_pv[0]).next         # assign the next function to a viable to speed up
        v2n = iter(ctrl_pv[1]).next

        chrom_max_len = len(treat_pv[0])+len(ctrl_pv[0])

        ret = [ np.zeros( chrom_max_len, dtype="int32" ),
                np.zeros( chrom_max_len, dtype="float32" ),
                np.zeros( chrom_max_len, dtype="float32" ) ] # p, t, c
            
        pre_p = 0
        index_ret = 0
        
        try:
            p1 = p1n()
            v1 = v1n()
            p2 = p2n()
            v2 = v2n()
            #if v2 == 0:
            #    print p2, v2
            while True:
                if p1 < p2:
                    # clip a region from pre_p to p1, then set pre_p as p1.
                    ret[0][index_ret] = p1
                    ret[1][index_ret] = v1
                    ret[2][index_ret] = v2                    
                    pre_p = p1
                    index_ret += 1
                    # call for the next p1 and v1
                    p1 = p1n()
                    v1 = v1n()
                elif p2 < p1:
                    # clip a region from pre_p to p2, then set pre_p as p2.
                    ret[0][index_ret] = p2
                    ret[1][index_ret] = v1
                    ret[2][index_ret] = v2                                        
                    pre_p = p2
                    index_ret += 1
                    # call for the next p2 and v2
                    p2 = p2n()
                    v2 = v2n()
                elif p1 == p2:
                    # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                    ret[0][index_ret] = p1
                    ret[1][index_ret] = v1
                    ret[2][index_ret] = v2
                    pre_p = p1
                    index_ret += 1
                    # call for the next p1, v1, p2, v2.
                    p1 = p1n()
                    v1 = v1n()
                    p2 = p2n()
                    v2 = v2n()
        except StopIteration:
            # meet the end of either bedGraphTrackI, simply exit
            pass
        ret[0].resize( index_ret, refcheck=False)
        ret[1].resize( index_ret, refcheck=False)
        ret[2].resize( index_ret, refcheck=False)
        return ret

    cdef np.ndarray __cal_score ( self, array1, array2, cal_func ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = cal_func( array1[i], array2[i] )
        return s

    cdef dict __cal_pvalue_qvalue_table ( self ):
        """After this function is called, self.pqtable is built. All
        chromosomes will be iterated. So it will take some time.
        
        """
        cdef:
            str chrom
            np.ndarray pos_array, treat_array, ctrl_array
            dict pvalue_stat = {}
            long n, pre_p, this_p, length, j, pre_l, l, i
            double this_v, pre_v, v, q, pre_q
            long N, k
            double f
            list unique_values
        
        for chrom in self.chromosomes:
            prev_p = 0
            #logging.info("#3 pileup %s" % chrom)
            self.__pileup_treat_ctrl_a_chromosome( chrom )
            [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

            #logging.info("#3 add up pvalue regions for %s" % chrom)
            for i in range(pos_array.shape[0]):
                #if ctrl_array[i] == 0:
                #    print "cal c p:", chrom, i, pos_array[i], treat_array[i], ctrl_array[i] 
                this_v = get_pscore( int(treat_array[i]), ctrl_array[i] )
                if pvalue_stat.has_key(this_v):
                    pvalue_stat[ this_v ] += pos_array[ i ] - prev_p 
                else:
                    pvalue_stat[ this_v ] = pos_array[ i ] - prev_p 
                #print pos_array[ i ], prev_p, pos_array[ i ] - prev_p
                prev_p = pos_array[ i ]
    
        N = sum(pvalue_stat.values()) # total length
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value

        self.pqtable = {}
        unique_values = sorted(pvalue_stat.keys(), reverse=True) #sorted(unique_values,reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = pvalue_stat[v]
            q = v + (log10(k) + f)
            q = max(0,min(pre_q,q))           # make q-score monotonic
            self.pqtable[ v ] = q
            pre_v = v
            pre_q = q
            k+=l
        
        return self.pqtable

    cpdef call_peaks ( self, list scoring_function_symbols, list score_cutoff_s, int min_length = 200, 
                       int max_gap = 50, bool call_summits = False ):
        """Call peaks for all chromosomes. Return a PeakIO object.
        
        scoring_function_s: symbols of functions to calculate score. 'p' for pscore, 'q' for qscore, 'f' for fold change, 's' for subtraction. for example: ['p', 'q']
        score_cutoff_s    : cutoff values corresponding to scoring functions
        min_length        : minimum length of peak
        max_gap           : maximum gap of 'insignificant' regions within a peak
        call_summits      : boolean. Whether or not call sub-peaks.
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        cdef:
            str chrom
            str s

        peaks = PeakIO()

        # prepare p-q table
        if not self.pqtable:
            logging.info("#3 Pre-compute pvalue-qvalue table...")
            self.__cal_pvalue_qvalue_table()

        # prepare bedGraph file
        if self.save_bedGraph:

            self.bedGraph_treat = open( self.bedGraph_filename_prefix + "_treat_pileup.bdg", "w" )
            self.bedGraph_ctrl = open( self.bedGraph_filename_prefix + "_control_lambda.bdg", "w" )
            logging.info ("#3 In the peak calling step, the following will be performed simultaneously:")
            logging.info ("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... %s" % self.bedGraph_filename_prefix + "_treat_pileup.bdg")
            logging.info ("#3   Write bedGraph files for control lambda (after scaling if necessary)... %s" % self.bedGraph_filename_prefix + "_control_lambda.bdg")

            if self.trackline:
                # this line is REQUIRED by the wiggle format for UCSC browser
                self.bedGraph_treat.write( "track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix )
                self.bedGraph_ctrl.write ( "track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix )

        logging.info("#3 Call peaks for each chromosome...")
        for chrom in self.chromosomes:
            self.__chrom_call_peak_using_certain_criteria ( peaks, chrom, scoring_function_symbols, score_cutoff_s, min_length, max_gap, call_summits, self.save_bedGraph )

        # close bedGraph file
        if self.save_bedGraph:
            self.bedGraph_treat.close()
            self.bedGraph_ctrl.close()
            self.save_bedGraph = False

        return peaks

    cdef __chrom_call_peak_using_certain_criteria ( self, peaks, str chrom, list scoring_function_s, list score_cutoff_s, int min_length, 
                                                   int max_gap, bool call_summits, bool save_bedGraph):
        """ Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        cdef:
            int i
            str s
            np.ndarray above_cutoff, above_cutoff_endpos, above_cutoff_startpos
            np.ndarray pos_array, treat_array, ctrl_array
            np.ndarray above_cutoff_index_array
            list score_array_s          # list to keep different types of scores
            list peak_content           #  to store information for a chunk in a peak region, it contains lists of: 1. left position; 2. right position; 3. treatment value; 4. control value; 5. list of scores at this chunk

        assert len(scoring_function_s) == len(score_cutoff_s), "number of functions and cutoffs should be the same!"
        
        peak_content = []           # to store points above cutoff

        # first, build pileup, self.chr_pos_treat_ctrl
        self.__pileup_treat_ctrl_a_chromosome( chrom )
        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        # while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
        if save_bedGraph:
            self.__write_bedGraph_for_a_chromosome ( chrom )

        # keep all types of scores needed
        score_array_s = []
        for i in range(len(scoring_function_s)):
            s = scoring_function_s[i]
            if s == 'p':
                score_array_s.append( self.__cal_pscore( treat_array, ctrl_array ) )
            elif s == 'q':
                score_array_s.append( self.__cal_qscore( treat_array, ctrl_array ) )
            elif s == 'f':
                score_array_s.append( self.__cal_FE( treat_array, ctrl_array ) )
            elif s == 's':
                score_array_s.append( self.__cal_subtraction( treat_array, ctrl_array ) )

        # get the regions with scores above cutoffs
        above_cutoff = np.nonzero( apply_multiple_cutoffs(score_array_s,score_cutoff_s) )[0] # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff_index_array = np.arange(pos_array.shape[0])[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return peaks

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], treat_array[above_cutoff_index_array[0]], ctrl_array[above_cutoff_index_array[0]], get_from_multiple_scores( score_array_s, above_cutoff_index_array[i])) )
        for i in range( 1,above_cutoff_startpos.size ):
            if above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                # append
                peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]], get_from_multiple_scores( score_array_s, above_cutoff_index_array[i]) ) )
            else:
                # close
                if call_summits:
                    self.__close_peak_with_subpeaks (peak_content, peaks, min_length, chrom, max_gap/2 )
                else:
                    self.__close_peak_wo_subpeaks   (peak_content, peaks, min_length, chrom, max_gap/2 )
                peak_content = [ (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]], get_from_multiple_scores( score_array_s, above_cutoff_index_array[i]) ), ]
            
        # save the last peak
        if not peak_content:
            return peaks
        else:
            if call_summits:
                self.__close_peak_with_subpeaks (peak_content, peaks, min_length, chrom, max_gap/2, score_cutoff_s ) # smooth length is 1/2 max-gap
            else:
                self.__close_peak_wo_subpeaks   (peak_content, peaks, min_length, chrom, max_gap/2, score_cutoff_s ) # smooth length is 1/2 max-gap

        return peaks

    cdef bool __close_peak_wo_subpeaks (self, list peak_content, peaks, int min_length,
                                          str chrom, int smoothlen=0, list score_cutoff_s=[]):
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peak_content contains [start, end, treat_p, ctrl_p, list_scores]

        peaks: a PeakIO object

        """
        cdef:
            int summit_pos, tstart, tend, tmpindex, summit_index, i, midindex
            double treat_v, ctrl_v, tsummitvalue, ttreat_p, tctrl_p, tscore, summit_treat, summit_ctrl, summit_p_score, summit_q_score
            list tlist_scores_p

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tsummit = []
            summit_pos   = 0
            summit_value = 0
            for i in range(len(peak_content)):
                (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
                tscore = ttreat_p #self.pqtable[ get_pscore(int(ttreat_p), tctrl_p) ] # use qscore as general score to find summit
                if not summit_value or summit_value < tscore:
                    tsummit = [(tend + tstart) / 2, ]
                    tsummit_index = [ i, ]
                    summit_value = tscore
                elif summit_value == tscore:
                    # remember continuous summit values
                    tsummit.append(int((tend + tstart) / 2))
                    tsummit_index.append( i )
            # the middle of all highest points in peak region is defined as summit
            midindex = int((len(tsummit) + 1) / 2) - 1
            summit_pos    = tsummit[ midindex ]
            summit_index  = tsummit_index[ midindex ]

            summit_treat = peak_content[ summit_index ][ 2 ]
            summit_ctrl = peak_content[ summit_index ][ 3 ]            
            
            # this is a double-check to see if the summit can pass cutoff values.
            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > peak_content[ summit_index ][ 4 ][i]:
                    return False # not passed, then disgard this peak.

            summit_p_score = get_pscore( int(summit_treat), summit_ctrl )
            summit_q_score = self.pqtable[ summit_p_score ]

            peaks.add( chrom,           # chromosome
                       peak_content[0][0], # start
                       peak_content[-1][1], # end
                       summit      = summit_pos, # summit position
                       peak_score  = summit_q_score, # score at summit
                       pileup      = summit_treat, # pileup
                       pscore      = summit_p_score, # pvalue
                       fold_change = float ( summit_treat + self.pseudocount ) / ( summit_ctrl + self.pseudocount ), # fold change
                       qscore      = summit_q_score # qvalue
                       )
            # start a new peak
            return True

    cdef bool __close_peak_with_subpeaks (self, list peak_content, peaks, int min_length,
                                         str chrom, int smoothlen=51, list score_cutoff_s=[],
                                         float min_valley = 0.9 ):
        """Algorithm implemented by Ben, to profile the pileup signals
        within a peak region then find subpeak summits. This method is
        highly recommended for TFBS or DNAase I sites.
        
        """
        cdef:
            int summit_pos, tstart, tend, tmpindex, summit_index, summit_offset
            int start, end, i, j, start_boundary, m, n
            double summit_value, tvalue, tsummitvalue
            np.ndarray[np.float32_t, ndim=1] peakdata
            np.ndarray[np.int32_t, ndim=1] peakindices, summit_offsets
            double ttreat_p, tctrl_p, tscore, summit_treat, summit_ctrl, summit_p_score, summit_q_score
            list tlist_scores_p
            
        # Add 10 bp padding to peak region so that we can get true minima
        end = peak_content[ -1 ][ 1 ] + 10
        start = peak_content[ 0 ][ 0 ] - 10
        if start < 0:
            start_boundary = 10 + start
            start = 0
        else:
            start_boundary = 10
        peak_length = end - start
        if end - start < min_length: return # if the region is too small, reject it

        peakdata = np.zeros(end - start, dtype='float32') # save the scores (qscore) for each position in this region
        peakindices = np.zeros(end - start, dtype='int32') # save the indices for each position in this region
        for i in range(len(peak_content)):
            (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
            #tscore = self.pqtable[ get_pscore(int(ttreat_p), tctrl_p) ] # use qscore as general score to find summit
            tscore = ttreat_p # use pileup as general score to find summit
            m = tstart - start + start_boundary
            n = tend - start + start_boundary
            peakdata[m:n] = tscore
            peakindices[m:n] = i

        summit_offsets = maxima(peakdata, smoothlen) # offsets are the indices for summits in peakdata/peakindices array.
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak_wo_subpeaks(peak_content, peaks, min_length, chrom)
        else:
            # remove maxima that occurred in padding
            m = np.searchsorted(summit_offsets, start_boundary)
            n = np.searchsorted(summit_offsets, peak_length + start_boundary, 'right')
            summit_offsets = summit_offsets[m:n]
        
        summit_offsets = enforce_peakyness(peakdata, summit_offsets)
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak_wo_subpeaks(peak_content, peaks, min_length, chrom)
        
        summit_indices = peakindices[summit_offsets] # indices are those point to peak_content
        summit_offsets -= start_boundary

        #peak_scores  = peakdata[ summit_offsets ]
        #if not (peak_scores > self.cutoff).all(): # fall back to old method
        #    return self.__close_peak_wo_subpeaks(peak_content, peaks, min_length, chrom)
        for summit_offset, summit_index in zip(summit_offsets, summit_indices):

            summit_treat = peak_content[ summit_index ][ 2 ]
            summit_ctrl = peak_content[ summit_index ][ 3 ]            

            summit_p_score = get_pscore( int(summit_treat), summit_ctrl )
            summit_q_score = self.pqtable[ summit_p_score ]

            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > peak_content[ summit_index ][ 4 ][i]:
                    return False # not passed, then disgard this summit.

            peaks.add( chrom,
                       start,
                       end,
                       summit      = start + summit_offset,
                       peak_score  = summit_q_score,
                       pileup      = summit_treat,
                       pscore      = summit_p_score,
                       fold_change = float ( summit_treat + self.pseudocount ) / ( summit_ctrl + self.pseudocount ), # fold change
                       qscore      = summit_q_score
                       )
        # start a new peak
        return True

    cdef np.ndarray __cal_pscore ( self, array1, array2 ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = get_pscore( int(array1[i]), array2[i] )
        return s

    cdef np.ndarray __cal_qscore ( self, array1, array2 ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = self.pqtable[ get_pscore( int(array1[i]), array2[i] ) ]
        return s

    cdef np.ndarray __cal_logLR ( self, array1, array2 ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = logLR( array1[i] + self.pseudocount, array2[i] + self.pseudocount ) 
        return s

    cdef np.ndarray __cal_logFE ( self, array1, array2 ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = get_logFE( array1[i] + self.pseudocount, array2[i] + self.pseudocount ) 
        return s

    cdef np.ndarray __cal_FE ( self, array1, array2 ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = (array1[i] + self.pseudocount) / ( array2[i] + self.pseudocount ) 
        return s

    cdef np.ndarray __cal_subtraction ( self, array1, array2 ):
        cdef:
            long i
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = array1[i] - array2[i]
        return s


    cdef bool __write_bedGraph_for_a_chromosome ( self, str chrom ):
        """Write treat/control values for a certain chromosome into a
        specified file handler.

        """
        cdef:
            np.ndarray pos_array, treat_array, ctrl_array
            int l, i
            int p, pre_p_t, pre_p_c # current position, previous position for treat, previous position for control
            float pre_v_t, pre_v_c, v_t, v_c # previous value for treat, for control, current value for treat, for control

        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        l = pos_array.shape[ 0 ]

        if l == 0:              # if there is no data, return
            return False

        # for treat
        t_write_func = self.bedGraph_treat.write
        c_write_func = self.bedGraph_ctrl.write
        
        pre_p_t = 0
        pre_p_c = 0
        pre_v_t = treat_array[ 0 ]
        pre_v_c = ctrl_array [ 0 ]
        
        for i in range( 1, l ):
            v_t = treat_array[ i ]
            v_c = ctrl_array [ i ]
            p   = pos_array  [ i-1 ]

            if abs(pre_v_t - v_t) > 1e-5: # precision is 5 digits
                t_write_func( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre_p_t, p, pre_v_t ) )
                pre_v_t = v_t
                pre_p_t = p

            if abs(pre_v_c - v_c) > 1e-5: # precision is 5 digits
                c_write_func( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre_p_c, p, pre_v_c ) )
                pre_v_c = v_c
                pre_p_c = p

        p = pos_array[ -1 ]
        # last one
        t_write_func( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre_p_t, p, pre_v_t ) )
        c_write_func( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre_p_c, p, pre_v_c ) )
            
        return True

    cpdef tuple call_broadpeaks (self, list scoring_function_symbols, list lvl1_cutoff_s, list lvl2_cutoff_s, int min_length=200, int lvl1_max_gap=50, int lvl2_max_gap=400):
        """This function try to find enriched regions within which,
        scores are continuously higher than a given cutoff for level
        1, and link them using the gap above level 2 cutoff with a
        maximum length of lvl2_max_gap.

        scoring_function_s: symbols of functions to calculate score. 'p' for pscore, 'q' for qscore, 'f' for fold change, 's' for subtraction. for example: ['p', 'q']

        lvl1_cutoff_s:  list of cutoffs at highly enriched regions, corresponding to scoring functions.
        lvl2_cutoff_s:  list of cutoffs at less enriched regions, corresponding to scoring functions.
        min_length :  minimum peak length, default 200.
        lvl1_max_gap   :  maximum gap to merge nearby enriched peaks, default 50.
        lvl2_max_gap   :  maximum length of linkage regions, default 400.        

        Return both general PeakIO object for highly enriched regions
        and gapped broad regions in BroadPeakIO.
        """
        cdef:
            int i
            str chrom
            object lvl1peaks, lvl1peakschrom, lvl1
            object lvl2peaks, lvl2peakschrom, lvl2
            object broadpeaks
            list chrs, tmppeakset

        lvl1peaks = PeakIO()
        lvl2peaks = PeakIO()

        # prepare p-q table
        if not self.pqtable:
            logging.info("#3 Pre-compute pvalue-qvalue table...")
            self.__cal_pvalue_qvalue_table()

        # prepare bedGraph file
        if self.save_bedGraph:

            self.bedGraph_treat = open( self.bedGraph_filename_prefix + "_treat_pileup.bdg", "w" )
            self.bedGraph_ctrl = open( self.bedGraph_filename_prefix + "_control_lambda.bdg", "w" )
            logging.info ("#3 In the peak calling step, the following will be performed simultaneously:")
            logging.info ("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... %s" % self.bedGraph_filename_prefix + "_treat_pileup.bdg")
            logging.info ("#3   Write bedGraph files for control lambda (after scaling if necessary)... %s" % self.bedGraph_filename_prefix + "_control_lambda.bdg")

            if self.trackline:
                # this line is REQUIRED by the wiggle format for UCSC browser
                self.bedGraph_treat.write( "track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix )
                self.bedGraph_ctrl.write ( "track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix )

        logging.info("#3 Call peaks for each chromosome...")
        for chrom in self.chromosomes:
            self.__chrom_call_broadpeak_using_certain_criteria ( lvl1peaks, lvl2peaks, chrom, scoring_function_symbols, lvl1_cutoff_s, lvl2_cutoff_s, min_length, lvl1_max_gap, lvl2_max_gap, self.save_bedGraph )

        # close bedGraph file
        if self.save_bedGraph:
            self.bedGraph_treat.close()
            self.bedGraph_ctrl.close()
            self.save_bedGraph = False

        # now combine lvl1 and lvl2 peaks
        chrs = lvl1peaks.get_chr_names()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in chrs:
            lvl1peakschrom = lvl1peaks.get_data_from_chrom(chrom)
            lvl2peakschrom = lvl2peaks.get_data_from_chrom(chrom)
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
        return lvl1peaks, broadpeaks

    cdef __chrom_call_broadpeak_using_certain_criteria ( self, lvl1peaks, lvl2peaks, str chrom, list scoring_function_s, list lvl1_cutoff_s, list lvl2_cutoff_s,
                                                         int min_length, int lvl1_max_gap, int lvl2_max_gap, bool save_bedGraph):
        """ Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        cdef:
            int i
            str s
            np.ndarray above_cutoff, above_cutoff_endpos, above_cutoff_startpos
            np.ndarray pos_array, treat_array, ctrl_array
            np.ndarray above_cutoff_index_array
            list score_array_s          # list to keep different types of scores
            list peak_content

        assert len(scoring_function_s) == len(lvl1_cutoff_s), "number of functions and cutoffs should be the same!"
        assert len(scoring_function_s) == len(lvl2_cutoff_s), "number of functions and cutoffs should be the same!"
        
        peak_content = []           # to store points above cutoff

        # first, build pileup, self.chr_pos_treat_ctrl
        self.__pileup_treat_ctrl_a_chromosome( chrom )
        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        # while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
        if save_bedGraph:
            self.__write_bedGraph_for_a_chromosome ( chrom )

        # keep all types of scores needed
        score_array_s = []
        for i in range(len(scoring_function_s)):
            s = scoring_function_s[i]
            if s == 'p':
                score_array_s.append( self.__cal_pscore( treat_array, ctrl_array ) )
            elif s == 'q':
                score_array_s.append( self.__cal_qscore( treat_array, ctrl_array ) )
            elif s == 'f':
                score_array_s.append( self.__cal_FE( treat_array, ctrl_array ) )
            elif s == 's':
                score_array_s.append( self.__cal_subtraction( treat_array, ctrl_array ) )

        # lvl1 : strong peaks
        # get the regions with scores above cutoffs
        above_cutoff = np.nonzero( apply_multiple_cutoffs(score_array_s,lvl1_cutoff_s) )[0] # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff_index_array = np.arange(pos_array.shape[0])[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return 

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], treat_array[above_cutoff_index_array[0]], ctrl_array[above_cutoff_index_array[0]]) )
        for i in range( 1,above_cutoff_startpos.size ):
            if above_cutoff_startpos[i] - peak_content[-1][1] <= lvl1_max_gap:
                # append
                peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]] ) )
            else:
                # close
                self.__close_peak_wo_subpeaks   (peak_content, lvl1peaks, min_length, chrom, lvl1_max_gap/2 )
                peak_content = [ (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]] ), ]
            
        # save the last peak
        if peak_content:
            self.__close_peak_wo_subpeaks   (peak_content, lvl1peaks, min_length, chrom, lvl1_max_gap/2 )            

        # lvl2 : weak peaks
        # get the regions with scores above cutoffs
        above_cutoff = np.nonzero( apply_multiple_cutoffs(score_array_s,lvl2_cutoff_s) )[0] # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff_index_array = np.arange(pos_array.shape[0])[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], treat_array[above_cutoff_index_array[0]], ctrl_array[above_cutoff_index_array[0]]) )
        for i in range( 1,above_cutoff_startpos.size ):
            if above_cutoff_startpos[i] - peak_content[-1][1] <= lvl2_max_gap:
                # append
                peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]] ) )
            else:
                # close
                self.__close_peak_wo_subpeaks   (peak_content, lvl2peaks, min_length, chrom, lvl2_max_gap/2 )
                peak_content = [ (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]] ), ]
            
        # save the last peak
        if peak_content:
            self.__close_peak_wo_subpeaks   (peak_content, lvl2peaks, min_length, chrom, lvl2_max_gap/2 )  

        return

    cdef __add_broadpeak (self, bpeaks, str chrom, object lvl2peak, list lvl1peakset):
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
        blockSizes = ",".join(map(str,map(itemgetter("length"),lvl1peakset))) #join( map(lambda x:str(x["length"]),lvl1peakset) )
        blockStarts = ",".join(getitem_then_subtract(lvl1peakset, start))     #join( map(lambda x:str(x["start"]-start),lvl1peakset) )
        #print blockSizes, blockStarts
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


cdef inline list getitem_then_subtract ( list peakset, int start ):
    cdef:
        list a
    
    a = map(itemgetter("start"), peakset)
    for i in range(len(a)):
        a[i] = str(a[i] - start)

    return a



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
