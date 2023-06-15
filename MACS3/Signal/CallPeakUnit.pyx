# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2022-09-15 17:06:17 Tao Liu>

"""Module for Calculate Scores.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

from collections import Counter
from copy import copy
from time import time as ttime
import _pickle as cPickle
from tempfile import mkstemp
import os

import logging
import MACS3.Utilities.Logger

logger = logging.getLogger(__name__)
debug   = logger.debug
info    = logger.info
# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
from cpython cimport bool
from cykhash import PyObjectMap, Float32to32Map

# ------------------------------------
# C lib
# ------------------------------------
from libc.stdio cimport *
from libc.math cimport exp,log,log10, M_LN10, log1p, erf, sqrt, floor, ceil

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.SignalProcessing import maxima, enforce_valleys, enforce_peakyness
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO, parse_peakname
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.Signal.Prob import poisson_cdf
# --------------------------------------------
# cached pscore function and LR_asym functions
# --------------------------------------------
pscore_dict = PyObjectMap()
logLR_dict = PyObjectMap()

cdef float32_t get_pscore ( tuple t ):
    """t: tuple of ( lambda, observation )
    """
    cdef:
        float32_t val
    if t in pscore_dict:
        return pscore_dict[ t ]
    else:
        # calculate and cache
        val = -1.0 * poisson_cdf ( t[0], t[1], False, True )
        pscore_dict[ t ] = val
        return val

cdef float32_t get_logLR_asym ( tuple t ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Set minus sign for depletion.
    """
    cdef:
        float32_t val
        float32_t x
        float32_t y
    if t in logLR_dict:
        return logLR_dict[ t ]
    else:
        x = t[0]
        y = t[1]
        # calculate and cache
        if x > y:
            val = (x*(log10(x)-log10(y))+y-x)
        elif x < y:
            val = (x*(-log10(x)+log10(y))-y+x)
        else:
            val = 0
        logLR_dict[ t ] = val
        return val

# ------------------------------------
# constants
# ------------------------------------
__version__ = "CallPeakUnit $Revision$"
__author__ = "Tao Liu <vladimir.liu@gmail.com>"
__doc__ = "CallPeakUnit"

LOG10_E = 0.43429448190325176

# ------------------------------------
# Misc functions
# ------------------------------------

cdef void clean_up_ndarray ( np.ndarray x ):
    # clean numpy ndarray in two steps
    cdef:
        int64_t i
    i = x.shape[0] // 2
    x.resize( 100000 if i > 100000 else i, refcheck=False)
    x.resize( 0, refcheck=False)
    return

cdef inline float32_t chi2_k1_cdf ( float32_t x ):
    return erf( sqrt(x/2) )

cdef inline float32_t log10_chi2_k1_cdf ( float32_t x ):
  return log10( erf( sqrt(x/2) ) )

cdef inline float32_t chi2_k2_cdf ( float32_t x ):
  return 1 - exp( -x/2 )

cdef inline float32_t log10_chi2_k2_cdf ( float32_t x ):
  return log1p( - exp( -x/2 ) ) * LOG10_E

cdef inline float32_t chi2_k4_cdf ( float32_t x ):
  return 1 - exp( -x/2 ) * ( 1 + x/2 )

cdef inline float32_t log10_chi2_k4_CDF ( float32_t x ):
  return log1p( - exp( -x/2 ) * ( 1 + x/2 ) ) * LOG10_E

cdef inline np.ndarray apply_multiple_cutoffs ( list multiple_score_arrays, list multiple_cutoffs ):
    cdef:
        int32_t i
        np.ndarray ret

    ret = multiple_score_arrays[0] > multiple_cutoffs[0]

    for i in range(1,len(multiple_score_arrays)):
        ret += multiple_score_arrays[i] > multiple_cutoffs[i]

    return ret

cdef inline list get_from_multiple_scores ( list multiple_score_arrays, int32_t index ):
    cdef:
        list ret = []
        int32_t i

    for i in range(len(multiple_score_arrays)):
        ret.append(multiple_score_arrays[i][index])
    return ret


cdef inline float32_t get_logFE ( float32_t x, float32_t y ):
    """ return 100* log10 fold enrichment with +1 pseudocount.
    """
    return log10( x/y )

cdef inline float32_t get_subtraction ( float32_t x, float32_t y):
    """ return subtraction.
    """
    return x - y

cdef inline list getitem_then_subtract ( list peakset, int32_t start ):
    cdef:
        list a

    a = [x["start"] for x in peakset]
    for i in range(len(a)):
        a[i] = a[i] - start
    return a

cdef inline int32_t left_sum ( data, int32_t pos, int32_t width ):
    """
    """
    return sum([data[x] for x in data if x <= pos and x >= pos - width])

cdef inline int32_t right_sum ( data, int32_t pos, int32_t width ):
    """
    """
    return sum([data[x] for x in data if x >= pos and x <= pos + width])

cdef inline int32_t left_forward ( data, int32_t pos, int32_t window_size ):
    return data.get(pos,0) - data.get(pos-window_size, 0)

cdef inline int32_t right_forward ( data, int32_t pos, int32_t window_size ):
    return data.get(pos + window_size, 0) - data.get(pos, 0)

cdef float32_t median_from_value_length ( np.ndarray[np.float32_t, ndim=1] value, list length ):
    """
    """
    cdef:
        list tmp
        int32_t c, tmp_l
        float32_t tmp_v, mid_l

    c = 0
    tmp = sorted(list(zip( value, length )))
    mid_l = sum( length )/2
    for (tmp_v, tmp_l) in tmp:
        c += tmp_l
        if c > mid_l:
            return tmp_v

cdef float32_t mean_from_value_length ( np.ndarray[np.float32_t, ndim=1] value, list length ):
    """take list of values and list of corresponding lengths, calculate the mean.
    An important function for bedGraph type of data.
    """
    cdef:
        int32_t i
        int32_t tmp_l, l
        float64_t tmp_v, sum_v, tmp_sum   #try to solve precision issue
        float32_t ret

    sum_v = 0
    l = 0

    for i in range( len(length) ):
        tmp_l = length[ i ]
        tmp_v = <float64_t>value[ i ]
        tmp_sum = tmp_v * tmp_l
        sum_v = tmp_sum + sum_v
        l += tmp_l

    ret = <float32_t>(sum_v/l)

    return ret


cdef tuple find_optimal_cutoff( list x, list y ):
    """Return the best cutoff x and y.

    We assume that total peak length increase exponentially while
    decreasing cutoff value. But while cutoff decreases to a point
    that background noises are captured, total length increases much
    faster. So we fit a linear model by taking the first 10 points,
    then look for the largest cutoff that

    """
    cdef:
        np.ndarray npx, npy, npA
        float32_t optimal_x, optimal_y
        int64_t l, i
        float32_t m, c # slop and intercept
        float32_t sst # sum of squared total
        float32_t sse # sum of squared error
        float32_t rsq # R-squared

    l = len(x)
    assert l == len(y)
    npx = np.array( x )
    npy = np.log10( np.array( y ) )
    npA = np.vstack( [npx, np.ones(len(npx))] ).T

    for i in range( 10, l ):
        # at least the largest 10 points
        m, c = np.linalg.lstsq( npA[:i], npy[:i], rcond=None )[ 0 ]
        sst = sum( ( npy[:i] - np.mean( npy[:i] ) ) ** 2 )
        sse = sum( ( npy[:i] - m*npx[:i] - c ) ** 2 )
        rsq = 1 - sse/sst
        #print i, x[i], y[i], m, c, rsq
    return ( 1.0, 1.0 )



# ------------------------------------
# Classes
# ------------------------------------
cdef class CallerFromAlignments:
    """A unit to calculate scores and call peaks from alignments --
    FWTrack or PETrack objects.

    It will compute for each chromosome separately in order to save
    memory usage.
    """
    cdef:
        object treat            # FWTrack or PETrackI object for ChIP
        object ctrl             # FWTrack or PETrackI object for Control

        int32_t  d                           # extension size for ChIP
        list ctrl_d_s                    # extension sizes for Control. Can be multiple values
        float32_t treat_scaling_factor       # scaling factor for ChIP
        list ctrl_scaling_factor_s       # scaling factor for Control, corresponding to each extension size.
        float32_t lambda_bg                  # minimum local bias to fill missing values
        list chromosomes                 # name of common chromosomes in ChIP and Control data
        float64_t pseudocount                # the pseudocount used to calcuate logLR, FE or logFE
        bytes bedGraph_filename_prefix     # prefix will be added to _pileup.bdg for treatment and _lambda.bdg for control

        int32_t  end_shift                   # shift of cutting ends before extension
        bool trackline                   # whether trackline should be saved in bedGraph
        bool save_bedGraph               # whether to save pileup and local bias in bedGraph files
        bool save_SPMR                   # whether to save pileup normalized by sequencing depth in million reads
        bool no_lambda_flag              # whether ignore local bias, and to use global bias instead
        bool PE_mode                     # whether it's in PE mode, will be detected during initiation

        # temporary data buffer
        list chr_pos_treat_ctrl          # temporary [position, treat_pileup, ctrl_pileup] for a given chromosome
        bytes bedGraph_treat_filename
        bytes bedGraph_control_filename
        FILE * bedGraph_treat_f
        FILE * bedGraph_ctrl_f

        # data needed to be pre-computed before peak calling
        object pqtable          # remember pvalue->qvalue convertion; saved in cykhash Float32to32Map
        bool pvalue_all_done             # whether the pvalue of whole genome is all calculated. If yes, it's OK to calculate q-value.

        dict pvalue_npeaks               # record for each pvalue cutoff, how many peaks can be called
        dict pvalue_length               # record for each pvalue cutoff, the total length of called peaks
        float32_t optimal_p_cutoff           # automatically decide the p-value cutoff ( can be translated into qvalue cutoff ) based
                                         # on p-value to total peak length analysis.
        bytes cutoff_analysis_filename     # file to save the pvalue-npeaks-totallength table

        dict pileup_data_files           # Record the names of temporary files for storing pileup values of each chromosome


    def __init__ (self, treat, ctrl,
                  int32_t d = 200, list ctrl_d_s = [200, 1000, 10000],
                  float32_t treat_scaling_factor = 1.0, list ctrl_scaling_factor_s = [1.0, 0.2, 0.02],
                  bool stderr_on = False,
                  float32_t pseudocount = 1,
                  int32_t end_shift = 0,
                  float32_t lambda_bg = 0,
                  bool save_bedGraph = False,
                  str  bedGraph_filename_prefix = "PREFIX",
                  str bedGraph_treat_filename = "TREAT.bdg",
                  str bedGraph_control_filename = "CTRL.bdg",
                  str cutoff_analysis_filename = "TMP.txt",
                  bool save_SPMR = False ):
        """Initialize.

        A calculator is unique to each comparison of treat and
        control. Treat_depth and ctrl_depth should not be changed
        during calculation.

        treat and ctrl are either FWTrack or PETrackI objects.

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
            int32_t i
            char * tmp
            bytes tmp_bytes
            float32_t p
        # decide PE mode
        if isinstance(treat, FWTrack):
            self.PE_mode = False
        elif isinstance(treat, PETrackI):
            self.PE_mode = True
        else:
            raise Exception("Should be FWTrack or PETrackI object!")
        # decide if there is control
        self.treat = treat
        if ctrl:
            self.ctrl = ctrl
        else:                   # while there is no control
            self.ctrl = treat
        self.trackline = False
        self.d = d              # note, self.d doesn't make sense in PE mode
        self.ctrl_d_s = ctrl_d_s# note, self.d doesn't make sense in PE mode
        self.treat_scaling_factor = treat_scaling_factor
        self.ctrl_scaling_factor_s= ctrl_scaling_factor_s
        self.end_shift = end_shift
        self.lambda_bg = lambda_bg
        self.pqtable = Float32to32Map( for_int = False ) # Float32 -> Float32 map
        self.save_bedGraph = save_bedGraph
        self.save_SPMR = save_SPMR
        self.bedGraph_filename_prefix = bedGraph_filename_prefix.encode()
        self.bedGraph_treat_filename = bedGraph_treat_filename.encode()
        self.bedGraph_control_filename = bedGraph_control_filename.encode()
        if not self.ctrl_d_s or not self.ctrl_scaling_factor_s:
            self.no_lambda_flag = True
        else:
            self.no_lambda_flag = False
        self.pseudocount = pseudocount
        # get the common chromosome names from both treatment and control
        chr1 = set(self.treat.get_chr_names())
        chr2 = set(self.ctrl.get_chr_names())
        self.chromosomes = sorted(list(chr1.intersection(chr2)))

        self.pileup_data_files = {}
        self.pvalue_length = {}
        self.pvalue_npeaks = {}
        for p in np.arange( 0.3, 10, 0.3 ): # step for optimal cutoff is 0.3 in -log10pvalue, we try from pvalue 1E-10 (-10logp=10) to 0.5 (-10logp=0.3)
            self.pvalue_length[ p ] = 0
            self.pvalue_npeaks[ p ] = 0
        self.optimal_p_cutoff = 0
        self.cutoff_analysis_filename = cutoff_analysis_filename.encode()

    cpdef destroy ( self ):
        """Remove temporary files for pileup values of each chromosome.

        Note: This function MUST be called if the class object won't
        be used anymore.

        """
        cdef:
            bytes f

        for f in self.pileup_data_files.values():
            if os.path.isfile( f ):
                os.unlink( f )
        return

    cpdef set_pseudocount( self, float32_t pseudocount ):
        self.pseudocount = pseudocount

    cpdef enable_trackline( self ):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True

    cdef __pileup_treat_ctrl_a_chromosome ( self, bytes chrom ):
        """After this function is called, self.chr_pos_treat_ctrl will
        be reset and assigned to the pileup values of the given
        chromosome.

        """
        cdef:
            list treat_pv, ctrl_pv
            int64_t i
            float32_t t
            object f
            str temp_filename

        assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom

        # check backup file of pileup values. If not exists, create
        # it. Otherwise, load them instead of calculating new pileup
        # values.
        if chrom in self.pileup_data_files:
            try:
                f = open( self.pileup_data_files[ chrom ],"rb" )
                self.chr_pos_treat_ctrl = cPickle.load( f )
                f.close()
                return
            except:
                temp_fd, temp_filename = mkstemp()
                os.close(temp_fd)
                self.pileup_data_files[ chrom ] = temp_filename
        else:
            temp_fd, temp_filename = mkstemp()
            os.close(temp_fd)
            self.pileup_data_files[ chrom ] = temp_filename.encode()

        # reset or clean existing self.chr_pos_treat_ctrl
        if self.chr_pos_treat_ctrl:     # not a beautiful way to clean
            clean_up_ndarray( self.chr_pos_treat_ctrl[0] )
            clean_up_ndarray( self.chr_pos_treat_ctrl[1] )
            clean_up_ndarray( self.chr_pos_treat_ctrl[2] )

        if self.PE_mode:
            treat_pv = self.treat.pileup_a_chromosome ( chrom, [self.treat_scaling_factor,], baseline_value = 0.0 )
        else:
            treat_pv = self.treat.pileup_a_chromosome( chrom, [self.d,], [self.treat_scaling_factor,], baseline_value = 0.0,
                                                       directional = True,
                                                       end_shift = self.end_shift )

        if not self.no_lambda_flag:
            if self.PE_mode:
                # note, we pileup up PE control as SE control because
                # we assume the bias only can be captured at the
                # surrounding regions of cutting sites from control experiments.
                ctrl_pv = self.ctrl.pileup_a_chromosome_c( chrom, self.ctrl_d_s, self.ctrl_scaling_factor_s, baseline_value = self.lambda_bg )
            else:
                ctrl_pv = self.ctrl.pileup_a_chromosome( chrom, self.ctrl_d_s, self.ctrl_scaling_factor_s,
                                                         baseline_value = self.lambda_bg,
                                                         directional = False )
        else:
            ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg,], dtype="float32")] # set a global lambda

        self.chr_pos_treat_ctrl = self.__chrom_pair_treat_ctrl( treat_pv, ctrl_pv)

        # clean treat_pv and ctrl_pv
        treat_pv = []
        ctrl_pv  = []

        # save data to temporary file
        try:
            f = open(self.pileup_data_files[ chrom ],"wb")
            cPickle.dump( self.chr_pos_treat_ctrl, f , protocol=2 )
            f.close()
        except:
            # fail to write then remove the key in pileup_data_files
            self.pileup_data_files.pop(chrom)
        return

    cdef list __chrom_pair_treat_ctrl ( self, treat_pv, ctrl_pv ):
        """*private* Pair treat and ctrl pileup for each region.

        treat_pv and ctrl_pv are [np.ndarray, np.ndarray].

        return [p, t, c] list, each element is a numpy array.
        """
        cdef:
            list ret
            int64_t pre_p, index_ret, it, ic, lt, lc
            np.ndarray[np.int32_t, ndim=1] t_p, c_p, ret_p
            np.ndarray[np.float32_t, ndim=1] t_v, c_v, ret_t, ret_c

            int32_t * t_p_ptr
            int32_t * c_p_ptr
            int32_t * ret_p_ptr

            float32_t * t_v_ptr
            float32_t * c_v_ptr
            float32_t * ret_t_ptr
            float32_t * ret_c_ptr

        [ t_p, t_v ] = treat_pv
        [ c_p, c_v ] = ctrl_pv

        lt = t_p.shape[0]
        lc = c_p.shape[0]

        chrom_max_len = lt + lc

        ret_p = np.zeros( chrom_max_len, dtype="int32" ) # position
        ret_t = np.zeros( chrom_max_len, dtype="float32" ) # value from treatment
        ret_c = np.zeros( chrom_max_len, dtype="float32" ) # value from control

        t_p_ptr = <int32_t *> t_p.data
        t_v_ptr = <float32_t *> t_v.data
        c_p_ptr = <int32_t *> c_p.data
        c_v_ptr = <float32_t *> c_v.data
        ret_p_ptr = <int32_t *> ret_p.data
        ret_t_ptr = <float32_t *>ret_t.data
        ret_c_ptr = <float32_t *>ret_c.data

        pre_p = 0
        index_ret = 0
        it = 0
        ic = 0

        while it < lt and ic < lc:
            if t_p_ptr[0] < c_p_ptr[0]:
                # clip a region from pre_p to p1, then set pre_p as p1.
                ret_p_ptr[0] = t_p_ptr[0]
                ret_t_ptr[0] = t_v_ptr[0]
                ret_c_ptr[0] = c_v_ptr[0]
                ret_p_ptr += 1
                ret_t_ptr += 1
                ret_c_ptr += 1
                pre_p = t_p_ptr[0]
                index_ret += 1
                # call for the next p1 and v1
                it += 1
                t_p_ptr += 1
                t_v_ptr += 1
            elif t_p_ptr[0] > c_p_ptr[0]:
                # clip a region from pre_p to p2, then set pre_p as p2.
                ret_p_ptr[0] = c_p_ptr[0]
                ret_t_ptr[0] = t_v_ptr[0]
                ret_c_ptr[0] = c_v_ptr[0]
                ret_p_ptr += 1
                ret_t_ptr += 1
                ret_c_ptr += 1
                pre_p = c_p_ptr[0]
                index_ret += 1
                # call for the next p2 and v2
                ic += 1
                c_p_ptr += 1
                c_v_ptr += 1
            else:
                # from pre_p to p1 or p2, then set pre_p as p1 or p2.
                ret_p_ptr[0] = t_p_ptr[0]
                ret_t_ptr[0] = t_v_ptr[0]
                ret_c_ptr[0] = c_v_ptr[0]
                ret_p_ptr += 1
                ret_t_ptr += 1
                ret_c_ptr += 1
                pre_p = t_p_ptr[0]
                index_ret += 1
                # call for the next p1, v1, p2, v2.
                it += 1
                ic += 1
                t_p_ptr += 1
                t_v_ptr += 1
                c_p_ptr += 1
                c_v_ptr += 1

        ret_p.resize( index_ret, refcheck=False)
        ret_t.resize( index_ret, refcheck=False)
        ret_c.resize( index_ret, refcheck=False)
        return [ret_p, ret_t, ret_c]

    cdef np.ndarray __cal_score ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2, cal_func ):
        cdef:
            int64_t i
            np.ndarray[np.float32_t, ndim=1] s
        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")
        for i in range(array1.shape[0]):
            s[i] = cal_func( array1[i], array2[i] )
        return s

    cdef void __cal_pvalue_qvalue_table ( self ):
        """After this function is called, self.pqtable is built. All
        chromosomes will be iterated. So it will take some time.

        """
        cdef:
            bytes chrom
            np.ndarray pos_array, treat_array, ctrl_array, score_array
            dict pscore_stat
            int64_t n, pre_p, length, pre_l, l, i, j
            float32_t this_v, pre_v, v, q, pre_q
            int64_t N, k, this_l
            float32_t f
            list unique_values
            int32_t * pos_ptr
            float32_t * treat_value_ptr
            float32_t * ctrl_value_ptr

        debug ( "Start to calculate pvalue stat..." )

        pscore_stat = {} #dict()
        for i in range( len( self.chromosomes ) ):
            chrom = self.chromosomes[ i ]
            pre_p = 0

            self.__pileup_treat_ctrl_a_chromosome( chrom )
            [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

            pos_ptr = <int32_t *> pos_array.data
            treat_value_ptr = <float32_t *> treat_array.data
            ctrl_value_ptr = <float32_t *> ctrl_array.data

            for j in range(pos_array.shape[0]):
                this_v = get_pscore( (<int32_t>(treat_value_ptr[0]), ctrl_value_ptr[0] ) )
                this_l = pos_ptr[0] - pre_p
                if this_v in pscore_stat:
                    pscore_stat[ this_v ] += this_l
                else:
                    pscore_stat[ this_v ] = this_l
                pre_p = pos_ptr[0]
                pos_ptr += 1
                treat_value_ptr += 1
                ctrl_value_ptr += 1

        N = sum(pscore_stat.values()) # total length
        k = 1                         # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647      # save the previous q-value

        self.pqtable = Float32to32Map( for_int = False )
        unique_values = sorted(list(pscore_stat.keys()), reverse=True) 
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = pscore_stat[v]
            q = v + (log10(k) + f)
            if q > pre_q:
                q = pre_q
            if q <= 0:
                q = 0
                break
            #q = max(0,min(pre_q,q))           # make q-score monotonic
            self.pqtable[ v ] = q
            pre_q = q
            k += l
        # bottom rank pscores all have qscores 0
        for j in range(i, len(unique_values) ):
            v = unique_values[ j ]
            self.pqtable[ v ] = 0
        return

    cdef void __pre_computes ( self, int32_t max_gap = 50, int32_t min_length = 200 ):
        """After this function is called, self.pqtable and self.pvalue_length is built. All
        chromosomes will be iterated. So it will take some time.

        """
        cdef:
            bytes chrom
            np.ndarray pos_array, treat_array, ctrl_array, score_array
            dict pscore_stat
            int64_t n, pre_p, this_p, length, j, pre_l, l, i
            float32_t q, pre_q, this_t, this_c
            float32_t this_v, pre_v, v, cutoff
            int64_t N, k, this_l
            float32_t f
            list unique_values
            float64_t t0, t1, t

            np.ndarray above_cutoff, above_cutoff_endpos, above_cutoff_startpos
            list peak_content
            int64_t peak_length, total_l, total_p

            list tmplist

            int32_t * acs_ptr   # above cutoff start position pointer
            int32_t * ace_ptr   # above cutoff end position pointer
            int32_t * pos_array_ptr # position array pointer
            float32_t * score_array_ptr # score array pointer

        debug ( "Start to calculate pvalue stat..." )

        # tmplist contains a list of log pvalue cutoffs from 0.3 to 10
        tmplist = [round(x,5) for x in sorted( list(np.arange(0.3, 10.0, 0.3)), reverse = True )]

        pscore_stat = {} #dict()
        #print (list(pscore_stat.keys()))
        #print (list(self.pvalue_length.keys()))
        #print (list(self.pvalue_npeaks.keys()))
        for i in range( len( self.chromosomes ) ):
            chrom = self.chromosomes[ i ]
            self.__pileup_treat_ctrl_a_chromosome( chrom )
            [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

            score_array = self.__cal_pscore( treat_array, ctrl_array )

            for n in range( len( tmplist ) ):
                cutoff = tmplist[ n ]
                total_l = 0           # total length in potential peak
                total_p = 0

                # get the regions with scores above cutoffs
                above_cutoff = np.nonzero( score_array > cutoff )[0]# this is not an optimized method. It would be better to store score array in a 2-D ndarray?
                above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
                above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

                if above_cutoff_endpos.size == 0:
                    continue

                # first bit of region above cutoff
                acs_ptr = <int32_t *> above_cutoff_startpos.data
                ace_ptr = <int32_t *> above_cutoff_endpos.data

                peak_content = [( acs_ptr[ 0 ], ace_ptr[ 0 ] ), ]
                lastp = ace_ptr[ 0 ]
                acs_ptr += 1
                ace_ptr += 1

                for i in range( 1, above_cutoff_startpos.size ):
                    tl = acs_ptr[ 0 ] - lastp
                    if tl <= max_gap:
                        peak_content.append( ( acs_ptr[ 0 ], ace_ptr[ 0 ] ) )
                    else:
                        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                        if peak_length >= min_length: # if the peak is too small, reject it
                            total_l +=  peak_length
                            total_p += 1
                        peak_content = [ ( acs_ptr[ 0 ], ace_ptr[ 0 ] ), ]
                    lastp = ace_ptr[ 0 ]
                    acs_ptr += 1
                    ace_ptr += 1

                if peak_content:
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        total_l +=  peak_length
                        total_p += 1
                self.pvalue_length[ cutoff ] = self.pvalue_length.get( cutoff, 0 ) + total_l
                self.pvalue_npeaks[ cutoff ] = self.pvalue_npeaks.get( cutoff, 0 ) + total_p

            pos_array_ptr = <int32_t *> pos_array.data
            score_array_ptr = <float32_t *> score_array.data

            pre_p = 0
            for i in range(pos_array.shape[0]):
                this_p = pos_array_ptr[ 0 ]
                this_l = this_p - pre_p
                this_v = score_array_ptr[ 0 ]
                if this_v in pscore_stat:
                    pscore_stat[ this_v ] += this_l
                else:
                    pscore_stat[ this_v ] = this_l
                pre_p = this_p #pos_array[ i ]
                pos_array_ptr += 1
                score_array_ptr += 1

        #debug ( "make pscore_stat cost %.5f seconds" % t )

        # add all pvalue cutoffs from cutoff-analysis part. So that we
        # can get the corresponding qvalues for them.
        for cutoff in tmplist:
            if cutoff not in pscore_stat:
                pscore_stat[ cutoff ] = 0

        nhval = 0

        N = sum(pscore_stat.values()) # total length
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value

        self.pqtable = Float32to32Map( for_int = False ) #{}
        unique_values = sorted(list(pscore_stat.keys()), reverse=True) #sorted(unique_values,reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = pscore_stat[v]
            q = v + (log10(k) + f)
            if q > pre_q:
                q = pre_q
            if q <= 0:
                q = 0
                break
            #q = max(0,min(pre_q,q))           # make q-score monotonic
            self.pqtable[ v ] = q
            pre_v = v
            pre_q = q
            k+=l
        for j in range(i, len(unique_values) ):
            v = unique_values[ j ]
            self.pqtable[ v ] = 0

        # write pvalue and total length of predicted peaks
        # this is the output from cutoff-analysis
        fhd = open( self.cutoff_analysis_filename, "w" )
        fhd.write( "pscore\tqscore\tnpeaks\tlpeaks\tavelpeak\n" )
        x = []
        y = []
        for cutoff in tmplist:
            if self.pvalue_npeaks[ cutoff ] > 0:
                fhd.write( "%.2f\t%.2f\t%d\t%d\t%.2f\n" % ( cutoff, self.pqtable[ cutoff ], self.pvalue_npeaks[ cutoff ], self.pvalue_length[ cutoff ], self.pvalue_length[ cutoff ]/self.pvalue_npeaks[ cutoff ] ) )
                x.append( cutoff )
                y.append( self.pvalue_length[ cutoff ] )
        fhd.close()
        info( "#3 Analysis of cutoff vs num of peaks or total length has been saved in %s" % self.cutoff_analysis_filename )
        #info( "#3 Suggest a cutoff..." )
        #optimal_cutoff, optimal_length = find_optimal_cutoff( x, y )
        #info( "#3 -10log10pvalue cutoff %.2f will call approximately %.0f bps regions as significant regions" % ( optimal_cutoff, optimal_length ) )
        #print (list(pqtable.keys()))
        #print (list(self.pvalue_length.keys()))
        #print (list(self.pvalue_npeaks.keys()))
        return

    cpdef call_peaks ( self, list scoring_function_symbols, list score_cutoff_s, int32_t min_length = 200,
                       int32_t max_gap = 50, bool call_summits = False, bool cutoff_analysis = False ):
        """Call peaks for all chromosomes. Return a PeakIO object.

        scoring_function_s: symbols of functions to calculate score. 'p' for pscore, 'q' for qscore, 'f' for fold change, 's' for subtraction. for example: ['p', 'q']
        score_cutoff_s    : cutoff values corresponding to scoring functions
        min_length        : minimum length of peak
        max_gap           : maximum gap of 'insignificant' regions within a peak. Note, for PE_mode, max_gap and max_length are both set as fragment length.
        call_summits      : boolean. Whether or not call sub-peaks.
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        cdef:
            bytes chrom
            bytes tmp_bytes

        peaks = PeakIO()

        # prepare p-q table
        if len( self.pqtable ) == 0:
            info("#3 Pre-compute pvalue-qvalue table...")
            if cutoff_analysis:
                info("#3 Cutoff vs peaks called will be analyzed!")
                self.__pre_computes( max_gap = max_gap, min_length = min_length )
            else:
                self.__cal_pvalue_qvalue_table()


        # prepare bedGraph file
        if self.save_bedGraph:
            self.bedGraph_treat_f = fopen( self.bedGraph_treat_filename, "w" )
            self.bedGraph_ctrl_f = fopen( self.bedGraph_control_filename, "w" )

            info ("#3 In the peak calling step, the following will be performed simultaneously:")
            info ("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... %s" % self.bedGraph_filename_prefix.decode() + "_treat_pileup.bdg")
            info ("#3   Write bedGraph files for control lambda (after scaling if necessary)... %s" % self.bedGraph_filename_prefix.decode() + "_control_lambda.bdg")

            if self.save_SPMR:
                info ( "#3   --SPMR is requested, so pileup will be normalized by sequencing depth in million reads." )
            elif self.treat_scaling_factor == 1:
                info ( "#3   Pileup will be based on sequencing depth in treatment." )
            else:
                info ( "#3   Pileup will be based on sequencing depth in control." )

            if self.trackline:
                # this line is REQUIRED by the wiggle format for UCSC browser
                tmp_bytes = ("track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf( self.bedGraph_treat_f, tmp_bytes )
                tmp_bytes = ("track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf( self.bedGraph_ctrl_f, tmp_bytes )

        info("#3 Call peaks for each chromosome...")
        for chrom in self.chromosomes:
            # treat/control bedGraph will be saved if requested by user.
            self.__chrom_call_peak_using_certain_criteria ( peaks, chrom, scoring_function_symbols, score_cutoff_s, min_length, max_gap, call_summits, self.save_bedGraph )

        # close bedGraph file
        if self.save_bedGraph:
            fclose(self.bedGraph_treat_f)
            fclose(self.bedGraph_ctrl_f)
            self.save_bedGraph = False

        return peaks

    cdef void __chrom_call_peak_using_certain_criteria ( self, peaks, bytes chrom, list scoring_function_s, list score_cutoff_s, int32_t min_length,
                                                   int32_t max_gap, bool call_summits, bool save_bedGraph ):
        """ Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object, the return value of this function
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        cdef:
            float64_t t0
            int32_t i, n
            str s
            np.ndarray above_cutoff
            np.ndarray[np.int32_t, ndim=1] above_cutoff_endpos, above_cutoff_startpos, pos_array, above_cutoff_index_array

            np.ndarray[np.float32_t, ndim=1] treat_array, ctrl_array
            list score_array_s          # list to keep different types of scores
            list peak_content           #  to store information for a
                                        #  chunk in a peak region, it
                                        #  contains lists of: 1. left
                                        #  position; 2. right
                                        #  position; 3. treatment
                                        #  value; 4. control value;
                                        #  5. list of scores at this
                                        #  chunk
            int64_t tl, lastp, ts, te, ti
            float32_t tp, cp
            int32_t * acs_ptr
            int32_t * ace_ptr
            int32_t * acia_ptr
            float32_t * treat_array_ptr
            float32_t * ctrl_array_ptr


        assert len(scoring_function_s) == len(score_cutoff_s), "number of functions and cutoffs should be the same!"

        peak_content = []           # to store points above cutoff

        # first, build pileup, self.chr_pos_treat_ctrl
        # this step will be speeped up if pqtable is pre-computed.
        self.__pileup_treat_ctrl_a_chromosome( chrom )
        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        # while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
        if save_bedGraph:
            self.__write_bedGraph_for_a_chromosome ( chrom )

        # keep all types of scores needed
        #t0 = ttime()
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
        above_cutoff_index_array = np.arange(pos_array.shape[0],dtype="int32")[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        #print "apply cutoff -- chrom:",chrom,"  time:", ttime() - t0
        # start to build peak regions
        #t0 = ttime()

        # first bit of region above cutoff
        acs_ptr = <int32_t *>above_cutoff_startpos.data
        ace_ptr = <int32_t *>above_cutoff_endpos.data
        acia_ptr= <int32_t *>above_cutoff_index_array.data
        treat_array_ptr = <float32_t *> treat_array.data
        ctrl_array_ptr = <float32_t *> ctrl_array.data

        ts = acs_ptr[ 0 ]
        te = ace_ptr[ 0 ]
        ti = acia_ptr[ 0 ]
        tp = treat_array_ptr[ ti ]
        cp = ctrl_array_ptr[ ti ]

        peak_content.append( ( ts, te, tp, cp, ti ) )
        lastp = te
        acs_ptr += 1
        ace_ptr += 1
        acia_ptr+= 1

        for i in range( 1, above_cutoff_startpos.shape[0] ):
            ts = acs_ptr[ 0 ]
            te = ace_ptr[ 0 ]
            ti = acia_ptr[ 0 ]
            acs_ptr += 1
            ace_ptr += 1
            acia_ptr+= 1
            tp = treat_array_ptr[ ti ]
            cp = ctrl_array_ptr[ ti ]
            tl = ts - lastp
            if tl <= max_gap:
                # append.
                peak_content.append( ( ts, te, tp, cp, ti ) )
                lastp = te #above_cutoff_endpos[i]
            else:
                # close
                if call_summits:
                    self.__close_peak_with_subpeaks (peak_content, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s = score_cutoff_s ) # smooth length is min_length, i.e. fragment size 'd'
                else:
                    self.__close_peak_wo_subpeaks   (peak_content, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s = score_cutoff_s ) # smooth length is min_length, i.e. fragment size 'd'
                peak_content = [ ( ts, te, tp, cp, ti ), ]
                lastp = te #above_cutoff_endpos[i]
        # save the last peak
        if not peak_content:
            return
        else:
            if call_summits:
                self.__close_peak_with_subpeaks (peak_content, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s = score_cutoff_s ) # smooth length is min_length, i.e. fragment size 'd'
            else:
                self.__close_peak_wo_subpeaks   (peak_content, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s = score_cutoff_s ) # smooth length is min_length, i.e. fragment size 'd'

        #print "close peaks -- chrom:",chrom,"  time:", ttime() - t0
        return

    cdef bool __close_peak_wo_subpeaks (self, list peak_content, peaks, int32_t min_length,
                                          bytes chrom, int32_t smoothlen, list score_array_s, list score_cutoff_s=[]):
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peak_content contains [start, end, treat_p, ctrl_p, index_in_score_array]

        peaks: a PeakIO object

        """
        cdef:
            int32_t summit_pos, tstart, tend, tmpindex, summit_index, i, midindex
            float64_t treat_v, ctrl_v, tsummitvalue, ttreat_p, tctrl_p, tscore, summit_treat, summit_ctrl, summit_p_score, summit_q_score
            int32_t tlist_scores_p

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tsummit = []
            summit_pos   = 0
            summit_value = 0
            for i in range(len(peak_content)):
                (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
                tscore = ttreat_p # use pscore as general score to find summit
                if not summit_value or summit_value < tscore:
                    tsummit = [(tend + tstart) // 2, ]
                    tsummit_index = [ i, ]
                    summit_value = tscore
                elif summit_value == tscore:
                    # remember continuous summit values
                    tsummit.append((tend + tstart) // 2)
                    tsummit_index.append( i )
            # the middle of all highest points in peak region is defined as summit
            midindex = (len(tsummit) + 1) // 2 - 1
            summit_pos    = tsummit[ midindex ]
            summit_index  = tsummit_index[ midindex ]

            summit_treat = peak_content[ summit_index ][ 2 ]
            summit_ctrl = peak_content[ summit_index ][ 3 ]

            # this is a double-check to see if the summit can pass cutoff values.
            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > score_array_s[ i ][ peak_content[ summit_index ][ 4 ] ]:
                    return False # not passed, then disgard this peak.

            summit_p_score = pscore_dict[ ( <int32_t>(summit_treat), summit_ctrl ) ] #get_pscore(( <int32_t>(summit_treat), summit_ctrl ) )
            summit_q_score = self.pqtable[ summit_p_score ]

            peaks.add( chrom,           # chromosome
                       peak_content[0][0], # start
                       peak_content[-1][1], # end
                       summit      = summit_pos, # summit position
                       peak_score  = summit_q_score, # score at summit
                       pileup      = summit_treat, # pileup
                       pscore      = summit_p_score, # pvalue
                       fold_change = ( summit_treat + self.pseudocount ) / ( summit_ctrl + self.pseudocount ), # fold change
                       qscore      = summit_q_score # qvalue
                       )
            # start a new peak
            return True

    cdef bool __close_peak_with_subpeaks (self, list peak_content, peaks, int32_t min_length,
                                         bytes chrom, int32_t smoothlen, list score_array_s, list score_cutoff_s=[],
                                         float32_t min_valley = 0.9 ):
        """Algorithm implemented by Ben, to profile the pileup signals
        within a peak region then find subpeak summits. This method is
        highly recommended for TFBS or DNAase I sites.

        """
        cdef:
            int32_t summit_pos, tstart, tend, tmpindex, summit_index, summit_offset
            int32_t start, end, i, j, start_boundary, m, n, l
            float64_t summit_value, tvalue, tsummitvalue, ttreat_p, tctrl_p, tscore, summit_treat, summit_ctrl, summit_p_score, summit_q_score
            np.ndarray[np.float32_t, ndim=1] peakdata
            np.ndarray[np.int32_t, ndim=1] peakindices, summit_offsets
            int32_t tlist_scores_p

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]

        if peak_length < min_length: return  # if the region is too small, reject it

        # Add 10 bp padding to peak region so that we can get true minima
        end = peak_content[ -1 ][ 1 ] + 10
        start = peak_content[ 0 ][ 0 ] - 10
        if start < 0:
            start_boundary = 10 + start # this is the offset of original peak boundary in peakdata list.
            start = 0
        else:
            start_boundary = 10 # this is the offset of original peak boundary in peakdata list.

        peakdata = np.zeros(end - start, dtype='float32') # save the scores (qscore) for each position in this region
        peakindices = np.zeros(end - start, dtype='int32') # save the indices for each position in this region
        for i in range(len(peak_content)):
            (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
            tscore = ttreat_p # use pileup as general score to find summit
            m = tstart - start + start_boundary
            n = tend - start + start_boundary
            peakdata[m:n] = tscore
            peakindices[m:n] = i

        summit_offsets = maxima(peakdata, smoothlen) # offsets are the indices for summits in peakdata/peakindices array.

        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak_wo_subpeaks(peak_content, peaks, min_length, chrom, smoothlen, score_array_s, score_cutoff_s)
        else:
            # remove maxima that occurred in padding
            m = np.searchsorted(summit_offsets, start_boundary)
            n = np.searchsorted(summit_offsets, peak_length + start_boundary, 'right')
            summit_offsets = summit_offsets[m:n]

        summit_offsets = enforce_peakyness(peakdata, summit_offsets)

        #print "enforced:",summit_offsets
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak_wo_subpeaks(peak_content, peaks, min_length, chrom, smoothlen, score_array_s, score_cutoff_s)

        summit_indices = peakindices[summit_offsets] # indices are those point to peak_content
        summit_offsets -= start_boundary

        for summit_offset, summit_index in list(zip(summit_offsets, summit_indices)):

            summit_treat = peak_content[ summit_index ][ 2 ]
            summit_ctrl = peak_content[ summit_index ][ 3 ]

            summit_p_score = pscore_dict[ ( <int32_t>(summit_treat), summit_ctrl ) ] # get_pscore(( <int32_t>(summit_treat), summit_ctrl ) )
            summit_q_score = self.pqtable[ summit_p_score ]

            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > score_array_s[ i ][ peak_content[ summit_index ][ 4 ] ]:
                    return False # not passed, then disgard this summit.

            peaks.add( chrom,
                       peak_content[ 0 ][ 0 ],
                       peak_content[ -1 ][ 1 ],
                       summit      = start + summit_offset,
                       peak_score  = summit_q_score,
                       pileup      = summit_treat,
                       pscore      = summit_p_score,
                       fold_change = (summit_treat + self.pseudocount ) / ( summit_ctrl + self.pseudocount ), # fold change
                       qscore      = summit_q_score
                       )
        # start a new peak
        return True

    cdef np.ndarray __cal_pscore ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2 ):
        cdef:
            int64_t i, array1_size
            np.ndarray[np.float32_t, ndim=1] s
            float32_t * a1_ptr
            float32_t * a2_ptr
            float32_t * s_ptr

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")

        a1_ptr = <float32_t *> array1.data
        a2_ptr = <float32_t *> array2.data
        s_ptr = <float32_t *> s.data

        array1_size = array1.shape[0]

        for i in range(array1_size):
            s_ptr[0] = get_pscore(( <int32_t>(a1_ptr[0]), a2_ptr[0] ))
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    cdef np.ndarray __cal_qscore ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2 ):
        cdef:
            int64_t i, array1_size
            np.ndarray[np.float32_t, ndim=1] s
            float32_t * a1_ptr
            float32_t * a2_ptr
            float32_t * s_ptr

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")

        a1_ptr = <float32_t *> array1.data
        a2_ptr = <float32_t *> array2.data
        s_ptr = <float32_t *> s.data

        for i in range(array1.shape[0]):
            s_ptr[0] = self.pqtable[ get_pscore(( <int32_t>(a1_ptr[0]), a2_ptr[0] )) ]
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    cdef np.ndarray __cal_logLR ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2 ):
        cdef:
            int64_t i, array1_size
            np.ndarray[np.float32_t, ndim=1] s
            float32_t * a1_ptr
            float32_t * a2_ptr
            float32_t * s_ptr

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")

        a1_ptr = <float32_t *> array1.data
        a2_ptr = <float32_t *> array2.data
        s_ptr = <float32_t *> s.data

        for i in range(array1.shape[0]):
            s_ptr[0] = get_logLR_asym( (a1_ptr[0] + self.pseudocount, a2_ptr[0] + self.pseudocount ) )
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    cdef np.ndarray __cal_logFE ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2 ):
        cdef:
            int64_t i, array1_size
            np.ndarray[np.float32_t, ndim=1] s
            float32_t * a1_ptr
            float32_t * a2_ptr
            float32_t * s_ptr

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")

        a1_ptr = <float32_t *> array1.data
        a2_ptr = <float32_t *> array2.data
        s_ptr = <float32_t *> s.data

        for i in range(array1.shape[0]):
            s_ptr[0] = get_logFE( a1_ptr[0] + self.pseudocount, a2_ptr[0] + self.pseudocount )
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    cdef np.ndarray __cal_FE ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2 ):
        cdef:
            int64_t i, array1_size
            np.ndarray[np.float32_t, ndim=1] s
            float32_t * a1_ptr
            float32_t * a2_ptr
            float32_t * s_ptr

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")

        a1_ptr = <float32_t *> array1.data
        a2_ptr = <float32_t *> array2.data
        s_ptr = <float32_t *> s.data

        for i in range(array1.shape[0]):
            s_ptr[0] = (a1_ptr[0] + self.pseudocount) / ( a2_ptr[0] + self.pseudocount )
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    cdef np.ndarray __cal_subtraction ( self, np.ndarray[np.float32_t, ndim=1] array1, np.ndarray[np.float32_t, ndim=1] array2 ):
        cdef:
            int64_t i, array1_size
            np.ndarray[np.float32_t, ndim=1] s
            float32_t * a1_ptr
            float32_t * a2_ptr
            float32_t * s_ptr

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="float32")

        a1_ptr = <float32_t *> array1.data
        a2_ptr = <float32_t *> array2.data
        s_ptr = <float32_t *> s.data

        for i in range(array1.shape[0]):
            s_ptr[0] = a1_ptr[0] - a2_ptr[0]
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s


    cdef bool __write_bedGraph_for_a_chromosome ( self, bytes chrom ):
        """Write treat/control values for a certain chromosome into a
        specified file handler.

        """
        cdef:
            np.ndarray[np.int32_t, ndim=1] pos_array
            np.ndarray[np.float32_t, ndim=1] treat_array, ctrl_array
            int32_t * pos_array_ptr
            float32_t * treat_array_ptr
            float32_t * ctrl_array_ptr
            int32_t l, i
            int32_t p, pre_p_t, pre_p_c # current position, previous position for treat, previous position for control
            float32_t pre_v_t, pre_v_c, v_t, v_c # previous value for treat, for control, current value for treat, for control
            float32_t denominator # 1 if save_SPMR is false, or depth in million if save_SPMR is true. Note, while piling up and calling peaks, treatment and control have been scaled to the same depth, so we need to find what this 'depth' is.
            FILE * ft
            FILE * fc
            basestring tmp_bytes

        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl
        pos_array_ptr = <int32_t *> pos_array.data
        treat_array_ptr = <float32_t *> treat_array.data
        ctrl_array_ptr = <float32_t *> ctrl_array.data

        if self.save_SPMR:
            if self.treat_scaling_factor == 1:
                # in this case, control has been asked to be scaled to depth of treatment
                denominator  = self.treat.total/1e6
            else:
                # in this case, treatment has been asked to be scaled to depth of control
                denominator  = self.ctrl.total/1e6
        else:
            denominator = 1.0

        l = pos_array.shape[ 0 ]

        if l == 0:              # if there is no data, return
            return False

        ft = self.bedGraph_treat_f
        fc = self.bedGraph_ctrl_f
        #t_write_func = self.bedGraph_treat.write
        #c_write_func = self.bedGraph_ctrl.write

        pre_p_t = 0
        pre_p_c = 0
        pre_v_t = treat_array_ptr[ 0 ]/denominator
        pre_v_c = ctrl_array_ptr [ 0 ]/denominator
        treat_array_ptr += 1
        ctrl_array_ptr += 1

        for i in range( 1, l ):
            v_t = treat_array_ptr[ 0 ]/denominator
            v_c = ctrl_array_ptr [ 0 ]/denominator
            p   = pos_array_ptr  [ 0 ]
            pos_array_ptr += 1
            treat_array_ptr += 1
            ctrl_array_ptr += 1

            if abs(pre_v_t - v_t) > 1e-5: # precision is 5 digits
                fprintf( ft, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_t, p, pre_v_t )
                pre_v_t = v_t
                pre_p_t = p

            if abs(pre_v_c - v_c) > 1e-5: # precision is 5 digits
                fprintf( fc, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_c, p, pre_v_c )
                pre_v_c = v_c
                pre_p_c = p

        p = pos_array_ptr[ 0 ]
        # last one
        fprintf( ft, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_t, p, pre_v_t )
        fprintf( fc, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_c, p, pre_v_c )

        return True

    cpdef call_broadpeaks (self, list scoring_function_symbols, list lvl1_cutoff_s, list lvl2_cutoff_s, int32_t min_length=200, int32_t lvl1_max_gap=50, int32_t lvl2_max_gap=400, bool cutoff_analysis = False):
        """This function try to find enriched regions within which,
        scores are continuously higher than a given cutoff for level
        1, and link them using the gap above level 2 cutoff with a
        maximum length of lvl2_max_gap.

        scoring_function_s: symbols of functions to calculate score. 'p' for pscore, 'q' for qscore, 'f' for fold change, 's' for subtraction. for example: ['p', 'q']

        lvl1_cutoff_s:  list of cutoffs at highly enriched regions, corresponding to scoring functions.
        lvl2_cutoff_s:  list of cutoffs at less enriched regions, corresponding to scoring functions.
        min_length :    minimum peak length, default 200.
        lvl1_max_gap   :  maximum gap to merge nearby enriched peaks, default 50.
        lvl2_max_gap   :  maximum length of linkage regions, default 400.

        Return both general PeakIO object for highly enriched regions
        and gapped broad regions in BroadPeakIO.
        """
        cdef:
            int32_t i, j
            bytes chrom
            object lvl1peaks, lvl1peakschrom, lvl1
            object lvl2peaks, lvl2peakschrom, lvl2
            object broadpeaks
            set chrs
            list tmppeakset

        lvl1peaks = PeakIO()
        lvl2peaks = PeakIO()

        # prepare p-q table
        if len( self.pqtable ) == 0:
            info("#3 Pre-compute pvalue-qvalue table...")
            if cutoff_analysis:
                info("#3 Cutoff value vs broad region calls will be analyzed!")
                self.__pre_computes( max_gap = lvl2_max_gap, min_length = min_length )
            else:
                self.__cal_pvalue_qvalue_table()

        # prepare bedGraph file
        if self.save_bedGraph:

            self.bedGraph_treat_f = fopen( self.bedGraph_treat_filename, "w" )
            self.bedGraph_ctrl_f = fopen( self.bedGraph_control_filename, "w" )
            info ("#3 In the peak calling step, the following will be performed simultaneously:")
            info ("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... %s" % self.bedGraph_filename_prefix.decode() + "_treat_pileup.bdg")
            info ("#3   Write bedGraph files for control lambda (after scaling if necessary)... %s" % self.bedGraph_filename_prefix.decode() + "_control_lambda.bdg")

            if self.trackline:
                # this line is REQUIRED by the wiggle format for UCSC browser
                tmp_bytes = ("track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf( self.bedGraph_treat_f, tmp_bytes )
                tmp_bytes = ("track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf( self.bedGraph_ctrl_f, tmp_bytes )


        info("#3 Call peaks for each chromosome...")
        for chrom in self.chromosomes:
            self.__chrom_call_broadpeak_using_certain_criteria ( lvl1peaks, lvl2peaks, chrom, scoring_function_symbols, lvl1_cutoff_s, lvl2_cutoff_s, min_length, lvl1_max_gap, lvl2_max_gap, self.save_bedGraph )

        # close bedGraph file
        if self.save_bedGraph:
            fclose( self.bedGraph_treat_f )
            fclose( self.bedGraph_ctrl_f )
            #self.bedGraph_ctrl.close()
            self.save_bedGraph = False

        # now combine lvl1 and lvl2 peaks
        chrs = lvl1peaks.get_chr_names()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in sorted(chrs):
            lvl1peakschrom = lvl1peaks.get_data_from_chrom(chrom)
            lvl2peakschrom = lvl2peaks.get_data_from_chrom(chrom)
            lvl1peakschrom_next = iter(lvl1peakschrom).__next__
            tmppeakset = []             # to temporarily store lvl1 region inside a lvl2 region
            # our assumption is lvl1 regions should be included in lvl2 regions
            try:
                lvl1 = lvl1peakschrom_next()
                for i in range( len(lvl2peakschrom) ):
                    # for each lvl2 peak, find all lvl1 peaks inside
                    # I assume lvl1 peaks can be ALL covered by lvl2 peaks.
                    lvl2 = lvl2peakschrom[i]

                    while True:
                        if lvl2["start"] <= lvl1["start"]  and lvl1["end"] <= lvl2["end"]:
                            tmppeakset.append(lvl1)
                            lvl1 = lvl1peakschrom_next()
                        else:
                            # make a hierarchical broad peak
                            #print lvl2["start"], lvl2["end"], lvl2["score"]
                            self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
                            tmppeakset = []
                            break
            except StopIteration:
                # no more strong (aka lvl1) peaks left
                self.__add_broadpeak ( broadpeaks, chrom, lvl2, tmppeakset)
                tmppeakset = []
                # add the rest lvl2 peaks
                for j in range( i+1, len(lvl2peakschrom) ):
                    self.__add_broadpeak( broadpeaks, chrom, lvl2peakschrom[j], tmppeakset )

        return broadpeaks

    cdef void  __chrom_call_broadpeak_using_certain_criteria ( self, lvl1peaks, lvl2peaks, bytes chrom, list scoring_function_s, list lvl1_cutoff_s, list lvl2_cutoff_s,
                                                         int32_t min_length, int32_t lvl1_max_gap, int32_t lvl2_max_gap, bool save_bedGraph):
        """ Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        cdef:
            int32_t i
            str s
            np.ndarray above_cutoff, above_cutoff_endpos, above_cutoff_startpos
            np.ndarray pos_array, treat_array, ctrl_array
            np.ndarray above_cutoff_index_array
            list score_array_s          # list to keep different types of scores
            list peak_content
            int32_t * acs_ptr
            int32_t * ace_ptr
            int32_t * acia_ptr
            float32_t * treat_array_ptr
            float32_t * ctrl_array_ptr

        assert len(scoring_function_s) == len(lvl1_cutoff_s), "number of functions and cutoffs should be the same!"
        assert len(scoring_function_s) == len(lvl2_cutoff_s), "number of functions and cutoffs should be the same!"

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
        peak_content = []           # to store points above cutoff

        # get the regions with scores above cutoffs
        above_cutoff = np.nonzero( apply_multiple_cutoffs(score_array_s,lvl1_cutoff_s) )[0] # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff_index_array = np.arange(pos_array.shape[0],dtype="int32")[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        acs_ptr = <int32_t *>above_cutoff_startpos.data
        ace_ptr = <int32_t *>above_cutoff_endpos.data
        acia_ptr= <int32_t *>above_cutoff_index_array.data
        treat_array_ptr = <float32_t *> treat_array.data
        ctrl_array_ptr = <float32_t *> ctrl_array.data

        ts = acs_ptr[ 0 ]
        te = ace_ptr[ 0 ]
        ti = acia_ptr[ 0 ]
        tp = treat_array_ptr[ ti ]
        cp = ctrl_array_ptr[ ti ]

        peak_content.append( ( ts, te, tp, cp, ti ) )
        acs_ptr += 1 # move ptr
        ace_ptr += 1
        acia_ptr+= 1
        lastp = te

        #peak_content.append( (above_cutoff_startpos[0], above_cutoff_endpos[0], treat_array[above_cutoff_index_array[0]], ctrl_array[above_cutoff_index_array[0]], score_array_s, above_cutoff_index_array[0] ) )
        for i in range( 1, above_cutoff_startpos.size ):
            ts = acs_ptr[ 0 ]
            te = ace_ptr[ 0 ]
            ti = acia_ptr[ 0 ]
            acs_ptr += 1
            ace_ptr += 1
            acia_ptr+= 1
            tp = treat_array_ptr[ ti ]
            cp = ctrl_array_ptr[ ti ]
            tl = ts - lastp
            if tl <= lvl1_max_gap:
                # append
                #peak_content.append( (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]], score_array_s, above_cutoff_index_array[i] ) )
                peak_content.append( ( ts, te, tp, cp, ti ) )
                lastp = te
            else:
                # close
                self.__close_peak_for_broad_region (peak_content, lvl1peaks, min_length, chrom, lvl1_max_gap//2, score_array_s )
                #peak_content = [ (above_cutoff_startpos[i], above_cutoff_endpos[i], treat_array[above_cutoff_index_array[i]], ctrl_array[above_cutoff_index_array[i]], score_array_s, above_cutoff_index_array[i]) , ]
                peak_content = [ ( ts, te, tp, cp, ti ), ]
                lastp = te #above_cutoff_endpos[i]

        # save the last peak
        if peak_content:
            self.__close_peak_for_broad_region (peak_content, lvl1peaks, min_length, chrom, lvl1_max_gap//2, score_array_s )

        # lvl2 : weak peaks
        peak_content = []           # to store points above cutoff

        # get the regions with scores above cutoffs
        above_cutoff = np.nonzero( apply_multiple_cutoffs(score_array_s,lvl2_cutoff_s) )[0] # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff_index_array = np.arange(pos_array.shape[0],dtype="int32")[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff] # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1] # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        acs_ptr = <int32_t *>above_cutoff_startpos.data
        ace_ptr = <int32_t *>above_cutoff_endpos.data
        acia_ptr= <int32_t *>above_cutoff_index_array.data
        treat_array_ptr = <float32_t *> treat_array.data
        ctrl_array_ptr = <float32_t *> ctrl_array.data

        ts = acs_ptr[ 0 ]
        te = ace_ptr[ 0 ]
        ti = acia_ptr[ 0 ]
        tp = treat_array_ptr[ ti ]
        cp = ctrl_array_ptr[ ti ]
        peak_content.append( ( ts, te, tp, cp, ti ) )
        acs_ptr += 1 # move ptr
        ace_ptr += 1
        acia_ptr+= 1

        lastp = te
        for i in range( 1, above_cutoff_startpos.size ):
            # for everything above cutoff
            ts = acs_ptr[ 0 ] # get the start
            te = ace_ptr[ 0 ] # get the end
            ti = acia_ptr[ 0 ]# get the index

            acs_ptr += 1 # move ptr
            ace_ptr += 1
            acia_ptr+= 1
            tp = treat_array_ptr[ ti ] # get the treatment pileup
            cp = ctrl_array_ptr[ ti ]  # get the control pileup
            tl = ts - lastp # get the distance from the current point to last position of existing peak_content

            if tl <= lvl2_max_gap:
                # append
                peak_content.append( ( ts, te, tp, cp, ti ) )
                lastp = te
            else:
                # close
                self.__close_peak_for_broad_region (peak_content, lvl2peaks, min_length, chrom, lvl2_max_gap//2, score_array_s )

                peak_content = [ ( ts, te, tp, cp, ti ), ]
                lastp = te

        # save the last peak
        if peak_content:
            self.__close_peak_for_broad_region (peak_content, lvl2peaks, min_length, chrom, lvl2_max_gap//2, score_array_s )

        return

    cdef bool __close_peak_for_broad_region (self, list peak_content, peaks, int32_t min_length,
                                             bytes chrom, int32_t smoothlen, list score_array_s, list score_cutoff_s=[]):
        """Close the broad peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peak_content contains [start, end, treat_p, ctrl_p, list_scores]

        peaks: a BroadPeakIO object

        """
        cdef:
            int32_t summit_pos, tstart, tend, tmpindex, summit_index, i, midindex
            float64_t treat_v, ctrl_v, tsummitvalue, ttreat_p, tctrl_p, tscore, summit_treat, summit_ctrl, summit_p_score, summit_q_score
            list tlist_pileup, tlist_control, tlist_length
            int32_t tlist_scores_p
            np.ndarray tarray_pileup, tarray_control, tarray_pscore, tarray_qscore, tarray_fc

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tlist_pileup = []
            tlist_control= []
            tlist_length = []
            for i in range(len(peak_content)): # each position in broad peak
                (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
                tlist_pileup.append( ttreat_p )
                tlist_control.append( tctrl_p )
                tlist_length.append( tend - tstart )

            tarray_pileup = np.array( tlist_pileup, dtype="float32")
            tarray_control = np.array( tlist_control, dtype="float32")
            tarray_pscore = self.__cal_pscore( tarray_pileup, tarray_control )
            tarray_qscore = self.__cal_qscore( tarray_pileup, tarray_control )
            tarray_fc     = self.__cal_FE    ( tarray_pileup, tarray_control )

            peaks.add( chrom,           # chromosome
                       peak_content[0][0], # start
                       peak_content[-1][1], # end
                       summit = 0,
                       peak_score  = mean_from_value_length( tarray_qscore, tlist_length ),
                       pileup      = mean_from_value_length( tarray_pileup, tlist_length ),
                       pscore      = mean_from_value_length( tarray_pscore, tlist_length ),
                       fold_change = mean_from_value_length( tarray_fc, tlist_length ),
                       qscore      = mean_from_value_length( tarray_qscore, tlist_length ),
                       )
            #if chrom == "chr1" and  peak_content[0][0] == 237643 and peak_content[-1][1] == 237935:
            #    print tarray_qscore, tlist_length
            # start a new peak
            return True

    cdef __add_broadpeak (self, bpeaks, bytes chrom, object lvl2peak, list lvl1peakset):
        """Internal function to create broad peak.

        *Note* lvl1peakset/strong_regions might be empty
        """

        cdef:
            int32_t blockNum, start, end
            bytes blockSizes, blockStarts, thickStart, thickEnd,

        start      = lvl2peak["start"]
        end        = lvl2peak["end"]

        if not lvl1peakset:
            # will complement by adding 1bps start and end to this region
            # may change in the future if gappedPeak format was improved.
            bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=(b"%d" % start), thickEnd=(b"%d" % end),
                       blockNum = 2, blockSizes = b"1,1", blockStarts = (b"0,%d" % (end-start-1)), pileup = lvl2peak["pileup"],
                       pscore = lvl2peak["pscore"], fold_change = lvl2peak["fc"],
                       qscore = lvl2peak["qscore"] )
            return bpeaks

        thickStart = b"%d" % (lvl1peakset[0]["start"])
        thickEnd   = b"%d" % (lvl1peakset[-1]["end"])
        blockNum   = len(lvl1peakset)
        blockSizes = b",".join([b"%d" % y for y in [x["length"] for x in lvl1peakset]])
        blockStarts = b",".join([b"%d" % x for x in getitem_then_subtract(lvl1peakset, start)])

        # add 1bp left and/or right block if necessary
        if int(thickStart) != start:
            # add 1bp left block
            thickStart = b"%d" % start
            blockNum += 1
            blockSizes = b"1,"+blockSizes
            blockStarts = b"0,"+blockStarts
        if int(thickEnd) != end:
            # add 1bp right block
            thickEnd = b"%d" % end
            blockNum += 1
            blockSizes = blockSizes + b",1"
            blockStarts = blockStarts + b"," + (b"%d" % (end-start-1))

        bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=thickStart, thickEnd=thickEnd,
                   blockNum = blockNum, blockSizes = blockSizes, blockStarts = blockStarts, pileup = lvl2peak["pileup"],
                   pscore = lvl2peak["pscore"], fold_change = lvl2peak["fc"],
                   qscore = lvl2peak["qscore"] )
        return bpeaks

 
