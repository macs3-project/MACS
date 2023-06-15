# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-09-15 17:24:53 Tao Liu>

"""Module for Feature IO classes.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from copy import copy
from functools import reduce

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.SignalProcessing import maxima, enforce_valleys, enforce_peakyness
from MACS3.Signal.Prob import poisson_cdf
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO, parse_peakname

# ------------------------------------
# Other modules
# ------------------------------------
cimport cython
import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
from cpython cimport bool
from cykhash import PyObjectMap, Float32to32Map

# ------------------------------------
# C lib
# ------------------------------------
from libc.math cimport log10,log, floor, ceil

# ------------------------------------
# constants
# ------------------------------------
__version__ = "scoreTrack $Revision$"
__author__ = "Tao Liu <vladimir.liu@gmail.com>"
__doc__ = "scoreTrack classes"

# ------------------------------------
# Misc functions
# ------------------------------------
cdef inline int32_t int_max(int32_t a, int32_t b): return a if a >= b else b
cdef inline int32_t int_min(int32_t a, int32_t b): return a if a <= b else b

LOG10_E = 0.43429448190325176

pscore_dict = PyObjectMap()

cdef float32_t get_pscore ( int32_t observed, float32_t expectation ):
    """Get p-value score from Poisson test. First check existing
    table, if failed, call poisson_cdf function, then store the result
    in table.

    """
    cdef:
        float64_t score

    try:
        return pscore_dict[(observed, expectation)]
    except KeyError:
        score = -1*poisson_cdf(observed,expectation,False,True)
        pscore_dict[(observed, expectation)] = score
        return score

asym_logLR_dict = PyObjectMap()

cdef float32_t logLR_asym ( float32_t x, float32_t y ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    chromatin bias ). Set minus sign for depletion.

    *asymmetric version*

    """
    cdef:
        float32_t s

    if (x,y) in asym_logLR_dict:
        return asym_logLR_dict[ ( x, y ) ]
    else:
        if x > y:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
        elif x < y:
            s = (x*(-log(x)+log(y))-y+x)*LOG10_E
        else:
            s = 0
        asym_logLR_dict[ ( x, y ) ] = s
        return s

sym_logLR_dict = PyObjectMap()

cdef float32_t logLR_sym ( float32_t x, float32_t y ):
    """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
    another enriched ). Set minus sign for H0>H1.

    * symmetric version *

    """
    cdef:
        float32_t s

    if (x,y) in sym_logLR_dict:
        return sym_logLR_dict[ ( x, y ) ]
    else:
        if x > y:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
        elif y > x:
            s = (y*(log(x)-log(y))+y-x)*LOG10_E
        else:
            s = 0
        sym_logLR_dict[ ( x, y ) ] = s
        return s

cdef float32_t get_logFE ( float32_t x, float32_t y ):
    """ return 100* log10 fold enrichment with +1 pseudocount.
    """
    return log10( x/y )

cdef float32_t get_subtraction ( float32_t x, float32_t y):
    """ return subtraction.
    """
    return x - y

# ------------------------------------
# Classes
# ------------------------------------

cdef class ScoreTrackII:
    """Class for a container to keep signals of each genomic position,
    including 1. score, 2. treatment and 2. control pileup.

    It also contains scoring methods and call_peak functions.
    """
    cdef:
        dict data                       # dictionary for data of each chromosome
        dict datalength                 # length of data array of each chromosome
        bool trackline                  # whether trackline should be saved in bedGraph
        float32_t treat_edm             # seq depth in million of treatment
        float32_t ctrl_edm              # seq depth in million of control
        char scoring_method             # method for calculating scores.
        char normalization_method       # scale to control? scale to treatment? both scale to 1million reads?
        float32_t pseudocount           # the pseudocount used to calcuate logLR, FE or logFE
        float32_t cutoff
        dict pvalue_stat                # save pvalue<->length dictionary


    def __init__ (self, float32_t treat_depth, float32_t ctrl_depth, float32_t pseudocount = 1.0 ):
        """Initialize.

        treat_depth and ctrl_depth are effective depth in million:
                                    sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        pseudocount: a pseudocount used to calculate logLR, FE or
                     logFE. Please note this value will not be changed
                     with normalization method. So if you really want
                     to set pseudocount 1 per million reads, set it
                     after you normalize treat and control by million
                     reads by `change_normalizetion_method(ord('M'))`.

        """
        self.data = {}           # for each chromosome, there is a l*4
                                 # matrix. First column: end position
                                 # of a region; Second: treatment
                                 # pileup; third: control pileup ;
                                 # forth: score ( can be
                                 # p/q-value/likelihood
                                 # ratio/fold-enrichment/subtraction
                                 # depending on -c setting)
        self.datalength = {}
        self.trackline = False
        self.treat_edm = treat_depth
        self.ctrl_edm = ctrl_depth
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
        self.pvalue_stat = {}

    cpdef set_pseudocount( self, float32_t pseudocount ):
        self.pseudocount = pseudocount

    cpdef enable_trackline( self ):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True

    cpdef add_chromosome ( self, bytes chrom, int32_t chrom_max_len ):
        """
        chrom: chromosome name
        chrom_max_len: maximum number of data points in this chromosome

        """
        if chrom not in self.data:
            #self.data[chrom] = np.zeros( ( chrom_max_len, 4 ), dtype="int32" ) # remember col #2-4 is actual value * 100, I use integer here.
            self.data[chrom] = [ np.zeros( chrom_max_len, dtype="int32" ), # pos
                                 np.zeros( chrom_max_len, dtype="float32" ), # pileup at each interval, in float32 format
                                 np.zeros( chrom_max_len, dtype="float32" ), # control at each interval, in float32 format
                                 np.zeros( chrom_max_len, dtype="float32" ) ] # score at each interval, in float32 format
            self.datalength[chrom] = 0

    cpdef add (self, bytes chromosome, int32_t endpos, float32_t chip, float32_t control):
        """Add a chr-endpos-sample-control block into data
        dictionary.

        chromosome: chromosome name in string
        endpos    : end position of each interval in integer
        chip      : ChIP pileup value of each interval in float
        control   : Control pileup value of each interval in float

        *Warning* Need to add regions continuously.
        """
        cdef int32_t i
        i = self.datalength[chromosome]
        c = self.data[chromosome]
        c[0][ i ] = endpos
        c[1][ i ] = chip
        c[2][ i ] = control
        self.datalength[chromosome] += 1

    cpdef finalize ( self ):
        """
        Adjust array size of each chromosome.

        """
        cdef:
            bytes chrom, k
            int32_t l

        for chrom in sorted(self.data.keys()):
            d = self.data[chrom]
            l = self.datalength[chrom]
            d[0].resize( l, refcheck = False )
            d[1].resize( l, refcheck = False )
            d[2].resize( l, refcheck = False )
            d[3].resize( l, refcheck = False )
        return

    cpdef get_data_by_chr (self, bytes chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if chromosome in self.data:
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
        if normalization_method == ord('T'):
            if self.normalization_method == ord('T'): # do nothing
                pass
            elif self.normalization_method == ord('C'):
                self.normalize( self.treat_edm/self.ctrl_edm, self.treat_edm/self.ctrl_edm )
            elif  self.normalization_method == ord('M'):
                self.normalize( self.treat_edm, self.treat_edm )
            elif self.normalization_method == ord('N'):
                self.normalize( 1, self.treat_edm/self.ctrl_edm )
            else:
                raise NotImplemented
            self.normalization_method = ord('T')
        elif normalization_method == ord('C'):
            if self.normalization_method == ord('T'):
                self.normalize( self.ctrl_edm/self.treat_edm, self.ctrl_edm/self.treat_edm )
            elif self.normalization_method == ord('C'): # do nothing
                pass
            elif  self.normalization_method == ord('M'):
                self.normalize( self.ctrl_edm, self.ctrl_edm )
            elif self.normalization_method == ord('N'):
                self.normalize( self.ctrl_edm/self.treat_edm, 1 )
            else:
                raise NotImplemented
            self.normalization_method = ord('C')
        elif normalization_method == ord('M'):
            if self.normalization_method == ord('T'):
                self.normalize( 1/self.treat_edm, 1/self.treat_edm )
            elif self.normalization_method == ord('C'):
                self.normalize( 1/self.ctrl_edm, 1/self.ctrl_edm )
            elif  self.normalization_method == ord('M'): # do nothing
                pass
            elif self.normalization_method == ord('N'):
                self.normalize( 1/self.treat_edm, 1/self.ctrl_edm )
            else:
                raise NotImplemented
            self.normalization_method = ord('M')
        elif normalization_method == ord('N'):
            if self.normalization_method == ord('T'):
                self.normalize( self.treat_edm, self.treat_edm )
            elif self.normalization_method == ord('C'):
                self.normalize( self.ctrl_edm, self.ctrl_edm )
            elif  self.normalization_method == ord('M'):
                self.normalize( self.treat_edm, self.ctrl_edm )
            elif self.normalization_method == ord('N'): # do nothing
                pass
            else:
                raise NotImplemented
            self.normalization_method = ord('N')

    cdef normalize ( self, float32_t treat_scale, float32_t control_scale ):
        cdef:
            np.ndarray p, c
            int64_t l, i

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            l = self.datalength[chrom]
            for i in range(l):
                p[ i ] *= treat_scale
                c[ i ] *= control_scale
        return

    cpdef change_score_method (self, char scoring_method):
        """
        scoring_method:  p: -log10 pvalue;
                         q: -log10 qvalue;
                         l: log10 likelihood ratio ( minus for depletion )
			 s: symmetric log10 likelihood ratio ( for comparing two ChIPs )
                         f: log10 fold enrichment
                         F: linear fold enrichment
                         d: subtraction
                         M: maximum
                         m: fragment pileup per million reads
        """
        if scoring_method == ord('p'):
            self.compute_pvalue()
        elif scoring_method == ord('q'):
            #if not already calculated p, compute pvalue first
            if self.scoring_method != ord('p'):
                self.compute_pvalue()
            self.compute_qvalue()
        elif scoring_method == ord('l'):
            self.compute_likelihood()
        elif scoring_method == ord('s'):
            self.compute_sym_likelihood()
        elif scoring_method == ord('f'):
            self.compute_logFE()
        elif scoring_method == ord('F'):
            self.compute_foldenrichment()
        elif scoring_method == ord('d'):
            self.compute_subtraction()
        elif scoring_method == ord('m'):
            self.compute_SPMR()
        elif scoring_method == ord('M'):
            self.compute_max()
        else:
            raise NotImplemented

    cdef compute_pvalue ( self ):
        """Compute -log_{10}(pvalue)
        """
        cdef:
            np.ndarray[np.float32_t] p, c, v
            np.ndarray[np.int32_t] pos
            int64_t l, i, prev_pos
            bytes chrom

        for chrom in sorted(self.data.keys()):
            prev_pos = 0
            pos = self.data[chrom][0]
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  get_pscore( <int32_t>(p[ i ]  + self.pseudocount) , c[ i ]  + self.pseudocount )
                try:
                    self.pvalue_stat[v[ i ]] += pos[ i ] - prev_pos
                except:
                    self.pvalue_stat[v[ i ]] = pos[ i ] - prev_pos
                prev_pos = pos[ i ]

        self.scoring_method = ord('p')
        return

    cdef compute_qvalue ( self ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            object pqtable
            int64_t i,l,j
            bytes chrom
            np.ndarray p, c, v

        # pvalue should be computed first!
        assert self.scoring_method == ord('p')
        # make pqtable
        pqtable = self.make_pq_table()

        # convert p to q
        for chrom in sorted(self.data.keys()):
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] = pqtable[ v[ i ] ]
                #v [ i ] =  g( v[ i ])

        self.scoring_method = ord('q')
        return

    cpdef object make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-log10pvalue:(-log10qvalue,rank,basepairs)} relationships.
        """
        cdef:
            int64_t n, pre_p, this_p, length, pre_l, l, i, j
            float32_t this_v, pre_v, v, q, pre_q # store the p and q scores
            int64_t N, k
            float32_t f
            bytes chrom
            np.ndarray v_chrom, pos_chrom
            object pvalue2qvalue
            dict pvalue_stat
            list unique_values

        assert self.scoring_method == ord('p')

        pvalue_stat = self.pvalue_stat

        N = sum(pvalue_stat.values())
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value

        pvalue2qvalue = Float32to32Map( for_int = False )
        unique_values = sorted(list(pvalue_stat.keys()), reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = pvalue_stat[v]
            q = v + (log10(k) + f)
            if q > pre_q:
                q = pre_q
            if q <= 0:
                q = 0
                break
            pvalue2qvalue[ v ] = q
            pre_q = q
            k+=l
        # bottom rank pscores all have qscores 0
        for j in range(i, len(unique_values) ):
            v = unique_values[ j ]
            pvalue2qvalue[ v ] = 0
        return pvalue2qvalue

    cdef compute_likelihood ( self ):
        """Calculate log10 likelihood.

        """
        cdef:
            #np.ndarray v, p, c
            int64_t l, i
            bytes chrom
            float32_t v1, v2
            float32_t pseudocount

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][ 1 ].flat.__next__ # pileup in treatment
            c = self.data[chrom][ 2 ].flat.__next__ # pileup in control
            v = self.data[chrom][ 3 ]               # score
            l = self.datalength[chrom]
            v1 = 2
            v2 = 1
            for i in range(l):
                v1 = p()
                v2 = c()
                v[ i ] =  logLR_asym( v1 + pseudocount, v2 + pseudocount )  #logLR( d[ i, 1]/100.0, d[ i, 2]/100.0 )
                #print v1, v2, v[i]
        self.scoring_method = ord('l')
        return

    cdef compute_sym_likelihood ( self ):
        """Calculate symmetric log10 likelihood.

        """
        cdef:
            #np.ndarray v, p, c
            int64_t l, i
            bytes chrom
            float32_t v1, v2
            float32_t pseudocount

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][ 1 ].flat.__next__
            c = self.data[chrom][ 2 ].flat.__next__
            v = self.data[chrom][ 3 ]
            l = self.datalength[chrom]
            v1 = 2
            v2 = 1
            for i in range(l):
                v1 = p()
                v2 = c()
                v[ i ] =  logLR_sym( v1 + pseudocount, v2 + pseudocount )  #logLR( d[ i, 1]/100.0, d[ i, 2]/100.0 )
        self.scoring_method = ord('s')
        return

    cdef compute_logFE ( self ):
        """Calculate log10 fold enrichment ( with 1 pseudocount ).

        """
        cdef:
            np.ndarray p, c, v
            int64_t l, i
            float32_t pseudocount

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] = get_logFE ( p[ i ] + pseudocount, c[ i ] + pseudocount)
        self.scoring_method = ord('f')
        return

    cdef compute_foldenrichment ( self ):
        """Calculate linear scale fold enrichment ( with 1 pseudocount ).

        """
        cdef:
            np.ndarray p, c, v
            int64_t l, i
            float32_t pseudocount

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  ( p[ i ] + pseudocount )/( c[ i ] + pseudocount )
        self.scoring_method = ord('F')
        return

    cdef compute_subtraction ( self ):
        cdef:
            np.ndarray p, c, v
            int64_t l, i

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] = p[ i ] - c[ i ]
        self.scoring_method = ord('d')
        return

    cdef compute_SPMR ( self ):
        cdef:
            np.ndarray p, v
            int64_t l, i
            float32_t scale
        if self.normalization_method == ord('T') or self.normalization_method == ord('N'):
            scale = self.treat_edm
        elif self.normalization_method == ord('C'):
            scale = self.ctrl_edm
        elif self.normalization_method == ord('M'):
            scale = 1

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  p[ i ] / scale # two digit precision may not be enough...
        self.scoring_method = ord('m')
        return

    cdef compute_max ( self ):
        cdef:
            np.ndarray p, c, v
            int64_t l, i

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] = max(p[ i ],c[ i ])
        self.scoring_method = ord('M')
        return

    cpdef write_bedGraph ( self, fhd, str name, str description, short column = 3):
        """Write all data to fhd in bedGraph Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 1: chip, 2: control, 3: score

        """
        cdef:
            bytes chrom
            int32_t l, pre, i, p
            float32_t pre_v, v
            set chrs
            np.ndarray pos, value

        assert column in range( 1, 4 ), "column should be between 1, 2 or 3."

        write = fhd.write

        if self.trackline:
            # this line is REQUIRED by the wiggle format for UCSC browser
            write( "track type=bedGraph name=\"%s\" description=\"%s\"\n" % ( name.decode(), description ) )

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            pos = self.data[ chrom ][ 0 ]
            value = self.data[ chrom ][ column ]
            l = self.datalength[ chrom ]
            pre = 0
            if pos.shape[ 0 ] == 0: continue # skip if there's no data
            pre_v = value[ 0 ]
            for i in range( 1, l ):
                v = value[ i ]
                p = pos[ i-1 ]
                #if ('%.5f' % pre_v) != ('%.5f' % v):
                if abs(pre_v - v) > 1e-5: # precision is 5 digits
                    write( "%s\t%d\t%d\t%.5f\n" % ( chrom.decode(), pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            write( "%s\t%d\t%d\t%.5f\n" % ( chrom.decode(), pre, p, pre_v ) )

        return True

    cpdef call_peaks (self, float32_t cutoff=5.0, int32_t min_length=200, int32_t max_gap=50, bool call_summits=False):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 5. For -log10pvalue, it means 10^-5.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        ptrack:  an optional track for pileup heights. If it's not None, use it to find summits. Otherwise, use self/scoreTrack.
        """
        cdef:
            int32_t i
            bytes chrom
            np.ndarray pos, sample, control, value, above_cutoff, above_cutoff_v, above_cutoff_endpos, above_cutoff_startpos, above_cutoff_sv
            list peak_content

        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        self.cutoff = cutoff
        for chrom in sorted(chrs):
            peak_content = []           # to store points above cutoff

            pos = self.data[chrom][ 0 ]
            sample = self.data[chrom][ 1 ]
            control = self.data[chrom][ 2 ]
            value = self.data[chrom][ 3 ]

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
                    if call_summits:
                        self.__close_peak2(peak_content, peaks, min_length, chrom, max_gap//2 )
                    else:
                        self.__close_peak(peak_content, peaks, min_length, chrom, max_gap//2 )
                    peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]

            # save the last peak
            if not peak_content:
                continue
            else:
                if call_summits:
                    self.__close_peak2(peak_content, peaks, min_length, chrom, max_gap//2 )
                else:
                    self.__close_peak(peak_content, peaks, min_length, chrom, max_gap//2 )

        return peaks

    cdef bool __close_peak (self, list peak_content, peaks, int32_t min_length,
                            bytes chrom, int32_t smoothlen=0):
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peaks: a PeakIO object

        """
        cdef:
            int32_t summit_pos, tstart, tend, tmpindex, summit_index, i, midindex
            float32_t summit_value, tvalue, tsummitvalue

        peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
        if peak_length >= min_length: # if the peak is too small, reject it
            tsummit = []
            summit_pos   = 0
            summit_value = 0
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
            if self.scoring_method == ord('q'):
                qscore = self.data[chrom][3][ summit_index ]
            else:
                # if q value is not computed, use -1
                qscore = -1

            peaks.add( chrom,
                       peak_content[0][0],
                       peak_content[-1][1],
                       summit      = summit_pos,
                       peak_score  = self.data[chrom][ 3 ][ summit_index ],
                       pileup      = self.data[chrom][ 1 ][ summit_index ], # should be the same as summit_value
                       pscore      = get_pscore(self.data[chrom][ 1 ][ summit_index ], self.data[chrom][ 2 ][ summit_index ]),
                       fold_change = ( self.data[chrom][ 1 ][ summit_index ] + self.pseudocount ) / ( self.data[chrom][ 2 ][ summit_index ] + self.pseudocount ),
                       qscore      = qscore,
                       )
            # start a new peak
            return True

    cdef bool __close_peak2 (self, list peak_content, peaks, int32_t min_length,
                             bytes chrom, int32_t smoothlen=51,
                             float32_t min_valley = 0.9):
        cdef:
            int32_t summit_pos, tstart, tend, tmpindex, summit_index, summit_offset
            int32_t start, end, i, j, start_boundary
            float32_t summit_value, tvalue, tsummitvalue
#            np.ndarray[np.float32_t, ndim=1] w
            np.ndarray[np.float32_t, ndim=1] peakdata
            np.ndarray[np.int32_t, ndim=1] peakindices, summit_offsets

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

        peakdata = np.zeros(end - start, dtype='float32')
        peakindices = np.zeros(end - start, dtype='int32')
        for (tstart,tend,tvalue,tsvalue, tmpindex) in peak_content:
            i = tstart - start + start_boundary
            j = tend - start + start_boundary
            peakdata[i:j] = tsvalue
            peakindices[i:j] = tmpindex
        summit_offsets = maxima(peakdata, smoothlen)
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak(peak_content, peaks, min_length, chrom)
        else:
            # remove maxima that occurred in padding
            i = np.searchsorted(summit_offsets, start_boundary)
            j = np.searchsorted(summit_offsets, peak_length + start_boundary, 'right')
            summit_offsets = summit_offsets[i:j]

        summit_offsets = enforce_peakyness(peakdata, summit_offsets)
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak(peak_content, peaks, min_length, chrom)

        summit_indices = peakindices[summit_offsets]
        summit_offsets -= start_boundary

        peak_scores  = self.data[chrom][3][ summit_indices ]
        if not (peak_scores > self.cutoff).all():
            return self.__close_peak(peak_content, peaks, min_length, chrom)
        for summit_offset, summit_index in zip(summit_offsets, summit_indices):
            if self.scoring_method == ord('q'):
                qscore = self.data[chrom][3][ summit_index ]
            else:
                # if q value is not computed, use -1
                qscore = -1
            peaks.add( chrom,
                       start,
                       end,
                       summit      = start + summit_offset,
                       peak_score  = self.data[chrom][3][ summit_index ],
                       pileup      = self.data[chrom][1][ summit_index ], # should be the same as summit_value
                       pscore      = get_pscore(self.data[chrom][ 1 ][ summit_index ], self.data[chrom][ 2 ][ summit_index ]),
                       fold_change = ( self.data[chrom][ 1 ][ summit_index ] + self.pseudocount ) / ( self.data[chrom][ 2 ][ summit_index ] + self.pseudocount ),
                       qscore      = qscore,
                       )
        # start a new peak
        return True

    cdef int64_t total ( self ):
        """Return the number of regions in this object.

        """
        cdef:
            int64_t t
            bytes chrom

        t = 0
        for chrom in sorted(self.data.keys()):
            t += self.datalength[chrom]
        return t

    cpdef call_broadpeaks (self, float32_t lvl1_cutoff=5.0, float32_t lvl2_cutoff=1.0, int32_t min_length=200, int32_t lvl1_max_gap=50, int32_t lvl2_max_gap=400):
        """This function try to find enriched regions within which,
        scores are continuously higher than a given cutoff for level
        1, and link them using the gap above level 2 cutoff with a
        maximum length of lvl2_max_gap.

        lvl1_cutoff:  cutoff of value at enriched regions, default 5.0.
        lvl2_cutoff:  cutoff of value at linkage regions, default 1.0.
        min_length :  minimum peak length, default 200.
        lvl1_max_gap   :  maximum gap to merge nearby enriched peaks, default 50.
        lvl2_max_gap   :  maximum length of linkage regions, default 400.

        Return both general PeakIO object for highly enriched regions
        and gapped broad regions in BroadPeakIO.
        """
        cdef:
            int32_t i
            bytes chrom

        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."
        lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff, min_length=min_length, max_gap=lvl1_max_gap)
        lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff, min_length=min_length, max_gap=lvl2_max_gap)
        chrs = lvl1_peaks.peaks.keys()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in sorted(chrs):
            lvl1peakschrom = lvl1_peaks.peaks[chrom]
            lvl2peakschrom = lvl2_peaks.peaks[chrom]
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

    def __add_broadpeak (self, bpeaks, bytes chrom, dict lvl2peak, list lvl1peakset):
        """Internal function to create broad peak.
        """

        cdef:
            int32_t blockNum, thickStart, thickEnd, start, end
            bytes blockSizes, blockStarts

        start      = lvl2peak["start"]
        end        = lvl2peak["end"]

        # the following code will add those broad/lvl2 peaks with no strong/lvl1 peaks inside
        if not lvl1peakset:
            # will complement by adding 1bps start and end to this region
            # may change in the future if gappedPeak format was improved.
            bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=(b"%d" % start), thickEnd=(b"%d" % end),
                       blockNum = 2, blockSizes = b"1,1", blockStarts = (b"0,%d" % (end-start-1)), pileup = lvl2peak["pileup"],
                       pscore = lvl2peak["pscore"], fold_change = lvl2peak["fc"],
                       qscore = lvl2peak["qscore"] )
            return bpeaks

        thickStart = b"%d" % lvl1peakset[0]["start"]
        thickEnd   = b"%d" % lvl1peakset[-1]["end"]
        blockNum   = int(len(lvl1peakset))
        blockSizes = b",".join( [b"%d" % x["length"] for x in lvl1peakset] )
        blockStarts = b",".join( [b"%d" % (x["start"]-start) for x in lvl1peakset] )

        if lvl2peak["start"] != thickStart:
            # add 1bp mark for the start of lvl2 peak
            thickStart = b"%d" % start
            blockNum += 1
            blockSizes = b"1,"+blockSizes
            blockStarts = b"0,"+blockStarts
        if lvl2peak["end"] != thickEnd:
            # add 1bp mark for the end of lvl2 peak
            thickEnd = b"%d" % end
            blockNum += 1
            blockSizes = blockSizes+b",1"
            blockStarts = blockStarts + b"," + (b"%d" % (end-start-1))

        # add to BroadPeakIO object
        bpeaks.add(chrom, start, end, score=lvl2peak["score"], thickStart=thickStart, thickEnd=thickEnd,
                   blockNum = blockNum, blockSizes = blockSizes, blockStarts = blockStarts,  pileup = lvl2peak["pileup"],
                   pscore = lvl2peak["pscore"], fold_change = lvl2peak["fc"],
                   qscore = lvl2peak["qscore"] )
        return bpeaks

cdef class TwoConditionScores:
    """Class for saving two condition comparison scores.
    """
    cdef:
        dict data                       # dictionary for data of each chromosome
        dict datalength                 # length of data array of each chromosome
        float32_t cond1_factor              # factor to apply to cond1 pileup values
        float32_t cond2_factor              # factor to apply to cond2 pileup values
        float32_t pseudocount               # the pseudocount used to calcuate LLR
        float32_t cutoff
        object t1bdg, c1bdg, t2bdg, c2bdg
        dict pvalue_stat1, pvalue_stat2, pvalue_stat3

    def __init__ (self, t1bdg, c1bdg, t2bdg, c2bdg, float32_t cond1_factor = 1.0, float32_t cond2_factor = 1.0, float32_t pseudocount = 0.01, float32_t proportion_background_empirical_distribution = 0.99999 ):
        """
        t1bdg: a bedGraphTrackI object for treat 1
        c1bdg: a bedGraphTrackI object for control 1
        t2bdg: a bedGraphTrackI object for treat 2
        c2bdg: a bedGraphTrackI object for control 2

        cond1_factor: this will be multiplied to values in t1bdg and c1bdg
        cond2_factor: this will be multiplied to values in t2bdg and c2bdg

        pseudocount: pseudocount, by default 0.01.

        proportion_background_empirical_distribution: proportion of genome as the background to build empirical distribution

        """

        self.data = {}           # for each chromosome, there is a l*4
                                 # matrix. First column: end position
                                 # of a region; Second: treatment
                                 # pileup; third: control pileup ;
                                 # forth: score ( can be
                                 # p/q-value/likelihood
                                 # ratio/fold-enrichment/subtraction
                                 # depending on -c setting)
        self.datalength = {}
        self.cond1_factor = cond1_factor
        self.cond2_factor = cond2_factor
        self.pseudocount = pseudocount
        self.pvalue_stat1 = {}
        self.pvalue_stat2 = {}
        self.t1bdg = t1bdg
        self.c1bdg = c1bdg
        self.t2bdg = t2bdg
        self.c2bdg = c2bdg

        #self.empirical_distr_llr = [] # save all values in histogram

    cpdef set_pseudocount( self, float32_t pseudocount ):
        self.pseudocount = pseudocount

    cpdef build ( self ):
        """Compute scores from 3 types of comparisons and store them in self.data.

        """
        cdef:
            set common_chrs
            bytes chrname
            int32_t chrom_max_len
        # common chromosome names
        common_chrs = self.get_common_chrs()
        for chrname in common_chrs:
            (cond1_treat_ps, cond1_treat_vs) = self.t1bdg.get_data_by_chr(chrname)
            (cond1_control_ps, cond1_control_vs) = self.c1bdg.get_data_by_chr(chrname)
            (cond2_treat_ps, cond2_treat_vs) = self.t2bdg.get_data_by_chr(chrname)
            (cond2_control_ps, cond2_control_vs) = self.c2bdg.get_data_by_chr(chrname)
            chrom_max_len = len(cond1_treat_ps) + len(cond1_control_ps) +\
                            len(cond2_treat_ps) + len(cond2_control_ps)
            self.add_chromosome( chrname, chrom_max_len )
            self.build_chromosome( chrname,
                                   cond1_treat_ps, cond1_control_ps,
                                   cond2_treat_ps, cond2_control_ps,
                                   cond1_treat_vs, cond1_control_vs,
                                   cond2_treat_vs, cond2_control_vs )


    cdef build_chromosome( self, chrname,
                           cond1_treat_ps, cond1_control_ps,
                           cond2_treat_ps, cond2_control_ps,
                           cond1_treat_vs, cond1_control_vs,
                           cond2_treat_vs, cond2_control_vs ):
        """Internal function to calculate scores for three types of comparisons.

        cond1_treat_ps, cond1_control_ps: position of treat and control of condition 1
        cond2_treat_ps, cond2_control_ps: position of treat and control of condition 2
        cond1_treat_vs, cond1_control_vs: value of treat and control of condition 1
        cond2_treat_vs, cond2_control_vs: value of treat and control of condition 2

        """
        cdef:
            int32_t c1tp, c1cp, c2tp, c2cp, minp, pre_p
            float32_t c1tv, c1cv, c2tv, c2cv 
        c1tpn = iter(cond1_treat_ps).__next__
        c1cpn = iter(cond1_control_ps).__next__
        c2tpn = iter(cond2_treat_ps).__next__
        c2cpn = iter(cond2_control_ps).__next__
        c1tvn = iter(cond1_treat_vs).__next__
        c1cvn = iter(cond1_control_vs).__next__
        c2tvn = iter(cond2_treat_vs).__next__
        c2cvn = iter(cond2_control_vs).__next__

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
        return

    cdef set get_common_chrs ( self ):
        cdef:
            set t1chrs, c1chrs, t2chrs, c2chrs, common
        t1chrs = self.t1bdg.get_chr_names()
        c1chrs = self.c1bdg.get_chr_names()
        t2chrs = self.t2bdg.get_chr_names()
        c2chrs = self.c2bdg.get_chr_names()
        common = reduce(lambda x,y:x.intersection(y), (t1chrs,c1chrs,t2chrs,c2chrs))
        return common

    cdef add_chromosome ( self, bytes chrom, int32_t chrom_max_len ):
        """
        chrom: chromosome name
        chrom_max_len: maximum number of data points in this chromosome

        """
        if chrom not in self.data:
            self.data[chrom] = [ np.zeros( chrom_max_len, dtype="int32" ), # pos
                                 np.zeros( chrom_max_len, dtype="float32" ), # LLR t1 vs c1
                                 np.zeros( chrom_max_len, dtype="float32" ), # LLR t2 vs c2
                                 np.zeros( chrom_max_len, dtype="float32" )] # LLR t1 vs t2
            self.datalength[chrom] = 0

    cdef add (self, bytes chromosome, int32_t endpos, float32_t t1, float32_t c1, float32_t t2, float32_t c2):
        """Take chr-endpos-sample1-control1-sample2-control2 and
        compute logLR for t1 vs c1, t2 vs c2, and t1 vs t2, then save
        values.

        chromosome: chromosome name in string
        endpos    : end position of each interval in integer
        t1        : Sample 1 ChIP pileup value of each interval in float
        c1        : Sample 1 Control pileup value of each interval in float
        t2        : Sample 2 ChIP pileup value of each interval in float
        c2        : Sample 2 Control pileup value of each interval in float

        *Warning* Need to add regions continuously.
        """
        cdef:
            int32_t i
            list c
        i = self.datalength[chromosome]
        c = self.data[chromosome]
        c[0][ i ] = endpos
        c[1][ i ] = logLR_asym( (t1+self.pseudocount)*self.cond1_factor, (c1+self.pseudocount)*self.cond1_factor )
        c[2][ i ] = logLR_asym( (t2+self.pseudocount)*self.cond2_factor, (c2+self.pseudocount)*self.cond2_factor )
        c[3][ i ] = logLR_sym( (t1+self.pseudocount)*self.cond1_factor, (t2+self.pseudocount)*self.cond2_factor )
        self.datalength[chromosome] += 1
        return

    cpdef finalize ( self ):
        """
        Adjust array size of each chromosome.

        """
        cdef:
            bytes chrom
            int32_t l
            list d

        for chrom in sorted(self.data.keys()):
            d = self.data[chrom]
            l = self.datalength[chrom]
            d[0].resize( l, refcheck = False )
            d[1].resize( l, refcheck = False )
            d[2].resize( l, refcheck = False )
            d[3].resize( l, refcheck = False )
        return

    cpdef get_data_by_chr (self, bytes chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if chromosome in self.data:
            return self.data[chromosome]
        else:
            return None

    cpdef get_chr_names (self):
        """Return all the chromosome names stored.

        """
        l = set(self.data.keys())
        return l

    cpdef write_bedGraph ( self, fhd, str name, str description, int32_t column = 3):
        """Write all data to fhd in bedGraph Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 1: cond1 chip vs cond1 ctrl, 2: cond2 chip vs cond2 ctrl, 3: cond1 chip vs cond2 chip

        """
        cdef:
            bytes chrom
            int32_t l, pre, i, p
            float32_t pre_v, v
            np.ndarray pos, value

        assert column in range( 1, 4 ), "column should be between 1, 2 or 3."

        write = fhd.write

        #if self.trackline:
        #    # this line is REQUIRED by the wiggle format for UCSC browser
        #    write( "track type=bedGraph name=\"%s\" description=\"%s\"\n" % ( name.decode(), description ) )

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
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
                    write( "%s\t%d\t%d\t%.5f\n" % ( chrom.decode(), pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            write( "%s\t%d\t%d\t%.5f\n" % ( chrom.decode(), pre, p, pre_v ) )

        return True

    cpdef write_matrix ( self, fhd, str name, str description ):
        """Write all data to fhd into five columns Format:

        col1: chr_start_end
        col2: t1 vs c1
        col3: t2 vs c2
        col4: t1 vs t2

        fhd: a filehandler to save the matrix.

        """
        cdef:
            bytes chrom
            int32_t l, pre, i, p
            float32_t v1, v2, v3
            np.ndarray pos, value1, value2, value3

        write = fhd.write

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            [ pos, value1, value2, value3 ] = self.data[ chrom ]
            l = self.datalength[ chrom ]
            pre = 0
            if pos.shape[ 0 ] == 0: continue # skip if there's no data
            for i in range( 0, l ):
                v1 = value1[ i ]
                v2 = value2[ i ]
                v3 = value3[ i ]
                p = pos[ i ]
                write( "%s:%d_%d\t%.5f\t%.5f\t%.5f\n" % ( chrom.decode(), pre, p, v1, v2, v3 ) )
                pre = p

        return True

    cpdef tuple call_peaks (self, float32_t cutoff=3, int32_t min_length=200, int32_t max_gap = 100,
                      bool call_summits=False):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        For bdgdiff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 3. For log10 LR, it means 1000 or -1000.
        min_length :  minimum peak length, default 200.
        max_gap   :  maximum gap to merge nearby peaks, default 100.
        ptrack:  an optional track for pileup heights. If it's not None, use it to find summits. Otherwise, use self/scoreTrack.
        """
        cdef:
            int32_t i
            bytes chrom
            np.ndarray pos, t1_vs_c1, t2_vs_c2, t1_vs_t2, \
                       cond1_over_cond2, cond2_over_cond1, cond1_equal_cond2, \
                       cond1_sig, cond2_sig,\
                       cat1, cat2, cat3, \
                       cat1_startpos, cat1_endpos, cat2_startpos, cat2_endpos, \
                       cat3_startpos, cat3_endpos
        chrs  = self.get_chr_names()
        cat1_peaks = PeakIO()       # dictionary to save peaks significant at condition 1
        cat2_peaks = PeakIO()       # dictionary to save peaks significant at condition 2
        cat3_peaks = PeakIO()       # dictionary to save peaks significant in both conditions

        self.cutoff = cutoff

        for chrom in sorted(chrs):
            pos = self.data[chrom][ 0 ]
            t1_vs_c1 = self.data[chrom][ 1 ]
            t2_vs_c2 = self.data[chrom][ 2 ]
            t1_vs_t2 = self.data[chrom][ 3 ]
            and_ = np.logical_and
            cond1_over_cond2 = t1_vs_t2 >= cutoff # regions with stronger cond1 signals
            cond2_over_cond1 = t1_vs_t2 <= -1*cutoff # regions with stronger cond2 signals
            cond1_equal_cond2= and_( t1_vs_t2 >= -1*cutoff, t1_vs_t2 <= cutoff )
            cond1_sig = t1_vs_c1 >= cutoff # enriched regions in condition 1
            cond2_sig = t2_vs_c2 >= cutoff # enriched regions in condition 2
            # indices where score is above cutoff
            cat1 = np.where( and_( cond1_sig, cond1_over_cond2 ) )[ 0 ] # cond1 stronger than cond2, the indices
            cat2 = np.where( and_( cond2_over_cond1, cond2_sig ) )[ 0 ] # cond2 stronger than cond1, the indices
            cat3 = np.where( and_( and_( cond1_sig, cond2_sig ), # cond1 and cond2 are equal, the indices
                                   cond1_equal_cond2 ) ) [ 0 ]

            cat1_endpos = pos[cat1] # end positions of regions where score is above cutoff
            cat1_startpos = pos[cat1-1] # start positions of regions where score is above cutoff
            cat2_endpos = pos[cat2] # end positions of regions where score is above cutoff
            cat2_startpos = pos[cat2-1] # start positions of regions where score is above cutoff
            cat3_endpos = pos[cat3] # end positions of regions where score is above cutoff
            cat3_startpos = pos[cat3-1] # start positions of regions where score is above cutoff

            # for cat1: condition 1 stronger regions
            self.__add_a_peak ( cat1_peaks, chrom, cat1, cat1_startpos, cat1_endpos, t1_vs_t2, max_gap, min_length )
            # for cat2: condition 2 stronger regions
            self.__add_a_peak ( cat2_peaks, chrom, cat2, cat2_startpos, cat2_endpos, -1 * t1_vs_t2, max_gap, min_length )
            # for cat3: commonly strong regions
            self.__add_a_peak ( cat3_peaks, chrom, cat3, cat3_startpos, cat3_endpos, abs(t1_vs_t2), max_gap, min_length )

        return (cat1_peaks, cat2_peaks, cat3_peaks)

    cdef object __add_a_peak ( self, object peaks, bytes chrom, np.ndarray indices, np.ndarray startpos, np.ndarray endpos,
                               np.ndarray score, int32_t max_gap, int32_t min_length ):
         """For a given chromosome, merge nearby significant regions,
         filter out smaller regions, then add regions to PeakIO
         object.

         """
         cdef:
             int32_t i
             list peak_content
             float32_t mean_logLR

         if startpos.size > 0:
             # if it is not empty
             peak_content = []
             if indices[0] == 0:
                 # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
                 startpos[0] = 0
             # first bit of region above cutoff
             peak_content.append( (startpos[0], endpos[0], score[indices[ 0 ]]) )
             for i in range( 1, startpos.size ):
                 if startpos[i] - peak_content[-1][1] <= max_gap:
                     # append
                     peak_content.append( ( startpos[i], endpos[i], score[indices[ i ]] ) )
                 else:
                     # close
                     if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:
                         mean_logLR = self.mean_from_peakcontent( peak_content )
                         #if peak_content[0][0] == 22414956:
                         #    print(f"{peak_content} {mean_logLR}")
                         peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                    summit = -1, peak_score  = mean_logLR, pileup = 0, pscore = 0,
                                    fold_change = 0, qscore = 0,
                                    )
                     peak_content = [(startpos[i], endpos[i], score[ indices[ i ] ]),]

             # save the last peak
             if peak_content:
                 if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:
                     mean_logLR = self.mean_from_peakcontent( peak_content )
                     peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                summit = -1, peak_score  = mean_logLR, pileup = 0, pscore = 0,
                                fold_change = 0, qscore = 0,
                                )

         return

    cdef float32_t mean_from_peakcontent ( self, list peakcontent ):
        """

        """
        cdef:
            int32_t tmp_s, tmp_e
            int32_t l
            float64_t tmp_v, sum_v        #for better precision
            float32_t r
            int32_t i

        l = 0
        sum_v = 0                         #initialize sum_v as 0
        for i in range( len(peakcontent) ):
            tmp_s = peakcontent[i][0]
            tmp_e = peakcontent[i][1]
            tmp_v = peakcontent[i][2]
            sum_v += tmp_v * ( tmp_e - tmp_s )
            l +=  tmp_e - tmp_s

        r = <float32_t>( sum_v / l )
        return r


    cdef int64_t total ( self ):
        """Return the number of regions in this object.

        """
        cdef:
            int64_t t
            bytes chrom

        t = 0
        for chrom in sorted(self.data.keys()):
            t += self.datalength[chrom]
        return t


