# cython: profile=True
# Time-stamp: <2012-06-22 14:39:05 Tao Liu>

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
from MACS2.cSignal import maxima, enforce_valleys, enforce_peakyness
#np_convolve = np.convolve

from libc.math cimport log10,log

from MACS2.Constants import *
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO, BroadPeakIO

from MACS2.hashtable import Int64HashTable, Float64HashTable

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
    #if pscore_dict.has_key(key_value):
    #    return pscore_dict[key_value]
    #else:
    #    score = int(-100*poisson_cdf(observed,expectation,False,True))
    #    pscore_dict[(observed,expectation)] = score
    #return score

logLR_khashtable = Int64HashTable()

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
        if x > y:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
        elif x < y:
            s = (x*(-log(x)+log(y))-y+x)*LOG10_E
        else:
            s = 0
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
            np.ndarray[np.int32_t, ndim=1] pos
            np.ndarray[np.float32_t, ndim=1] value
        
        if colname not in ['V1','V2']:
            raise Exception("%s not supported!" % colname)
        chrs = self.get_chr_names()
        write = fhd.write
        for chrom in chrs:
            d = self.data[chrom]
            l = self.pointer[chrom]
            pos   = d['pos']
            value = d[colname]
            pre = 0
            for i in range( l ):
                write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,pos[i],value[i]))
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
        float pseudocount                # the pseudocount used to calcuate logLR, FE or logFE
        float cutoff
        dict pvalue_stat                 # save pvalue<->length dictionary

    
    def __init__ (self, float treat_depth, float ctrl_depth, bool stderr_on = False, float pseudocount = 1.0 ):
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

    cpdef set_pseudocount( self, float pseudocount ):
        self.pseudocount = pseudocount
        
    cpdef enable_trackline( self ):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True

    cpdef add_chromosome ( self, str chrom, int chrom_max_len ):
        """
        chrom: chromosome name
        chrom_max_len: maximum number of data points in this chromosome
        
        """
        if not self.data.has_key(chrom):
            #self.data[chrom] = np.zeros( ( chrom_max_len, 4 ), dtype="int32" ) # remember col #2-4 is actual value * 100, I use integer here.
            self.data[chrom] = [ np.zeros( chrom_max_len, dtype="int32" ), # pos
                                 np.zeros( chrom_max_len, dtype="float32" ), # pileup at each interval, in float format
                                 np.zeros( chrom_max_len, dtype="float32" ), # control at each interval, in float format
                                 np.zeros( chrom_max_len, dtype="float32" ) ] # score at each interval, in float format
            self.datalength[chrom] = 0

    cpdef add (self, str chromosome, int endpos, float chip, float control):
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
            str chrom, k
            int l

        for chrom in self.data.keys():
            d = self.data[chrom]
            l = self.datalength[chrom]
            d[0].resize( l, refcheck = False )
            d[1].resize( l, refcheck = False )
            d[2].resize( l, refcheck = False )
            d[3].resize( l, refcheck = False )            
        return

    # cpdef sort ( self, int column = 1 ):
    #     """ Sort data for each chromosome, by certain column.

    #     column: 1: position, 2: sample, 3: control, 4: score

    #     Default: sort by positions.
    #     """
    #     cdef:
    #         str chrom
            
        
    #     for chrom in self.data.keys():
    #         d = self.data[chrom]
    #         d.view('int32,int32,int32,int32').sort(axis=0,order=column-1)
    #     return

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
            np.ndarray p, c
            long l, i
        
        for chrom in self.data.keys():
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
                         f: log10 fold enrichment
                         F: linear fold enrichment
                         d: subtraction
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
            self.compute_logFE()
        elif scoring_method == 'F':
            self.compute_foldenrichment()
        elif scoring_method == 'd':
            self.compute_subtraction()
        elif scoring_method == 'm':
            self.compute_SPMR()
        else:
            raise NotImplemented
            
    cdef compute_pvalue ( self ):
        """Compute -log_{10}(pvalue)
        """
        cdef:
            np.ndarray[np.float32_t] p, c, v
            np.ndarray[np.int32_t] pos
            long l, i, prev_pos
            str chrom
        
        for chrom in self.data.keys():
            prev_pos = 0
            pos = self.data[chrom][0]
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  get_pscore( int(p[ i ]) , c[ i ] )
                try:
                    self.pvalue_stat[v[ i ]] += pos[ i ] - prev_pos
                except:
                    self.pvalue_stat[v[ i ]] = pos[ i ] - prev_pos
                prev_pos = pos[ i ]
                    
        self.scoring_method = 'p'
        return 

    cdef compute_qvalue ( self ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            dict pqtable
            long i,l,j
            double k
            str chrom
            np.ndarray p, c, v
            
        # pvalue should be computed first!
        assert self.scoring_method == 'p'
        # make pqtable
        pqtable = self.make_pq_table()
        
        # convert p to q

        # convert pvalue2qvalue to a simple dict based on khash
        # khash has big advantage while checking keys for millions of times.
        s_p2q = Float64HashTable()
        for k in pqtable.keys():
            s_p2q.set_item(k,pqtable[k])

        g = s_p2q.get_item
        
        for chrom in self.data.keys():
            #p = self.data[chrom][1]
            #c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  g( v[ i ])
        
        self.scoring_method = 'q'
        return

    cdef dict make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-log10pvalue:(-log10qvalue,rank,basepairs)} relationships.
        """
        cdef:
            long n, pre_p, this_p, length, j, pre_l, l, i
            double this_v, pre_v, v, q, pre_q
            long N, k
            double f
            str chrom
            np.ndarray v_chrom, pos_chrom
            dict pvalue2qvalue
            dict value_dict
            list unique_values

        assert self.scoring_method == 'p'

        # value_dict = {} #Float64HashTable()
        # #unique_values = pyarray(FBYTE4,[])
        # # this is a table of how many positions each p value occurs at
        # chroms = self.data.keys()
        # for chrom in self.data.keys():
        #     # for each chromosome
        #     pre_p  = 0
        #     #d_chrom = self.data[chrom]
        #     pos_chrom = self.data[chrom][0]
        #     v_chrom = self.data[chrom][3]            
        #     length = self.datalength[chrom]
        #     for j in range(length):
        #         this_p = pos_chrom[ j ]
        #         this_v = v_chrom[ j ]
        #         #assert this_v == this_v, "NaN at %d" % pos
        #         if value_dict.has_key(this_v):
        #             #value_dict.set_item(this_v, value_dict.get_item(this_v) + this_p - pre_p)
        #             value_dict[this_v] += this_p - pre_p
        #         else:
        #             #value_dict.set_item(this_v, this_p - pre_p)
        #             value_dict[this_v] = this_p - pre_p
        #             #unique_values.append(this_v)
        #         pre_p = this_p
        value_dict = self.pvalue_stat
        #logging.info("####test#### 2")
        N = sum(value_dict.values())
        #for i in range(len(unique_values)):
            #N += value_dict.get_item(unique_values[i])
        k = 1                           # rank
        f = -log10(N)
        pre_v = -2147483647
        pre_l = 0
        pre_q = 2147483647              # save the previous q-value
        #pvalue2qvalue = {pre_v:[0,k,0]}              # pvalue:[qvalue,rank,bp_with_this_pvalue]
        pvalue2qvalue = {}#Float64HashTable()
        unique_values = sorted(value_dict.keys(), reverse=True) #sorted(unique_values,reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            #l = value_dict.get_item(v)
            l = value_dict[v]
            q = v + (log10(k) + f)
            q = max(0,min(pre_q,q))           # make q-score monotonic
            #pvalue2qvalue[v] = [q, k, 0]
            #pvalue2qvalue.set_item( v, q )
            pvalue2qvalue[ v ] = q
            #pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
            pre_v = v
            pre_q = q
            k+=l
        #pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
        # pop the first -1e100 one
        #pvalue2qvalue.pop(-2147483647)
        #logging.info("####test#### 3")
        return pvalue2qvalue

    cdef compute_likelihood ( self ):
        """Calculate log10 likelihood.
        
        """
        cdef:
            #np.ndarray v, p, c
            long l, i
            str chrom
            float v1, v2
            float pseudocount

        pseudocount = self.pseudocount
        
        for chrom in self.data.keys():
            p = self.data[chrom][ 1 ].flat.next
            c = self.data[chrom][ 2 ].flat.next
            v = self.data[chrom][ 3 ]
            l = self.datalength[chrom]
            v1 = 2
            v2 = 1
            for i in range(l):
                v1 = p() 
                v2 = c()
                v[ i ] =  logLR( v1 + pseudocount, v2 + pseudocount )  #logLR( d[ i, 1]/100.0, d[ i, 2]/100.0 )
        self.scoring_method = 'l'
        return 

    cdef compute_logFE ( self ):
        """Calculate log10 fold enrichment ( with 1 pseudocount ).

        """
        cdef:
            np.ndarray p, c, v
            long l, i
            float pseudocount

        pseudocount = self.pseudocount
        
        for chrom in self.data.keys():
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] = get_logFE ( p[ i ] + pseudocount, c[ i ] + pseudocount)
        self.scoring_method = 'f'
        return

    cdef compute_foldenrichment ( self ):
        """Calculate linear scale fold enrichment ( with 1 pseudocount ).

        """
        cdef:
            np.ndarray p, c, v
            long l, i
            float pseudocount

        pseudocount = self.pseudocount
        
        for chrom in self.data.keys():
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  ( p[ i ] + pseudocount )/( c[ i ] + pseudocount )
        self.scoring_method = 'F'
        return

    cdef compute_subtraction ( self ):
        cdef:
            np.ndarray p, c, v
            long l, i
        
        for chrom in self.data.keys():
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] = p[ i ] - c[ i ]
        self.scoring_method = 'd'
        return

    cdef compute_SPMR ( self ):
        cdef:
            np.ndarray p, v
            long l, i
            float scale
        if self.normalization_method == 'T' or self.normalization_method == 'N':
            scale = self.treat_edm
        elif self.normalization_method == 'C':
            scale = self.ctrl_edm
        elif self.normalization_method == 'M':
            scale = 1
        
        for chrom in self.data.keys():
            p = self.data[chrom][1]
            v = self.data[chrom][3]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  p[ i ] / scale # two digit precision may not be enough...
        self.scoring_method = 'm'
        return

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
                if pre_v != v: 
                    write( "%s\t%d\t%d\t%.2f\n" % ( chrom, pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            write( "%s\t%d\t%d\t%.2f\n" % ( chrom, pre, p, pre_v ) )
            
        return True

    cpdef call_peaks (self, float cutoff=5.0, int min_length=200, int max_gap=50, bool call_summits=False):
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
            int i
            str chrom
            np.ndarray pos, sample, control, value, above_cutoff, above_cutoff_v, above_cutoff_endpos, above_cutoff_startpos, above_cutoff_sv
            list peak_content
        
        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        self.cutoff = cutoff
        for chrom in chrs:
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
                        self.__close_peak2(peak_content, peaks, min_length, chrom, max_gap/2 )
                    else:
                        self.__close_peak(peak_content, peaks, min_length, chrom, max_gap/2 )
                    peak_content = [(above_cutoff_startpos[i], above_cutoff_endpos[i], above_cutoff_v[i], above_cutoff_sv[i], above_cutoff[i]),]
            
            # save the last peak
            if not peak_content:
                continue
            else:
                if call_summits:
                    self.__close_peak2(peak_content, peaks, min_length, chrom, max_gap/2 )
                else:
                    self.__close_peak(peak_content, peaks, min_length, chrom, max_gap/2 )

        return peaks

    cdef bool __close_peak (self, list peak_content, peaks, int min_length,
                            str chrom, int smoothlen=0):
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
                       fold_change = float ( self.data[chrom][ 1 ][ summit_index ] + self.pseudocount ) / ( self.data[chrom][ 2 ][ summit_index ] + self.pseudocount ),
                       qscore      = qscore,
                       )
            # start a new peak
            return True

    cdef bool __close_peak2 (self, list peak_content, peaks, int min_length,
                             str chrom, int smoothlen=51,
                             float min_valley = 0.9):
        cdef:
            int summit_pos, tstart, tend, tmpindex, summit_index, summit_offset
            int start, end, i, j, start_boundary
            double summit_value, tvalue, tsummitvalue
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
        # apply smoothing window of smoothlen
#        w = np.ones(smoothlen, dtype='float32') / smoothlen
#        if smoothlen > 0:
#            smoothdata = np_convolve(w, peakdata, mode='same')
#        else:
#            smoothdata = peakdata.copy()
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
        
#        summit_offsets = enforce_valleys(peakdata, summit_offsets, min_valley = min_valley)
        summit_indices = peakindices[summit_offsets]
        summit_offsets -= start_boundary

        peak_scores  = self.data[chrom][3][ summit_indices ]
        if not (peak_scores > self.cutoff).all():
            return self.__close_peak(peak_content, peaks, min_length, chrom)
        for summit_offset, summit_index in zip(summit_offsets, summit_indices):
            if self.scoring_method == 'q':
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
                       fold_change = float ( self.data[chrom][ 1 ][ summit_index ] + self.pseudocount ) / ( self.data[chrom][ 2 ][ summit_index ] + self.pseudocount ),
                       qscore      = qscore,
                       )
        # start a new peak
        return True
        
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

    cpdef tuple call_broadpeaks (self, float lvl1_cutoff=5.0, float lvl2_cutoff=1.0, int min_length=200, int lvl1_max_gap=50, int lvl2_max_gap=400):
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
            int i
            str chrom
        
        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."        
        lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff, min_length=min_length, max_gap=lvl1_max_gap)
        lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff, min_length=min_length, max_gap=lvl2_max_gap)
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

    def __add_broadpeak (self, bpeaks, str chrom, dict lvl2peak, list lvl1peakset):
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
