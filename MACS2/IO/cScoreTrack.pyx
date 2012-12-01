# Time-stamp: <2012-10-18 15:46:03 Tao Liu>

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

from operator import itemgetter

from cpython cimport bool
#from scipy.signal import fftconvolve
from MACS2.cSignal import maxima, enforce_valleys, enforce_peakyness
#np_convolve = np.convolve

# Experimental
from scipy.stats import chi2

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
        if x==y:
            s = 0
        elif y > x:
            y, x = x, y
            s = -(x*(log(x)-log(y))+y-x)*LOG10_E
        else:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
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
                write("%s\t%d\t%d\t%.5f\n" % (chrom,pre,pos[i],value[i]))
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

    cpdef dict make_pq_table ( self ):
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
                #if ('%.5f' % pre_v) != ('%.5f' % v):
                if abs(pre_v - v) > 1e-5: # precision is 5 digits
                    write( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre, p, pre_v ) )
                    pre_v = v
                    pre = p
            p = pos[ -1 ]
            # last one
            write( "%s\t%d\t%d\t%.5f\n" % ( chrom, pre, p, pre_v ) )
            
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

cdef class TwoConditionScores:
    cdef:
        dict data                       # dictionary for data of each chromosome
        dict datalength                 # length of data array of each chromosome
        float cond1_depth                 # seq depth in million of treatment
        float cond2_depth                  # seq depth in million of control
        float pseudocount                # the pseudocount used to calcuate LLR
        float cutoff
        object t1bdg, c1bdg, t2bdg, c2bdg
        dict pvalue_stat1, pvalue_stat2, pvalue_stat3
    
    def __init__ (self, t1bdg, c1bdg, t2bdg, c2bdg, float cond1_depth = 1.0, float cond2_depth = 1.0, float pseudocount = 0.01 ):
        self.data = {}           # for each chromosome, there is a l*4
                                 # matrix. First column: end position
                                 # of a region; Second: treatment
                                 # pileup; third: control pileup ;
                                 # forth: score ( can be
                                 # p/q-value/likelihood
                                 # ratio/fold-enrichment/subtraction
                                 # depending on -c setting)
        self.datalength = {}
        self.cond1_depth = cond1_depth
        self.cond2_depth = cond2_depth
        self.pseudocount = pseudocount
        self.pvalue_stat1 = {}
        self.pvalue_stat2 = {}
        self.t1bdg = t1bdg
        self.c1bdg = c1bdg
        self.t2bdg = t2bdg
        self.c2bdg = c2bdg

    cpdef set_pseudocount( self, float pseudocount ):
        self.pseudocount = pseudocount
        
    cpdef build ( self ):
        cdef:
            set common_chrs
            str chrname
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

    def get_common_chrs ( self ):
        t1chrs = self.t1bdg.get_chr_names()
        c1chrs = self.c1bdg.get_chr_names()
        t2chrs = self.t2bdg.get_chr_names()
        c2chrs = self.c2bdg.get_chr_names()        
        common = reduce(lambda x,y:x.intersection(y), (t1chrs,c1chrs,t2chrs,c2chrs))
        return common

    cdef add_chromosome ( self, str chrom, int chrom_max_len ):
        """
        chrom: chromosome name
        chrom_max_len: maximum number of data points in this chromosome
        
        """
        if not self.data.has_key(chrom):
            self.data[chrom] = [ np.zeros( chrom_max_len, dtype="int32" ), # pos
                                 np.zeros( chrom_max_len, dtype="float32" ), # LLR t1 vs c1
                                 np.zeros( chrom_max_len, dtype="float32" ), # LLR t2 vs c2
                                 np.zeros( chrom_max_len, dtype="float32" )] # LLR t1 vs t2
                                 #np.zeros( chrom_max_len, dtype="float32" )] # LLR t2 vs t1
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
        c = self.data[chromosome]
        c[0][ i ] = endpos
        c[1][ i ] = logLR( (t1+self.pseudocount)*self.cond1_depth, (c1+self.pseudocount)*self.cond1_depth )
        c[2][ i ] = logLR( (t2+self.pseudocount)*self.cond2_depth, (c2+self.pseudocount)*self.cond2_depth )
        c[3][ i ] = logLR( (t1+self.pseudocount)*self.cond1_depth, (t2+self.pseudocount)*self.cond2_depth )
        #c[4][ i ] = logLR( t2+self.pseudocount, t1+self.pseudocount )
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
    
        
    cpdef compute_all_pvalues( self ):
        """Compute -log_{10}(pvalue)
        """
        cdef:
            np.ndarray[np.float32_t] v1, v2, v3
            np.ndarray[np.int32_t] pos
            long l, i, prev_pos
            str chrom
        rv = chi2(1) # a chi-squared distribution with one degree of freedom
        sf = rv.sf
        log10 = np.log10
        for chrom in self.data.keys():
            prev_pos = 0
            pos = self.data[chrom][0]
            l = self.datalength[chrom]
            self.data[chrom][1] = -log10(sf(2 * self.data[chrom][1] / LOG10_E)).astype('float32')
            self.data[chrom][2] = -log10(sf(2 * self.data[chrom][2] / LOG10_E)).astype('float32')
            self.data[chrom][3] = -log10(sf(2 * np.abs(self.data[chrom][3]) / LOG10_E)).astype('float32')
            v1 = self.data[chrom][1]
            v2 = self.data[chrom][2]
            v3 = self.data[chrom][3]
            for i in range(l):
                try:
                    self.pvalue_stat1[v1[ i ]] += pos[ i ] - prev_pos
                except:
                    self.pvalue_stat1[v1[ i ]] = pos[ i ] - prev_pos
                try:
                    self.pvalue_stat2[v2[ i ]] += pos[ i ] - prev_pos
                except:
                    self.pvalue_stat2[v2[ i ]] = pos[ i ] - prev_pos
                prev_pos = pos[ i ]
            
    cpdef compute_track_qvalues ( self ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            dict pqtable
            tuple pqtables
            long i,l,j
            double k
            str chrom
            np.ndarray p, c, v
            
        # make pqtables
        # not sure if we should pool all the p-values or treat each comparison
        # separately, but doing the latter for now
        pqtables = self.make_pq_tables()
        
        # convert p to q

        # convert pvalue2qvalue to a simple dict based on khash
        # khash has big advantage while checking keys for millions of times.
        for table_index in range(2):
            pqtable = pqtables[table_index]
            s_p2q = Float64HashTable()
            for k in pqtable.keys():
                s_p2q.set_item(k,pqtable[k])
    
            g = s_p2q.get_item
            
            for chrom in self.data.keys():
                v = self.data[chrom][table_index + 1]
                l = self.datalength[chrom]
                for i in range(l):
                    v[ i ] =  g( v[ i ])
        
        return

    cpdef compute_category_qvalues ( self, dict category ):
        """Compute -log_{10}(qvalue)
        """
        cdef:
            dict value_dict, pqtable
            list unique_values
            long i,l,j
            np.ndarray p, c, v, data
            long pre_l
            double pre_v, unique_v, q, pre_q
            long N, k
            double f, pvalue
        
        # make pqtable for category
        value_dict = {}
        for chrom in self.data.keys():
            if category[chrom].size == 0: continue

            pos = self.data[chrom][0]
            data = self.data[chrom][3]
            for j in range(category[chrom].size):
                i = category[chrom][j]
                try:
                    value_dict[data[i]] += pos[i + 1] - pos[i]
                except IndexError:
                    if not value_dict.has_key(data[i]):
                        value_dict[data[i]] = 0
                except KeyError:
                    value_dict[data[i]] = pos[i + 1] - pos[i]

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
        
        qvalues = {}
        for chrom in self.data.keys():
            v = self.data[chrom][3][category[chrom]]
            qvalues[chrom] = v.copy()
            for i in range(v.size):
                qvalues[chrom][ i ] =  g( v[ i ])
        
        return qvalues

    # borrowed from CombinedTwoTrack
    cdef tuple make_pq_tables ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-log10pvalue:(-log10qvalue,rank,basepairs)} relationships.
        """
        cdef:
            long pre_l, l, i
            double pre_v, v, q, pre_q
            long N, k
            double f
            str chrom
            dict pvalue2qvalue1, pvalue2qvalue2
            dict value_dict
            list unique_values

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

    cpdef write_matrix ( self, fhd, str name, str description ):
        """Write all data to fhd into five columns Format:

        col1: chr_start_end
        col2: t1 vs c1
        col3: t2 vs c2
        col4: t1 vs t2
        col5: t2 vs t1

        fhd: a filehandler to save the matrix.

        """
        cdef:
            str chrom
            int l, pre, i, p 
            float v1, v2, v3, v4
            np.ndarray pos, value

        write = fhd.write

        chrs = self.get_chr_names()
        for chrom in chrs:
            pos = self.data[ chrom ][ 0 ]
            value = self.data[ chrom ][ column ]
            l = self.datalength[ chrom ]
            pre = 0
            if pos.shape[ 0 ] == 0: continue # skip if there's no data
            for i in range( 0, l ):
                v1 = self.data[ i ][ 1 ]
                v2 = self.data[ i ][ 2 ]
                v3 = self.data[ i ][ 3 ]
                v4 = self.data[ i ][ 4 ]                
                p = pos[ i ]
                write( "%s:%d_%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % ( chrom, pre, p, v1, v2, v3, v4 ) )
                pre = p
            
        return True

    cpdef call_peaks (self, float cutoff=0.05, int min_length=200,
                      int max_gap=0, bool call_summits=False):
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
            int i
            str chrom
            np.ndarray pos, sample, control, value, above_cutoff, \
                       above_cutoff_v, above_cutoff_endpos, \
                       above_cutoff_startpos, above_cutoff_sv
            np.ndarray[np.float32_t] cat1_qvalues, cat2_qvalues, cat3_qvalues
            list peak_content
        
        chrs  = self.get_chr_names()
        cat1_peaks = PeakIO()       # dictionary to save peaks significant at condition 1
        cat2_peaks = PeakIO()       # dictionary to save peaks significant at condition 2
        cat3_peaks = PeakIO()       # dictionary to save peaks significant in both conditions
        #cat4_peaks = PeakIO()       # dictionary to save peaks        

        self.cutoff = cutoff
        logcutoff = -log10(cutoff)
        cat1 = {}
        cat2 = {}
        cat3 = {}
        cat1_qvalues_by_chrom = {}
        cat2_qvalues_by_chrom = {}
        cat3_qvalues_by_chrom = {}
        for chrom in chrs:
            t1vsc1 = self.data[chrom][ 1 ]
            t2vsc2 = self.data[chrom][ 2 ]            
            and_ = np.logical_and
            not_ = np.logical_not
            cond1_sig = t1vsc1 >= logcutoff
            cond2_sig = t2vsc2 >= logcutoff
            # indices where score is above cutoff
            cat1[chrom] = np.where(and_(cond1_sig,       not_(cond2_sig)))[0]
            cat2[chrom] = np.where(and_(not_(cond1_sig), cond2_sig))[0]
            cat3[chrom] = np.where(and_(cond1_sig,       cond2_sig))[0]

        # qvalue conversion for each category to improve statistical power
        cat1_qvalues_by_chrom = self.compute_category_qvalues(cat1)
        cat2_qvalues_by_chrom = self.compute_category_qvalues(cat2)
        cat3_qvalues_by_chrom = self.compute_category_qvalues(cat3)

        for chrom in chrs:
            peak_content = []           # to store points above cutoff

            pos = self.data[chrom][ 0 ]
            t1vst2 = self.data[chrom][ 3 ]

            cat1_qvalues = cat1_qvalues_by_chrom[chrom]
            cat2_qvalues = cat2_qvalues_by_chrom[chrom]
            cat3_qvalues = cat3_qvalues_by_chrom[chrom]


            where_cat1_above_cutoff = cat1_qvalues >= logcutoff
            where_cat2_above_cutoff = cat2_qvalues >= logcutoff
            where_cat3_above_cutoff = cat3_qvalues >= logcutoff
            cat1_above_cutoff = cat1[chrom][where_cat1_above_cutoff]
            cat2_above_cutoff = cat2[chrom][where_cat2_above_cutoff]
            cat3_above_cutoff = cat3[chrom][where_cat2_above_cutoff]
            cat1_above_cutoff_endpos = pos[cat1_above_cutoff] # end positions of regions where score is above cutoff
            cat1_above_cutoff_startpos = pos[cat1_above_cutoff-1] # start positions of regions where score is above cutoff
            cat2_above_cutoff_endpos = pos[cat2_above_cutoff] # end positions of regions where score is above cutoff
            cat2_above_cutoff_startpos = pos[cat2_above_cutoff-1] # start positions of regions where score is above cutoff
            cat3_above_cutoff_endpos = pos[cat3_above_cutoff] # end positions of regions where score is above cutoff
            cat3_above_cutoff_startpos = pos[cat3_above_cutoff-1] # start positions of regions where score is above cutoff
            #cat4_above_cutoff_endpos = pos[cat4_above_cutoff] # end positions of regions where score is above cutoff
            #cat4_above_cutoff_startpos = pos[cat4_above_cutoff-1] # start positions of regions where score is above cutoff            

            # should we actually save these statistics???

            if cat1_above_cutoff_startpos.size > 0:
                qscores = t1vst2[cat1_above_cutoff]
                peak_content = []
                # cat 1 is not empty
                if cat1_above_cutoff[0] == 0:
                    # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
                    cat1_above_cutoff_startpos[0] = 0
                # first bit of region above cutoff
                peak_content.append( (cat1_above_cutoff_startpos[0], cat1_above_cutoff_endpos[0], qscores[0]) )
                for i in range( 1, cat1_above_cutoff_startpos.size ):
                    if cat1_above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                        # append
                        peak_content.append( (cat1_above_cutoff_startpos[i], cat1_above_cutoff_endpos[i], qscores[i] ) )
                    else:
                        # close
                        if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:
                            best_qscore = min(map(itemgetter(2), peak_content))
                            cat1_peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                            summit = -1, peak_score  = 0, pileup = 0, pscore = 0, 
                                            fold_change = 0, qscore = best_qscore,
                                            )
                        peak_content = [(cat1_above_cutoff_startpos[i], cat1_above_cutoff_endpos[i], qscores[i]),]
                # save the last peak
                if not peak_content:
                    continue
                else:
                    if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:                    
                        best_qscore = min(map(itemgetter(2), peak_content))
                        cat1_peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                        summit = -1, peak_score  = 0, pileup = 0, pscore = 0, 
                                        fold_change = 0, qscore = best_qscore, 
                                        )
            if cat2_above_cutoff_startpos.size > 0:
                qscores = t1vst2[cat2_above_cutoff]
                peak_content = []
                # cat 1 is not empty
                if cat2_above_cutoff[0] == 0:
                    # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
                    cat2_above_cutoff_startpos[0] = 0
                # first bit of region above cutoff
                peak_content.append( (cat2_above_cutoff_startpos[0], cat2_above_cutoff_endpos[0], qscores[0]) )
                for i in range( 1, cat2_above_cutoff_startpos.size ):
                    if cat2_above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                        # append
                        peak_content.append( (cat2_above_cutoff_startpos[i], cat2_above_cutoff_endpos[i], qscores[i] ) )
                    else:
                        # close
                        if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:                        
                            best_qscore = min(map(itemgetter(2), peak_content))
                            cat2_peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                            summit = -1, peak_score  = 0, pileup = 0, pscore = 0, 
                                            fold_change = 0, qscore = best_qscore, 
                                            )
                        peak_content = [(cat2_above_cutoff_startpos[i], cat2_above_cutoff_endpos[i], qscores[i]),]
                # save the last peak
                if not peak_content:
                    continue
                else:
                    if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:                    
                        best_qscore = min(map(itemgetter(2), peak_content))
                        cat2_peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                        summit = -1, peak_score  = 0, pileup = 0, pscore = 0, 
                                        fold_change = 0, qscore = best_qscore, 
                                        )

            if cat3_above_cutoff_startpos.size > 0:
                qscores = t1vst2[cat3_above_cutoff]
                peak_content = []
                # cat 1 is not empty
                if cat3_above_cutoff[0] == 0:
                    # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
                    cat3_above_cutoff_startpos[0] = 0
                # first bit of region above cutoff
                peak_content.append( (cat3_above_cutoff_startpos[0], cat3_above_cutoff_endpos[0], qscores[0]) )
                for i in range( 1, cat3_above_cutoff_startpos.size ):
                    if cat3_above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                        # append
                        peak_content.append( (cat3_above_cutoff_startpos[i], cat3_above_cutoff_endpos[i], qscores[i] ) )
                    else:
                        # close
                        if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:                        
                            best_qscore = min(map(itemgetter(2), peak_content))
                            cat3_peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                            summit = -1, peak_score  = 0, pileup = 0, pscore = 0, 
                                            fold_change = 0, qscore = best_qscore, 
                                            )
                        peak_content = [(cat3_above_cutoff_startpos[i], cat3_above_cutoff_endpos[i], qscores[i]),]
                # save the last peak
                if not peak_content:
                    continue
                else:
                    if peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ] >= min_length:
                        best_qscore = min(map(itemgetter(2), peak_content))
                        cat3_peaks.add( chrom, peak_content[0][0], peak_content[-1][1],
                                        summit = -1, peak_score  = 0, pileup = 0, pscore = 0, 
                                        fold_change = 0, qscore = best_qscore, 
                                        )

        return cat1_peaks, cat2_peaks, cat3_peaks

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

cdef class DiffScoreTrackI:
    cdef:
        dict pos, t1, c1, t2, c2, tvsc1, tvsc2, t1vs2, tlogLR
        dict in_peaks, where_peaks, diff_peaks
        dict datalength                 # length of data array of each chromosome
        float cond1_depth                 # seq depth in million of treatment
        float cond2_depth                  # seq depth in million of control
        float pseudocount                # the pseudocount used to calcuate LLR
        float cutoff
        dict pvalue_stat1, pvalue_stat2
        str track_scoring_method
    
    def __init__ (self, t1bdg, c1bdg, t2bdg, c2bdg, float cond1_depth = 1.0, float cond2_depth = 1.0, float pseudocount = 0.01 ):
        self.pos = {}
        self.t1 = {}
        self.c1 = {}
        self.tvsc1 = {}
        self.t2 = {}
        self.c2 = {}
        self.tvsc2 = {}
        self.t1vs2 = {}
        self.tlogLR = {}
        self.in_peaks = {}
        self.where_peaks = {}
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
            
    cpdef set_pseudocount( self, float pseudocount ):
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
            self.t1vs2[chrom] = np.ones( chrom_max_len, dtype="float32" ) # c1
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
        self.t1[chromosome][ i ] = t1 *self.cond1_depth
        self.c1[chromosome][ i ] = c2 *self.cond1_depth
        self.t2[chromosome][ i ] = t2 *self.cond2_depth
        self.c2[chromosome][ i ] = c2 *self.cond2_depth
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
        
        for chrom in self.pos.keys():
            prev_pos = 0
            pos = self.pos[chrom]
            p = self.t1[chrom]
            c = self.c1[chrom]
            v = self.tvsc1[chrom]
            l = self.datalength[chrom]
            for i in range(l):
                v[ i ] =  get_pscore( int(p[ i ]) , c[ i ] )
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
                v[ i ] =  get_pscore( int(p[ i ]) , c[ i ] )
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
            above_cutoff_startpos = pos[above_cutoff-1] # start positions of regions where score is above cutoff

            print "Regions > cutoff: %d" % above_cutoff.size
            first_i = -1
            if above_cutoff.size > 1:
                if above_cutoff_startpos[0] > above_cutoff_endpos[0]: above_cutoff_startpos[0] = 0
                first_i = 0
                for i in range(1, above_cutoff_startpos.size):
                    if first_i == -1:
                        first_i = i
                    if not above_cutoff_startpos[i] - above_cutoff_endpos[i - 1] <= max_gap:
                        length = above_cutoff_endpos[i - 1] - above_cutoff_startpos[first_i]
                        if length >= min_length:
                            in_peaks[first_i:i] = True
                        first_i = -1
            
            if not first_i == -1:
                length = above_cutoff_endpos[i - 1] - above_cutoff_startpos[first_i]
                if length >= min_length:
                    in_peaks[first_i:(i+1)] = True
            
            value = self.tvsc2[chrom]
            above_cutoff = np.nonzero( value >= cutoff )[0] # indices where score is above cutoff
            above_cutoff_endpos = pos[above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = pos[above_cutoff-1] # start positions of regions where score is above cutoff

            first_i = -1
            if above_cutoff.size > 1:
                if above_cutoff_startpos[0] > above_cutoff_endpos[0]: above_cutoff_startpos[0] = 0
                first_i = 0
                for i in range(1, above_cutoff_startpos.size):
                    if first_i == -1:
                        first_i = i
                    if not above_cutoff_startpos[i] - above_cutoff_endpos[i - 1] <= max_gap:
                        length = above_cutoff_endpos[i - 1] - above_cutoff_startpos[first_i]
                        if length >= min_length:
                            in_peaks[first_i:i] = True # not including i
                        first_i = -1
            
            if not first_i == -1:
                length = above_cutoff_endpos[i - 1] - above_cutoff_startpos[first_i]
                if length >= min_length:
                    in_peaks[first_i:(i+1)] = True # including i

            self.in_peaks[chrom] = in_peaks
            self.where_peaks[chrom] = np.where(self.in_peaks[chrom])[0]
        
        return

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
            int i
            float cond1_d, cond2_d, pseudocount
            np.ndarray[np.float32_t] t1, t2, v
        rv = chi2(1) # a chi-squared distribution with one degree of freedom
        sf = rv.sf
        log10 = np.log10
        cond1_d = self.cond1_depth
        cond2_d = self.cond2_depth
        pseudocount = self.pseudocount
        for chrom in self.pos.keys():
            t1 = self.t1[chrom]
            t2 = self.t2[chrom]
            v = self.tlogLR[chrom]
            for i in range(self.pos[chrom].size): 
                v[i] = logLR(cond1_d * (t1[i] + pseudocount),
                             cond2_d * (t2[i] + pseudocount))
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
                    if (prev_i + 1) == i:
                        value_dict[stat[i]] += pos[i] - pos[prev_i]
                    else:
                        value_dict[stat[i]] += pos[i] - pos[i - 1]
                except IndexError:
                    if not value_dict.has_key(stat[i]):
                        value_dict[stat[i]] = 0
                except KeyError:
                    if (prev_i + 1) == i:
                        value_dict[stat[i]] = pos[i] - pos[prev_i]
                    else:
                        value_dict[stat[i]] += pos[i] - pos[i - 1]
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
        
        qvalues = {}
        for chrom in self.t1vs2.keys():
            v = self.t1vs2[chrom][self.where_peaks[chrom]]
            qvalues[chrom] = v.copy()
            for i in range(v.size):
                qvalues[chrom][ i ] =  g(v[ i ])
        
        return qvalues

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

    cpdef write_matrix ( self, fhd, str name, str description ):
        """Write all data to fhd into five columns Format:

        col1: chr_start_end
        col2: t1 vs c1
        col3: t2 vs c2
        col4: t1 vs t2
        col5: t2 vs t1

        fhd: a filehandler to save the matrix.

        """
        cdef:
            str chrom
            int l, pre, i, p 
            float v1, v2, v3, v4
            np.ndarray pos, value

        write = fhd.write

        chrs = self.get_chr_names()
        for chrom in chrs:
            pos = self.data[ chrom ][ 0 ]
            value = self.data[ chrom ][ column ]
            l = self.datalength[ chrom ]
            pre = 0
            if pos.shape[ 0 ] == 0: continue # skip if there's no data
            for i in range( 0, l ):
                v1 = self.data[ i ][ 1 ]
                v2 = self.data[ i ][ 2 ]
                v3 = self.data[ i ][ 3 ]
                v4 = self.data[ i ][ 4 ]                
                p = pos[ i ]
                write( "%s:%d_%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % ( chrom, pre, p, v1, v2, v3, v4 ) )
                pre = p
            
        return True

    cpdef call_diff_peaks (self, float cutoff=2.0, int min_length=50,
                      int max_gap=20, str score_method='q'):
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
            np.ndarray pos, sample, control, value, above_cutoff, \
                       above_cutoff_v, above_cutoff_endpos, \
                       above_cutoff_startpos, above_cutoff_sv
            np.ndarray[np.int64_t, ndim=2] diff_peaks
            list peak_content
        
        chrs  = self.get_chr_names()

        self.cutoff = cutoff

        # qvalue conversion for each category to improve statistical power
        diff_qvalues_by_chrom = self.compute_diff_pvalues()
        diff_qvalues_by_chrom = self.compute_diff_qvalues()
        if not (score_method == 'q' or score_method == 'p'):
            raise NotImplementedError

        for chrom in chrs:
            peak_content = []           # to store points above cutoff

            pos = self.pos[chrom]
            qpos = self.where_peaks[chrom]
            t1vst2 = self.t1vs2[chrom]

            if score_method == 'q':
                diff_qvalues = diff_qvalues_by_chrom[chrom]
            elif score_method == 'p':
                diff_qvalues = self.t1vs2[chrom][self.in_peaks[chrom]]
            else:
                raise NotImplementedError

            is_above_cutoff = diff_qvalues >= cutoff
            above_cutoff = qpos[np.nonzero(is_above_cutoff)[0]]
            above_cutoff_endpos = pos[above_cutoff] # end positions of regions where score is above cutoff
            above_cutoff_startpos = pos[above_cutoff-1] # start positions of regions where score is above cutoff

            # save original indices so that we can output everything
            first_i = -1
            n_diff_peaks = 0
            diff_peaks = np.ndarray((above_cutoff.size, 2), dtype='int64')
            print above_cutoff.size
            if above_cutoff.size > 1:
                if above_cutoff_startpos[0] > above_cutoff_startpos[1]:
                    above_cutoff_startpos[0] = 0
                first_i = 0
                for i in range(1, above_cutoff_startpos.size):
                    if first_i == -1:
                        first_i = i
                    if not above_cutoff_startpos[i] - above_cutoff_endpos[i - 1] <= max_gap:
                        length = above_cutoff_endpos[i - 1] - above_cutoff_startpos[first_i]
                        if length >= min_length:
                            diff_peaks[n_diff_peaks,0] = first_i
                            diff_peaks[n_diff_peaks,1] = i - 1
                            n_diff_peaks += 1
                        first_i = -1
           
            if not first_i == -1:
                length = above_cutoff_endpos[i - 1] - above_cutoff_startpos[first_i]
                if length >= min_length:
                    diff_peaks[n_diff_peaks,0] = first_i
                    diff_peaks[n_diff_peaks,1] = i
                    n_diff_peaks += 1

            diff_peaks.resize((n_diff_peaks, 2), refcheck=False)
            self.diff_peaks[chrom] = diff_peaks
            print chrom, n_diff_peaks
        return
    
    cpdef print_some_peaks(self, int max = 10):
        """for testing only
        """
        cdef:
            int start, end, i
            str chrom
        for chrom in self.diff_peaks.keys():
            i = 0
            for end in self.where_peaks[chrom]:
                i += 1
                if i == max: break
                print '%s\t%d' % (chrom, self.pos[chrom][end])
        return
    
    cpdef print_diff_peaks(self):
        """ for testing only
        """
        cdef:
            int start, end
            str chrom
        for chrom in self.diff_peaks.keys():
            for start,end in self.diff_peaks[chrom]:
                print '%s\t%d\t%d' % (chrom, start, end)
        return

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

