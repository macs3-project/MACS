# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-12-04 23:10:35 Tao Liu>

"""Module for SAPPER PosReadsInfo class.

Copyright (c) 2017 Tao Liu <tliu4@buffalo.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
from MACS3.Signal.VariantStat import CalModel_Homo, CalModel_Heter_noAS, CalModel_Heter_AS, calculate_GQ, calculate_GQ_heterASsig
from MACS3.Signal.Prob import binomial_cdf
from MACS3.Signal.PeakVariants import Variant

from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t, float32_t

LN10 = 2.3025850929940458

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    void free(void *ptr)
    int strcmp(char *a, char *b)
    char * strcpy(char *a, char *b)
    long atol(char *bytes)
    int atoi(char *bytes)

# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <tliu4@buffalo.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class PosReadsInfo:
    cdef:
        long ref_pos
        bytes ref_allele
        bytes alt_allele
        bool filterout          # if true, do not output

        dict bq_set_T     #{A:[], C:[], G:[], T:[], N:[]} for treatment
        dict bq_set_C
        dict n_reads_T    #{A:[], C:[], G:[], T:[], N:[]} for treatment
        dict n_reads_C
        dict n_reads

        list n_strand #[{A:[], C:[], G:[], T:[], N:[]},{A:[], C:[], G:[], T:[], N:[]}] for total appearance on plus strand and minus strand for ChIP sample only
        dict n_tips # count of nt appearing at tips

        bytes top1allele
        bytes top2allele
        float top12alleles_ratio

        double lnL_homo_major,lnL_heter_AS,lnL_heter_noAS,lnL_homo_minor
        double BIC_homo_major,BIC_heter_AS,BIC_heter_noAS,BIC_homo_minor
        double PL_00, PL_01, PL_11
        double deltaBIC
        int heter_noAS_kc, heter_noAS_ki
        int heter_AS_kc, heter_AS_ki
        double heter_AS_alleleratio

        int GQ_homo_major,GQ_heter_noAS,GQ_heter_AS  #phred scale of prob by standard formular
        int GQ_heter_ASsig #phred scale of prob, to measure the difference between AS and noAS

        double GQ
        
        str GT
        str type
        str mutation_type       # SNV or Insertion or Deletion

        bool hasfermiinfor #if no fermi bam overlap in the position, false; if fermi bam in the position GT: N, false; if anyone of top2allele is not in fermi GT NTs, false;
        bytearray fermiNTs # 

    def __cinit__ ( self ):
        self.filterout = False
        self.GQ = 0
        self.GT = "unsure"
        self.alt_allele = b'.'
        
    def __init__ ( self, long ref_pos, bytes ref_allele ):
        self.ref_pos = ref_pos
        self.ref_allele = ref_allele
        self.bq_set_T = { ref_allele:[],b'A':[], b'C':[], b'G':[], b'T':[], b'N':[], b'*':[] }
        self.bq_set_C = { ref_allele:[],b'A':[], b'C':[], b'G':[], b'T':[], b'N':[], b'*':[] }
        self.n_reads_T = { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0, b'*':0 }
        self.n_reads_C = { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0, b'*':0 }
        self.n_reads =  { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0, b'*':0 }
        self.n_strand = [ { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0, b'*':0 }, { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0, b'*':0 } ]
        self.n_tips = { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0, b'*':0 }
        

    #cpdef void merge ( self, PosReadsInfo PRI2 ):
    #    """Merge two PRIs. No check available.
    #
    #    """
    #    assert self.ref_pos == PRI2.ref_pos
    #    assert self.ref_allele == PRI2.ref_allele
    #    for b in set( self.n_reads.keys() ).union( set( PRI2.n_reads.keys() ) ):
    #        self.bq_set_T[ b ] = self.bq_set_T.get( b, []).extend( PRI2.bq_set_T.get( b, [] ) )
    #        self.bq_set_C[ b ] = self.bq_set_C.get( b, []).extend( PRI2.bq_set_C.get( b, [] ) )
    #        self.n_reads_T[ b ] = self.n_reads_T.get( b, 0) + PRI2.n_reads_T.get( b, 0 )
    #        self.n_reads_C[ b ] = self.n_reads_C.get( b, 0) + PRI2.n_reads_C.get( b, 0 )
    #        self.n_reads[ b ] = self.n_reads.get( b, 0) + PRI2.n_reads.get( b, 0 )
    #    return

    def __getstate__ ( self ):
        return ( self.ref_pos, self.ref_allele, self.alt_allele, self.filterout,
                 self.bq_set_T, self.bq_set_C, self.n_reads_T, self.n_reads_C, self.n_reads, self.n_strand, self.n_tips,
                 self.top1allele, self.top2allele, self.top12alleles_ratio,
                 self.lnL_homo_major, self.lnL_heter_AS, self.lnL_heter_noAS, self.lnL_homo_minor,
                 self.BIC_homo_major, self.BIC_heter_AS, self.BIC_heter_noAS, self.BIC_homo_minor,
                 self.heter_noAS_kc, self.heter_noAS_ki,
                 self.heter_AS_kc, self.heter_AS_ki,
                 self.heter_AS_alleleratio,
                 self.GQ_homo_major, self.GQ_heter_noAS, self.GQ_heter_AS,
                 self.GQ_heter_ASsig,
                 self.GQ,
                 self.GT,
                 self.type,
                 self.hasfermiinfor,
                 self.fermiNTs )

    def __setstate__ ( self, state ):
        ( self.ref_pos, self.ref_allele, self.alt_allele, self.filterout,
          self.bq_set_T, self.bq_set_C, self.n_reads_T, self.n_reads_C, self.n_reads, self.n_strand, self.n_tips,
          self.top1allele, self.top2allele, self.top12alleles_ratio,
          self.lnL_homo_major, self.lnL_heter_AS, self.lnL_heter_noAS, self.lnL_homo_minor,
          self.BIC_homo_major, self.BIC_heter_AS, self.BIC_heter_noAS, self.BIC_homo_minor,
          self.heter_noAS_kc, self.heter_noAS_ki,
          self.heter_AS_kc, self.heter_AS_ki,
          self.heter_AS_alleleratio,
          self.GQ_homo_major, self.GQ_heter_noAS, self.GQ_heter_AS,
          self.GQ_heter_ASsig,
          self.GQ,
          self.GT,
          self.type,
          self.hasfermiinfor,
          self.fermiNTs ) = state

    cpdef bool filterflag ( self ):
        return self.filterout

    cpdef void apply_GQ_cutoff ( self, int min_homo_GQ = 50, int min_heter_GQ = 100 ):
        if self.filterout:
            return
        if self.type.startswith('homo') and self.GQ < min_homo_GQ:
            self.filterout = True
        elif self.type.startswith('heter') and self.GQ < min_heter_GQ:
            self.filterout = True
        return

    cpdef void apply_deltaBIC_cutoff ( self, float min_delta_BIC = 10 ):
        if self.filterout:
            return
        if self.deltaBIC < min_delta_BIC:
            self.filterout = True
        return

    cpdef void add_T ( self, int read_index, bytes read_allele, int read_bq, int strand, bool tip, int Q=20 ):
        """ Strand 0: plus, 1: minus

        Q is the quality cutoff. By default, only consider Q20 or read_bq > 20.
        """
        if read_bq <= Q:
            return
        if not self.n_reads.has_key( read_allele ):
            self.bq_set_T[read_allele] = []
            self.bq_set_C[read_allele] = []
            self.n_reads_T[read_allele] = 0
            self.n_reads_C[read_allele] = 0
            self.n_reads[read_allele] = 0
            self.n_strand[ 0 ][ read_allele ] = 0
            self.n_strand[ 1 ][ read_allele ] = 0
            self.n_tips[read_allele] = 0
        self.bq_set_T[read_allele].append( read_bq )
        self.n_reads_T[ read_allele ] += 1
        self.n_reads[ read_allele ] += 1
        self.n_strand[ strand ][ read_allele ] += 1
        if tip: self.n_tips[ read_allele ] += 1

    cpdef void add_C ( self, int read_index, bytes read_allele, int read_bq, int strand, int Q=20 ):
        if read_bq <= Q:
            return
        if not self.n_reads.has_key( read_allele ):
            self.bq_set_T[read_allele] = []
            self.bq_set_C[read_allele] = []
            self.n_reads_T[read_allele] = 0
            self.n_reads_C[read_allele] = 0
            self.n_reads[read_allele] = 0
            self.n_strand[ 0 ][ read_allele ] = 0
            self.n_strand[ 1 ][ read_allele ] = 0
            self.n_tips[read_allele] = 0
        self.bq_set_C[read_allele].append( read_bq )
        self.n_reads_C[ read_allele ] += 1
        self.n_reads[ read_allele ] += 1
        #self.n_strand[ strand ][ read_allele ] += 1

    cpdef int raw_read_depth ( self, str opt = "all" ):
        if opt == "all":
            return sum( self.n_reads.values() )
        elif opt == "T":
            return sum( self.n_reads_T.values() )
        elif opt == "C":
            return sum( self.n_reads_C.values() )
        else:
            raise Exception( "opt should be either 'all', 'T' or 'C'." )

    cpdef void update_top_alleles ( self, float min_top12alleles_ratio = 0.8, int min_altallele_count = 2, float max_allowed_ar = 0.95 ):
    #cpdef update_top_alleles ( self, float min_top12alleles_ratio = 0.8 ):
        """Identify top1 and top2 NT.  the ratio of (top1+top2)/total
        """
        cdef:
            float r

        [self.top1allele, self.top2allele] = sorted(self.n_reads, key=self.n_reads_T.get, reverse=True)[:2]

        # if top2 allele count in ChIP is lower than
        # min_altallele_count, or when allele ratio top1/(top1+top2)
        # is larger than max_allowed_ar in ChIP, we won't consider
        # this allele at all.  we set values of top2 allele in
        # dictionaries to [] and ignore top2 allele entirely.
 
        # max(self.n_strand[ 0 ][ self.top2allele ], self.n_strand[ 1 ][ self.top2allele ]) < min_altallele_count
        #if self.ref_pos == 52608504:
        #    print self.ref_pos, self.n_reads_T[ self.top1allele ], self.n_reads_T[ self.top2allele ], self.n_reads_C[ self.top1allele ], self.n_reads_C[ self.top2allele ]
        if self.n_reads_T[ self.top1allele ] + self.n_reads_T[ self.top2allele ] == 0:
            self.filterout = True
            return

        if (len(self.top1allele)==1 and len(self.top2allele)==1) and ( self.top2allele != self.ref_allele and ( ( self.n_reads_T[ self.top2allele ] - self.n_tips[ self.top2allele ] ) < min_altallele_count ) or \
                self.n_reads_T[ self.top1allele ]/(self.n_reads_T[ self.top1allele ] + self.n_reads_T[ self.top2allele ]) > max_allowed_ar ):
            self.bq_set_T[ self.top2allele ] = []
            self.bq_set_C[ self.top2allele ] = []
            self.n_reads_T[ self.top2allele ] = 0
            self.n_reads_C[ self.top2allele ] = 0
            self.n_reads[ self.top2allele ] = 0
            self.n_tips[ self.top2allele ] = 0
            if ( self.top1allele != self.ref_allele and ( self.n_reads_T[ self.top1allele ] - self.n_tips[ self.top1allele ] ) < min_altallele_count ):
                self.bq_set_T[ self.top1allele ] = []
                self.bq_set_C[ self.top1allele ] = []
                self.n_reads_T[ self.top1allele ] = 0
                self.n_reads_C[ self.top1allele ] = 0
                self.n_reads[ self.top1allele ] = 0
                self.n_tips[ self.top1allele ] = 0

        if self.n_reads_T[ self.top1allele ] + self.n_reads_T[ self.top2allele ] == 0:
            self.filterout = True
            return

        self.top12alleles_ratio = ( self.n_reads[ self.top1allele ] + self.n_reads[ self.top2allele ] ) /  sum( self.n_reads.values() )
        if self.top12alleles_ratio < min_top12alleles_ratio:
            self.filterout = True
            return

        if self.top1allele == self.ref_allele and self.n_reads[ self.top2allele ] == 0:
            # This means this position only contains top1allele which is the ref_allele. So the GT must be 0/0
            self.type = "homo_ref"
            self.filterout = True
            return
        return

    cpdef void top12alleles ( self ):
        print ( self.ref_pos, self.ref_allele)
        print ("Top1allele",self.top1allele, "Treatment", self.bq_set_T[self.top1allele], "Control", self.bq_set_C[self.top1allele])
        print ("Top2allele",self.top2allele, "Treatment", self.bq_set_T[self.top2allele], "Control", self.bq_set_C[self.top2allele])
    
    cpdef void call_GT ( self, float max_allowed_ar = 0.99 ):
        """Require update_top_alleles being called.
        """
        cdef:
            np.ndarray[np.int32_t, ndim=1] top1_bq_T
            np.ndarray[np.int32_t, ndim=1] top2_bq_T
            np.ndarray[np.int32_t, ndim=1] top1_bq_C
            np.ndarray[np.int32_t, ndim=1] top2_bq_C
            int i
            list top1_bq_T_l
            list top2_bq_T_l
            list top1_bq_C_l
            list top2_bq_C_l
            list tmp_mutation_type
            bytes tmp_alt

        if self.filterout:
            return
        
        top1_bq_T = np.array( self.bq_set_T[ self.top1allele ], dtype="int32" )
        top2_bq_T = np.array( self.bq_set_T[ self.top2allele ], dtype="int32" )
        top1_bq_C = np.array( self.bq_set_C[ self.top1allele ], dtype="int32" )
        top2_bq_C = np.array( self.bq_set_C[ self.top2allele ], dtype="int32" )
        (self.lnL_homo_major, self.BIC_homo_major) = CalModel_Homo( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_homo_minor, self.BIC_homo_minor) = CalModel_Homo( top2_bq_T, top2_bq_C, top1_bq_T, top1_bq_C )
        (self.lnL_heter_noAS, self.BIC_heter_noAS) = CalModel_Heter_noAS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_heter_AS, self.BIC_heter_AS)     = CalModel_Heter_AS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C, max_allowed_ar )

        #if self.ref_pos == 71078525:
        #    print "---"
        #    print len( top1_bq_T ), len( top1_bq_C ), len( top2_bq_T ), len( top2_bq_C )
        #    print self.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS, self.lnL_heter_AS
        #    print self.BIC_homo_major, self.BIC_homo_minor, self.BIC_heter_noAS, self.BIC_heter_AS
             
        if self.top1allele != self.ref_allele and self.n_reads[ self.top2allele ] == 0:
            # in this case, there is no top2 nt (or socalled minor
            # allele) in either treatment or control, we should assume
            # it's a 1/1 genotype. We will take 1/1 if it passes BIC
            # test (deltaBIC >=2), and will skip this loci if it can't
            # pass the test.
            
            self.deltaBIC = min( self.BIC_heter_noAS, self.BIC_heter_AS, self.BIC_homo_minor ) - self.BIC_homo_major
            if self.deltaBIC < 2:
               self.filterout = True
               return

            self.type = "homo"
            self.GT = "1/1"

            self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
            self.PL_01 = -10.0 * max( self.lnL_heter_noAS, self.lnL_heter_AS ) / LN10
            self.PL_11 = -10.0 * self.lnL_homo_major / LN10

            self.PL_00 = max( 0, self.PL_00 - self.PL_11 )
            self.PL_01 = max( 0, self.PL_01 - self.PL_11 )
            self.PL_11 = 0

            self.GQ = min( self.PL_00, self.PL_01 )
            self.alt_allele = self.top1allele
        else:
            # assign GQ, GT, and type
            if self.ref_allele != self.top1allele and self.BIC_homo_major + 2 <= self.BIC_homo_minor and self.BIC_homo_major + 2 <= self.BIC_heter_noAS and self.BIC_homo_major + 2 <= self.BIC_heter_AS:
                self.type = "homo"
                self.deltaBIC = min( self.BIC_heter_noAS, self.BIC_heter_AS, self.BIC_homo_minor ) - self.BIC_homo_major
                self.GT = "1/1"
                self.alt_allele = self.top1allele

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * max( self.lnL_heter_noAS, self.lnL_heter_AS ) / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_11
                self.PL_01 = self.PL_01 - self.PL_11
                self.PL_11 = 0

                self.GQ = min( self.PL_00, self.PL_01 )
                
            elif self.BIC_heter_noAS + 2 <= self.BIC_homo_major and self.BIC_heter_noAS + 2 <= self.BIC_homo_minor and self.BIC_heter_noAS + 2 <= self.BIC_heter_AS :
                self.type = "heter_noAS"
                self.deltaBIC = min( self.BIC_homo_major, self.BIC_homo_minor ) - self.BIC_heter_noAS

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * self.lnL_heter_noAS / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_01
                self.PL_11 = self.PL_11 - self.PL_01
                self.PL_01 = 0

                self.GQ = min( self.PL_00, self.PL_11 )
                
            elif self.BIC_heter_AS + 2 <= self.BIC_homo_major and self.BIC_heter_AS + 2 <= self.BIC_homo_minor and self.BIC_heter_AS + 2 <= self.BIC_heter_noAS:
                self.type = "heter_AS"
                self.deltaBIC = min( self.BIC_homo_major, self.BIC_homo_minor ) - self.BIC_heter_AS

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * self.lnL_heter_AS / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_01
                self.PL_11 = self.PL_11 - self.PL_01
                self.PL_01 = 0

                self.GQ = min( self.PL_00, self.PL_11 )

            elif self.BIC_heter_AS + 2 <= self.BIC_homo_major and self.BIC_heter_AS + 2 <= self.BIC_homo_minor:
                # can't decide if it's noAS or AS
                self.type = "heter_unsure"
                self.deltaBIC = min( self.BIC_homo_major, self.BIC_homo_minor ) - max( self.BIC_heter_AS, self.BIC_heter_noAS )

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * max( self.lnL_heter_noAS, self.lnL_heter_AS ) / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_01
                self.PL_11 = self.PL_11 - self.PL_01
                self.PL_01 = 0

                self.GQ = min( self.PL_00, self.PL_11 )
                
            elif self.ref_allele == self.top1allele and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
                self.type = "homo_ref"
                # we do not calculate GQ if type is homo_ref
                self.GT = "0/0"
                self.filterout = True
            else:
                self.type="unsure"
                self.filterout = True

            if self.type.startswith( "heter" ):
                if self.ref_allele == self.top1allele:
                    self.alt_allele = self.top2allele
                    self.GT = "0/1"
                elif self.ref_allele == self.top2allele:
                    self.alt_allele = self.top1allele
                    self.GT = "0/1"
                else:
                    self.alt_allele = self.top1allele+b','+self.top2allele
                    self.GT = "1/2"
                # strand bias filter, uncomment following if wish to debug
                # calculate SB score
                #print "calculate SB score for ", self.ref_pos, "a/b/c/d:", self.n_strand[ 0 ][ self.top1allele ], self.n_strand[ 0 ][ self.top2allele ], self.n_strand[ 1 ][ self.top1allele ], self.n_strand[ 1 ][ self.top2allele ]
                #SBscore = self.SB_score_ChIP( self.n_strand[ 0 ][ self.top1allele ], self.n_strand[ 0 ][ self.top2allele ], self.n_strand[ 1 ][ self.top1allele ], self.n_strand[ 1 ][ self.top2allele ] )
                #SBscore = 0
                #if SBscore >= 1:
                #    print "disgard variant at", self.ref_pos, "type", self.type
                #    self.filterout = True
                
                # if self.ref_allele == self.top1allele：
                #     self.n_strand[ 0 ][ self.top1allele ] + self.n_strand[ 1 ][ self.top1allele ]
                #     if and self.n_strand[ 0 ][ self.top2allele ] == 0 or self.n_strand[ 1 ][ self.top2allele ] == 0:
                #         self.filterout = True
                #         print self.ref_pos


            # self.deltaBIC = self.deltaBIC

        tmp_mutation_type = []
        for tmp_alt in self.alt_allele.split(b','):
            if tmp_alt == b'*':
                tmp_mutation_type.append( "Deletion" )
            elif len( tmp_alt ) > 1:
                tmp_mutation_type.append( "Insertion" )
            else:
                tmp_mutation_type.append( "SNV" )
        self.mutation_type = ",".join( tmp_mutation_type )
        return

    cdef float SB_score_ChIP( self, int a, int b, int c, int d ):
        """ calculate score for filtering variants with strange strand biases.

        a: top1/major allele plus strand
        b: top2/minor allele plus strand
        c: top1/major allele minus strand
        d: top2/minor allele minus strand

        Return a float value so that if this value >= 1, the variant will be filtered out.
        """
        cdef:
            float score
            double p
            double p1_l, p1_r
            double p2_l, p2_r
            double top2_sb, top1_sb

        if a+b == 0 or c+d == 0:
            # if major allele and minor allele both bias to the same strand, allow it
            return 0.0

        # Rule:
        # if there is bias in top2 allele then bias in top1 allele should not be significantly smaller than it.
        # or there is no significant bias (0.5) in top2 allele.
        
        #print a, b, c, d
        p1_l = binomial_cdf( a, (a+c), 0.5, lower=True )      # alternative: less than 0.5
        p1_r = binomial_cdf( c, (a+c), 0.5, lower=True )   #              greater than 0.5
        p2_l = binomial_cdf( b, (b+d), 0.5, lower=True )      # alternative: less than 0.5
        p2_r = binomial_cdf( d, (b+d), 0.5, lower=True )   #              greater than 0.5
        #print p1_l, p1_r, p2_l, p2_r

        if (p1_l < 0.05 and p2_r < 0.05) or (p1_r < 0.05 and p2_l < 0.05):
            # we reject loci where the significant biases are inconsistent between top1 and top2 alleles.
            return 1.0
        else:
            # if b<=2 and d=0 or b=0 and d<=2 -- highly possible FPs
            #if ( b<=2 and d==0 or b==0 and d<=2 ):
            #    return 1
            # can't decide
            return 0.0

    cdef float SB_score_ATAC( self, int a, int b, int c, int d ):
        """ calculate score for filtering variants with strange strand biases.

        ATAC-seq version

        a: top1/major allele plus strand
        b: top2/minor allele plus strand
        c: top1/major allele minus strand
        d: top2/minor allele minus strand

        Return a float value so that if this value >= 1, the variant will be filtered out.
        """
        cdef:
            float score
            double p
            double p1_l, p1_r
            double p2_l, p2_r
            double top2_sb, top1_sb

        if a+b == 0 or c+d == 0:
            # if major allele and minor allele both bias to the same strand, allow it
            return 0.0

        # Rule:
        # if there is bias in top2 allele then bias in top1 allele should not be significantly smaller than it.
        # or there is no significant bias (0.5) in top2 allele.
        
        #print a, b, c, d
        p1_l = binomial_cdf( a, (a+c), 0.5, lower=True )      # alternative: less than 0.5
        p1_r = binomial_cdf( c, (a+c), 0.5, lower=True )   #              greater than 0.5
        p2_l = binomial_cdf( b, (b+d), 0.5, lower=True )      # alternative: less than 0.5
        p2_r = binomial_cdf( d, (b+d), 0.5, lower=True )   #              greater than 0.5
        #print p1_l, p1_r, p2_l, p2_r

        if (p1_l < 0.05 and p2_r < 0.05) or (p1_r < 0.05 and p2_l < 0.05):
            # we reject loci where the significant biases are inconsistent between top1 and top2 alleles.
            return 1.0
        else:
            # can't decide
            return 0.0

    cpdef str to_vcf ( self ):
        """Output REF,ALT,QUAL,FILTER,INFO,FORMAT, SAMPLE columns.
        """
        cdef:
            str vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample

        vcf_ref = self.ref_allele.decode()
        vcf_alt = self.alt_allele.decode()
        vcf_qual = "%d" % self.GQ
        vcf_filter = "."
        vcf_info = (b"M=%s;MT=%s;DPT=%d;DPC=%d;DP1T=%d%s;DP2T=%d%s;DP1C=%d%s;DP2C=%d%s;SB=%d,%d,%d,%d;DBIC=%.2f;BICHOMOMAJOR=%.2f;BICHOMOMINOR=%.2f;BICHETERNOAS=%.2f;BICHETERAS=%.2f;AR=%.2f" % \
            (self.type.encode(), self.mutation_type.encode(), sum( self.n_reads_T.values() ), sum( self.n_reads_C.values() ), 
             self.n_reads_T[self.top1allele], self.top1allele, self.n_reads_T[self.top2allele], self.top2allele,
             self.n_reads_C[self.top1allele], self.top1allele, self.n_reads_C[self.top2allele], self.top2allele,
             self.n_strand[ 0 ][ self.top1allele ], self.n_strand[ 0 ][ self.top2allele ], self.n_strand[ 1 ][ self.top1allele ], self.n_strand[ 1 ][ self.top2allele ],
             self.deltaBIC,
             self.BIC_homo_major, self.BIC_homo_minor, self.BIC_heter_noAS,self.BIC_heter_AS,
             self.n_reads_T[self.top1allele]/(self.n_reads_T[self.top1allele]+self.n_reads_T[self.top2allele])
             )).decode()
        vcf_format = "GT:DP:GQ:PL"
        vcf_sample = "%s:%d:%d:%d,%d,%d" % (self.GT, self.raw_read_depth( opt = "all" ), self.GQ, self.PL_00, self.PL_01, self.PL_11)
        return "\t".join( ( vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample ) )

    cpdef toVariant ( self ):
        cdef:
            object v
        v = Variant( 
                     self.ref_allele.decode(),
                     self.alt_allele.decode(),
                     self.GQ,
                     '.',
                     self.type,
                     self.mutation_type,
                     self.top1allele.decode(),
                     self.top2allele.decode(),
                     sum( self.n_reads_T.values() ),
                     sum( self.n_reads_C.values() ),
                     self.n_reads_T[self.top1allele],
                     self.n_reads_T[self.top2allele],
                     self.n_reads_C[self.top1allele],
                     self.n_reads_C[self.top2allele],
                     self.n_strand[ 0 ][ self.top1allele ],
                     self.n_strand[ 0 ][ self.top2allele ],
                     self.n_strand[ 1 ][ self.top1allele ],
                     self.n_strand[ 1 ][ self.top2allele ],
                     self.deltaBIC,
                     self.BIC_homo_major,
                     self.BIC_homo_minor,
                     self.BIC_heter_noAS,
                     self.BIC_heter_AS,
                     self.n_reads_T[self.top1allele]/(self.n_reads_T[self.top1allele]+self.n_reads_T[self.top2allele]),
                     self.GT,
                     self.raw_read_depth( opt = "all" ),
                     self.PL_00,
                     self.PL_01,
                     self.PL_11 )
        return v
