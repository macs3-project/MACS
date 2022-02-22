# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-02-18 11:44:57 Tao Liu>

"""Module for SAPPER BAMParser class

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
from collections import Counter
from operator import itemgetter
from copy import copy

from MACS3.Signal.ReadAlignment import ReadAlignment
from MACS3.Signal.PosReadsInfo import PosReadsInfo
from MACS3.IO.PeakIO import PeakIO

from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t, int64_t

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

__DNACOMPLEMENT__ = b'\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQRSAUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff' # A trans table to convert A to T, C to G, G to C, and T to A.

__CIGARCODE__ = "MIDNSHP=X"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class UnitigRAs:
    """
    """
    cdef:
        list RAlists            # [RAlists_T, RAlists_C]
        bytes seq
        bytes unitig_aln
        bytes reference_aln
        bytes chrom
        long lpos
        long rpos
        long unitig_length
        long reference_length
        long aln_length

    def __init__ ( self, bytes chrom, long lpos, long rpos, bytes unitig_aln, bytes reference_aln, list RAlists ):
        assert len( unitig_aln )==len( reference_aln ), Exception("aln on unitig and reference should be the same length!")
        self.chrom = chrom
        self.lpos = lpos
        self.rpos = rpos
        self.unitig_aln = unitig_aln
        self.reference_aln = reference_aln
        self.RAlists = RAlists
        # fill in other information
        self.seq = self.unitig_aln.replace(b'-',b'')
        self.unitig_length = len( self.seq )
        self.reference_length = rpos - lpos
        self.aln_length = len( unitig_aln )

    def __getitem__ ( self, keyname ):
        if keyname == "chrom":
            return self.chrom
        elif keyname == "lpos":
            return self.lpos
        elif keyname == "rpos":
            return self.rpos
        elif keyname == "seq":
            return self.seq
        elif keyname == "unitig_aln":
            return self.unitig_aln
        elif keyname == "reference_aln":
            return self.reference_aln	    
        elif keyname == "unitig_length":
            return self.unitig_length
        elif keyname == "reference_length":
            return self.reference_length
        elif keyname == "aln_length":
            return self.aln_length
        elif keyname == "count":
            return len( self.RAlists[0] ) + len( self.RAlists[1] )
        else:
            raise KeyError("Unavailable key:", keyname)

    def __getstate__ ( self ):
        return (self.RAlists, self.seq, self.unitig_aln, self.reference_aln, self.chrom, self.lpos, self.rpos, self.unitig_length, self.reference_length, self.aln_length )
        
    def __setstate__ ( self, state ):
        (self.RAlists, self.seq, self.unitig_aln, self.reference_aln, self.chrom, self.lpos, self.rpos, self.unitig_length, self.reference_length, self.aln_length ) = state


    cpdef tuple get_variant_bq_by_ref_pos( self, long ref_pos ):
        """
        
        return ( s, bq_list_t, bq_list_c, strand_list_t, strand_list_c ) 
        """
        cdef:
            long i
            long index_aln
            long index_unitig
            long residue
            object ra
            bytes s
            list bq_list_t = []
            list bq_list_c = []
            list strand_list_t = []
            list strand_list_c = []
            list tip_list_t = []
            list pos_list_t = []
            list pos_list_c = []
            bytes ra_seq
            long ra_pos
            int p_seq
            int l_read

        #  b'TTATTAGAAAAAAT' find = 2
        #          b'AAAAATCCCACAGG'
        #b'TTTTATTAGAAAAAATCCCACAGGCAGCCACTAGGTGGCAGTAACAGGCTTTTGCCAGCGGCTCCAGTCAGCATGGCTTGACTGTGTGCT'
        #b'TTTTATTACAAAAA-TCCCACAGGCAGCCACTAGGTGGCAGTAACAGGCTTTTGCCAGCGGCTCCAGTCAGCATGGCTTGACTGTGTGCT' lpos=100
        #          |    |       |
        #genome   108  113     120
        #aln       8    13      21
        #unitig    8    13      21
        #ref       8    13      20
        #read1     6    11
        #read2           3      11
        # find the position
        residue = ref_pos - self.lpos + 1
        index_aln  = 0
        for i in range( self.aln_length ):
            if self.reference_aln[ i ] != 45: # 45 means b'-'
                residue -= 1
            if residue == 0:
                break
            index_aln += 1

        # index_aln should be the position on aln
        s = self.unitig_aln[ index_aln:index_aln+1 ]
        # find the index on unitig
        index_unitig = len( self.unitig_aln[:index_aln+1].replace(b'-',b'') )

        if s == b'-':                     #deletion
            for ra in self.RAlists[ 0 ]:
                ra_seq = ra["SEQ"]
                l_read = ra["l"]
                ra_pos = index_unitig - self.seq.find( ra_seq ) - 1
                if ra_pos == 0 or ra_pos == l_read -1:
                    tip_list_t.append( True )
                else:
                    tip_list_t.append( False )
                bq_list_t.append(  93 )
                strand_list_t.append( ra["strand"] )
                pos_list_t.append( ra_pos )
            for ra in self.RAlists[ 1 ]:
                ra_seq = ra["SEQ"]
                ra_pos = index_unitig - self.seq.find( ra_seq ) - 1
                bq_list_c.append(  93 )
                strand_list_c.append( ra["strand"] )
                pos_list_c.append( ra_pos )
            return ( bytearray(b'*'), bq_list_t, bq_list_c, strand_list_t, strand_list_c, tip_list_t, pos_list_t, pos_list_c )

        if index_aln < self.aln_length - 1:
            for i in range( index_aln + 1, self.aln_length ):
                if self.reference_aln[ i ] == 45:   #insertion detected, 45 means b'-'
                    s += self.unitig_aln[ i:i+1 ] # we extend the s string to contain the inserted seq
                else:
                    break

        for ra in self.RAlists[0]:        #treatment
            ra_seq = ra["SEQ"]
            l_read = ra["l"]
            ra_pos = index_unitig - self.seq.find( ra_seq ) - 1
            if ra_pos < l_read and ra_pos >= 0:
                pos_list_t.append( ra_pos )
                if ra_pos == 0 or ra_pos == l_read -1:
                    tip_list_t.append( True )
                else:
                    tip_list_t.append( False )
                bq_list_t.append( ra["binaryqual"][ra_pos] )
                strand_list_t.append( ra["strand"] )

        for ra in self.RAlists[1]:        #control
            ra_seq = ra["SEQ"]
            l_read = ra["l"]
            ra_pos = index_unitig - self.seq.find( ra_seq ) - 1
            if ra_pos < l_read and ra_pos >= 0:
                pos_list_c.append( ra_pos )
                bq_list_c.append( ra["binaryqual"][ra_pos] )                
                strand_list_c.append( ra["strand"] )

        return (bytearray(s), bq_list_t, bq_list_c, strand_list_t, strand_list_c, tip_list_t, pos_list_t, pos_list_c )

cdef class UnitigCollection:
    """A collection of ReadAlignment objects and the corresponding
    PeakIO.

    """
    cdef:
        bytes chrom
        object peak             # A PeakIO object
        list URAs_list
        long left               # left position of peak
        long right              # right position of peak
        long length             # length of peak
        long URAs_left          # left position of all RAs in the collection
        long URAs_right         # right position of all RAs in the collection
        bool sorted             # if sorted by lpos

    def __init__ ( self, chrom, peak, URAs_list=[] ):
        self.chrom = chrom
        self.peak = peak
        self.URAs_list = URAs_list
        self.left = peak["start"]
        self.right = peak["end"]
        self.length =  self.right - self.left
        self.URAs_left = URAs_list[ 0 ]["lpos"] # initial assignment of RAs_left
        self.URAs_right = URAs_list[-1]["rpos"] # initial assignment of RAs_right
        self.sort()                           # it will set self.sorted = True
        # check RAs_left and RAs_right
        for ura in URAs_list:
            if ura[ "lpos" ] < self.URAs_left:
                self.URAs_left = ura[ "lpos" ]
            if ura[ "rpos" ] > self.URAs_right:
                self.URAs_right = ura[ "rpos" ]

    def __getitem__ ( self, keyname ):
        if keyname == "chrom":
            return self.chrom
        elif keyname == "left":
            return self.left
        elif keyname == "right":
            return self.right
        elif keyname == "URAs_left":
            return self.URAs_left
        elif keyname == "URAs_right":
            return self.URAs_right
        elif keyname == "length":
            return self.length
        elif keyname == "count":
            return len( self.URAs_list )
        elif keyname == "URAs_list":
            return self.URAs_list
        else:
            raise KeyError("Unavailable key:", keyname)

    def __getstate__ ( self ):
        return (self.chrom, self.peak, self.URAs_list, self.left, self.right, self.length, self.URAs_left, self.URAs_right, self.sorted)
        
    def __setstate__ ( self, state ):
        (self.chrom, self.peak, self.URAs_list, self.left, self.right, self.length, self.URAs_left, self.URAs_right, self.sorted) = state
        
    cpdef sort ( self ):
        """Sort RAs according to lpos. Should be used after realignment.

        """
        self.URAs_list.sort(key=itemgetter("lpos"))
        self.sorted = True
        return
        
    cpdef object get_PosReadsInfo_ref_pos ( self, long ref_pos, bytes ref_nt, int Q=20 ):
        """Generate a PosReadsInfo object for a given reference genome
        position.

        Return a PosReadsInfo object.

        """
        cdef:
            bytearray s, bq
            list bq_list_t, bq_list_c, strand_list_t, strand_list_c, tip_list_t, pos_list_t, pos_list_c
            object ura
            int i

        posreadsinfo_p = PosReadsInfo( ref_pos, ref_nt )
        for i in range( len( self.URAs_list ) ):
            ura = self.URAs_list[ i ]
            if ura[ "lpos" ] <= ref_pos and ura[ "rpos" ] > ref_pos:
                ( s, bq_list_t, bq_list_c, strand_list_t, strand_list_c, tip_list_t, pos_list_t, pos_list_c ) = ura.get_variant_bq_by_ref_pos( ref_pos )
                for i in range( len(bq_list_t) ):
                    posreadsinfo_p.add_T( i, bytes(s), bq_list_t[i], strand_list_t[i], tip_list_t[i], Q=Q )
                for i in range( len(bq_list_c) ):
                    posreadsinfo_p.add_C( i, bytes(s), bq_list_c[i], strand_list_c[i], Q=Q )

        return posreadsinfo_p
