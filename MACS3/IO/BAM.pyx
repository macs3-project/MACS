# cython: language_level=3
# cython: profile=True
# Time-stamp: <2021-01-22 17:12:44 Tao Liu>

"""Module for SAPPER BAMParser class. It's derived from BAMParser in
Parser.pyx. The difference is that, in this new BAMParser class, we
will keep the alignment information and read sequence. Also, we do not
differentiate single-end vs paired-end data, and intend to read BAM
always as single-end.

Copyright (c) 2021 Tao Liu <vladimir.liu@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

"""

# ------------------------------------
# python modules
# ------------------------------------
import logging
from logging import info, debug
import struct
from struct import unpack
import gzip
import io
import sys

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Signal.ReadAlignment import ReadAlignment

# ------------------------------------
# Other modules
# ------------------------------------
from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t, int64_t

is_le = sys.byteorder == "little"

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
__version__ = "BAM $Revision$"
__author__ = "Tao Liu <vladimir.liu@gmail.com>"
__doc__ = "SAPPER BAMParser class"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
class StrandFormatError( Exception ):
    """Exception about strand format error.

    Example:
    raise StrandFormatError('Must be F or R','X')
    """
    def __init__ ( self, string, strand ):
        self.strand = strand
        self.string = string
        
    def __str__ ( self ):
        return repr( "Strand information can not be recognized in this line: \"%s\",\"%s\"" % ( self.string, self.strand ) )

class MDTagMissingError( Exception ):
    """Exception about missing MD tag

    Example:
    raise MDTagMissingError( name, aux )

    aux is the auxiliary data part.
    """
    def __init__ ( self, name, aux ):
        self.name = name
        self.aux = aux
        
    def __str__ ( self ):
        return repr( "MD tag is missing! Please use \"samtools calmd\" command to add MD tags!\nName of sequence:%s\nCurrent auxiliary data section: %s" % (self.name.decode(), self.aux.decode() ) )

cdef class BAMParser:
    """File Parser Class for BAM File.

    File is gzip-compatible and binary.
    Information available is the same that is in SAM format.
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    2048	supplementary alignment
    """
    cdef str filename
    cdef bool gzipped
    cdef object fhd             # BAM file handler
    cdef list references        # name of references/chromosomes, contain the order of chromosomes
    cdef dict rlengths          # lengths of references/chromosomes
    cdef bool sorted            # if BAM is sorted
    cdef long SOA               # position in file for the Start of Alignment blocks

    def __init__ ( self, str filename ):
        """Open input file. Determine whether it's a gzipped file.

        'filename' must be a string object.

        This function initialize the following attributes:

        1. self.filename: the filename for input file.
        2. self.gzipped: a boolean indicating whether input file is gzipped.
        3. self.fhd: buffered I/O stream of input file
        """
        self.filename = filename
        self.gzipped = True

        # try gzip first
        f = gzip.open( filename )
        try:
            f.read( 10 )
        except IOError:
            # not a gzipped file
            self.gzipped = False
        f.close()
        if self.gzipped:
            # open with gzip.open, then wrap it with BufferedReader!
            self.fhd = io.BufferedReader( gzip.open( filename, mode='rb' ), buffer_size = READ_BUFFER_SIZE )
        else:
            self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!
        #self.get_tsize()
        self.parse_header()
        # if not sorted, raise expection and terminate the program!
        if not self.sorted:
            raise Exception("BAM should be sorted by coordinates!")

    cpdef close ( self ):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.fhd.close()


    cpdef reset ( self ):
        self.fhd.seek( self.SOA )

    cpdef sniff( self ):
        """Check the first 3 bytes of BAM file. If it's 'BAM', check
        is success.

        """
        magic_header = self.fhd.read( 3 )
        if magic_header == b"BAM":
            tsize  = self.tsize()
            if tsize > 0:
                self.fhd.seek( 0 )
                return True
            else:
                self.fhd.seek( 0 )
                raise Exception( "File is not of a valid BAM format! %d" % tsize )
        else:
            self.fhd.seek( 0 )
            return False

    cdef parse_header( self ):
        """Parse the HEADER of BAM. It will do the following thing:

        1. read in references from BAM header, and save in
        self.references, and self.rlengths.  

        2. invoke self.check_sorted to see if BAM is sorted by
        coordinates, and set self.sorted.

        3. remember the file pointer at the end of HEADER or the start
        of the alignment blocks, and set self.SOA.
        
        """
        cdef:
            int header_len, x, nc, nlength
            bytes refname
            list references = []
            dict rlengths = {}
            bytes header_text
            
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something, the first 4 bytes is called magic string
        fseek(4)
        # next 4bytes is length of header 
        header_len =  unpack( '<i', fread( 4 ) )[ 0 ]
        header_text = unpack( '<%ds' % header_len, fread( header_len ) )[0]

        # get the number of chromosome
        nc = unpack( '<i', fread( 4 ) )[ 0 ]
        for x in range( nc ):
            # read each chromosome name
            nlength = unpack( '<i', fread( 4 ) )[ 0 ]
            refname = fread( nlength )[ :-1 ]
            references.append( refname )
            # don't jump over chromosome size
            # we can use it to avoid falling of chrom ends during peak calling
            rlengths[refname] = unpack( '<i', fread( 4 ) )[ 0 ]
        self.references = references
        self.rlengths = rlengths

        # read sorting information
        self.check_sorted( header_text )
        
        # remember the Start of Alignment blocks
        self.SOA = self.fhd.tell()
        return

    cdef check_sorted( self, bytes header ):
        """Check if the BAM is sorted by looking at the HEADER text
        """
        cdef:
            bytes l
            int i

        self.sorted = False
        # HD line is the first line
        l = header.split(b"\n")[0]
        if l.startswith(b"@HD"):
            # locate HD line
            i = l.find(b"SO:")
            if i != -1:
                # I will only check the first 5 characters
                if l[i+3:i+8] == b"coord":
                    self.sorted = True
        return
    
    cpdef get_chromosomes ( self ):
        """Get chromosomes in order of their appearance in BAM HEADER.

        """
        return self.references

    cpdef get_rlengths ( self ):
        """Get chromosomes in order of their appearance in BAM HEADER.

        """
        return self.rlengths

    cpdef get_reads_in_region ( self, bytes chrom, int left, int right, int maxDuplicate = 1 ):
        """Get reads in a given region. Initial call will start at
        'self.SOA', but will keep incrementing.

        Return: 1. list of ReadAlignment; 2. current position in file
        """
        cdef:
            int i = 0
            int m = 0
            int entrylength, fpos, strand, chrid
            list readslist
            int cur_duplicates = 0
            object read
            object previous_read
        
        readslist = []
        fread = self.fhd.read
        fseek = self.fhd.seek
        ftell = self.fhd.tell
        cur_duplicates = 0
        previous_read = None
        while True:
            try:
                entrylength = unpack( '<i', fread( 4 ) )[ 0 ]
            except struct.error:
                # if reaching the EOF, this will trigger
                break
            read = self.__fw_binary_parse( fread( entrylength ) )
            #print( "Got a record", left, right, "\n" )            
            if read != None:
                if read["chrom"] == chrom and read["lpos"] < right and read["rpos"] > left:
                    #print( "Found an overlapped read", read["chrom"],read["lpos"],read["rpos"], "\n" )
                    # an overlap is found
                    if previous_read != None and previous_read["lpos"] == read["lpos"] and previous_read["rpos"] == read["rpos"] \
                            and previous_read["strand"] ==  read["strand"] and previous_read["cigar"] == read["cigar"]:
                        cur_duplicates += 1
                    else:
                        cur_duplicates = 1
                    if cur_duplicates <= maxDuplicate:
                        readslist.append( read )
                    previous_read = read
                elif read["chrom"] != chrom or read["lpos"] > right:
                    # pass the region, rewind, then trigger 'break'
                    fseek( ftell()-entrylength-4 )
                    #print( "rewind:", ftell() )
                    break
                #else:
                #    print( "NOT an overlapped read", read["chrom"],read["lpos"],read["rpos"], "\n" )

        #print( "end", ftell() )
        #print( "Get # of reads:", len(readslist),"\n")
        return readslist

    cdef object __fw_binary_parse (self, data, min_MAPQ=1 ):
        """ Read information from an alignment block. Discard 
        alignment below a minimum MAPQ score, default: 1.

        Return: ReadAlignment object with only selected information,
        including refname (chromosome), leftmost alignment location,
        strand information, sequence, quality, CIGAR and MD tag.

        Note: We do require MD tag exists in BAM file.
        """
        cdef:
            int ref, leftmost, rightmost, i, j, l_seq, strand
            short bwflag, l_read_name, n_cigar_op, MAPQ
            short bin_bam
            bytes read_name
            bytes seq #note: for each byte, 1st base in the highest 4bit; 2nd in the lowest 4bit. "=ACMGRSVTWYHKDBN" -> [0,15]
            bytes qual
            tuple cigar_op  # op_len<<4|op, op: "MIDNSHP=X" -> 012345678
            bytes tag
            bytes MD
        
        if not data: return None

        # read number of CIGAR operators, Bitwise FLAG
        (n_cigar_op,  bwflag ) = unpack( '<HH' , data[ 12:16 ] )

        # first, we will discard problematic reads
        if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
            return None       #unmapped sequence or bad sequence or  secondary or supplementary alignment 
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # Different with MACS, both reads will be kept.
            if not bwflag & 2:
                return None   # not a proper pair
            if bwflag & 8:
                return None   # the mate is unmapped
        
        # read length of readname, MAPQ, and bin information
        (l_read_name, MAPQ, bin_bam) = unpack( '<BBH', data[ 8:12 ] )

        # we will also discard reads having MAPQ lower than min_MAPQ. MAPQ=255 means no alignment
        if MAPQ < min_MAPQ or MAPQ == 255:
            return None

        # Now we read other information

        # read ref name
        ref = unpack( '<i', data[ 0:4 ] )[ 0 ]
        # read leftmost position of alignment
        leftmost = unpack( '<i', data[ 4:8 ] )[ 0 ]
        # read length of query sequence
        l_seq = unpack( '<i', data[ 16:20 ]  )[0]
        # readname , skip next_refID, next_pos, tlen, which is 12 bytes or from the 32th byte
        read_name = unpack( '<%ds' % (l_read_name), data[ 32: 32+l_read_name ] )[0][:-1] # last byte is \x00
        # cigar_op, find the index i first.
        i = 32 + l_read_name
        cigar_op = unpack( '<%dI' % (n_cigar_op) , data[ i : i + n_cigar_op*4 ] )
        # read sequence information
        i += n_cigar_op*4
        seq = unpack( '<%ds' % int((l_seq+1)/2), data[ i: i + int((l_seq+1)/2)] )[0]
        # read quality information. Note: the value in BAM is the acutal Phred score, there is no +33!
        i += int((l_seq+1)/2)
        qual = unpack( '<%ds' % (l_seq), data[ i : i + l_seq])[0]

        rightmost = leftmost
        for j in cigar_op:
            if j & 15 in [ 0, 2, 3, 7, 8 ]:   # they are CIGAR op M/D/N/=/X, no need to add I, S, H, or P
                rightmost += j >> 4
                    
        # strand information
        if bwflag & 16:
            # reverse strand
            strand = 1
        else:
            strand = 0

        # MD tag
        MD = b''
        i += l_seq
        tag = data[ i : ]
        j = tag.find(b'MDZ')
        if j == -1: raise MDTagMissingError( data[ 32: 32+l_read_name], tag )
        MD = tag[ j+3 : tag[j:].find(b"\0") + j ]

        # construct a ReadAlignment object and return
        return ReadAlignment( read_name, self.references[ref], leftmost, rightmost, strand, seq, qual, cigar_op, MD )

