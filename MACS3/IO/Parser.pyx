# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2022-02-02 13:24:26 Tao Liu>

"""Module for all MACS Parser classes for input. Please note that the
parsers are for reading the alignment files ONLY.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).

"""

# ------------------------------------
# python modules
# ------------------------------------
import struct
from struct import unpack
from re import findall
import gzip
import io
import sys

import logging
import MACS3.Utilities.Logger

logger = logging.getLogger(__name__)
debug   = logger.debug
info    = logger.info
# ------------------------------------
# MACS3 modules
# ------------------------------------

from MACS3.Utilities.Constants import *
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI

# ------------------------------------
# Other modules
# ------------------------------------
from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t

is_le = sys.byteorder == "little"

cdef extern from "stdlib.h":
    ctypedef uint32_t size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    void free(void *ptr)
    int32_t strcmp(char *a, char *b)
    char * strcpy(char *a, char *b)
    int64_t atol(char *str)
    int32_t atoi(char *str)

# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

cpdef guess_parser ( fname, int64_t buffer_size = 100000 ):
    # Note: BAMPE and BEDPE can't be automatically detected.
    ordered_parser_dict = {"BAM":BAMParser,
                           "BED":BEDParser,
                           "ELAND":ELANDResultParser,
                           "ELANDMULTI":ELANDMultiParser,
                           "ELANDEXPORT":ELANDExportParser,
                           "SAM":SAMParser,
                           "BOWTIE":BowtieParser}

    for f in ordered_parser_dict:
        p = ordered_parser_dict[ f ]
        t_parser = p( fname, buffer_size = buffer_size )
        debug( "Testing format %s" % f )
        s = t_parser.sniff()
        if s:
            info( "Detected format is: %s" % ( f ) )
            if t_parser.is_gzipped():
                info( "* Input file is gzipped." )
            return t_parser
        else:
            t_parser.close()
    raise Exception( "Can't detect format!" )

cdef tuple __bam_fw_binary_parse_le ( const unsigned char * data ):
    """Parse a BAM SE entry in little endian system
    """
    cdef:
        int32_t thisref, thisstart, thisstrand
        uint16_t bwflag
        uint8_t l_read_name
        uint16_t n_cigar_op
        int32_t cigar_code
        uint8_t *ui8
        int32_t *i32
        uint16_t *ui16
        uint32_t *ui32

    # we skip lot of the available information in data (i.e. tag name, quality etc etc)
    # no data, return, does it really happen without raising struct.error?
    if not data: return ( -1, -1, -1 )

    ui16 = <uint16_t *>data
    bwflag = ui16[7]

    # we filter out unmapped sequence or bad sequence or  secondary or supplementary alignment
    # we filter out 2nd mate, not a proper pair, mate is unmapped
    if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)): return ( -1, -1, -1 )

    i32 = <int32_t *>data
    thisref = i32[0]
    thisstart = i32[1]
    n_cigar_op = ui16[6]
    # In case of paired-end we have now skipped all possible "bad" pairs
    # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
    # we can treat it as a single read, so just check the strand and calculate its
    # start position... hope I'm right!
    if bwflag & 16:
        # read mapped to minus strand; then we have to compute cigar to find the rightmost position
        ui8 = <uint8_t *>data
        l_read_name = ui8[8]
        # need to decipher CIGAR string
        ui32 = <uint32_t *>(data + 32 + l_read_name) # move pointer to cigar_code
        for cigar_code in ui32[:n_cigar_op]:#unpack( '<%dI' % (n_cigar_op) , data[ 32 + l_read_name : 32 + l_read_name + n_cigar_op*4 ] ):
            if cigar_code & 15 in [ 0, 2, 3, 7, 8 ]:   # they are CIGAR op M/D/N/=/X
                thisstart += cigar_code >> 4
        thisstrand = 1
    else:
        thisstrand = 0

    return ( thisref, thisstart, thisstrand )

cdef tuple __bam_fw_binary_parse_be ( const unsigned char * data ):
    """Big endian version. We need byte swap.
    """
    cdef:
        int32_t thisref, thisstart, thisstrand
        uint16_t bwflag
        uint8_t l_read_name
        uint16_t n_cigar_op
        int32_t cigar_code
        uint8_t *ui8                      # we will only cast 1 byte at a time
        int32_t i
        uint32_t shift0, shift


    # we skip lot of the available information in data (i.e. tag name, quality etc etc)
    # no data, return, does it really happen without raising struct.error?
    if not data: return ( -1, -1, -1 )

    ui8 = <uint8_t *>data
    bwflag = ui8[15] << 8 | ui8[14]       # it works as bwflag = ui16[7] in little-endian

    # we filter out unmapped sequence or bad sequence or  secondary or supplementary alignment
    # we filter out 2nd mate, not a proper pair, mate is unmapped
    if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)): return ( -1, -1, -1 )

    # the following three lins are for little-endian
    #thisref = i32[0]
    #thisstart = i32[1]
    #n_cigar_op = ui16[6]

    # to simplify the byte swap, we pretend all original numbers (thisref, pos, nextpos) positive
    thisref = ui8[3] << 24 | ui8[2] << 16 | ui8[1] << 8 | ui8[0]
    thisstart = ui8[7] << 24 | ui8[6] << 16 | ui8[5] << 8 | ui8[4]
    n_cigar_op = ui8[13] << 8 | i8[12]

    # In case of paired-end we have now skipped all possible "bad" pairs
    # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
    # we can treat it as a single read, so just check the strand and calculate its
    # start position... hope I'm right!
    if bwflag & 16:
        # read mapped to minus strand; then we have to compute cigar to find the rightmost position
        l_read_name = ui8[8]
        # need to decipher CIGAR string
        # move pointer to cigar_code
        shift0 = 32 + l_read_name
        for i in range(n_cigar_op):
            shift = shift0 + i*4          # move 32bit at a time
            cigar_code = ui8[shift0+3] << 24 | ui8[shift0+2] << 16 | ui8[shift0+1] << 8 | ui8[shift0] # it works like cigar_code = ui32[...] in little-endian
            if cigar_code & 15 in [ 0, 2, 3, 7, 8 ]:   # they are CIGAR op M/D/N/=/X
                thisstart += cigar_code >> 4
        thisstrand = 1
    else:
        thisstrand = 0

    return ( thisref, thisstart, thisstrand )

cdef tuple __bampe_pe_binary_parse_le (const unsigned char * data):
    """Parse a BAMPE record in little-endian system.
    """
    cdef:
        int32_t thisref, thisstart, thistlen
        int32_t nextpos, pos
        uint16_t bwflag
        uint8_t *ui8
        int32_t *i32
        uint16_t *ui16

    # we skip lot of the available information in data (i.e. tag name, quality etc etc)
    if not data:
        return ( -1, -1, -1 )

    ui16 = <uint16_t *>data
    bwflag = ui16[7]
    # we filter out unmapped, bad sequence, secondary/supplementary alignment
    # we filter out other mate of paired reads, not a proper pair, or mate is unmapped
    if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)):
        return ( -1, -1, -1 )

    i32 = <int32_t *>data
    ui8 = <uint8_t *>data
    thisref = i32[0]

    pos = i32[1]
    nextpos = i32[6]
    thistlen = i32[7]
    thisstart = min(pos, nextpos) # we keep only the leftmost
                                  # position which means this must
                                  # be at + strand. So we don't
                                  # need to decipher CIGAR string.
    thistlen = abs( thistlen )                    # Actually, if
    #                                             # the value
    #                                             # unpacked is
    #                                             # negative, then
    #                                             # nextpos is the
    #                                             # leftmost
    #                                             # position.

    return ( thisref, thisstart, thistlen )

cdef tuple __bampe_pe_binary_parse_be (const unsigned char * data):
    """Parse a BAMPE record in big-endian system. And we need byte swap.
    """
    cdef:
        int32_t thisref, thisstart, thistlen
        uint32_t tmp_thistlen
        int32_t nextpos, pos
        uint16_t bwflag
        uint8_t *ui8                      # we will only cast 1 byte at a time

    # we skip lot of the available information in data (i.e. tag name, quality etc etc)
    if not data:
        return ( -1, -1, -1 )

    ui8 = <uint8_t *>data

    bwflag = ui8[15] << 8 | ui8[14]       # as le: bwflag = ui16[7]
    # we filter out unmapped, bad sequence, secondary/supplementary alignment
    # we filter out other mate of paired reads, not a proper pair, or mate is unmapped
    if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)):
        return ( -1, -1, -1 )

    i8 = <int8_t *>data

    # the following three lins are for little-endian



    # to simplify the byte swap, we pretend all original numbers (thisref, pos, nextpos) positive
    thisref = ui8[3] << 24 | ui8[2] << 16 | ui8[1] << 8 | ui8[0]   # as le:thisref = i32[0]
    pos = ui8[7] << 24 | ui8[6] << 16 | ui8[5] << 8 | ui8[4]       # as le:pos = i32[1]
    nextpos = ui8[27] << 24 | ui8[26] << 16 | ui8[25] << 8 | ui8[24] # as le:nextpos = i32[6]

    # thistlen can be negative, so we byte swap it then convert to int32_t then take abs (maybe there is more effecient way?)
    tmp_thistlen = ui8[31] << 24 | ui8[30] << 16 | ui8[29] << 8 | ui8[28] # as le:tmp_thistlen = ui32[7]
    thistlen = abs(<int32_t> tmp_thistlen)

    # position which means this must
    # be at + strand. So we don't
    # need to decipher CIGAR string.
    thisstart = pos if nextpos > pos else nextpos #min(pos, nextpos) # we keep only the leftmost

    return ( thisref, thisstart, thistlen )


# choose a parser according to endian
if is_le:
    bam_se_entry_parser = __bam_fw_binary_parse_le
    bampe_pe_entry_parser = __bampe_pe_binary_parse_le
else:
    bam_se_entry_parser = __bam_fw_binary_parse_be
    bampe_pe_entry_parser = __bampe_pe_binary_parse_be

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

cdef class GenericParser:
    """Generic Parser class.

    Inherit this class to write your own parser. In most cases, you need to override:

    1. __tlen_parse_line which returns tag length of a line
    2.  __fw_parse_line which returns tuple of ( chromosome, 5'position, strand )
    """
    cdef:
        str filename
        bool gzipped
        int32_t tag_size
        object fhd
        int64_t buffer_size

    def __init__ ( self, str filename, int64_t buffer_size = 100000 ):
        """Open input file. Determine whether it's a gzipped file.

        'filename' must be a string object.

        This function initialize the following attributes:

        1. self.filename: the filename for input file.
        2. self.gzipped: a boolean indicating whether input file is gzipped.
        3. self.fhd: buffered I/O stream of input file
        """
        self.filename = filename
        self.gzipped = True
        self.tag_size = -1
        self.buffer_size = buffer_size
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
            self.fhd = io.BufferedReader( gzip.open( filename, mode='rb' ), buffer_size = READ_BUFFER_SIZE ) # buffersize set to 10M
        else:
            self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!
        self.__skip_first_commentlines()

    cdef void __skip_first_commentlines ( self ):
        """Some parser needs to skip the first several comment lines.

        Redefine this if it's necessary!
        """
        return

    cpdef int32_t tsize( self ):
        """General function to detect tag size.

        * Although it can be used by most parsers, it must be
          rewritten by BAMParser!
        """
        cdef:
            int32_t s = 0
            int32_t n = 0     # number of successful/valid read alignments
            int32_t m = 0     # number of trials
            int32_t this_taglength
            bytes thisline

        if self.tag_size != -1:
            # if we have already calculated tag size (!= -1),  return it.
            return self.tag_size

        # try 10k times or retrieve 10 successfule alignments
        while n < 10 and m < 10000:
            m += 1
            thisline = self.fhd.readline()
            this_taglength = self.__tlen_parse_line( thisline )
            if this_taglength > 0:
                # this_taglength == 0 means this line doesn't contain
                # successful alignment.
                s += this_taglength
                n += 1
        # done
        self.fhd.seek( 0 )
        self.__skip_first_commentlines()
        if n != 0:              # else tsize = -1
            self.tag_size = <int32_t>(s/n)
        return self.tag_size

    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Abstract function to detect tag length.

        """
        raise NotImplemented

    cpdef build_fwtrack ( self ):
        """Generic function to build FWTrack object. Create a new
        FWTrack object. If you want to append new records to an
        existing FWTrack object, try append_fwtrack function.

        * BAMParser for binary BAM format should have a different one.
        """
        cdef:
            int64_t i, fpos, strand
            bytes chromosome
            bytes tmp = b""

        fwtrack = FWTrack( buffer_size = self.buffer_size )
        i = 0
        while True:
            # for each block of input
            tmp += self.fhd.read( READ_BUFFER_SIZE )
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[ -1 ]
            for thisline in lines[ :-1 ]:
                ( chromosome, fpos, strand ) = self.__fw_parse_line( thisline )
                if fpos < 0 or not chromosome:
                    # normally __fw_parse_line will return -1 if the line
                    # contains no successful alignment.
                    continue
                i += 1
                if i % 1000000 == 0:
                    info( " %d reads parsed" % i )
                fwtrack.add_loc( chromosome, fpos, strand )
        # last one
        if tmp:
            ( chromosome, fpos, strand ) = self.__fw_parse_line( tmp )
            if fpos >= 0 and chromosome:
                i += 1
                fwtrack.add_loc( chromosome, fpos, strand )
        # close file stream.
        self.close()
        return fwtrack

    cpdef append_fwtrack ( self, fwtrack ):
        """Add more records to an existing FWTrack object.

        """
        cdef:
            int64_t i, fpos, strand
            bytes chromosome
            bytes tmp = b""
        i = 0
        while True:
            # for each block of input
            tmp += self.fhd.read( READ_BUFFER_SIZE )
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[ -1 ]
            for thisline in lines[ :-1 ]:
                ( chromosome, fpos, strand ) = self.__fw_parse_line( thisline )
                if fpos < 0 or not chromosome:
                    # normally __fw_parse_line will return -1 if the line
                    # contains no successful alignment.
                    continue
                i += 1
                if i % 1000000 == 0:
                    info( " %d reads parsed" % i )
                fwtrack.add_loc( chromosome, fpos, strand )

        # last one
        if tmp:
            ( chromosome, fpos, strand ) = self.__fw_parse_line( tmp )
            if fpos >= 0 and chromosome:
                i += 1
                fwtrack.add_loc( chromosome, fpos, strand )
        # close file stream.
        self.close()
        return fwtrack

    cdef tuple __fw_parse_line ( self, bytes thisline ):
        """Abstract function to parse chromosome, 5' end position and
        strand.

        """
        cdef bytes chromosome = b""
        cdef int32_t fpos = -1
        cdef int32_t strand = -1
        return ( chromosome, fpos, strand )

    cpdef sniff ( self ):
        """Detect whether this parser is the correct parser for input
        file.

        Rule: try to find the tag size using this parser, if error
        occurs or tag size is too small or too big, check is failed.

        * BAMParser has a different sniff function.
        """
        cdef int32_t t

        t = self.tsize()
        if t <= 10 or t >= 10000: # tsize too small or too big
            self.fhd.seek( 0 )
            return False
        else:
            self.fhd.seek( 0 )
            self.__skip_first_commentlines()
            return True

    cpdef close ( self ):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.fhd.close()

    cpdef bool is_gzipped ( self ):
        return self.gzipped

cdef class BEDParser( GenericParser ):
    """File Parser Class for BED File.

    """
    cdef void __skip_first_commentlines ( self ):
        """BEDParser needs to skip the first several comment lines.
        """
        cdef:
            int32_t l_line
            bytes this_line
        for thisline in self.fhd:
            l_line = len( thisline )
            if thisline and ( thisline[ :5 ] != b"track" ) \
               and ( thisline[ :7 ] != b"browser" ) \
               and ( thisline[ 0 ] != 35 ): # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek( -l_line, 1 )
        return

    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Parse 5' and 3' position, then calculate frag length.

        """
        thisline = thisline.rstrip()
        if not thisline:
            return 0

        thisfields = thisline.split( b'\t' )
        return atoi( thisfields[ 2 ] )-atoi( thisfields[ 1 ] )

    cdef tuple __fw_parse_line ( self, bytes thisline ):
        #cdef list thisfields
        cdef:
            bytes chromname
            list thisfields

        thisline = thisline.rstrip()
        thisfields = thisline.split( b'\t' )
        chromname = thisfields[ 0 ]
        try:
            #if not strcmp(thisfields[ 5 ],b"+"):
            if thisfields[5] == b"+":
                return ( chromname,
                         atoi( thisfields[ 1 ] ),
                         0 )
            #elif not strcmp(thisfields[ 5 ], b"-"):
            elif thisfields[5] == b"-":
                return ( chromname,
                         atoi( thisfields[ 2 ] ),
                         1 )
            else:
                raise StrandFormatError( thisline, thisfields[ 5 ] )
        except IndexError:
            # default pos strand if no strand
            # info can be found
            return ( chromname,
                     atoi( thisfields[ 1 ] ),
                     0 )

cdef class BEDPEParser(GenericParser):
    """Parser for BED format file containing PE information, and also
    can be used for the cases when users predefine the fragment
    locations by shifting/extending by themselves.

    Format:

    chromosome_name	frag_leftend	frag_rightend


    Note: Only the first columns are used!
    """

    cdef public int32_t n
    cdef public float32_t d

    cdef void __skip_first_commentlines ( self ):
        """BEDPEParser needs to skip the first several comment lines.
        """
        cdef:
            int32_t l_line
            bytes this_line
        for thisline in self.fhd:
            l_line = len( thisline )
            if thisline and ( thisline[ :5 ] != b"track" ) \
               and ( thisline[ :7 ] != b"browser" ) \
               and ( thisline[ 0 ] != 35 ): # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek( -l_line, 1 )
        return

    cdef __pe_parse_line ( self, bytes thisline ):
        """ Parse each line, and return chromosome, left and right positions

        """
        cdef:
            list thisfields

        thisline = thisline.rstrip()

        # still only support tabular as delimiter.
        thisfields = thisline.split( b'\t' )
        try:
            return ( thisfields[ 0 ],
                     atoi( thisfields[ 1 ] ),
                     atoi( thisfields[ 2 ] ) )
        except IndexError:
            raise Exception("Less than 3 columns found at this line: %s\n" % thisline)

    cpdef build_petrack ( self ):
        """Build PETrackI from all lines.

        """
        cdef:
            bytes chromname
            int32_t left_pos
            int32_t right_pos
            int64_t i = 0          # number of fragments
            int64_t m = 0          # sum of fragment lengths
            bytes tmp = b""

        petrack = PETrackI( buffer_size = self.buffer_size )
        add_loc = petrack.add_loc

        while True:
            # for each block of input
            tmp += self.fhd.read( READ_BUFFER_SIZE )
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[ -1 ]
            for thisline in lines[ :-1 ]:
                ( chromosome, left_pos, right_pos ) = self.__pe_parse_line( thisline )
                if left_pos < 0 or not chromosome:
                    continue
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                m += right_pos - left_pos
                i += 1
                if i % 1000000 == 0:
                    info( " %d fragments parsed" % i )
                add_loc( chromosome, left_pos, right_pos )
        # last one
        if tmp:
            ( chromosome, left_pos, right_pos ) = self.__pe_parse_line( thisline )
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                add_loc( chromosome, left_pos, right_pos )
                
        self.d = <float32_t>( m ) / i
        self.n = i
        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"

        self.close()
        petrack.set_rlengths( {"DUMMYCHROM":0} )
        return petrack

    cpdef append_petrack (self, petrack):
        """Build PETrackI from all lines, return a PETrackI object.
        """
        cdef:
            bytes chromname
            int32_t left_pos
            int32_t right_pos
            int64_t i = 0          # number of fragments
            int64_t m = 0          # sum of fragment lengths
            bytes tmp = b""

        add_loc = petrack.add_loc
        while True:
            # for each block of input
            tmp += self.fhd.read( READ_BUFFER_SIZE )
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[ -1 ]
            for thisline in lines[ :-1 ]:
                ( chromosome, left_pos, right_pos ) = self.__pe_parse_line( thisline )
                if left_pos < 0 or not chromosome:
                    continue
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                m += right_pos - left_pos
                i += 1
                if i % 1000000 == 0:
                    info( " %d fragments parsed" % i )
                add_loc( chromosome, left_pos, right_pos )
        # last one
        if tmp:
            ( chromosome, left_pos, right_pos ) = self.__pe_parse_line( thisline )
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                add_loc( chromosome, left_pos, right_pos )

        self.d = ( self.d * self.n + m ) / ( self.n + i )
        self.n += i

        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"
        self.close()
        petrack.set_rlengths( {"DUMMYCHROM":0} )
        return petrack

cdef class ELANDResultParser( GenericParser ):
    """File Parser Class for tabular File.

    """
    cdef void __skip_first_commentlines ( self ):
        """ELANDResultParser needs to skip the first several comment lines.
        """
        cdef:
            int32_t l_line
            bytes this_line
        for thisline in self.fhd:
            l_line = len( thisline )
            if thisline and thisline[ 0 ] != 35: # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek( -l_line, 1 )
        return

    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( b'\t' )
        if thisfields[1].isdigit():
            return 0
        else:
            return len( thisfields[ 1 ] )

    cdef tuple __fw_parse_line ( self, bytes thisline ):
        cdef:
            bytes chromname, strand
            int32_t thistaglength
        thisline = thisline.rstrip()
        if not thisline: return ( b"", -1, -1 )

        thisfields = thisline.split( b'\t' )
        thistaglength = strlen( thisfields[ 1 ] )

        if len( thisfields ) <= 6:
            return ( b"", -1, -1 )

        try:
            chromname = thisfields[ 6 ]
            chromname = chromname[ :chromname.rindex( b".fa" ) ]
        except ValueError:
            pass

        if thisfields[ 2 ] == b"U0" or thisfields[ 2 ] == b"U1" or thisfields[ 2 ] == b"U2":
            # allow up to 2 mismatches...
            strand = thisfields[ 8 ]
            if strand == b"F":
                return ( chromname,
                         atoi( thisfields[ 7 ] ) - 1,
                         0 )
            elif strand == b"R":
                return ( chromname,
                         atoi( thisfields[ 7 ] ) + thistaglength - 1,
                         1 )
            else:
                raise StrandFormatError( thisline, strand )
        else:
            return ( b"", -1, -1 )

cdef class ELANDMultiParser( GenericParser ):
    """File Parser Class for ELAND multi File.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields:
    1. Sequence name
    2. Sequence
    3. Either NM, QC, RM (as described above) or the following:
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches
    found
    5. Blank, if no matches found or if too many matches found, or the following:
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction
    starting at position 160322 with one error, one in the forward direction starting at
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    cdef void __skip_first_commentlines ( self ):
        """ELANDResultParser needs to skip the first several comment lines.
        """
        cdef:
            int32_t l_line
            bytes this_line
        for thisline in self.fhd:
            l_line = len( thisline )
            if thisline and thisline[ 0 ] != 35: # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek( -l_line, 1 )
        return

    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( b'\t' )
        if thisfields[1].isdigit():
            return 0
        else:
            return len( thisfields[ 1 ] )

    cdef tuple __fw_parse_line ( self, bytes thisline ):
        cdef:
            list thisfields
            bytes thistagname, pos, strand
            int32_t thistaglength, thistaghits

        if not thisline: return ( b"", -1, -1 )
        thisline = thisline.rstrip()
        if not thisline: return ( b"", -1, -1 )

        thisfields = thisline.split( b'\t' )
        thistagname = thisfields[ 0 ]        # name of tag
        thistaglength = len( thisfields[ 1 ] ) # length of tag

        if len( thisfields ) < 4:
            return ( b"", -1, -1 )
        else:
            thistaghits = sum( [<int32_t>(x) for x in thisfields[ 2 ].split( b':' ) ] )
            if thistaghits > 1:
                # multiple hits
                return ( b"", -1, -1 )
            else:
                ( chromname, pos ) = thisfields[ 3 ].split( b':' )

                try:
                    chromname = chromname[ :chromname.rindex( b".fa" ) ]
                except ValueError:
                    pass

                strand  = pos[ -2 ]
                if strand == b"F":
                    return ( chromname,
                             <int32_t>( pos[ :-2 ] )-1,
                             0 )
                elif strand == b"R":
                    return ( chromname,
                             <int32_t>( pos[ :-2 ] ) + thistaglength - 1,
                             1 )
                else:
                    raise StrandFormatError( thisline,strand )


cdef class ELANDExportParser( GenericParser ):
    """File Parser Class for ELAND Export File.

    """
    cdef void __skip_first_commentlines ( self ):
        """ELANDResultParser needs to skip the first several comment lines.
        """
        cdef:
            int32_t l_line
            bytes this_line
        for thisline in self.fhd:
            l_line = len( thisline )
            if thisline and thisline[ 0 ] != 35: # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek( -l_line, 1 )
        return

    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( b'\t' )
        if len( thisfields ) > 12 and thisfields[ 12 ]:
            # a successful alignment has over 12 columns
            return len( thisfields[ 8 ] )
        else:
            return 0

    cdef tuple __fw_parse_line ( self, bytes thisline ):
        cdef:
            list thisfields
            bytes thisname, strand
            int32_t thistaglength

        thisline = thisline.rstrip()
        if not thisline: return ( b"", -1, -1 )

        thisfields = thisline.split( b"\t" )

        if len(thisfields) > 12 and thisfields[ 12 ]:
            thisname = b":".join( thisfields[ 0:6 ] )
            thistaglength = len( thisfields[ 8 ] )
            strand = thisfields[ 13 ]
            if strand == b"F":
                return ( thisfields[ 10 ], atoi( thisfields[ 12 ] ) - 1, 0 )
            elif strand == b"R":
                return ( thisfields[ 10 ], atoi( thisfields[ 12 ] ) + thistaglength - 1, 1 )
            else:
                raise StrandFormatError( thisline, strand )
        else:
            return ( b"", -1, -1 )

### Contributed by Davide, modified by Tao
cdef class SAMParser( GenericParser ):
    """File Parser Class for SAM File.

    Each line of the output file contains at least:
    1. Sequence name
    2. Bitwise flag
    3. Reference name
    4. 1-based leftmost position fo clipped alignment
    5. Mapping quality
    6. CIGAR string
    7. Mate Reference Name
    8. 1-based leftmost Mate Position
    9. Inferred insert size
    10. Query sequence on the same strand as the reference
    11. Query quality

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

    cdef void __skip_first_commentlines ( self ):
        """SAMParser needs to skip the first several comment lines.
        """
        cdef:
            int32_t l_line
            bytes this_line
        for thisline in self.fhd:
            l_line = len( thisline )
            if thisline and thisline[ 0 ] != 64: # 64 is b"@"
                break

        # rewind from SEEK_CUR
        self.fhd.seek( -l_line, 1 )
        return

    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Parse tag sequence, then tag length.

        """
        cdef:
            list thisfields
            int32_t bwflag

        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( b'\t' )
        bwflag = atoi( thisfields[ 1 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
            return 0       #unmapped sequence or bad sequence or 2nd or sup alignment
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag!
            if not bwflag & 2:
                return 0   # not a proper pair
            if bwflag & 8:
                return 0   # the mate is unmapped
            # From Benjamin Schiller https://github.com/benjschiller
            if bwflag & 128:
                # this is not the first read in a pair
                return 0
        return len( thisfields[ 9 ] )

    cdef tuple __fw_parse_line ( self, bytes thisline ):
        cdef:
            list thisfields
            bytes thistagname, thisref
            int32_t bwflag, thisstrand, thisstart

        thisline = thisline.rstrip()
        if not thisline: return ( b"", -1, -1 )
        thisfields = thisline.split( b'\t' )
        thistagname = thisfields[ 0 ]         # name of tag
        thisref = thisfields[ 2 ]
        bwflag = atoi( thisfields[ 1 ] )
        CIGAR = thisfields[ 5 ]

        if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)): return ( b"", -1, -1 )

        #if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
        #    return ( b"", -1, -1 )       #unmapped sequence or bad sequence or 2nd or sup alignment
        #if bwflag & 1:
        #    # paired read. We should only keep sequence if the mate is mapped
        #    # and if this is the left mate, all is within  the flag!
        #    if not bwflag & 2:
        #        return ( b"", -1, -1 )   # not a proper pair
        #    if bwflag & 8:
        #        return ( b"", -1, -1 )   # the mate is unmapped
        #    # From Benjamin Schiller https://github.com/benjschiller
        #    if bwflag & 128:
        #        # this is not the first read in a pair
        #        return ( b"", -1, -1 )
        #    # end of the patch
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            # minus strand, we have to decipher CIGAR string
            thisstrand = 1
            thisstart = atoi( thisfields[ 3 ] ) - 1 + sum( [ <int32_t>(x) for x in findall(b"(\d+)[MDNX=]",CIGAR) ] )	#reverse strand should be shifted alen bp
        else:
            thisstrand = 0
            thisstart = atoi( thisfields[ 3 ] ) - 1

        try:
            thisref = thisref[ :thisref.rindex( b".fa" ) ]
        except ValueError:
            pass

        return ( thisref, thisstart, thisstrand )

cdef class BAMParser( GenericParser ):
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
    def __init__ ( self, str filename, int64_t buffer_size = 100000 ):
        """Open input file. Determine whether it's a gzipped file.

        'filename' must be a string object.

        This function initialize the following attributes:

        1. self.filename: the filename for input file.
        2. self.gzipped: a boolean indicating whether input file is gzipped.
        3. self.fhd: buffered I/O stream of input file
        """
        self.filename = filename
        self.gzipped = True
        self.tag_size = -1
        self.buffer_size = buffer_size
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
            self.fhd = io.BufferedReader( gzip.open( filename, mode='rb' ), buffer_size = READ_BUFFER_SIZE) # buffersize set to 1M
        else:
            self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!

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

    cpdef int32_t tsize ( self ):
        """Get tag size from BAM file -- read l_seq field.

        Refer to: http://samtools.sourceforge.net/SAM1.pdf

        * This may not work for BAM file from bedToBAM (bedtools),
        since the l_seq field seems to be 0.
        """

        cdef:
            int32_t x, header_len, nc, nlength
            int32_t n = 0                   # successful read of tag size
            float64_t s = 0                # sum of tag sizes

        if self.tag_size != -1:
            # if we have already calculated tag size (!= -1),  return it.
            return self.tag_size

        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek( 4 )
        header_len =  unpack( '<i', fread( 4 ) )[ 0 ]
        fseek( header_len + ftell() )
        # get the number of chromosome
        nc = unpack( '<i', fread( 4 ) )[ 0 ]
        for x in range( nc ):
            # read each chromosome name
            nlength = unpack( '<i' , fread( 4 ) )[ 0 ]
            # jump over chromosome size, we don't need it
            fread( nlength )
            fseek( ftell() + 4 )
        while n < 10:
            entrylength = unpack( '<i', fread( 4 ) )[ 0 ]
            data = fread( entrylength )
            a = unpack( '<i', data[16:20] )[ 0 ]
            s += a
            n += 1
        fseek( 0 )
        self.tag_size = <int32_t>(s/n)
        return self.tag_size

    cpdef tuple get_references( self ):
        """
        read in references from BAM header

        return a tuple (references (list of names),
                        rlengths (dict of lengths)
        """
        cdef:
            int32_t header_len, x, nc, nlength
            bytes refname
            list references = []
            dict rlengths = {}

        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len =  unpack( '<i', fread( 4 ) )[ 0 ]
        fseek( header_len + ftell() )
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
        return (references, rlengths)

    cpdef build_fwtrack ( self ):
        """Build FWTrack from all lines, return a FWTrack object.

        Note only the unique match for a tag is kept.
        """
        cdef:
            int64_t i = 0                           #number of reads kept
            int32_t entrylength, fpos, strand, chrid
            list references
            dict rlengths

        fwtrack = FWTrack( buffer_size = self.buffer_size )
        references, rlengths = self.get_references() # after this, ptr at list of alignments
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell

        while True:
            try:
                entrylength = unpack( "<i", fread( 4 ) )[0]
            except struct.error:
                break
            ( chrid, fpos, strand ) = bam_se_entry_parser( fread( entrylength ) )
            if chrid == -1: continue
            fwtrack.add_loc( references[ chrid ], fpos, strand )
            i += 1
            if i % 1000000 == 0:
                info( " %d reads parsed" % i )

        #print( f"{references[chrid]:},{fpos:},{strand:}" )
        info( "%d reads have been read." % i )
        self.fhd.close()
        fwtrack.set_rlengths( rlengths )
        return fwtrack

    cpdef append_fwtrack ( self, fwtrack ):
        """Build FWTrack from all lines, return a FWTrack object.

        Note only the unique match for a tag is kept.
        """
        cdef:
            int64_t i = 0                     #number of reads kept
            int32_t entrylength, fpos, strand, chrid
            list references
            dict rlengths

        references, rlengths = self.get_references()
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell

        while True:
            try:
                entrylength = unpack( '<i', fread( 4 ) )[ 0 ]
            except struct.error:
                break
            ( chrid, fpos, strand ) = bam_se_entry_parser( fread( entrylength ) )
            if chrid == -1: continue
            fwtrack.add_loc( references[ chrid ], fpos, strand )
            i += 1
            if i % 1000000 == 0:
                info( " %d reads parsed" % i )

        info( "%d reads have been read." % i )
        self.fhd.close()
        #fwtrack.finalize()
        # this is the problematic part. If fwtrack is finalized, then it's impossible to increase the length of it in a step of buffer_size for multiple input files.
        fwtrack.set_rlengths( rlengths )
        return fwtrack

cdef class BAMPEParser(BAMParser):
    """File Parser Class for BAM File containing paired-end reads
    Only counts valid pairs, discards everything else
    Uses the midpoint of every read and calculates the average fragment size
    on the fly instead of modeling it

    File is gzip-compatible and binary.
    Information available is the same that is in SAM format.

    The bitwise flag is made like this:
    dec    meaning
    ---    -------
    1    paired read
    2    proper pair
    4    query unmapped
    8    mate unmapped
    16    strand of the query (1 -> reverse)
    32    strand of the mate
    64    first read in pair
    128    second read in pair
    256    alignment is not primary
    512    does not pass quality check
    1024    PCR or optical duplicate
    2048    supplementary alignment
    """
    cdef public int32_t n           # total number of fragments
    cdef public float32_t d         # the average length of fragments

    cpdef object build_petrack ( self ):
        """Build PETrackI from all lines, return a FWTrack object.
        """
        cdef:
            int64_t i          # number of fragments kept
            int64_t m          # sum of fragment lengths
            int32_t entrylength, fpos, chrid, tlen
            list references
            dict rlengths

        i = 0
        m = 0

        petrack = PETrackI( buffer_size = self.buffer_size )
        references, rlengths = self.get_references()
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell

        # for convenience, only count valid pairs
        add_loc = petrack.add_loc
        err = struct.error

        while True:
            try:
                entrylength = unpack( '<i', fread(4) )[0]
            except err:
                #e1 += 1
                break
            ( chrid, fpos, tlen ) = bampe_pe_entry_parser( fread(entrylength) )
            if chrid == -1:
                #e2 += 1
                continue
            add_loc(references[ chrid ], fpos, fpos + tlen)
            m += tlen
            i += 1
            if i % 1000000 == 0:
                info( " %d fragments parsed" % i )

        #print( f"{references[chrid]:},{fpos:},{tlen:}" )
        info( "%d fragments have been read." % i )
        #debug( f" {e1} Can't identify the length of entry, it may be the end of file, stop looping..." )
        #debug( f" {e2} Chromosome name can't be found which means this entry is skipped ..." )
        #assert i > 0, "Something went wrong, no fragment has been read! Check input file!"
        self.d = m / i
        self.n = i
        #assert self.d >= 0, "Something went wrong (mean fragment size was negative: %d = %d / %d)" % (self.d, m, i)
        self.fhd.close()
        petrack.set_rlengths( rlengths )
        return petrack

    cpdef append_petrack (self, petrack):
        """Build PETrackI from all lines, return a PETrackI object.
        """
        cdef:
            int64_t i = 0          # number of fragments kept
            int64_t m = 0          # sum of fragment lengths
            int32_t entrylength, fpos, chrid, tlen
            list references
            dict rlengths

        references, rlengths = self.get_references()
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell

        # for convenience, only count valid pairs
        add_loc = petrack.add_loc
        err = struct.error

        while True:
            try:
                entrylength = unpack('<i', fread(4))[0]
            except err:
                break
            ( chrid, fpos, tlen ) = bampe_pe_entry_parser( fread(entrylength) )
            if chrid == -1:
                continue
            add_loc(references[ chrid ], fpos, fpos + tlen)
            m += tlen
            i += 1
            if i % 1000000 == 0:
                info(" %d fragments parsed" % i)

        info( "%d fragments have been read." % i )
        self.d = ( self.d * self.n + m ) / ( self.n + i )
        self.n += i
        #assert self.d >= 0, "Something went wrong (mean fragment size was negative: %d = %d / %d)" % (self.d, m, i)
        self.fhd.close()
        # this is the problematic part. If fwtrack is finalized, then it's impossible to increase the length of it in a step of buffer_size for multiple input files.
        # petrack.finalize()
        petrack.set_rlengths( rlengths )
        return petrack

cdef class BowtieParser( GenericParser ):
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    cdef int32_t __tlen_parse_line ( self, bytes thisline ):
        """Parse tag sequence, then tag length.

        """
        cdef list thisfields

        thisline = thisline.rstrip()
        if not thisline: return ( b"", -1, -1 )
        if thisline[ 0 ] == b"#": return ( b"", -1 , -1 ) # comment line is skipped
        thisfields = thisline.split( b'\t' ) # I hope it will never bring me more trouble
        return len( thisfields[ 4 ] )

    cdef tuple __fw_parse_line (self, bytes thisline ):
        """
        The following definition comes from bowtie website:

        The bowtie aligner outputs each alignment on a separate
        line. Each line is a collection of 8 fields separated by tabs;
        from left to right, the fields are:

        1. Name of read that aligned

        2. Orientation of read in the alignment, - for reverse
        complement, + otherwise

        3. Name of reference sequence where alignment occurs, or
        ordinal ID if no name was provided

        4. 0-based offset into the forward reference strand where
        leftmost character of the alignment occurs

        5. Read sequence (reverse-complemented if orientation is -)

        6. ASCII-encoded read qualities (reversed if orientation is
        -). The encoded quality values are on the Phred scale and the
        encoding is ASCII-offset by 33 (ASCII char !).

        7. Number of other instances where the same read aligns
        against the same reference characters as were aligned against
        in this alignment. This is not the number of other places the
        read aligns with the same number of mismatches. The number in
        this column is generally not a good proxy for that number
        (e.g., the number in this column may be '0' while the number
        of other alignments with the same number of mismatches might
        be large). This column was previously described as 'Reserved'.

        8. Comma-separated list of mismatch descriptors. If there are
        no mismatches in the alignment, this field is empty. A single
        descriptor has the format offset:reference-base>read-base. The
        offset is expressed as a 0-based offset from the high-quality
        (5') end of the read.

        """
        cdef:
            list thisfields
            bytes chromname

        thisline = thisline.rstrip()
        if not thisline: return ( b"", -1, -1 )
        if thisline[ 0 ] == b"#": return ( b"", -1, -1 ) # comment line is skipped
        thisfields = thisline.split( b'\t' ) # I hope it will never bring me more trouble

        chromname = thisfields[ 2 ]
        try:
            chromname = chromname[ :chromname.rindex( b".fa" ) ]
        except ValueError:
            pass

            if thisfields[ 1 ] == b"+":
                return ( chromname,
                         atoi( thisfields[ 3 ] ),
                         0 )
            elif thisfields[ 1 ] == b"-":
                return ( chromname,
                         atoi( thisfields[ 3 ] ) + strlen( thisfields[ 4 ] ),
                         1 )
            else:
                raise StrandFormatError( thisline, thisfields[ 1 ] )

