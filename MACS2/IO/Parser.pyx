# Time-stamp: <2015-03-11 15:04:10 Tao Liu>

"""Module for all MACS Parser classes for input.

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import logging
import struct
from struct import unpack
from re import findall
import gzip
import io
from MACS2.Constants import *
from MACS2.IO.FixWidthTrack import FWTrack
from MACS2.IO.PairedEndTrack import PETrackI

cdef public bint HAS_PYSAM

try:
    import pysam
    HAS_PYSAM = True
except:
    HAS_PYSAM = False

HAS_PYSAM=False

from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    void free(void *ptr)
    int strcmp(char *a, char *b)
    char * strcpy(char *a, char *b)
    long atol(char *str)
    int atoi(char *str)

# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

cpdef guess_parser ( fhd, long buffer_size = 100000 ):
    order_list = ("BAM",
                  "BED",
                  "ELAND",
                  "ELANDMULTI",
                  "ELANDEXPORT",
                  "SAM",
                  "BOWTIE",
                  )

    for f in order_list:
        if f == 'BED':
            p = BEDParser( fhd, buffer_size = buffer_size )
        elif f == "ELAND":
            p = ELANDResultParser( fhd, buffer_size = buffer_size )
        elif f ==  "ELANDMULTI":
            p = ELANDMultiParser( fhd, buffer_size = buffer_size )
        elif f == "ELANDEXPORT":
            p = ELANDExportParser( fhd, buffer_size = buffer_size )
        elif f == "SAM":
            p = SAMParser( fhd, buffer_size = buffer_size )
        elif f == "BAM":
            p = BAMParser( fhd, buffer_size = buffer_size )
        elif f == "BOWTIE":
            p = BowtieParser( fhd, buffer_size = buffer_size )
        logging.debug( "Testing format %s" % f )
        s = p.sniff()
        if s:
            logging.info( "Detected format is: %s" % ( f ) )
            if p.gzipped:
                logging.info( "* Input file is gzipped." )
            return p
        else:
            p.close()
    raise Exception( "Can't detect format!" )

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
    cdef str filename
    cdef bool gzipped
    cdef int tag_size
    cdef object fhd
    cdef long buffer_size
    
    def __init__ ( self, str filename, long buffer_size = 100000 ):
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
            self.fhd = io.BufferedReader( gzip.open( filename, mode='rb' ) )
        else:
            self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!

    cpdef int tsize( self ):
        """General function to detect tag size.

        * Although it can be used by most parsers, it must be
          rewritten by BAMParser!
        """
        cdef:
            int s = 0
            int n = 0     # number of successful/valid read alignments
            int m = 0     # number of trials
            int this_taglength
        
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
        self.tag_size = s/n
        return self.tag_size

    cdef __tlen_parse_line ( self, str thisline ):
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
            long i, m, fpos, strand
            str chromosome
        
        fwtrack = FWTrack( buffer_size = self.buffer_size )
        i = 0
        m = 0
        for thisline in self.fhd:
            ( chromosome, fpos, strand ) = self.__fw_parse_line( thisline )
            i+=1
            if fpos < 0 or not chromosome:
                # normally __fw_parse_line will return -1 if the line
                # contains no successful alignment.
                continue
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i=0
            fwtrack.add_loc( chromosome, fpos, strand )

        # close fwtrack and sort
        fwtrack.finalize()
        # close file stream.
        self.close()
        return fwtrack

    cpdef append_fwtrack ( self, fwtrack ):
        """Add more records to an existing FWTrack object. 

        """
        i = 0
        m = 0
        for thisline in self.fhd:
            ( chromosome, fpos, strand ) = self.__fw_parse_line( thisline )
            i+=1
            if fpos < 0 or not chromosome:
                # normally __fw_parse_line will return -1 if the line
                # contains no successful alignment.
                continue
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i=0
            fwtrack.add_loc( chromosome, fpos, strand )

        # close fwtrack and sort
        fwtrack.finalize()
        self.close()
        return fwtrack
        
    cdef __fw_parse_line ( self, str thisline ):
        """Abstract function to parse chromosome, 5' end position and
        strand.
        
        """
        cdef str chromosome = ""
        cdef int fpos = -1
        cdef int strand = -1
        return ( chromosome, fpos, strand )

    cpdef sniff ( self ):
        """Detect whether this parser is the correct parser for input
        file.

        Rule: try to find the tag size using this parser, if error
        occurs or tag size is too small or too big, check is failed.

        * BAMParser has a different sniff function.
        """
        cdef int t
        
        try:
            t = self.tsize()
        except:
            self.fhd.seek( 0 )
            return False
        else:
            if t <= 10 or t >= 10000:
                self.fhd.seek( 0 )
                return False
            else:
                self.fhd.seek( 0 )
                return True
            
    cpdef close ( self ):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.fhd.close()

cdef class BEDParser( GenericParser ):
    """File Parser Class for tabular File.

    """
    cdef __tlen_parse_line ( self, str thisline ):
        """Parse 5' and 3' position, then calculate tag length.

        """
        thisline = thisline.rstrip()
        if not thisline \
           or thisline[ :5 ] == "track" \
           or thisline[ :7 ] == "browser"\
           or thisline[ 0 ] == "#":
            return 0

        thisfields = thisline.split( '\t' )
        return atoi( thisfields[ 2 ] )-atoi( thisfields[ 1 ] )
    
    cdef __fw_parse_line ( self, str thisline ):
        #cdef list thisfields
        cdef char * chromname
        
        thisline = thisline.rstrip()

        if not thisline or thisline[ :5 ] == "track" \
            or thisline[ :7 ] == "browser" \
            or thisline[ 0 ] == "#":
             return ( "", -1, -1 )

        thisfields = thisline.split( '\t' )
        chromname = thisfields[ 0 ]
        #try:
        ##    chromname = chromname[ :chromname.rindex( ".fa" ) ]
        #except ValueError:
        #    pass

        try:
            if not strcmp(thisfields[ 5 ],"+"):
                return ( chromname,
                         atoi( thisfields[ 1 ] ),
                         0 )
            elif not strcmp(thisfields[ 5 ], "-"):
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
            
cdef class ELANDResultParser( GenericParser ):
    """File Parser Class for tabular File.

    """
    cdef __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( '\t' )
        if thisfields[1].isdigit():
            return 0
        else:
            return len( thisfields[ 1 ] )

    cdef __fw_parse_line ( self, str thisline ):
        cdef:
            str chromname, strand
            int thistaglength
        #if thisline.startswith("#") or thisline.startswith("track") or thisline.startswith("browser"): return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )

        thisfields = thisline.split( '\t' )
        thistaglength = strlen( thisfields[ 1 ] )

        if len( thisfields ) <= 6:
            return ( "", -1, -1 )

        try:
            chromname = thisfields[ 6 ]
            chromname = chromname[ :chromname.rindex( ".fa" ) ]
        except ValueError:
            pass

        if thisfields[ 2 ] == "U0" or thisfields[ 2 ] == "U1" or thisfields[ 2 ] == "U2":
            # allow up to 2 mismatches...
            strand = thisfields[ 8 ]
            if strand == "F":
                return ( chromname,
                         atoi( thisfields[ 7 ] ) - 1,
                         0 )
            elif strand == "R":
                return ( chromname,
                         atoi( thisfields[ 7 ] ) + thistaglength - 1,
                         1 )
            else:
                raise StrandFormatError( thisline, strand )
        else:
            return ( "", -1, -1 )

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
    cdef __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( '\t' )
        if thisfields[1].isdigit():
            return 0
        else:
            return len( thisfields[ 1 ] )

    cdef __fw_parse_line ( self, str thisline ):
        cdef:
            list thisfields
            str thistagname, pos, strand
            int thistaglength, thistaghits
        
        if not thisline: return ( "", -1, -1 )
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )

        #if thisline[ 0 ] == "#": return ( "", -1, -1 ) # comment line is skipped
        
        thisfields = thisline.split( '\t' )
        thistagname = thisfields[ 0 ]        # name of tag
        thistaglength = len( thisfields[ 1 ] ) # length of tag

        if len( thisfields ) < 4:
            return ( "", -1, -1 )
        else:
            thistaghits = sum( map( int, thisfields[ 2 ].split( ':' ) ) )
            if thistaghits > 1:
                # multiple hits
                return ( "", -1, -1 )
            else:
                ( chromname, pos ) = thisfields[ 3 ].split( ':' )

                try:
                    chromname = chromname[ :chromname.rindex( ".fa" ) ]
                except ValueError:
                    pass
                
                strand  = pos[ -2 ]
                if strand == "F":
                    return ( chromname,
                             int( pos[ :-2 ] )-1,
                             0 )
                elif strand == "R":
                    return ( chromname,
                             int( pos[ :-2 ] ) + thistaglength - 1,
                             1 )
                else:
                    raise StrandFormatError( thisline,strand )


cdef class ELANDExportParser( GenericParser ):
    """File Parser Class for ELAND Export File.

    """
    cdef __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( '\t' )
        if len( thisfields ) > 12 and thisfields[ 12 ]:
            # a successful alignment has over 12 columns
            return len( thisfields[ 8 ] )
        else:
            return 0
        
    cdef __fw_parse_line ( self, str thisline ):
        cdef:
            list thisfields
            str thisname, strand
            int thistaglength
        
        #if thisline.startswith("#") : return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )
    
        thisfields = thisline.split( "\t" )

        if len(thisfields) > 12 and thisfields[ 12 ]:
            thisname = ":".join( thisfields[ 0:6 ] )
            thistaglength = len( thisfields[ 8 ] )
            strand = thisfields[ 13 ]
            if strand == "F":
                return ( thisfields[ 10 ], atoi( thisfields[ 12 ] ) - 1, 0 )
            elif strand == "R":
                return ( thisfields[ 10 ], atoi( thisfields[ 12 ] ) + thistaglength - 1, 1 )
            else:
                raise StrandFormatError( thisline, strand )
        else:
            return ( -1, -1, -1 )

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
    """

    cdef __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        cdef:
            list thisfields
            int bwflag
        
        thisline = thisline.rstrip()
        if not thisline: return 0
        if thisline[ 0 ] == "@": return 0 # header line started with '@' is skipped
        thisfields = thisline.split( '\t' )
        bwflag = atoi( thisfields[ 1 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 1024 or bwflag & 256 or bwflag & 2048:
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

    cdef __fw_parse_line ( self, str thisline ):
        cdef:
            list thisfields
            str thistagname, thisref
            int bwflag, thisstrand, thisstart

        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )
        if thisline[ 0 ] == "@": return ( "", -1, -1 ) # header line started with '@' is skipped
        thisfields = thisline.split( '\t' )
        thistagname = thisfields[ 0 ]         # name of tag
        thisref = thisfields[ 2 ]
        bwflag = atoi( thisfields[ 1 ] )
        CIGAR = thisfields[ 5 ]
        if bwflag & 4 or bwflag & 512 or bwflag & 1024 or bwflag & 256 or bwflag & 2048:
            return ( "", -1, -1 )       #unmapped sequence or bad sequence or 2nd or sup alignment
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return ( "", -1, -1 )   # not a proper pair
            if bwflag & 8:
                return ( "", -1, -1 )   # the mate is unmapped
            # From Benjamin Schiller https://github.com/benjschiller
            if bwflag & 128:
                # this is not the first read in a pair
                return ( "", -1, -1 )
            # end of the patch
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            # minus strand, we have to decipher CIGAR string
            
            thisstrand = 1
            thisstart = atoi( thisfields[ 3 ] ) - 1 + sum(map(int, findall("(\d+)[MDNX=]",CIGAR)))	#reverse strand should be shifted alen bp 
        else:
            thisstrand = 0
            thisstart = atoi( thisfields[ 3 ] ) - 1

        try:
            thisref = thisref[ :thisref.rindex( ".fa" ) ]
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
    """
    def __init__ ( self, str filename, long buffer_size = 100000 ):
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
        if HAS_PYSAM:
            self.fhd = pysam.Samfile(filename)
        else:
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
                self.fhd = io.BufferedReader( gzip.open( filename, mode='rb' ) )
            else:
                self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!

    cpdef sniff( self ):
        """Check the first 3 bytes of BAM file. If it's 'BAM', check
        is success.

        """
        if HAS_PYSAM:
            try:
                self.fhd.tell()
            except:
                return False
            else:
                return True
        else:
            magic_header = self.fhd.read( 3 )
            if magic_header == "BAM":
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

    cpdef int tsize ( self ):
        """Get tag size from BAM file -- read l_seq field.

        Refer to: http://samtools.sourceforge.net/SAM1.pdf

        * This may not work for BAM file from bedToBAM (bedtools),
        since the l_seq field seems to be 0.
        """
        
        cdef:
            int x, header_len, nc, nlength
            int n = 0                   # successful read of tag size
            double s = 0                # sum of tag sizes

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
        self.tag_size = int( s/n )
        return self.tag_size

    cpdef tuple get_references( self ):
        if HAS_PYSAM:
            return self.__get_references_w_pysam()
        else:
            return self.__get_references_wo_pysam()

    cdef tuple __get_references_w_pysam( self ):
        """
        read in references from BAM header
        
        return a tuple (references (list of names),
                        rlengths (dict of lengths)
        """
        cdef:
            list references = []
            dict rlengths = {}

        references = list(self.fhd.references)
        rlengths = dict( zip( references, self.fhd.lengths ) )
        return (references, rlengths)
        
    cdef tuple __get_references_wo_pysam( self ):
        """
        read in references from BAM header
        
        return a tuple (references (list of names),
                        rlengths (dict of lengths)
        """
        cdef:
            int header_len, x, nc, nlength
            str refname
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
        try:
            if HAS_PYSAM:
                return self.__build_fwtrack_w_pysam()
            else:
                return self.__build_fwtrack_wo_pysam()
        except IOError:
            logging.error( "BAM files might be corrupted!" )
            sys.exit(0)
            

    cdef __build_fwtrack_w_pysam ( self ):
        cdef:
            int i = 0
            int m = 0
            int fpos, strand, chrid
            list references
            dict rlengths
        
        fwtrack = FWTrack( buffer_size = self.buffer_size )
        self.fhd.reset()
        references, rlengths = self.get_references()
        while True:
            try:
                a  =  self.fhd.next()
            except StopIteration:
                break
            chrid = a.tid
            if a.is_unmapped or (a.is_paired and (not a.is_proper_pair or a.is_read2)):
                fpos = -1
            else:
                if a.is_reverse:
                    strand = 1              # minus strand
                    fpos = a.aend     # rightmost position
                else:
                    strand = 0              # plus strand
                    fpos = a.pos
            i+=1
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i = 0
            if fpos >= 0:
                fwtrack.add_loc( references[ chrid ], fpos, strand )
        self.fhd.close()
        fwtrack.finalize()
        fwtrack.set_rlengths( rlengths )
        return fwtrack

    cdef __build_fwtrack_wo_pysam ( self ):
        """Build FWTrack from all lines, return a FWTrack object.

        Note only the unique match for a tag is kept.
        """
        cdef:
            int i = 0
            int m = 0
            int entrylength, fpos, strand, chrid
            list references
            dict rlengths
        
        fwtrack = FWTrack( buffer_size = self.buffer_size )
        references, rlengths = self.get_references()
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        
        while True:
            try:
                entrylength = unpack( '<i', fread( 4 ) )[ 0 ]
            except struct.error:
                break
            ( chrid, fpos, strand ) = self.__fw_binary_parse_wo_pysam( fread( entrylength ) )
            i+=1
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i = 0
            if fpos >= 0:
                fwtrack.add_loc( references[ chrid ], fpos, strand )
        self.fhd.close()
        fwtrack.finalize()
        fwtrack.set_rlengths( rlengths )
        return fwtrack

    cpdef append_fwtrack ( self, fwtrack ):
        if HAS_PYSAM:
            return self.__append_fwtrack_w_pysam( fwtrack )
        else:
            return self.__append_fwtrack_wo_pysam( fwtrack )

    cdef __append_fwtrack_w_pysam ( self, fwtrack ):
        cdef:
            int i = 0
            int m = 0
            int fpos, strand, chrid
            list references
            dict rlengths
        
        self.fhd.reset()
        references, rlengths = self.get_references()
        while True:
            try:
                a  =  self.fhd.next()
            except StopIteration:
                break
            chrid = a.tid
            if a.is_unmapped or (a.is_paired and (not a.is_proper_pair or a.is_read2)):
                fpos = -1
            else:
                if a.is_reverse:
                    strand = 1              # minus strand
                    fpos = a.aend           # rightmost position
                else:
                    strand = 0              # plus strand
                    fpos = a.pos
            i+=1
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i = 0
            if fpos >= 0:
                fwtrack.add_loc( references[ chrid ], fpos, strand )
        self.fhd.close()
        fwtrack.finalize()
        if isinstance( fwtrack.rlengths, dict ):
            fwtrack.set_rlengths(rlengths)
        return fwtrack

    cdef __append_fwtrack_wo_pysam ( self, fwtrack ):
        """Build FWTrack from all lines, return a FWTrack object.

        Note only the unique match for a tag is kept.
        """
        cdef:
            int i = 0
            int m = 0
            int entrylength, fpos, strand, chrid
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
            ( chrid, fpos, strand ) = self.__fw_binary_parse_wo_pysam( fread( entrylength ) )
            i+=1
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i = 0
            if fpos >= 0:
                fwtrack.add_loc( references[ chrid ], fpos, strand )
        self.fhd.close()
        fwtrack.finalize()
        fwtrack.set_rlengths( rlengths )
        return fwtrack
    
    cdef tuple __fw_binary_parse_wo_pysam (self, data ):
        cdef:
            int thisref, thisstart, thisstrand, i
            short bwflag, l_read_name, n_cigar_op
            int cigar_code
        
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return ( -1, -1, -1 )

        thisref = unpack( '<i', data[ 0:4 ] )[ 0 ]
        thisstart = unpack( '<i', data[ 4:8 ] )[ 0 ]
        (n_cigar_op,  bwflag ) = unpack( '<HH' , data[ 12:16 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 1024 or bwflag & 256 or bwflag & 2048:
            return ( -1, -1, -1 )       #unmapped sequence or bad sequence or  secondary or supplementary alignment 
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return ( -1, -1, -1 )   # not a proper pair
            if bwflag & 8:
                return ( -1, -1, -1 )   # the mate is unmapped
            # From Benjamin Schiller https://github.com/benjschiller
            if bwflag & 128:
                # this is not the first read in a pair
                return ( -1, -1, -1 )
            # end of the patch
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            # read mapped to minus strand
            l_read_name = unpack( '<B', data[ 8:9 ] )[ 0 ]
            # need to decipher CIGAR string
            for cigar_code in unpack( '<%dI' % (n_cigar_op) , data[ 32 + l_read_name : 32 + l_read_name + n_cigar_op*4 ] ):
                if cigar_code & 15 in [ 0, 2, 3, 7, 8 ]:   # they are CIGAR op M/D/N/=/X
                    thisstart += cigar_code >> 4
            thisstrand = 1
        else:
            thisstrand = 0

        return ( thisref, thisstart, thisstrand )

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
    """
    cdef public int n
    cdef public int d

    cpdef build_petrack ( self ):
        if HAS_PYSAM:
            return self.__build_petrack_w_pysam()
        else:
            return self.__build_petrack_wo_pysam()

    cdef __build_petrack_w_pysam ( self ):
        return
#         cdef:
#             int i = 0
#             int m = 0
#             int entrylength, fpos, chrid, tlen
#             int *asint
#             list references
#             dict rlengths
#             float d = 0.0
#             char *rawread, *rawentrylength
#             _BAMPEParsed read
        
#         petrack = PETrackI()
#         self.fhd.reset()

#         references, rlengths = self.get_references()
#         while True:
#             try:
#                 a = self.fhd.next()
#             except StopIteration:
#                 break
#             chrid = a.tid
#             if a.is_paired:
#                 # not from paired-end sequencing? Pass
                
#             if a.is_unmapped or (a.is_paired and (not a.is_proper_pair or a.is_read2)):
#                 fpos = -1
#             else:
#                 if a.is_reverse:
#                     strand = 1              # minus strand
#                     fpos = a.pos + a.rlen     # rightmost position
#                 else:
#                     strand = 0              # plus strand
#                     fpos = a.pos
#             i+=1
#             if i == 1000000:
#                 m += 1
#                 logging.info( " %d" % ( m*1000000 ) )
#                 i = 0
#             if fpos >= 0:
#                 fwtrack.add_loc( references[ chrid ], fpos, strand )
#         self.fhd.close()
#         fwtrack.finalize()
#         fwtrack.set_rlengths( rlengths )
#         return fwtrack
            
            
#         # for convenience, only count valid pairs
#         add_loc = petrack.add_loc
#         info = logging.info
#         err = struct.error
#         while True:
#             try: entrylength = unpack('<i', fread(4))[0]
#             except err: break
#             rawread = <bytes>fread(32)
# #            rawread = <bytes>fread(entrylength)
#             read = self.__pe_binary_parse(rawread)
#             fseek(entrylength - 32, 1)
#             if read.ref == -1: continue
#             d = (d * i + abs(read.tlen)) / (i + 1) # keep track of avg fragment size
#             i+=1
            
#             if i == 1000000:
#                 m += 1
#                 info(" %d" % (m*1000000))
#                 i=0
#             petrack.add_loc(references[read.ref], read.start, read.start + read.tlen)
#         self.n = m * 1000000 + i
#         self.d = int(d)
#         assert d >= 0, "Something went wrong (mean fragment size was negative)"
#         self.fhd.close()
#         petrack.finalize()
#         petrack.set_rlengths( rlengths )
#         return petrack


    cdef __build_petrack_wo_pysam ( self ):
        """Build PETrackI from all lines, return a FWTrack object.
        """
        cdef:
            int i = 0
            int m = 0
            int entrylength, fpos, chrid, tlen
            int *asint
            list references
            dict rlengths
            float d = 0.0
            str rawread
            str rawentrylength
            _BAMPEParsed read
        
        petrack = PETrackI( buffer_size = self.buffer_size )

        references, rlengths = self.get_references()
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        
        # for convenience, only count valid pairs
        add_loc = petrack.add_loc
        info = logging.info
        err = struct.error
        while True:
            try: entrylength = unpack('<i', fread(4))[0]
            except err: break
            rawread = fread(32)
#            rawread = <bytes>fread(entrylength)
            read = self.__pe_binary_parse(rawread)
            fseek(entrylength - 32, 1)
            if read.ref == -1: continue
            d = (d * i + abs(read.tlen)) / (i + 1) # keep track of avg fragment size
            i+=1
            
            if i == 1000000:
                m += 1
                info(" %d" % (m*1000000))
                i=0
            add_loc(references[read.ref], read.start, read.start + read.tlen)
        self.n = m * 1000000 + i
        self.d = int(d)
        assert d >= 0, "Something went wrong (mean fragment size was negative)"
        self.fhd.close()
        petrack.finalize()
        petrack.set_rlengths( rlengths )
        return petrack

    cpdef append_petrack (self, petrack):
        """Build PETrackI from all lines, return a PETrackI object.
        """
        cdef:
            int i = 0
            int m = 0
            int entrylength, fpos, chrid, tlen
            int *asint
            list references
            dict rlengths
            float d = 0.0
            str rawread
            str rawentrylength
            _BAMPEParsed read
        
        references, rlengths = self.get_references()
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        
        # for convenience, only count valid pairs
        add_loc = petrack.add_loc
        info = logging.info
        err = struct.error
        while True:
            try: entrylength = unpack('<i', fread(4))[0]
            except err: break
            rawread = fread(32)
#            rawread = <bytes>fread(entrylength)
            read = self.__pe_binary_parse(rawread)
            fseek(entrylength - 32, 1)
            if read.ref == -1: continue
            d = (d * i + abs(read.tlen)) / (i + 1) # keep track of avg fragment size
            i+=1
            if i == 1000000:
                m += 1
                info(" %d" % (m*1000000))
                i=0
            add_loc(references[read.ref], read.start, read.start + read.tlen)
        self.n = m * 1000000 + i
        self.d = int(d)
        assert d >= 0, "Something went wrong (mean fragment size was negative)"
        self.fhd.close()
        petrack.finalize()
        petrack.set_rlengths( rlengths )
        return petrack
        
    cdef _BAMPEParsed __pe_binary_parse (self, str data):
        cdef:
            int nextpos, pos, cigar_op_len, i
            short bwflag, l_read_name, n_cigar_op, cigar_op
            _BAMPEParsed ret
#            int *asint = <int*>data
#            short *asshort = <short *>data
#            int thisref = asint[0]
#            int pos = asint[1]
#            short bwflag = asshort[7]
#            int nextpos = asint[6]
#            int tlen = asint[7]
        
        ret.ref = -1
        ret.start = -1
        ret.tlen = 0
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return ret

        (n_cigar_op,  bwflag ) = unpack( '<HH' , data[ 12:16 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 1024 or bwflag & 256 or bwflag & 2048:
            return ret       #unmapped sequence or bad sequence or 2nd or sup alignment
        if bwflag & 256 or bwflag & 2048:
            return ret          # secondary or supplementary alignment
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return ret  # not a proper pair
            if bwflag & 8:
                return ret  # the mate is unmapped
            if bwflag & 128:
                # this is not the first read in a pair
                return ret
                       
        ret.ref = unpack('<i', data[0:4])[0]
        pos = unpack('<i', data[4:8])[0]
        nextpos = unpack('<i', data[24:28])[0]
        ret.start = min(pos, nextpos) # we keep only the leftmost
                                      # position which means this must
                                      # be at + strand. So we don't
                                      # need to decipher CIGAR string.
        ret.tlen = abs(unpack('<i', data[28:32])[0]) # Actually, if
                                                     # the value
                                                     # unpacked is
                                                     # negative, then
                                                     # nextpos is the
                                                     # leftmost
                                                     # position.
        return ret

cdef struct _BAMPEParsed:
    int ref
    int start
    int tlen

### End ###

cdef class BowtieParser( GenericParser ):
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    cdef __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        cdef list thisfields
        
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )
        if thisline[ 0 ]=="#": return ( "", -1 , -1 ) # comment line is skipped
        thisfields = thisline.split( '\t' ) # I hope it will never bring me more trouble
        return len( thisfields[ 4 ] )

    cdef __fw_parse_line (self, str thisline ):
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
            str chromname
        
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )
        if thisline[ 0 ]=="#": return ( "", -1, -1 ) # comment line is skipped
        thisfields = thisline.split( '\t' ) # I hope it will never bring me more trouble

        chromname = thisfields[ 2 ]
        try:
            chromname = chromname[ :chromname.rindex( ".fa" ) ]
        except ValueError:
            pass

            if thisfields[ 1 ] == "+":
                return ( chromname,
                         atoi( thisfields[ 3 ] ),
                         0 )
            elif thisfields[ 1 ] == "-":
                return ( chromname,
                         atoi( thisfields[ 3 ] ) + strlen( thisfields[ 4 ] ),
                         1 )
            else:
                raise StrandFormatError( thisline, thisfields[ 1 ] )

cdef class PySAMParser:
    """Parser using PySAM to parse SAM or BAM

    """
    cdef str filename
    cdef bool gzipped
    cdef int tag_size
    cdef object fhd
    
    def __init__ ( self, str filename, long buffer_size = 100000 ):
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
            self.fhd = io.BufferedReader( gzip.open( filename, mode='rb' ) )
        else:
            self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!

    def tsize( self ):
        """General function to detect tag size.

        * Although it can be used by most parsers, it must be
          rewritten by BAMParser!
        """
        cdef:
            int s = 0
            int n = 0     # number of successful/valid read alignments
            int m = 0     # number of trials
            int this_taglength
        
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
        self.tag_size = s/n
        return self.tag_size

    cdef __tlen_parse_line ( self, str thisline ):
        """Abstract function to detect tag length.
        
        """
        raise NotImplemented
    
    def build_fwtrack ( self ):
        """Generic function to build FWTrack object. Create a new
        FWTrack object. If you want to append new records to an
        existing FWTrack object, try append_fwtrack function.

        * BAMParser for binary BAM format should have a different one.
        """
        cdef:
            long i, m, fpos, strand
            str chromosome
        
        fwtrack = FWTrack( buffer_size = self.buffer_size )
        i = 0
        m = 0
        for thisline in self.fhd:
            ( chromosome, fpos, strand ) = self.__fw_parse_line( thisline )
            i+=1
            if fpos < 0 or not chromosome:
                # normally __fw_parse_line will return -1 if the line
                # contains no successful alignment.
                continue
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i=0
            fwtrack.add_loc( chromosome, fpos, strand )

        # close fwtrack and sort
        fwtrack.finalize()
        # close file stream.
        self.close()
        return fwtrack

    def append_fwtrack ( self, fwtrack ):
        """Add more records to an existing FWTrack object. 

        """
        i = 0
        m = 0
        for thisline in self.fhd:
            ( chromosome, fpos, strand ) = self.__fw_parse_line( thisline )
            i+=1
            if fpos < 0 or not chromosome:
                # normally __fw_parse_line will return -1 if the line
                # contains no successful alignment.
                continue
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i=0
            fwtrack.add_loc( chromosome, fpos, strand )

        # close fwtrack and sort
        fwtrack.finalize()
        self.close()
        return fwtrack
        


    cdef __fw_parse_line ( self, str thisline ):
        """Abstract function to parse chromosome, 5' end position and
        strand.
        
        """
        cdef str chromosome = ""
        cdef int fpos = -1
        cdef int strand = -1
        return ( chromosome, fpos, strand )

    cpdef sniff ( self ):
        """Detect whether this parser is the correct parser for input
        file.

        Rule: try to find the tag size using this parser, if error
        occurs or tag size is too small or too big, check is failed.

        * BAMParser has a different sniff function.
        """
        cdef int t
        
        try:
            t = self.tsize()
        except:
            self.fhd.seek( 0 )
            return False
        else:
            if t <= 10 or t >= 10000:
                self.fhd.seek( 0 )
                return False
            else:
                self.fhd.seek( 0 )
                return True
            
    def close ( self ):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.fhd.close()
