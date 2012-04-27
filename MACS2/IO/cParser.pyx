# cython: profile=True
# Time-stamp: <2012-04-27 03:46:07 Tao Liu>

"""Module for all MACS Parser classes for input.

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
import logging
import struct
from struct import unpack
import gzip
import io
from MACS2.Constants import *
from MACS2.IO.cFixWidthTrack import FWTrackIII

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

def guess_parser ( fhd ):
    parser_dict = {"BED":BEDParser,
                   "ELAND":ELANDResultParser,
                   "ELANDMULTI":ELANDMultiParser,
                   "ELANDEXPORT":ELANDExportParser,
                   "SAM":SAMParser,
                   "BAM":BAMParser,
                   "BAMPE": BAMPEParser,
                   "BOWTIE":BowtieParser
                   }
    order_list = ("BAM",
                  "BED",
                  "ELAND",
                  "ELANDMULTI",
                  "ELANDEXPORT",
                  "SAM",
                  "BOWTIE",
                  )
    
    for f in order_list:
        p = parser_dict[ f ]( fhd )
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
        
class GenericParser:
    """Generic Parser class.

    Inherit this class to write your own parser. In most cases, you need to override:

    1. __tlen_parse_line which returns tag length of a line
    2.  __fw_parse_line which returns tuple of ( chromosome, 5'position, strand )
    """
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
        self.tag_size = -1
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
        cdef int s, n, m, this_taglength
        
        if self.tag_size != -1:
            # if we have already calculated tag size (!= -1),  return it.
            return self.tag_size
        
        s = 0
        n = 0                           # number of successful/valid read alignments
        m = 0                           # number of trials

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

    def __tlen_parse_line ( self, str thisline ):
        """Abstract function to detect tag length.
        
        """
        return 0
    
    def build_fwtrack ( self ):
        """Generic function to build FWTrackIII object. 

        * BAMParser for binary BAM format should have a different one.
        """
        cdef long i, m, fpos, strand
        cdef str chromosome
        
        fwtrack = FWTrackIII()
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

        fwtrack.finalize()
        # close file stream.
        self.close()
        return fwtrack

    def __fw_parse_line ( self, str thisline ):
        """Abstract function to parse chromosome, 5' end position and
        strand.
        
        """
        cdef str chromosome = ""
        cdef int fpos = -1
        cdef int strand = -1
        return ( chromosome, fpos, strand )

    def sniff ( self ):
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

class BEDParser( GenericParser ):
    """File Parser Class for tabular File.

    """
    def __tlen_parse_line ( self, str thisline ):
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
    
    def __fw_parse_line ( self, str thisline ):
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
            

class ELANDResultParser( GenericParser ):
    """File Parser Class for tabular File.

    """
    def __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( '\t' )
        return len( thisfields[ 1 ] )

    def __fw_parse_line ( self, str thisline ):
        cdef str chromname, strand
        cdef int thistaglength
        #if thisline.startswith("#") or thisline.startswith("track") or thisline.startswith("browser"): return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )

        thisfields = thisline.split( '\t' )
        thistaglength = len( thisfields[ 1 ] )

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

class ELANDMultiParser( GenericParser ):
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
    def __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline: return 0
        thisfields = thisline.split( '\t' )
        return len( thisfields[ 1 ] )

    def __fw_parse_line ( self, str thisline ):
        cdef list thisfields
        cdef str thistagname, pos, strand
        cdef int thistaglength, thistaghits
        
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


class ELANDExportParser( GenericParser ):
    """File Parser Class for ELAND Export File.

    """
    def __tlen_parse_line ( self, str thisline ):
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
        
    def __fw_parse_line ( self, str thisline ):
        cdef list thisfields
        cdef str thisname, strand
        cdef int thistaglength
        
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
class SAMParser( GenericParser ):
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

    def __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        cdef list thisfields
        cdef int bwflag
        
        thisline = thisline.rstrip()
        if not thisline: return 0
        if thisline[ 0 ] == "@": return 0 # header line started with '@' is skipped
        thisfields = thisline.split( '\t' )
        bwflag = atoi( thisfields[ 1 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return 0       #unmapped sequence or bad sequence
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

    def __fw_parse_line ( self, thisline ):
        cdef list thisfields
        cdef str thistagname, thisref
        cdef int bwflag, thisstrand, thisstart

        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )
        if thisline[ 0 ] == "@": return ( "", -1, -1 ) # header line started with '@' is skipped
        thisfields = thisline.split( '\t' )
        thistagname = thisfields[ 0 ]         # name of tag
        thisref = thisfields[ 2 ]
        bwflag = atoi( thisfields[ 1 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return ( "", -1, -1 )       #unmapped sequence or bad sequence
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
            thisstrand = 1
            thisstart = atoi( thisfields[ 3 ] ) - 1 + atoi( thisfields[ 9 ] )    #reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0
            thisstart = atoi( thisfields[ 3 ] ) - 1    

        try:
            thisref = thisref[ :thisref.rindex( ".fa" ) ]
        except ValueError:
            pass
        return ( thisref, thisstart, thisstrand )

class BAMParser( GenericParser ):
    """File Parser Class for BAM File.

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

    def sniff( self ):
        """Check the first 3 bytes of BAM file. If it's 'BAM', check
        is success.

        """
        magic_header = self.fhd.read( 3 )
        if magic_header == "BAM":
            tsize  = self.tsize()
            if tsize > 0:
                self.fhd.seek( 0 )
                return True
            else:
                self.fhd.seek( 0 )
                raise Exception( "File is not of a valid BAM format! %d" % tsize )
            return False
        else:
            self.fhd.seek( 0 )
            return False
            
    def tsize( self ):
        """Get tag size from BAM file -- read l_seq field.

        Refer to: http://samtools.sourceforge.net/SAM1.pdf

        * This may not work for BAM file from bedToBAM (bedtools),
        since the l_seq field seems to be 0.
        """
        
        cdef int x, header_len, nc, nlength, n
        cdef double s

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
        s = 0
        n = 0
        while n < 10:
            entrylength = unpack( '<i', fread( 4 ) )[ 0 ]
            data = fread( entrylength )
            a = unpack( '<i', data[16:20] )[ 0 ]
            s += a
            n += 1
        fseek( 0 )
        self.tag_size = int( s/n )
        return self.tag_size

    def build_fwtrack ( self ):
        """Build FWTrackIII from all lines, return a FWTrackIII object.

        Note only the unique match for a tag is kept.
        """
        cdef int i, m, header_len, nc, x, nlength
        cdef int entrylength, fpos, strand, chrid
        cdef list references
        
        fwtrack = FWTrackIII()
        i = 0
        m = 0
        references = []
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
            references.append( fread( nlength )[ :-1 ] )
            # jump over chromosome size, we don't need it
            fseek( ftell() + 4 )
        
        while True:
            try:
                entrylength = unpack( '<i', fread( 4 ) )[ 0 ]
            except struct.error:
                break
            ( chrid, fpos, strand ) = self.__fw_binary_parse( fread( entrylength ) )
            i+=1
            if i == 1000000:
                m += 1
                logging.info( " %d" % ( m*1000000 ) )
                i = 0
            if fpos >= 0:
                fwtrack.add_loc( references[ chrid ], fpos, strand )
        self.fhd.close()
        fwtrack.finalize()
        return fwtrack
    
    def __fw_binary_parse (self, data ):
        cdef int thisref, thisstart, thisstrand
        cdef short cigar, bwflag
        
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return ( -1, -1, -1 )

        thisref = unpack( '<i', data[ 0:4 ] )[ 0 ]
        thisstart = unpack( '<i', data[ 4:8 ] )[ 0 ]
        ( cigar, bwflag ) = unpack( '<hh' , data[ 12:16 ] )
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return ( -1, -1, -1 )       #unmapped sequence or bad sequence
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
        l = unpack( '<i', data[ 16:20 ] )[ 0 ]
        if bwflag & 16:
            thisstrand = 1
            thisstart = thisstart + unpack( '<i', data[ 16:20 ] )[ 0 ]    #reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0

        return ( thisref, thisstart, thisstrand )

### End ###
class BAMPEParser(BAMParser):
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
    def tsize(self):
        if hasattr(self, 'd') and self.d is not None:
            pass
        else:
            self.build_fwtrack()
        return self.d

    def __build_fwtrack (self):
        """Build FWTrackIII from all lines, return a FWTrackIII object.

        Note only the unique match for a tag is kept.
        """
        cdef int i, m, header_len, nc, x, nlength
        cdef int entrylength, fpos, strand, chrid
        cdef list references
        cdef float d
        
        fwtrack = FWTrackIII()
        i = 0
        m = 0
        references = []
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            references.append(fread(nlength)[:-1])
            # jump over chromosome size, we don't need it
            fseek(ftell() + 4)
        
        d = 0
        # for convenience, only count valid pairs
        while True:
            try:
                entrylength = struct.unpack('<i', fread(4))[0]
            except struct.error:
                break
            (chrid,fpos,strand,tlen) = self.__fw_binary_parse(fread(entrylength))
            if chrid < 0: continue
            d = (d * i + abs(tlen)) / (i + 1) # keep track of avg fragment size
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            fwtrack.add_loc(references[chrid],fpos,strand)
        self.d = int(d)
        assert d >= 0, "Something went wrong (average fragment size was negative)"
        self.fhd.close()
        self.fwtrack = fwtrack
        
    def build_fwtrack (self):
        if hasattr(self, 'fwtrack') and self.fwtrack is not None:
            pass
        else:
            self.__build_fwtrack()
        return self.fwtrack
    
    def __fw_binary_parse (self, data ):
        cdef int thisref, thisstart, thisstrand
        cdef short cigar, bwflag
        
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return (-1,-1,-1,0)

        thisref = struct.unpack('<i', data[0:4])[0]
        thisstart = struct.unpack('<i', data[4:8])[0]
        (cigar, bwflag) = struct.unpack('<hh', data[12:16])
        tlen = struct.unpack('<i', data[28:32])[0]
        midpoint = thisstart + tlen / 2
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (-1, -1, -1, 0)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return (-1, -1, -1, 0)  # not a proper pair
            if bwflag & 8:
                return (-1, -1, -1, 0)  # the mate is unmapped
            if bwflag & 128:
                # this is not the first read in a pair
                return (-1, -1, -1, 0)
                
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        thisstrand = bool(bwflag & 16)

        return (thisref, midpoint, thisstrand, tlen)

### End ###

class BowtieParser( GenericParser ):
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    def __tlen_parse_line ( self, str thisline ):
        """Parse tag sequence, then tag length.

        """
        cdef list thisfields
        
        thisline = thisline.rstrip()
        if not thisline: return ( "", -1, -1 )
        if thisline[ 0 ]=="#": return ( "", -1 , -1 ) # comment line is skipped
        thisfields = thisline.split( '\t' ) # I hope it will never bring me more trouble
        return len( thisfields[ 4 ] )

    def __fw_parse_line (self, str thisline ):
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
        cdef list thisfields
        cdef str chromname
        
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
                         atoi( thisfields[ 3 ] ) + atoi( thisfields[ 4 ] ),
                         1 )
            else:
                raise StrandFormatError( thisline, thisfields[ 1 ] )

