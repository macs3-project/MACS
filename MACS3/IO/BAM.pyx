# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2022-02-02 00:04:16 Tao Liu>

"""

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import struct
from struct import unpack
import gzip
import io
import os
import sys
from zlib import decompress, MAX_WBITS

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
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t

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
cdef list get_bins_by_region ( uint32_t beg, uint32_t end ):
    """ Get the possible bins by given a region
    """
    cdef:
        list bins = [ 0 ]
        uint32_t k
    # for different levels
    for k in range(1 + (beg>>26), 2 + (end>>26) ):
        bins.append( k )
    for k in range(9 + (beg>>23), 10 + (end>>23) ):
        bins.append( k )
    for k in range(73 + (beg>>20), 74 + (end>>20) ):
        bins.append( k )
    for k in range(585 + (beg>>17), 586 + (end>>17) ):
        bins.append( k )
    for k in range(4681 + (beg>>14), 4682 + (end>>14) ):
        bins.append( k )
    return bins

cdef list reg2bins(rbeg, rend):
    """Generates bin ids which overlap the specified region.
    Args:
        rbeg (int): inclusive beginning position of region
        rend (int): exclusive end position of region
    Yields:
        (int): bin IDs for overlapping bins of region
    Raises:
        AssertionError (Exception): if the range is malformed or invalid
    """
    # Based off the algorithm presented in:
    # https://samtools.github.io/hts-specs/SAMv1.pdf

    # Bin calculation constants.
    BIN_ID_STARTS = (0, 1, 9, 73, 585, 4681)

    # Maximum range supported by specifications.
    MAX_RNG = (2 ** 29) - 1

    assert 0 <= rbeg <= rend <= MAX_RNG, 'Invalid region {}, {}'.format(rbeg, rend)

    l = []
    
    for start, shift in zip(BIN_ID_STARTS, range(29, 13, -3)):
        i = rbeg >> shift if rbeg > 0 else 0
        j = rend >> shift if rend < MAX_RNG else MAX_RNG >> shift

        for bin_id_offset in range(i, j + 1):
            #yield start + bin_id_offset
            l.append( start + bin_id_offset )
    return l
            
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


cdef class BAIFile:
    """BAI File Class for BAI (index of BAM) File.

    While initiating the object, BAI file will be loaded and the
    information of bins and chunks will be saved in the class object.

    Please refer to https://samtools.github.io/hts-specs/SAMv1.pdf for
    detail definition of BAI file.
    """
    cdef:
        str filename            # filename
        object fhd              # file handler
        bytes magic             # magic code for this file
        uint32_t n_ref          # number of refseqs/chromosomes
        dict metadata           # metadata for ref seqs
        uint64_t n_bins         # total number of bins
        uint64_t n_chunks       # total number of chunks
        uint64_t n_mapped       # total mapped reads
        uint64_t n_unmapped     # total unmapped reads
        list bins

    def __init__ ( self, str filename ):
        """Open input file. Determine whether it's a gzipped file.

        'filename' must be a string object.

        This function initialize the following attributes:

        1. self.filename: the filename for input file.
        2. self.gzipped: a boolean indicating whether input file is gzipped.
        3. self.fhd: buffered I/O stream of input file
        """
        self.filename = filename
        self.fhd = io.open( filename, mode='rb' ) # binary mode! I don't expect unicode here!
        self.magic = self.fhd.read( 4 )
        if self.magic != b"BAI\1":
            raise Exception(f"Not a BAI file. The first 4 bytes are \'{self.magic}\'")
        self.n_ref = self.__read_n_ref()
        self.bins = list( range( self.n_ref ) )
        self.metadata = dict.fromkeys( self.bins )
        self.__load_bins()
        self.fhd.close()
        return

    def __cinit__ ( self ):
        self.n_ref = 0
        self.n_bins = 0
        self.n_chunks = 0
        self.n_mapped = 0
        self.n_unmapped = 0

    def __str__ ( self ):
        tmp = list(self.bins[0].keys())[:10]
        return f"""
Summary of BAI File:
filename: {self.filename}
magic: {self.magic}
number of reference sequences/chromosomes: {self.n_ref}
number of total bins: {self.n_bins}
number of total chunks: {self.n_chunks}
number of total mapped reads: {self.n_mapped}
number of total unmapped reads: {self.n_unmapped}
Example of bins: ref 0, {tmp} ..
Example of metadata: ref 0, {self.metadata[ 0 ]}
"""

    cdef uint32_t __read_n_ref ( self ):
        cdef:
            uint32_t ret_val
        self.fhd.seek( 4 )      # skip magic
        ret_val = unpack( '<I', self.fhd.read( 4 ) )[ 0 ]
        return ret_val

    cdef __load_bins ( self ):
        cdef:
            uint32_t i, j, m
            uint32_t this_n_bin, this_n_intv
            uint32_t this_bin
            uint32_t this_n_chunk
            dict this_bin_dict
            list this_chunks
            list pseudobin
        # skip the magic and n_ref
        self.fhd.seek( 8 )
        # for each ref seq/chromosome
        for i in range( self.n_ref ):
            # get number of bins
            this_n_bin= unpack( '<I', self.fhd.read( 4 ) )[ 0 ]
            # must be less than 37451
            assert this_n_bin <= 37451
            # increment the total n_bins counter
            self.n_bins += this_n_bin
            # initiate this_bin_dict for this ref seq
            this_bin_dict = {}
            # for each bin
            for j in range( this_n_bin ):
                ( this_bin, this_n_chunk ) = unpack( '<II', self.fhd.read( 8 ) )
                # increment the total n_chunks counter
                self.n_chunks += this_n_chunk
                this_chunks = []
                for m in range( this_n_chunk ):
                    this_chunks.append( unpack( '<qq', self.fhd.read( 16 ) ) )
                # put the list of chunks in the this_bin_dict
                this_bin_dict[ this_bin ] = this_chunks
            # put the this_bin_dict in the list
            self.bins[ i ] = this_bin_dict
            # we will skip the next linear index part -- since I don't
            # think we need them in regular
            # ChIP-seq/ATAC-seq/other-seq analysis. Let me know if I am wrong!
            this_n_intv = unpack( '<I', self.fhd.read( 4 ) )[ 0 ]
            # skip each 64bits
            self.fhd.seek( 8 * this_n_intv, 1 )
        # we will also skip n_no_coor, the Number of unplaced unmapped reads
        # now extract and clean up the pseudo-bins
        # if empty, just skip
        for i in range( self.n_ref ):
            if self.bins[ i ]:
                pseudobin = self.bins[ i ].pop( 37450 )
                self.metadata[ i ] = {"ref_beg": pseudobin[0][0], "ref_end": pseudobin[0][1],
                                    "n_mapped": pseudobin[1][0], "n_unmapped": pseudobin[1][1] }
                self.n_mapped += pseudobin[1][0]
                self.n_unmapped += pseudobin[1][1]
        return

    cpdef list get_chunks_by_bin ( self, uint32_t ref_n, uint32_t bin_n ):
        """Get the chunks by bin number, for a given ref seq (not the name,
        but the index).

        The chunks will be sorted using default python sorted
        function. Therefore the result will be the order of the offset
        of beginning of each chunks.

        """
        return sorted( self.bins[ ref_n ].get( bin_n, [] ) )

    cpdef list get_chunks_by_list_of_bins ( self, uint32_t ref_n, list bins ):
        """Similar to get_chunks_by_bin, but accept a list of bins.

        Note: The redudant bins in the list will be removed.

        """
        cdef:
            uint32_t bin_n
            list chunks = []
            set bin_set
        bin_set = set( bins )
        for bin_n in bin_set:
            chunks.extend( self.bins[ ref_n ].get( bin_n, [] ) )
        return sorted( chunks )

    cpdef dict get_metadata_by_refseq ( self, uint32_t ref_n ):
        return self.metadata[ ref_n ]

    cpdef list get_chunks_by_region ( self, uint32_t ref_n, uint32_t beg, uint32_t end ):
        """Get the chunks by given a region in a given refseq (not the name,
        but the index).

        """
        cdef:
            list bins
            list chunks
        bins = get_bins_by_region( beg, end )
        chunks = self.get_chunks_by_list_of_bins( ref_n, bins )
        return chunks

    cpdef list get_chunks_by_list_of_regions ( self, uint32_t ref_n, list regions ):
        """ Similar to get_chunks_by_region, but accept a list of regions
        """
        cdef:
            int i
            list temp_bins
            list bins = []
            uint32_t beg, end

        for i in range( len(regions) ):
            beg = regions[ i ][ 0 ]
            end = regions[ i ][ 1 ]
            temp_bins = get_bins_by_region( beg, end )
            bins.extend( temp_bins )
        return self.get_chunks_by_list_of_bins( ref_n, bins )

    cpdef uint64_t get_coffset_by_region ( self, uint32_t ref_n, uint32_t beg, uint32_t end ):
        """Find the offset to access the BGZF block which contains the
        leftmost reads overlapping with the given genomic region.
        """
        cdef:
            uint64_t voffset_tmp, coffset_tmp
            list chunks
            tuple chunk
            int i
            uint64_t coffset

        chunks = self.get_chunks_by_region( ref_n, beg, end )
        if not chunks:
            return 0
        coffset = chunks[0][0] >> 16
        for i in range( 1, len(chunks) ):
            voffset_tmp = chunks[i][0]
            coffset_tmp = voffset_tmp >> 16
            if coffset_tmp < coffset:
                coffset = coffset_tmp
        return coffset

    cpdef uint64_t get_coffsets_by_list_of_regions ( self, uint32_t ref_n, list regions ):
        """Find the offset to access the BGZF block which contains the
        leftmost reads overlapping with the given genomic region.  
        """
        cdef:
            uint32_t beg, end
            int i
            uint64_t coffset
            list coffset_list

        coffset_list = []
        for i in range( len(regions) ):
            beg = regions[ i ][ 0 ]
            end = regions[ i ][ 1 ]            
            coffset = get_coffset_by_region( ref_n, beg, end )
            coffset_list.append( coffset )
        return coffset_list

cdef class BAMaccessor:
    """This BAM file reader is designed to access certain alignment
    records that overlap with givin genomic coordinates.

    The BAMaccessor needs a matching BAM index file (.bai) for fast
    access. Please refer to SAM format specification:

    https://samtools.github.io/hts-specs/SAMv1.pdf

    Note, the header of BAM will still be read through 'traditional'
    way -- using gzip to read through. But for accessing specific
    alignments, we will read BAM as binary, access the compressed
    chunk then use zlib.decompress.

    """
    cdef:
        # all private
        str bam_filename        # BAM filename
        str bai_filename        # BAI filename
        object bamfile          # BAM file handler "rb" mode
        BAIFile baifile         # Matching BAI file
        list references         # name of references/chromosomes, contain the order of chromosomes
        dict rlengths           # lengths of references/chromosomes
        bytes bgzf_block_cache  # cache of decompressed bgzf_block
        uint64_t coffset_cache  # coffset of the cached bgzf_block
        uint64_t noffset_cache  # coffset of the next block of the cached bgzf_block

    def __init__ ( self, str BAM_filename ):
        """Open input file. Determine whether it's a gzipped file.

        'filename' must be a string object.

        It will call __parse_header to check if the file is BAM format
        (through magic string); check if it's sorted by coordinates
        (SO). It will then check if BAI is available.

        """
        self.bam_filename = BAM_filename
        # parse the header
        self.__parse_header()
        # check if BAI is available
        self.bai_filename = self.bam_filename + ".bai"
        if os.path.exists( self.bai_filename ):
            self.baifile = BAIFile( self.bai_filename ) # The baifile is not for IO, it already contains all content.
        else:
            raise Exception(f"BAI is not available! Please make sure the `{self.bai_filename}` file exists in the same path")
        self.bamfile = io.BufferedReader( io.open( self.bam_filename, mode='rb' ), buffer_size = READ_BUFFER_SIZE )  # binary mode to read the raw BGZF file
        self.bgzf_block_cache = b""
        self.coffset_cache = 0
        self.noffset_cache = 0

    cpdef close ( self ):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.bamfile.close()

    cdef __parse_header( self ):
        """Parse the HEADER of BAM. It will do the following thing:

        1. read in references from BAM header, and save in
        self.references, and self.rlengths.  

        2. invoke self.check_sorted to see if BAM is sorted by
        coordinates, and set self.sorted.

        3. remember the file pointer at the end of HEADER or the start
        of the alignment blocks, and set self.SOA.
        
        """
        cdef:
            object fhd          # file handler from gzip.open; temporary use
            int header_len, x, nc, nlength
            bytes magic, refname
            list references = []
            dict rlengths = {}
            bytes header_text

        # we use traditional way to read the header -- through gzip
        fhd = io.BufferedReader( gzip.open( self.bam_filename, mode='rb' ), buffer_size = READ_BUFFER_SIZE )
        # get the first 4 bytes
        magic = fhd.read( 4 )
        if magic != b"BAM\1":
            raise Exception( f"File \"{self.filenmame}\" is not a BAM file. The magic string is not \"BAM\\1\", but \"{magic}\" " )
        # next 4bytes is length of header 
        header_len =  unpack( '<i', fhd.read( 4 ) )[ 0 ]
        header_text = unpack( '<%ds' % header_len, fhd.read( header_len ) )[0]
        # get the number of chromosome
        nc = unpack( '<i', fhd.read( 4 ) )[ 0 ]
        for x in range( nc ):
            # read each chromosome name
            nlength = unpack( '<i', fhd.read( 4 ) )[ 0 ]
            refname = fhd.read( nlength )[ :-1 ]
            references.append( refname )
            # don't jump over chromosome size
            # we can use it to avoid falling of chrom ends during peak calling
            rlengths[refname] = unpack( '<i', fhd.read( 4 ) )[ 0 ]
        self.references = references
        self.rlengths = rlengths
        # read sorting information
        if not self.__check_sorted( header_text ):
            raise Exception( f"BAM should be sorted by coordinates!" )
        # close
        fhd.close()
        return

    cdef bool __check_sorted( self, bytes header ):
        """Check if the BAM is sorted by looking at the HEADER text
        """
        cdef:
            bytes l
            int i

        # HD line is the first line
        l = header.split(b"\n")[0]
        if l.startswith(b"@HD"):
            # locate HD line
            i = l.find(b"SO:")
            if i != -1:
                # I will only check the first 5 characters
                if l[i+3:i+8] == b"coord":
                    return True
        return False
    
    cpdef list get_chromosomes ( self ):
        """Get chromosomes in order of their appearance in BAM HEADER.

        """
        return self.references

    cpdef dict get_rlengths ( self ):
        """Get chromosomes in order of their appearance in BAM HEADER.

        """
        return self.rlengths

    cdef tuple __decode_voffset ( self, uint64_t voffset ):
        """
        """
        cdef:
            uint64_t coffset
            uint32_t uoffset
        coffset = voffset >> 16
        uoffset = voffset ^ (coffset << 16)
        return ( coffset, uoffset )

    cdef bool __seek ( self, uint64_t offset ):
        """Move to given position
        """
        self.bamfile.seek( offset, 0 )
        return True
    
    cdef bool __retrieve_cdata_from_bgzf_block ( self ):
        """Retrieve uncompressed CDATA from BGZF block
        """
        cdef:
            uint16_t xlen, bsize
            bytes extra, cdata

        self.bamfile.seek( 10, 1 ) # skip 10 bytes
        # get XLEN
        xlen = unpack("<H", self.bamfile.read(2) )[0]
        # get extra subfields
        extra = self.bamfile.read(xlen)
        bsize = unpack( '<H', extra[extra.index(b'BC\x02\x00') + 4 : ])[0]
        cdata = self.bamfile.read( bsize - xlen - 19 )
        self.bgzf_block_cache = decompress( cdata, -MAX_WBITS )
        self.bamfile.seek( 8, 1 ) # move to the end of the block
        return True

    cpdef list get_reads_in_region ( self, bytes chrom, int left, int right, int maxDuplicate = 1 ):
        """Get reads in a given region.

        Return: list of ReadAlignment
        """
        cdef:
            int entrylength, fpos, strand, chrid
            list readslist
            int cur_duplicates = 0
            object read
            object previous_read
            int chrom_index
            bytes dcdata
            bool flag_end_searching = False
            uint64_t coffset

        readslist = []
        cur_duplicates = 0

        previous_read = None
        chrom_index = self.references.index( chrom )

        # get the coffset from BAI
        coffset = self.baifile.get_coffset_by_region ( chrom_index, left, right )
        if coffset == 0:
            return readslist
        if coffset != self.coffset_cache:
            self.coffset_cache = coffset
            self.__seek( coffset )
            self.__retrieve_cdata_from_bgzf_block ( ) # note: self.cached_bgzf_block will be updated.
        dcdata = self.bgzf_block_cache

        while True:
            dcdata_length = len(dcdata)
            # now parse the reads in dcdata
            i_bytes = 0
            while i_bytes < dcdata_length:
                entrylength = unpack( '<I', dcdata[ i_bytes : i_bytes + 4 ] )[ 0 ]
                i_bytes += 4
                # I am not sure if i_bytes+entrylength can be outside of dcdata_length...
                read = self.__fw_binary_parse( dcdata[ i_bytes : i_bytes + entrylength ] )
                i_bytes += entrylength

                if read == None:
                    # no mapping information, skip to next one
                    continue
                if read["lpos"] > right:
                    # outside of region, end searching
                    flag_end_searching = True
                    break
                elif read["rpos"] > left:
                    # found overlap
                    # check if redundant
                    if previous_read != None and previous_read["lpos"] == read["lpos"] and previous_read["rpos"] == read["rpos"] and previous_read["strand"] ==  read["strand"] and previous_read["cigar"] == read["cigar"]:
                        cur_duplicates += 1
                    else:
                        cur_duplicates = 1
                    if cur_duplicates <= maxDuplicate:
                        readslist.append( read )
                    previous_read = read
            if flag_end_searching:
                break
            # if not, get the next block, search again
            self.coffset_cache = self.bamfile.tell()
            self.__retrieve_cdata_from_bgzf_block ()
            dcdata = self.bgzf_block_cache
            if not dcdata:
                # empty means EOF is reached.
                break
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
            # Different with other MACS subcommands, both reads will be kept.
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


