# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2024-10-07 16:09:06 Tao Liu>

"""

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from struct import unpack
import gzip
import io
import os
import sys
from zlib import decompress, MAX_WBITS

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Utilities.Constants import READ_BUFFER_SIZE
from MACS3.Signal.ReadAlignment import ReadAlignment

# ------------------------------------
# Other modules
# ------------------------------------
import cython
from cython.cimports.cpython import bool

# ------------------------------------
# constants
# ------------------------------------
if sys.byteorder == "little":
    endian_prefix = "<"
elif sys.byteorder == "big":
    endian_prefix = ">"
else:
    raise Exception("Byteorder should be either big or little-endian.")

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.cfunc
def get_bins_by_region(beg: cython.uint, end: cython.uint) -> list:
    """ Get the possible bins by given a region
    """
    bins: list = [0]
    k: cython.uint

    # for different levels
    for k in range(1 + (beg >> 26), 2 + (end >> 26)):
        bins.append(k)
    for k in range(9 + (beg >> 23), 10 + (end >> 23)):
        bins.append(k)
    for k in range(73 + (beg >> 20), 74 + (end >> 20)):
        bins.append(k)
    for k in range(585 + (beg >> 17), 586 + (end >> 17)):
        bins.append(k)
    for k in range(4681 + (beg >> 14), 4682 + (end >> 14)):
        bins.append(k)
    return bins


@cython.cfunc
def reg2bins(rbeg, rend) -> list:
    """Generates bin ids which overlap the specified region.
    Args:
        rbeg (int): inclusive beginning position of region
        rend (int): exclusive end position of region
    Yields:
        (int): bin IDs for overlapping bins of region
    Raises:
        AssertionError (Exception): if the range is malformed or invalid
    """
    l_bins: list

    # Based off the algorithm presented in:
    # https://samtools.github.io/hts-specs/SAMv1.pdf

    # Bin calculation constants.
    BIN_ID_STARTS = (0, 1, 9, 73, 585, 4681)

    # Maximum range supported by specifications.
    MAX_RNG = (2 ** 29) - 1

    assert 0 <= rbeg <= rend <= MAX_RNG, 'Invalid region {}, {}'.format(rbeg, rend)

    l_bins = []

    for start, shift in zip(BIN_ID_STARTS, range(29, 13, -3)):
        i = rbeg >> shift if rbeg > 0 else 0
        j = rend >> shift if rend < MAX_RNG else MAX_RNG >> shift

        for bin_id_offset in range(i, j + 1):
            # yield start + bin_id_offset
            l_bins.append(start + bin_id_offset)
    return l_bins

# ------------------------------------
# Classes
# ------------------------------------


class StrandFormatError(Exception):
    """Exception about strand format error.

    Example:
    raise StrandFormatError('Must be F or R','X')
    """
    def __init__(self, string, strand):
        self.strand = strand
        self.string = string

    def __str__(self):
        return repr("Strand information can not be recognized in this line: \"%s\",\"%s\"" % (self.string, self.strand))


class MDTagMissingError(Exception):
    """Exception about missing MD tag

    Example:
    raise MDTagMissingError(name, aux)

    aux is the auxiliary data part.
    """
    def __init__(self, name, aux):
        self.name = name
        self.aux = aux

    def __str__(self):
        return repr("MD tag is missing! Please use \"samtools calmd\" command to add MD tags!\nName of sequence:%s\nCurrent auxiliary data section: %s" % (self.name.decode(), self.aux.decode()))


@cython.cclass
class BAIFile:
    """BAI File Class for BAI (index of BAM) File.

    While initiating the object, BAI file will be loaded and the
    information of bins and chunks will be saved in the class object.

    Please refer to https://samtools.github.io/hts-specs/SAMv1.pdf for
    detail definition of BAI file.
    """
    filename: str               # filename
    fhd: object                 # file handler
    magic: bytes                # magic code for this file
    n_ref: cython.uint             # number of refseqs/chromosomes
    metadata: dict              # metadata for ref seqs
    n_bins: cython.ulong            # total number of bins
    n_chunks: cython.ulong          # total number of chunks
    n_mapped: cython.ulong          # total mapped reads
    n_unmapped: cython.ulong        # total unmapped reads
    bins: list

    def __init__(self, filename: str):
        """Open input file. Determine whether it's a gzipped file.

        'filename' must be a string object.

        This function initialize the following attributes:

        1. self.filename: the filename for input file.
        2. self.gzipped: a boolean indicating whether input file is gzipped.
        3. self.fhd: buffered I/O stream of input file
        """
        self.filename = filename
        self.fhd = io.open(filename, mode='rb')  # binary mode! I don't expect unicode here!
        self.magic = self.fhd.read(4)
        if self.magic != b"BAI\1":
            raise Exception(f"Not a BAI file. The first 4 bytes are \'{self.magic}\'")
        self.n_ref = self.__read_n_ref()
        self.bins = list(range(self.n_ref))
        self.metadata = dict.fromkeys(self.bins)
        self.__load_bins()
        self.fhd.close()
        return

    def __cinit__(self):
        self.n_ref = 0
        self.n_bins = 0
        self.n_chunks = 0
        self.n_mapped = 0
        self.n_unmapped = 0

    def __str__(self):
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
Example of metadata: ref 0, {self.metadata[0]}
"""

    @cython.cfunc
    def __read_n_ref(self) -> cython.uint:
        ret_val: cython.uint

        self.fhd.seek(4)      # skip magic
        ret_val = unpack(endian_prefix+'I', self.fhd.read(4))[0]
        return ret_val

    @cython.cfunc
    def __load_bins(self):
        i: cython.uint
        j: cython.uint
        m: cython.uint
        this_n_bin: cython.uint
        this_n_intv: cython.uint
        this_bin: cython.uint
        this_n_chunk: cython.uint
        this_bin_dict: dict
        this_chunks: list
        pseudobin: list

        # skip the magic and n_ref
        self.fhd.seek(8)
        # for each ref seq/chromosome
        for i in range(self.n_ref):
            # get number of bins
            this_n_bin = unpack(endian_prefix+'I', self.fhd.read(4))[0]
            # must be less than 37451
            assert this_n_bin <= 37451
            # increment the total n_bins counter
            self.n_bins += this_n_bin
            # initiate this_bin_dict for this ref seq
            this_bin_dict = {}
            # for each bin
            for j in range(this_n_bin):
                (this_bin, this_n_chunk) = unpack(endian_prefix+'II', self.fhd.read(8))
                # increment the total n_chunks counter
                self.n_chunks += this_n_chunk
                this_chunks = []
                for m in range(this_n_chunk):
                    this_chunks.append(unpack(endian_prefix+'qq', self.fhd.read(16)))
                # put the list of chunks in the this_bin_dict
                this_bin_dict[this_bin] = this_chunks
            # put the this_bin_dict in the list
            self.bins[i] = this_bin_dict
            # we will skip the next linear index part -- since I don't
            # think we need them in regular
            # ChIP-seq/ATAC-seq/other-seq analysis. Let me know if I am wrong!
            this_n_intv = unpack(endian_prefix+'I', self.fhd.read(4))[0]
            # skip each 64bits
            self.fhd.seek(8 * this_n_intv, 1)
        # we will also skip n_no_coor, the Number of unplaced unmapped reads
        # now extract and clean up the pseudo-bins
        # if empty, just skip
        for i in range(self.n_ref):
            if self.bins[i]:
                pseudobin = self.bins[i].pop(37450)
                self.metadata[i] = {"ref_beg": pseudobin[0][0], "ref_end": pseudobin[0][1],
                                    "n_mapped": pseudobin[1][0], "n_unmapped": pseudobin[1][1]}
                self.n_mapped += pseudobin[1][0]
                self.n_unmapped += pseudobin[1][1]
        return

    @cython.ccall
    def get_chunks_by_bin(self, ref_n: cython.uint, bin_n: cython.uint) -> list:
        """Get the chunks by bin number, for a given ref seq (not the name,
        but the index).

        The chunks will be sorted using default python sorted
        function. Therefore the result will be the order of the offset
        of beginning of each chunks.

        """
        return sorted(self.bins[ref_n].get(bin_n, []))

    @cython.ccall
    def get_chunks_by_list_of_bins(self, ref_n: cython.uint, bins: list) -> list:
        """Similar to get_chunks_by_bin, but accept a list of bins.

        Note: The redudant bins in the list will be removed.

        """
        bin_n: cython.uint
        chunks: list = []
        bin_set: set

        bin_set = set(bins)
        for bin_n in bin_set:
            chunks.extend(self.bins[ref_n].get(bin_n, []))
        return sorted(chunks)

    @cython.ccall
    def get_metadata_by_refseq(self, ref_n: cython.uint) -> dict:
        return self.metadata[ref_n]

    @cython.ccall
    def get_chunks_by_region(self, ref_n: cython.uint, beg: cython.uint, end: cython.uint) -> list:
        """Get the chunks by given a region in a given refseq (not the name,
        but the index).

        """
        bins: list
        chunks: list

        bins = get_bins_by_region(beg, end)
        chunks = self.get_chunks_by_list_of_bins(ref_n, bins)
        return chunks

    @cython.ccall
    def get_chunks_by_list_of_regions(self, ref_n: cython.uint, regions: list) -> list:
        """ Similar to get_chunks_by_region, but accept a list of regions
        """
        i: int
        temp_bins: list
        bins: list = []
        beg: cython.uint
        end: cython.uint

        for i in range(len(regions)):
            beg = regions[i][0]
            end = regions[i][1]
            temp_bins = get_bins_by_region(beg, end)
            bins.extend(temp_bins)
        return self.get_chunks_by_list_of_bins(ref_n, bins)

    @cython.ccall
    def get_coffset_by_region(self, ref_n: cython.uint, beg: cython.uint, end: cython.uint) -> cython.ulong:
        """Find the offset to access the BGZF block which contains the
        leftmost reads overlapping with the given genomic region.
        """
        voffset_tmp: cython.ulong
        coffset_tmp: cython.ulong
        chunks: list
        i: cython.int
        coffset: cython.ulong

        chunks = self.get_chunks_by_region(ref_n, beg, end)
        if not chunks:
            return 0
        coffset = chunks[0][0] >> 16
        for i in range(1, len(chunks)):
            voffset_tmp = chunks[i][0]
            coffset_tmp = voffset_tmp >> 16
            if coffset_tmp < coffset:
                coffset = coffset_tmp
        return coffset

    @cython.ccall
    def get_coffsets_by_list_of_regions(self, ref_n: cython.uint, regions: list) -> cython.ulong:
        """Find the offset to access the BGZF block which contains the
        leftmost reads overlapping with the given genomic region.
        """
        beg: cython.uint
        end: cython.uint
        i: cython.int
        coffset: cython.ulong
        coffset_list: list

        coffset_list = []
        for i in range(len(regions)):
            beg = regions[i][0]
            end = regions[i][1]
            coffset = self.get_coffset_by_region(ref_n, beg, end)
            coffset_list.append(coffset)
        return coffset_list


@cython.cclass
class BAMaccessor:
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
    # all private
    bam_filename: str           # BAM filename
    bai_filename: str            # BAI filename
    bamfile: object             # BAM file handler "rb" mode
    baifile: BAIFile            # Matching BAI file
    # name of references/chromosomes, contain the order of chromosomes
    references: list
    rlengths: dict    # lengths of references/chromosomes
    bgzf_block_cache: bytes     # cache of decompressed bgzf_block
    coffset_cache: cython.ulong     # coffset of the cached bgzf_block
    # coffset of the next block of the cached bgzf_block
    noffset_cache: cython.ulong

    def __init__(self, BAM_filename: str):
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
        if os.path.exists(self.bai_filename):
            # The baifile is not for IO, it already contains all content.
            self.baifile = BAIFile(self.bai_filename)
        else:
            raise Exception(f"BAI is not available! Please make sure the `{self.bai_filename}` file exists in the same path")
        # binary mode to read the raw BGZF file
        self.bamfile = io.BufferedReader(io.open(self.bam_filename, mode='rb'),
                                         buffer_size=READ_BUFFER_SIZE)
        self.bgzf_block_cache = b""
        self.coffset_cache = 0
        self.noffset_cache = 0

    @cython.ccall
    def close(self):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.bamfile.close()

    @cython.cfunc
    def __parse_header(self):
        """Parse the HEADER of BAM. It will do the following thing:

        1. read in references from BAM header, and save in
        self.references, and self.rlengths.

        2. invoke self.check_sorted to see if BAM is sorted by
        coordinates, and set self.sorted.

        3. remember the file pointer at the end of HEADER or the start
        of the alignment blocks, and set self.SOA.

        """
        fhd: object          # file handler from gzip.open; temporary use
        header_len: cython.int
        x: cython.int
        nc: cython.int
        nlength: cython.int
        magic: bytes
        refname: bytes
        header_text: bytes
        references: list = []
        rlengths: dict = {}

        # we use traditional way to read the header -- through gzip
        fhd = io.BufferedReader(gzip.open(self.bam_filename, mode='rb'),
                                buffer_size=READ_BUFFER_SIZE)
        # get the first 4 bytes
        magic = fhd.read(4)
        if magic != b"BAM\1":
            raise Exception(f"File \"{self.filenmame}\" is not a BAM file. The magic string is not \"BAM\\1\", but \"{magic}\" ")
        # next 4bytes is length of header
        header_len = unpack(endian_prefix+'i', fhd.read(4))[0]
        header_text = unpack(endian_prefix+'%ds' % header_len, fhd.read(header_len))[0]
        # get the number of chromosome
        nc = unpack(endian_prefix+'i', fhd.read(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = unpack(endian_prefix+'i', fhd.read(4))[0]
            refname = fhd.read(nlength)[:-1]
            references.append(refname)
            # don't jump over chromosome size
            # we can use it to avoid falling of chrom ends during peak calling
            rlengths[refname] = unpack(endian_prefix+'i', fhd.read(4))[0]
        self.references = references
        self.rlengths = rlengths
        # read sorting information
        if not self.__check_sorted(header_text):
            raise Exception("BAM should be sorted by coordinates!")
        # close
        fhd.close()
        return

    @cython.cfunc
    def __check_sorted(self, header: bytes) -> bool:
        """Check if the BAM is sorted by looking at the HEADER text
        """
        tl: bytes
        i: cython.int

        # HD line is the first line
        tl = header.split(b"\n")[0]
        if tl.startswith(b"@HD"):
            # locate HD line
            i = tl.find(b"SO:")
            if i != -1:
                # I will only check the first 5 characters
                if tl[i+3:i+8] == b"coord":
                    return True
        return False

    @cython.ccall
    def get_chromosomes(self) -> list:
        """Get chromosomes in order of their appearance in BAM HEADER.

        """
        return self.references

    @cython.ccall
    def get_rlengths(self) -> dict:
        """Get chromosomes in order of their appearance in BAM HEADER.

        """
        return self.rlengths

    @cython.ccall
    def __decode_voffset(self, voffset: cython.ulong) -> tuple:
        """
        """
        coffset: cython.ulong
        uoffset: cython.uint

        coffset = voffset >> 16
        uoffset = voffset ^ (coffset << 16)
        return (coffset, uoffset)

    @cython.ccall
    def __seek(self, offset: cython.ulong) -> bool:
        """Move to given position
        """
        self.bamfile.seek(offset, 0)
        return True

    @cython.ccall
    def __retrieve_cdata_from_bgzf_block(self) -> bool:
        """Retrieve uncompressed CDATA from BGZF block
        """
        xlen: cython.ushort
        bsize: cython.ushort
        extra: bytes
        cdata: bytes

        self.bamfile.seek(10, 1)  # skip 10 bytes
        # get XLEN
        xlen = unpack(endian_prefix+"H", self.bamfile.read(2))[0]
        # get extra subfields
        extra = self.bamfile.read(xlen)
        bsize = unpack(endian_prefix+'H', extra[extra.index(b'BC\x02\x00') + 4:])[0]
        cdata = self.bamfile.read(bsize - xlen - 19)
        self.bgzf_block_cache = decompress(cdata, -MAX_WBITS)
        self.bamfile.seek(8, 1)  # move to the end of the block
        return True

    @cython.ccall
    def get_reads_in_region(self, chrom: bytes, left: cython.int,
                            right: cython.int, maxDuplicate: cython.int = 1) -> list:
        """Get reads in a given region.

        Return: list of ReadAlignment
        """
        entrylength: cython.int
        chrom_index: cython.int
        cur_duplicates: cython.int = 0
        readslist: list = []
        read: ReadAlignment
        previous_read: ReadAlignment
        dcdata: bytes
        flag_end_searching: bool = False
        coffset: cython.ulong

        previous_read = None
        chrom_index = self.references.index(chrom)

        # get the coffset from BAI
        coffset = self.baifile.get_coffset_by_region(chrom_index, left, right)
        if coffset == 0:
            return readslist
        if coffset != self.coffset_cache:
            self.coffset_cache = coffset
            self.__seek(coffset)
            self.__retrieve_cdata_from_bgzf_block()  # note: self.cached_bgzf_block will be updated.
        dcdata = self.bgzf_block_cache

        while True:
            dcdata_length = len(dcdata)
            # now parse the reads in dcdata
            i_bytes = 0
            while i_bytes < dcdata_length:
                entrylength = unpack(endian_prefix+'I', dcdata[i_bytes: i_bytes + 4])[0]
                i_bytes += 4
                # I am not sure if i_bytes+entrylength can be outside of dcdata_length...
                read = self.__fw_binary_parse(dcdata[i_bytes: i_bytes + entrylength])
                i_bytes += entrylength

                if read is None:
                    # no mapping information, skip to next one
                    continue
                if read["lpos"] > right:
                    # outside of region, end searching
                    flag_end_searching = True
                    break
                elif read["rpos"] > left:
                    # found overlap
                    # check if redundant
                    if previous_read is not None and previous_read["lpos"] == read["lpos"] and previous_read["rpos"] == read["rpos"] and previous_read["strand"] == read["strand"] and previous_read["cigar"] == read["cigar"]:
                        cur_duplicates += 1
                    else:
                        cur_duplicates = 1
                    if cur_duplicates <= maxDuplicate:
                        readslist.append(read)
                    previous_read = read
            if flag_end_searching:
                break
            # if not, get the next block, search again
            self.coffset_cache = self.bamfile.tell()
            self.__retrieve_cdata_from_bgzf_block()
            dcdata = self.bgzf_block_cache
            if not dcdata:
                # empty means EOF is reached.
                break
        return readslist

    @cython.ccall
    def __fw_binary_parse(self, data, min_MAPQ=1):
        """ Read information from an alignment block. Discard
        alignment below a minimum MAPQ score, default: 1.

        Return: ReadAlignment object with only selected information,
        including refname (chromosome), leftmost alignment location,
        strand information, sequence, quality, CIGAR and MD tag.

        Note: We do require MD tag exists in BAM file.
        """
        ref: cython.int
        leftmost: cython.int
        rightmost: cython.int
        i: cython.int
        j: cython.int
        l_seq: cython.int
        strand: cython.int
        bwflag: cython.ushort
        l_read_name: cython.ushort
        n_cigar_op: cython.ushort
        MAPQ: cython.ushort
        bin_bam: cython.ushort
        cigar_op: tuple  # op_len<<4|op, op: "MIDNSHP=X" -> 012345678
        read_name: bytes
        # note: for each byte, 1st base in the highest 4bit; 2nd in
        # the lowest 4bit. "=ACMGRSVTWYHKDBN" -> [0,15]
        seq: bytes
        qual: bytes
        tag: bytes
        MD: bytes

        if not data:
            return None

        # read number of CIGAR operators, Bitwise FLAG
        (n_cigar_op,  bwflag) = unpack(endian_prefix+'HH', data[12:16])

        # first, we will discard problematic reads
        if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
            # unmapped sequence or bad sequence or secondary or
            # supplementary alignment
            return None
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # Different with other MACS subcommands, both reads will be kept.
            if not bwflag & 2:
                return None   # not a proper pair
            if bwflag & 8:
                return None   # the mate is unmapped

        # read length of readname, MAPQ, and bin information
        (l_read_name, MAPQ, bin_bam) = unpack(endian_prefix+'BBH', data[8:12])

        # we will also discard reads having MAPQ lower than
        # min_MAPQ. MAPQ=255 means no alignment
        if MAPQ < min_MAPQ or MAPQ == 255:
            return None

        # Now we read other information

        # read ref name
        ref = unpack(endian_prefix+'i', data[0:4])[0]
        # read leftmost position of alignment
        leftmost = unpack(endian_prefix+'i', data[4:8])[0]
        # read length of query sequence
        l_seq = unpack(endian_prefix+'i', data[16:20])[0]
        # readname , skip next_refID, next_pos, tlen, which is 12
        # bytes or from the 32th byte
        read_name = unpack(endian_prefix+'%ds' % (l_read_name),
                           data[32: 32+l_read_name])[0][:-1]  # last byte is \x00
        # cigar_op, find the index i first.
        i = 32 + l_read_name
        cigar_op = unpack(endian_prefix+'%dI' % (n_cigar_op), data[i: i + n_cigar_op*4])
        # read sequence information
        i += n_cigar_op*4
        seq = unpack(endian_prefix+'%ds' % int((l_seq+1)/2), data[i: i + int((l_seq+1)/2)])[0]
        # read quality information. Note: the value in BAM is the
        # acutal Phred score, there is no +33!
        i += int((l_seq+1)/2)
        qual = unpack(endian_prefix+'%ds' % (l_seq), data[i: i + l_seq])[0]

        rightmost = leftmost
        for j in cigar_op:
            # they are CIGAR op M/D/N/=/X, no need to add I, S, H, or P
            if j & 15 in [0, 2, 3, 7, 8]:
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
        tag = data[i:]
        j = tag.find(b'MDZ')
        if j == -1:
            raise MDTagMissingError(data[32: 32+l_read_name], tag)
        MD = tag[j+3: tag[j:].find(b"\0") + j]

        # construct a ReadAlignment object and return
        return ReadAlignment(read_name, self.references[ref],
                             leftmost, rightmost, strand, seq, qual,
                             cigar_op, MD)
