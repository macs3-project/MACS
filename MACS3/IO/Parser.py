# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2024-10-16 00:09:32 Tao Liu>

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

import cython
from cython.cimports.cpython import bool

from cython.cimports.libc.stdlib import atoi

from MACS3.Utilities.Constants import READ_BUFFER_SIZE
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Utilities.Logger import logging

logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

# ------------------------------------
# Other modules
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


@cython.ccall
def guess_parser(fname, buffer_size: cython.long = 100000):
    # Note: BAMPE and BEDPE can't be automatically detected.
    ordered_parser_dict = {"BAM": BAMParser,
                           "BED": BEDParser,
                           "ELAND": ELANDResultParser,
                           "ELANDMULTI": ELANDMultiParser,
                           "ELANDEXPORT": ELANDExportParser,
                           "SAM": SAMParser,
                           "BOWTIE": BowtieParser}

    for f in ordered_parser_dict:
        p = ordered_parser_dict[f]
        t_parser = p(fname, buffer_size=buffer_size)
        debug("Testing format %s" % f)
        s = t_parser.sniff()
        if s:
            info("Detected format is: %s" % (f))
            if t_parser.is_gzipped():
                info("* Input file is gzipped.")
            return t_parser
        else:
            t_parser.close()
    raise Exception("Can't detect format!")


@cython.cfunc
def bam_fw_binary_parse(data: cython.pointer[cython.const[cython.uchar]],
                        endian_prefix=endian_prefix) -> tuple:
    """Parse a BAM SE entry.

    If in little endian system, prefix should "<", otherwise,
    ">".

    """
    thisref: cython.int
    thisstart: cython.int
    thisstrand: cython.int

    i: cython.int

    cigar_code: cython.int
    cigar_op: tuple
    n_cigar_op: cython.ushort

    bwflag: cython.ushort
    l_read_name: cython.uchar   # length for the read name

    # we skip lot of the available information in data (i.e. tag name,
    # quality etc etc) no data, return, does it really happen without
    # raising struct.error?
    if not data:
        return (-1, -1, -1)

    (n_cigar_op,  bwflag) = unpack(endian_prefix+'HH', data[12:16])

    # we filter out unmapped sequence or bad sequence or secondary or
    # supplementary alignment. we filter out 2nd mate, not a proper
    # pair, mate is unmapped
    if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)):
        return (-1, -1, -1)

    (thisref, thisstart) = unpack(endian_prefix+'ii', data[0:8])

    # In case of paired-end we have now skipped all possible "bad"
    # pairs in case of proper pair we have skipped the rightmost
    # one... if the leftmost pair comes we can treat it as a single
    # read, so just check the strand and calculate its start
    # position... hope I'm right!
    if bwflag & 16:
        # read mapped to minus strand; then we have to compute cigar
        # to find the rightmost position

        (l_read_name, MAPQ) = unpack(endian_prefix+'BB', data[8:10])

        # need to decipher CIGAR string
        i = 32 + l_read_name
        cigar_op = unpack(endian_prefix+'%dI' % (n_cigar_op),
                          data[i: i + n_cigar_op*4])

        for cigar_code in cigar_op:
            if cigar_code & 15 in [0, 2, 3, 7, 8]:
                # they are CIGAR op M/D/N/=/X
                thisstart += cigar_code >> 4
        thisstrand = 1
    else:
        thisstrand = 0

    return (thisref, thisstart, thisstrand)


@cython.cfunc
def bampe_pe_binary_parse(data: cython.pointer[cython.const[cython.uchar]],
                          endian_prefix=endian_prefix) -> tuple:
    """Parse a BAMPE record.

    If in little endian system, prefix should "<", otherwise,
    ">".
    """
    thisref: cython.int
    thisstart: cython.int
    thistlen: cython.int
    nextpos: cython.int
    pos: cython.int
    bwflag: cython.ushort

    # we skip lot of the available information in data (i.e. tag name,
    # quality etc etc)
    if not data:
        return (-1, -1, -1)

    bwflag = unpack(endian_prefix+'H', data[14:16])[0]

    # we filter out unmapped, bad sequence, secondary/supplementary
    # alignment we filter out other mate of paired reads, not a proper
    # pair, or mate is unmapped
    if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)):
        return (-1, -1, -1)

    (thisref, pos) = unpack(endian_prefix+'ii', data[0:8])

    (nextpos, thistlen) = unpack(endian_prefix+'ii', data[24:32])

    # we keep only the leftmost position which means this must be at +
    # strand. So we don't need to decipher CIGAR string.
    thisstart = min(pos, nextpos)
    # Actually, if the value unpacked is negative, then nextpos is the
    # leftmost position.
    thistlen = abs(thistlen)

    return (thisref, thisstart, thistlen)

# ------------------------------------
# Classes
# ------------------------------------


class StrandFormatError(BaseException):
    """Exception about strand format error.

    Example:
    raise StrandFormatError('Must be F or R','X')
    """
    def __init__(self, string, strand):
        self.strand = strand
        self.string = string

    def __str__(self):
        return repr("Strand information can not be recognized in this line: \"%s\",\"%s\"" % (self.string, self.strand))


@cython.cclass
class GenericParser:
    """Generic Parser class.

    Inherit this class to write your own parser. In most cases, you
    need to override:

    1. tlen_parse_line which returns tag length of a line
    2.  fw_parse_line which returns tuple of (chromosome, 5'position, strand)

    """
    filename: str
    gzipped: bool
    tag_size: cython.int
    fhd: object
    buffer_size: cython.long

    def __init__(self, filename: str, buffer_size: cython.long = 100000):
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
        f = gzip.open(filename)
        try:
            f.read(10)
        except IOError:
            # not a gzipped file
            self.gzipped = False
        f.close()
        if self.gzipped:
            # open with gzip.open, then wrap it with BufferedReader!
            # buffersize set to 10M by default.
            self.fhd = io.BufferedReader(gzip.open(filename, mode='rb'),
                                         buffer_size=READ_BUFFER_SIZE)
        else:
            # binary mode! I don't expect unicode here!            
            self.fhd = io.open(filename, mode='rb')
        self.skip_first_commentlines()

    @cython.cfunc
    def skip_first_commentlines(self):
        """Some parser needs to skip the first several comment lines.

        Redefine this if it's necessary!
        """
        return

    @cython.ccall
    def tsize(self) -> cython.int:
        """General function to detect tag size.

        * Although it can be used by most parsers, it must be
          rewritten by BAMParser!
        """
        s: cython.int = 0
        n: cython.int = 0  # number of successful/valid read alignments
        m: cython.int = 0  # number of trials
        this_taglength: cython.int
        thisline: bytes

        if self.tag_size != -1:
            # if we have already calculated tag size (!= -1),  return it.
            return self.tag_size

        # try 10k times or retrieve 10 successfule alignments
        while n < 10 and m < 10000:
            m += 1
            thisline = self.fhd.readline()
            this_taglength = self.tlen_parse_line(thisline)
            if this_taglength > 0:
                # this_taglength == 0 means this line doesn't contain
                # successful alignment.
                s += this_taglength
                n += 1
        # done
        self.fhd.seek(0)
        self.skip_first_commentlines()
        if n != 0:              # else tsize = -1
            self.tag_size = cython.cast(cython.int, (s/n))
        return self.tag_size

    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Abstract function to detect tag length.

        """
        raise NotImplementedError

    @cython.ccall
    def build_fwtrack(self):
        """Generic function to build FWTrack object. Create a new
        FWTrack object. If you want to append new records to an
        existing FWTrack object, try append_fwtrack function.

        * BAMParser for binary BAM format should have a different one.
        """
        i: cython.long
        fpos: cython.long
        strand: cython.long
        chromosome: bytes
        tmp: bytes = b""

        fwtrack = FWTrack(buffer_size=self.buffer_size)
        i = 0
        while True:
            # for each block of input
            tmp += self.fhd.read(READ_BUFFER_SIZE)
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[-1]
            for thisline in lines[:-1]:
                (chromosome, fpos, strand) = self.fw_parse_line(thisline)
                if fpos < 0 or not chromosome:
                    # normally fw_parse_line will return -1 if the line
                    # contains no successful alignment.
                    continue
                i += 1
                if i % 1000000 == 0:
                    info(" %d reads parsed" % i)
                fwtrack.add_loc(chromosome, fpos, strand)
        # last one
        if tmp:
            (chromosome, fpos, strand) = self.fw_parse_line(tmp)
            if fpos >= 0 and chromosome:
                i += 1
                fwtrack.add_loc(chromosome, fpos, strand)
        # close file stream.
        self.close()
        return fwtrack

    @cython.ccall
    def append_fwtrack(self, fwtrack):
        """Add more records to an existing FWTrack object.

        """
        i: cython.long
        fpos: cython.long
        strand: cython.long
        chromosome: bytes
        tmp: bytes = b""

        i = 0
        while True:
            # for each block of input
            tmp += self.fhd.read(READ_BUFFER_SIZE)
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[-1]
            for thisline in lines[:-1]:
                (chromosome, fpos, strand) = self.fw_parse_line(thisline)
                if fpos < 0 or not chromosome:
                    # normally fw_parse_line will return -1 if the line
                    # contains no successful alignment.
                    continue
                i += 1
                if i % 1000000 == 0:
                    info(" %d reads parsed" % i)
                fwtrack.add_loc(chromosome, fpos, strand)

        # last one
        if tmp:
            (chromosome, fpos, strand) = self.fw_parse_line(tmp)
            if fpos >= 0 and chromosome:
                i += 1
                fwtrack.add_loc(chromosome, fpos, strand)
        # close file stream.
        self.close()
        return fwtrack

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        """Abstract function to parse chromosome, 5' end position and
        strand.

        """
        chromosome: bytes = b""
        fpos: cython.int = -1
        strand: cython.int = -1
        return (chromosome, fpos, strand)

    @cython.ccall
    def sniff(self):
        """Detect whether this parser is the correct parser for input
        file.

        Rule: try to find the tag size using this parser, if error
        occurs or tag size is too small or too big, check is failed.

        * BAMParser has a different sniff function.
        """
        t: cython.int

        t = self.tsize()
        if t <= 10 or t >= 10000:  # tsize too small or too big
            self.fhd.seek(0)
            return False
        else:
            self.fhd.seek(0)
            self.skip_first_commentlines()
            return True

    @cython.ccall
    def close(self):
        """Run this when this Parser will be never used.

        Close file I/O stream.
        """
        self.fhd.close()

    @cython.ccall
    def is_gzipped(self) -> bool:
        return self.gzipped


@cython.cclass
class BEDParser(GenericParser):
    """File Parser Class for BED File.

    """

    @cython.cfunc
    def skip_first_commentlines(self):
        """BEDParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and (thisline[:5] != b"track") \
               and (thisline[:7] != b"browser") \
               and (thisline[0] != 35):  # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Parse 5' and 3' position, then calculate frag length.

        """
        thisline = thisline.rstrip()
        if not thisline:
            return 0

        thisfields = thisline.split(b'\t')
        return atoi(thisfields[2]) - atoi(thisfields[1])

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        # cdef list thisfields
        chromname: bytes
        thisfields: list

        thisline = thisline.rstrip()
        thisfields = thisline.split(b'\t')
        chromname = thisfields[0]
        try:
            if thisfields[5] == b"+":
                return (chromname,
                        atoi(thisfields[1]),
                        0)
            elif thisfields[5] == b"-":
                return (chromname,
                        atoi(thisfields[2]),
                        1)
            else:
                raise StrandFormatError(thisline, thisfields[5])
        except IndexError:
            # default pos strand if no strand
            # info can be found
            return (chromname,
                    atoi(thisfields[1]),
                    0)


@cython.cclass
class BEDPEParser(GenericParser):
    """Parser for BED format file containing PE information, and also
    can be used for the cases when users predefine the fragment
    locations by shifting/extending by themselves.

    Format:

    chromosome_name	frag_leftend	frag_rightend


    Note: Only the first columns are used!

    """
    n = cython.declare(cython.int, visibility='public')
    d = cython.declare(cython.float, visibility='public')

    @cython.cfunc
    def skip_first_commentlines(self):
        """BEDPEParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and (thisline[:5] != b"track") \
               and (thisline[:7] != b"browser") \
               and (thisline[0] != 35):  # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def pe_parse_line(self, thisline: bytes):
        """ Parse each line, and return chromosome, left and right positions

        """
        thisfields: list

        thisline = thisline.rstrip()

        # still only support tabular as delimiter.
        thisfields = thisline.split(b'\t')
        try:
            return (thisfields[0],
                    atoi(thisfields[1]),
                    atoi(thisfields[2]))
        except IndexError:
            raise Exception("Less than 3 columns found at this line: %s\n" % thisline)

    @cython.ccall
    def build_petrack(self):
        """Build PETrackI from all lines.

        """
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths
        tmp: bytes = b""

        petrack = PETrackI(buffer_size=self.buffer_size)
        add_loc = petrack.add_loc

        while True:
            # for each block of input
            tmp += self.fhd.read(READ_BUFFER_SIZE)
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[-1]
            for thisline in lines[:-1]:
                (chromosome, left_pos, right_pos) = self.pe_parse_line(thisline)
                if left_pos < 0 or not chromosome:
                    continue
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                m += right_pos - left_pos
                i += 1
                if i % 1000000 == 0:
                    info(" %d fragments parsed" % i)
                add_loc(chromosome, left_pos, right_pos)
        # last one
        if tmp:
            (chromosome, left_pos, right_pos) = self.pe_parse_line(thisline)
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                add_loc(chromosome, left_pos, right_pos)

        self.d = cython.cast(cython.float, m) / i
        self.n = i
        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"

        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack

    @cython.ccall
    def append_petrack(self, petrack):
        """Build PETrackI from all lines, return a PETrackI object.
        """
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths
        tmp: bytes = b""

        add_loc = petrack.add_loc
        while True:
            # for each block of input
            tmp += self.fhd.read(READ_BUFFER_SIZE)
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[-1]
            for thisline in lines[:-1]:
                (chromosome, left_pos, right_pos) = self.pe_parse_line(thisline)
                if left_pos < 0 or not chromosome:
                    continue
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                m += right_pos - left_pos
                i += 1
                if i % 1000000 == 0:
                    info(" %d fragments parsed" % i)
                add_loc(chromosome, left_pos, right_pos)
        # last one
        if tmp:
            (chromosome, left_pos, right_pos) = self.pe_parse_line(thisline)
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                add_loc(chromosome, left_pos, right_pos)

        self.d = (self.d * self.n + m) / (self.n + i)
        self.n += i

        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"
        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack


@cython.cclass
class ELANDResultParser(GenericParser):
    """File Parser Class for tabular File.

    """

    @cython.cfunc
    def skip_first_commentlines(self):
        """ELANDResultParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and thisline[0] != 35:  # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Parse tag sequence, then tag length.

        """
        thisfields: list

        thisline = thisline.rstrip()
        if not thisline:
            return 0
        thisfields = thisline.split(b'\t')
        if thisfields[1].isdigit():
            return 0
        else:
            return len(thisfields[1])

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        chromname: bytes
        strand: bytes
        thistaglength: cython.int
        thisfields: list

        thisline = thisline.rstrip()
        if not thisline:
            return (b"", -1, -1)

        thisfields = thisline.split(b'\t')
        thistaglength = len(thisfields[1])

        if len(thisfields) <= 6:
            return (b"", -1, -1)

        try:
            chromname = thisfields[6]
            chromname = chromname[:chromname.rindex(b".fa")]
        except ValueError:
            pass

        if thisfields[2] == b"U0" or thisfields[2] == b"U1" or thisfields[2] == b"U2":
            # allow up to 2 mismatches...
            strand = thisfields[8]
            if strand == b"F":
                return (chromname,
                        atoi(thisfields[7]) - 1,
                        0)
            elif strand == b"R":
                return (chromname,
                        atoi(thisfields[7]) + thistaglength - 1,
                        1)
            else:
                raise StrandFormatError(thisline, strand)
        else:
            return (b"", -1, -1)


@cython.cclass
class ELANDMultiParser(GenericParser):
    """File Parser Class for ELAND multi File.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields:
    1. Sequence name
    2. Sequence
    3. Either NM, QC, RM (as described above) or the following:
    4. x:y:z where x, y, and z are the number of exact, single-error, and
    2-error matches found
    5. Blank, if no matches found or if too many matches found, or the
    following:

    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1

    This says there are two matches to BAC_plus_vector.fa: one in the
    reverse direction starting at position 160322 with one error, one
    in the forward direction starting at position 170128 with two
    errors. There is also a single-error match to E_coli.fa.

    """

    @cython.cfunc
    def skip_first_commentlines(self):
        """ELANDResultParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and thisline[0] != 35:  # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline:
            return 0
        thisfields = thisline.split(b'\t')
        if thisfields[1].isdigit():
            return 0
        else:
            return len(thisfields[1])

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        # thistagname: bytes
        pos: bytes
        strand: bytes
        thisfields: list
        thistaglength: cython.int
        thistaghits: cython.int

        if not thisline:
            return (b"", -1, -1)
        thisline = thisline.rstrip()
        if not thisline:
            return (b"", -1, -1)

        thisfields = thisline.split(b'\t')
        # thistagname = thisfields[0]        # name of tag
        thistaglength = len(thisfields[1])  # length of tag

        if len(thisfields) < 4:
            return (b"", -1, -1)
        else:
            thistaghits = sum([cython.cast(cython.int, x) for x in thisfields[2].split(b':')])
            if thistaghits > 1:
                # multiple hits
                return (b"", -1, -1)
            else:
                (chromname, pos) = thisfields[3].split(b':')

                try:
                    chromname = chromname[:chromname.rindex(b".fa")]
                except ValueError:
                    pass

                strand = pos[-2]
                if strand == b"F":
                    return (chromname,
                            cython.cast(cython.int, pos[:-2]) - 1,
                            0)
                elif strand == b"R":
                    return (chromname,
                            cython.cast(cython.int, pos[:-2]) + thistaglength - 1,
                            1)
                else:
                    raise StrandFormatError(thisline, strand)


@cython.cclass
class ELANDExportParser(GenericParser):
    """File Parser Class for ELAND Export File.

    """
    @cython.cfunc
    def skip_first_commentlines(self):
        """ELANDResultParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and thisline[0] != 35:  # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Parse tag sequence, then tag length.

        """
        thisline = thisline.rstrip()
        if not thisline:
            return 0
        thisfields = thisline.split(b'\t')
        if len(thisfields) > 12 and thisfields[12]:
            # a successful alignment has over 12 columns
            return len(thisfields[8])
        else:
            return 0

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        # thisname: bytes
        strand: bytes
        thistaglength: cython.int
        thisfields: list

        thisline = thisline.rstrip()
        if not thisline:
            return (b"", -1, -1)

        thisfields = thisline.split(b"\t")

        if len(thisfields) > 12 and thisfields[12]:
            # thisname = b":".join(thisfields[0:6])
            thistaglength = len(thisfields[8])
            strand = thisfields[13]
            if strand == b"F":
                return (thisfields[10], atoi(thisfields[12]) - 1, 0)
            elif strand == b"R":
                return (thisfields[10], atoi(thisfields[12]) + thistaglength - 1, 1)
            else:
                raise StrandFormatError(thisline, strand)
        else:
            return (b"", -1, -1)


# Contributed by Davide, modified by Tao

@cython.cclass
class SAMParser(GenericParser):
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

    @cython.cfunc
    def skip_first_commentlines(self):
        """SAMParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and thisline[0] != 64:  # 64 is b"@"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Parse tag sequence, then tag length.

        """
        thisfields: list
        bwflag: cython.int

        thisline = thisline.rstrip()
        if not thisline:
            return 0
        thisfields = thisline.split(b'\t')
        bwflag = atoi(thisfields[1])
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
        return len(thisfields[9])

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        # thistagname: bytes
        thisref: bytes
        thisfields: list
        bwflag: cython.int
        thisstrand: cython.int
        thisstart: cython.int

        thisline = thisline.rstrip()
        if not thisline:
            return (b"", -1, -1)
        thisfields = thisline.split(b'\t')
        # thistagname = thisfields[0]         # name of tag
        thisref = thisfields[2]
        bwflag = atoi(thisfields[1])
        CIGAR = thisfields[5]

        if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)):
            return (b"", -1, -1)

        # if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
        #    return (b"", -1, -1)       #unmapped sequence or bad sequence or 2nd or sup alignment
        # if bwflag & 1:
        #    # paired read. We should only keep sequence if the mate is mapped
        #    # and if this is the left mate, all is within  the flag!
        #    if not bwflag & 2:
        #        return (b"", -1, -1)   # not a proper pair
        #    if bwflag & 8:
        #        return (b"", -1, -1)   # the mate is unmapped
        #    # From Benjamin Schiller https://github.com/benjschiller
        #    if bwflag & 128:
        #        # this is not the first read in a pair
        #        return (b"", -1, -1)
        #    # end of the patch
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            # minus strand, we have to decipher CIGAR string
            thisstrand = 1
            thisstart = atoi(thisfields[3]) - 1 + sum([cython.cast(cython.int, x) for x in findall(b"(\d+)[MDNX=]", CIGAR)])  #reverse strand should be shifted alen bp
        else:
            thisstrand = 0
            thisstart = atoi(thisfields[3]) - 1

        try:
            thisref = thisref[:thisref.rindex(b".fa")]
        except ValueError:
            pass

        return (thisref, thisstart, thisstrand)


@cython.cclass
class BAMParser(GenericParser):
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
    def __init__(self, filename: str,
                 buffer_size: cython.long = 100000):
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
        f = gzip.open(filename)
        try:
            f.read(10)
        except IOError:
            # not a gzipped file
            self.gzipped = False
        f.close()
        if self.gzipped:
            # open with gzip.open, then wrap it with BufferedReader!
            self.fhd = io.BufferedReader(gzip.open(filename, mode='rb'), buffer_size=READ_BUFFER_SIZE)  # buffersize set to 1M
        else:
            self.fhd = io.open(filename, mode='rb')  # binary mode! I don't expect unicode here!

    @cython.ccall
    def sniff(self):
        """Check the first 3 bytes of BAM file. If it's 'BAM', check
        is success.

        """
        magic_header: bytes
        tsize: cython.int

        magic_header = self.fhd.read(3)
        if magic_header == b"BAM":
            tsize = self.tsize()
            if tsize > 0:
                self.fhd.seek(0)
                return True
            else:
                self.fhd.seek(0)
                raise Exception("File is not of a valid BAM format! %d" % tsize)
        else:
            self.fhd.seek(0)
            return False

    @cython.ccall
    def tsize(self) -> cython.int:
        """Get tag size from BAM file -- read l_seq field.

        Refer to: http://samtools.sourceforge.net/SAM1.pdf

        * This may not work for BAM file from bedToBAM (bedtools),
        since the l_seq field seems to be 0.
        """
        x: cython.int
        header_len: cython.int
        nc: cython.int
        nlength: cython.int
        n: cython.int = 0          # successful read of tag size
        s: cython.double = 0        # sum of tag sizes

        if self.tag_size != -1:
            # if we have already calculated tag size (!= -1),  return it.
            return self.tag_size

        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len = unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = unpack('<i', fread(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = unpack('<i', fread(4))[0]
            # jump over chromosome size, we don't need it
            fread(nlength)
            fseek(ftell() + 4)
        while n < 10:
            entrylength = unpack('<i', fread(4))[0]
            data = fread(entrylength)
            a = unpack('<i', data[16:20])[0]
            s += a
            n += 1
        fseek(0)
        self.tag_size = cython.cast(cython.int, (s/n))
        return self.tag_size

    @cython.ccall
    def get_references(self) -> tuple:
        """
        read in references from BAM header

        return a tuple (references (list of names),
                        rlengths (dict of lengths)
        """
        header_len: cython.int
        x: cython.int
        nc: cython.int
        nlength: cython.int
        refname: bytes
        references: list = []
        rlengths: dict = {}

        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len = unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = unpack('<i', fread(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = unpack('<i', fread(4))[0]
            refname = fread(nlength)[:-1]
            references.append(refname)
            # don't jump over chromosome size
            # we can use it to avoid falling of chrom ends during peak calling
            rlengths[refname] = unpack('<i', fread(4))[0]
        return (references, rlengths)

    @cython.ccall
    def build_fwtrack(self):
        """Build FWTrack from all lines, return a FWTrack object.

        Note only the unique match for a tag is kept.
        """
        i: cython.long = 0          # number of reads kept
        entrylength: cython.int
        fpos: cython.int
        strand: cython.int
        chrid: cython.int
        references: list
        rlengths: dict

        fwtrack = FWTrack(buffer_size=self.buffer_size)
        references, rlengths = self.get_references()  # after this, ptr at list of alignments
        # fseek = self.fhd.seek
        fread = self.fhd.read
        # ftell = self.fhd.tell

        while True:
            try:
                entrylength = unpack("<i", fread(4))[0]
            except struct.error:
                break
            (chrid, fpos, strand) = bam_fw_binary_parse(fread(entrylength))
            if chrid == -1:
                continue
            fwtrack.add_loc(references[chrid], fpos, strand)
            i += 1
            if i % 1000000 == 0:
                info(" %d reads parsed" % i)

        # print(f"{references[chrid]:},{fpos:},{strand:}")
        info("%d reads have been read." % i)
        self.fhd.close()
        fwtrack.set_rlengths(rlengths)
        return fwtrack

    @cython.ccall
    def append_fwtrack(self, fwtrack):
        """Build FWTrack from all lines, return a FWTrack object.

        Note only the unique match for a tag is kept.
        """
        i: cython.long = 0          # number of reads kept
        entrylength: cython.int
        fpos: cython.int
        strand: cython.int
        chrid: cython.int
        references: list
        rlengths: dict

        references, rlengths = self.get_references()
        # fseek = self.fhd.seek
        fread = self.fhd.read
        # ftell = self.fhd.tell

        while True:
            try:
                entrylength = unpack('<i', fread(4))[0]
            except struct.error:
                break
            (chrid, fpos, strand) = bam_fw_binary_parse(fread(entrylength))
            if chrid == -1:
                continue
            fwtrack.add_loc(references[chrid], fpos, strand)
            i += 1
            if i % 1000000 == 0:
                info(" %d reads parsed" % i)

        info("%d reads have been read." % i)
        self.fhd.close()
        # fwtrack.finalize()
        # this is the problematic part. If fwtrack is finalized, then it's impossible to increase the length of it in a step of buffer_size for multiple input files.
        fwtrack.set_rlengths(rlengths)
        return fwtrack


@cython.cclass
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
    2048    supplementary alignment
    """
    # total number of fragments
    n = cython.declare(cython.int, visibility='public')
    # the average length of fragments
    d = cython.declare(cython.float, visibility='public')

    @cython.ccall
    def build_petrack(self):
        """Build PETrackI from all lines, return a FWTrack object.
        """
        i: cython.long = 0          # number of fragments kept
        m: cython.long = 0          # sum of fragment lengths
        entrylength: cython.int
        fpos: cython.int
        chrid: cython.int
        tlen: cython.int
        references: list
        rlengths: dict

        petrack = PETrackI(buffer_size=self.buffer_size)
        references, rlengths = self.get_references()
        # fseek = self.fhd.seek
        fread = self.fhd.read
        # ftell = self.fhd.tell

        # for convenience, only count valid pairs
        add_loc = petrack.add_loc
        err = struct.error

        while True:
            try:
                entrylength = unpack('<i', fread(4))[0]
            except err:
                # e1 += 1
                break
            (chrid, fpos, tlen) = bampe_pe_binary_parse(fread(entrylength))
            if chrid == -1:
                # e2 += 1
                continue
            add_loc(references[chrid], fpos, fpos + tlen)
            m += tlen
            i += 1
            if i % 1000000 == 0:
                info(" %d fragments parsed" % i)

        # print(f"{references[chrid]:},{fpos:},{tlen:}")
        info("%d fragments have been read." % i)
        # debug(f" {e1} Can't identify the length of entry, it may be the end of file, stop looping...")
        # debug(f" {e2} Chromosome name can't be found which means this entry is skipped ...")
        # assert i > 0, "Something went wrong, no fragment has been read! Check input file!"
        self.d = m / i
        self.n = i
        # assert self.d >= 0, "Something went wrong (mean fragment size was negative: %d = %d / %d)" % (self.d, m, i)
        self.fhd.close()
        petrack.set_rlengths(rlengths)
        return petrack

    @cython.ccall
    def append_petrack(self, petrack):
        """Build PETrackI from all lines, return a PETrackI object.
        """
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths
        entrylength: cython.int
        fpos: cython.int
        chrid: cython.int
        tlen: cython.int
        references: list
        rlengths: dict

        references, rlengths = self.get_references()
        # fseek = self.fhd.seek
        fread = self.fhd.read
        # ftell = self.fhd.tell

        # for convenience, only count valid pairs
        add_loc = petrack.add_loc
        err = struct.error

        while True:
            try:
                entrylength = unpack('<i', fread(4))[0]
            except err:
                break
            (chrid, fpos, tlen) = bampe_pe_binary_parse(fread(entrylength))
            if chrid == -1:
                continue
            add_loc(references[chrid], fpos, fpos + tlen)
            m += tlen
            i += 1
            if i % 1000000 == 0:
                info(" %d fragments parsed" % i)

        info("%d fragments have been read." % i)
        self.d = (self.d * self.n + m) / (self.n + i)
        self.n += i
        # assert self.d >= 0, "Something went wrong (mean fragment size was negative: %d = %d / %d)" % (self.d, m, i)
        self.fhd.close()
        # this is the problematic part. If fwtrack is finalized, then it's impossible to increase the length of it in a step of buffer_size for multiple input files.
        # petrack.finalize()
        petrack.set_rlengths(rlengths)
        return petrack


@cython.cclass
class BowtieParser(GenericParser):
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Parse tag sequence, then tag length.

        """
        thisfields: list

        thisline = thisline.rstrip()
        if not thisline:
            return (b"", -1, -1)
        if thisline[0] == b"#":
            return (b"", -1, -1)  # comment line is skipped
        thisfields = thisline.split(b'\t')  # I hope it will never bring me more trouble
        return len(thisfields[4])

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
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
        thisfields: list
        chromname: bytes

        thisline = thisline.rstrip()
        if not thisline:
            return (b"", -1, -1)
        if thisline[0] == b"#":
            return (b"", -1, -1)  # comment line is skipped
        # I hope it will never bring me more trouble
        thisfields = thisline.split(b'\t')

        chromname = thisfields[2]
        try:
            chromname = chromname[:chromname.rindex(b".fa")]
        except ValueError:
            pass

            if thisfields[1] == b"+":
                return (chromname,
                        atoi(thisfields[3]),
                        0)
            elif thisfields[1] == b"-":
                return (chromname,
                        atoi(thisfields[3]) + len(thisfields[4]),
                        1)
            else:
                raise StrandFormatError(thisline, thisfields[1])


@cython.cclass
class FragParser(GenericParser):
    """Parser for Fragment file containing scATAC-seq information.

    Format:

    chromosome frag_leftend frag_rightend barcode count

    Note: Only the first five columns are used!

    """
    n = cython.declare(cython.int, visibility='public')
    d = cython.declare(cython.float, visibility='public')

    @cython.cfunc
    def skip_first_commentlines(self):
        """BEDPEParser needs to skip the first several comment lines.
        """
        l_line: cython.int
        thisline: bytes

        for thisline in self.fhd:
            l_line = len(thisline)
            if thisline and (thisline[:5] != b"track") \
               and (thisline[:7] != b"browser") \
               and (thisline[0] != 35):  # 35 is b"#"
                break

        # rewind from SEEK_CUR
        self.fhd.seek(-l_line, 1)
        return

    @cython.cfunc
    def pe_parse_line(self, thisline: bytes):
        """Parse each line, and return chromosome, left and right
        positions, barcode and count.

        """
        thisfields: list

        thisline = thisline.rstrip()

        # still only support tabular as delimiter.
        thisfields = thisline.split(b'\t')
        try:
            return (thisfields[0],
                    atoi(thisfields[1]),
                    atoi(thisfields[2]),
                    thisfields[3],
                    atoi(thisfields[4]))
        except IndexError:
            raise Exception("Less than 5 columns found at this line: %s\n" % thisline)

    @cython.ccall
    def build_petrack2(self):
        """Build PETrackII from all lines.

        """
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        barcode: bytes
        count: cython.uchar
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths
        tmp: bytes = b""

        petrack = PETrackII(buffer_size=self.buffer_size)
        add_loc = petrack.add_loc

        while True:
            # for each block of input
            tmp += self.fhd.read(READ_BUFFER_SIZE)
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[-1]
            for thisline in lines[:-1]:
                (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
                if left_pos < 0 or not chromosome:
                    continue
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                m += right_pos - left_pos
                i += 1
                if i % 1000000 == 0:
                    info(" %d fragments parsed" % i)
                add_loc(chromosome, left_pos, right_pos, barcode, count)
        # last one
        if tmp:
            (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                add_loc(chromosome, left_pos, right_pos, barcode, count)

        self.d = cython.cast(cython.float, m) / i
        self.n = i
        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"

        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack

    @cython.ccall
    def append_petrack(self, petrack):
        """Build PETrackI from all lines, return a PETrackI object.
        """
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        barcode: bytes
        count: cython.uchar
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths
        tmp: bytes = b""

        add_loc = petrack.add_loc
        while True:
            # for each block of input
            tmp += self.fhd.read(READ_BUFFER_SIZE)
            if not tmp:
                break
            lines = tmp.split(b"\n")
            tmp = lines[-1]
            for thisline in lines[:-1]:
                (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
                if left_pos < 0 or not chromosome:
                    continue
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                m += right_pos - left_pos
                i += 1
                if i % 1000000 == 0:
                    info(" %d fragments parsed" % i)
                add_loc(chromosome, left_pos, right_pos, barcode, count)
        # last one
        if tmp:
            (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                add_loc(chromosome, left_pos, right_pos, barcode, count)

        self.d = (self.d * self.n + m) / (self.n + i)
        self.n += i

        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"
        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack
