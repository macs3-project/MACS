# cython: language_level=3
# cython: profile=True
# cython: linetrace=True
# Time-stamp: <2025-11-12 22:14:22 Tao Liu>

"""Input parsers used across MACS3 for reading alignment-like formats.

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
warn = logger.warn
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
    """Return the first parser that recognises ``fname`` as a supported format.

    Args:
        fname: Path to inspect.
        buffer_size: Buffer size handed to candidate parser constructors.

    Returns:
        Parser instance configured for the detected format.

    Raises:
        Exception: If no parser can identify the file structure.
    """
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
    """Exception raised when strand annotations cannot be interpreted."""
    def __init__(self, string, strand):
        """Capture the offending line and strand token."""
        self.strand = strand
        self.string = string

    def __str__(self):
        """Return a descriptive error message for the malformed strand."""
        return repr("Strand information can not be recognized in this line: \"%s\",\"%s\"" %
                    (self.string, self.strand))


@cython.cclass
class GenericParser:
    """Base parser with helpers for streaming alignment-like text files.

    Attributes:
        filename: Path to the input file.
        gzipped: Whether the input stream is gzipped. (bool)
        tag_size: tag size.
        fhd: Open file handle for the input stream.
        buffer_size: Buffer size for streaming reads.
    """
    filename: str
    gzipped: bool
    tag_size: cython.int
    fhd: object
    buffer_size: cython.long

    def __init__(self, filename: str, buffer_size: cython.long = 100000):
        """Prepare the parser and open the target file.

        Args:
            filename: Path to the input alignment file.
            buffer_size: Chunk size used when reading the stream.
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
        """Advance the stream past any leading comment or header lines."""
        return

    @cython.ccall
    def tsize(self) -> cython.int:
        """Estimate tag length from a sample of valid alignments."""
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
        """Return the inferred tag length for ``thisline`` or ``0`` if invalid."""
        raise NotImplementedError

    @cython.ccall
    def build_fwtrack(self):
        """Create a new ``FWTrack`` populated from the underlying stream."""
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
        """Append parsed locations to an existing ``FWTrack``."""
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
        """Return ``(chromosome, position, strand)`` parsed from ``thisline``."""
        chromosome: bytes = b""
        fpos: cython.int = -1
        strand: cython.int = -1
        return (chromosome, fpos, strand)

    @cython.ccall
    def sniff(self):
        """Return ``True`` when the input appears compatible with this parser."""
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
        """Close the underlying file handle."""
        self.fhd.close()

    @cython.ccall
    def is_gzipped(self) -> bool:
        """Report whether the underlying input stream is gzip compressed."""
        return self.gzipped


@cython.cclass
class BEDParser(GenericParser):
    """Parser for standard BED records with optional strand column."""

    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip ``track``/``browser``/``#`` lines at the top of BED files."""
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
        """Return fragment length encoded in a BED line or ``0`` if invalid."""
        thisline = thisline.rstrip()
        if not thisline:
            return 0

        thisfields = thisline.split(b'\t')
        return atoi(thisfields[2]) - atoi(thisfields[1])

    @cython.cfunc
    def fw_parse_line(self, thisline: bytes) -> tuple:
        """Parse a BED entry into ``(chromosome, position, strand)``.

        Args:
            thisline: Raw line from the BED file.

        Returns:
            Tuple containing the chromosome name, 5' coordinate for the
            strand, and strand flag (0 for ``+``, 1 for ``-``).
        """
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
    """Parser for three-column BEDPE-style fragments (chrom, left, right)."""
    n = cython.declare(cython.int, visibility='public')
    d = cython.declare(cython.float, visibility='public')

    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip ``track``/``browser``/``#`` lines at the top of BEDPE files."""
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
        """Parse a fragment line into ``(chrom, left, right)`` integers."""
        thisfields: list

        thisline = thisline.rstrip()

        # still only support tabular as delimiter.
        thisfields = thisline.split(b'\t')
        try:
            return (thisfields[0],
                    atoi(thisfields[1]),
                    atoi(thisfields[2]))
        except IndexError:
            raise Exception("Less than 3 columns found at this line: %s\n" %
                            thisline)

    @cython.ccall
    def build_petrack(self):
        """Return a ``PETrackI`` constructed from the entire stream.

        Returns:
            PETrackI: Paired-end track populated from the input stream.

        Examples:
            .. code-block:: python

                from MACS3.IO.Parser import BEDPEParser
                parser = BEDPEParser("fragments.bedpe")
                petrack = parser.build_petrack()
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
        """Append fragments from the stream to an existing ``PETrackI``."""
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
    """Parser for single-end ELAND result tables."""

    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip lines beginning with ``#`` before data rows."""
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
        """Return tag length for ELAND single-end entries or ``0`` if invalid."""
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
        """Parse an ELAND result line into location and strand tuple.

        Args:
            thisline: Raw ELAND text line.

        Returns:
            ``(chromosome, position, strand)`` with ``position`` set to
            ``-1`` when the line does not encode a valid alignment.
        """
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
    """Parser for ELAND multi-hit reports (``s_N_eland_multi`` format)."""

    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip lines beginning with ``#`` before data rows."""
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
        """Return tag length for ELAND multi entries or ``0`` if ambiguous."""
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
        """Parse ELAND multi-format line into a single-hit location tuple.

        Args:
            thisline: Raw ELAND multi line.

        Returns:
            ``(chromosome, position, strand)`` or ``(b"", -1, -1)`` if the
            entry is not uniquely mappable.
        """
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
    """Parser for ELAND export tab-delimited files."""
    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip lines beginning with ``#`` before data rows."""
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
        """Return tag length for ELAND export entries or ``0`` if invalid."""
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
        """Parse an ELAND export entry to ``(chromosome, position, strand)``.

        Args:
            thisline: Raw ELAND export line.

        Returns:
            Location tuple when the line contains a successful alignment;
            otherwise ``(b"", -1, -1)``.
        """
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
    """Parser for SAM alignment text files with standard SAM flags."""

    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip SAM header lines beginning with ``@``."""
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
        """Return read length for valid SAM records or ``0`` if filtered."""
        thisfields: list
        bwflag: cython.int

        thisline = thisline.rstrip()
        if not thisline:
            return 0
        thisfields = thisline.split(b'\t')
        bwflag = atoi(thisfields[1])
        if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
            # unmapped sequence or bad sequence or 2nd or sup alignment
            return 0

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
        """Parse a SAM alignment into ``(chromosome, position, strand)``.

        Args:
            thisline: Raw SAM alignment line.

        Returns:
            Tuple identifying the leftmost coordinate for forward strands or
            the rightmost coordinate for reverse strands. Returns
            ``(b"", -1, -1)`` for reads that do not satisfy mapping filters.
        """
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
            thisstart = atoi(thisfields[3]) - 1 + sum([cython.cast(cython.int, x) for x in findall(b"(\\d+)[MDNX=]", CIGAR)])  #reverse strand should be shifted alen bp
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
    """Parser for BAM binary alignment files."""
    def __init__(self, filename: str,
                 buffer_size: cython.long = 100000):
        """Prepare the parser and open the BAM file.

        Args:
            filename: Path to the BAM file.
            buffer_size: Chunk size used when reading BGZF blocks.
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
            self.fhd = io.BufferedReader(gzip.open(filename, mode='rb'),
                                         buffer_size=READ_BUFFER_SIZE)
        else:
            # binary mode! I don't expect unicode here!
            self.fhd = io.open(filename, mode='rb')

    @cython.ccall
    def sniff(self):
        """Return ``True`` if the file begins with the BAM magic string."""
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
                raise Exception("File is not of a valid BAM format! %d" %
                                tsize)
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
        """Return ``(references, lengths)`` extracted from the BAM header."""
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
        """Append uniquely mapped reads to an existing ``FWTrack``."""
        i: cython.long = 0          # number of reads kept
        entrylength: cython.int
        fpos: cython.int
        strand: cython.int
        chrid: cython.int
        references: list
        rlengths: dict

        fwtrack = FWTrack(buffer_size=self.buffer_size)
        # after this, ptr at list of alignments
        references, rlengths = self.get_references()  
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
        """Append uniquely mapped reads to an existing ``FWTrack``."""
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
        # this is the problematic part. If fwtrack is finalized, then
        # it's impossible to increase the length of it in a step of
        # buffer_size for multiple input files.
        fwtrack.set_rlengths(rlengths)
        return fwtrack


@cython.cclass
class BAMPEParser(BAMParser):
    """Parser for paired-end BAM files that yields fragment spans.

    Attributes:
        n: Total number of fragments.
        d: Average fragment length.
        
    """
    # total number of fragments
    n = cython.declare(cython.int, visibility='public')
    # the average length of fragments
    d = cython.declare(cython.float, visibility='public')

    @cython.ccall
    def build_petrack(self):
        """Return a ``PETrackI`` populated with inferred fragments.

        Returns:
            PETrackI: Paired-end track populated from BAM pairs.

        Examples:
            .. code-block:: python

                from MACS3.IO.Parser import BAMPEParser
                parser = BAMPEParser("reads.bam")
                petrack = parser.build_petrack()
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

        info("%d fragments have been read." % i)
        self.d = m / i
        self.n = i
        self.fhd.close()
        petrack.set_rlengths(rlengths)
        return petrack

    @cython.ccall
    def append_petrack(self, petrack):
        """Append inferred fragments to an existing ``PETrackI``.

        Args:
            petrack: Existing paired-end track to append to.

        Returns:
            PETrackI: The updated track instance.

        Examples:
            .. code-block:: python

                from MACS3.IO.Parser import BAMPEParser
                parser = BAMPEParser("reads.bam")
                petrack = parser.build_petrack()
                # Later, append more fragments from another file:
                parser2 = BAMPEParser("more_reads.bam")
                petrack = parser2.append_petrack(petrack)
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
        fread = self.fhd.read

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
        self.fhd.close()
        petrack.set_rlengths(rlengths)
        return petrack


@cython.cclass
class BowtieParser(GenericParser):
    """Parser for Bowtie or Maqview single-end map files."""
    @cython.cfunc
    def tlen_parse_line(self, thisline: bytes) -> cython.int:
        """Return read length for Bowtie map entries or ``0`` if invalid."""
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
    """Parser for scATAC fragment TSV files with barcode counts."""
    n = cython.declare(cython.int, visibility='public')
    d = cython.declare(cython.float, visibility='public')

    @cython.cfunc
    def skip_first_commentlines(self):
        """Skip ``track``/``browser``/``#`` lines at the top of fragment files."""
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
        """Parse a fragment line into ``(chrom, left, right, barcode, count)``."""
        thisfields: list
        thiscount: cython.ushort

        thisline = thisline.rstrip()

        # still only support tabular as delimiter.
        thisfields = thisline.split(b'\t')
        try:
            try:
                thiscount = atoi(thisfields[4])
            except OverflowError:
                thiscount = 65535
                warn(f"The count in this line is over 65535, and will be capped at 65535: {thisline}")
            return (thisfields[0],
                    atoi(thisfields[1]),
                    atoi(thisfields[2]),
                    thisfields[3],
                    thiscount)
        except IndexError:
            raise Exception("Less than 5 columns found at this line: %s\n" %
                            thisline)

    @cython.ccall
    def build_petrack(self, max_count=0):
        """Return a ``PETrackII`` populated with fragments and barcodes.

        Args:
            max_count: Optional cap applied to per-fragment count values.

        Returns:
            PETrackII: Paired-end track with barcode and count metadata.

        Examples:
            .. code-block:: python

                from MACS3.IO.Parser import FragParser
                parser = FragParser("fragments.tsv.gz")
                petrack = parser.build_petrack(max_count=10)
        """
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        barcode: bytes
        count: cython.ushort
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
                if max_count:
                    count = min(count, max_count)
                add_loc(chromosome, left_pos, right_pos, barcode, count)
        # last one
        if tmp:
            (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                if max_count:
                    count = min(count, max_count)                
                add_loc(chromosome, left_pos, right_pos, barcode, count)

        self.d = cython.cast(cython.float, m) / i
        self.n = i
        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"

        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack

    @cython.ccall
    def append_petrack(self, petrack, max_count=0):
        """Append barcode-aware fragments to an existing ``PETrackI``.

        Args:
            petrack: Existing paired-end track to append to.
            max_count: Optional cap applied to per-fragment count values.

        Returns:
            PETrackII: The updated track instance.

        Examples:
            .. code-block:: python

                from MACS3.IO.Parser import FragParser
                parser = FragParser("fragments.tsv.gz")
                petrack = parser.build_petrack()
                parser2 = FragParser("more_fragments.tsv.gz")
                petrack = parser2.append_petrack(petrack, max_count=10)
        """
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        barcode: bytes
        count: cython.ushort
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
                if max_count:
                    count = min(count, max_count)                    
                add_loc(chromosome, left_pos, right_pos, barcode, count)
        # last one
        if tmp:
            (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
            if left_pos >= 0 and chromosome:
                assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
                i += 1
                m += right_pos - left_pos
                if max_count:
                    count = min(count, max_count)                
                add_loc(chromosome, left_pos, right_pos, barcode, count)

        self.d = (self.d * self.n + m) / (self.n + i)
        self.n += i

        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"
        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack
