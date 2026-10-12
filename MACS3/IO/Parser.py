# cython: language_level=3
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
import os
import sys

import cython
from cython.cimports.cpython import bool

from cython.cimports.libc.stdlib import atoi, malloc, calloc, realloc, free
from cython.cimports.libc.string import memcpy, memmove, memset, memchr
from cython.cimports.cpython.bytes import (PyBytes_FromStringAndSize,
                                         PyBytes_AS_STRING)
from cython.cimports.cpython.exc import PyErr_CheckSignals
from cython.cimports.cpython.bytearray import PyByteArray_AS_STRING
from cython.cimports.MACS3.IO.czlib import (z_stream, inflateInit2, inflate,
                                           inflateReset, inflateEnd, Z_OK,
                                           Z_STREAM_END, Z_BUF_ERROR,
                                           Z_NO_FLUSH)
from cython.cimports.MACS3.IO.cbgzf_mt import (bgzf_job_t, bgzf_pool_t,
                                              bgzf_split, bgzf_pool_new,
                                              bgzf_pool_submit,
                                              bgzf_pool_wait,
                                              bgzf_pool_free, bgzf_pool_run,
                                              bgzf_task_fn)
import numpy as np

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


# ------------------------------------
# BAM alignment records, read in large decompressed blocks
# ------------------------------------

# Compressed BAM input is read _BAM_IN_SIZE bytes at a time and inflated
# into a buffer of _BAM_OUT_SIZE bytes (grown if one record is larger),
# whose records are then walked in C.
_BAM_IN_SIZE = cython.declare(cython.Py_ssize_t, 1048576)
_BAM_OUT_SIZE = cython.declare(cython.Py_ssize_t, 4194304)
# bytes of an alignment record before read_name, after block_size
_BAM_CORE_SIZE = cython.declare(cython.Py_ssize_t, 32)

# BGZF input is inflated on min(_BAM_MAX_THREADS, cores) threads when
# more than one core is available: _BAM_PAR_IN compressed bytes at a time
# are split into blocks, which are inflated into _BAM_PAR_OUT bytes
# following _BAM_PAR_HEAD bytes kept for the unconsumed tail of the
# previous batch.
_BAM_MAX_THREADS = cython.declare(cython.int, 8)
_BAM_PAR_IN = cython.declare(cython.Py_ssize_t, 4194304)
_BAM_PAR_OUT = cython.declare(cython.Py_ssize_t, 16777216)
_BAM_PAR_HEAD = cython.declare(cython.Py_ssize_t, 262144)
_BAM_PAR_MAXJOBS = cython.declare(cython.Py_ssize_t, 8192)

# one of the two batches of the parallel path
_BAMBatch = cython.struct(
    buf=cython.p_uchar,              # decompressed output, from offset _BAM_PAR_HEAD
    cap=cython.Py_ssize_t,           # bytes allocated at buf
    inb=cython.p_uchar,              # compressed input
    in_len=cython.Py_ssize_t,        # bytes read into inb
    jobs=cython.pointer(bgzf_job_t),  # one per block
    njobs=cython.Py_ssize_t,
    out_len=cython.Py_ssize_t)       # bytes the blocks inflate to


def _inflate_threads() -> int:
    """Return how many threads inflate a BGZF file (BAM, or a BGZF
    fragment file): min(8, the cores this process may run on)."""
    if hasattr(os, "sched_getaffinity"):
        return min(_BAM_MAX_THREADS, len(os.sched_getaffinity(0)))
    return min(_BAM_MAX_THREADS, os.cpu_count() or 1)


def _starts_with_bgzf(filename: str) -> bool:
    """Return whether ``filename`` begins with a BGZF block header: a
    deflate gzip member whose extra field starts with the "BC" subfield.

    This only decides whether to start inflate threads for a text file;
    ``_BAMStream`` checks every block itself and reads anything that is
    not well-formed BGZF serially.
    """
    with io.open(filename, mode='rb') as f:
        h = f.read(16)
    return (len(h) == 16 and h[0] == 31 and h[1] == 139 and h[2] == 8 and
            (h[3] & 4) != 0 and h[12] == 66 and h[13] == 67 and
            h[14] == 2 and h[15] == 0)


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
def _le_uint32(p: cython.p_uchar) -> cython.uint:
    """Little-endian uint32 at ``p`` (BAM fields are little-endian)."""
    return (cython.cast(cython.uint, p[0]) |
            (cython.cast(cython.uint, p[1]) << 8) |
            (cython.cast(cython.uint, p[2]) << 16) |
            (cython.cast(cython.uint, p[3]) << 24))


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
def _le_int32(p: cython.p_uchar) -> cython.int:
    """Little-endian int32 at ``p``."""
    return cython.cast(cython.int, _le_uint32(p))


@cython.cclass
class _BAMStream:
    """The decompressed bytes of a BAM file, produced in large blocks.

    Compressed input is read ``_BAM_IN_SIZE`` bytes at a time and
    inflated by zlib straight into ``buf``. zlib parses the gzip framing
    of every BGZF block (or gzip member) and checks its CRC32 and
    length, as Python's gzip module does, and zero bytes between members
    are skipped, as ``gzip.GzipFile`` does. A file that is not gzipped
    is copied through unchanged.

    When more than one core is available, a gzipped file is read by the
    parallel path first: ``_BAM_PAR_IN`` compressed bytes at a time are
    split into BGZF blocks by their headers, and the blocks are inflated
    by libdeflate (``bgzf_mt.c``), each checked against its CRC32 and
    ISIZE, on ``_inflate_threads()`` threads, the next batch
    while the caller walks the records of the current one, and the
    compressed input of the batch after it is read while the threads
    inflate. From the
    first block that is not a well-formed BGZF block (another kind of
    gzip member, gzip flags other than FEXTRA alone, zero padding, a
    block that is corrupt or truncated, or that does not inflate to
    exactly its ISIZE), the rest of the file is
    read by the serial path, so the bytes delivered and the errors raised
    are the serial path's. With ``parallel`` False, only the serial path
    is used.

    With ``window`` set as well (the FRAG reader), ``refill`` leaves in
    ``buf[0:cap]`` exactly the bytes the serial path would leave there:
    the unconsumed bytes and what follows them, ``cap`` bytes in all
    (``_BAM_OUT_SIZE``, or more after a ``need`` larger than that). When
    those bytes lie in one batch, ``buf`` points into it; otherwise they
    are copied into a buffer of the stream's own, and from the first
    block the parallel path does not take, the serial path fills that
    buffer on. So every ``refill`` gives the same bytes, and raises the
    same error, as with ``parallel`` False.

    ``buf[start:end]`` holds the bytes not consumed yet: callers walk
    records there, advance ``start`` and call ``refill`` for more.
    """
    fh: object
    inbuf: bytearray
    inptr: cython.p_uchar
    strm: z_stream
    gzipped: cython.bint
    zready: cython.bint
    in_member: cython.bint
    eof: cython.bint
    buf: cython.p_uchar
    cap: cython.Py_ssize_t
    start: cython.Py_ssize_t
    end: cython.Py_ssize_t
    # the parallel path; buf is batches[cur].buf while par is set
    par: cython.bint
    pool: cython.pointer(bgzf_pool_t)
    batches: cython.pointer(_BAMBatch)
    cur: cython.int
    inflight: cython.bint       # batches[1 - cur] is being inflated
    par_eof: cython.bint        # the parallel path has read the whole file
    inbufs: list                # the bytearrays behind batches[i].inb
    inviews: list
    carry: cython.p_uchar       # compressed bytes the last batch did not take
    carry_len: cython.Py_ssize_t
    staged: cython.bint         # par_refill has read the next batch's input
    joined: bytearray           # the serial path's first input, when two parts
    # window mode: buf is pbuf, or points into batches[cur]
    window: cython.bint
    pbuf: cython.p_uchar        # the stream's own buffer
    pcap: cython.Py_ssize_t
    src_lo: cython.p_uchar      # the bytes of batches[cur]: src_lo[0:],
    src_pos: cython.p_uchar     # those after buf[0:end] from src_pos,
    src_end: cython.p_uchar     # up to src_end
    win_stop: cython.bint       # after src_end, the serial path reads on
    stop_src: cython.p_uchar    # from stop_src[0:stop_n], then
    stop_n: cython.Py_ssize_t   # stop_src2[0:stop_n2]
    stop_src2: cython.p_uchar
    stop_n2: cython.Py_ssize_t

    def __cinit__(self):
        self.buf = cython.NULL
        self.zready = False
        self.par = False
        self.pool = cython.NULL
        self.batches = cython.NULL
        self.inflight = False
        self.window = False

    def __init__(self, filename: str, gzipped: cython.bint,
                 parallel: cython.bint = True, window: cython.bint = False):
        i: cython.int
        nthreads: cython.int

        memset(cython.address(self.strm), 0, cython.sizeof(z_stream))
        self.gzipped = gzipped
        self.in_member = False
        self.eof = False
        self.start = 0
        self.end = 0
        self.inbuf = bytearray(_BAM_IN_SIZE)
        self.inptr = cython.cast(cython.p_uchar,
                                 PyByteArray_AS_STRING(self.inbuf))
        if gzipped:
            # windowBits 15 + 16: 32 kB window, gzip framing
            if inflateInit2(cython.address(self.strm), 15 + 16) != Z_OK:
                raise MemoryError("zlib could not be initialized")
            self.zready = True
            nthreads = _inflate_threads() if parallel else 1
            if nthreads > 1:
                self.pool = bgzf_pool_new(nthreads)
        if self.pool != cython.NULL:
            self.batches = cython.cast(cython.pointer(_BAMBatch),
                                       calloc(2, cython.sizeof(_BAMBatch)))
            if self.batches == cython.NULL:
                raise MemoryError()
            self.inbufs = []
            self.inviews = []
            for i in range(2):
                self.batches[i].cap = _BAM_PAR_HEAD + _BAM_PAR_OUT
                self.batches[i].buf = cython.cast(cython.p_uchar,
                                                  malloc(self.batches[i].cap))
                self.batches[i].jobs = cython.cast(
                    cython.pointer(bgzf_job_t),
                    malloc(_BAM_PAR_MAXJOBS * cython.sizeof(bgzf_job_t)))
                if (self.batches[i].buf == cython.NULL or
                        self.batches[i].jobs == cython.NULL):
                    raise MemoryError()
                self.inbufs.append(bytearray(_BAM_PAR_IN))
                self.inviews.append(memoryview(self.inbufs[i]))
                self.batches[i].inb = cython.cast(
                    cython.p_uchar, PyByteArray_AS_STRING(self.inbufs[i]))
            self.par = True
            self.par_eof = False
            self.cur = 0
            self.carry = self.batches[1].inb
            self.carry_len = 0
            self.staged = False
            if window:
                self.window = True
                self.src_lo = cython.NULL
                self.src_pos = cython.NULL
                self.src_end = cython.NULL
                self.win_stop = False
                self.cap = _BAM_OUT_SIZE
                self.pcap = _BAM_OUT_SIZE
                self.pbuf = cython.cast(cython.p_uchar, malloc(self.pcap))
                if self.pbuf == cython.NULL:
                    raise MemoryError()
                self.buf = self.pbuf
            else:
                self.buf = self.batches[0].buf
                self.cap = self.batches[0].cap
        else:
            self.cap = _BAM_OUT_SIZE
            self.buf = cython.cast(cython.p_uchar, malloc(self.cap))
            if self.buf == cython.NULL:
                raise MemoryError()
        self.fh = io.open(filename, mode='rb', buffering=0)

    def __dealloc__(self):
        i: cython.int

        self.stop_threads()
        if self.zready:
            inflateEnd(cython.address(self.strm))
        if self.batches != cython.NULL:
            # buf is one of the batches' while par is set, or pbuf
            for i in range(2):
                free(self.batches[i].buf)
                free(self.batches[i].jobs)
            free(self.batches)
        if self.par and self.window:
            free(self.pbuf)
        if not self.par and self.buf != cython.NULL:
            free(self.buf)

    @cython.cfunc
    def stop_threads(self):
        """Finish the batch being inflated, if any, and join the threads."""
        if self.pool != cython.NULL:
            if self.inflight:
                with cython.nogil:
                    bgzf_pool_wait(self.pool)
                self.inflight = False
            bgzf_pool_free(self.pool)
            self.pool = cython.NULL

    @cython.ccall
    def close(self):
        """Join the inflating threads and close the underlying file."""
        self.stop_threads()
        if self.fh is not None:
            self.fh.close()
            self.fh = None

    @cython.cfunc
    def par_read(self, t: cython.int):
        """Put the compressed bytes the last batch did not take at the
        start of ``batches[t].inb`` and read the file after them, until
        it holds ``_BAM_PAR_IN`` bytes or the file ends."""
        b: cython.pointer(_BAMBatch) = self.batches + t
        view: object = self.inviews[t]
        ilen: cython.Py_ssize_t = self.carry_len
        n: cython.Py_ssize_t

        if ilen > 0 and self.carry != b.inb:
            memmove(b.inb, self.carry, ilen)
        while ilen < _BAM_PAR_IN and not self.par_eof:
            n = self.fh.readinto(view[ilen:])
            if n == 0:
                self.par_eof = True
            ilen += n
        b.in_len = ilen

    @cython.cfunc
    def par_submit(self) -> cython.bint:
        """Read compressed input after the bytes the last batch did not
        take, unless ``par_refill`` has read it already, split it into
        BGZF blocks and start inflating them into ``batches[1 - cur]``.
        Returns False, starting nothing, when the input does not begin
        with a complete BGZF block."""
        b: cython.pointer(_BAMBatch) = self.batches + (1 - self.cur)
        ilen: cython.Py_ssize_t
        in_used: cython.size_t
        out_used: cython.size_t
        reason: cython.int

        if not self.staged:
            self.par_read(1 - self.cur)
        self.staged = False
        ilen = b.in_len
        b.njobs = bgzf_split(b.inb, ilen, b.buf + _BAM_PAR_HEAD,
                             b.cap - _BAM_PAR_HEAD, b.jobs, _BAM_PAR_MAXJOBS,
                             cython.address(in_used), cython.address(out_used),
                             cython.address(reason))
        b.out_len = out_used
        self.carry = b.inb + in_used
        self.carry_len = ilen - in_used
        if b.njobs == 0:
            return False
        bgzf_pool_submit(self.pool, b.jobs, b.njobs)
        self.inflight = True
        return True

    @cython.cfunc
    def to_serial(self, src: cython.p_uchar, n: cython.Py_ssize_t,
                  src2: cython.p_uchar, n2: cython.Py_ssize_t):
        """Leave the parallel path: the serial path inflates ``src[0:n]``,
        then ``src2[0:n2]``, then the rest of the file, into the current
        buffer."""
        p: cython.p_uchar

        self.stop_threads()
        self.par = False
        if self.window:
            # buf is pbuf, the serial path's from now on
            free(self.batches[self.cur].buf)
            self.pbuf = cython.NULL
        # buf, the current batch's, is the serial path's from now on
        self.batches[self.cur].buf = cython.NULL
        free(self.batches[1 - self.cur].buf)
        self.batches[1 - self.cur].buf = cython.NULL
        if n2 > 0:
            # input read ahead of src's batch: join the two parts
            self.joined = bytearray(n + n2)
            p = cython.cast(cython.p_uchar, PyByteArray_AS_STRING(self.joined))
            memcpy(p, src, n)
            memcpy(p + n, src2, n2)
            src = p
            n += n2
        self.strm.next_in = src
        self.strm.avail_in = n
        self.in_member = False

    @cython.cfunc
    def par_refill(self, need: cython.Py_ssize_t) -> cython.bint:
        """``refill`` on the parallel path. Returns False when it has
        handed the rest of the file to the serial path."""
        b: cython.pointer(_BAMBatch)
        c: cython.pointer(_BAMBatch)
        nj: cython.Py_ssize_t
        g: cython.Py_ssize_t
        d: cython.Py_ssize_t
        L: cython.Py_ssize_t
        newcap: cython.Py_ssize_t
        newbuf: cython.p_uchar
        s: cython.pointer(_BAMBatch)
        cl: cython.Py_ssize_t = 0
        src2: cython.p_uchar
        n2: cython.Py_ssize_t

        if self.end - self.start >= need:
            return True
        while True:
            if not self.inflight and not self.par_submit():
                self.to_serial(self.carry, self.carry_len, cython.NULL, 0)
                return False
            b = self.batches + (1 - self.cur)
            s = self.batches + self.cur
            L = self.end - self.start
            if L <= _BAM_PAR_HEAD:
                # b becomes the current batch below, so the next batch is
                # s, whose blocks are all inflated: read its input now,
                # while the threads inflate b
                cl = self.carry_len
                self.par_read(self.cur)
                self.staged = True
            with cython.nogil:
                bgzf_pool_wait(self.pool)
            self.inflight = False
            # the blocks before the first one the threads did not accept
            nj = b.njobs
            g = 0
            while g < nj and b.jobs[g].ok:
                g += 1
            if g < nj:
                d = b.jobs[g].dst - (b.buf + _BAM_PAR_HEAD)
            else:
                d = b.out_len
            if L <= _BAM_PAR_HEAD:
                # move the unconsumed bytes in front of the new batch
                memcpy(b.buf + _BAM_PAR_HEAD - L, self.buf + self.start, L)
                self.cur = 1 - self.cur
                self.buf = b.buf
                self.cap = b.cap
                self.start = _BAM_PAR_HEAD - L
                self.end = _BAM_PAR_HEAD + d
            else:
                # a long record: append the new batch to the current one
                c = self.batches + self.cur
                memmove(self.buf, self.buf + self.start, L)
                if L + d > c.cap:
                    newcap = c.cap
                    while newcap < L + d:
                        newcap *= 2
                    newbuf = cython.cast(cython.p_uchar,
                                         realloc(c.buf, newcap))
                    if newbuf == cython.NULL:
                        raise MemoryError()
                    c.buf = newbuf
                    c.cap = newcap
                    self.buf = newbuf
                    self.cap = newcap
                memcpy(self.buf + L, b.buf + _BAM_PAR_HEAD, d)
                self.start = 0
                self.end = L + d
            if g < nj:
                # read on from the rejected block serially, then from the
                # input read ahead into s past its first cl bytes, which
                # are b's input after its last block
                src2 = cython.NULL
                n2 = 0
                if self.staged:
                    self.staged = False
                    src2 = s.inb + cl
                    n2 = s.in_len - cl
                self.to_serial(cython.cast(cython.p_uchar, b.jobs[g].src),
                               b.in_len -
                               (cython.cast(cython.p_uchar, b.jobs[g].src) -
                                b.inb),
                               src2, n2)
                return False
            # inflate the next batch while the caller walks this one
            if not self.par_submit():
                self.to_serial(self.carry, self.carry_len, cython.NULL, 0)
                return False
            if self.end - self.start >= need:
                return True

    @cython.cfunc
    def win_refill(self, need: cython.Py_ssize_t) -> cython.bint:
        """``refill`` in window mode. Returns False, with the unconsumed
        bytes at ``buf[0]`` (``buf`` is ``pbuf``), ``cap`` at least
        ``need``, and ``buf`` filled as far as the parallel path went,
        when it has handed the rest of the file to the serial path."""
        L: cython.Py_ssize_t = self.end - self.start
        r: cython.p_uchar
        newcap: cython.Py_ssize_t
        newbuf: cython.p_uchar

        # the serial path's buffer size
        if need > self.cap:
            newcap = self.cap
            while newcap < need:
                newcap *= 2
            self.cap = newcap
        # the bytes after the unconsumed ones are at src_pos, in
        # batches[cur]; if the unconsumed ones are there too, and the
        # cap bytes from them, buf points at them
        if self.src_lo != cython.NULL and self.src_pos - self.src_lo >= L:
            r = self.src_pos - L
            if self.src_end - r >= self.cap:
                self.buf = r
                self.start = 0
                self.end = self.cap
                self.src_pos = r + self.cap
                return True
        # otherwise the unconsumed bytes go to pbuf and the rest is copied
        if self.buf == self.pbuf:
            if self.start > 0:
                memmove(self.pbuf, self.pbuf + self.start, L)
            if self.cap > self.pcap:
                newbuf = cython.cast(cython.p_uchar,
                                     realloc(self.pbuf, self.cap))
                if newbuf == cython.NULL:
                    raise MemoryError()
                self.pbuf = newbuf
                self.pcap = self.cap
        else:
            if self.cap > self.pcap:
                newbuf = cython.cast(cython.p_uchar,
                                     realloc(self.pbuf, self.cap))
                if newbuf == cython.NULL:
                    raise MemoryError()
                self.pbuf = newbuf
                self.pcap = self.cap
            memcpy(self.pbuf, self.buf + self.start, L)
        self.buf = self.pbuf
        self.start = 0
        self.end = L
        return self.win_fill()

    @cython.cfunc
    def win_fill(self) -> cython.bint:
        """Copy the parallel path's bytes into ``buf``, which is
        ``pbuf``, until it holds ``cap`` bytes. Returns False, ``buf``
        filled as far as the parallel path went, when it has handed the
        rest of the file to the serial path."""
        b: cython.pointer(_BAMBatch)
        s: cython.pointer(_BAMBatch)
        k: cython.Py_ssize_t
        nj: cython.Py_ssize_t
        g: cython.Py_ssize_t
        d: cython.Py_ssize_t
        cl: cython.Py_ssize_t

        while self.end < self.cap:
            k = self.src_end - self.src_pos
            if k > 0:
                if k > self.cap - self.end:
                    k = self.cap - self.end
                memcpy(self.buf + self.end, self.src_pos, k)
                self.src_pos += k
                self.end += k
                continue
            if self.win_stop:
                self.to_serial(self.stop_src, self.stop_n,
                               self.stop_src2, self.stop_n2)
                return False
            if not self.inflight and not self.par_submit():
                self.to_serial(self.carry, self.carry_len, cython.NULL, 0)
                return False
            b = self.batches + (1 - self.cur)
            s = self.batches + self.cur
            # every byte of s is copied: read the input of the batch
            # after b into it while the threads inflate b
            cl = self.carry_len
            self.par_read(self.cur)
            self.staged = True
            with cython.nogil:
                bgzf_pool_wait(self.pool)
            self.inflight = False
            # the blocks before the first one the threads did not accept
            nj = b.njobs
            g = 0
            while g < nj and b.jobs[g].ok:
                g += 1
            if g < nj:
                d = b.jobs[g].dst - (b.buf + _BAM_PAR_HEAD)
            else:
                d = b.out_len
            self.cur = 1 - self.cur
            self.src_lo = b.buf + _BAM_PAR_HEAD
            self.src_pos = self.src_lo
            self.src_end = self.src_lo + d
            if g < nj:
                # after b's bytes, read on serially from the rejected
                # block, then from the input read ahead into s past its
                # first cl bytes, which are b's input after its last block
                self.win_stop = True
                self.staged = False
                self.stop_src = cython.cast(cython.p_uchar, b.jobs[g].src)
                self.stop_n = b.in_len - (self.stop_src - b.inb)
                self.stop_src2 = s.inb + cl
                self.stop_n2 = s.in_len - cl
            elif not self.par_submit():
                # the input after b's blocks is not a complete BGZF block
                self.win_stop = True
                self.stop_src = self.carry
                self.stop_n = self.carry_len
                self.stop_src2 = cython.NULL
                self.stop_n2 = 0
        return True

    @cython.cfunc
    def refill(self, need: cython.Py_ssize_t) -> cython.bint:
        """Keep the unconsumed bytes, then read and decompress until the
        buffer is full or the file ends.

        Returns whether at least ``need`` bytes are available.
        """
        n: cython.Py_ssize_t
        ret: cython.int
        newcap: cython.Py_ssize_t
        newbuf: cython.p_uchar
        msg: bytes

        if self.par:
            if self.window:
                if self.win_refill(need):
                    return self.end - self.start >= need
            elif self.par_refill(need):
                return True
        if self.start > 0:
            memmove(self.buf, self.buf + self.start, self.end - self.start)
            self.end -= self.start
            self.start = 0
        if need > self.cap:
            newcap = self.cap
            while newcap < need:
                newcap *= 2
            newbuf = cython.cast(cython.p_uchar, realloc(self.buf, newcap))
            if newbuf == cython.NULL:
                raise MemoryError()
            self.buf = newbuf
            self.cap = newcap

        while self.end < self.cap and not self.eof:
            if self.strm.avail_in == 0:
                n = self.fh.readinto(self.inbuf)
                if n == 0:
                    if self.in_member:
                        raise EOFError("Compressed file ended before the "
                                       "end-of-stream marker was reached")
                    self.eof = True
                    break
                self.strm.next_in = self.inptr
                self.strm.avail_in = n
            if not self.gzipped:
                n = min(cython.cast(cython.Py_ssize_t, self.strm.avail_in),
                        self.cap - self.end)
                memcpy(self.buf + self.end, self.strm.next_in, n)
                self.strm.next_in += n
                self.strm.avail_in -= n
                self.end += n
                continue
            if not self.in_member:
                # zero padding between gzip members is skipped
                while self.strm.avail_in > 0 and self.strm.next_in[0] == 0:
                    self.strm.next_in += 1
                    self.strm.avail_in -= 1
                if self.strm.avail_in == 0:
                    continue
                inflateReset(cython.address(self.strm))
                self.in_member = True
            self.strm.next_out = self.buf + self.end
            self.strm.avail_out = self.cap - self.end
            ret = inflate(cython.address(self.strm), Z_NO_FLUSH)
            self.end = self.cap - self.strm.avail_out
            if ret == Z_STREAM_END:
                self.in_member = False
            elif ret != Z_OK and (ret != Z_BUF_ERROR or
                                  (self.strm.avail_in > 0 and
                                   self.strm.avail_out > 0)):
                msg = b""
                if self.strm.msg != cython.NULL:
                    msg = self.strm.msg
                raise gzip.BadGzipFile("Error %d while decompressing data: %s"
                                       % (ret, msg.decode(errors="replace")))
        return self.end - self.start >= need

    @cython.cfunc
    def skip(self, n: cython.Py_ssize_t):
        """Discard the next ``n`` decompressed bytes."""
        k: cython.Py_ssize_t

        while n > 0:
            if self.start == self.end and not self.refill(1):
                raise EOFError("BAM file ended inside its header")
            k = min(n, self.end - self.start)
            self.start += k
            n -= k


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def _bam_load_fwtrack(stream: _BAMStream, references: list,
                      fwtrack) -> cython.long:
    """Add every kept alignment's 5' end in ``stream`` to ``fwtrack``.

    The reads kept, their positions and their order are those of parsing
    each record and calling ``fwtrack.add_loc(references[refID], pos,
    strand)``: unmapped, QC-failed, secondary and supplementary alignments
    are dropped, and of paired reads the second mate, improper pairs and
    pairs with an unmapped mate; a minus-strand read's position is moved
    past its CIGAR M/D/N/=/X operations; records with refID -1 are
    dropped. Kept reads are gathered over one decompressed block at a
    time, grouped by chromosome (in order of first appearance) and strand
    with their order kept, and appended with ``fwtrack.add_loc_arrays``.

    Returns the number of reads added.
    """
    nref: cython.Py_ssize_t = len(references)
    total: cython.long = 0
    nextlog: cython.long = 1000000
    epoch: cython.int = 0
    need: cython.Py_ssize_t = 4
    kcap: cython.Py_ssize_t = 0
    maxrec: cython.Py_ssize_t
    nk: cython.Py_ssize_t
    nd: cython.Py_ssize_t
    k: cython.Py_ssize_t
    j: cython.Py_ssize_t
    r: cython.Py_ssize_t
    nc: cython.Py_ssize_t
    lrn: cython.Py_ssize_t
    bs: cython.Py_ssize_t
    cur: cython.int
    p: cython.p_uchar
    q: cython.p_uchar
    rec: cython.p_uchar
    cig: cython.p_uchar
    flag: cython.uint
    c: cython.uint
    upos: cython.uint
    refid: cython.int
    fpos: cython.int
    strand: cython.uchar
    # per reference: the block it was last seen in, its read count and
    # the next free slot in `grouped`, for each strand
    seen: cython.int[::1] = np.full(max(nref, 1), -1, dtype=np.int32)
    c0: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    c1: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    o0: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    o1: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    order: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    kref: cython.int[::1]
    kpos: cython.int[::1]
    kstr: cython.uchar[::1]
    grouped: cython.int[::1]

    while stream.refill(need):
        p = stream.buf + stream.start
        q = stream.buf + stream.end
        maxrec = (q - p) // (4 + _BAM_CORE_SIZE) + 1
        if maxrec > kcap:
            kcap = maxrec
            kref = np.empty(kcap, dtype=np.int32)
            kpos = np.empty(kcap, dtype=np.int32)
            kstr = np.empty(kcap, dtype=np.uint8)
            grouped_a = np.empty(kcap, dtype=np.int32)
            grouped = grouped_a

        # walk the complete records in the buffer
        nk = 0
        while q - p >= 4:
            bs = _le_int32(p)
            if bs < _BAM_CORE_SIZE:
                raise Exception("Invalid BAM record: block_size is %d" % bs)
            if q - p - 4 < bs:
                break
            rec = p + 4
            flag = rec[14] | (rec[15] << 8)
            if (flag & 2820) == 0 and \
               ((flag & 1) == 0 or ((flag & 136) == 0 and (flag & 2) != 0)):
                refid = _le_int32(rec)
                fpos = _le_int32(rec + 4)
                strand = 0
                if flag & 16:
                    # minus strand: move past the CIGAR M/D/N/=/X lengths
                    lrn = rec[8]
                    nc = rec[12] | (rec[13] << 8)
                    if _BAM_CORE_SIZE + lrn + 4 * nc > bs:
                        raise Exception("Invalid BAM record: CIGAR runs past the record")
                    cig = rec + _BAM_CORE_SIZE + lrn
                    upos = cython.cast(cython.uint, fpos)
                    for j in range(nc):
                        c = _le_uint32(cig + 4 * j)
                        if c > 0x7FFFFFFF:
                            raise OverflowError("value too large to convert to int")
                        if (0x18D >> (c & 15)) & 1:
                            upos += c >> 4
                    fpos = cython.cast(cython.int, upos)
                    strand = 1
                if refid != -1:
                    r = refid
                    if r < 0:
                        r += nref
                    if r < 0 or r >= nref:
                        raise IndexError("list index out of range")
                    kref[nk] = r
                    kpos[nk] = fpos
                    kstr[nk] = strand
                    nk += 1
            p += 4 + bs
        stream.start = p - stream.buf
        need = 4
        if q - p >= 4:
            need = 4 + max(_le_int32(p), 0)

        if nk == 0:
            continue
        # group this block's reads by chromosome and strand
        epoch += 1
        nd = 0
        for k in range(nk):
            r = kref[k]
            if seen[r] != epoch:
                seen[r] = epoch
                order[nd] = r
                nd += 1
                c0[r] = 0
                c1[r] = 0
            if kstr[k]:
                c1[r] += 1
            else:
                c0[r] += 1
        cur = 0
        for j in range(nd):
            r = order[j]
            o0[r] = cur
            cur += c0[r]
            o1[r] = cur
            cur += c1[r]
        for k in range(nk):
            r = kref[k]
            if kstr[k]:
                grouped[o1[r]] = kpos[k]
                o1[r] += 1
            else:
                grouped[o0[r]] = kpos[k]
                o0[r] += 1
        for j in range(nd):
            r = order[j]
            fwtrack.add_loc_arrays(references[r],
                                   grouped_a[o0[r] - c0[r]:o0[r]],
                                   grouped_a[o1[r] - c1[r]:o1[r]])
        total += nk
        while total >= nextlog:
            info(" %d reads parsed" % nextlog)
            nextlog += 1000000
    return total


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def _bam_load_petrack(stream: _BAMStream, references: list, petrack,
                      tlen_sum: cython.pointer(cython.long)) -> cython.long:
    """Add every kept alignment in ``stream`` to ``petrack`` as a fragment.

    The fragments kept, their coordinates and their order are those of
    parsing each record and calling ``petrack.add_loc(references[refID],
    start, start + abs(tlen))`` with ``start = min(pos, mate pos)``, with
    the same flag filter as ``_bam_load_fwtrack``; records with refID -1
    are dropped. Fragments are gathered over one decompressed block at a
    time, grouped by chromosome (in order of first appearance) with their
    order kept, and appended with ``petrack.add_loc_arrays``.

    Returns the number of fragments added, and adds the sum of their
    lengths to ``tlen_sum[0]``.
    """
    nref: cython.Py_ssize_t = len(references)
    total: cython.long = 0
    m: cython.long = 0
    nextlog: cython.long = 1000000
    epoch: cython.int = 0
    need: cython.Py_ssize_t = 4
    kcap: cython.Py_ssize_t = 0
    maxrec: cython.Py_ssize_t
    nk: cython.Py_ssize_t
    nd: cython.Py_ssize_t
    k: cython.Py_ssize_t
    j: cython.Py_ssize_t
    r: cython.Py_ssize_t
    bs: cython.Py_ssize_t
    cur: cython.int
    p: cython.p_uchar
    q: cython.p_uchar
    rec: cython.p_uchar
    flag: cython.uint
    refid: cython.int
    pos: cython.int
    nextpos: cython.int
    start: cython.int
    tlen: cython.int
    # per reference: the block it was last seen in, its fragment count
    # and the next free slot in `gl`/`gr`
    seen: cython.int[::1] = np.full(max(nref, 1), -1, dtype=np.int32)
    c0: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    o0: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    order: cython.int[::1] = np.zeros(max(nref, 1), dtype=np.int32)
    kref: cython.int[::1]
    kl: cython.int[::1]
    kr: cython.int[::1]
    gl: cython.int[::1]
    gr: cython.int[::1]

    while stream.refill(need):
        p = stream.buf + stream.start
        q = stream.buf + stream.end
        maxrec = (q - p) // (4 + _BAM_CORE_SIZE) + 1
        if maxrec > kcap:
            kcap = maxrec
            kref = np.empty(kcap, dtype=np.int32)
            kl = np.empty(kcap, dtype=np.int32)
            kr = np.empty(kcap, dtype=np.int32)
            gl_a = np.empty(kcap, dtype=np.int32)
            gr_a = np.empty(kcap, dtype=np.int32)
            gl = gl_a
            gr = gr_a

        # walk the complete records in the buffer
        nk = 0
        while q - p >= 4:
            bs = _le_int32(p)
            if bs < _BAM_CORE_SIZE:
                raise Exception("Invalid BAM record: block_size is %d" % bs)
            if q - p - 4 < bs:
                break
            rec = p + 4
            flag = rec[14] | (rec[15] << 8)
            if (flag & 2820) == 0 and \
               ((flag & 1) == 0 or ((flag & 136) == 0 and (flag & 2) != 0)):
                refid = _le_int32(rec)
                pos = _le_int32(rec + 4)
                nextpos = _le_int32(rec + 24)
                tlen = _le_int32(rec + 28)
                # the leftmost end, so no CIGAR is needed
                start = min(pos, nextpos)
                if tlen < 0:
                    tlen = cython.cast(cython.int,
                                       0 - cython.cast(cython.uint, tlen))
                if refid != -1:
                    r = refid
                    if r < 0:
                        r += nref
                    if r < 0 or r >= nref:
                        raise IndexError("list index out of range")
                    kref[nk] = r
                    kl[nk] = start
                    kr[nk] = cython.cast(cython.int,
                                         cython.cast(cython.uint, start) +
                                         cython.cast(cython.uint, tlen))
                    m += tlen
                    nk += 1
            p += 4 + bs
        stream.start = p - stream.buf
        need = 4
        if q - p >= 4:
            need = 4 + max(_le_int32(p), 0)

        if nk == 0:
            continue
        # group this block's fragments by chromosome
        epoch += 1
        nd = 0
        for k in range(nk):
            r = kref[k]
            if seen[r] != epoch:
                seen[r] = epoch
                order[nd] = r
                nd += 1
                c0[r] = 0
            c0[r] += 1
        cur = 0
        for j in range(nd):
            r = order[j]
            o0[r] = cur
            cur += c0[r]
        for k in range(nk):
            r = kref[k]
            gl[o0[r]] = kl[k]
            gr[o0[r]] = kr[k]
            o0[r] += 1
        for j in range(nd):
            r = order[j]
            petrack.add_loc_arrays(references[r],
                                   gl_a[o0[r] - c0[r]:o0[r]],
                                   gr_a[o0[r] - c0[r]:o0[r]])
        total += nk
        while total >= nextlog:
            info(" %d fragments parsed" % nextlog)
            nextlog += 1000000
    tlen_sum[0] += m
    return total

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

    @cython.cfunc
    def alignment_stream(self) -> _BAMStream:
        """Return a ``_BAMStream`` at the first alignment record.

        Call after ``get_references``, which leaves ``self.fhd`` there;
        ``self.fhd`` is closed.
        """
        offset: cython.Py_ssize_t = self.fhd.tell()
        stream: _BAMStream

        self.fhd.close()
        stream = _BAMStream(self.filename, self.gzipped)
        stream.skip(offset)
        return stream

    @cython.ccall
    def build_fwtrack(self):
        """Append uniquely mapped reads to an existing ``FWTrack``."""
        i: cython.long = 0          # number of reads kept
        references: list
        rlengths: dict

        fwtrack = FWTrack(buffer_size=self.buffer_size)
        # after this, ptr at list of alignments
        references, rlengths = self.get_references()
        stream = self.alignment_stream()
        i = _bam_load_fwtrack(stream, references, fwtrack)
        stream.close()

        info("%d reads have been read." % i)
        fwtrack.set_rlengths(rlengths)
        return fwtrack

    @cython.ccall
    def append_fwtrack(self, fwtrack):
        """Append uniquely mapped reads to an existing ``FWTrack``."""
        i: cython.long = 0          # number of reads kept
        references: list
        rlengths: dict

        references, rlengths = self.get_references()
        stream = self.alignment_stream()
        i = _bam_load_fwtrack(stream, references, fwtrack)
        stream.close()

        info("%d reads have been read." % i)
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
        references: list
        rlengths: dict

        petrack = PETrackI(buffer_size=self.buffer_size)
        references, rlengths = self.get_references()
        stream = self.alignment_stream()
        # for convenience, only count valid pairs
        i = _bam_load_petrack(stream, references, petrack, cython.address(m))
        stream.close()

        info("%d fragments have been read." % i)
        self.d = m / i
        self.n = i
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
        references: list
        rlengths: dict

        references, rlengths = self.get_references()
        stream = self.alignment_stream()
        # for convenience, only count valid pairs
        i = _bam_load_petrack(stream, references, petrack, cython.address(m))
        stream.close()

        info("%d fragments have been read." % i)
        self.d = (self.d * self.n + m) / (self.n + i)
        self.n += i
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


# ------------------------------------
# Fragment (FRAG) files, read in large decompressed blocks
# ------------------------------------

# fragments are gathered in batches of this many before they are
# appended to the track
_FRAG_BATCH = cython.declare(cython.Py_ssize_t, 262144)
# With the inflate pool running (a BGZF file on two or more cores), a
# window of at least _FRAG_PAR_MIN bytes is cut at line ends into
# _FRAG_TASKS_PER_THREAD chunks per thread, whose lines the pool's threads
# parse (_frag_chunk)
_FRAG_PAR_MIN = cython.declare(cython.Py_ssize_t, 65536)
_FRAG_TASKS_PER_THREAD = cython.declare(cython.int, 2)
# per chunk: at most _FRAG_NLC chromosomes not known before the window,
# and runs of one chromosome, up to _FRAG_NRUN, each appended in one
# add_loc_arrays call (a chunk of more runs goes through _FragBatch)
_FRAG_NLC = cython.declare(cython.int, 8)
_FRAG_NRUN = cython.declare(cython.int, 8)
# the fewest bytes of a line _frag_chunk parses: "c\t0\t1\t\t1\n"
_FRAG_MIN_LINE = cython.declare(cython.Py_ssize_t, 9)
# multipliers of _hash_bytes
_HASH_K0 = cython.declare(cython.ulonglong, 0x9E3779B97F4A7C15)
_HASH_K1 = cython.declare(cython.ulonglong, 0xBF58476D1CE4E5B9)
_HASH_K2 = cython.declare(cython.ulonglong, 0x94D049BB133111EB)
_HASH_K3 = cython.declare(cython.ulonglong, 0xD6E8FEB86659FD93)
# bytes 0x01, 0x09, 0x0a and 0x80 in every byte of a word, for finding
# tabs and newlines eight bytes at a time; byte-order dependent, so the
# word-at-a-time scan runs only on little-endian machines
_W_ONES = cython.declare(cython.ulonglong, 0x0101010101010101)
_W_TABS = cython.declare(cython.ulonglong, 0x0909090909090909)
_W_NLS = cython.declare(cython.ulonglong, 0x0A0A0A0A0A0A0A0A)
_W_HIGH = cython.declare(cython.ulonglong, 0x8080808080808080)
_W_SPACES = cython.declare(cython.ulonglong, 0x2020202020202020)
_W_INDEX = cython.declare(cython.ulonglong, 0x0001020304050607)
_W_ZEROS = cython.declare(cython.ulonglong, 0x3030303030303030)
_W_SEVENTY6 = cython.declare(cython.ulonglong, 0x7676767676767676)
_W_LOW1 = cython.declare(cython.ulonglong, 0x000000FF000000FF)
_W_MUL1 = cython.declare(cython.ulonglong, 0x000F424000000064)
_W_MUL2 = cython.declare(cython.ulonglong, 0x0000271000000001)
_LITTLE_ENDIAN = cython.declare(cython.bint, sys.byteorder == "little")


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _hash_round(h: cython.ulonglong, w: cython.ulonglong) -> cython.ulonglong:
    """One word of ``_hash_bytes``."""
    h = (h ^ w) * _HASH_K1
    return h ^ (h >> 31)


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _hash_words(k0: cython.ulonglong, k1: cython.ulonglong,
                k2: cython.ulonglong, n: cython.Py_ssize_t) -> cython.ulonglong:
    """``_hash_bytes`` of a key of ``n`` <= 24 bytes whose three
    zero-padded eight-byte words are ``k0``, ``k1`` and ``k2``: a sum of
    products, whose high bits, which ``_Interner`` uses, depend on every
    bit of the key."""
    return ((k0 * _HASH_K1) ^ (k1 * _HASH_K2) ^ (k2 * _HASH_K3) ^
            (cython.cast(cython.ulonglong, n) * _HASH_K0))


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _hash_bytes(p: cython.p_uchar, n: cython.Py_ssize_t) -> cython.ulonglong:
    """A 64-bit hash of ``p[0:n]``, for ``_Interner``, whose high bits
    depend on every byte: of its eight-byte words, the last one padded
    with zero bytes (``_hash_words`` for a key of at most 24 bytes)."""
    h: cython.ulonglong
    w: cython.ulonglong
    k: cython.Py_ssize_t = 0
    kb: cython.uchar[24]
    k0: cython.ulonglong
    k1: cython.ulonglong
    k2: cython.ulonglong

    if n <= 24:
        memset(cython.address(kb[0]), 0, 24)
        memcpy(cython.address(kb[0]), p, n)
        memcpy(cython.address(k0), cython.address(kb[0]), 8)
        memcpy(cython.address(k1), cython.address(kb[8]), 8)
        memcpy(cython.address(k2), cython.address(kb[16]), 8)
        return _hash_words(k0, k1, k2, n)
    h = cython.cast(cython.ulonglong, n) * _HASH_K0
    while k + 8 <= n:
        memcpy(cython.address(w), p + k, 8)
        h = _hash_round(h, w)
        k += 8
    if k < n:
        w = 0
        memcpy(cython.address(w), p + k, n - k)
        h = _hash_round(h, w)
    return h * _HASH_K2


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _tab_or_newline(w: cython.ulonglong) -> cython.ulonglong:
    """For a word read from memory on a little-endian machine: 0x80 in
    the lowest byte of ``w`` that is a tab or a newline (and maybe in
    bytes above it), or 0 when no byte is. A borrow can mark a byte
    only above a tab or newline, so the lowest mark is exact."""
    x: cython.ulonglong = w ^ _W_TABS
    y: cython.ulonglong = w ^ _W_NLS

    return (((x - _W_ONES) & ~x) | ((y - _W_ONES) & ~y)) & _W_HIGH


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _below_space(w: cython.ulonglong) -> cython.ulonglong:
    """For a word read from memory on a little-endian machine: 0x80 in
    the lowest byte of ``w`` below 0x20 (and maybe in bytes above it),
    or 0 when no byte is."""
    return (w - _W_SPACES) & ~w & _W_HIGH


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _low_bytes(n: cython.Py_ssize_t) -> cython.ulonglong:
    """A word whose lowest ``n`` bytes (none when ``n`` <= 0, all when
    ``n`` >= 8) are 0xff and the others 0."""
    if n <= 0:
        return 0
    if n >= 8:
        return ~cython.cast(cython.ulonglong, 0)
    return (cython.cast(cython.ulonglong, 1) << (8 * n)) - 1


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _lowest_mark(m: cython.ulonglong) -> cython.Py_ssize_t:
    """The index of the lowest byte of ``m`` (not 0) with 0x80 set."""
    return cython.cast(cython.Py_ssize_t,
                       (((m & (~m + 1)) >> 7) * _W_INDEX) >> 56)


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _same_bytes(a: cython.p_uchar, b: cython.p_uchar,
                n: cython.Py_ssize_t) -> cython.bint:
    """Whether ``a[0:n]`` and ``b[0:n]`` are equal, compared inline a
    word at a time, the last word overlapping the ones before."""
    x: cython.ulonglong
    y: cython.ulonglong
    x4: cython.uint
    y4: cython.uint
    k: cython.Py_ssize_t = 0

    if n >= 8:
        while k + 8 < n:
            memcpy(cython.address(x), a + k, 8)
            memcpy(cython.address(y), b + k, 8)
            if x != y:
                return False
            k += 8
        memcpy(cython.address(x), a + n - 8, 8)
        memcpy(cython.address(y), b + n - 8, 8)
        return x == y
    if n >= 4:
        memcpy(cython.address(x4), a, 4)
        memcpy(cython.address(y4), b, 4)
        if x4 != y4:
            return False
        memcpy(cython.address(x4), a + n - 4, 4)
        memcpy(cython.address(y4), b + n - 4, 4)
        return x4 == y4
    while k < n:
        if a[k] != b[k]:
            return False
        k += 1
    return True


# One slot of _Interner's table, a 64-byte cache line: the key's hash,
# its length + 1 (0: the slot is empty), its offset in the arena (keys
# longer than _ISLOT_KEY bytes), its value, and the key itself when it
# is at most _ISLOT_KEY bytes long
_ISlot = cython.struct(h=cython.ulonglong, n1=cython.Py_ssize_t,
                       off=cython.Py_ssize_t, value=cython.int,
                       key=cython.uchar[36])
_ISLOT_KEY = cython.declare(cython.Py_ssize_t, 36)


@cython.final
@cython.cclass
class _Interner:
    """A table from byte strings, given as a pointer and a length, to
    non-negative ints, so that looking up a string needs no bytes
    object. Open addressing with linear probing, at most a quarter of
    the slots used; a slot holds the key's hash, value and length and,
    for a key of up to ``_ISLOT_KEY`` bytes, the key, so that most
    lookups read one cache line. Longer keys are copied into an arena."""
    slots: cython.pointer(_ISlot)   # 64-byte aligned, in raw
    raw: cython.p_void
    mask: cython.Py_ssize_t         # number of slots - 1
    shift: cython.int               # a key's slot is its hash >> shift
    n: cython.Py_ssize_t
    arena: cython.p_uchar
    alen: cython.Py_ssize_t
    acap: cython.Py_ssize_t

    def __cinit__(self):
        self.raw = cython.NULL
        self.slots = cython.NULL
        self.mask = 0
        self.shift = 64
        self.n = 0
        self.alen = 0
        self.acap = 0
        self.arena = cython.NULL
        self.new_slots(1024)

    def __dealloc__(self):
        free(self.raw)
        free(self.arena)

    @cython.cfunc
    @cython.exceptval(-1)
    def new_slots(self, nslots: cython.Py_ssize_t) -> cython.int:
        """Make the table ``nslots`` empty slots, a power of two, and put
        the keys of the old table back in it, in the old slots' order.
        Returns 0."""
        raw: cython.p_void
        slots: cython.pointer(_ISlot)
        old: cython.pointer(_ISlot) = self.slots
        oldraw: cython.p_void = self.raw
        oldn: cython.Py_ssize_t = self.mask + 1 if old != cython.NULL else 0
        mask: cython.Py_ssize_t = nslots - 1
        shift: cython.int = 64
        k: cython.Py_ssize_t
        j: cython.Py_ssize_t

        raw = calloc(nslots + 1, cython.sizeof(_ISlot))
        if raw == cython.NULL:
            raise MemoryError()
        slots = cython.cast(cython.pointer(_ISlot),
                            (cython.cast(cython.size_t, raw) + 63) &
                            ~cython.cast(cython.size_t, 63))
        k = nslots
        while k > 1:
            shift -= 1
            k >>= 1
        for k in range(oldn):
            if old[k].n1 != 0:
                j = old[k].h >> shift
                while slots[j].n1 != 0:
                    j = (j + 1) & mask
                slots[j] = old[k]
        self.slots = slots
        self.raw = raw
        self.mask = mask
        self.shift = shift
        free(oldraw)
        return 0

    @cython.cfunc
    @cython.inline
    @cython.exceptval(check=False)
    @cython.profile(False)
    def find(self, p: cython.p_uchar, n: cython.Py_ssize_t,
             h: cython.ulonglong) -> cython.int:
        """The value of ``p[0:n]``, whose ``_hash_bytes`` is ``h``, or -1."""
        j: cython.Py_ssize_t = h >> self.shift
        s: cython.pointer(_ISlot)

        while True:
            s = self.slots + j
            if s.n1 == 0:
                return -1
            if s.h == h and s.n1 == n + 1:
                if n <= _ISLOT_KEY:
                    if _same_bytes(cython.address(s.key[0]), p, n):
                        return s.value
                elif _same_bytes(self.arena + s.off, p, n):
                    return s.value
            j = (j + 1) & self.mask

    @cython.cfunc
    @cython.inline
    @cython.exceptval(check=False)
    @cython.profile(False)
    def find_words(self, k0: cython.ulonglong, k1: cython.ulonglong,
                   k2: cython.ulonglong, n: cython.Py_ssize_t,
                   h: cython.ulonglong) -> cython.int:
        """``find`` for a key of ``n`` <= 24 bytes given as its three
        zero-padded words (``_hash_words``). A slot's key bytes past the
        key's length are zero, so its three words are compared."""
        j: cython.Py_ssize_t = h >> self.shift
        s: cython.pointer(_ISlot)
        x0: cython.ulonglong
        x1: cython.ulonglong
        x2: cython.ulonglong

        while True:
            s = self.slots + j
            if s.n1 == 0:
                return -1
            if s.h == h and s.n1 == n + 1:
                memcpy(cython.address(x0), cython.address(s.key[0]), 8)
                memcpy(cython.address(x1), cython.address(s.key[8]), 8)
                memcpy(cython.address(x2), cython.address(s.key[16]), 8)
                if x0 == k0 and x1 == k1 and x2 == k2:
                    return s.value
            j = (j + 1) & self.mask

    @cython.cfunc
    @cython.exceptval(-1)
    def add(self, p: cython.p_uchar, n: cython.Py_ssize_t,
            h: cython.ulonglong, value: cython.int) -> cython.int:
        """Add ``p[0:n]``, which is not in the table and whose
        ``_hash_bytes`` is ``h``, with ``value``. Returns 0."""
        newcap: cython.Py_ssize_t
        j: cython.Py_ssize_t
        s: cython.pointer(_ISlot)
        tmp: cython.p_void

        if 4 * (self.n + 1) > self.mask + 1:
            # keep the table at most a quarter full
            self.new_slots(2 * (self.mask + 1))
        j = h >> self.shift
        while self.slots[j].n1 != 0:
            j = (j + 1) & self.mask
        s = self.slots + j
        if n <= _ISLOT_KEY:
            memcpy(cython.address(s.key[0]), p, n)
        else:
            if self.alen + n > self.acap:
                newcap = self.acap if self.acap > 0 else 16384
                while self.alen + n > newcap:
                    newcap *= 2
                tmp = realloc(self.arena, newcap)
                if tmp == cython.NULL:
                    raise MemoryError()
                self.arena = cython.cast(cython.p_uchar, tmp)
                self.acap = newcap
            memcpy(self.arena + self.alen, p, n)
            s.off = self.alen
            self.alen += n
        s.h = h
        s.value = value
        s.n1 = n + 1
        self.n += 1
        return 0


# A read-only view of an _Interner's table, for lookups without the GIL
# (the FRAG reader's threads), made by _itab_of while nothing adds keys
_ITab = cython.struct(slots=cython.pointer(_ISlot), mask=cython.Py_ssize_t,
                      shift=cython.int, arena=cython.p_uchar)


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
def _itab_of(t: cython.pointer(_ITab), it: _Interner) -> cython.void:
    """Point ``t`` at the table of ``it`` as it is now."""
    t.slots = it.slots
    t.mask = it.mask
    t.shift = it.shift
    t.arena = it.arena


@cython.cfunc
@cython.inline
@cython.nogil
@cython.exceptval(check=False)
@cython.profile(False)
def _itab_find(t: cython.pointer(_ITab), p: cython.p_uchar,
               n: cython.Py_ssize_t, h: cython.ulonglong) -> cython.int:
    """``_Interner.find`` on the view ``t``."""
    j: cython.Py_ssize_t = h >> t.shift
    s: cython.pointer(_ISlot)

    while True:
        s = t.slots + j
        if s.n1 == 0:
            return -1
        if s.h == h and s.n1 == n + 1:
            if n <= _ISLOT_KEY:
                if _same_bytes(cython.address(s.key[0]), p, n):
                    return s.value
            elif _same_bytes(t.arena + s.off, p, n):
                return s.value
        j = (j + 1) & t.mask


@cython.cfunc
@cython.inline
@cython.nogil
@cython.exceptval(check=False)
@cython.profile(False)
def _itab_find_words(t: cython.pointer(_ITab), k0: cython.ulonglong,
                     k1: cython.ulonglong, k2: cython.ulonglong,
                     n: cython.Py_ssize_t, h: cython.ulonglong) -> cython.int:
    """``_Interner.find_words`` on the view ``t``."""
    j: cython.Py_ssize_t = h >> t.shift
    s: cython.pointer(_ISlot)
    x0: cython.ulonglong
    x1: cython.ulonglong
    x2: cython.ulonglong

    while True:
        s = t.slots + j
        if s.n1 == 0:
            return -1
        if s.h == h and s.n1 == n + 1:
            memcpy(cython.address(x0), cython.address(s.key[0]), 8)
            memcpy(cython.address(x1), cython.address(s.key[8]), 8)
            memcpy(cython.address(x2), cython.address(s.key[16]), 8)
            if x0 == k0 and x1 == k1 and x2 == k2:
                return s.value
        j = (j + 1) & t.mask


@cython.cclass
class _FragBatch:
    """Fragments waiting to be appended to a ``PETrackII``, in file order:
    a chromosome index (into the list given to ``flush``), the two ends,
    the count and the barcode id of each."""
    n: cython.Py_ssize_t
    cap: cython.Py_ssize_t
    chrom: cython.int[::1]
    left: cython.int[::1]
    right: cython.int[::1]
    count: cython.ushort[::1]
    bc: cython.int[::1]
    left_a: object
    right_a: object
    count_a: object
    bc_a: object
    # the same, grouped by chromosome
    gleft: cython.int[::1]
    gright: cython.int[::1]
    gcount: cython.ushort[::1]
    gbc: cython.int[::1]
    gleft_a: object
    gright_a: object
    gcount_a: object
    gbc_a: object
    # per chromosome: the flush it was last seen in, its fragment count
    # and next free slot in this flush, and the order of first appearance
    seen: cython.int[::1]
    ccount: cython.int[::1]
    coff: cython.int[::1]
    order: cython.int[::1]
    epoch: cython.int

    def __init__(self, cap: cython.Py_ssize_t):
        self.n = 0
        self.cap = cap
        self.chrom = np.empty(cap, dtype=np.int32)
        self.left_a = np.empty(cap, dtype=np.int32)
        self.right_a = np.empty(cap, dtype=np.int32)
        self.count_a = np.empty(cap, dtype=np.uint16)
        self.bc_a = np.empty(cap, dtype=np.int32)
        self.left = self.left_a
        self.right = self.right_a
        self.count = self.count_a
        self.bc = self.bc_a
        self.gleft_a = np.empty(cap, dtype=np.int32)
        self.gright_a = np.empty(cap, dtype=np.int32)
        self.gcount_a = np.empty(cap, dtype=np.uint16)
        self.gbc_a = np.empty(cap, dtype=np.int32)
        self.gleft = self.gleft_a
        self.gright = self.gright_a
        self.gcount = self.gcount_a
        self.gbc = self.gbc_a
        self.size_chroms(64)
        self.epoch = 0

    @cython.cfunc
    def size_chroms(self, nc: cython.Py_ssize_t):
        """Make the per-chromosome arrays hold ``nc`` chromosomes."""
        self.seen = np.full(nc, -1, dtype=np.int32)
        self.ccount = np.zeros(nc, dtype=np.int32)
        self.coff = np.zeros(nc, dtype=np.int32)
        self.order = np.zeros(nc, dtype=np.int32)

    @cython.cfunc
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    def flush(self, petrack, chroms: list):
        """Append the waiting fragments to ``petrack``, one
        ``add_loc_arrays`` call per chromosome in order of first
        appearance, each chromosome's fragments in file order."""
        n: cython.Py_ssize_t = self.n
        nd: cython.Py_ssize_t = 0
        k: cython.Py_ssize_t
        j: cython.Py_ssize_t
        r: cython.int
        o: cython.int
        cur: cython.int

        if n == 0:
            return
        if self.seen.shape[0] < len(chroms):
            self.size_chroms(2 * len(chroms))
        self.epoch += 1
        for k in range(n):
            r = self.chrom[k]
            if self.seen[r] != self.epoch:
                self.seen[r] = self.epoch
                self.order[nd] = r
                nd += 1
                self.ccount[r] = 0
            self.ccount[r] += 1
        self.n = 0
        if nd == 1:
            petrack.add_loc_arrays(chroms[self.order[0]], self.left_a[:n],
                                   self.right_a[:n], self.count_a[:n],
                                   self.bc_a[:n])
            return
        cur = 0
        for j in range(nd):
            r = self.order[j]
            self.coff[r] = cur
            cur += self.ccount[r]
        for k in range(n):
            r = self.chrom[k]
            o = self.coff[r]
            self.gleft[o] = self.left[k]
            self.gright[o] = self.right[k]
            self.gcount[o] = self.count[k]
            self.gbc[o] = self.bc[k]
            self.coff[r] = o + 1
        for j in range(nd):
            r = self.order[j]
            o = self.coff[r]
            cur = o - self.ccount[r]
            petrack.add_loc_arrays(chroms[r], self.gleft_a[cur:o],
                                   self.gright_a[cur:o], self.gcount_a[cur:o],
                                   self.gbc_a[cur:o])


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
def _frag_atoi(a: cython.p_uchar, b: cython.p_uchar) -> cython.int:
    """``atoi`` of the field ``a[0:b - a]``, as ``atoi`` of that field as
    a bytes object gives it. ``b`` must point into a writable buffer."""
    n: cython.Py_ssize_t = b - a
    v: cython.int = 0
    d: cython.uint
    k: cython.Py_ssize_t
    c: cython.uchar

    if 0 < n <= 9:
        # up to nine decimal digits: the value, below 2^31
        for k in range(n):
            d = cython.cast(cython.uint, a[k]) - 48
            if d > 9:
                break
            v = v * 10 + cython.cast(cython.int, d)
        else:
            return v
    # anything else (signs, spaces, long or empty fields): libc's atoi on
    # the field alone, terminated where its bytes object would be
    c = b[0]
    b[0] = 0
    v = atoi(cython.cast(cython.p_char, a))
    b[0] = c
    return v


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _frag_digits(s: cython.p_uchar, v: cython.p_int) -> cython.p_uchar:
    """The end of the run of ASCII digits at ``s``, with its value in
    ``v[0]``, when the run is one to nine digits long; NULL otherwise.
    A byte that is not a digit must follow the run."""
    a: cython.p_uchar = s
    x: cython.uint = 0
    d: cython.uint = cython.cast(cython.uint, s[0]) - 48

    while d <= 9:
        x = x * 10 + d
        s += 1
        d = cython.cast(cython.uint, s[0]) - 48
    if s == a or s - a > 9:
        return cython.NULL
    v[0] = cython.cast(cython.int, x)
    return s


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _frag_digits_w(s: cython.p_uchar, v: cython.p_int) -> cython.p_uchar:
    """``_frag_digits`` eight bytes at a time, on a little-endian
    machine; ``s[0:8]`` must be readable. The digits' values are
    ``s[k] - 48``; the lowest byte that is not a digit is found exactly
    (a borrow or carry can only mark bytes above it), and up to eight
    digits are converted together."""
    w: cython.ulonglong
    x: cython.ulonglong
    m: cython.ulonglong
    n: cython.Py_ssize_t
    d: cython.uint

    memcpy(cython.address(w), s, 8)
    x = w - _W_ZEROS
    m = ((x + _W_SEVENTY6) | x) & _W_HIGH
    if m != 0:
        n = _lowest_mark(m)
        if n == 0:
            return cython.NULL
        # the n digits to the top of the word, zeros before them
        x <<= 8 * (8 - n)
    else:
        n = 8
    x = x * 10 + (x >> 8)
    x = (((x & _W_LOW1) * _W_MUL1) + (((x >> 16) & _W_LOW1) * _W_MUL2)) >> 32
    if n == 8:
        d = cython.cast(cython.uint, s[8]) - 48
        if d <= 9:
            # a ninth digit; a tenth makes the run too long
            if cython.cast(cython.uint, s[9]) - 48 <= 9:
                return cython.NULL
            x = x * 10 + d
            n = 9
    v[0] = cython.cast(cython.int, x)
    return s + n


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
@cython.profile(False)
@cython.nogil
def _frag_digits_tab(s: cython.p_uchar, n: cython.Py_ssize_t,
                     v: cython.p_int) -> cython.p_uchar:
    """``s + n``, with the value of ``s[0:n]`` in ``v[0]``, when
    ``s[0:n]`` are digits and ``s[n]`` a tab, for a guessed ``n`` from 1
    to 9; NULL otherwise. On a little-endian machine; ``s[0:10]`` must
    be readable. Where the field ends is the guess, which the checks
    confirm, so the next field's address does not wait for them."""
    w: cython.ulonglong
    x: cython.ulonglong
    m: cython.ulonglong
    d: cython.uint

    if s[n] != 9:
        return cython.NULL
    memcpy(cython.address(w), s, 8)
    x = w - _W_ZEROS
    # the lowest byte that is not a digit is marked exactly
    m = ((x + _W_SEVENTY6) | x) & _W_HIGH
    if n < 8:
        if (m & ((cython.cast(cython.ulonglong, 1) << (8 * n)) - 1)) != 0:
            return cython.NULL
        x <<= 8 * (8 - n)
    elif m != 0:
        return cython.NULL
    x = x * 10 + (x >> 8)
    x = (((x & _W_LOW1) * _W_MUL1) + (((x >> 16) & _W_LOW1) * _W_MUL2)) >> 32
    if n == 9:
        d = cython.cast(cython.uint, s[8]) - 48
        if d > 9:
            return cython.NULL
        x = x * 10 + d
    v[0] = cython.cast(cython.int, x)
    return s + n


# One window of the threaded line walk: the chunks, chunk j being the
# lines in cuts[j]:cuts[j + 1], and what _frag_chunk leaves for each.
# Chunk j's records are at base[j]: base[j] + n[j] of the record arrays
# (chromosome id, ends, count, barcode id). An id below 0 is chunk j's
# chromosome -1 - id of those not in ctab (lcp/lcn[j * _FRAG_NLC + id]),
# or, for a barcode not in btab, -1 with the barcode in bptr/blen.
_FragJob = cython.struct(
    stop=cython.p_uchar,           # one past the window's last newline
    cuts=cython.pointer(cython.p_uchar),
    base=cython.pointer(cython.Py_ssize_t),
    capped=cython.bint,
    cap=cython.int,
    ctab=_ITab,
    btab=_ITab,
    chrom=cython.p_int,
    left=cython.p_int,
    right=cython.p_int,
    count=cython.pointer(cython.ushort),
    bc=cython.p_int,
    bptr=cython.pointer(cython.p_uchar),
    blen=cython.p_int,
    # per chunk
    n=cython.pointer(cython.Py_ssize_t),
    msum=cython.p_long,            # sum of right - left
    nmiss=cython.pointer(cython.Py_ssize_t),
    fail=cython.pointer(cython.p_uchar),   # the first line not parsed
    nlc=cython.p_int,
    lcp=cython.pointer(cython.p_uchar),
    lcn=cython.pointer(cython.Py_ssize_t),
    nrun=cython.p_int,             # above _FRAG_NRUN: too many to list
    runc=cython.p_int,             # [j * _FRAG_NRUN + k]: run k's id
    runend=cython.pointer(cython.Py_ssize_t))   # and its end


@cython.cfunc
@cython.nogil
@cython.exceptval(check=False)
def _frag_chunk(arg: cython.p_void, j: cython.size_t) -> cython.void:
    """Parse the lines of chunk ``j`` of the window described by
    ``arg`` (a ``_FragJob``) as ``FragParser.load_petrack``'s fast path
    parses them, until the first line it does not take, which, and every
    line after it in the chunk, is left to the walk on the calling
    thread. Runs on the inflate pool's threads, without the GIL: it
    reads the window and the tables, and writes only chunk ``j``'s
    records and results."""
    job: cython.pointer(_FragJob) = cython.cast(cython.pointer(_FragJob), arg)
    p: cython.p_uchar = job.cuts[j]
    end: cython.p_uchar = job.cuts[j + 1]
    stop: cython.p_uchar = job.stop
    r: cython.Py_ssize_t = job.base[j]
    nl: cython.p_uchar
    t: cython.p_uchar
    s: cython.p_uchar
    bstart: cython.p_uchar
    left: cython.int
    right: cython.int
    count: cython.int
    clen: cython.Py_ssize_t
    blen: cython.Py_ssize_t
    ok: cython.bint
    words: cython.bint
    csame: cython.bint
    found: cython.bint
    mark: cython.ulonglong
    k0: cython.ulonglong
    k1: cython.ulonglong
    k2: cython.ulonglong
    h: cython.ulonglong
    cw: cython.ulonglong
    ci: cython.int
    bi: cython.int
    k: cython.int
    havelast: cython.bint = False
    lastc: cython.int = 0
    lastp: cython.p_uchar = cython.NULL
    lastn: cython.Py_ssize_t = 8
    lastw: cython.ulonglong = 0
    lastm: cython.ulonglong = 0
    gl: cython.Py_ssize_t = 8
    gr: cython.Py_ssize_t = 8
    gb: cython.Py_ssize_t = 18
    bm0: cython.ulonglong = _low_bytes(18)
    bm1: cython.ulonglong = _low_bytes(10)
    bm2: cython.ulonglong = _low_bytes(2)
    nlc: cython.int = 0
    nrun: cython.int = 0
    curc: cython.int = 0
    nmiss: cython.Py_ssize_t = 0
    msum: cython.long = 0
    lcp: cython.pointer(cython.p_uchar) = job.lcp + j * _FRAG_NLC
    lcn: cython.pointer(cython.Py_ssize_t) = job.lcn + j * _FRAG_NLC
    runc: cython.p_int = job.runc + j * _FRAG_NRUN
    runend: cython.pointer(cython.Py_ssize_t) = job.runend + j * _FRAG_NRUN

    while p < end:
        # the fast path of load_petrack, line for line
        ok = False
        words = False
        csame = False
        t = p
        if (_LITTLE_ENDIAN and lastn < 8 and p + 8 <= stop and
                p[lastn] == 9):
            memcpy(cython.address(cw), p, 8)
            if (cw & lastm) == lastw:
                t = p + lastn
                csame = True
        if not csame:
            while t[0] != 9 and t[0] != 10:
                t += 1
        clen = t - p
        if t[0] == 9 and clen > 0:
            s = t + 1
            t = cython.NULL
            if _LITTLE_ENDIAN and s + 10 <= stop:
                t = _frag_digits_tab(s, gl, cython.address(left))
                if t == cython.NULL:
                    t = _frag_digits_w(s, cython.address(left))
                    if t != cython.NULL:
                        gl = t - s
            else:
                t = _frag_digits(s, cython.address(left))
            if t != cython.NULL and t[0] == 9:
                s = t + 1
                t = cython.NULL
                if _LITTLE_ENDIAN and s + 10 <= stop:
                    t = _frag_digits_tab(s, gr, cython.address(right))
                    if t == cython.NULL:
                        t = _frag_digits_w(s, cython.address(right))
                        if t != cython.NULL:
                            gr = t - s
                else:
                    t = _frag_digits(s, cython.address(right))
                if t != cython.NULL and t[0] == 9 and right > left:
                    bstart = t + 1
                    t = bstart
                    if (_LITTLE_ENDIAN and bstart + 24 <= stop and
                            bstart[gb] == 9):
                        memcpy(cython.address(k0), bstart, 8)
                        memcpy(cython.address(k1), bstart + 8, 8)
                        memcpy(cython.address(k2), bstart + 16, 8)
                        if (((_below_space(k0) & bm0) |
                             (_below_space(k1) & bm1) |
                             (_below_space(k2) & bm2)) == 0):
                            k0 &= bm0
                            k1 &= bm1
                            k2 &= bm2
                            blen = gb
                            words = True
                    if not words and _LITTLE_ENDIAN and bstart + 24 <= stop:
                        memcpy(cython.address(k0), bstart, 8)
                        mark = _tab_or_newline(k0)
                        if mark != 0:
                            blen = _lowest_mark(mark)
                            k0 &= (cython.cast(cython.ulonglong, 1) << (8 * blen)) - 1
                            k1 = 0
                            k2 = 0
                            words = True
                        else:
                            memcpy(cython.address(k1), bstart + 8, 8)
                            mark = _tab_or_newline(k1)
                            if mark != 0:
                                blen = _lowest_mark(mark)
                                k1 &= (cython.cast(cython.ulonglong, 1) << (8 * blen)) - 1
                                k2 = 0
                                blen += 8
                                words = True
                            else:
                                memcpy(cython.address(k2), bstart + 16, 8)
                                mark = _tab_or_newline(k2)
                                if mark != 0:
                                    blen = _lowest_mark(mark)
                                    k2 &= (cython.cast(cython.ulonglong, 1) << (8 * blen)) - 1
                                    blen += 16
                                    words = True
                                else:
                                    t = bstart + 24
                        if words:
                            gb = blen
                            bm0 = _low_bytes(blen)
                            bm1 = _low_bytes(blen - 8)
                            bm2 = _low_bytes(blen - 16)
                    if words:
                        t = bstart + blen
                    while t[0] != 9 and t[0] != 10:
                        t += 1
                    if t[0] == 9:
                        blen = t - bstart
                        nl = _frag_digits(t + 1, cython.address(count))
                        if nl != cython.NULL:
                            while nl[0] != 10:
                                nl += 1
                            ok = True
        if not ok:
            break
        count = cython.cast(cython.ushort, count)
        if job.capped and job.cap < count:
            count = job.cap
        # the chromosome: the last line's, one in ctab, or one of this
        # chunk's new ones
        if csame or (havelast and lastn == clen and
                     _same_bytes(lastp, p, clen)):
            ci = lastc
        else:
            h = _hash_bytes(p, clen)
            ci = _itab_find(cython.address(job.ctab), p, clen, h)
            if ci < 0:
                found = False
                for k in range(nlc):
                    if lcn[k] == clen and _same_bytes(lcp[k], p, clen):
                        ci = -1 - k
                        found = True
                        break
                if not found:
                    if nlc == _FRAG_NLC:
                        break
                    lcp[nlc] = p
                    lcn[nlc] = clen
                    ci = -1 - nlc
                    nlc += 1
            havelast = True
            lastc = ci
            lastp = p
            lastn = clen
            if clen < 8:
                lastw = 0
                memcpy(cython.address(lastw), p, clen)
                lastm = _low_bytes(clen)
        if words:
            h = _hash_words(k0, k1, k2, blen)
            bi = _itab_find_words(cython.address(job.btab), k0, k1, k2, blen, h)
        else:
            h = _hash_bytes(bstart, blen)
            bi = _itab_find(cython.address(job.btab), bstart, blen, h)
        if bi < 0:
            job.bptr[r] = bstart
            job.blen[r] = cython.cast(cython.int, blen)
            nmiss += 1
        if nrun == 0 or ci != curc:
            if 0 < nrun <= _FRAG_NRUN:
                runend[nrun - 1] = r
            if nrun < _FRAG_NRUN:
                runc[nrun] = ci
            if nrun <= _FRAG_NRUN:
                nrun += 1
            curc = ci
        job.chrom[r] = ci
        job.left[r] = left
        job.right[r] = right
        job.count[r] = cython.cast(cython.ushort, count)
        job.bc[r] = bi
        r += 1
        msum += right - left
        p = nl + 1
    if 0 < nrun <= _FRAG_NRUN:
        runend[nrun - 1] = r
    job.fail[j] = p
    job.n[j] = r - job.base[j]
    job.msum[j] = msum
    job.nmiss[j] = nmiss
    job.nlc[j] = nlc
    job.nrun[j] = nrun


@cython.final
@cython.cclass
class _FragChunks:
    """The buffers of the threaded line walk: a ``_FragJob`` for up to
    ``nmax`` chunks and records for one window."""
    job: _FragJob
    nmax: cython.Py_ssize_t
    rcap: cython.Py_ssize_t
    chrom_a: object
    left_a: object
    right_a: object
    count_a: object
    bc_a: object
    blen_a: object
    bptr_b: bytearray

    def __cinit__(self):
        memset(cython.address(self.job), 0, cython.sizeof(_FragJob))

    def __init__(self, nmax: cython.Py_ssize_t):
        self.nmax = nmax
        self.job.cuts = cython.cast(cython.pointer(cython.p_uchar),
                                    calloc(nmax + 1, cython.sizeof(cython.p_uchar)))
        self.job.base = cython.cast(cython.pointer(cython.Py_ssize_t),
                                    calloc(nmax + 1, cython.sizeof(cython.Py_ssize_t)))
        self.job.n = cython.cast(cython.pointer(cython.Py_ssize_t),
                                 calloc(nmax, cython.sizeof(cython.Py_ssize_t)))
        self.job.msum = cython.cast(cython.p_long,
                                    calloc(nmax, cython.sizeof(cython.long)))
        self.job.nmiss = cython.cast(cython.pointer(cython.Py_ssize_t),
                                     calloc(nmax, cython.sizeof(cython.Py_ssize_t)))
        self.job.fail = cython.cast(cython.pointer(cython.p_uchar),
                                    calloc(nmax, cython.sizeof(cython.p_uchar)))
        self.job.nlc = cython.cast(cython.p_int,
                                   calloc(nmax, cython.sizeof(cython.int)))
        self.job.lcp = cython.cast(cython.pointer(cython.p_uchar),
                                   calloc(nmax * _FRAG_NLC, cython.sizeof(cython.p_uchar)))
        self.job.lcn = cython.cast(cython.pointer(cython.Py_ssize_t),
                                   calloc(nmax * _FRAG_NLC, cython.sizeof(cython.Py_ssize_t)))
        self.job.nrun = cython.cast(cython.p_int,
                                    calloc(nmax, cython.sizeof(cython.int)))
        self.job.runc = cython.cast(cython.p_int,
                                    calloc(nmax * _FRAG_NRUN, cython.sizeof(cython.int)))
        self.job.runend = cython.cast(cython.pointer(cython.Py_ssize_t),
                                      calloc(nmax * _FRAG_NRUN, cython.sizeof(cython.Py_ssize_t)))
        if (self.job.cuts == cython.NULL or self.job.base == cython.NULL or
                self.job.n == cython.NULL or self.job.msum == cython.NULL or
                self.job.nmiss == cython.NULL or self.job.fail == cython.NULL or
                self.job.nlc == cython.NULL or self.job.lcp == cython.NULL or
                self.job.lcn == cython.NULL or self.job.nrun == cython.NULL or
                self.job.runc == cython.NULL or self.job.runend == cython.NULL):
            raise MemoryError()
        self.rcap = 0

    def __dealloc__(self):
        free(self.job.cuts)
        free(self.job.base)
        free(self.job.n)
        free(self.job.msum)
        free(self.job.nmiss)
        free(self.job.fail)
        free(self.job.nlc)
        free(self.job.lcp)
        free(self.job.lcn)
        free(self.job.nrun)
        free(self.job.runc)
        free(self.job.runend)

    @cython.cfunc
    def cut(self, p: cython.p_uchar, stop: cython.p_uchar,
            nch: cython.Py_ssize_t):
        """Cut the lines in ``p:stop`` (``stop`` one past a newline) into
        ``nch`` chunks of about equal size at line ends, place each
        chunk's records and make the record arrays hold them."""
        j: cython.Py_ssize_t
        x: cython.p_uchar
        need: cython.Py_ssize_t
        size: cython.Py_ssize_t = stop - p
        chrom_v: cython.int[::1]
        left_v: cython.int[::1]
        right_v: cython.int[::1]
        count_v: cython.ushort[::1]
        bc_v: cython.int[::1]
        blen_v: cython.int[::1]

        self.job.cuts[0] = p
        for j in range(1, nch):
            x = p + size * j // nch
            if x < self.job.cuts[j - 1]:
                x = self.job.cuts[j - 1]
            elif x > p:
                # the line that holds x - 1 ends the chunk before
                x = cython.cast(cython.p_uchar, memchr(x - 1, 10, stop - (x - 1))) + 1
            self.job.cuts[j] = x
        self.job.cuts[nch] = stop
        # a chunk of b bytes holds at most b // _FRAG_MIN_LINE records
        self.job.base[0] = 0
        for j in range(nch):
            self.job.base[j + 1] = (self.job.base[j] + 1 +
                                    (self.job.cuts[j + 1] - self.job.cuts[j]) //
                                    _FRAG_MIN_LINE)
        need = self.job.base[nch]
        if need > self.rcap:
            self.rcap = max(need, 2 * self.rcap)
            self.chrom_a = np.empty(self.rcap, dtype=np.int32)
            self.left_a = np.empty(self.rcap, dtype=np.int32)
            self.right_a = np.empty(self.rcap, dtype=np.int32)
            self.count_a = np.empty(self.rcap, dtype=np.uint16)
            self.bc_a = np.empty(self.rcap, dtype=np.int32)
            self.blen_a = np.empty(self.rcap, dtype=np.int32)
            self.bptr_b = bytearray(self.rcap * cython.sizeof(cython.p_uchar))
            chrom_v = self.chrom_a
            left_v = self.left_a
            right_v = self.right_a
            count_v = self.count_a
            bc_v = self.bc_a
            blen_v = self.blen_a
            self.job.chrom = cython.address(chrom_v[0])
            self.job.left = cython.address(left_v[0])
            self.job.right = cython.address(right_v[0])
            self.job.count = cython.address(count_v[0])
            self.job.bc = cython.address(bc_v[0])
            self.job.blen = cython.address(blen_v[0])
            self.job.bptr = cython.cast(cython.pointer(cython.p_uchar),
                                        PyByteArray_AS_STRING(self.bptr_b))
        self.job.stop = stop

    @cython.cfunc
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def emit(self, j: cython.Py_ssize_t, petrack, chroms: list,
             ctab: _Interner, btab: _Interner,
             batch: _FragBatch) -> cython.Py_ssize_t:
        """Append the records of chunk ``j`` to ``petrack`` after giving
        its new chromosomes and barcodes their ids, in file order, as the
        walk on the calling thread does; returns how many there are."""
        job: cython.pointer(_FragJob) = cython.address(self.job)
        b0: cython.Py_ssize_t = job.base[j]
        n: cython.Py_ssize_t = job.n[j]
        r: cython.Py_ssize_t
        a: cython.Py_ssize_t
        e: cython.Py_ssize_t
        k: cython.Py_ssize_t
        nr: cython.int
        lp: cython.p_uchar
        ln: cython.Py_ssize_t
        h: cython.ulonglong
        ci: cython.int
        bi: cython.int
        lmap: cython.int[8]

        # the chunk's chromosomes not in ctab, in order of appearance
        for k in range(job.nlc[j]):
            lp = job.lcp[j * _FRAG_NLC + k]
            ln = job.lcn[j * _FRAG_NLC + k]
            h = _hash_bytes(lp, ln)
            ci = ctab.find(lp, ln, h)
            if ci < 0:
                ci = len(chroms)
                chroms.append(PyBytes_FromStringAndSize(
                    cython.cast(cython.p_char, lp), ln))
                ctab.add(lp, ln, h, ci)
            lmap[k] = ci
        # its barcodes not in btab, in file order
        if job.nmiss[j] > 0:
            for r in range(b0, b0 + n):
                if job.bc[r] < 0:
                    lp = job.bptr[r]
                    ln = job.blen[r]
                    h = _hash_bytes(lp, ln)
                    bi = btab.find(lp, ln, h)
                    if bi < 0:
                        bi = petrack.barcode_id(PyBytes_FromStringAndSize(
                            cython.cast(cython.p_char, lp), ln))
                        btab.add(lp, ln, h, bi)
                    job.bc[r] = bi
        if n == 0:
            return 0
        nr = job.nrun[j]
        if nr <= _FRAG_NRUN:
            # each run of one chromosome in one call, after the fragments
            # before it
            batch.flush(petrack, chroms)
            a = b0
            for k in range(nr):
                ci = job.runc[j * _FRAG_NRUN + k]
                if ci < 0:
                    ci = lmap[-1 - ci]
                e = job.runend[j * _FRAG_NRUN + k]
                petrack.add_loc_arrays(chroms[ci], self.left_a[a:e],
                                       self.right_a[a:e], self.count_a[a:e],
                                       self.bc_a[a:e])
                a = e
        else:
            # many runs: through the batch, which groups them
            for r in range(b0, b0 + n):
                ci = job.chrom[r]
                if ci < 0:
                    ci = lmap[-1 - ci]
                k = batch.n
                batch.chrom[k] = ci
                batch.left[k] = job.left[r]
                batch.right[k] = job.right[r]
                batch.count[k] = job.count[r]
                batch.bc[k] = job.bc[r]
                batch.n = k + 1
                if batch.n == batch.cap:
                    batch.flush(petrack, chroms)
        return n


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

    @cython.cfunc
    def add_line(self, add_loc, thisline: bytes, max_count,
                 flen: cython.pointer(cython.long)) -> cython.bint:
        """Add the fragment of one line with ``add_loc`` as
        ``build_petrack`` does line by line: returns whether it was
        added, with its length in ``flen[0]``, or raises."""
        chromosome: bytes
        left_pos: cython.int
        right_pos: cython.int
        barcode: bytes
        count: cython.ushort

        (chromosome, left_pos, right_pos, barcode, count) = self.pe_parse_line(thisline)
        if left_pos < 0 or not chromosome:
            return False
        assert right_pos > left_pos, "Right position must be larger than left position, check your BED file at line: %s" % thisline
        flen[0] = right_pos - left_pos
        if max_count:
            count = min(count, max_count)
        add_loc(chromosome, left_pos, right_pos, barcode, count)
        return True

    @cython.cfunc
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    def load_petrack(self, petrack, max_count,
                     msum: cython.pointer(cython.long)) -> cython.long:
        """Add the fragment of every line after the comment lines to
        ``petrack``; returns how many were added, and the sum of their
        lengths in ``msum[0]``.

        The result is that of parsing each line with ``pe_parse_line``,
        skipping it when its chromosome is empty or its left end
        negative, capping its count at ``max_count`` when that is set,
        and calling ``petrack.add_loc``, as earlier versions did line by
        line. The file is decompressed in large blocks (``_BAMStream``,
        on ``_inflate_threads()`` threads when the file starts with a
        BGZF block, as tabix-indexed fragment files do, and serially
        otherwise; the blocks are the serial path's either way), whose
        lines are parsed in C; the fragments of a block
        are appended with ``petrack.add_loc_arrays`` and their barcodes
        given ids with ``petrack.barcode_id``, in file order. While the
        inflate threads run, a block's lines are cut into chunks parsed
        on those threads (``_frag_chunk``), and each chunk's fragments
        are then appended, and its new chromosomes and barcodes given
        ids, on this thread in file order, so the track is the same.
        A line this
        does not parse in C, any line when ``petrack`` is not a
        ``PETrackII`` with a positive buffer size, and any line when
        ``max_count`` is not None or a non-negative int, is handled by
        ``add_line`` in its place in the file.
        """
        offset: cython.Py_ssize_t
        stream: _BAMStream
        fast: cython.bint
        capped: cython.bint = False
        cap: cython.long = 0
        total: cython.long = 0
        m: cython.long = 0
        flen: cython.long = 0
        nextlog: cython.long = 1000000
        need: cython.Py_ssize_t = 1
        more: cython.bint
        p: cython.p_uchar
        q: cython.p_uchar
        nl: cython.p_uchar
        stop: cython.p_uchar
        e: cython.p_uchar
        t: cython.p_uchar
        tabs: cython.p_uchar[5]
        nt: cython.int
        f4end: cython.p_uchar
        left: cython.int
        right: cython.int
        count: cython.int
        clen: cython.Py_ssize_t
        bstart: cython.p_uchar
        blen: cython.Py_ssize_t
        ok: cython.bint
        words: cython.bint
        mark: cython.ulonglong
        k0: cython.ulonglong
        k1: cython.ulonglong
        k2: cython.ulonglong
        h: cython.ulonglong
        ci: cython.int
        bi: cython.int
        lastc: cython.int = -1
        lastp: cython.p_uchar = cython.NULL
        lastn: cython.Py_ssize_t = 8
        lastw: cython.ulonglong = 0
        lastm: cython.ulonglong = 0
        cw: cython.ulonglong
        csame: cython.bint
        s: cython.p_uchar
        gl: cython.Py_ssize_t = 8      # guessed lengths of the two ends
        gr: cython.Py_ssize_t = 8
        gb: cython.Py_ssize_t = 18     # and of the barcode, with its
        bm0: cython.ulonglong = _low_bytes(18)   # three word masks
        bm1: cython.ulonglong = _low_bytes(10)
        bm2: cython.ulonglong = _low_bytes(2)
        k: cython.Py_ssize_t
        chroms: list = []
        ctab: _Interner
        btab: _Interner
        batch: _FragBatch
        chunks: _FragChunks = None
        tasks: cython.Py_ssize_t = 0
        nch: cython.Py_ssize_t = 0     # chunks of the window, the next one
        jc: cython.Py_ssize_t = 0
        segend: cython.p_uchar

        # as line by line: a track without add_loc fails before any read
        add_loc = petrack.add_loc
        fast = (type(petrack) is PETrackII and petrack.buffer_size > 0 and
                (max_count is None or
                 (type(max_count) is int and max_count >= 0)))
        if fast and max_count is not None and max_count > 0:
            capped = True
            cap = min(max_count, 65536)
        ctab = _Interner()
        btab = _Interner()
        batch = _FragBatch(_FRAG_BATCH)

        # the stream starts where skip_first_commentlines left self.fhd;
        # a BGZF file is inflated on several threads, in window mode
        offset = self.fhd.tell()
        self.fhd.close()
        stream = _BAMStream(self.filename, self.gzipped,
                            self.gzipped and _starts_with_bgzf(self.filename),
                            True)
        stream.skip(offset)
        if fast and stream.pool != cython.NULL:
            tasks = _inflate_threads() * _FRAG_TASKS_PER_THREAD
            chunks = _FragChunks(tasks)
            chunks.job.capped = capped
            chunks.job.cap = cython.cast(cython.int, cap)

        while True:
            if jc < nch:
                # the next chunk the threads parsed: its fragments, then
                # its lines from the first one they did not take
                total += chunks.emit(jc, petrack, chroms, ctab, btab, batch)
                m += chunks.job.msum[jc]
                p = chunks.job.fail[jc]
                segend = chunks.job.cuts[jc + 1]
                jc += 1
            else:
                more = stream.refill(need)
                p = stream.buf + stream.start
                q = stream.buf + stream.end
                # the lines are those before stop, one past the last newline
                # in the window; every scan of them ends at a newline
                stop = p
                if p < q and memchr(p, 10, q - p) != cython.NULL:
                    stop = q
                    while (stop - 1)[0] != 10:
                        stop -= 1
                segend = stop
                nch = 0
                jc = 0
                if (tasks > 0 and stream.pool != cython.NULL and
                        stop - p >= _FRAG_PAR_MIN):
                    # the window's lines in chunks, parsed on the inflate
                    # pool's threads with the tables as they are now
                    chunks.cut(p, stop, tasks)
                    _itab_of(cython.address(chunks.job.ctab), ctab)
                    _itab_of(cython.address(chunks.job.btab), btab)
                    with cython.nogil:
                        bgzf_pool_run(stream.pool, _frag_chunk,
                                      cython.address(chunks.job), tasks)
                    nch = tasks
                    continue
            while p < segend:
                ok = False
                words = False
                csame = False
                if fast:
                    # The usual line in one pass: a chromosome that is not
                    # empty, a tab, one to nine digits, a tab, one to nine
                    # digits that are more than the first, a tab, a
                    # barcode, a tab and one to nine digits. Its fields are
                    # the ones the general walk below finds (the digits
                    # after the fourth tab put it before the rstripped
                    # end, and atoi stops where the digits do), so a line
                    # it takes gives the same fragment.
                    t = p
                    if (_LITTLE_ENDIAN and lastn < 8 and p + 8 <= stop and
                            p[lastn] == 9):
                        # the last line's chromosome, a name of up to seven
                        # bytes, followed by a tab
                        memcpy(cython.address(cw), p, 8)
                        if (cw & lastm) == lastw:
                            t = p + lastn
                            csame = True
                    if not csame:
                        while t[0] != 9 and t[0] != 10:
                            t += 1
                    clen = t - p
                    if t[0] == 9 and clen > 0:
                        # Each end is first read as being as long as the
                        # same end of the line before, which in a sorted
                        # file it nearly always is, and only scanned for
                        # its end when it is not.
                        s = t + 1
                        t = cython.NULL
                        if _LITTLE_ENDIAN and s + 10 <= stop:
                            t = _frag_digits_tab(s, gl, cython.address(left))
                            if t == cython.NULL:
                                t = _frag_digits_w(s, cython.address(left))
                                if t != cython.NULL:
                                    gl = t - s
                        else:
                            t = _frag_digits(s, cython.address(left))
                        if t != cython.NULL and t[0] == 9:
                            s = t + 1
                            t = cython.NULL
                            if _LITTLE_ENDIAN and s + 10 <= stop:
                                t = _frag_digits_tab(s, gr, cython.address(right))
                                if t == cython.NULL:
                                    t = _frag_digits_w(s, cython.address(right))
                                    if t != cython.NULL:
                                        gr = t - s
                            else:
                                t = _frag_digits(s, cython.address(right))
                            if t != cython.NULL and t[0] == 9 and right > left:
                                bstart = t + 1
                                t = bstart
                                if (_LITTLE_ENDIAN and bstart + 24 <= stop and
                                        bstart[gb] == 9):
                                    # the barcode as long as the last one:
                                    # no control byte (so no tab or
                                    # newline) before its end
                                    memcpy(cython.address(k0), bstart, 8)
                                    memcpy(cython.address(k1), bstart + 8, 8)
                                    memcpy(cython.address(k2), bstart + 16, 8)
                                    if (((_below_space(k0) & bm0) |
                                         (_below_space(k1) & bm1) |
                                         (_below_space(k2) & bm2)) == 0):
                                        k0 &= bm0
                                        k1 &= bm1
                                        k2 &= bm2
                                        blen = gb
                                        words = True
                                if not words and _LITTLE_ENDIAN and bstart + 24 <= stop:
                                    # the first tab or newline in the next
                                    # 24 bytes, a word at a time, and the
                                    # barcode's words up to it
                                    memcpy(cython.address(k0), bstart, 8)
                                    mark = _tab_or_newline(k0)
                                    if mark != 0:
                                        blen = _lowest_mark(mark)
                                        k0 &= (cython.cast(cython.ulonglong, 1) << (8 * blen)) - 1
                                        k1 = 0
                                        k2 = 0
                                        words = True
                                    else:
                                        memcpy(cython.address(k1), bstart + 8, 8)
                                        mark = _tab_or_newline(k1)
                                        if mark != 0:
                                            blen = _lowest_mark(mark)
                                            k1 &= (cython.cast(cython.ulonglong, 1) << (8 * blen)) - 1
                                            k2 = 0
                                            blen += 8
                                            words = True
                                        else:
                                            memcpy(cython.address(k2), bstart + 16, 8)
                                            mark = _tab_or_newline(k2)
                                            if mark != 0:
                                                blen = _lowest_mark(mark)
                                                k2 &= (cython.cast(cython.ulonglong, 1) << (8 * blen)) - 1
                                                blen += 16
                                                words = True
                                            else:
                                                t = bstart + 24
                                    if words:
                                        # the next barcode is guessed as long
                                        gb = blen
                                        bm0 = _low_bytes(blen)
                                        bm1 = _low_bytes(blen - 8)
                                        bm2 = _low_bytes(blen - 16)
                                if words:
                                    t = bstart + blen
                                while t[0] != 9 and t[0] != 10:
                                    t += 1
                                if t[0] == 9:
                                    blen = t - bstart
                                    nl = _frag_digits(t + 1, cython.address(count))
                                    if nl != cython.NULL:
                                        while nl[0] != 10:
                                            nl += 1
                                        ok = True
                if not ok:
                    words = False
                    csame = False
                    nl = cython.cast(cython.p_uchar, memchr(p, 10, stop - p))
                    if fast:
                        # bytes.rstrip(): drop trailing ASCII whitespace
                        e = nl
                        while e > p and ((e - 1)[0] == 32 or 9 <= (e - 1)[0] <= 13):
                            e -= 1
                        # the first five tabs, as bytes.split(b'\t') would cut
                        nt = 0
                        t = p
                        while t < e:
                            if t[0] == 9:
                                tabs[nt] = t
                                nt += 1
                                if nt == 5:
                                    break
                            t += 1
                        if nt >= 4:
                            f4end = tabs[4] if nt == 5 else e
                            clen = tabs[0] - p
                            left = _frag_atoi(tabs[0] + 1, tabs[1])
                            if clen == 0 or left < 0:
                                # skipped, as add_line skips it
                                p = nl + 1
                                continue
                            right = _frag_atoi(tabs[1] + 1, tabs[2])
                            if right > left:
                                count = _frag_atoi(tabs[3] + 1, f4end)
                                bstart = tabs[2] + 1
                                blen = tabs[3] - bstart
                                ok = True
                if ok:
                    # pe_parse_line's C int to unsigned short
                    count = cython.cast(cython.ushort, count)
                    if capped and cap < count:
                        count = cap
                    if csame or (lastc >= 0 and lastn == clen and
                                 _same_bytes(lastp, p, clen)):
                        ci = lastc
                    else:
                        h = _hash_bytes(p, clen)
                        ci = ctab.find(p, clen, h)
                        if ci < 0:
                            ci = len(chroms)
                            chroms.append(PyBytes_FromStringAndSize(
                                cython.cast(cython.p_char, p), clen))
                            ctab.add(p, clen, h, ci)
                        lastc = ci
                        # the name's bytes, held by chroms
                        lastp = cython.cast(cython.p_uchar,
                                            PyBytes_AS_STRING(chroms[ci]))
                        lastn = clen
                        # and, when shorter than eight bytes, as a word
                        if clen < 8:
                            lastw = 0
                            memcpy(cython.address(lastw), lastp, clen)
                            lastm = _low_bytes(clen)
                    if words:
                        h = _hash_words(k0, k1, k2, blen)
                        bi = btab.find_words(k0, k1, k2, blen, h)
                    else:
                        h = _hash_bytes(bstart, blen)
                        bi = btab.find(bstart, blen, h)
                    if bi < 0:
                        bi = petrack.barcode_id(PyBytes_FromStringAndSize(
                            cython.cast(cython.p_char, bstart), blen))
                        btab.add(bstart, blen, h, bi)
                    k = batch.n
                    batch.chrom[k] = ci
                    batch.left[k] = left
                    batch.right[k] = right
                    batch.count[k] = count
                    batch.bc[k] = bi
                    batch.n = k + 1
                    total += 1
                    m += right - left
                    if batch.n == batch.cap:
                        batch.flush(petrack, chroms)
                    p = nl + 1
                    continue
                # this line, in its place in the file, as line by line
                batch.flush(petrack, chroms)
                while total >= nextlog:
                    info(" %d fragments parsed" % nextlog)
                    nextlog += 1000000
                if self.add_line(add_loc, PyBytes_FromStringAndSize(
                        cython.cast(cython.p_char, p), nl - p), max_count,
                        cython.address(flen)):
                    total += 1
                    m += flen
                p = nl + 1
            if jc < nch:
                continue
            stream.start = p - stream.buf
            batch.flush(petrack, chroms)
            while total >= nextlog:
                info(" %d fragments parsed" % nextlog)
                nextlog += 1000000
            if not more:
                break
            # the unconsumed bytes hold no newline: read on past them
            need = stream.end - stream.start + 1

        if stream.end > stream.start:
            # The last line has no newline. The line-by-line reader this
            # replaces never returned on such a file: its loop re-split
            # that line forever without parsing it. Kept as it was.
            while True:
                PyErr_CheckSignals()
        stream.close()
        msum[0] = m
        return total

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
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths

        petrack = PETrackII(buffer_size=self.buffer_size)
        i = self.load_petrack(petrack, max_count, cython.address(m))

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
        i: cython.long = 0          # number of fragments
        m: cython.long = 0          # sum of fragment lengths

        i = self.load_petrack(petrack, max_count, cython.address(m))

        self.d = (self.d * self.n + m) / (self.n + i)
        self.n += i

        assert self.d >= 0, "Something went wrong (mean fragment size was negative)"
        self.close()
        petrack.set_rlengths({"DUMMYCHROM": 0})
        return petrack
