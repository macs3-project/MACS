# cython: language_level=3
"""Internal BAM/BGZF byte streaming and native-pool ownership.

Parser retains format decoding, FRAG interning and track assembly. This
module depends only on native IO declarations and Performance defaults;
it must not import Parser or tracks. Its Cython interface is internal.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import gzip
import io

import cython
from cython.cimports.libc.stdlib import malloc, calloc, realloc, free
from cython.cimports.libc.string import memcpy, memmove, memset
from cython.cimports.cpython.bytearray import PyByteArray_AS_STRING
from cython.cimports.MACS3.IO.czlib import (
    z_stream, inflateInit2, inflate, inflateReset, inflateEnd,
    Z_OK, Z_STREAM_END, Z_BUF_ERROR, Z_NO_FLUSH)
from cython.cimports.MACS3.IO.cbgzf_mt import (
    bgzf_job_t, bgzf_split, bgzf_pool_new, bgzf_pool_submit,
    bgzf_pool_wait, bgzf_pool_free, bgzf_pool_run, bgzf_task_fn)
from MACS3.Utilities.Performance import DEFAULTS


# Compressed BAM input is read _BAM_IN_SIZE bytes at a time and inflated
# into a buffer of _BAM_OUT_SIZE bytes (grown if one record is larger),
# whose records are then walked in C.
_BAM_IN_SIZE = cython.declare(cython.Py_ssize_t, DEFAULTS.bam_serial_input_bytes)
_BAM_OUT_SIZE = cython.declare(cython.Py_ssize_t, DEFAULTS.bam_serial_output_bytes)

# BGZF inflation uses the invocation thread budget: _BAM_PAR_IN compressed bytes at a time
# are split into blocks, which are inflated into _BAM_PAR_OUT bytes
# following _BAM_PAR_HEAD bytes kept for the unconsumed tail of the
# previous batch.
_BAM_PAR_IN = cython.declare(cython.Py_ssize_t, DEFAULTS.bgzf_batch_input_bytes)
_BAM_PAR_OUT = cython.declare(cython.Py_ssize_t, DEFAULTS.bgzf_batch_output_bytes)
_BAM_PAR_HEAD = cython.declare(cython.Py_ssize_t, DEFAULTS.bgzf_tail_reserve_bytes)
_BAM_PAR_MAXJOBS = cython.declare(cython.Py_ssize_t, DEFAULTS.bgzf_max_batch_jobs)


@cython.cfunc
def _starts_with_bgzf(filename: str) -> cython.bint:
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
    ISIZE, on the invocation's budgeted threads, the next batch
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
    records there, advance ``start`` and call ``refill`` for more. The
    pointer is borrowed until the next refill, fallback or destruction;
    refill may move, resize or replace it. False at EOF may leave a short
    unread tail for the parser to validate.

    The stream owns the file, zlib state, pool, batch/job allocations and
    Python input buffers. Jobs borrow batch storage until completion;
    carry, staged input and zlib pointers borrow the retained bytearrays
    or joined fallback input. Normal parallel output aliases a batch;
    window output aliases a batch slice or the owned pbuf. Serial fallback
    transfers the active output allocation and clears its former owning
    alias, while retaining compressed input owners through decompression.

    Failed pool creation uses the serial path. Partial construction is
    safe to destroy. close joins pending inflation and closes the file;
    destruction then releases zlib and native buffers. Neither may free
    storage while workers still borrow it. C fields/methods are declared
    in BAMStream.pxd solely for internal compiled consumers.
    """
    def __cinit__(self):
        self.buf = cython.NULL
        self.zready = False
        self.par = False
        self.pool = cython.NULL
        self.batches = cython.NULL
        self.inflight = False
        self.window = False

    def __init__(self, filename: str, gzipped: cython.bint,
                 parallel: cython.bint = True, window: cython.bint = False,
                 threads: cython.int = 1):
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
            nthreads = threads if parallel else 1
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
    @cython.exceptval(check=False)
    def has_pool(self) -> cython.bint:
        """Whether native tasks can use the stream's pool, before close."""
        return self.pool != cython.NULL

    @cython.cfunc
    @cython.nogil
    @cython.exceptval(check=False)
    def run_tasks(self, fn: bgzf_task_fn, arg: cython.p_void,
                  ntasks: cython.size_t) -> cython.void:
        """Run native tasks on the live pool, synchronously.

        Requires has_pool(). fn and arg must remain valid until return;
        callbacks run without the GIL and may overlap pending inflation.
        The caller must keep task output disjoint from inflate output and
        must not refill or close this stream until all tasks return.
        """
        bgzf_pool_run(self.pool, fn, arg, ntasks)

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
