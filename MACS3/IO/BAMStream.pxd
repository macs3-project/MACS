# cython: language_level=3
# Internal Cython contract; not a supported public extension API.
from MACS3.IO.czlib cimport z_stream
from MACS3.IO.cbgzf_mt cimport bgzf_job_t, bgzf_pool_t, bgzf_task_fn


cdef struct _BAMBatch:
    unsigned char *buf  # owned output, following _BAM_PAR_HEAD reserve
    Py_ssize_t cap
    unsigned char *inb  # borrowed from inbufs; pinned by inviews
    Py_ssize_t in_len
    bgzf_job_t *jobs  # owned; input/output borrowed until pool completion
    Py_ssize_t njobs
    Py_ssize_t out_len


cdef bint _starts_with_bgzf(str filename) except? -1


cdef class _BAMStream:
    cdef object fh
    cdef bytearray inbuf
    cdef unsigned char *inptr
    cdef z_stream strm
    cdef bint gzipped
    cdef bint zready
    cdef bint in_member
    cdef bint eof
    cdef unsigned char *buf
    cdef Py_ssize_t cap
    cdef Py_ssize_t start
    cdef Py_ssize_t end
    # the parallel path; buf is batches[cur].buf while par is set
    cdef bint par
    cdef bgzf_pool_t *pool
    cdef _BAMBatch *batches
    cdef int cur
    cdef bint inflight  # batches[1 - cur] is being inflated
    cdef bint par_eof  # the parallel path has read the whole file
    cdef list inbufs  # the bytearrays behind batches[i].inb
    cdef list inviews
    cdef unsigned char *carry  # compressed bytes the last batch did not take
    cdef Py_ssize_t carry_len
    cdef bint staged  # par_refill has read the next batch's input
    cdef bytearray joined  # the serial path's first input, when two parts
    # window mode: buf is pbuf, or points into batches[cur]
    cdef bint window
    cdef unsigned char *pbuf  # the stream's own buffer
    cdef Py_ssize_t pcap
    cdef unsigned char *src_lo  # the bytes of batches[cur]: src_lo[0:],
    cdef unsigned char *src_pos  # those after buf[0:end] from src_pos,
    cdef unsigned char *src_end  # up to src_end
    cdef bint win_stop  # after src_end, the serial path reads on
    cdef unsigned char *stop_src  # from stop_src[0:stop_n], then
    cdef Py_ssize_t stop_n  # stop_src2[0:stop_n2]
    cdef unsigned char *stop_src2
    cdef Py_ssize_t stop_n2

    cdef object stop_threads(self)
    cpdef close(self)
    cdef bint has_pool(self) noexcept
    cdef void run_tasks(self, bgzf_task_fn fn, void *arg, size_t ntasks) noexcept nogil
    cdef object par_read(self, int t)
    cdef bint par_submit(self) except *
    cdef object to_serial(self, unsigned char *src, Py_ssize_t n,
                          unsigned char *src2, Py_ssize_t n2)
    cdef bint par_refill(self, Py_ssize_t need) except *
    cdef bint win_refill(self, Py_ssize_t need) except *
    cdef bint win_fill(self) except *
    cdef bint refill(self, Py_ssize_t need) except *
    cdef object skip(self, Py_ssize_t n)
