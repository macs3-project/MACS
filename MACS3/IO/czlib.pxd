# cython: language_level=3

# The part of zlib's inflate API used by MACS3.IO.Parser to decompress
# BAM (BGZF) files in large blocks. Link the extension with libz.

cdef extern from "zlib.h" nogil:
    ctypedef struct z_stream:
        unsigned char *next_in
        unsigned int avail_in
        unsigned long total_in
        unsigned char *next_out
        unsigned int avail_out
        unsigned long total_out
        char *msg

    int Z_OK
    int Z_STREAM_END
    int Z_NEED_DICT
    int Z_BUF_ERROR
    int Z_NO_FLUSH

    int inflateInit2(z_stream *strm, int windowBits)
    int inflate(z_stream *strm, int flush)
    int inflateReset(z_stream *strm)
    int inflateEnd(z_stream *strm)
