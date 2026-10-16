# Compile-time FRAG layout bounds. Keep arrays, allocation strides and
# comparison limits coupled; none of these is a Python runtime setting.
cdef enum:
    FRAG_LOCAL_CHROMOSOMES = 8
    FRAG_LOCAL_RUNS = 8
    INTERNER_INLINE_KEY_BYTES = 36
    INTERNER_ALIGNMENT_BYTES = 64

# Cython 3.0 cannot use imported enum bounds in pure-Python array syntax.
# Declare the bounded arrays here so every supported compiler uses the
# same enums for storage, allocation strides and indexing limits.
ctypedef int FragLocalMap[FRAG_LOCAL_CHROMOSOMES]

cdef struct InternerSlot:
    unsigned long long h
    Py_ssize_t n1
    Py_ssize_t off
    int value
    unsigned char key[INTERNER_INLINE_KEY_BYTES]
