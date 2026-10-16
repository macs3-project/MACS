# Compile-time FRAG layout bounds. Keep arrays, allocation strides and
# comparison limits coupled; none of these is a Python runtime setting.
cdef enum:
    FRAG_LOCAL_CHROMOSOMES = 8
    FRAG_LOCAL_RUNS = 8
    INTERNER_INLINE_KEY_BYTES = 36
    INTERNER_ALIGNMENT_BYTES = 64
