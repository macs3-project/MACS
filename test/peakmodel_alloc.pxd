from libc.stddef cimport size_t

cdef extern from "peakmodel_alloc.h" nogil:
    void *pm_realloc(void *ptr, size_t size) noexcept
    void pm_free(void *ptr) noexcept
    void pm_reset(size_t target, int worker_only) noexcept
    void pm_stats(size_t *out) noexcept
