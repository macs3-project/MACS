# cython: language_level=3

# BGZF inflation and native tasks on the stream's budgeted pool.
# bgzf_mt.c is compiled only into BAMStream; Parser submits its FRAG
# callback through the stream's Cython interface, without owning the pool.

cdef extern from "bgzf_mt.h" nogil:
    ctypedef struct bgzf_job_t:
        const unsigned char *src
        unsigned char *dst
        unsigned int src_len
        unsigned int dst_len
        int ok

    ctypedef struct bgzf_pool_t:
        pass

    int BGZF_STOP_MORE
    int BGZF_STOP_FULL
    int BGZF_STOP_OTHER
    int BGZF_MAX_ISIZE

    size_t bgzf_split(const unsigned char *src, size_t in_len,
                      unsigned char *out, size_t out_cap,
                      bgzf_job_t *jobs, size_t max_jobs,
                      size_t *in_used, size_t *out_used, int *reason)
    bgzf_pool_t *bgzf_pool_new(int nthreads)
    void bgzf_pool_submit(bgzf_pool_t *pool, bgzf_job_t *jobs, size_t njobs)
    void bgzf_pool_wait(bgzf_pool_t *pool)
    void bgzf_pool_free(bgzf_pool_t *pool)

    ctypedef void (*bgzf_task_fn)(void *arg, size_t i) noexcept nogil
    void bgzf_pool_run(bgzf_pool_t *pool, bgzf_task_fn fn, void *arg,
                       size_t ntasks)
