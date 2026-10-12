/* Inflate the blocks of a BGZF file (BAM) on several threads, and run
 * other work on the same threads (bgzf_pool_run).
 *
 * Every BGZF block is a complete gzip member whose total size is stored
 * in its header (the BSIZE field of the "BC" extra subfield) and whose
 * uncompressed size is its last four bytes (ISIZE), so blocks can be
 * located without inflating them and inflated independently, each into
 * its own place in the output.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the BSD License (see the file LICENSE included with
 * the distribution).
 */
#ifndef MACS3_BGZF_MT_H
#define MACS3_BGZF_MT_H

#include <stddef.h>

/* One block to inflate: src_len compressed bytes at src, inflated into
 * exactly dst_len bytes at dst. ok is set to 1 when the gzip member's
 * flags are FEXTRA only and its DEFLATE data, from the end of the extra
 * field to the 8-byte trailer, inflated to exactly dst_len bytes whose
 * CRC32 is the trailer's, and to 0 otherwise. dst is written only when
 * ok is 1. */
typedef struct {
	const unsigned char *src;
	unsigned char *dst;
	unsigned int src_len;
	unsigned int dst_len;
	int ok;
} bgzf_job_t;

/* Why bgzf_split stopped. */
#define BGZF_STOP_MORE 0	/* the next block is not complete in the input */
#define BGZF_STOP_FULL 1	/* the next block does not fit in the output, or max_jobs */
#define BGZF_STOP_OTHER 2	/* the next bytes are not a BGZF block header */

/* The largest ISIZE a BGZF block may have. */
#define BGZF_MAX_ISIZE 65536

/* Split in[0:in_len] into consecutive BGZF blocks, laying their outputs
 * end to end from out, for as long as each block is complete in the input
 * and its output fits in out_cap bytes. Fills jobs[0:n] and returns n;
 * *in_used and *out_used are the bytes the n blocks take, and *reason
 * says why the next block was not taken. Nothing is inflated. */
size_t bgzf_split(const unsigned char *in, size_t in_len,
		  unsigned char *out, size_t out_cap,
		  bgzf_job_t *jobs, size_t max_jobs,
		  size_t *in_used, size_t *out_used, int *reason);

typedef struct bgzf_pool_s bgzf_pool_t;

/* A pool of nthreads - 1 worker threads, the caller of bgzf_pool_wait
 * being the nthreads-th. Returns NULL if nthreads < 2 or if a
 * decompressor or a thread could not be set up. */
bgzf_pool_t *bgzf_pool_new(int nthreads);

/* Start inflating jobs[0:njobs] on the workers and return. jobs and the
 * memory they point to must stay untouched until bgzf_pool_wait returns.
 * Call only when no batch is in progress. */
void bgzf_pool_submit(bgzf_pool_t *pool, bgzf_job_t *jobs, size_t njobs);

/* Inflate jobs of the current batch on the calling thread too, then wait
 * until every job of the batch is done. */
void bgzf_pool_wait(bgzf_pool_t *pool);

/* Stop and join every worker thread and free the pool. Jobs of a batch in
 * progress may be left undone, but no thread touches them after this. */
void bgzf_pool_free(bgzf_pool_t *pool);

/* A task of bgzf_pool_run: called once for each i of a batch. */
typedef void (*bgzf_task_fn)(void *arg, size_t i);

/* Call fn(arg, i) for every i in [0, ntasks) on the pool's threads and
 * the calling thread, and return when every call has returned. The
 * calls run ahead of the jobs of an inflate batch in progress: a worker
 * inflating a block finishes that block first. They must not depend on
 * one another or on the order they run in, and must not use the pool.
 * Call only from the thread that submits and waits for inflate batches. */
void bgzf_pool_run(bgzf_pool_t *pool, bgzf_task_fn fn, void *arg,
		   size_t ntasks);

#endif
