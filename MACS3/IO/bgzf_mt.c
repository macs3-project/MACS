/* Inflate the blocks of a BGZF file (BAM) on several threads; see
 * bgzf_mt.h.
 *
 * Each block is inflated by libdeflate (vendored in libdeflate/: release
 * v1.26, commit 92e6a0db9fa848d742f9eb286c92afc60f2c3dda of
 * https://github.com/ebiggers/libdeflate, the files decompression and
 * CRC32 need, unmodified; MIT licence in libdeflate/COPYING), which
 * picks its instructions for this CPU at run time. A block is ok only
 * when its gzip flags are exactly FEXTRA, as BGZF writers set them, its
 * DEFLATE data inflates to exactly ISIZE bytes and ends exactly at the
 * 8-byte trailer, and the CRC32 of those bytes is the trailer's. Any
 * other block is marked not ok; the caller then reads from that block on
 * serially with zlib, which reproduces the serial reader's behaviour and
 * errors.
 *
 * Each thread inflates a block into a block-sized buffer of its own and
 * then copies it to its place in the output. Measured on a 2-socket EPYC
 * (medians of 10 runs), a thread inflating straight into the output ran
 * at 1.56-1.70 GB/s on the CCD of the thread that walks the output,
 * 1.24-1.44 GB/s on another CCD of the same socket and 0.98-1.12 GB/s on
 * the other socket; through the buffer, 1.54-1.68 GB/s on all three. The
 * cause is inferred, not measured, since perf counters were unavailable:
 * the stores go into cache lines the walking thread's cache holds, and
 * libdeflate's stores, interleaved with decoding, would wait for each
 * line in turn, where the copy fetches many at once.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the BSD License (see the file LICENSE included with
 * the distribution).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>

#include "bgzf_mt.h"
#include "libdeflate/libdeflate.h"

static unsigned int le16(const unsigned char *p)
{
	return (unsigned int)p[0] | ((unsigned int)p[1] << 8);
}

static uint32_t le32(const unsigned char *p)
{
	return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
		((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

size_t bgzf_split(const unsigned char *in, size_t in_len,
		  unsigned char *out, size_t out_cap,
		  bgzf_job_t *jobs, size_t max_jobs,
		  size_t *in_used, size_t *out_used, int *reason)
{
	size_t pos = 0, opos = 0, n = 0;
	size_t rem, xend, q, bsize, slen;
	uint32_t isize;
	const unsigned char *b;

	*reason = BGZF_STOP_FULL;
	while (n < max_jobs) {
		rem = in_len - pos;
		b = in + pos;
		if (rem < 12) {
			*reason = BGZF_STOP_MORE;
			break;
		}
		/* gzip magic, deflate, FEXTRA */
		if (b[0] != 31 || b[1] != 139 || b[2] != 8 || !(b[3] & 4)) {
			*reason = BGZF_STOP_OTHER;
			break;
		}
		xend = 12 + (size_t)le16(b + 10);
		if (rem < xend) {
			*reason = BGZF_STOP_MORE;
			break;
		}
		/* the "BC" subfield holds BSIZE, the block size minus 1 */
		bsize = 0;
		q = 12;
		while (q + 4 <= xend) {
			slen = le16(b + q + 2);
			if (b[q] == 66 && b[q + 1] == 67 && slen == 2) {
				if (q + 6 <= xend)
					bsize = (size_t)le16(b + q + 4) + 1;
				break;
			}
			q += 4 + slen;
		}
		if (bsize < xend + 8) {
			/* no BC subfield, or too small to hold the trailer */
			*reason = BGZF_STOP_OTHER;
			break;
		}
		if (rem < bsize) {
			*reason = BGZF_STOP_MORE;
			break;
		}
		isize = le32(b + bsize - 4);
		if (isize > BGZF_MAX_ISIZE) {
			*reason = BGZF_STOP_OTHER;
			break;
		}
		if (out_cap - opos < isize) {
			*reason = BGZF_STOP_FULL;
			break;
		}
		jobs[n].src = b;
		jobs[n].dst = out + opos;
		jobs[n].src_len = (unsigned int)bsize;
		jobs[n].dst_len = isize;
		jobs[n].ok = 0;
		n++;
		pos += bsize;
		opos += isize;
	}
	*in_used = pos;
	*out_used = opos;
	return n;
}

/* The BGZF_MAX_ISIZE-byte buffer the calling thread inflates into: its
 * worker's, set by bgzf_worker, or the pool's, set by bgzf_pool_wait. */
static __thread unsigned char *bgzf_scratch;

/* Inflate one block, a whole gzip member with flags FEXTRA only, whose
 * DEFLATE data runs from the end of the extra field to the trailer, into
 * bgzf_scratch, and copy it to j->dst if it is ok. */
static void bgzf_inflate_one(struct libdeflate_decompressor *d, bgzf_job_t *j)
{
	const unsigned char *s = j->src;
	unsigned char *out = bgzf_scratch;
	size_t xend, clen, in_used, out_used;

	j->ok = 0;
	if (s[3] != 4)
		return;
	/* bgzf_split took the block only if src_len >= xend + 8 */
	xend = 12 + (size_t)le16(s + 10);
	clen = j->src_len - xend - 8;
	if (libdeflate_deflate_decompress_ex(d, s + xend, clen, out,
					     j->dst_len, &in_used, &out_used)
	    != LIBDEFLATE_SUCCESS)
		return;
	j->ok = in_used == clen && out_used == j->dst_len &&
		libdeflate_crc32(0, out, j->dst_len) ==
		le32(s + j->src_len - 8);
	if (j->ok)
		memcpy(j->dst, out, j->dst_len);
}

typedef struct {
	bgzf_pool_t *pool;
	struct libdeflate_decompressor *d;
	unsigned char *scratch;	/* BGZF_MAX_ISIZE bytes */
	pthread_t tid;
} bgzf_worker_t;

struct bgzf_pool_s {
	pthread_mutex_t mu;
	pthread_cond_t work_cv;	/* a new batch, or quit */
	pthread_cond_t done_cv;	/* every job of the batch is done */
	bgzf_job_t *jobs;	/* the current batch, under mu */
	size_t njobs;
	size_t ndone;
	unsigned int gen;	/* batch number, under mu */
	int quit;
	/* (gen << 32) | index of the next job to take; a thread takes jobs
	 * only of the batch it was woken for */
	uint64_t claim;
	struct libdeflate_decompressor *d;	/* the caller's */
	unsigned char *scratch;	/* the caller's, BGZF_MAX_ISIZE bytes */
	int nworkers;
	bgzf_worker_t *workers;
	/* the current batch of bgzf_pool_run, under mu */
	bgzf_task_fn tfn;
	void *targ;
	size_t ntasks;
	size_t tnext;		/* the next task to take */
	size_t tdone;
	pthread_cond_t tdone_cv;	/* every task of the batch is done */
	/* tnext < ntasks; set and cleared under mu, read without it between
	 * inflate jobs */
	int tpending;
};

/* Take and run tasks of the current task batch until none is left. */
static void bgzf_tasks(bgzf_pool_t *p)
{
	bgzf_task_fn fn;
	void *arg;
	size_t i;

	pthread_mutex_lock(&p->mu);
	while (p->tnext < p->ntasks) {
		i = p->tnext++;
		if (p->tnext == p->ntasks)
			__atomic_store_n(&p->tpending, 0, __ATOMIC_RELAXED);
		fn = p->tfn;
		arg = p->targ;
		pthread_mutex_unlock(&p->mu);
		fn(arg, i);
		pthread_mutex_lock(&p->mu);
		/* the batch cannot change before its tasks are counted */
		if (++p->tdone == p->ntasks)
			pthread_cond_signal(&p->tdone_cv);
	}
	pthread_mutex_unlock(&p->mu);
}

static long bgzf_claim(bgzf_pool_t *p, unsigned int gen, size_t njobs)
{
	uint64_t v = __atomic_load_n(&p->claim, __ATOMIC_ACQUIRE);

	for (;;) {
		if ((unsigned int)(v >> 32) != gen || (v & 0xffffffffu) >= njobs)
			return -1;
		if (__atomic_compare_exchange_n(&p->claim, &v, v + 1, 1,
						__ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE))
			return (long)(v & 0xffffffffu);
	}
}

/* Inflate jobs of batch gen until none is left; return how many.
 * Tasks of bgzf_pool_run waiting to be taken are run first. */
static size_t bgzf_run(bgzf_pool_t *p, struct libdeflate_decompressor *d,
		       unsigned int gen, bgzf_job_t *jobs, size_t njobs)
{
	size_t mine = 0;
	long i;

	for (;;) {
		if (__atomic_load_n(&p->tpending, __ATOMIC_RELAXED))
			bgzf_tasks(p);
		if ((i = bgzf_claim(p, gen, njobs)) < 0)
			break;
		bgzf_inflate_one(d, &jobs[i]);
		mine++;
	}
	return mine;
}

static void *bgzf_worker(void *arg)
{
	bgzf_worker_t *w = (bgzf_worker_t *)arg;
	bgzf_pool_t *p = w->pool;
	unsigned int seen = 0, gen;
	bgzf_job_t *jobs;
	size_t njobs, mine;

	bgzf_scratch = w->scratch;
	pthread_mutex_lock(&p->mu);
	for (;;) {
		while (!p->quit && p->gen == seen && p->tnext >= p->ntasks)
			pthread_cond_wait(&p->work_cv, &p->mu);
		if (p->quit)
			break;
		if (p->tnext < p->ntasks) {
			pthread_mutex_unlock(&p->mu);
			bgzf_tasks(p);
			pthread_mutex_lock(&p->mu);
			continue;
		}
		gen = seen = p->gen;
		jobs = p->jobs;
		njobs = p->njobs;
		pthread_mutex_unlock(&p->mu);
		mine = bgzf_run(p, w->d, gen, jobs, njobs);
		pthread_mutex_lock(&p->mu);
		if (mine > 0) {
			/* the batch cannot change before its jobs are counted */
			p->ndone += mine;
			if (p->ndone == p->njobs)
				pthread_cond_signal(&p->done_cv);
		}
	}
	pthread_mutex_unlock(&p->mu);
	return NULL;
}

bgzf_pool_t *bgzf_pool_new(int nthreads)
{
	bgzf_pool_t *p;
	sigset_t all, old;
	int i, ok = 1;

	if (nthreads < 2)
		return NULL;
	p = (bgzf_pool_t *)calloc(1, sizeof(bgzf_pool_t));
	if (p == NULL)
		return NULL;
	p->workers = (bgzf_worker_t *)calloc(nthreads - 1, sizeof(bgzf_worker_t));
	if (p->workers == NULL) {
		free(p);
		return NULL;
	}
	p->d = libdeflate_alloc_decompressor();
	p->scratch = (unsigned char *)malloc(BGZF_MAX_ISIZE);
	if (p->d == NULL || p->scratch == NULL) {
		if (p->d != NULL)
			libdeflate_free_decompressor(p->d);
		free(p->scratch);
		free(p->workers);
		free(p);
		return NULL;
	}
	pthread_mutex_init(&p->mu, NULL);
	pthread_cond_init(&p->work_cv, NULL);
	pthread_cond_init(&p->done_cv, NULL);
	pthread_cond_init(&p->tdone_cv, NULL);
	/* workers take no signals; the caller's thread handles them */
	sigfillset(&all);
	pthread_sigmask(SIG_SETMASK, &all, &old);
	for (i = 0; i < nthreads - 1; i++) {
		p->workers[i].pool = p;
		p->workers[i].d = libdeflate_alloc_decompressor();
		p->workers[i].scratch = (unsigned char *)malloc(BGZF_MAX_ISIZE);
		if (p->workers[i].d == NULL || p->workers[i].scratch == NULL ||
		    pthread_create(&p->workers[i].tid, NULL, bgzf_worker,
				   &p->workers[i]) != 0) {
			if (p->workers[i].d != NULL)
				libdeflate_free_decompressor(p->workers[i].d);
			free(p->workers[i].scratch);
			ok = 0;
			break;
		}
		p->nworkers++;
	}
	pthread_sigmask(SIG_SETMASK, &old, NULL);
	if (!ok) {
		bgzf_pool_free(p);
		return NULL;
	}
	return p;
}

void bgzf_pool_submit(bgzf_pool_t *p, bgzf_job_t *jobs, size_t njobs)
{
	pthread_mutex_lock(&p->mu);
	p->jobs = jobs;
	p->njobs = njobs;
	p->ndone = 0;
	p->gen++;
	__atomic_store_n(&p->claim, (uint64_t)p->gen << 32, __ATOMIC_RELEASE);
	pthread_cond_broadcast(&p->work_cv);
	pthread_mutex_unlock(&p->mu);
}

void bgzf_pool_wait(bgzf_pool_t *p)
{
	unsigned int gen;
	bgzf_job_t *jobs;
	size_t njobs, mine;

	pthread_mutex_lock(&p->mu);
	gen = p->gen;
	jobs = p->jobs;
	njobs = p->njobs;
	pthread_mutex_unlock(&p->mu);
	bgzf_scratch = p->scratch;
	mine = bgzf_run(p, p->d, gen, jobs, njobs);
	pthread_mutex_lock(&p->mu);
	p->ndone += mine;
	while (p->ndone < p->njobs)
		pthread_cond_wait(&p->done_cv, &p->mu);
	pthread_mutex_unlock(&p->mu);
}

void bgzf_pool_run(bgzf_pool_t *p, bgzf_task_fn fn, void *arg,
		   size_t ntasks)
{
	if (ntasks == 0)
		return;
	pthread_mutex_lock(&p->mu);
	p->tfn = fn;
	p->targ = arg;
	p->ntasks = ntasks;
	p->tnext = 0;
	p->tdone = 0;
	__atomic_store_n(&p->tpending, 1, __ATOMIC_RELAXED);
	pthread_cond_broadcast(&p->work_cv);
	pthread_mutex_unlock(&p->mu);
	bgzf_tasks(p);
	pthread_mutex_lock(&p->mu);
	while (p->tdone < p->ntasks)
		pthread_cond_wait(&p->tdone_cv, &p->mu);
	pthread_mutex_unlock(&p->mu);
}

void bgzf_pool_free(bgzf_pool_t *p)
{
	int i;

	if (p == NULL)
		return;
	pthread_mutex_lock(&p->mu);
	p->quit = 1;
	pthread_cond_broadcast(&p->work_cv);
	pthread_mutex_unlock(&p->mu);
	for (i = 0; i < p->nworkers; i++) {
		pthread_join(p->workers[i].tid, NULL);
		libdeflate_free_decompressor(p->workers[i].d);
		free(p->workers[i].scratch);
	}
	libdeflate_free_decompressor(p->d);
	free(p->scratch);
	pthread_cond_destroy(&p->done_cv);
	pthread_cond_destroy(&p->tdone_cv);
	pthread_cond_destroy(&p->work_cv);
	pthread_mutex_destroy(&p->mu);
	free(p->workers);
	free(p);
}
