#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "rld0.h"

#define RLD_IBITS_PLUS 4

#define rld_file_size(e) ((4 + (e)->asize) * 8 + (e)->n_bytes + 8 * (e)->n_frames * ((e)->asize + 1))

#ifndef xcalloc
#define xcalloc(n, s) calloc(n, s)
#endif
#ifndef xmalloc
#define xmalloc(s) malloc(s)
#endif

/******************
 * Delta encoding *
 ******************/

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	register uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

static inline int ilog2(uint64_t v)
{
	return v>>32? 32 + ilog2_32(v>>32) : ilog2_32(v);
}

static inline int64_t rld_delta_enc1(uint64_t x, int *width)
{
	int y = ilog2(x);
	int z = ilog2_32(y + 1);
	*width = (z<<1) + 1 + y;
	return (x ^ (uint64_t)1<<y) | (uint64_t)(y+1)<<y;
}

/***********************************
 * Initialization and deallocation *
 ***********************************/

rld_t *rld_init(int asize, int bbits)
{
	rld_t *e;
	e = xcalloc(1, sizeof(rld_t));
	e->n = 1;
	e->z = xmalloc(sizeof(void*));
	e->z[0] = xcalloc(RLD_LSIZE, 8);
	e->ssize = 1<<bbits;
	e->cnt = xcalloc(asize + 1, 8);
	e->mcnt = xcalloc(asize + 1, 8);
	e->abits = ilog2(asize) + 1;
	e->asize = asize;
	e->sbits = bbits;
	e->asize1 = asize + 1;
	e->offset0[0] = (e->asize1*16+63)/64;
	e->offset0[1] = (e->asize1*32+63)/64;
	e->offset0[2] = e->asize1;
	return e;
}

void rld_destroy(rld_t *e)
{
	int i = 0;
	if (e == 0) return;
	if (e->mem) {
		close(e->fd);
		munmap(e->mem, rld_file_size(e));
	} else {
		for (i = 0; i < e->n; ++i) free(e->z[i]);
		free(e->frame);
	}
	free(e->z); free(e->cnt); free(e->mcnt); free(e);
}

void rld_itr_init(const rld_t *e, rlditr_t *itr, uint64_t k)
{
	itr->i = e->z + (k >> RLD_LBITS);
	itr->shead = *itr->i + k%RLD_LSIZE;
	itr->stail = rld_get_stail(e, itr);
	itr->p = itr->shead + e->offset0[rld_block_type(*itr->shead)];
	itr->q = (uint8_t*)itr->p;
	itr->r = 64;
	itr->c = -1;
	itr->l = 0;
}

/************
 * Encoding *
 ************/

static inline void enc_next_block(rld_t *e, rlditr_t *itr)
{
	int i, type;
	if (itr->stail + 2 - *itr->i == RLD_LSIZE) {
		++e->n;
		e->z = realloc(e->z, e->n * sizeof(void*));
		itr->i = e->z + e->n - 1;
		itr->shead = *itr->i = xcalloc(RLD_LSIZE, 8);
	} else itr->shead += e->ssize;
	if (e->cnt[0] - e->mcnt[0] < 0x4000) {
		uint16_t *p = (uint16_t*)itr->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		type = 0;
	} else if (e->cnt[0] - e->mcnt[0] < 0x40000000) {
		uint32_t *p = (uint32_t*)itr->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		type = 1;
	} else {
		uint64_t *p = (uint64_t*)itr->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		type = 2;
	}
	*itr->shead |= (uint64_t)type<<62;
	itr->p = itr->shead + e->offset0[type];
	itr->stail = rld_get_stail(e, itr);
	itr->q = (uint8_t*)itr->p;
	itr->r = 64;
	for (i = 0; i <= e->asize; ++i) e->mcnt[i] = e->cnt[i];
}

static int rld_enc1(rld_t *e, rlditr_t *itr, int64_t l, uint8_t c)
{
	int w;
	uint64_t x = rld_delta_enc1(l, &w) << e->abits | c;
	w += e->abits;
	if (w >= itr->r && itr->p == itr->stail) enc_next_block(e, itr);
	if (w > itr->r) {
		w -= itr->r;
		*itr->p++ |= x >> w;
		*itr->p = x << (itr->r = 64 - w);
	} else itr->r -= w, *itr->p |= x << itr->r;
	e->cnt[0] += l;
	e->cnt[c + 1] += l;
	return 0;
}

int rld_enc(rld_t *e, rlditr_t *itr, int64_t l, uint8_t c)
{
	if (l == 0) return 0;
	if (itr->c != c) {
		if (itr->l) rld_enc1(e, itr, itr->l, itr->c);
		itr->l = l; itr->c = c;
	} else itr->l += l;
	return 0;
}

void rld_rank_index(rld_t *e)
{
	uint64_t last, n_blks, i, k, *cnt;
	int j;

	n_blks = e->n_bytes * 8 / 64 / e->ssize + 1;
	last = rld_last_blk(e);
	cnt = alloca(e->asize * 8);
	e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
	e->n_frames = ((e->mcnt[0] + (1ll<<e->ibits) - 1) >> e->ibits) + 1;
	e->frame = xcalloc(e->n_frames * e->asize1, 8);
	e->frame[0] = 0;
	for (j = 0; j < e->asize; ++j) cnt[j] = 0;
	for (i = e->ssize, k = 1; i <= last; i += e->ssize) {
		uint64_t sum, *p = rld_seek_blk(e, i);
		int type = rld_block_type(*p);
		if (type == 0) {
			uint16_t *q = (uint16_t*)p;
			for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j];
		} else if (type == 1) {
			uint32_t *q = (uint32_t*)p;
			for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j] & 0x3fffffff;
		} else {
			uint64_t *q = (uint64_t*)p;
			for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j];
		}
		for (j = 0, sum = 0; j < e->asize; ++j) sum += cnt[j];
		while (sum >= k<<e->ibits) ++k;
		if (k < e->n_frames) {
			uint64_t x = k * e->asize1;
			e->frame[x] = i;
			for (j = 0; j < e->asize; ++j) e->frame[x + j + 1] = cnt[j];
		}
	}
	assert(k >= e->n_frames - 1);
	for (k = 1; k < e->n_frames; ++k) { // fill zero cells
		uint64_t x = k * e->asize1;
		if (e->frame[x] == 0) {
			for (j = 0; j <= e->asize; ++j)
				e->frame[x + j] = e->frame[x - e->asize1 + j];
		}
	}
}

uint64_t rld_enc_finish(rld_t *e, rlditr_t *itr)
{
	int i;
	if (itr->l) rld_enc1(e, itr, itr->l, itr->c);
	enc_next_block(e, itr);
	e->n_bytes = (((uint64_t)(e->n - 1) * RLD_LSIZE) + (itr->p - *itr->i)) * 8;
	// recompute e->cnt as the accumulative count; e->mcnt[] keeps the marginal counts
	for (e->cnt[0] = 0, i = 1; i <= e->asize; ++i) e->cnt[i] += e->cnt[i - 1];
	rld_rank_index(e);
	return e->n_bytes;
}

/*****************
 * Save and load *
 *****************/

int rld_dump(const rld_t *e, const char *fn)
{
	uint64_t k = 0;
	int i;
	uint32_t a;
	FILE *fp;
	fp = strcmp(fn, "-")? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
	if (fp == 0) return -1;
	a = e->asize<<16 | e->sbits;
	fwrite("RLD\3", 1, 4, fp); // write magic
	fwrite(&a, 4, 1, fp); // write sbits and asize
	fwrite(&k, 8, 1, fp); // preserve 8 bytes for future uses
	fwrite(&e->n_bytes, 8, 1, fp); // n_bytes can always be divided by 8
	fwrite(&e->n_frames, 8, 1, fp); // number of frames
	fwrite(e->mcnt + 1, 8, e->asize, fp); // write the marginal counts
	for (i = 0, k = e->n_bytes / 8; i < e->n - 1; ++i, k -= RLD_LSIZE)
		fwrite(e->z[i], 8, RLD_LSIZE, fp);
	fwrite(e->z[i], 8, k, fp);
	fwrite(e->frame, 8 * e->asize1, e->n_frames, fp);
	fclose(fp);
	return 0;
}

static rld_t *rld_restore_header(const char *fn, FILE **_fp)
{
	FILE *fp;
	rld_t *e;
	char magic[4];
	uint64_t a[3];
	int32_t i, x;

	if (strcmp(fn, "-") == 0) *_fp = fp = stdin;
	else if ((*_fp = fp = fopen(fn, "rb")) == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, "RLD\3", 4)) return 0;
	fread(&x, 4, 1, fp);
	e = rld_init(x>>16, x&0xffff);
	fread(a, 8, 3, fp);
	e->n_bytes = a[1]; e->n_frames = a[2];
	fread(e->mcnt + 1, 8, e->asize, fp);
	for (i = 0; i <= e->asize; ++i) e->cnt[i] = e->mcnt[i];
	for (i = 1; i <= e->asize; ++i) e->cnt[i] += e->cnt[i - 1];
	e->mcnt[0] = e->cnt[e->asize];
	return e;
}

rld_t *rld_restore(const char *fn)
{
	FILE *fp;
	rld_t *e;
	uint64_t k, n_blks;
	int32_t i;

	if ((e = rld_restore_header(fn, &fp)) == 0) { // then load as plain DNA rle
		uint8_t *buf;
		int l;
		rlditr_t itr;
		buf = malloc(0x10000);
		e = rld_init(6, 3);
		rld_itr_init(e, &itr, 0);
		while ((l = fread(buf, 1, 0x10000, fp)) != 0)
			for (i = 0; i < l; ++i)
				if (buf[i]>>3) rld_enc(e, &itr, buf[i]>>3, buf[i]&7);
		fclose(fp);
		free(buf);
		rld_enc_finish(e, &itr);
		return e;
	}
	if (e->n_bytes / 8 > RLD_LSIZE) { // allocate enough memory
		e->n = (e->n_bytes / 8 + RLD_LSIZE - 1) / RLD_LSIZE;
		e->z = realloc(e->z, e->n * sizeof(void*));
		for (i = 1; i < e->n; ++i)
			e->z[i] = xcalloc(RLD_LSIZE, 8);
	}
	for (i = 0, k = e->n_bytes / 8; i < e->n - 1; ++i, k -= RLD_LSIZE)
		fread(e->z[i], 8, RLD_LSIZE, fp);
	fread(e->z[i], 8, k, fp);
	e->frame = xmalloc(e->n_frames * e->asize1 * 8);
	fread(e->frame, 8 * e->asize1, e->n_frames, fp);
	fclose(fp);
	n_blks = e->n_bytes * 8 / 64 / e->ssize + 1;
	e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
	return e;
}

rld_t *rld_restore_mmap(const char *fn)
{
	FILE *fp;
	rld_t *e;
	int i;
	int64_t n_blks;

	e = rld_restore_header(fn, &fp);
	fclose(fp);
	free(e->z[0]); free(e->z);
	e->n = (e->n_bytes / 8 + RLD_LSIZE - 1) / RLD_LSIZE;
	e->z = xcalloc(e->n, sizeof(void*));
	e->fd = open(fn, O_RDONLY);
	e->mem = (uint64_t*)mmap(0, rld_file_size(e), PROT_READ, MAP_PRIVATE, e->fd, 0);
	for (i = 0; i < e->n; ++i) e->z[i] = e->mem + (4 + e->asize) + (size_t)i * RLD_LSIZE;
	e->frame = e->mem + (4 + e->asize) + e->n_bytes/8;
	n_blks = e->n_bytes * 8 / 64 / e->ssize + 1;
	e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
	return e;
}

/******************
 * Computing rank *
 ******************/

#ifdef _DNA_ONLY
static inline int64_t rld_dec0_fast_dna(const rld_t *e, rlditr_t *itr, int *c)
{ // This is NOT a replacement of rld_dec0(). It does not do boundary check.
	uint64_t x = itr->r == 64? itr->p[0] : itr->p[0] << (64 - itr->r) | itr->p[1] >> itr->r;
	if (x>>63 == 0) {
		int64_t y;
		int l, w = 0x333333335555779bll>>(x>>59<<2)&0xf;
		l = (x >> (64 - w)) - 1;
		y = x << w >> (64 - l) | 1u << l;
		w += l;
		*c = x << w >> 61;
		w += 3;
		itr->r -= w;
		if (itr->r <= 0) ++itr->p, itr->r += 64;
		return y;
	} else {
		*c = x << 1 >> 61;
		itr->r -= 4;
		if (itr->r <= 0) ++itr->p, itr->r += 64;
		return 1;
	}
}
#endif

static inline uint64_t rld_locate_blk(const rld_t *e, rlditr_t *itr, uint64_t k, uint64_t *cnt, uint64_t *sum)
{
	int j;
	uint64_t c = 0, *q, *z = e->frame + (k>>e->ibits) * e->asize1;
	itr->i = e->z + (*z>>RLD_LBITS);
	q = itr->p = *itr->i + (*z&RLD_LMASK);
	for (j = 1, *sum = 0; j < e->asize1; ++j) *sum += (cnt[j-1] = z[j]);
	while (1) { // seek to the small block
		int type;
		q += e->ssize;
		if (q - *itr->i == RLD_LSIZE) q = *++itr->i;
		type = rld_block_type(*q);
		c = type == 2? *q&0x3fffffffffffffffULL : type == 1? *(uint32_t*)q : *(uint16_t*)q;
		if (*sum + c > k) break;
		if (type == 0) {
			uint16_t *p = (uint16_t*)q + 1;
#ifdef _DNA_ONLY
			cnt[0] += p[0]; cnt[1] += p[1]; cnt[2] += p[2]; cnt[3] += p[3]; cnt[4] += p[4]; cnt[5] += p[5];
#else
			for (j = 0; j < e->asize; ++j) cnt[j] += p[j];
#endif
		} else if (type == 1) {
			uint32_t *p = (uint32_t*)q + 1;
			for (j = 0; j < e->asize; ++j) cnt[j] += p[j] & 0x3fffffff;
		} else {
			uint64_t *p = (uint64_t*)q + 1;
			for (j = 0; j < e->asize; ++j) cnt[j] += p[j];
		}
		*sum += c;
		itr->p = q;
	}
	itr->shead = itr->p;
	itr->stail = rld_get_stail(e, itr);
	itr->p += e->offset0[rld_block_type(*itr->shead)];
	itr->q = (uint8_t*)itr->p;
	itr->r = 64;
	return c + *sum;
}

void rld_rank21(const rld_t *e, uint64_t k, uint64_t l, int c, uint64_t *ok, uint64_t *ol) // FIXME: can be faster
{
	*ok = rld_rank11(e, k, c);
	*ol = rld_rank11(e, l, c);
}

int rld_rank1a(const rld_t *e, uint64_t k, uint64_t *ok)
{
	uint64_t z, l;
	int a = -1;
	rlditr_t itr;
	if (k == 0) {
		for (a = 0; a < e->asize; ++a) ok[a] = 0;
		return -1;
	}
	rld_locate_blk(e, &itr, k-1, ok, &z);
	while (1) {
#ifdef _DNA_ONLY
		l = rld_dec0_fast_dna(e, &itr, &a);
#else
		l = rld_dec0(e, &itr, &a);
#endif
		if (z + l >= k) break;
		z += l; ok[a] += l;
	}
	ok[a] += k - z;
	return a;
}

uint64_t rld_rank11(const rld_t *e, uint64_t k, int c)
{
	uint64_t *ok;
	if (k == (uint64_t)-1) return 0;
	ok = alloca(e->asize1 * 8);
	rld_rank1a(e, k, ok);
	return ok[c];
}

void rld_rank2a(const rld_t *e, uint64_t k, uint64_t l, uint64_t *ok, uint64_t *ol)
{
	uint64_t z, y, len;
	rlditr_t itr;
	int a = -1;
	if (k == 0) {
		for (a = 0; a < e->asize; ++a) ok[a] = 0;
		rld_rank1a(e, l, ol);
		return;
	}
	y = rld_locate_blk(e, &itr, k-1, ok, &z); // locate the block bracketing k
	while (1) { // compute ok[]
#ifdef _DNA_ONLY
		len = rld_dec0_fast_dna(e, &itr, &a);
#else
		len = rld_dec0(e, &itr, &a);
#endif
		if (z + len >= k) break;
		z += len; ok[a] += len;
	}
	if (y > l) { // we do not need to decode other blocks
		int b;
		for (b = 0; b < e->asize; ++b) ol[b] = ok[b]; // copy ok[] to ol[]
		ok[a] += k - z; // finalize ok[a]
		if (z + len < l) { // we need to decode the next run
			z += len; ol[a] += len;
			while (1) {
				len = rld_dec0(e, &itr, &a);
				if (z + len >= l) break;
				z += len; ol[a] += len;
			}
		}
		ol[a] += l - z;
	} else { // we have to decode other blocks
		ok[a] += k - z;
		rld_rank1a(e, l, ol);
	}
}

int rld_extend(const rld_t *e, const rldintv_t *ik, rldintv_t ok[6], int is_back)
{ // TODO: this can be accelerated a little by using rld_rank1a() when ik.x[2]==1
	uint64_t tk[6], tl[6];
	int i;
	rld_rank2a(e, ik->x[!is_back], ik->x[!is_back] + ik->x[2], tk, tl);
	for (i = 0; i < 6; ++i) {
		ok[i].x[!is_back] = e->cnt[i] + tk[i];
		ok[i].x[2] = (tl[i] -= tk[i]);
	}
	ok[0].x[is_back] = ik->x[is_back];
	ok[4].x[is_back] = ok[0].x[is_back] + tl[0];
	ok[3].x[is_back] = ok[4].x[is_back] + tl[4];
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[5].x[is_back] = ok[1].x[is_back] + tl[1];
	return 0;
}
