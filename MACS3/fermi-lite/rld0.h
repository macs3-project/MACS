#ifndef RLDELTA0_H
#define RLDELTA0_H

#define _DNA_ONLY

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define RLD_LBITS 23
#define RLD_LSIZE (1<<RLD_LBITS)
#define RLD_LMASK (RLD_LSIZE - 1)

typedef struct {
	int r, c; // $r: bits remained in the last 64-bit integer; $c: pending symbol
	int64_t l; // $l: pending length
	uint64_t *p, *shead, *stail, **i;
	uint8_t *q;
} rlditr_t;

typedef struct rld_t {
	// initialized in the constructor
	uint8_t asize, asize1; // alphabet size; asize1=asize+1
	int8_t abits; // bits required to store a symbol
	int8_t sbits; // bits per small block
	int8_t ibits; // modified during indexing; here for a better alignment
	int8_t offset0[3]; // 0 for 16-bit blocks; 1 for 32-bit blocks; 2 for 64-bit blocks
	int ssize; // ssize = 1<<sbits
	// modified during encoding
	int n; // number of blocks (unchanged in decoding)
	uint64_t n_bytes; // total number of bits (unchanged in decoding)
	uint64_t **z; // the actual data (unchanged in decoding)
	uint64_t *cnt, *mcnt; // after enc_finish, cnt keeps the accumulative count and mcnt keeps the marginal
	// modified during indexing
	uint64_t n_frames;
	uint64_t *frame;
	//
	int fd;
	uint64_t *mem; // only used for memory mapped file
} rld_t;

typedef struct {
	uint64_t x[3]; // 0: start of the interval, backward; 1: forward; 2: size of the interval
	uint64_t info;
} rldintv_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_init(int asize, int bbits);
	void rld_destroy(rld_t *e);
	int rld_dump(const rld_t *e, const char *fn);
	rld_t *rld_restore(const char *fn);
	rld_t *rld_restore_mmap(const char *fn);

	void rld_itr_init(const rld_t *e, rlditr_t *itr, uint64_t k);
	int rld_enc(rld_t *e, rlditr_t *itr, int64_t l, uint8_t c);
	uint64_t rld_enc_finish(rld_t *e, rlditr_t *itr);

	uint64_t rld_rank11(const rld_t *e, uint64_t k, int c);
	int rld_rank1a(const rld_t *e, uint64_t k, uint64_t *ok);
	void rld_rank21(const rld_t *e, uint64_t k, uint64_t l, int c, uint64_t *ok, uint64_t *ol);
	void rld_rank2a(const rld_t *e, uint64_t k, uint64_t l, uint64_t *ok, uint64_t *ol);

	int rld_extend(const rld_t *e, const rldintv_t *ik, rldintv_t ok[6], int is_back);

#ifdef __cplusplus
}
#endif

#define rld_last_blk(e) ((e)->n_bytes>>3>>(e)->sbits<<(e)->sbits)
#define rld_seek_blk(e, k) ((e)->z[(k)>>RLD_LBITS] + ((k)&RLD_LMASK))
#define rld_get_stail(e, itr) ((itr)->shead + (e)->ssize - ((itr)->shead + (e)->ssize - *(itr)->i == RLD_LSIZE? 2 : 1))

#define rld_block_type(x) ((uint64_t)(x)>>62)

static inline int64_t rld_dec0(const rld_t *e, rlditr_t *itr, int *c)
{
	int w;
	uint64_t x;
	int64_t l, y = 0;
	x = itr->p[0] << (64 - itr->r) | (itr->p != itr->stail && itr->r != 64? itr->p[1] >> itr->r : 0);
	if (x>>63 == 0) {
		if ((w = 0x333333335555779bll>>(x>>59<<2)&0xf) == 0xb && x>>58 == 0) return 0;
		l = (x >> (64 - w)) - 1;
		y = x << w >> (64 - l) | 1u << l;
		w += l;
	} else w = y = 1;
	*c = x << w >> (64 - e->abits);
	w += e->abits;
	if (itr->r > w) itr->r -= w;
	else ++itr->p, itr->r = 64 + itr->r - w;
	return y;
}

static inline int64_t rld_dec(const rld_t *e, rlditr_t *itr, int *_c, int is_free)
{
	int64_t l = rld_dec0(e, itr, _c);
	if (l == 0 || *_c > e->asize) {
		uint64_t last = rld_last_blk(e);
		if (itr->p - *itr->i > RLD_LSIZE - e->ssize) {
			if (is_free) {
				free(*itr->i); *itr->i = 0;
			}
			itr->shead = *++itr->i;
		} else itr->shead += e->ssize;
		if (itr->shead == rld_seek_blk(e, last)) return -1;
		itr->p = itr->shead + e->offset0[rld_block_type(*itr->shead)];
		itr->q = (uint8_t*)itr->p;
		itr->stail = rld_get_stail(e, itr);
		itr->r = 64;
		return rld_dec0(e, itr, _c);
	} else return l;
}

// take k symbols from e0 and write it to e
static inline void rld_dec_enc(rld_t *e, rlditr_t *itr, const rld_t *e0, rlditr_t *itr0, int64_t k)
{
	if (itr0->l >= k) { // there are more pending symbols
		rld_enc(e, itr, k, itr0->c);
		itr0->l -= k; // l - k symbols remains
	} else { // use up all pending symbols
		int c = -1; // to please gcc
		int64_t l;
		rld_enc(e, itr, itr0->l, itr0->c); // write all pending symbols
		k -= itr0->l;
		for (; k > 0; k -= l) { // we always go into this loop because l0<k
			l = rld_dec(e0, itr0, &c, 1);
			rld_enc(e, itr, k < l? k : l, c);
		}
		itr0->l = -k; itr0->c = c;
	}
}

#endif
