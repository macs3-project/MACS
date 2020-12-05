#include <assert.h>
#include <string.h>
#include <math.h>
#include "kvec.h"
#include "kstring.h"
#include "rld0.h"
#include "mag.h"
#include "internal.h"

/******************
 *** From fermi ***
 ******************/

typedef struct { size_t n, m; int32_t *a; } fm32s_v;
typedef struct { size_t n, m; rldintv_t *a; } rldintv_v;

static uint64_t utg_primes[] = { 123457, 234571, 345679, 456791, 567899, 0 };

#define fm6_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))
#define fm6_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(int)(c)], (ik).x[2] = (e)->cnt[(int)(c)+1] - (e)->cnt[(int)(c)], (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

int rld_extend0(const rld_t *e, const rldintv_t *ik, rldintv_t *ok0, int is_back)
{ // FIXME: this can be accelerated a little by using rld_rank1a() when ik.x[2]==1
	uint64_t tk[6], tl[6];
	rld_rank2a(e, ik->x[!is_back], ik->x[!is_back] + ik->x[2], tk, tl);
	ok0->x[!is_back] = tk[0];
	ok0->x[is_back] = ik->x[is_back];
	ok0->x[2] = tl[0] - tk[0];
	return 0;
}

uint64_t fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s, rldintv_t *k2, int *contained)
{
	uint64_t k = x, ok[6];
	rldintv_t ok2[6];
	s->l = 0; *contained = 0;
	while (1) {
		int c = rld_rank1a(e, k + 1, ok);
		k = e->cnt[c] + ok[c] - 1;
		if (c == 0) break;
		if (s->l > 0) {
			if (k2->x[2] == 1) k2->x[0] = k;
			else {
				rld_extend(e, k2, ok2, 1);
				*k2 = ok2[c];
			}
		} else fm6_set_intv(e, c, *k2);
		kputc(c, s);
	}
	if (k2->x[2] != 1) {
		rld_extend(e, k2, ok2, 1);
		if (ok2[0].x[2] != k2->x[2]) *contained |= 1; // left contained
		*k2 = ok2[0];
	} else k2->x[0] = k;
	rld_extend(e, k2, ok2, 0);
	if (ok2[0].x[2] != k2->x[2]) *contained |= 2; // right contained
	*k2 = ok2[0];
	return k;
}

/*****************
 *** Main body ***
 *****************/

#define info_lt(a, b) ((a).info < (b).info)

#include "ksort.h"
KSORT_INIT(infocmp, rldintv_t, info_lt)

static inline void set_bit(uint64_t *bits, uint64_t x)
{
	uint64_t *p = bits + (x>>6);
	uint64_t z = 1LLU<<(x&0x3f);
	__sync_fetch_and_or(p, z);
}

static inline void set_bits(uint64_t *bits, const rldintv_t *p)
{
	uint64_t k;
	for (k = 0; k < p->x[2]; ++k) {
		set_bit(bits, p->x[0] + k);
		set_bit(bits, p->x[1] + k);
	}
}

static rldintv_t overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, rldintv_v *p, int inc_sentinel)
{ // requirement: seq[j] matches the end of a read
	int c, depth, dir, end;
	rldintv_t ik, ok[6];
	p->n = 0;
	dir = at5? 1 : -1; // at5 is true iff we start from the 5'-end of a read
	end = at5? len : -1;
	c = seq[j];
	fm6_set_intv(e, c, ik);
	for (depth = 1, j += dir; j != end; j += dir, ++depth) {
		c = at5? fm6_comp(seq[j]) : seq[j];
		rld_extend(e, &ik, ok, !at5);
		if (!ok[c].x[2]) break; // cannot be extended
		if (depth >= min && ok[0].x[2]) {
			if (inc_sentinel) {
				ok[0].info = j - dir;
				kv_push(rldintv_t, *p, ok[0]);
			} else {
				ik.info = j - dir;
				kv_push(rldintv_t, *p, ik);
			}
		}
		ik = ok[c];
	}
	kv_reverse(rldintv_t, *p, 0); // reverse the array such that the smallest interval comes first
	return ik;
}

typedef struct {
	const rld_t *e;
	int min_match, min_merge_len;
	rldintv_v a[2], nei;
	fm32s_v cat;
	uint64_t *used, *bend;
	kstring_t str;
	uint64_t n, sum, sum2, unpaired;
} aux_t;

int fm6_is_contained(const rld_t *e, int min_match, const kstring_t *s, rldintv_t *intv, rldintv_v *ovlp)
{ // for s is a sequence in e, test if s is contained in other sequences in e; return intervals right overlapping with s
	rldintv_t ik, ok[6];
	int ret = 0;
	assert(s->l > min_match);
	ovlp->n = 0;
	ik = overlap_intv(e, s->l, (uint8_t*)s->s, min_match, s->l - 1, 0, ovlp, 0);
	rld_extend(e, &ik, ok, 1); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -1; // the sequence is left contained
	ik = ok[0];
	rld_extend(e, &ik, ok, 0); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -1; // the sequence is right contained
	*intv = ok[0];
	return ret;
}

int fm6_get_nei(const rld_t *e, int min_match, int beg, kstring_t *s, rldintv_v *nei, // input and output variables
				rldintv_v *prev, rldintv_v *curr, fm32s_v *cat,                        // temporary arrays
				uint64_t *used)                                                      // optional info
{
	int ori_l = s->l, j, i, c, rbeg, is_forked = 0;
	rldintv_v *swap;
	rldintv_t ok[6], ok0;

	curr->n = nei->n = cat->n = 0;
	if (prev->n == 0) { // when this routine is called for the seed, prev may filled by fm6_is_contained()
		overlap_intv(e, s->l - beg, (uint8_t*)s->s + beg, min_match, s->l - beg - 1, 0, prev, 0);
		if (prev->n == 0) return -1; // no overlap
		for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
	}
	kv_resize(int, *cat, prev->m);
	for (j = 0; j < prev->n; ++j) cat->a[j] = 0; // only one interval; all point to 0
	while (prev->n) {
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			rldintv_t *p = &prev->a[j];
			if (cat->a[j] < 0) continue;
			rld_extend(e, p, ok, 0); // forward extension
			if (ok[0].x[2] && ori_l != s->l) { // some (partial) reads end here
				rld_extend0(e, &ok[0], &ok0, 1); // backward extension to look for sentinels
				if (ok0.x[2]) { // the match is bounded by sentinels - a full-length match
					if (ok[0].x[2] == p->x[2] && p->x[2] == ok0.x[2]) { // never consider a read contained in another read
						int cat0 = cat->a[j]; // a category approximately corresponds to one neighbor, though not always
						assert(j == 0 || cat->a[j] > cat->a[j-1]); // otherwise not irreducible
						ok0.info = ori_l - (p->info&0xffffffffU);
						for (i = j; i < prev->n && cat->a[i] == cat0; ++i) cat->a[i] = -1; // mask out other intervals of the same cat
						kv_push(rldintv_t, *nei, ok0); // keep in the neighbor vector
						continue; // no need to go through for(c); do NOT set "used" as this neighbor may be rejected later
					} else if (used) set_bits(used, &ok0); // the read is contained in another read; mark it as used
				}
			} // ~if(ok[0].x[2])
			if (cat->a[j] < 0) continue; // no need to proceed if we have finished this path
			for (c = 1; c < 5; ++c) // collect extensible intervals
				if (ok[c].x[2]) {
					rld_extend0(e, &ok[c], &ok0, 1);
					if (ok0.x[2]) { // do not extend intervals whose left end is not bounded by a sentinel
						ok[c].info = (p->info&0xfffffff0ffffffffLLU) | (uint64_t)c<<32;
						kv_push(rldintv_t, *curr, ok[c]);
					}
				}
		} // ~for(j)
		if (curr->n) { // update category
			uint32_t last, cat0;
			kv_resize(int, *cat, curr->m);
			c = curr->a[0].info>>32&0xf;
			kputc(fm6_comp(c), s);
			ks_introsort(infocmp, curr->n, curr->a);
			last = curr->a[0].info >> 32;
			cat->a[0] = 0;
			curr->a[0].info &= 0xffffffff;
			for (j = 1, cat0 = 0; j < curr->n; ++j) { // this loop recalculate cat
				if (curr->a[j].info>>32 != last)
					last = curr->a[j].info>>32, cat0 = j;
				cat->a[j] = cat0;
				curr->a[j].info = (curr->a[j].info&0xffffffff) | (uint64_t)cat0<<36;
			}
			if (cat0 != 0) is_forked = 1;
		}
		swap = curr; curr = prev; prev = swap; // swap curr and prev
	} // ~while(prev->n)
	if (nei->n == 0) return -1; // no overlap
	rbeg = ori_l - (uint32_t)nei->a[0].info;
	if (nei->n == 1 && is_forked) { // this may happen if there are contained reads; fix this
		fm6_set_intv(e, 0, ok0);
		for (i = rbeg; i < ori_l; ++i) {
			rld_extend(e, &ok0, ok, 0);
			ok0 = ok[fm6_comp(s->s[i])];
		}
		for (i = ori_l; i < s->l; ++i) {
			int c0 = -1;
			rld_extend(e, &ok0, ok, 0);
			for (c = 1, j = 0; c < 5; ++c)
				if (ok[c].x[2] && ok[c].x[0] <= nei->a[0].x[0] && ok[c].x[0] + ok[c].x[2] >= nei->a[0].x[0] + nei->a[0].x[2])
					++j, c0 = c;
			if (j == 0 && ok[0].x[2]) break;
			assert(j == 1);
			s->s[i] = fm6_comp(c0);
			ok0 = ok[c0];
		}
		s->l = i; s->s[s->l] = 0;
	}
	if (nei->n > 1) s->l = ori_l, s->s[s->l] = 0;
	return rbeg;
}

static int try_right(aux_t *a, int beg, kstring_t *s)
{
	return fm6_get_nei(a->e, a->min_match, beg, s, &a->nei, &a->a[0], &a->a[1], &a->cat, a->used);
}

static int check_left_simple(aux_t *a, int beg, int rbeg, const kstring_t *s)
{
	rldintv_t ok[6];
	rldintv_v *prev = &a->a[0], *curr = &a->a[1], *swap;
	int i, j;

	overlap_intv(a->e, s->l, (uint8_t*)s->s, a->min_match, rbeg, 1, prev, 1);
	for (i = rbeg - 1; i >= beg; --i) {
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			rldintv_t *p = &prev->a[j];
			rld_extend(a->e, p, ok, 1);
			if (ok[0].x[2]) set_bits(a->used, &ok[0]); // some reads end here; they must be contained in a longer read
			if (ok[0].x[2] + ok[(int)s->s[i]].x[2] != p->x[2]) return -1; // potential backward bifurcation
			kv_push(rldintv_t, *curr, ok[(int)s->s[i]]);
		}
		swap = curr; curr = prev; prev = swap;
	} // ~for(i)
	return 0;
}

static int check_left(aux_t *a, int beg, int rbeg, const kstring_t *s)
{
	int i, ret;
	rldintv_t tmp;
	assert(a->nei.n == 1);
	ret = check_left_simple(a, beg, rbeg, s);
	if (ret == 0) return 0;
	// when ret<0, the back fork may be caused by a contained read. we have to do more to confirm this.
	tmp = a->nei.a[0]; // backup the neighbour as it will be overwritten by try_right()
	a->a[0].n = a->a[1].n = a->nei.n = 0;
	ks_resize(&a->str, s->l - rbeg + 1);
	for (i = s->l - 1, a->str.l = 0; i >= rbeg; --i)
		a->str.s[a->str.l++] = fm6_comp(s->s[i]);
	a->str.s[a->str.l] = 0;
	try_right(a, 0, &a->str);
	assert(a->nei.n >= 1);
	ret = a->nei.n > 1? -1 : 0;
	a->nei.n = 1; a->nei.a[0] = tmp; // recover the original neighbour
	return ret;
}

static int unitig_unidir(aux_t *a, kstring_t *s, kstring_t *cov, int beg0, uint64_t k0, uint64_t *end, int *is_loop)
{
	int i, beg = beg0, rbeg, ori_l = s->l, n_reads = 0;
	*is_loop = 0;
	while ((rbeg = try_right(a, beg, s)) >= 0) { // loop if there is at least one overlap
		uint64_t k;
		if (a->nei.n > 1) { // forward bifurcation
			set_bit(a->bend, *end);
			break;
		}
		if ((k = a->nei.a[0].x[0]) == *end) break; // a loop like b>>c>>a><a; keep the link but stop extension
		if (((a->bend[k>>6]>>(k&0x3f)&1) || check_left(a, beg, rbeg, s) < 0)) { // backward bifurcation
			set_bit(a->bend, k);
			break;
		}
		if (k == k0) { // a loop like a>>b>>c>>a
			*is_loop = 1;
			break;
		}
		if (a->nei.a[0].x[1] == *end) { // a loop like b>>c>>a>>a; cut the last link
			a->nei.n = 0;
			break;
		}
		if ((int)a->nei.a[0].info < a->min_merge_len) break; // the overlap is not long enough
		*end = a->nei.a[0].x[1];
		set_bits(a->used, &a->nei.a[0]); // successful extension
		++n_reads;
		if (cov->m < s->m) ks_resize(cov, s->m);
		cov->l = s->l; cov->s[cov->l] = 0;
		for (i = rbeg; i < ori_l; ++i) // update the coverage string
			if (cov->s[i] != '~') ++cov->s[i];
		for (i = ori_l; i < s->l; ++i) cov->s[i] = '"';
		beg = rbeg; ori_l = s->l; a->a[0].n = a->a[1].n = 0; // prepare for the next round of loop
	}
	cov->l = s->l = ori_l; s->s[ori_l] = cov->s[ori_l] = 0;
	return n_reads;
}

static void copy_nei(ku128_v *dst, const rldintv_v *src)
{
	int i;
	for (i = 0; i < src->n; ++i) {
		ku128_t z;
		z.x = src->a[i].x[0]; z.y = src->a[i].info;
		kv_push(ku128_t, *dst, z);
	}
}

static int unitig1(aux_t *a, int64_t seed, kstring_t *s, kstring_t *cov, uint64_t end[2], ku128_v nei[2], int *n_reads)
{
	rldintv_t intv0;
	int seed_len, ret, is_loop, contained;
	int64_t k;
	size_t i;

	*n_reads = nei[0].n = nei[1].n = 0;
	if (a->used[seed>>6]>>(seed&0x3f)&1) return -2; // used
	// retrieve the sequence pointed by seed
	k = fm6_retrieve(a->e, seed, s, &intv0, &contained);
	seq_reverse(s->l, (uint8_t*)s->s);
	seed_len = s->l;
	// check contained status
	if (intv0.x[2] > 1 && k != intv0.x[0]) return -3; // duplicated, but not the first
	set_bits(a->used, &intv0);
	if (contained) return -3; // contained
	// check length, containment and if used before
	if (s->l <= a->min_match) return -1; // too short
	ret = fm6_is_contained(a->e, a->min_match, s, &intv0, &a->a[0]);
	*n_reads = 1;
	// initialize the coverage string
	if (cov->m < s->m) ks_resize(cov, s->m);
	cov->l = s->l; cov->s[cov->l] = 0;
	for (i = 0; i < cov->l; ++i) cov->s[i] = '"';
	// left-wards extension
	end[0] = intv0.x[1]; end[1] = intv0.x[0];
	if (a->a[0].n) { // no need to extend to the right if there is no overlap
		*n_reads += unitig_unidir(a, s, cov, 0, intv0.x[0], &end[0], &is_loop);
		copy_nei(&nei[0], &a->nei);
		if (is_loop) {
			ku128_t z;
			z.x = end[0]; z.y = a->nei.a[0].info;
			kv_push(ku128_t, nei[1], z);
			return 0;
		}
	}
	// right-wards extension
	a->a[0].n = a->a[1].n = a->nei.n = 0;
	seq_revcomp6(s->l, (uint8_t*)s->s); // reverse complement for extension in the other direction
	seq_reverse(cov->l, (uint8_t*)cov->s); // reverse the coverage
	*n_reads += unitig_unidir(a, s, cov, s->l - seed_len, intv0.x[1], &end[1], &is_loop);
	copy_nei(&nei[1], &a->nei);
	return 0;
}

typedef struct {
	long max_l;
	aux_t a;
	kstring_t str, cov;
	magv_t z;
	magv_v v;
} thrdat_t;

typedef struct {
	uint64_t prime, *used, *bend, *visited;
	const rld_t *e;
	thrdat_t *d;
} worker_t;

static void worker(void *data, long _i, int tid)
{
	worker_t *w = (worker_t*)data;
	thrdat_t *d = &w->d[tid];
	uint64_t i = (w->prime * _i) % w->e->mcnt[1];
	if (unitig1(&d->a, i, &d->str, &d->cov, d->z.k, d->z.nei, &d->z.nsr) >= 0) { // then we keep the unitig
		uint64_t *p[2], x[2];
		magv_t *q;
		p[0] = w->visited + (d->z.k[0]>>6); x[0] = 1LLU<<(d->z.k[0]&0x3f);
		p[1] = w->visited + (d->z.k[1]>>6); x[1] = 1LLU<<(d->z.k[1]&0x3f);
		if ((__sync_fetch_and_or(p[0], x[0])&x[0]) || (__sync_fetch_and_or(p[1], x[1])&x[1])) return;
		d->z.len = d->str.l;
		if (d->max_l < d->str.m) {
			d->max_l = d->str.m;
			d->z.seq = realloc(d->z.seq, d->max_l);
			d->z.cov = realloc(d->z.cov, d->max_l);
		}
		memcpy(d->z.seq, d->str.s, d->z.len);
		memcpy(d->z.cov, d->cov.s, d->z.len + 1);
		kv_pushp(magv_t, d->v, &q);
		mag_v_copy_to_empty(q, &d->z);
	}
}

mag_t *fml_fmi2mag_core(const rld_t *e, int min_match, int min_merge_len, int n_threads)
{
	extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
	worker_t w;
	int j;
	mag_t *g;

	w.used    = (uint64_t*)calloc((e->mcnt[1] + 63)/64, 8);
	w.bend    = (uint64_t*)calloc((e->mcnt[1] + 63)/64, 8);
	w.visited = (uint64_t*)calloc((e->mcnt[1] + 63)/64, 8);
	w.e       = e;
	assert(e->mcnt[1] >= n_threads * 2);
	w.d = calloc(n_threads, sizeof(thrdat_t));
	w.prime = 0;
	for (j = 0; utg_primes[j] > 0; ++j)
		if (e->mcnt[1] % utg_primes[j] != 0) {
			w.prime = utg_primes[j];
			break;
		}
	assert(w.prime);
	for (j = 0; j < n_threads; ++j) {
		w.d[j].a.e = e; w.d[j].a.min_match = min_match; w.d[j].a.min_merge_len = min_merge_len;
		w.d[j].a.used = w.used; w.d[j].a.bend = w.bend;
	}
	kt_for(n_threads, worker, &w, e->mcnt[1]);
	g = (mag_t*)calloc(1, sizeof(mag_t));
	for (j = 0; j < n_threads; ++j) {
		kv_resize(magv_t, g->v, g->v.n + w.d[j].v.n);
		memcpy(g->v.a + g->v.n, w.d[j].v.a, w.d[j].v.n * sizeof(magv_t));
		g->v.n += w.d[j].v.n;
		free(w.d[j].v.a);
		free(w.d[j].a.a[0].a); free(w.d[j].a.a[1].a); free(w.d[j].a.nei.a); free(w.d[j].a.cat.a);
		free(w.d[j].z.nei[0].a); free(w.d[j].z.nei[1].a); free(w.d[j].z.seq); free(w.d[j].z.cov);
		free(w.d[j].a.str.s); free(w.d[j].str.s); free(w.d[j].cov.s);
	}
	free(w.d); free(w.used); free(w.bend); free(w.visited);

	mag_g_build_hash(g);
	mag_g_amend(g);
	g->rdist = mag_cal_rdist(g);
	return g;
}

mag_t *fml_fmi2mag(const fml_opt_t *opt, rld_t *e)
{
	mag_t *g;
	g = fml_fmi2mag_core(e, opt->min_asm_ovlp, opt->min_merge_len, opt->n_threads);
	rld_destroy(e);
	return g;
}
