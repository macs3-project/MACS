#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "htab.h"
#include "kmer.h"
#include "internal.h"
#include "fml.h"

/*******************
 *** BFC options ***
 *******************/

typedef struct {
	int n_threads, q, k, l_pre;
	int min_cov; // a k-mer is considered solid if the count is no less than this

	int max_end_ext;
	int win_multi_ec;
	float min_trim_frac;

	// these ec options cannot be changed on the command line
	int w_ec, w_ec_high, w_absent, w_absent_high;
	int max_path_diff, max_heap;
} bfc_opt_t;

void bfc_opt_init(bfc_opt_t *opt)
{
	memset(opt, 0, sizeof(bfc_opt_t));
	opt->n_threads = 1;
	opt->q = 20;
	opt->k = -1;
	opt->l_pre = -1;

	opt->min_cov = 4; // in BFC, this defaults to 3 because it has Bloom pre-filter
	opt->win_multi_ec = 10;
	opt->max_end_ext = 5;
	opt->min_trim_frac = .8;

	opt->w_ec = 1;
	opt->w_ec_high = 7;
	opt->w_absent = 3;
	opt->w_absent_high = 1;
	opt->max_path_diff = 15;
	opt->max_heap = 100;
}

/**********************
 *** K-mer counting ***
 **********************/

#define CNT_BUF_SIZE 256

typedef struct { // cache to reduce locking
	uint64_t y[2];
	int is_high;
} insbuf_t;

typedef struct {
	int k, q;
	int n_seqs;
	const bseq1_t *seqs;
	bfc_ch_t *ch;
	int *n_buf;
	insbuf_t **buf;
} cnt_step_t;

bfc_kmer_t bfc_kmer_null = {{0,0,0,0}};

static int bfc_kmer_bufclear(cnt_step_t *cs, int forced, int tid)
{
	int i, k, r;
	if (cs->ch == 0) return 0;
	for (i = k = 0; i < cs->n_buf[tid]; ++i) {
		r = bfc_ch_insert(cs->ch, cs->buf[tid][i].y, cs->buf[tid][i].is_high, forced);
		if (r < 0) cs->buf[tid][k++] = cs->buf[tid][i];
	}
	cs->n_buf[tid] = k;
	return k;
}

static void bfc_kmer_insert(cnt_step_t *cs, const bfc_kmer_t *x, int is_high, int tid)
{
	int k = cs->k;
	uint64_t y[2], hash;
	hash = bfc_kmer_hash(k, x->x, y);
	if (bfc_ch_insert(cs->ch, y, is_high, 0) < 0) {
		insbuf_t *p;
		if (bfc_kmer_bufclear(cs, 0, tid) == CNT_BUF_SIZE)
			bfc_kmer_bufclear(cs, 1, tid);
		p = &cs->buf[tid][cs->n_buf[tid]++];
		p->y[0] = y[0], p->y[1] = y[1], p->is_high = is_high;
	}
}

static void worker_count(void *_data, long k, int tid)
{
	cnt_step_t *cs = (cnt_step_t*)_data;
	const bseq1_t *s = &cs->seqs[k];
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	uint64_t qmer = 0, mask = (1ULL<<cs->k) - 1;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(cs->k, x.x, c);
			qmer = (qmer<<1 | (s->qual == 0 || s->qual[i] - 33 >= cs->q)) & mask;
			if (++l >= cs->k) bfc_kmer_insert(cs, &x, (qmer == mask), tid);
		} else l = 0, qmer = 0, x = bfc_kmer_null;
	}
}

struct bfc_ch_s *fml_count(int n, const bseq1_t *seq, int k, int q, int l_pre, int n_threads)
{
	int i;
	cnt_step_t cs;
	cs.n_seqs = n, cs.seqs = seq, cs.k = k, cs.q = q;
	cs.ch = bfc_ch_init(cs.k, l_pre);
	cs.n_buf = calloc(n_threads, sizeof(int));
	cs.buf = calloc(n_threads, sizeof(void*));
	for (i = 0; i < n_threads; ++i)
		cs.buf[i] = malloc(CNT_BUF_SIZE * sizeof(insbuf_t));
	kt_for(n_threads, worker_count, &cs, cs.n_seqs);
	for (i = 0; i < n_threads; ++i) free(cs.buf[i]);
	free(cs.buf); free(cs.n_buf);
	return cs.ch;
}

/***************
 *** Correct ***
 ***************/

#define BFC_MAX_KMER     63
#define BFC_MAX_BF_SHIFT 37

#define BFC_MAX_PATHS 4
#define BFC_EC_HIST 5
#define BFC_EC_HIST_HIGH 2

#define BFC_EC_MIN_COV_COEF .1

/**************************
 * Sequence struct for ec *
 **************************/

#include "kvec.h"

typedef struct { // NOTE: unaligned memory
	uint8_t b:3, q:1, ob:3, oq:1;
	uint8_t dummy;
	uint16_t lcov:6, hcov:6, solid_end:1, high_end:1, ec:1, absent:1;
	int i;
} ecbase_t;

typedef kvec_t(ecbase_t) ecseq_t;

static int bfc_seq_conv(const char *s, const char *q, int qthres, ecseq_t *seq)
{
	int i, l;
	l = strlen(s);
	kv_resize(ecbase_t, *seq, l);
	seq->n = l;
	for (i = 0; i < l; ++i) {
		ecbase_t *c = &seq->a[i];
		c->b = c->ob = seq_nt6_table[(int)s[i]] - 1;
		c->q = c->oq = !q? 1 : q[i] - 33 >= qthres? 1 : 0;
		if (c->b > 3) c->q = c->oq = 0;
		c->i = i;
	}
	return l;
}

static inline ecbase_t ecbase_comp(const ecbase_t *b)
{
	ecbase_t r = *b;
	r.b = b->b < 4? 3 - b->b : 4;
	r.ob = b->ob < 4? 3 - b->ob : 4;
	return r;
}

static void bfc_seq_revcomp(ecseq_t *seq)
{
	int i;
	for (i = 0; i < seq->n>>1; ++i) {
		ecbase_t tmp;
		tmp = ecbase_comp(&seq->a[i]);
		seq->a[i] = ecbase_comp(&seq->a[seq->n - 1 - i]);
		seq->a[seq->n - 1 - i] = tmp;
	}
	if (seq->n&1) seq->a[i] = ecbase_comp(&seq->a[i]);
}

/***************************
 * Independent ec routines *
 ***************************/

int bfc_ec_greedy_k(int k, int mode, const bfc_kmer_t *x, const bfc_ch_t *ch)
{
	int i, j, max = 0, max_ec = -1, max2 = 0;
	for (i = 0; i < k; ++i) {
		int c = (x->x[1]>>i&1)<<1 | (x->x[0]>>i&1);
		for (j = 0; j < 4; ++j) {
			bfc_kmer_t y = *x;
			int ret;
			if (j == c) continue;
			bfc_kmer_change(k, y.x, i, j);
			ret = bfc_ch_kmer_occ(ch, &y);
			if (ret < 0) continue;
			if ((max&0xff) < (ret&0xff)) max2 = max, max = ret, max_ec = i<<2 | j;
			else if ((max2&0xff) < (ret&0xff)) max2 = ret;
		}
	}
	return (max&0xff) * 3 > mode && (max2&0xff) < 3? max_ec : -1;
}

int bfc_ec_first_kmer(int k, const ecseq_t *s, int start, bfc_kmer_t *x)
{
	int i, l;
	*x = bfc_kmer_null;
	for (i = start, l = 0; i < s->n; ++i) {
		ecbase_t *c = &s->a[i];
		if (c->b < 4) {
			bfc_kmer_append(k, x->x, c->b);
			if (++l == k) break;
		} else l = 0, *x = bfc_kmer_null;
	}
	return i;
}

void bfc_ec_kcov(int k, int min_occ, ecseq_t *s, const bfc_ch_t *ch)
{
	int i, l, r, j;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->n; ++i) {
		ecbase_t *c = &s->a[i];
		c->high_end = c->solid_end = c->lcov = c->hcov = 0;
		if (c->b < 4) {
			bfc_kmer_append(k, x.x, c->b);
			if (++l >= k) {
				if ((r = bfc_ch_kmer_occ(ch, &x)) >= 0) {
					if ((r>>8&0x3f) >= min_occ+1) c->high_end = 1;
					if ((r&0xff) >= min_occ) {
						c->solid_end = 1;
						for (j = i - k + 1; j <= i; ++j)
							++s->a[j].lcov, s->a[j].hcov += c->high_end;
					}
				}
			}
		} else l = 0, x = bfc_kmer_null;
	}
}

uint64_t bfc_ec_best_island(int k, const ecseq_t *s)
{ // IMPORTANT: call bfc_ec_kcov() before calling this function!
	int i, l, max, max_i;
	for (i = k - 1, max = l = 0, max_i = -1; i < s->n; ++i) {
		if (!s->a[i].solid_end) {
			if (l > max) max = l, max_i = i;
			l = 0;
		} else ++l;
	}
	if (l > max) max = l, max_i = i;
	return max > 0? (uint64_t)(max_i - max - k + 1) << 32 | max_i : 0;
}

/********************
 * Correct one read *
 ********************/

#include "ksort.h"

#define ECCODE_MISC      1
#define ECCODE_MANY_N    2
#define ECCODE_NO_SOLID  3
#define ECCODE_UNCORR_N  4
#define ECCODE_MANY_FAIL 5

typedef struct {
	uint32_t ec_code:3, brute:1, n_ec:14, n_ec_high:14;
	uint32_t n_absent:24, max_heap:8;
} ecstat_t;

typedef struct {
	uint8_t ec:1, ec_high:1, absent:1, absent_high:1, b:4;
} bfc_penalty_t;

typedef struct {
	int tot_pen;
	int i; // base position
	int k; // position in the stack
	int32_t ecpos_high[BFC_EC_HIST_HIGH];
	int32_t ecpos[BFC_EC_HIST];
	bfc_kmer_t x;
} echeap1_t;

typedef struct {
	int parent, i, tot_pen;
	uint8_t b;
	bfc_penalty_t pen;
	uint16_t cnt;
} ecstack1_t;

typedef struct {
	const bfc_opt_t *opt;
	const bfc_ch_t *ch;
	kvec_t(echeap1_t) heap;
	kvec_t(ecstack1_t) stack;
	ecseq_t seq, tmp, ec[2];
	int mode;
	ecstat_t ori_st;
} bfc_ec1buf_t;

#define heap_lt(a, b) ((a).tot_pen > (b).tot_pen)
KSORT_INIT(ec, echeap1_t, heap_lt)

static bfc_ec1buf_t *ec1buf_init(const bfc_opt_t *opt, const bfc_ch_t *ch)
{
	bfc_ec1buf_t *e;
	e = calloc(1, sizeof(bfc_ec1buf_t));
	e->opt = opt, e->ch = ch;
	return e;
}

static void ec1buf_destroy(bfc_ec1buf_t *e)
{	
	free(e->heap.a); free(e->stack.a); free(e->seq.a); free(e->tmp.a); free(e->ec[0].a); free(e->ec[1].a);
	free(e);
}

#define weighted_penalty(o, p) ((o)->w_ec * (p).ec + (o)->w_ec_high * (p).ec_high + (o)->w_absent * (p).absent + (o)->w_absent_high * (p).absent_high)

static void buf_update(bfc_ec1buf_t *e, const echeap1_t *prev, bfc_penalty_t pen, int cnt)
{
	ecstack1_t *q;
	echeap1_t *r;
	const bfc_opt_t *o = e->opt;
	int b = pen.b;
	// update stack
	kv_pushp(ecstack1_t, e->stack, &q);
	q->parent = prev->k;
	q->i = prev->i;
	q->b = b;
	q->pen = pen;
	q->cnt = cnt > 0? cnt&0xff : 0;
	q->tot_pen = prev->tot_pen + weighted_penalty(o, pen);
	// update heap
	kv_pushp(echeap1_t, e->heap, &r);
	r->i = prev->i + 1;
	r->k = e->stack.n - 1;
	r->x = prev->x;
	if (pen.ec_high) {
		memcpy(r->ecpos_high + 1, prev->ecpos_high, (BFC_EC_HIST_HIGH - 1) * 4);
		r->ecpos_high[0] = prev->i;
	} else memcpy(r->ecpos_high, prev->ecpos_high, BFC_EC_HIST_HIGH * 4);
	if (pen.ec) {
		memcpy(r->ecpos + 1, prev->ecpos, (BFC_EC_HIST - 1) * 4);
		r->ecpos[0] = prev->i;
	} else memcpy(r->ecpos, prev->ecpos, BFC_EC_HIST * 4);
	r->tot_pen = q->tot_pen;
	bfc_kmer_append(e->opt->k, r->x.x, b);
	ks_heapup_ec(e->heap.n, e->heap.a);
}

static int buf_backtrack(ecstack1_t *s, int end, const ecseq_t *seq, ecseq_t *path)
{
	int i, n_absent = 0;
	kv_resize(ecbase_t, *path, seq->n);
	path->n = seq->n;
	while (end >= 0) {
		if ((i = s[end].i) < seq->n) {
			path->a[i].b = s[end].b;
			path->a[i].ec = s[end].pen.ec;
			path->a[i].absent = s[end].pen.absent;
			n_absent += s[end].pen.absent;
		}
		end = s[end].parent;
	}
	return n_absent;
}

static int bfc_ec1dir(bfc_ec1buf_t *e, const ecseq_t *seq, ecseq_t *ec, int start, int end, int *max_heap)
{
	echeap1_t z;
	int i, l, rv = -1, path[BFC_MAX_PATHS], n_paths = 0, min_path = -1, min_path_pen = INT_MAX, n_failures = 0;
	assert(end <= seq->n && end - start >= e->opt->k);
	e->heap.n = e->stack.n = 0;
	*max_heap = 0;
	memset(&z, 0, sizeof(echeap1_t));
	kv_resize(ecbase_t, *ec, seq->n);
	ec->n = seq->n;
	for (z.i = start, l = 0; z.i < end; ++z.i) {
		int c = seq->a[z.i].b;
		if (c < 4) {
			if (++l == e->opt->k) break;
			bfc_kmer_append(e->opt->k, z.x.x, c);
		} else l = 0, z.x = bfc_kmer_null;
	}
	assert(z.i < end); // before calling this function, there must be at least one solid k-mer
	z.k = -1;
	for (i = 0; i < BFC_EC_HIST; ++i) z.ecpos[i] = -1;
	for (i = 0; i < BFC_EC_HIST_HIGH; ++i) z.ecpos_high[i] = -1;
	kv_push(echeap1_t, e->heap, z);
	for (i = 0; i < seq->n; ++i) ec->a[i].b = seq->a[i].b, ec->a[i].ob = seq->a[i].ob;
	// exhaustive error correction
	while (1) {
		int stop = 0;
		*max_heap = *max_heap > 255? 255 : *max_heap > e->heap.n? *max_heap : e->heap.n;
		if (e->heap.n == 0) { // may happen when there is an uncorrectable "N"
			rv = -2;
			break;
		}
		z = e->heap.a[0];
		e->heap.a[0] = kv_pop(e->heap);
		ks_heapdown_ec(0, e->heap.n, e->heap.a);
		if (min_path >= 0 && z.tot_pen > min_path_pen + e->opt->max_path_diff) break;
		if (z.i - end > e->opt->max_end_ext) stop = 1;
		if (!stop) {
			ecbase_t *c = z.i < seq->n? &seq->a[z.i] : 0;
			int b, os = -1, fixed = 0, other_ext = 0, n_added = 0, added_cnt[4];
			bfc_penalty_t added[4];
			// test if the read extension alone is enough
			if (z.i > end) fixed = 1;
			if (c && c->b < 4) { // A, C, G or T
				bfc_kmer_t x = z.x;
				bfc_kmer_append(e->opt->k, x.x, c->b);
				os = bfc_ch_kmer_occ(e->ch, &x);
				if (c->q && (os&0xff) >= e->opt->min_cov + 1 && c->lcov >= e->opt->min_cov + 1) fixed = 1;
				else if (c->hcov > e->opt->k * .75) fixed = 1;
			}
			// extension
			for (b = 0; b < 4; ++b) {
				bfc_penalty_t pen;
				if (fixed && c && b != c->b) continue;
				if (c == 0 || b != c->b) {
					int s;
					bfc_kmer_t x = z.x;
					pen.ec = 0, pen.ec_high = 0, pen.absent = 0, pen.absent_high = 0, pen.b = b;
					if (c) { // not over the end
						if (c->q && z.ecpos_high[BFC_EC_HIST_HIGH-1] >= 0 && z.i - z.ecpos_high[BFC_EC_HIST_HIGH-1] < e->opt->win_multi_ec) continue; // no close highQ corrections
						if (z.ecpos[BFC_EC_HIST-1] >= 0 && z.i - z.ecpos[BFC_EC_HIST-1] < e->opt->win_multi_ec) continue; // no clustered corrections
					}
					bfc_kmer_append(e->opt->k, x.x, b);
					s = bfc_ch_kmer_occ(e->ch, &x);
					if (s < 0 || (s&0xff) < e->opt->min_cov) continue; // not solid
					//if (os >= 0 && (s&0xff) - (os&0xff) < 2) continue; // not sufficiently better than the read path
					pen.ec = c && c->b < 4? 1 : 0;
					pen.ec_high = pen.ec? c->oq : 0;
					pen.absent = 0;
					pen.absent_high = ((s>>8&0xff) < e->opt->min_cov);
					pen.b = b;
					added_cnt[n_added] = s;
					added[n_added++] = pen;
					++other_ext;
				} else {
					pen.ec = pen.ec_high = 0;
					pen.absent = (os < 0 || (os&0xff) < e->opt->min_cov);
					pen.absent_high = (os < 0 || (os>>8&0xff) < e->opt->min_cov);
					pen.b = b;
					added_cnt[n_added] = os;
					added[n_added++] = pen;
				}
			} // ~for(b)
			if (fixed == 0 && other_ext == 0) ++n_failures;
			if (n_failures > seq->n * 2) {
				rv = -3;
				break;
			}
			if (c || n_added == 1) {
				if (n_added > 1 && e->heap.n > e->opt->max_heap) { // to prevent heap explosion
					int min_b = -1, min = INT_MAX;
					for (b = 0; b < n_added; ++b) {
						int t = weighted_penalty(e->opt, added[b]);
						if (min > t) min = t, min_b = b;
					}
					buf_update(e, &z, added[min_b], added_cnt[min_b]);
				} else {
					for (b = 0; b < n_added; ++b)
						buf_update(e, &z, added[b], added_cnt[b]);
				}
			} else {
				if (n_added == 0)
					e->stack.a[z.k].tot_pen += e->opt->w_absent * (e->opt->max_end_ext - (z.i - end));
				stop = 1;
			}
		} // ~if(!stop)
		if (stop) {
			if (e->stack.a[z.k].tot_pen < min_path_pen)
				min_path_pen = e->stack.a[z.k].tot_pen, min_path = n_paths;
			path[n_paths++] = z.k;
			if (n_paths == BFC_MAX_PATHS) break;
		}
	} // ~while(1)
	// backtrack
	if (n_paths == 0) return rv;
	assert(min_path >= 0 && min_path < n_paths && e->stack.a[path[min_path]].tot_pen == min_path_pen);
	rv = buf_backtrack(e->stack.a, path[min_path], seq, ec);
	for (i = 0; i < ec->n; ++i) // mask out uncorrected regions
		if (i < start + e->opt->k || i >= end) ec->a[i].b = 4;
	return rv;
}

ecstat_t bfc_ec1(bfc_ec1buf_t *e, char *seq, char *qual)
{
	int i, start = 0, end = 0, n_n = 0, rv[2], max_heap[2];
	uint64_t r;
	ecstat_t s;

	s.ec_code = ECCODE_MISC, s.brute = 0, s.n_ec = s.n_ec_high = 0, s.n_absent = s.max_heap = 0;
	bfc_seq_conv(seq, qual, e->opt->q, &e->seq);
	for (i = 0; i < e->seq.n; ++i)
		if (e->seq.a[i].ob > 3) ++n_n;
	if (n_n > e->seq.n * .05) {
		s.ec_code = ECCODE_MANY_N;
		return s;
	}
	bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->ch);
	r = bfc_ec_best_island(e->opt->k, &e->seq);
	if (r == 0) { // no solid k-mer
		bfc_kmer_t x;
		int ec = -1;
		while ((end = bfc_ec_first_kmer(e->opt->k, &e->seq, start, &x)) < e->seq.n) {
			ec = bfc_ec_greedy_k(e->opt->k, e->mode, &x, e->ch);
			if (ec >= 0) break;
			if (end + (e->opt->k>>1) >= e->seq.n) break;
			start = end - (e->opt->k>>1);
		}
		if (ec >= 0) {
			e->seq.a[end - (ec>>2)].b = ec&3;
			++end; start = end - e->opt->k;
			s.brute = 1;
		} else {
			s.ec_code = ECCODE_NO_SOLID;
			return s;
		}
	} else start = r>>32, end = (uint32_t)r;
	if ((rv[0] = bfc_ec1dir(e, &e->seq, &e->ec[0], start, e->seq.n, &max_heap[0])) < 0) {
		s.ec_code = rv[0] == -2? ECCODE_UNCORR_N : rv[0] == -3? ECCODE_MANY_FAIL : ECCODE_MISC;
		return s;
	}
	bfc_seq_revcomp(&e->seq);
	if ((rv[1] = bfc_ec1dir(e, &e->seq, &e->ec[1], e->seq.n - end, e->seq.n, &max_heap[1])) < 0) {
		s.ec_code = rv[1] == -2? ECCODE_UNCORR_N : rv[1] == -3? ECCODE_MANY_FAIL : ECCODE_MISC;
		return s;
	}
	s.max_heap = max_heap[0] > max_heap[1]? max_heap[0] : max_heap[1];
	s.ec_code = 0, s.n_absent = rv[0] + rv[1];
	bfc_seq_revcomp(&e->ec[1]);
	bfc_seq_revcomp(&e->seq);
	for (i = 0; i < e->seq.n; ++i) {
		ecbase_t *c = &e->seq.a[i];
		if (e->ec[0].a[i].b == e->ec[1].a[i].b)
			c->b = e->ec[0].a[i].b > 3? e->seq.a[i].b : e->ec[0].a[i].b;
		else if (e->ec[1].a[i].b > 3) c->b = e->ec[0].a[i].b;
		else if (e->ec[0].a[i].b > 3) c->b = e->ec[1].a[i].b;
		else c->b = e->seq.a[i].ob;
	}
	for (i = 0; i < e->seq.n; ++i) {
		int is_diff = !(e->seq.a[i].b == e->seq.a[i].ob);
		if (is_diff) {
			++s.n_ec;
			if (e->seq.a[i].q) ++s.n_ec_high;
		}
		seq[i] = (is_diff? "acgtn" : "ACGTN")[e->seq.a[i].b];
		if (qual) qual[i] = is_diff? 34 + e->seq.a[i].ob : "+?"[e->seq.a[i].q];
	}
	return s;
}

/********************
 * Error correction *
 ********************/

typedef struct {
	const bfc_opt_t *opt;
	const bfc_ch_t *ch;
	bfc_ec1buf_t **e;
	int64_t n_processed;
	int n_seqs, flt_uniq;
	bseq1_t *seqs;
} ec_step_t;

static uint64_t max_streak(int k, const bfc_ch_t *ch, const bseq1_t *s)
{
	int i, l;
	uint64_t max = 0, t = 0;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) { // not an ambiguous base
			bfc_kmer_append(k, x.x, c);
			if (++l >= k) { // ok, we have a k-mer now
				if (bfc_ch_kmer_occ(ch, &x) > 0) t += 1ULL<<32;
				else t = i + 1;
			} else t = i + 1;
		} else l = 0, x = bfc_kmer_null, t = i + 1;
		max = max > t? max : t;
	}
	return max;
}

static void worker_ec(void *_data, long k, int tid)
{
	ec_step_t *es = (ec_step_t*)_data;
	bseq1_t *s = &es->seqs[k];
	if (es->flt_uniq) {
		uint64_t max;
		//printf("meme\n");
		max = max_streak(es->opt->k, es->ch, s);
		//printf("meme1\n");
		if (max>>32 && (double)((max>>32) + es->opt->k - 1) / s->l_seq > es->opt->min_trim_frac) {
		  //printf("meme1.1\n");
			int start = (uint32_t)max, end = start + (max>>32);
			start -= es->opt->k - 1;
			assert(start >= 0 && end <= s->l_seq);
			memmove(s->seq, s->seq + start, end - start);
			s->l_seq = end - start;
			s->seq[s->l_seq] = 0;
			if (s->qual) {
				memmove(s->qual, s->qual + start, s->l_seq);
				s->qual[s->l_seq] = 0;
			}
		} else {
		  //printf("meme1.2\n");
		        free(s->seq); free(s->qual);
			s->l_seq = 0, s->seq = s->qual = 0;
		}
	} else bfc_ec1(es->e[tid], s->seq, s->qual);
	//printf("meme2\n");
}

float fml_correct_core(const fml_opt_t *opt, int flt_uniq, int n, bseq1_t *seq)
{
	bfc_ch_t *ch;
	int i, mode;
	uint64_t hist[256], hist_high[64], tot_len = 0, sum_k = 0, tot_k = 0;
	ec_step_t es;
	bfc_opt_t bfc_opt;
	float kcov;

	// initialize BFC options
	bfc_opt_init(&bfc_opt);
	bfc_opt.n_threads = opt->n_threads; // copy from FML options
	bfc_opt.k = flt_uniq? opt->min_asm_ovlp : opt->ec_k;
	for (i = 0; i < n; ++i) tot_len += seq[i].l_seq; // compute total length
	bfc_opt.l_pre = tot_len - 8 < 20? tot_len - 8 : 20;
	//printf("hoho!\n");
	memset(&es, 0, sizeof(ec_step_t));
	es.opt = &bfc_opt, es.n_seqs = n, es.seqs = seq, es.flt_uniq = flt_uniq;
	//printf("hoho!\n");
	es.ch = ch = fml_count(n, seq, bfc_opt.k, bfc_opt.q, bfc_opt.l_pre, bfc_opt.n_threads);
	mode = bfc_ch_hist(ch, hist, hist_high);
	for (i = opt->min_cnt; i < 256; ++i)
		sum_k += hist[i], tot_k += i * hist[i];
	kcov = (float)tot_k / sum_k;
	bfc_opt.min_cov = (int)(BFC_EC_MIN_COV_COEF * kcov + .499);
	bfc_opt.min_cov = bfc_opt.min_cov < opt->max_cnt? bfc_opt.min_cov : opt->max_cnt;
	bfc_opt.min_cov = bfc_opt.min_cov > opt->min_cnt? bfc_opt.min_cov : opt->min_cnt;
	//printf("hoho3!\n");
	es.e = calloc(es.opt->n_threads, sizeof(void*));
	//printf("hoho3.1!\n");
	for (i = 0; i < es.opt->n_threads; ++i)
		es.e[i] = ec1buf_init(es.opt, ch), es.e[i]->mode = mode;
	//printf("hoho3.2!\n");
	kt_for(es.opt->n_threads, worker_ec, &es, es.n_seqs);
	//printf("hoho3.3!\n");
	for (i = 0; i < es.opt->n_threads; ++i)
		ec1buf_destroy(es.e[i]);
	//printf("hoho4!\n");
	free(es.e);
	bfc_ch_destroy(ch);
	return kcov;
}

float fml_correct(const fml_opt_t *opt, int n, bseq1_t *seq)
{
	return fml_correct_core(opt, 0, n, seq);
}

float fml_fltuniq(const fml_opt_t *opt, int n, bseq1_t *seq)
{
	return fml_correct_core(opt, 1, n, seq);
}
