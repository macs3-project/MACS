#include <limits.h>
#include <stdio.h>
#include "mag.h"
#include "kvec.h"
#include "ksw.h"
#include "internal.h"
#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

typedef khash_t(64) hash64_t;

#define MAX_N_DIFF 2.01 // for evaluating alignment after SW
#define MAX_R_DIFF 0.1
#define L_DIFF_COEF 0.2 // n_diff=|l_0 - l_1|*L_DIFF_COEF

#define edge_mark_del(_x) ((_x).x = (uint64_t)-2, (_x).y = 0)
#define edge_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0)

/******************
 * Closed bubbles *
 ******************/

typedef struct {
	uint64_t id;
	int cnt[2];
	int n[2][2], d[2][2];
	uint64_t v[2][2];
} trinfo_t;

const trinfo_t g_trinull = {-1, {0, 0}, {{INT_MIN, INT_MIN}, {INT_MIN, INT_MIN}}, {{INT_MIN, INT_MIN}, {INT_MIN, INT_MIN}}, {{-1, -1}, {-1, -1}}};

typedef struct {
	int n, m;
	trinfo_t **buf;
} tipool_t;

struct mogb_aux {
	tipool_t pool;
	ku64_v stack;
	hash64_t *h;
};

mogb_aux_t *mag_b_initaux(void)
{
	mogb_aux_t *aux = calloc(1, sizeof(mogb_aux_t));
	aux->h = kh_init(64);
	return aux;
}

void mag_b_destroyaux(mogb_aux_t *b)
{
	int i;
	for (i = 0; i < b->pool.m; ++i)
		free(b->pool.buf[i]);
	free(b->pool.buf); free(b->stack.a);
	kh_destroy(64, b->h);
	free(b);
}

#define tiptr(p) ((trinfo_t*)(p)->ptr)

static inline trinfo_t *tip_alloc(tipool_t *pool, uint32_t id)
{ // allocate an object from the memory pool
	trinfo_t *p;
	if (pool->n == pool->m) {
		int i, new_m = pool->m? pool->m<<1 : 256;
		pool->buf = realloc(pool->buf, new_m * sizeof(void*));
		for (i = pool->m; i < new_m; ++i)
			pool->buf[i] = malloc(sizeof(trinfo_t));
		pool->m = new_m;
	}
	p = pool->buf[pool->n++];
	*p = g_trinull;
	p->id = id;
	return p;
}

static void backtrace(mag_t *g, uint64_t end, uint64_t start, hash64_t *h)
{
	while (end>>32 != start) {
		int ret;
		kh_put(64, h, end>>33, &ret);
		end = tiptr(&g->v.a[end>>33])->v[(end>>32^1)&1][end&1];
	}
}

void mag_vh_simplify_bubble(mag_t *g, uint64_t idd, int max_vtx, int max_dist, mogb_aux_t *a)
{
	int i, n_pending = 0;
	magv_t *p, *q;

	p = &g->v.a[idd>>1];
	if (p->len < 0 || p->nei[idd&1].n < 2) return; // stop if p is deleted or it has 0 or 1 neighbor
	// reset aux data
	a->stack.n = a->pool.n = 0;
	if (kh_n_buckets(a->h) >= 64) {
		kh_destroy(64, a->h);
		a->h = kh_init(64);
	} else kh_clear(64, a->h);
	// add the initial vertex
	p->ptr = tip_alloc(&a->pool, idd>>1);
	tiptr(p)->d[(idd&1)^1][0] = -p->len;
	tiptr(p)->n[(idd&1)^1][0] = -p->nsr;
	kv_push(uint64_t, a->stack, idd^1);
	// essentially a topological sorting
	while (a->stack.n) {
		uint64_t x, y;
		ku128_v *r;
		if (a->stack.n == 1 && a->stack.a[0] != (idd^1) && n_pending == 0) break; // found the other end of the bubble
		x = kv_pop(a->stack);
		p = &g->v.a[x>>1];
		//printf("%lld:%lld\n", p->k[0], p->k[1]);
		r = &p->nei[(x&1)^1]; // we will look the the neighbors from the other end of the unitig
		if (a->pool.n > max_vtx || tiptr(p)->d[x&1][0] > max_dist || tiptr(p)->d[x&1][1] > max_dist || r->n == 0) break; // we failed
		// set the distance to p's neighbors
		for (i = 0; i < r->n; ++i) {
			int nsr, dist, which;
			if ((int64_t)r->a[i].x < 0) continue;
			y = mag_tid2idd(g->h, r->a[i].x);
			if (y == (idd^1)) { // there is a loop involving the initial vertex
				a->stack.n = 0;
				break; // not a bubble; stop; this will jump out of the while() loop
			}
			q = &g->v.a[y>>1];
			if (q->ptr == 0) { // has not been attempted
				q->ptr = tip_alloc(&a->pool, y>>1), ++n_pending;
				mag_v128_clean(&q->nei[y&1]); // make sure there are no deleted edges
			}
			nsr  = tiptr(p)->n[x&1][0] + p->nsr; which = 0;
			dist = tiptr(p)->d[x&1][0] + p->len - r->a[i].y;
			//printf("01 [%d]\t[%d,%d]\t[%d,%d]\n", i, tiptr(q)->n[y&1][0], tiptr(q)->n[y&1][1], tiptr(q)->d[y&1][0], tiptr(q)->d[y&1][1]);
			// test and possibly update the tentative distance
			if (nsr > tiptr(q)->n[y&1][0]) { // then move the best to the 2nd best and update the best
				tiptr(q)->n[y&1][1] = tiptr(q)->n[y&1][0]; tiptr(q)->n[y&1][0] = nsr;
				tiptr(q)->v[y&1][1] = tiptr(q)->v[y&1][0]; tiptr(q)->v[y&1][0] = (x^1)<<32|i<<1|which;
				tiptr(q)->d[y&1][1] = tiptr(q)->d[y&1][0]; tiptr(q)->d[y&1][0] = dist;
				nsr  = tiptr(p)->n[x&1][1] + p->nsr; which = 1; // now nsr is the 2nd best
				dist = tiptr(p)->d[x&1][1] + p->len - r->a[i].y;
			}
			if (nsr > tiptr(q)->n[y&1][1]) // update the 2nd best
				tiptr(q)->n[y&1][1] = nsr, tiptr(q)->v[y&1][1] = (x^1)<<32|i<<1|which, tiptr(q)->d[y&1][1] = dist;
			if (++tiptr(q)->cnt[y&1] == q->nei[y&1].n) { // all q's predecessors have been processed; then push
				kv_push(uint64_t, a->stack, y);
				--n_pending;
			}
		}
	}
	if (n_pending == 0 && a->stack.n == 1) { // found a bubble
		uint64_t x = a->stack.a[0];
		p = &g->v.a[x>>1];
		//printf("(%d,%d)\t(%d,%d)\n", tiptr(p)->n[x&1][0], tiptr(p)->n[x&1][1], tiptr(p)->d[x&1][0], tiptr(p)->d[x&1][1]);
		backtrace(g, tiptr(p)->v[x&1][0], idd, a->h);
		backtrace(g, tiptr(p)->v[x&1][1], idd, a->h);
	}
	for (i = 0; i < a->pool.n; ++i) // reset p->ptr
		g->v.a[a->pool.buf[i]->id].ptr = 0;
	if (kh_size(a->h)) { // bubble detected; then remove verticies not in the top two paths
		for (i = 1; i < a->pool.n; ++i) { // i=0 corresponds to the initial vertex which we want to exclude
			uint64_t id = a->pool.buf[i]->id;
			if (id != a->stack.a[0]>>1 && kh_get(64, a->h, id) == kh_end(a->h)) // not in the top two paths
				mag_v_del(g, &g->v.a[id]);
		}
	}
}

void mag_g_simplify_bubble(mag_t *g, int max_vtx, int max_dist)
{
	int64_t i;
	mogb_aux_t *a;
	a = mag_b_initaux();
	for (i = 0; i < g->v.n; ++i) {
		mag_vh_simplify_bubble(g, i<<1|0, max_vtx, max_dist, a);
		mag_vh_simplify_bubble(g, i<<1|1, max_vtx, max_dist, a);
	}
	mag_b_destroyaux(a);
	mag_g_merge(g, 0, 0);
}

int mag_vh_pop_simple(mag_t *g, uint64_t idd, float max_cov, float max_frac, int max_bdiff, int aggressive)
{
	magv_t *p = &g->v.a[idd>>1], *q[2];
	ku128_v *r;
	int i, j, k, dir[2], l[2], ret = -1;
	char *seq[2], *cov[2];
	float n_diff, r_diff, avg[2], max_n_diff = aggressive? MAX_N_DIFF * 2. : MAX_N_DIFF;

	if (p->len < 0 || p->nei[idd&1].n != 2) return ret; // deleted or no bubble
	r = &p->nei[idd&1];
	for (j = 0; j < 2; ++j) {
		uint64_t x;
		if ((int64_t)r->a[j].x < 0) return ret;
		x = mag_tid2idd(g->h, r->a[j].x);
		dir[j] = x&1;
		q[j] = &g->v.a[x>>1];
		if (q[j]->nei[0].n != 1 || q[j]->nei[1].n != 1) return ret; // no bubble
		l[j] = q[j]->len - (int)(q[j]->nei[0].a->y + q[j]->nei[1].a->y); // bubble length, excluding overlaps
	}
	if (q[0]->nei[dir[0]^1].a->x != q[1]->nei[dir[1]^1].a->x) return ret; // no bubble
	if (l[0] - l[1] > max_bdiff || l[1] - l[0] > max_bdiff) return 1; // huge bubble differences
	for (j = 0; j < 2; ++j) { // set seq[] and cov[], and compute avg[]
		if (l[j] > 0) {
			seq[j] = malloc(l[j]<<1);
			cov[j] = seq[j] + l[j];
			for (i = 0; i < l[j]; ++i) {
				seq[j][i] = q[j]->seq[i + q[j]->nei[0].a->y];
				cov[j][i] = q[j]->cov[i + q[j]->nei[0].a->y];
			}
			if (dir[j]) {
				seq_revcomp6(l[j], (uint8_t*)seq[j]);
				seq_reverse(l[j], (uint8_t*)cov[j]);
			}
			for (i = 0, avg[j] = 0.; i < l[j]; ++i) {
				--seq[j][i]; // change DNA6 encoding to DNA4 for SW below
				avg[j] += cov[j][i] - 33;
			}
			avg[j] /= l[j];
		} else { // l[j] <= 0; this may happen around a tandem repeat
			int beg, end;
			seq[j] = cov[j] = 0;
			beg = q[j]->nei[0].a->y; end = q[j]->len - q[j]->nei[1].a->y;
			if (beg > end) beg ^= end, end ^= beg, beg ^= end; // swap
			if (beg < end) {
				for (i = beg, avg[j] = 0.; i < end; ++i)
					avg[j] += q[j]->cov[i] - 33;
				avg[j] /= end - beg;
			} else avg[j] = q[j]->cov[beg] - 33; // FIXME: when q[j] is contained, weird thing may happen
		}
	}
	ret = 1;
	if (l[0] > 0 && l[1] > 0) { // then do SW to compute n_diff and r_diff
		int8_t mat[16];
		kswr_t aln;
		for (i = k = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j)
				mat[k++] = i == j? 5 : -4;
		aln = ksw_align(l[0], (uint8_t*)seq[0], l[1], (uint8_t*)seq[1], 4, mat, 5, 2, 0, 0);
		n_diff = ((l[0] < l[1]? l[0] : l[1]) * 5. - aln.score) / (5. + 4.); // 5: matching score; -4: mismatchig score
		r_diff = n_diff / ((l[0] + l[1]) / 2.);
		//fprintf(stderr, "===> %f %f <===\n", n_diff, r_diff); for (j = 0; j < 2; ++j) { for (i = 0; i < l[j]; ++i) fputc("ACGTN"[(int)seq[j][i]], stderr); fputc('\n', stderr); }
	} else {
		n_diff = abs(l[0] - l[1]) * L_DIFF_COEF;
		r_diff = 1.;
		//fprintf(stderr, "---> (%d,%d) <---\n", l[0], l[1]);
	}
	if (n_diff < max_n_diff || r_diff < MAX_R_DIFF) {
		j = avg[0] < avg[1]? 0 : 1;
		if (aggressive || (avg[j] < max_cov && avg[j] / (avg[j^1] + avg[j]) < max_frac)) {
			mag_v_del(g, q[j]);
			ret = 2;
		}
	}
	free(seq[0]); free(seq[1]);
	return ret;
}

void mag_g_pop_simple(mag_t *g, float max_cov, float max_frac, int min_merge_len, int max_bdiff, int aggressive)
{
	int64_t i, n_examined = 0, n_popped = 0;
	int ret;

	for (i = 0; i < g->v.n; ++i) {
		ret = mag_vh_pop_simple(g, i<<1|0, max_cov, max_frac, max_bdiff, aggressive);
		if (ret >= 1) ++n_examined;
		if (ret >= 2) ++n_popped;
		ret = mag_vh_pop_simple(g, i<<1|1, max_cov, max_frac, max_bdiff, aggressive);
		if (ret >= 1) ++n_examined;
		if (ret >= 2) ++n_popped;
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] examined %ld bubbles and popped %ld\n", __func__, (long)n_examined, (long)n_popped);
	mag_g_merge(g, 0, min_merge_len);
}

/****************
 * Open bubbles *
 ****************/

void mag_v_pop_open(mag_t *g, magv_t *p, int min_elen)
{
	int i, j, k, l, dir, max_l, l_qry;
	magv_t *q, *t;
	ku128_v *r, *s;
	uint8_t *seq;
	int8_t mat[16];

	if (p->len < 0 || p->len >= min_elen) return;
	//if (p->nei[0].n && p->nei[1].n) return; // FIXME: between this and the next line, which is better?
	if (p->nei[0].n + p->nei[1].n != 1) return;
	dir = p->nei[0].n? 0 : 1;
	// initialize the scoring system
	for (i = k = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? 5 : -4;
	
	s = &p->nei[dir];
	for (l = 0; l < s->n; ++l) { // if we use "if (p->nei[0].n + p->nei[1].n != 1)", s->n == 1
		uint64_t v;
		kswq_t *qry;
		if ((int64_t)s->a[l].x < 0) continue;
		v = mag_tid2idd(g->h, s->a[l].x);
		q = &g->v.a[v>>1];
		if (q == p || q->nei[v&1].n == 1) continue;
		// get the query ready
		max_l = (p->len - s->a[l].y) * 2;
		seq = malloc(max_l + 1);
		if (dir == 0) { // forward strand
			for (j = s->a[l].y, k = 0; j < p->len; ++j)
				seq[k++] = p->seq[j] - 1;
		} else { // reverse
			for (j = p->len - s->a[l].y - 1, k = 0; j >= 0; --j)
				seq[k++] = 4 - p->seq[j];
		}
		l_qry = k;
		qry = ksw_qinit(2, l_qry, seq, 4, mat);
		//fprintf(stderr, "===> %lld:%lld:%d[%d], %d, %ld <===\n", p->k[0], p->k[1], s->n, l, p->nsr, q->nei[v&1].n);
		//for (j = 0; j < k; ++j) fputc("ACGTN"[(int)seq[j]], stderr); fputc('\n', stderr);

		r = &q->nei[v&1];
		for (i = 0; i < r->n; ++i) {
			uint64_t w;
			kswr_t aln;
			if (r->a[i].x == p->k[dir] || (int64_t)r->a[i].x < 0) continue;
			w = mag_tid2idd(g->h, r->a[i].x);
			// get the target sequence
			t = &g->v.a[w>>1];
			if (w&1) { // reverse strand
				for (j = t->len - r->a[i].y - 1, k = 0; j >= 0 && k < max_l; --j)
					seq[k++] = 4 - t->seq[j];
			} else {
				for (j = r->a[i].y, k = 0; j < t->len && k < max_l; ++j)
					seq[k++] = t->seq[j] - 1;
			}
			aln = ksw_align(0, 0, k, seq, 4, mat, 5, 2, 0, &qry);
			//for (j = 0; j < k; ++j) fputc("ACGTN"[(int)seq[j]], stderr); fprintf(stderr, "\t%d\t%f\n", aln.score, (l_qry * 5. - aln.score) / (5. + 4.));
			if (aln.score >= l_qry * 5 / 2) {
				double r_diff, n_diff;
				n_diff = (l_qry * 5. - aln.score) / (5. + 4.); // 5: matching score; -4: mismatchig score
				r_diff = n_diff / l_qry;
				if (n_diff < MAX_N_DIFF || r_diff < MAX_R_DIFF) break;
			}
		}

		if (i != r->n) {
			// mark delete in p and delete in q
			edge_mark_del(s->a[l]);
			for (i = 0; i < r->n; ++i)
				if (r->a[i].x == p->k[dir])
					edge_mark_del(r->a[i]);
		}
		free(seq); free(qry);
	}

	for (i = 0; i < s->n; ++i)
		if (!edge_is_del(s->a[i])) break;
	if (i == s->n) mag_v_del(g, p); // p is not connected to any other vertices
}

void mag_g_pop_open(mag_t *g, int min_elen)
{
	int64_t i;
	for (i = 0; i < g->v.n; ++i)
		mag_v_pop_open(g, &g->v.a[i], min_elen);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M:%s] popped open bubbles\n", __func__);
	mag_g_merge(g, 0, 0);
}
