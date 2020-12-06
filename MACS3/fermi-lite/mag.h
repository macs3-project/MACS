#ifndef FM_MOG_H
#define FM_MOG_H

#include <stdint.h>
#include <stdlib.h>
#include "kstring.h"
#include "fml.h"

#ifndef KINT_DEF
#define KINT_DEF
typedef struct { uint64_t x, y; } ku128_t;
typedef struct { size_t n, m; uint64_t *a; } ku64_v;
typedef struct { size_t n, m; ku128_t *a; } ku128_v;
#endif

typedef struct {
	int len, nsr;    // length; number supporting reads
	uint32_t max_len;// allocated seq/cov size
	uint64_t k[2];   // bi-interval
	ku128_v nei[2];  // neighbors
	char *seq, *cov; // sequence and coverage
	void *ptr;       // additional information
} magv_t;

typedef struct { size_t n, m; magv_t *a; } magv_v;

typedef struct mag_t {
	magv_v v;
	float rdist;  // read distance
	int min_ovlp; // minimum overlap seen from the graph
	void *h;
} mag_t;

struct mogb_aux;
typedef struct mogb_aux mogb_aux_t;

#ifdef __cplusplus
extern "C" {
#endif

	void mag_init_opt(magopt_t *o);
	void mag_g_clean(mag_t *g, const magopt_t *opt);

	void mag_g_destroy(mag_t *g);
	void mag_g_amend(mag_t *g);
	void mag_g_build_hash(mag_t *g);
	void mag_g_print(const mag_t *g);
	int mag_g_rm_vext(mag_t *g, int min_len, int min_nsr);
	void mag_g_rm_edge(mag_t *g, int min_ovlp, double min_ratio, int min_len, int min_nsr);
	void mag_g_merge(mag_t *g, int rmdup, int min_merge_len);
	void mag_g_simplify_bubble(mag_t *g, int max_vtx, int max_dist);
	void mag_g_pop_simple(mag_t *g, float max_cov, float max_frac, int min_merge_len, int max_bdiff, int aggressive);
	void mag_g_pop_open(mag_t *g, int min_elen);
	void mag_g_trim_open(mag_t *g, const magopt_t *opt);

	void mag_v_copy_to_empty(magv_t *dst, const magv_t *src); // NB: memory leak if dst is allocated
	void mag_v_del(mag_t *g, magv_t *p);
	void mag_v_write(const magv_t *p, kstring_t *out);
	void mag_v_pop_open(mag_t *g, magv_t *p, int min_elen);

	uint64_t mag_tid2idd(void *h, uint64_t tid);
	void mag_v128_clean(ku128_v *r);
	double mag_cal_rdist(mag_t *g);

#ifdef __cplusplus
}
#endif

#endif
