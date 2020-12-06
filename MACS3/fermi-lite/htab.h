#ifndef BFC_HTAB_H
#define BFC_HTAB_H

#include <stdint.h>
#include "kmer.h"

#define BFC_CH_KEYBITS 50
#define BFC_CH_MAXPRE  20

struct bfc_ch_s;
typedef struct bfc_ch_s bfc_ch_t;

bfc_ch_t *bfc_ch_init(int k, int l_pre);
void bfc_ch_destroy(bfc_ch_t *ch);
int bfc_ch_insert(bfc_ch_t *ch, const uint64_t x[2], int is_high, int forced);
int bfc_ch_get(const bfc_ch_t *ch, const uint64_t x[2]);
uint64_t bfc_ch_count(const bfc_ch_t *ch);
int bfc_ch_hist(const bfc_ch_t *ch, uint64_t cnt[256], uint64_t high[64]);
int bfc_ch_get_k(const bfc_ch_t *ch);

int bfc_ch_kmer_occ(const bfc_ch_t *ch, const bfc_kmer_t *z);

#endif
