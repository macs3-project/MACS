#ifndef BFC_KMER_H
#define BFC_KMER_H

#include <stdint.h>

typedef struct {
	uint64_t x[4];
} bfc_kmer_t;

static inline void bfc_kmer_append(int k, uint64_t x[4], int c)
{ // IMPORTANT: 0 <= c < 4
	uint64_t mask = (1ULL<<k) - 1;
	x[0] = (x[0]<<1 | (c&1))  & mask;
	x[1] = (x[1]<<1 | (c>>1)) & mask;
	x[2] = x[2]>>1 | (1ULL^(c&1))<<(k-1);
	x[3] = x[3]>>1 | (1ULL^c>>1) <<(k-1);
}

static inline void bfc_kmer_change(int k, uint64_t x[4], int d, int c) // d-bp from the 3'-end of k-mer; 0<=d<k
{ // IMPORTANT: 0 <= c < 4
	uint64_t t = ~(1ULL<<d);
	x[0] = (uint64_t) (c&1)<<d | (x[0]&t);
	x[1] = (uint64_t)(c>>1)<<d | (x[1]&t);
	t = ~(1ULL<<(k-1-d));
	x[2] = (uint64_t)(1^(c&1))<<(k-1-d) | (x[2]&t);
	x[3] = (uint64_t)(1^ c>>1)<<(k-1-d) | (x[3]&t);
}

// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
static inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

static inline uint64_t bfc_hash_64_inv(uint64_t key, uint64_t mask)
{
	uint64_t tmp;
 
	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;
 
	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;
 
	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;
 
	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;
 
	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;
 
	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;
 
	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;
 
	return key;
}

static inline uint64_t bfc_kmer_hash(int k, const uint64_t x[4], uint64_t h[2])
{
	int t = k>>1, u = ((x[1]>>t&1) > (x[3]>>t&1)); // the middle base is always different
	uint64_t mask = (1ULL<<k) - 1, ret;
	h[0] = bfc_hash_64((x[u<<1|0] + x[u<<1|1]) & mask, mask);
	h[1] = bfc_hash_64(h[0] ^ x[u<<1|1], mask);
	ret = (h[0] ^ h[1]) << k | ((h[0] + h[1]) & mask);
	h[0] = (h[0] + h[1]) & mask;
	return ret;
}

static inline void bfc_kmer_hash_inv(int k, const uint64_t h[2], uint64_t y[2])
{
	uint64_t mask = (1ULL<<k) - 1, t = (h[0] - h[1]) & mask;
	y[1] = bfc_hash_64_inv(h[1], mask) ^ t;
	y[0] = (bfc_hash_64_inv(t, mask) - y[1]) & mask;
}

static inline char *bfc_kmer_2str(int k, const uint64_t y[2], char *buf)
{
	int l;
	for (l = 0; l < k; ++l)
		buf[k - 1 - l] = "ACGT"[(y[1]>>l&1)<<1 | (y[0]>>l&1)];
	buf[k] = 0;
	return buf;
}

#endif
