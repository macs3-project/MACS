#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fml.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

bseq1_t *bseq_read(const char *fn, int *n_)
{
	gzFile fp;
	bseq1_t *seqs;
	kseq_t *ks;
	int m, n;
	uint64_t size = 0;

	*n_ = 0;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = kseq_init(fp);

	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->l_seq = ks->seq.l;
		size += seqs[n++].l_seq;
	}
	*n_ = n;

	kseq_destroy(ks);
	gzclose(fp);
	return seqs;
}

void seq_reverse(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		s[l-1-i] = s[i]; s[i] = tmp;
	}
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}
