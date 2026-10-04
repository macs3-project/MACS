# cython: language_level=3

# Text for the narrow-peak writers in MACS3.IO.PeakIO, formatted in C.
#
# Each function writes the same characters the Python formatting it
# replaces produces, for every value it accepts, and returns -1 (having
# written nothing that is kept) for a value it does not, so the caller
# can fall back to the Python formatting:
#
#   m3_g6_f32(f, out)     "%.6g" % f, for a float32 value f. Python formats
#                         it with dtoa: correctly rounded, ties to even.
#                         Here the rounding is done in exact 64-bit
#                         integer arithmetic on f = m * 2**e. Refused:
#                         inf, nan, subnormals, and values the integers
#                         cannot hold (below about 1e-7, at or above 1e19).
#   m3_r2g6_f32(f, out)   "%.6g" % round(f, 2), for a float32 value f.
#                         round() rounds the exact value to two decimals
#                         (dtoa, ties to even) and converts that decimal R
#                         back to the nearest double; below 10000, R has at
#                         most six significant digits, so "%.6g" gives R
#                         itself with trailing zeros dropped. Refused: f
#                         negative or -0.0, inf, nan, subnormal, and R of
#                         10000 or more.
#   m3_put_i64(p, v)      str(v) for an integer.
#   m3_put_letters(p, i)  subpeak_letters(i) in PeakIO.
#
# The first two write at most 13 characters; the last two return the end
# of what they wrote.

# const char *, spelled so that Cython 3.0 also accepts it in PeakIO's
# pure-Python annotations (cython.p_const_char needs Cython 3.1).
ctypedef const char *const_char_p

cdef extern from *:
	"""
	#include <stdint.h>
	#include <string.h>

	static const uint64_t m3_p10[19] = {
		1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL,
		10000000ULL, 100000000ULL, 1000000000ULL, 10000000000ULL,
		100000000000ULL, 1000000000000ULL, 10000000000000ULL,
		100000000000000ULL, 1000000000000000ULL, 10000000000000000ULL,
		100000000000000000ULL, 1000000000000000000ULL};

	static inline char *m3_put_u64(char *p, uint64_t v)
	{
		char tmp[20];
		int n = 0;
		do {
			tmp[n++] = (char)('0' + (int)(v % 10u));
			v /= 10u;
		} while (v);
		while (n)
			*p++ = tmp[--n];
		return p;
	}

	static inline char *m3_put_i64(char *p, int64_t v)
	{
		if (v < 0) {
			*p++ = '-';
			return m3_put_u64(p, (uint64_t)0 - (uint64_t)v);
		}
		return m3_put_u64(p, (uint64_t)v);
	}

	static char *m3_put_letters(char *p, int64_t i)
	{
		if (i < 26) {
			*p++ = (char)(97 + i);
			return p;
		}
		p = m3_put_letters(p, i / 26);
		*p++ = (char)(97 + i % 26);
		return p;
	}

	/* m * 2**e * 10**q as the fraction A / B of two integers, with
	   B <= 2**62 so that 2 * (A % B) cannot overflow; 0 when they do not
	   fit. m < 2**24. */
	static inline int m3_scale(uint64_t m, int e, int q, uint64_t *A,
	                           uint64_t *B)
	{
		if (q >= 0) {
			uint64_t t;
			if (q > 12)
				return 0;
			t = m * m3_p10[q];             /* < 2**24 * 10**12 < 2**64 */
			if (e >= 0) {
				if (e > 62 || t > (UINT64_MAX >> e))
					return 0;
				*A = t << e;
				*B = 1;
			} else {
				if (-e > 62)
					return 0;
				*A = t;
				*B = (uint64_t)1 << -e;
			}
		} else {
			uint64_t D;
			if (-q > 18)
				return 0;
			D = m3_p10[-q];
			if (e >= 0) {
				if (e > 39)                /* m << e < 2**63 */
					return 0;
				*A = m << e;
				*B = D;
			} else {
				if (-e > 62 || D > ((UINT64_MAX >> 2) >> -e))
					return 0;
				*A = m;
				*B = D << -e;
			}
		}
		return 1;
	}

	static inline int m3_g6_f32(float f, char *out)
	{
		uint32_t bits, be, frac;
		uint64_t m, A, B, F, N, r;
		int e, k, t, i, it, nd, x;
		char d[6];
		char *p = out;

		memcpy(&bits, &f, sizeof bits);
		be = (bits >> 23) & 0xffu;
		frac = bits & 0x7fffffu;
		if (be == 0xffu)
			return -1;                      /* inf, nan */
		if (be == 0 && frac != 0)
			return -1;                      /* subnormal */
		if (bits >> 31)
			*p++ = '-';
		if (be == 0) {                      /* 0.0, -0.0 */
			*p++ = '0';
			return (int)(p - out);
		}
		m = frac | 0x800000u;
		e = (int)be - 150;                  /* |f| = m * 2**e */

		/* 2**(e+23) <= |f| < 2**(e+24), so floor(log10|f|) is
		   floor((e+23) * log10(2)) or one more; 78913 / 2**18 gives
		   that floor exactly for |e+23| < 1650 */
		t = (e + 23) * 78913;
		k = t >= 0 ? t >> 18 : -((-t + (1 << 18) - 1) >> 18);
		for (it = 0; ; it++) {
			if (it == 3 || !m3_scale(m, e, 5 - k, &A, &B))
				return -1;
			F = A / B;
			if (F >= 1000000u)
				k++;
			else if (F < 100000u)
				k--;
			else
				break;
		}
		/* |f| * 10**(5-k) = F + r / B, with 10**5 <= F < 10**6 */
		r = A - F * B;
		N = F;
		if (2 * r > B || (2 * r == B && (F & 1u)))
			N++;
		if (N == 1000000u) {
			N = 100000u;
			k++;
		}
		for (i = 5; i >= 0; i--) {
			d[i] = (char)('0' + (int)(N % 10u));
			N /= 10u;
		}
		nd = 6;
		while (nd > 1 && d[nd - 1] == '0')
			nd--;
		if (k < -4 || k >= 6) {
			*p++ = d[0];
			if (nd > 1) {
				*p++ = '.';
				for (i = 1; i < nd; i++)
					*p++ = d[i];
			}
			*p++ = 'e';
			x = k;
			if (x < 0) {
				*p++ = '-';
				x = -x;
			} else
				*p++ = '+';
			if (x < 10)
				*p++ = '0';
			p = m3_put_u64(p, (uint64_t)x);
		} else if (k >= 0) {
			for (i = 0; i <= k; i++)
				*p++ = d[i];
			if (nd > k + 1) {
				*p++ = '.';
				for (i = k + 1; i < nd; i++)
					*p++ = d[i];
			}
		} else {
			*p++ = '0';
			*p++ = '.';
			for (i = 0; i < -k - 1; i++)
				*p++ = '0';
			for (i = 0; i < nd; i++)
				*p++ = d[i];
		}
		return (int)(p - out);
	}

	static inline int m3_r2g6_f32(float f, char *out)
	{
		uint32_t bits, be, frac, fr;
		uint64_t m, N, A, F, r, h;
		int e, s;
		char *p = out;

		memcpy(&bits, &f, sizeof bits);
		be = (bits >> 23) & 0xffu;
		frac = bits & 0x7fffffu;
		if ((bits >> 31) || be == 0xffu)
			return -1;                      /* negative, -0.0, inf, nan */
		if (be == 0) {
			if (frac != 0)
				return -1;                  /* subnormal */
			*p++ = '0';
			return 1;
		}
		m = frac | 0x800000u;
		e = (int)be - 150;
		if (e >= 0) {
			if (e > 10)
				return -1;
			N = (m << e) * 100u;
		} else {
			s = -e;
			A = m * 100u;                   /* < 2**31 */
			if (s > 40)
				N = 0;                      /* f * 100 < 2**-10 */
			else {
				F = A >> s;
				r = A & (((uint64_t)1 << s) - 1);
				h = (uint64_t)1 << (s - 1);
				N = F;
				if (r > h || (r == h && (F & 1u)))
					N++;
			}
		}
		if (N >= 1000000u)
			return -1;
		p = m3_put_u64(p, N / 100u);
		fr = (uint32_t)(N % 100u);
		if (fr) {
			*p++ = '.';
			*p++ = (char)('0' + fr / 10u);
			if (fr % 10u)
				*p++ = (char)('0' + fr % 10u);
		}
		return (int)(p - out);
	}
	"""
	char *m3_put_i64(char *p, long long v) nogil
	char *m3_put_letters(char *p, long long i) nogil
	int m3_g6_f32(float f, char *out) nogil
	int m3_r2g6_f32(float f, char *out) nogil
