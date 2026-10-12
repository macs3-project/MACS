# cython: language_level=3

# A prefetch hint, for C loops that know which memory they will touch
# some iterations ahead (MACS3.Signal.CallPeakUnit fills a hash table
# far larger than the cache, one random slot per key). It changes no
# result; where the compiler has no __builtin_prefetch it does nothing.

cdef extern from *:
	"""
	#if defined(__GNUC__) || defined(__clang__)
	#define macs3_prefetch_w(p) __builtin_prefetch((p), 1, 1)
	#else
	#define macs3_prefetch_w(p) ((void)(p))
	#endif
	"""
	void macs3_prefetch_w(const void *p) nogil
