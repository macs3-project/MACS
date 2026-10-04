# cython: language_level=3

# Return the free memory held by malloc to the system. Used by
# MACS3.Signal.CallPeakUnit before it forks worker processes: fork()
# copies the page-table entry of every resident page of the parent, and
# each worker unmaps them all again when it exits, so memory that was
# freed but is still resident costs time twice per worker. Only glibc
# has malloc_trim; elsewhere this does nothing.

cdef extern from *:
	"""
	#include <stdlib.h>
	#if defined(__GLIBC__)
	#include <malloc.h>
	static int macs3_release_free_heap(void) { return malloc_trim(0); }
	#else
	static int macs3_release_free_heap(void) { return 0; }
	#endif
	"""
	int macs3_release_free_heap() nogil

# Ask the kernel to back the ``n`` bytes at ``p``, a fresh allocation
# not yet written, with transparent huge pages, as numpy does for its
# own arrays of 4 MB or more: fewer page faults when it is first
# written and fewer TLB misses when it is read at random. Only the
# whole pages inside the block are advised; a smaller block, a failed
# madvise and any system without MADV_HUGEPAGE leave it as it was.

cdef extern from *:
	"""
	#if defined(__linux__)
	#include <sys/mman.h>
	#include <unistd.h>
	#endif
	static void macs3_advise_huge_pages(void *p, size_t n) {
	#if defined(__linux__) && defined(MADV_HUGEPAGE)
		long page = sysconf(_SC_PAGESIZE);
		size_t start, end;
		if (p == NULL || n < ((size_t)1 << 22) || page <= 0)
			return;
		start = ((size_t)p + (size_t)page - 1) & ~((size_t)page - 1);
		end = ((size_t)p + n) & ~((size_t)page - 1);
		if (end > start)
			madvise((void *)start, end - start, MADV_HUGEPAGE);
	#else
		(void)p;
		(void)n;
	#endif
	}
	"""
	void macs3_advise_huge_pages(void *p, size_t n) nogil
