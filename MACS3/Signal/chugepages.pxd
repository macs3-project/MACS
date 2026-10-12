# cython: language_level=3

# A numpy memory handler (PyDataMem_SetHandler) for arrays that grow by
# resize, as the per-chromosome arrays of
# MACS3.Signal.PairedEndTrack.PETrackI do while a file is read. An
# array below one transparent huge page is a malloc block, as it would
# be with numpy's own handler. Once a resize takes it to a huge page
# or more, it moves (once) into an anonymous mapping of
# its own: aligned to a huge page, advised MADV_HUGEPAGE, and reserving
# more address space than it uses (a quarter and a huge page more), so
# that later growth stays in place and is faulted in a huge page at a
# time. Growth past the reservation moves the mapping with mremap,
# without copying. A resize to a smaller size unmaps the whole pages
# past the new end, which ends the reservation. Where address space is
# limited (RLIMIT_AS) and a reservation cannot be mapped, a block maps
# only what it needs.
#
# Why: the callpeak pools fork up to 8 workers one after another, and
# fork() copies one page-table entry per resident 4 KB page of the
# parent but one per huge page, so a track held in huge pages costs
# each fork a fraction of the time. Memory freed into malloc's heap and
# reused is already faulted in 4 KB pages, so only fresh memory advised
# before its first touch gets huge pages.
#
# Each block starts with a 64-byte header (kind, bytes reserved, bytes
# in use) that the data follows, so that realloc and free know what
# they hold. Fresh anonymous memory reads as zeros, so calloc of a
# mapped block maps as malloc does.
#
# macs3_hp_available() is 1 only on Linux with MADV_HUGEPAGE and
# transparent huge pages set to "always" or "madvise"; elsewhere the
# handler is never installed, and the other functions are not called.

cdef extern from *:
	"""
	#include <stddef.h>
	#include <stdint.h>
	#include <stdlib.h>
	#include <string.h>
	#include "numpy/ndarraytypes.h"
	#if defined(__linux__)
	#include <sys/mman.h>
	#include <fcntl.h>
	#include <unistd.h>
	#endif

	#if defined(__linux__) && defined(MADV_HUGEPAGE) && defined(MREMAP_MAYMOVE)
	#define MACS3_HP_HEAD ((size_t)64)
	#define MACS3_HP_SMALL ((size_t)0x6d61637333480001ULL)
	#define MACS3_HP_MAPPED ((size_t)0x6d61637333480002ULL)
	/* the largest huge page size this is used with */
	#define MACS3_HP_MAX_HUGE ((size_t)32 << 20)

	typedef struct {
		size_t kind;
		size_t reserved;   /* SMALL: bytes malloc'd; MAPPED: bytes mapped */
		size_t used;       /* bytes of data, the size numpy last asked for */
	} macs3_hp_header;

	static size_t macs3_hp_page = 0;
	static size_t macs3_hp_huge = 0;
	static int macs3_hp_state = -1;

	static size_t macs3_hp_read_sys(const char *path, char *buf, size_t cap)
	{
		ssize_t n;
		int fd = open(path, O_RDONLY);
		if (fd < 0)
			return 0;
		n = read(fd, buf, cap - 1);
		close(fd);
		if (n <= 0)
			return 0;
		buf[n] = 0;
		return (size_t)n;
	}

	static int macs3_hp_available(void)
	{
		char buf[256];
		long page;
		unsigned long huge = 0;

		if (macs3_hp_state >= 0)
			return macs3_hp_state;
		macs3_hp_state = 0;
		page = sysconf(_SC_PAGESIZE);
		if (!macs3_hp_read_sys("/sys/kernel/mm/transparent_hugepage/enabled",
				       buf, sizeof(buf)))
			return 0;
		if (!strstr(buf, "[always]") && !strstr(buf, "[madvise]"))
			return 0;
		if (!macs3_hp_read_sys("/sys/kernel/mm/transparent_hugepage/hpage_pmd_size",
				       buf, sizeof(buf)))
			return 0;
		huge = strtoul(buf, NULL, 10);
		if (page <= 0 || (page & (page - 1)) != 0 || huge < (unsigned long)page ||
		    (huge & (huge - 1)) != 0 || huge > MACS3_HP_MAX_HUGE)
			return 0;
		macs3_hp_page = (size_t)page;
		macs3_hp_huge = (size_t)huge;
		macs3_hp_state = 1;
		return 1;
	}

	/* sizes above this are refused, so no sum below can wrap */
	#define MACS3_HP_MAX_SIZE (SIZE_MAX / 4)

	static size_t macs3_hp_round(size_t n, size_t unit)
	{
		return (n + unit - 1) & ~(unit - 1);
	}

	static void *macs3_hp_small(size_t size, int zero)
	{
		macs3_hp_header *h;
		size_t n = size + MACS3_HP_HEAD;

		h = (macs3_hp_header *)(zero ? calloc(1, n) : malloc(n));
		if (h == NULL)
			return NULL;
		h->kind = MACS3_HP_SMALL;
		h->reserved = n;
		h->used = size;
		return (char *)h + MACS3_HP_HEAD;
	}

	/* a mapping of `reserve` bytes (a multiple of the page size),
	   advised huge pages, starting on a huge-page boundary unless
	   `align` is 0 */
	static char *macs3_hp_map(size_t reserve, int align)
	{
		size_t span = reserve + (align ? macs3_hp_huge : 0);
		char *raw, *p, *end;

		raw = (char *)mmap(NULL, span, PROT_READ | PROT_WRITE,
				   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
		if (raw == (char *)MAP_FAILED)
			return NULL;
		p = align ? (char *)macs3_hp_round((uintptr_t)raw, macs3_hp_huge) : raw;
		end = raw + span;
		if (p > raw)
			munmap(raw, (size_t)(p - raw));
		if (end > p + reserve)
			munmap(p + reserve, (size_t)(end - (p + reserve)));
		madvise(p, reserve, MADV_HUGEPAGE);
		return p;
	}

	/* address space to reserve for a block of `size` bytes: what it
	   needs, a quarter more and one huge page more, in whole huge
	   pages, so that growth by resize stays in place for a while and
	   the pages it touches lie wholly inside the mapping, which the
	   kernel needs to fault in a huge page */
	static size_t macs3_hp_reserve_for(size_t size)
	{
		size_t need = size + MACS3_HP_HEAD;
		return macs3_hp_round(need + need / 4 + macs3_hp_huge, macs3_hp_huge);
	}

	static void *macs3_hp_mapped(size_t size)
	{
		macs3_hp_header *h;
		size_t reserve = macs3_hp_reserve_for(size);

		h = (macs3_hp_header *)macs3_hp_map(reserve, 1);
		if (h == NULL) {
			/* where address space is limited (RLIMIT_AS), no more
			   than the block needs */
			reserve = macs3_hp_round(size + MACS3_HP_HEAD, macs3_hp_page);
			h = (macs3_hp_header *)macs3_hp_map(reserve, 0);
			if (h == NULL)
				return NULL;
		}
		h->kind = MACS3_HP_MAPPED;
		h->reserved = reserve;
		h->used = size;
		return (char *)h + MACS3_HP_HEAD;
	}

	static void *macs3_hp_new(size_t size, int zero)
	{
		if (size > MACS3_HP_MAX_SIZE)
			return NULL;
		if (size + MACS3_HP_HEAD < macs3_hp_huge)
			return macs3_hp_small(size, zero);
		return macs3_hp_mapped(size);
	}

	static void *macs3_hp_malloc(void *ctx, size_t size)
	{
		(void)ctx;
		return macs3_hp_new(size, 0);
	}

	static void *macs3_hp_calloc(void *ctx, size_t nelem, size_t elsize)
	{
		(void)ctx;
		if (elsize != 0 && nelem > SIZE_MAX / elsize)
			return NULL;
		return macs3_hp_new(nelem * elsize, 1);
	}

	static void macs3_hp_free(void *ctx, void *ptr, size_t size)
	{
		macs3_hp_header *h;

		(void)ctx;
		(void)size;
		if (ptr == NULL)
			return;
		h = (macs3_hp_header *)((char *)ptr - MACS3_HP_HEAD);
		if (h->kind == MACS3_HP_MAPPED)
			munmap((void *)h, h->reserved);
		else
			free(h);
	}

	static void *macs3_hp_realloc(void *ctx, void *ptr, size_t new_size)
	{
		macs3_hp_header *h, *g;
		char *q;
		size_t keep, reserve;

		(void)ctx;
		if (ptr == NULL)
			return macs3_hp_new(new_size, 0);
		if (new_size > MACS3_HP_MAX_SIZE)
			return NULL;
		h = (macs3_hp_header *)((char *)ptr - MACS3_HP_HEAD);
		if (h->kind == MACS3_HP_SMALL) {
			if (new_size + MACS3_HP_HEAD < macs3_hp_huge) {
				g = (macs3_hp_header *)realloc(h, new_size + MACS3_HP_HEAD);
				if (g == NULL)
					return NULL;
				g->reserved = new_size + MACS3_HP_HEAD;
				g->used = new_size;
				return (char *)g + MACS3_HP_HEAD;
			}
			/* grown to a huge page or more: into a mapping, once */
			q = (char *)macs3_hp_mapped(new_size);
			if (q == NULL)
				return NULL;
			memcpy(q, ptr, h->used < new_size ? h->used : new_size);
			free(h);
			return q;
		}
		if (new_size < h->used) {
			/* shrinking ends the reservation: unmap whole pages past
			   the new end */
			keep = macs3_hp_round(new_size + MACS3_HP_HEAD, macs3_hp_page);
			if (keep < h->reserved) {
				munmap((char *)h + keep, h->reserved - keep);
				h->reserved = keep;
			}
			h->used = new_size;
			return ptr;
		}
		if (new_size + MACS3_HP_HEAD <= h->reserved) {
			h->used = new_size;
			return ptr;
		}
		/* past the reservation: move the mapping, page tables and all,
		   to a new reservation, or to just what it needs where address
		   space is limited, or else map anew and copy */
		reserve = macs3_hp_reserve_for(new_size);
		g = (macs3_hp_header *)mremap((void *)h, h->reserved, reserve, MREMAP_MAYMOVE);
		if (g == (macs3_hp_header *)MAP_FAILED) {
			reserve = macs3_hp_round(new_size + MACS3_HP_HEAD, macs3_hp_page);
			g = (macs3_hp_header *)mremap((void *)h, h->reserved, reserve,
						      MREMAP_MAYMOVE);
		}
		if (g != (macs3_hp_header *)MAP_FAILED) {
			madvise((void *)g, reserve, MADV_HUGEPAGE);
			g->reserved = reserve;
			g->used = new_size;
			return (char *)g + MACS3_HP_HEAD;
		}
		q = (char *)macs3_hp_mapped(new_size);
		if (q == NULL)
			return NULL;
		memcpy(q, ptr, h->used);
		munmap((void *)h, h->reserved);
		return q;
	}

	static PyDataMem_Handler macs3_hp_handler = {
		"macs3_huge_pages",
		1,
		{NULL, macs3_hp_malloc, macs3_hp_calloc, macs3_hp_realloc, macs3_hp_free}
	};

	static PyObject *macs3_hp_capsule(void)
	{
		return PyCapsule_New((void *)&macs3_hp_handler, "mem_handler", NULL);
	}
	#else
	static int macs3_hp_available(void) { return 0; }
	static PyObject *macs3_hp_capsule(void) { Py_RETURN_NONE; }
	#endif
	"""
	int macs3_hp_available()
	object macs3_hp_capsule()
