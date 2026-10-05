# cython: language_level=3
# Time-stamp: <2025-11-14 19:08:47 Tao Liu>

"""Module for filter duplicate tags from paired-end data

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# Python modules
# ------------------------------------
import io
import sys
from array import array as pyarray
from collections import Counter,defaultdict
from operator import itemgetter
# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.Pileup import se_all_in_one_pileup_max3
from MACS3.Signal.BedGraph import (bedGraphTrackI,
                                   bedGraphTrackII)
from MACS3.Signal.PileupV2 import (LR_DTYPE,
                                   LRC_DTYPE,
                                   pileup_from_LR_hmmratac,
                                   pileup_from_LRC,
                                   pileup_from_LRC_as_list,
                                   pileup_from_LRC_centers_as_list,
                                   pileup_from_LR_as_list,
                                   pileup_from_PN_shifted,
                                   pileup_LRC_as_list_equal,
                                   pileup_LRC_centers_as_list_equal,
                                   over_two_pv_array)
from MACS3.Signal.Region import Regions
# ------------------------------------
# Other modules
# ------------------------------------
import cython
import numpy as np
import cython.cimports.numpy as cnp
from cython.cimports.cpython import bool
from cython.cimports.libc.stdint import INT32_MAX as INT_MAX
from cython.cimports.libc.stdint import UINT32_MAX
from cython.cimports.libc.string import memcpy
from cython.cimports.MACS3.Signal.crefcount import Py_REFCNT
from cython.cimports.MACS3.Signal.chugepages import (macs3_hp_available,
                                                     macs3_hp_capsule)

from MACS3.Utilities.Logger import logging
from MACS3.Utilities.Threads import MIN_THREADED_SIZE, map_chromosomes

logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

# finalize and filter_dup hand a chromosome of at least this many
# fragments to map_chromosomes, and run any other on the spot
_MIN_THREADED_SIZE = cython.declare(cython.long, MIN_THREADED_SIZE)
# the numpy memory handler of the per-chromosome arrays
# (chugepages.pxd): one that reaches a huge page moves into memory of
# its own advised to use huge pages; None, and numpy's own, where
# transparent huge pages are not available
_TRACK_MEMORY = macs3_hp_capsule() if macs3_hp_available() else None

# Let numpy enforce PE-ness using ndarray, gives bonus speedup when sorting
# PE data doesn't have strandedness

# ------------------------------------
# C fast paths for the track methods
# ------------------------------------
# Each runs only on a plain array, one whose layout a C loop can read
# through a pointer (is_plain_lr_array, is_plain_lrc_array). sort_lr_array
# and argsort_lrc check that themselves and PETrackI.filter_dup checks it,
# and that the array owns its data and nothing else holds it, before
# calling filter_dup_lr, which filters in place; any other array takes
# the original numpy
# or Python code, so it gives the original result or the original
# exception. LR_DTYPE and LRC_DTYPE, the record layouts of PETrackI and
# PETrackII chromosomes, come from PileupV2.


@cython.cfunc
def is_plain_lr_array(locs: cnp.ndarray) -> cython.bint:
    """Whether `locs` is a writable C-contiguous array of LR_DTYPE, the
    layout sort_lr_array and filter_dup_lr read and write through a
    pointer."""
    return (locs.dtype == LR_DTYPE and locs.flags.c_contiguous and
            locs.flags.writeable)


@cython.cfunc
@cython.nogil
@cython.inline
@cython.exceptval(check=False)
def lr_key(l: cython.int, r: cython.int) -> cython.ulonglong:
    """(l, r) as one uint64 key, (l ^ 2^31) * 2^32 + (r ^ 2^31) with l
    and r taken as uint32. Flipping the sign bit maps the signed int32
    order onto the unsigned order, so the keys of two records compare
    as the records do: by l, then by r, both ascending."""
    sign: cython.uint = 0x80000000
    return ((cython.cast(cython.ulonglong,
                         cython.cast(cython.uint, l) ^ sign) << 32) |
            (cython.cast(cython.uint, r) ^ sign))


@cython.cfunc
@cython.nogil
@cython.inline
@cython.exceptval(check=False)
def lr_key_l(key: cython.ulonglong) -> cython.int:
    """The l of an lr_key."""
    sign: cython.uint = 0x80000000
    return cython.cast(cython.int, cython.cast(cython.uint, key >> 32) ^ sign)


@cython.cfunc
@cython.nogil
@cython.inline
@cython.exceptval(check=False)
def lr_key_r(key: cython.ulonglong) -> cython.int:
    """The r of an lr_key."""
    sign: cython.uint = 0x80000000
    return cython.cast(cython.int, cython.cast(cython.uint, key) ^ sign)


@cython.cfunc
@cython.nogil
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def lr_records_to_keys(q: cython.pointer(cython.char),
                       n: cython.Py_ssize_t) -> cython.void:
    """Overwrite each of the n 8-byte (l, r) records at q with its
    lr_key, read and written through memcpy, so q needs no alignment
    and no access aliases another type."""
    i: cython.Py_ssize_t
    lr: cython.int[2]
    key: cython.ulonglong

    for i in range(n):
        memcpy(lr, q + 8 * i, 8)
        key = lr_key(lr[0], lr[1])
        memcpy(q + 8 * i, cython.address(key), 8)


@cython.cfunc
@cython.nogil
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def lr_keys_to_records(q: cython.pointer(cython.char),
                       n: cython.Py_ssize_t) -> cython.void:
    """The inverse of lr_records_to_keys: overwrite each of the n
    lr_keys at q with its (l, r) record."""
    i: cython.Py_ssize_t
    lr: cython.int[2]
    key: cython.ulonglong

    for i in range(n):
        memcpy(cython.address(key), q + 8 * i, 8)
        lr[0] = lr_key_l(key)
        lr[1] = lr_key_r(key)
        memcpy(q + 8 * i, lr, 8)


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def sort_lr_array(locs: cnp.ndarray):
    """Sort a PETrackI location array in place, as
    `locs.sort(order=['l', 'r'])` does, but through an integer sort of
    the records' lr_keys.

    A record holds nothing but (l, r), so records with equal keys are
    identical, and any sort of the keys gives the same bytes as the
    structured sort. An array that is not plain (is_plain_lr_array)
    takes the structured sort itself. A one-dimensional plain array
    is turned into its keys in place, an LR record being 8 bytes, and
    sorted as a uint64 view of itself, so nothing is allocated; any
    other plain array sorts a separate array of keys. The plain path
    runs without the GIL but for the view or the allocation.
    """
    n: cython.long
    i: cython.long
    keys: cnp.ndarray
    p: cython.pointer(cython.int)
    k: cython.pointer(cython.ulonglong)
    q: cython.pointer(cython.char)

    if not is_plain_lr_array(locs):
        locs.sort(order=['l', 'r'])
        return
    n = locs.shape[0]
    if n < 2:
        return
    if locs.ndim == 1:
        keys = locs.view(dtype=np.uint64, type=np.ndarray)
        q = locs.data
        with cython.nogil:
            lr_records_to_keys(q, n)
        try:
            keys.sort()
        finally:
            # back to records, whether or not the sort raised
            with cython.nogil:
                lr_keys_to_records(q, n)
        return
    p = cython.cast(cython.pointer(cython.int), locs.data)
    keys = np.empty(n, dtype=np.uint64)
    k = cython.cast(cython.pointer(cython.ulonglong), keys.data)
    with cython.nogil:
        for i in range(n):
            k[i] = lr_key(p[2 * i], p[2 * i + 1])
    keys.sort()
    with cython.nogil:
        for i in range(n):
            p[2 * i] = lr_key_l(k[i])
            p[2 * i + 1] = lr_key_r(k[i])
    return


@cython.cfunc
@cython.nogil
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def filter_dup_lr(p: cython.pointer(cython.int), size: cython.Py_ssize_t,
                  maxnum: cython.int,
                  kept: cython.pointer(cython.Py_ssize_t)) -> cython.ulonglong:
    """The loop of PETrackI.filter_dup, in C and in place, over the
    `size` sorted (l, r) records at `p`: drop each record after the
    first `maxnum` of a run of equal records, moving the kept ones
    forward in their order, set kept[0] to how many were kept, and
    return the sum of the dropped records' r - l, each taken as uint64.

    The first record of a run is always kept, even at `maxnum` 0, and
    the comparisons are the original loop's, on the same int32 values.
    Each write lands at or before the record just read, so no record is
    overwritten before it is read. The caller subtracts the sum from
    the unsigned `length`, which wraps exactly as subtracting the
    lengths one by one does. `size` must be at least 2.
    """
    i: cython.Py_ssize_t
    j: cython.Py_ssize_t = 1
    n: cython.int = 1
    loc_start: cython.int
    loc_end: cython.int
    current_loc_start: cython.int = p[0]
    current_loc_end: cython.int = p[1]
    removed_length: cython.ulonglong = 0

    for i in range(1, size):
        loc_start = p[2 * i]
        loc_end = p[2 * i + 1]
        if loc_start != current_loc_start or loc_end != current_loc_end:
            current_loc_start = loc_start
            current_loc_end = loc_end
            n = 1
        else:
            n += 1
            if n > maxnum:
                removed_length += cython.cast(cython.ulonglong, current_loc_end - current_loc_start)
                continue
        p[2 * j] = loc_start
        p[2 * j + 1] = loc_end
        j += 1
    kept[0] = j
    return removed_length


@cython.cfunc
def is_plain_lrc_array(locs: cnp.ndarray) -> cython.bint:
    """Whether `locs` is a one-dimensional C-contiguous array of
    LRC_DTYPE, the layout argsort_lrc reads through a pointer."""
    return (locs.ndim == 1 and locs.dtype == LRC_DTYPE and
            locs.flags.c_contiguous)


# One element of argsort_lrc's work array: the record's lr_key, its
# count, and its index in the array being sorted.
LRCItem = cython.struct(k=cython.ulonglong, c=cython.uint, i=cython.uint)


@cython.cfunc
@cython.nogil
@cython.inline
@cython.exceptval(check=False)
def lrc_lt(a: LRCItem, b: LRCItem) -> cython.bint:
    """Whether a's record compares below b's: by l, then r, then c.

    This is the comparison numpy's structured argsort makes for
    `order=['l', 'r']` on LRC_DTYPE, since numpy breaks ties on the
    fields left out of `order` in dtype order, here c.
    """
    return a.k < b.k or (a.k == b.k and a.c < b.c)


# lrc_aheapsort and lrc_aquicksort below are translations of NumPy's
# generic heapsort and introsort (numpy/_core/src/npysort/heapsort.cpp
# and quicksort_generic.cpp), used under NumPy's licence:
#
# Copyright (c) 2005-2025, NumPy Developers.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.
#
#     * Neither the name of the NumPy Developers nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@cython.cfunc
@cython.nogil
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def lrc_aheapsort(a: cython.pointer(LRCItem), n: cython.Py_ssize_t) -> cython.void:
    """numpy's generic `npy_aheapsort` (numpy/_core/src/npysort/
    heapsort.cpp, numpy 2.5.3) on the n items at a, with numpy's
    1-based a[x] written as a[x - 1]. numpy moves indices into the
    data; this moves the items, which carry their index along.
    """
    i: cython.Py_ssize_t
    j: cython.Py_ssize_t
    l: cython.Py_ssize_t
    tmp: LRCItem

    l = n >> 1
    while l > 0:
        tmp = a[l - 1]
        i = l
        j = l << 1
        while j <= n:
            if j < n and lrc_lt(a[j - 1], a[j]):
                j += 1
            if lrc_lt(tmp, a[j - 1]):
                a[i - 1] = a[j - 1]
                i = j
                j += j
            else:
                break
        a[i - 1] = tmp
        l -= 1

    while n > 1:
        tmp = a[n - 1]
        a[n - 1] = a[0]
        n -= 1
        i = 1
        j = 2
        while j <= n:
            if j < n and lrc_lt(a[j - 1], a[j]):
                j += 1
            if lrc_lt(tmp, a[j - 1]):
                a[i - 1] = a[j - 1]
                i = j
                j += j
            else:
                break
        a[i - 1] = tmp
    return


@cython.cfunc
@cython.nogil
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def lrc_aquicksort(v: cython.pointer(LRCItem), num: cython.Py_ssize_t) -> cython.void:
    """numpy's generic `npy_aquicksort_impl` (numpy/_core/src/npysort/
    quicksort_generic.cpp, numpy 2.5.3), the introsort numpy runs for
    an argsort of a structured array: median-of-three pivot, insertion
    sort below 16 items, heapsort past a depth of 2 * floor(log2(num)).

    numpy permutes an index array and compares the records the indices
    point to; this permutes the items, which carry the index along, so
    every comparison sees the same two records and the same branch is
    taken at every step. The indices end in numpy's order, ties
    included.
    """
    pl: cython.Py_ssize_t = 0
    pr: cython.Py_ssize_t = num - 1
    pm: cython.Py_ssize_t
    pi: cython.Py_ssize_t
    pj: cython.Py_ssize_t
    pk: cython.Py_ssize_t
    stack: cython.Py_ssize_t[128]
    sptr: cython.int = 0
    depth: cython.int[128]
    psdepth: cython.int = 0
    cdepth: cython.int = 0
    u: cython.size_t
    vp: LRCItem
    tmp: LRCItem

    # npy_get_msb(num) * 2
    u = cython.cast(cython.size_t, num) >> 1
    while u:
        cdepth += 1
        u >>= 1
    cdepth *= 2

    while True:
        if cdepth < 0:
            lrc_aheapsort(v + pl, pr - pl + 1)
        else:
            while pr - pl > 15:
                # quicksort partition
                pm = pl + ((pr - pl) >> 1)
                if lrc_lt(v[pm], v[pl]):
                    tmp = v[pm]
                    v[pm] = v[pl]
                    v[pl] = tmp
                if lrc_lt(v[pr], v[pm]):
                    tmp = v[pr]
                    v[pr] = v[pm]
                    v[pm] = tmp
                if lrc_lt(v[pm], v[pl]):
                    tmp = v[pm]
                    v[pm] = v[pl]
                    v[pl] = tmp
                vp = v[pm]
                pi = pl
                pj = pr - 1
                tmp = v[pm]
                v[pm] = v[pj]
                v[pj] = tmp
                while True:
                    pi += 1
                    while lrc_lt(v[pi], vp) and pi < pj:
                        pi += 1
                    pj -= 1
                    while lrc_lt(vp, v[pj]) and pi < pj:
                        pj -= 1
                    if pi >= pj:
                        break
                    tmp = v[pi]
                    v[pi] = v[pj]
                    v[pj] = tmp
                pk = pr - 1
                tmp = v[pi]
                v[pi] = v[pk]
                v[pk] = tmp
                # push largest partition on stack
                if pi - pl < pr - pi:
                    stack[sptr] = pi + 1
                    stack[sptr + 1] = pr
                    sptr += 2
                    pr = pi - 1
                else:
                    stack[sptr] = pl
                    stack[sptr + 1] = pi - 1
                    sptr += 2
                    pl = pi + 1
                cdepth -= 1
                depth[psdepth] = cdepth
                psdepth += 1

            # insertion sort
            pi = pl + 1
            while pi <= pr:
                vp = v[pi]
                pj = pi
                pk = pi - 1
                while pj > pl and lrc_lt(vp, v[pk]):
                    v[pj] = v[pk]
                    pj -= 1
                    pk -= 1
                v[pj] = vp
                pi += 1
        # stack_pop
        if sptr == 0:
            break
        sptr -= 2
        pl = stack[sptr]
        pr = stack[sptr + 1]
        psdepth -= 1
        cdepth = depth[psdepth]
    return


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def argsort_lrc(locs: cnp.ndarray) -> cnp.ndarray:
    """Return `np.argsort(locs, order=['l', 'r'])` for a PETrackII
    location array, the same indices in the same order, through a typed
    comparison instead of numpy's per-field structured one.

    Records equal in l, r and c are ties, and numpy's default argsort
    is not stable, so their order is whatever its introsort leaves;
    their barcodes can differ. lrc_aquicksort makes the comparisons
    numpy makes, in the same sequence, so it leaves them in the same
    order. A plain array (is_plain_lrc_array) of 2 to 2^32 - 1 records
    takes this path, without the GIL but for its allocations; anything
    else goes to numpy's argsort itself.
    """
    n: cython.Py_ssize_t
    i: cython.Py_ssize_t
    work: cnp.ndarray
    out: cnp.ndarray
    v: cython.pointer(LRCItem)
    o: cython.pointer(cnp.npy_intp)
    rec: cython.pointer(cython.char)
    l: cython.int
    r: cython.int
    c: cython.ushort

    if not (is_plain_lrc_array(locs) and
            2 <= locs.shape[0] <= UINT32_MAX):
        return np.argsort(locs, order=['l', 'r'])
    n = locs.shape[0]

    work = np.empty((n, 2), dtype=np.uint64)
    v = cython.cast(cython.pointer(LRCItem), work.data)
    rec = cython.cast(cython.pointer(cython.char), locs.data)
    with cython.nogil:
        for i in range(n):
            # records are 10 bytes, so the fields are read unaligned
            memcpy(cython.address(l), rec, 4)
            memcpy(cython.address(r), rec + 4, 4)
            memcpy(cython.address(c), rec + 8, 2)
            rec += 10
            v[i].k = lr_key(l, r)
            v[i].c = c
            v[i].i = cython.cast(cython.uint, i)

        lrc_aquicksort(v, n)

    out = np.empty(n, dtype=np.intp)
    o = cython.cast(cython.pointer(cnp.npy_intp), out.data)
    with cython.nogil:
        for i in range(n):
            o[i] = v[i].i
    return out


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def _lrc_append(locs, bcs, i: cython.Py_ssize_t, n: cython.Py_ssize_t,
                starts, ends, counts, barcode_ids,
                dlength: cython.pointer(cython.longlong)) -> cython.int:
    """Write the n fragments of the arrays into records i to i + n - 1
    of the PETrackII location array ``locs`` and its barcode array
    ``bcs`` in one C loop, putting the sum of their C-int ``(end -
    start) * count`` (each wrapped as C int arithmetic wraps it) in
    ``dlength[0]``, and return 1. This is what add_loc_arrays' numpy
    assignments and sum give. Returns 0, writing nothing, unless every
    array is 1-D and C-contiguous with the dtypes FragParser passes:
    locs LRC_DTYPE, bcs, starts, ends and barcode_ids int32, counts
    uint16, each of n elements, and locs and bcs hold at least i + n
    records."""
    k: cython.Py_ssize_t
    sv: cython.int[::1]
    ev: cython.int[::1]
    cv: cython.ushort[::1]
    bv: cython.int[::1]
    ov: cython.int[::1]
    la: cnp.ndarray
    rec: cython.pointer(cython.char)
    s: cython.int
    e: cython.int
    c: cython.ushort
    dl: cython.longlong = 0

    if not (type(locs) is np.ndarray and type(bcs) is np.ndarray and
            type(starts) is np.ndarray and type(ends) is np.ndarray and
            type(counts) is np.ndarray and type(barcode_ids) is np.ndarray):
        return 0
    if not (locs.ndim == 1 and locs.dtype == LRC_DTYPE and
            locs.flags.c_contiguous and locs.shape[0] >= i + n and
            bcs.ndim == 1 and bcs.dtype == np.int32 and
            bcs.flags.c_contiguous and bcs.shape[0] >= i + n and
            starts.ndim == 1 and starts.dtype == np.int32 and
            starts.flags.c_contiguous and starts.shape[0] == n and
            ends.ndim == 1 and ends.dtype == np.int32 and
            ends.flags.c_contiguous and ends.shape[0] == n and
            counts.ndim == 1 and counts.dtype == np.uint16 and
            counts.flags.c_contiguous and counts.shape[0] == n and
            barcode_ids.ndim == 1 and barcode_ids.dtype == np.int32 and
            barcode_ids.flags.c_contiguous and barcode_ids.shape[0] == n):
        return 0
    sv = starts
    ev = ends
    cv = counts
    bv = barcode_ids
    ov = bcs
    la = locs
    # records are 10 bytes, so the fields are written unaligned
    rec = cython.cast(cython.pointer(cython.char), la.data) + 10 * i
    for k in range(n):
        s = sv[k]
        e = ev[k]
        c = cv[k]
        memcpy(rec, cython.address(s), 4)
        memcpy(rec + 4, cython.address(e), 4)
        memcpy(rec + 8, cython.address(c), 2)
        rec += 10
        ov[i + k] = bv[k]
        dl += cython.cast(cython.int,
                          (cython.cast(cython.uint, e) -
                           cython.cast(cython.uint, s)) *
                          cython.cast(cython.uint, c))
    dlength[0] = dl
    return 1


@cython.cfunc
def new_track_array(n, dtype) -> cnp.ndarray:
    """np.zeros(n, dtype=dtype), with its memory, and that of every
    resize of it, from _TRACK_MEMORY where there is one. Such an array
    owns its data like any other.

    The callpeak pools fork their workers from the process holding
    the track, and fork() copies a page-table entry for each 4 KB page
    of it but one for each huge page.

    PETrackI's arrays only: PETrackII's are replaced by sorted copies
    in finalize, so the arrays it reads into are gone before any
    fork, and its reserve already over-allocates (4x), to which the
    handler's reservation would add a quarter more address space."""
    if _TRACK_MEMORY is None:
        return np.zeros(n, dtype=dtype)
    old = cnp.PyDataMem_SetHandler(_TRACK_MEMORY)
    try:
        a = np.zeros(n, dtype=dtype)
    finally:
        cnp.PyDataMem_SetHandler(old)
    return a


@cython.cclass
class PETrackI:
    """In-memory paired-end fragment container grouped by chromosome.
    
    The track exposes utilities for sorting, filtering, downsampling, and pileup
    generation on numpy structured arrays of left/right coordinates.
    """
    locations = cython.declare(dict, visibility="public")
    size = cython.declare(dict, visibility="public")
    buf_size = cython.declare(dict, visibility="public")
    is_sorted = cython.declare(bool, visibility="public")
    total = cython.declare(cython.ulong, visibility="public")
    annotation = cython.declare(str, visibility="public")
    # rlengths: reference chromosome lengths dictionary
    rlengths = cython.declare(dict, visibility="public")
    buffer_size = cython.declare(cython.long, visibility="public")
    length = cython.declare(cython.ulonglong, visibility="public")
    average_template_length = cython.declare(cython.float, visibility="public")
    is_destroyed: bool

    def __init__(self, anno: str = "", buffer_size: cython.long = 100000):
        """Initialize an empty paired-end track.
        
        Parameters
        ----------
        anno : str, optional
            Annotation label retained with the track metadata.
        buffer_size : int, optional
            Number of fragment slots allocated per growth chunk for each chromosome.
        """
        # dictionary with chrname as key, nparray with
        # [('l','i4'),('r','i4')] as value
        self.locations = {}
        # dictionary with chrname as key, size of the above nparray as value
        # size is to remember the number of the fragments added to this chromosome
        self.size = {}
        # dictionary with chrname as key, size of the above nparray as value
        self.buf_size = {}
        self.is_sorted = False
        self.total = 0           # total fragments
        self.annotation = anno   # need to be figured out
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.average_template_length = 0.0
        self.is_destroyed = False

    @cython.ccall
    def add_loc(self, chromosome: bytes,
                start: cython.int, end: cython.int):
        """Append a paired-end fragment to the track.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) that owns the fragment.
        start : int
            Zero-based start coordinate of the fragment (5' end).
        end : int
            Zero-based end coordinate of the fragment (3' end).
        
        Notes
        -----
        Fragments are stored in structured numpy arrays keyed by chromosome, and
        the running fragment count and total template length are updated in place.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.add_loc(b"chr1", 50, 90)
        """
        i: cython.int

        if chromosome not in self.locations:
            self.buf_size[chromosome] = self.buffer_size
            # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            self.locations[chromosome] = new_track_array(self.buffer_size,
                                                         [('l', 'i4'), ('r', 'i4')])
            self.locations[chromosome][0] = (start, end)
            self.size[chromosome] = 1
        else:
            i = self.size[chromosome]
            if self.buf_size[chromosome] == i:
                self.buf_size[chromosome] += self.buffer_size
                self.locations[chromosome].resize((self.buf_size[chromosome]),
                                                  refcheck=False)
            self.locations[chromosome][i] = (start, end)
            self.size[chromosome] = i + 1
        self.length += end - start
        return

    @cython.ccall
    def add_loc_arrays(self, chromosome: bytes, starts, ends):
        """Append many paired-end fragments of one chromosome to the track.

        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) that owns the fragments.
        starts : numpy.ndarray
            int32 leftmost ends of the fragments, in the order to append.
        ends : numpy.ndarray
            int32 rightmost ends of the fragments, same length as ``starts``.

        Notes
        -----
        Leaves the track exactly as calling ``add_loc(chromosome,
        starts[k], ends[k])`` for each ``k`` in turn would: the same array,
        size, buffer size and total ``length``.
        """
        i: cython.long
        n: cython.long
        b: cython.long
        dlength: cython.longlong

        n = len(starts)
        if n == 0:
            return
        if self.buffer_size <= 0:
            # add_loc's own behaviour, errors included
            for k in range(n):
                self.add_loc(chromosome, starts[k], ends[k])
            return

        if chromosome not in self.locations:
            self.buf_size[chromosome] = self.buffer_size
            # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            self.locations[chromosome] = new_track_array(self.buffer_size,
                                                         [('l', 'i4'), ('r', 'i4')])
            self.size[chromosome] = 0
        i = self.size[chromosome]
        b = self.buf_size[chromosome]
        if i + n > b:
            # grow in steps of buffer_size, as add_loc does
            while b < i + n:
                b += self.buffer_size
            self.locations[chromosome].resize((b), refcheck=False)
            self.buf_size[chromosome] = b
        self.locations[chromosome]['l'][i:i + n] = starts
        self.locations[chromosome]['r'][i:i + n] = ends
        self.size[chromosome] = i + n
        # add_loc adds each int32 difference end - start to length
        dlength = np.subtract(ends, starts, dtype=np.int32).sum(dtype=np.int64)
        self.length += dlength
        return

    @cython.ccall
    def destroy(self):
        """Release numpy buffers held by the track.
        
        All per-chromosome arrays are resized to zero so the memory footprint is
        returned to the allocator, and the track is marked as destroyed.
        """
        chrs: set
        chromosome: bytes

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.locations:
                self.locations[chromosome].resize(self.buffer_size,
                                                  refcheck=False)
                self.locations[chromosome].resize(0,
                                                  refcheck=False)
                self.locations[chromosome] = None
                self.locations.pop(chromosome)
        self.is_destroyed = True
        return

    @cython.ccall
    def set_rlengths(self, rlengths: dict) -> bool:
        """Attach reference chromosome lengths to the track.
        
        Parameters
        ----------
        rlengths : dict
            Mapping from chromosome name (bytes) to reference length.
        
        Returns
        -------
        bool
            True when the length mapping has been updated.
        
        Notes
        -----
        Any chromosome stored in the track but missing from ``rlengths`` is
        assigned ``INT_MAX`` so downstream bounds checks can succeed.
        """
        valid_chroms: set
        missed_chroms: set
        chrom: bytes

        valid_chroms = set(self.locations.keys()).intersection(rlengths.keys())
        for chrom in sorted(valid_chroms):
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.locations.keys()).difference(rlengths.keys())
        for chrom in sorted(missed_chroms):
            self.rlengths[chrom] = INT_MAX
        return True

    @cython.ccall
    def get_rlengths(self) -> dict:
        """Return the reference chromosome lengths associated with the track.
        
        Returns
        -------
        dict
            Mapping from chromosome name (bytes) to reference length. Chromosomes
            without a recorded length default to ``INT_MAX``.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

    @cython.ccall
    def finalize(self):
        """Shrink backing arrays and sort fragments in place.
        
        Each per-chromosome array is resized to the observed fragment count, sorted
        by the left and right coordinates, and the aggregate counters ``total`` and
        ``average_template_length`` are refreshed. Call this after loading data.
        """
        c: bytes
        chrnames: set
        big: list
        sizes: list

        self.total = 0

        chrnames = self.get_chr_names()

        # the sorts, each chromosome on its own: the large ones on
        # threads (map_chromosomes), any other here and now
        big = []
        sizes = []
        for c in chrnames:
            self.locations[c].resize((self.size[c]), refcheck=False)
            if self.size[c] == 0:
                if c in self.size:
                    del self.size[c]
                if c in self.locations:
                    del self.locations[c]
                if c in self.rlengths:
                    del self.rlengths[c]
                continue
            if self.size[c] >= _MIN_THREADED_SIZE:
                big.append(c)
                sizes.append(self.size[c])
            else:
                self._sort_chrom(c)
            self.total += self.size[c]
        map_chromosomes(self._sort_chrom, big, sizes)

        self.is_sorted = True
        self.average_template_length = cython.cast(cython.float, self.length) / self.total
        return

    @cython.ccall
    def _sort_chrom(self, c: bytes):
        """finalize's sort of chromosome `c`, which touches nothing
        shared, so chromosomes can run on separate threads."""
        sort_lr_array(self.locations[c])

    @cython.ccall
    def get_locations_by_chr(self, chromosome: bytes):
        """Return the fragment array for a chromosome.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name, provided as bytes.
        
        Returns
        -------
        numpy.ndarray
            Structured array with ``('l', 'i4')`` and ``('r', 'i4')`` fields.
        
        Raises
        ------
        Exception
            If the chromosome is not present in the track.
        """
        if chromosome in self.locations:
            return self.locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return the set of chromosome names stored in the track.
        
        Returns
        -------
        set
            Chromosome names (bytes) that currently have fragments.
        """
        return set(self.locations.keys())

    @cython.ccall
    def sort(self):
        """Sort fragments for each chromosome by genomic coordinate.
        
        Fragments are ordered first by their left coordinate and then by their right
        coordinate. The ``is_sorted`` flag is set to ``True`` when sorting completes.
        """
        c: bytes
        chrnames: set

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c].sort(order=['l', 'r'])  # sort by the leftmost location
        self.is_sorted = True
        return

    @cython.ccall
    def count_fraglengths(self) -> dict:
        """Count observed fragment lengths across the track.
        
        Returns
        -------
        dict
            Mapping from fragment length to observed count, useful for downstream
            models such as HMMRATAC.
        """
        sizes: cnp.ndarray(cnp.int32_t, ndim=1)
        s: cython.int
        locs: cnp.ndarray
        chrnames: list
        i: cython.int

        counter = Counter()
        chrnames = list(self.get_chr_names())
        for i in range(len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = locs['r'] - locs['l']
            for s in sizes:
                counter[s] += 1
        return dict(counter)

    @cython.ccall
    def fraglengths(self) -> cnp.ndarray:
        """Return all fragment lengths as a single array.
        
        Returns
        -------
        numpy.ndarray
            Concatenated array of ``end - start`` for every stored fragment across
            chromosomes.
        """
        sizes: cnp.ndarray(np.int32_t, ndim=1)
        locs: cnp.ndarray
        chrnames: list
        i: cython.int

        chrnames = list(self.get_chr_names())
        locs = self.locations[chrnames[0]]
        sizes = locs['r'] - locs['l']
        for i in range(1, len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = np.concatenate((sizes, locs['r'] - locs['l']))
        return sizes

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def exclude(self, regions):
        """Remove fragments that overlap the provided exclusion regions.
        
        Parameters
        ----------
        regions : MACS3.Signal.Region.Regions
            Sorted region collection whose intervals should be excluded.
            The default path merges overlapping or adjacent intervals on a copy.
        
        Notes
        -----
        The operation mutates the track in place and finishes by calling
        :meth:`finalize` to refresh cached statistics.
        """
        i: cython.ulong
        j: cython.ulong
        k: bytes
        locs: cnp.ndarray
        locs_size: cython.ulong
        chrnames: set
        regions_c: list
        selected_idx: cnp.ndarray
        regions_chrs: list
        r1: cnp.void            # this is the location in numpy.void -- like a tuple
        r2: tuple               # this is the region
        n_rl1: cython.long
        n_rl2: cython.long

        if not self.is_sorted:
            self.sort()

        assert isinstance(regions, Regions)
        regions.sort()
        regions_chrs = list(regions.regions.keys())

        chrnames = self.get_chr_names()

        for k in chrnames:      # for each chromosome
            locs = self.locations[k]
            locs_size = self.size[k]
            # let's check if k is in regions_chr
            if k not in regions_chrs:
                # do nothing and continue
                self.total += locs_size
                continue

            # discard overlapping reads and make a new locations[k]
            # initialize boolean array as all TRUE, or all being kept
            selected_idx = np.ones(locs_size, dtype=bool)

            regions_c = regions.regions[k]

            i = 0
            n_rl1 = locs_size   # the number of locations left
            n_rl2 = len(regions_c)  # the number of regions left
            rl1_k = iter(locs).__next__
            rl2_k = iter(regions_c).__next__
            r1 = rl1_k()        # take the first value
            n_rl1 -= 1          # remaining rl1
            r2 = rl2_k()
            n_rl2 -= 1          # remaining rl2
            while (True):
                # we do this until there is no r1 or r2 left.
                if r2[0] < r1[1] and r1[0] < r2[1]:
                    # since we found an overlap, r1 will be skipped/excluded
                    # and move to the next r1
                    # get rid of this one
                    n_rl1 -= 1
                    self.length -= cython.cast(cython.ulonglong, r1[1] - r1[0])
                    selected_idx[i] = False

                    if n_rl1 > 0:
                        r1 = rl1_k()  # take the next location
                        i += 1
                        continue
                    else:
                        break
                if r1[1] < r2[1]:
                    # in this case, we need to move to the next r1,
                    n_rl1 -= 1
                    if n_rl1 > 0:
                        r1 = rl1_k()  # take the next location
                        i += 1
                    else:
                        # no more r1 left
                        break
                else:
                    # in this case, we need to move the next r2
                    if n_rl2:
                        r2 = rl2_k()  # take the next region
                        n_rl2 -= 1
                    else:
                        # no more r2 left
                        break

            self.locations[k] = locs[selected_idx]
            self.size[k] = self.locations[k].shape[0]
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            selected_idx.resize(self.buffer_size, refcheck=False)
            selected_idx.resize(0, refcheck=False)
        self.finalize()         # use length to set average_template_length and total
        return

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def filter_dup(self, maxnum: cython.int = -1):
        """Limit the number of duplicate fragments at identical coordinates.
        
        Parameters
        ----------
        maxnum : int, optional
            Maximum number of fragments allowed per unique ``(start, end)`` pair.
            A negative value disables duplicate filtering.
        
        Notes
        -----
        Fragments exceeding ``maxnum`` for the same coordinates are dropped and the
        aggregate template length is adjusted accordingly.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.add_loc(b"chr1", 10, 40)
            pe.finalize()
            pe.filter_dup(maxnum=1)
        """
        k: bytes
        i: cython.Py_ssize_t
        big: list
        sizes: list
        removed: list

        if maxnum < 0:
            return              # condition to return if not filtering

        if not self.is_sorted:
            self.sort()

        self.total = 0
        # self.length = 0
        self.average_template_length = 0.0

        # each chromosome on its own: the large ones on threads
        # (map_chromosomes), any other here and now. length is
        # unsigned, so subtracting each chromosome's sum wraps as
        # subtracting each fragment's length does
        big = []
        sizes = []
        for k in self.get_chr_names():
            if self.size[k] >= _MIN_THREADED_SIZE:
                big.append(k)
                sizes.append(self.size[k])
            else:
                self.length -= self._filter_dup_chrom(k, maxnum)
                self.total += self.size[k]
        removed = map_chromosomes(self._filter_dup_chrom, big, sizes, maxnum)
        for i in range(len(big)):
            self.length -= cython.cast(cython.ulonglong, removed[i])
            self.total += self.size[big[i]]
        self.average_template_length = self.length / self.total
        return

    @cython.ccall
    def _filter_dup_chrom(self, k: bytes,
                          maxnum: cython.int) -> cython.ulonglong:
        """filter_dup's work on chromosome `k`: drop the duplicates
        past `maxnum` and set locations[k] and size[k], and return the
        sum of the dropped fragments' lengths, each taken as uint64, for
        the caller to subtract from length. Touches nothing shared but
        the existing key k of locations and size, so chromosomes can run
        on separate threads.

        A plain one-dimensional array that owns its data and is held by
        nothing but locations[k] is filtered in place (filter_dup_lr)
        and shrunk, and stays locations[k]; any other takes the
        original loop and a new array, leaving the old one as it was."""
        n: cython.int
        loc_start: cython.int
        loc_end: cython.int
        current_loc_start: cython.int
        current_loc_end: cython.int
        i: cython.ulong
        locs_size: cython.ulong
        locs: cnp.ndarray
        selected_idx: cnp.ndarray
        p: cython.pointer(cython.int)
        kept: cython.Py_ssize_t = 0
        removed: cython.ulonglong = 0

        locs = self.locations[k]
        locs_size = self.size[k]
        if locs_size == 1:
            # do nothing
            return removed
        # the references to locs are locations[k]'s and this
        # function's, so no view or other holder of the array sees it
        # change
        if (locs_size >= 2 and locs.shape[0] == locs_size and
                Py_REFCNT(locs) == 2 and
                locs.ndim == 1 and is_plain_lr_array(locs) and
                cnp.PyArray_CHKFLAGS(locs, cnp.NPY_ARRAY_OWNDATA) and
                cnp.PyArray_BASE(locs) == cython.NULL):
            # the loop below, in C, in place and without the GIL
            p = cython.cast(cython.pointer(cython.int), locs.data)
            with cython.nogil:
                removed = filter_dup_lr(p, locs_size, maxnum,
                                        cython.address(kept))
            locs.resize(kept, refcheck=False)
            self.size[k] = locs.shape[0]
            return removed
        # discard duplicate reads and make a new locations[k]
        # initialize boolean array as all TRUE, or all being kept
        selected_idx = np.ones(locs_size, dtype=bool)
        # get the first loc
        (current_loc_start, current_loc_end) = locs[0]
        i = 1               # index of new_locs
        n = 1  # the number of tags in the current genomic location
        for i in range(1, locs_size):
            (loc_start, loc_end) = locs[i]
            if loc_start != current_loc_start or loc_end != current_loc_end:
                # not the same, update currnet_loc_start/end/l, reset n
                current_loc_start = loc_start
                current_loc_end = loc_end
                n = 1
                continue
            else:
                # both ends are the same, add 1 to duplicate number n
                n += 1
                if n > maxnum:
                    # change the flag to False
                    selected_idx[i] = False
                    # for the caller to subtract from self.length
                    removed += cython.cast(cython.ulonglong, current_loc_end - current_loc_start)
        self.locations[k] = locs[selected_idx]
        self.size[k] = self.locations[k].shape[0]
        # free memory?
        # I know I should shrink it to 0 size directly,
        # however, on Mac OSX, it seems directly assigning 0
        # doesn't do a thing.
        selected_idx.resize(self.buffer_size, refcheck=False)
        selected_idx.resize(0, refcheck=False)
        return removed

    @cython.ccall
    def sample_percent(self,
                       percent: cython.float,
                       seed: cython.int = -1):
        """Down-sample fragments in place by a fixed percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of fragments to keep per chromosome between 0 and 1 (inclusive).
        seed : int, optional
            Deterministic seed for the RNG; a negative value uses NumPy's global state.
        
        Notes
        -----
        Sampling is performed independently for each chromosome by shuffling the
        fragments, resizing the arrays, and restoring coordinate order.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.add_loc(b"chr1", 50, 90)
            pe.finalize()
            pe.sample_percent(0.5, seed=123)
        """
        # num: number of reads allowed on a certain chromosome
        num: cython.uint
        k: bytes
        chrnames: set

        self.total = 0
        self.length = 0
        self.average_template_length = 0.0

        chrnames = self.get_chr_names()

        if seed >= 0:
            info(f"#   A random seed {seed} has been used")
            rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
            rs_shuffle = rs.shuffle
        else:
            rs_shuffle = np.random.shuffle

        for k in sorted(chrnames):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            num = cython.cast(cython.uint,
                              round(self.locations[k].shape[0] * percent, 5))
            rs_shuffle(self.locations[k])
            self.locations[k].resize(num, refcheck=False)
            self.locations[k].sort(order=['l', 'r'])  # sort by leftmost positions
            self.size[k] = self.locations[k].shape[0]
            self.length += (self.locations[k]['r'] - self.locations[k]['l']).sum()
            self.total += self.size[k]
        self.average_template_length = cython.cast(cython.float, self.length)/self.total
        return

    @cython.ccall
    def sample_percent_copy(self,
                            percent: cython.float,
                            seed: cython.int = -1):
        """Return a down-sampled copy of the track.
        
        Parameters
        ----------
        percent : float
            Fraction of fragments to retain per chromosome between 0 and 1.
        seed : int, optional
            Deterministic seed used when shuffling; a negative value disables seeding.
        
        Returns
        -------
        PETrackI
            New track containing the sampled fragments with metadata copied over.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.add_loc(b"chr1", 50, 90)
            pe.finalize()
            pe_subset = pe.sample_percent_copy(0.5, seed=123)
        """
        # num: number of reads allowed on a certain chromosome
        num: cython.uint
        k: bytes
        chrnames: set
        ret_petrackI: PETrackI
        loc: cnp.ndarray

        ret_petrackI = PETrackI(anno=self.annotation,
                                buffer_size=self.buffer_size)
        chrnames = self.get_chr_names()

        if seed >= 0:
            info(f"# A random seed {seed} has been used in the sampling function")
            rs = np.random.default_rng(seed)
        else:
            rs = np.random.default_rng()

        rs_shuffle = rs.shuffle

        # chrnames need to be sorted otherwise we can't assure reproducibility
        for k in sorted(chrnames):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            loc = np.copy(self.locations[k])
            num = cython.cast(cython.uint, round(loc.shape[0] * percent, 5))
            rs_shuffle(loc)
            loc.resize(num, refcheck=False)
            loc.sort(order=['l', 'r'])  # sort by leftmost positions
            ret_petrackI.locations[k] = loc
            ret_petrackI.size[k] = loc.shape[0]
            ret_petrackI.length += (loc['r'] - loc['l']).sum()
            ret_petrackI.total += ret_petrackI.size[k]
        ret_petrackI.average_template_length = cython.cast(cython.float, ret_petrackI.length)/ret_petrackI.total
        ret_petrackI.set_rlengths(self.get_rlengths())
        return ret_petrackI

    @cython.ccall
    def sample_num(self,
                   samplesize: cython.ulong,
                   seed: cython.int = -1):
        """Down-sample fragments in place to approximately ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target number of fragments across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent`.
        
        Notes
        -----
        The method converts ``samplesize`` into a sampling fraction using ``self.total``.
        Ensure :meth:`finalize` has been called so counts are up to date.
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize)/self.total
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def sample_num_copy(self,
                        samplesize: cython.ulong,
                        seed: cython.int = -1):
        """Return a down-sampled copy with approximately ``samplesize`` fragments.
        
        Parameters
        ----------
        samplesize : int
            Target number of fragments across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent_copy`.
        
        Returns
        -------
        PETrackI
            New track containing the sampled fragments.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.add_loc(b"chr1", 50, 90)
            pe.finalize()
            pe_subset = pe.sample_num_copy(samplesize=1, seed=123)
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize)/self.total
        return self.sample_percent_copy(percent, seed)

    @cython.ccall
    def print_to_bed(self, fhd=None):
        """Write fragments to a three-column BEDPE-style stream.
        
        Parameters
        ----------
        fhd : io.IOBase, optional
            Writable file-like object. Defaults to ``sys.stdout``.
        
        Notes
        -----
        Each fragment is emitted as ``chrom	start	end`` using decoded chromosome
        names and the stored integer coordinates.

        Examples
        --------
        .. code-block:: python

            import io
            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.finalize()
            buf = io.StringIO()
            pe.print_to_bed(buf)
        """
        i: cython.int
        s: cython.int
        e: cython.int
        k: bytes
        chrnames: set

        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, io.IOBase)

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            locs = self.locations[k]

            for i in range(locs.shape[0]):
                s, e = locs[i]
                fhd.write("%s\t%d\t%d\n" % (k.decode(), s, e))
        return

    @cython.ccall
    def pileup_a_chromosome(self,
                            chrom: bytes,
                            scale_factor: cython.float = 1.0,
                            baseline_value: cython.float = 0.0) -> list:
        """Compute a coverage pileup for a single chromosome.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        scale_factor : float, optional
            Value used to scale the resulting coverage.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` with numpy arrays describing
            the pileup breakpoints and scaled coverage.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.add_loc(b"chr1", 50, 90)
            pe.finalize()
            p, v = pe.pileup_a_chromosome(b"chr1")
        """
        tmp_pileup: list

        tmp_pileup = pileup_from_LR_as_list(self.locations[chrom],
                                            scale_factor,
                                            baseline_value)
        return tmp_pileup

    @cython.ccall
    def pileup_a_chromosome_c(self,
                              chrom: bytes,
                              ds: list,
                              scale_factor_s: list,
                              baseline_value: cython.float = 0.0) -> list:
        """Project paired-end fragments into pseudo single-end pileups.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        ds : list[int]
            Fragment lengths used to build the projections.
        scale_factor_s : list[float]
            Scale factors paired with each entry in ``ds``.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` representing the merged pileup
            with the maximum value taken across projections.
        """
        tmp_pileup: list
        prev_pileup: list
        five_shift_s: list
        scale_factor: cython.float
        d: cython.long
        five_shift: cython.long
        three_shift: cython.long
        rlength: cython.long = self.get_rlengths()[chrom]

        if not self.is_sorted:
            self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        # three windows (d, slocal and llocal): their pileups and the
        # maximum in one sweep, the same arrays as the loop below
        if len(ds) == 3:
            five_shift_s = []
            for i in range(3):
                d = ds[i]
                five_shift_s.append(d//2)
            prev_pileup = se_all_in_one_pileup_max3(self.locations[chrom]['l'],
                                                    self.locations[chrom]['r'],
                                                    five_shift_s,
                                                    five_shift_s,
                                                    rlength,
                                                    scale_factor_s,
                                                    baseline_value)
            if prev_pileup is not None:
                return prev_pileup

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            five_shift = d//2
            three_shift = d//2

            tmp_pileup = pileup_from_PN_shifted(self.locations[chrom]['l'],
                                                self.locations[chrom]['r'],
                                                five_shift,
                                                three_shift,
                                                rlength,
                                                scale_factor,
                                                baseline_value)

            if prev_pileup:
                prev_pileup = over_two_pv_array(prev_pileup,
                                                tmp_pileup,
                                                func="max")
            else:
                prev_pileup = tmp_pileup

        return prev_pileup

    @cython.ccall
    def pileup_bdg(self,
                   scale_factor: cython.float = 1.0,
                   baseline_value: cython.float = 0.0):
        """Build a ``bedGraphTrackI`` with pileups for every chromosome.
        
        Parameters
        ----------
        scale_factor : float, optional
            Value used to scale the coverage for each chromosome.
        baseline_value : float, optional
            Minimum value enforced on the coverage arrays.
        
        Returns
        -------
        bedGraphTrackI
            BedGraph track populated with per-chromosome pileup data.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.finalize()
            bdg = pe.pileup_bdg()
        """
        tmp_pileup: list
        chrom: bytes
        bdg: bedGraphTrackI

        bdg = bedGraphTrackI(baseline_value=baseline_value)

        for chrom in sorted(self.get_chr_names()):
            tmp_pileup = pileup_from_LR_as_list(self.locations[chrom],
                                                scale_factor,
                                                baseline_value)

            # save to bedGraph
            bdg.add_chrom_data(chrom,
                               pyarray('i', tmp_pileup[0]),
                               pyarray('f', tmp_pileup[1]))
        return bdg

    @cython.ccall
    def pileup_bdg_hmmr(self,
                        mapping: list,
                        baseline_value: cython.float = 0.0) -> list:
        """Generate HMMRATAC-style pileups for every chromosome.
        
        Parameters
        ----------
        mapping : list
            Weight mapping produced by HMMRATAC EM training describing the short,
            mono-, di-, and tri-nucleosomal signals.
        baseline_value : float, optional
            Reserved parameter for API compatibility; not currently applied.
        
        Returns
        -------
        list
            List of dictionaries mirroring ``mapping`` where each dictionary maps
            chromosome names to pileup arrays returned by
            :func:`pileup_from_LR_hmmratac`.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackI

            pe = PETrackI()
            pe.add_loc(b"chr1", 10, 40)
            pe.finalize()
            mapping = [{"short": 1.0, "mono": 0.0, "di": 0.0, "tri": 0.0}]
            pileups = pe.pileup_bdg_hmmr(mapping)
        """
        ret_pileup: list
        chroms: set
        chrom: bytes
        i: cython.int

        ret_pileup = []
        for i in range(len(mapping)):
            ret_pileup.append({})
        chroms = self.get_chr_names()
        for i in range(len(mapping)):
            for chrom in sorted(chroms):
                ret_pileup[i][chrom] = pileup_from_LR_hmmratac(self.locations[chrom], mapping[i])
        return ret_pileup


@cython.cclass
class PETrackII:
    """Paired-end track for single-cell ATAC fragments with barcode metadata.
    
    Each chromosome stores a structured array of fragment coordinates and counts
    alongside an integer-encoded barcode array to support barcode-aware analyses.

    Examples
    --------
    .. code-block:: python

        from MACS3.Signal.PairedEndTrack import PETrackII

        pe = PETrackII(anno="scatac")
        pe.add_loc(b"chr1", 10, 40, barcode=b"BC01", count=1)
        pe.add_loc(b"chr1", 50, 90, barcode=b"BC02", count=2)
        pe.finalize()
        pe.set_rlengths({b"chr1": 1000})
    """
    locations = cython.declare(dict, visibility="public")
    # add another dict for storing barcode for each fragment we will
    # first convert barcode into integer and remember them in the
    # barcode_dict, which will map the rule to numerize
    # key:bytes as value:4bytes_integer
    barcodes = cython.declare(dict, visibility="public")
    barcode_dict = cython.declare(dict, visibility="public")
    # the last number for barcodes, used to map barcode to integer
    barcode_last_n: cython.int

    size = cython.declare(dict, visibility="public")
    buf_size = cython.declare(dict, visibility="public")
    is_sorted = cython.declare(bool, visibility="public")
    total = cython.declare(cython.ulong, visibility="public")
    annotation = cython.declare(str, visibility="public")
    # rlengths: reference chromosome lengths dictionary
    rlengths = cython.declare(dict, visibility="public")
    buffer_size = cython.declare(cython.long, visibility="public")
    length = cython.declare(cython.ulonglong, visibility="public")  # total length of all fragments
    average_template_length = cython.declare(cython.float, visibility="public")
    is_destroyed: bool

    def __init__(self, anno: str = "", buffer_size: cython.long = 100000):
        # dictionary with chrname as key, nparray with
        # [('l','i4'),('r','i4'),('c','u2')] as value
        self.locations = {}
        # dictionary with chrname as key, size of the above nparray as value
        # size is to remember the size of the fragments added to this chromosome
        self.size = {}
        # dictionary with chrname as key, size of the above nparray as value
        self.buf_size = {}
        self.is_sorted = False
        self.total = 0           # total fragments
        self.annotation = anno   # need to be figured out
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.average_template_length = 0.0
        self.is_destroyed = False

        self.barcodes = {}
        self.barcode_dict = {}
        self.barcode_last_n = 0

    @cython.ccall
    def add_loc(self,
                chromosome: bytes,
                start: cython.int,
                end: cython.int,
                barcode: bytes,
                count: cython.ushort):
        """Append a fragment together with its barcode and count.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) for the fragment.
        start : int
            Zero-based start coordinate of the fragment.
        end : int
            Zero-based end coordinate of the fragment.
        barcode : bytes
            Raw barcode sequence associated with the fragment.
        count : int
            Number of occurrences represented by the fragment.
        
        Notes
        -----
        Barcodes are interned into integers via ``barcode_dict`` for compact storage
        and the accumulated template length is weighted by ``count``.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackII

            pe = PETrackII()
            pe.add_loc(b"chr1", 10, 40, barcode=b"BC01", count=1)
            pe.add_loc(b"chr1", 50, 90, barcode=b"BC02", count=2)
        """
        i: cython.int
        # bn: the integer in barcode_dict for this barcode
        bn: cython.int

        if barcode not in self.barcode_dict:
            self.barcode_dict[barcode] = self.barcode_last_n
            self.barcode_last_n += 1
        bn = self.barcode_dict[barcode]

        if chromosome not in self.locations:
            self.buf_size[chromosome] = self.buffer_size
            # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            # ['c'] is the count number of this fragment
            self.locations[chromosome] = np.zeros(shape=self.buffer_size,
                                                  dtype=[('l', 'i4'), ('r', 'i4'), ('c', 'u2')])
            self.barcodes[chromosome] = np.zeros(shape=self.buffer_size,
                                                 dtype='i4')
            self.locations[chromosome][0] = (start, end, count)
            self.barcodes[chromosome][0] = bn
            self.size[chromosome] = 1
        else:
            i = self.size[chromosome]
            if self.buf_size[chromosome] == i:
                self.buf_size[chromosome] += self.buffer_size
                self.locations[chromosome].resize((self.buf_size[chromosome]),
                                                  refcheck=False)
                self.barcodes[chromosome].resize((self.buf_size[chromosome]),
                                                 refcheck=False)                
            self.locations[chromosome][i] = (start, end, count)
            self.barcodes[chromosome][i] = bn
            self.size[chromosome] = i + 1
        self.length += (end - start) * count
        return

    @cython.ccall
    def barcode_id(self, barcode: bytes) -> cython.int:
        """Return the integer that encodes ``barcode`` in ``barcodes``.

        A barcode not seen before is given the next integer, as
        ``add_loc`` gives it, so ids follow the order of first appearance.
        """
        if barcode not in self.barcode_dict:
            self.barcode_dict[barcode] = self.barcode_last_n
            self.barcode_last_n += 1
        return self.barcode_dict[barcode]

    @cython.ccall
    def add_loc_arrays(self, chromosome: bytes, starts, ends, counts,
                       barcode_ids):
        """Append many fragments of one chromosome to the track.

        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) for the fragments.
        starts : numpy.ndarray
            int32 leftmost ends of the fragments, in the order to append.
        ends : numpy.ndarray
            int32 rightmost ends, same length as ``starts``.
        counts : numpy.ndarray
            uint16 counts, same length.
        barcode_ids : numpy.ndarray
            int32 barcode ids from ``barcode_id``, same length.

        Notes
        -----
        Leaves the track as calling ``add_loc(chromosome, starts[k],
        ends[k], barcode, counts[k])`` for each ``k`` in turn would,
        where ``barcode`` is the barcode whose id is
        ``barcode_ids[k]``: the same records, sizes, buffer sizes and
        ``length``. The arrays may be longer than ``buf_size``, with
        zeros past ``size`` (see ``reserve``); after
        ``trim_to_buf_size`` they are exactly add_loc's arrays. Needs a
        positive ``buffer_size``.
        """
        i: cython.long
        n: cython.long
        b: cython.long
        dlength: cython.longlong

        n = len(starts)
        if n == 0:
            return
        if self.buffer_size <= 0:
            raise ValueError("add_loc_arrays needs a positive buffer_size")

        if chromosome not in self.locations:
            self.buf_size[chromosome] = self.buffer_size
            self.locations[chromosome] = np.zeros(shape=self.buffer_size,
                                                  dtype=[('l', 'i4'), ('r', 'i4'), ('c', 'u2')])
            self.barcodes[chromosome] = np.zeros(shape=self.buffer_size,
                                                 dtype='i4')
            self.size[chromosome] = 0
        i = self.size[chromosome]
        b = self.buf_size[chromosome]
        if i + n > b:
            # buf_size grows in steps of buffer_size, as add_loc's does
            while b < i + n:
                b += self.buffer_size
            self.buf_size[chromosome] = b
            if (self.locations[chromosome].shape[0] < b or
                    self.barcodes[chromosome].shape[0] < b):
                self.reserve(chromosome, b)
        locs = self.locations[chromosome]
        bcs = self.barcodes[chromosome]
        if (_lrc_append(locs, bcs, i, n, starts, ends, counts, barcode_ids,
                        cython.address(dlength)) == 0):
            locs['l'][i:i + n] = starts
            locs['r'][i:i + n] = ends
            locs['c'][i:i + n] = counts
            bcs[i:i + n] = barcode_ids
            self.size[chromosome] = i + n
            # add_loc adds each C int (end - start) * count, wrapping as C
            # int arithmetic does, to length
            dlength = np.multiply(np.subtract(ends, starts, dtype=np.int32),
                                  counts, dtype=np.int32).sum(dtype=np.int64)
        else:
            self.size[chromosome] = i + n
        self.length += dlength
        return

    @cython.cfunc
    def reserve(self, chromosome: bytes, b: cython.long):
        """Give the arrays of ``chromosome`` room for at least ``b``
        records and four times their current length: new zeroed
        arrays (numpy asks for huge pages at 4 MB or more) into which
        the first ``size`` records are copied.

        Growing in place with ``resize`` in steps of ``buffer_size``,
        as add_loc does, zero-fills each step and so faults in every
        page of it on the spot, one 4 kB page at a time; that cost
        more than these copies. Memory past what is written is never
        touched, so the larger capacity costs address space only.
        """
        i: cython.long
        cap: cython.long

        locs = self.locations[chromosome]
        bcs = self.barcodes[chromosome]
        i = self.size[chromosome]
        cap = max(b, 4 * cython.cast(cython.long, locs.shape[0]))
        new_locs = np.zeros(cap, dtype=locs.dtype)
        new_bcs = np.zeros(cap, dtype=bcs.dtype)
        new_locs[:i] = locs[:i]
        new_bcs[:i] = bcs[:i]
        self.locations[chromosome] = new_locs
        self.barcodes[chromosome] = new_bcs

    @cython.ccall
    def trim_to_buf_size(self):
        """Shorten each chromosome's arrays that ``reserve`` made longer
        than ``buf_size`` to ``buf_size`` records, the length add_loc
        gives them. Called by ``set_rlengths``, which FragParser calls
        once a file is read."""
        c: bytes
        b: cython.long

        for c in self.locations:
            b = self.buf_size[c]
            if self.locations[c] is not None and self.locations[c].shape[0] > b:
                self.locations[c].resize((b), refcheck=False)
            if self.barcodes[c] is not None and self.barcodes[c].shape[0] > b:
                self.barcodes[c].resize((b), refcheck=False)

    @cython.ccall
    def destroy(self):
        """Release fragment and barcode arrays held by the track.
        
        All per-chromosome arrays are resized to zero, barcode mappings are cleared,
        and the track is marked as destroyed.
        """
        chrs: set
        chromosome: bytes

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.locations:
                self.locations[chromosome].resize(self.buffer_size,
                                                  refcheck=False)
                self.locations[chromosome].resize(0,
                                                  refcheck=False)
                self.locations[chromosome] = None
                self.locations.pop(chromosome)
                self.barcodes.resize(self.buffer_size,
                                     refcheck=False)
                self.barcodes.resize(0,
                                     refcheck=False)
                self.barcodes[chromosome] = None
                self.barcodes.pop(chromosome)
        self.barcode_dict = {}
        self.is_destroyed = True
        return

    @cython.ccall
    def set_rlengths(self, rlengths: dict) -> bool:
        """Attach reference chromosome lengths to the track.
        
        Parameters
        ----------
        rlengths : dict
            Mapping from chromosome name (bytes) to reference length.
        
        Returns
        -------
        bool
            True when the length mapping has been updated.
        
        Notes
        -----
        Any chromosome stored in the track but missing from ``rlengths`` is assigned
        ``INT_MAX`` so downstream bounds checks can succeed.

        FragParser calls this once a file is read, so it is also where
        arrays that ``add_loc_arrays`` made longer than ``buf_size`` are
        cut back to add_loc's length (``trim_to_buf_size``); on any other
        track that is a no-op.
        """
        valid_chroms: set
        missed_chroms: set
        chrom: bytes

        self.trim_to_buf_size()
        valid_chroms = set(self.locations.keys()).intersection(rlengths.keys())
        for chrom in sorted(valid_chroms):
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.locations.keys()).difference(rlengths.keys())
        for chrom in sorted(missed_chroms):
            self.rlengths[chrom] = INT_MAX
        return True

    @cython.ccall
    def get_rlengths(self) -> dict:
        """Return the reference chromosome lengths associated with the track.
        
        Returns
        -------
        dict
            Mapping from chromosome name (bytes) to reference length. Chromosomes
            without a recorded length default to ``INT_MAX``.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

    @cython.ccall
    def finalize(self):
        """Shrink arrays, sort fragments, and refresh aggregate counters.
        
        Each per-chromosome fragment array is resized to its observed length, sorted
        by ``('l', 'r')``, and the accompanying barcode array is reordered to match.
        The method updates ``total`` and ``average_template_length`` using count
        weights and marks the track as sorted.
        
        Raises
        ------
        AssertionError
            If no fragments are present when finalizing.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackII

            pe = PETrackII()
            pe.add_loc(b"chr1", 10, 40, barcode=b"BC01", count=1)
            pe.finalize()
        """
        c: bytes
        chrnames: set
        big: list
        sizes: list

        self.total = 0

        chrnames = self.get_chr_names()

        # the sorts, each chromosome on its own: the large ones on
        # threads (map_chromosomes), any other here and now
        big = []
        sizes = []
        for c in chrnames:
            self.locations[c].resize((self.size[c]), refcheck=False)
            if self.size[c] == 0:
                if c in self.size:
                    del self.size[c]
                if c in self.locations:
                    del self.locations[c]
                if c in self.barcodes:
                    del self.barcodes[c]
                if c in self.rlengths:
                    del self.rlengths[c]
                continue
            if self.size[c] >= _MIN_THREADED_SIZE:
                big.append(c)
                sizes.append(self.size[c])
            else:
                self.total += self._sort_chrom(c)  # self.size[c]
        for count in map_chromosomes(self._sort_chrom, big, sizes):
            self.total += count  # self.size[c]

        assert self.total > 0, "Error: no fragments in PETrackII"

        self.is_sorted = True
        self.average_template_length = cython.cast(cython.float,
                                                   self.length) / self.total
        return

    @cython.ccall
    def _sort_chrom(self, c: bytes):
        """finalize's sort of chromosome `c`, fragments and barcodes,
        returning the sum of its counts for the caller to add to total.
        Touches nothing shared but the existing key c of locations and
        barcodes, so chromosomes can run on separate threads."""
        indices: cnp.ndarray

        indices = argsort_lrc(self.locations[c])
        self.locations[c] = self.locations[c][indices]
        self.barcodes[c] = self.barcodes[c][indices]
        return np.sum(self.locations[c]['c'])

    @cython.ccall
    def get_locations_by_chr(self, chromosome: bytes):
        """Return the fragment array for a chromosome.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name, provided as bytes.
        
        Returns
        -------
        numpy.ndarray
            Structured array with ``('l', 'i4')``, ``('r', 'i4')``, and ``('c', 'u2')`` fields.
        
        Raises
        ------
        Exception
            If the chromosome is not present in the track.
        """
        if chromosome in self.locations:
            return self.locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return the set of chromosome names stored in the track.
        
        Returns
        -------
        set
            Chromosome names (bytes) that currently have fragments.
        """
        return set(self.locations.keys())

    @cython.ccall
    def sort(self):
        """Sort fragments and barcodes for each chromosome.
        
        Fragments are ordered first by their left coordinate and then by their right
        coordinate, and the barcode array is reordered alongside the fragment array.
        The ``is_sorted`` flag is set to ``True`` when sorting completes.
        """
        c: bytes
        chrnames: set
        indices: cnp.ndarray

        chrnames = self.get_chr_names()

        for c in chrnames:
            indices = argsort_lrc(self.locations[c])
            self.locations[c] = self.locations[c][indices]
            self.barcodes[c] = self.barcodes[c][indices]
        self.is_sorted = True
        return

    @cython.ccall
    def count_fraglengths(self) -> dict:
        """Count fragment lengths weighted by per-fragment counts.
        
        Returns
        -------
        dict
            Mapping from fragment length to the total count contributed by fragments
            of that length.
        """
        sizes: cnp.ndarray(cnp.int32_t, ndim=1)
        s: cython.int
        locs: cnp.ndarray
        chrnames: list
        i: cython.int
        j: cython.int

        counter = Counter()
        chrnames = list(self.get_chr_names())
        for i in range(len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = locs['r'] - locs['l']
            for j in range(len(sizes)):
                s = sizes[j]
                counter[s] += locs['c'][j]
        return dict(counter)

    @cython.ccall
    def fraglengths(self) -> cnp.ndarray:
        """Return all fragment lengths expanded by their counts.
        
        Returns
        -------
        numpy.ndarray
            Array of ``end - start`` values repeated according to the stored counts.
        """
        sizes: cnp.ndarray(np.int32_t, ndim=1)
        chrnames: list
        out: list
        chrom: bytes

        chrnames = list(self.get_chr_names())
        out = []
        for chrom in chrnames:
            locs = self.locations[chrom]
            sizes = locs['r'] - locs['l']
            counts = locs['c']
            out.append(np.repeat(sizes, counts))
        if out:
            return np.concatenate(out).astype(np.int32)
        else:
            return np.array([], dtype=np.int32)

    @cython.ccall
    def subset(self, selected_barcodes: set):
        """Build a new track containing only fragments from selected barcodes.
        
        Parameters
        ----------
        selected_barcodes : set
            Set of barcode byte strings to retain.
        
        Returns
        -------
        PETrackII
            New track restricted to the provided barcodes with metadata preserved.

        Examples
        --------
        .. code-block:: python

            from MACS3.Signal.PairedEndTrack import PETrackII

            pe = PETrackII()
            pe.add_loc(b"chr1", 10, 40, barcode=b"BC01", count=1)
            pe.add_loc(b"chr1", 50, 90, barcode=b"BC02", count=1)
            pe.finalize()
            subset = pe.subset({b"BC01"})
        """
        indices: cnp.ndarray
        chrs: set
        selected_barcodes_filtered: list
        selected_barcodes_n: list
        chromosome: bytes
        ret: PETrackII

        ret = PETrackII()
        chrs = self.get_chr_names()

        # first we need to convert barcodes into integers in our
        # barcode_dict
        selected_barcodes_filtered = [b
                                      for b in selected_barcodes
                                      if b in self.barcode_dict]
        ret.barcode_dict = {b: self.barcode_dict[b]
                            for b in selected_barcodes_filtered}
        selected_barcodes_n = [self.barcode_dict[b]
                               for b in selected_barcodes_filtered]
        ret.barcode_last_n = self.barcode_last_n

        # pass some values from self to ret
        ret.annotation = self.annotation
        ret.is_sorted = self.is_sorted
        ret.rlengths = self.rlengths
        ret.buffer_size = self.buffer_size
        # ret.total = 0
        ret.length = 0
        # ret.average_template_length = 0
        ret.is_destroyed = True

        for chromosome in sorted(chrs):
            indices = np.where(np.isin(self.barcodes[chromosome],
                                       list(selected_barcodes_n)))[0]
            ret.barcodes[chromosome] = self.barcodes[chromosome][indices]
            ret.locations[chromosome] = self.locations[chromosome][indices]
            ret.size[chromosome] = len(ret.locations[chromosome])
            ret.buf_size[chromosome] = ret.size[chromosome]
            # ret.total += np.sum(ret.locations[chromosome]['c'])
            ret.length += np.sum((ret.locations[chromosome]['r'] -
                                  ret.locations[chromosome]['l']) *
                                 ret.locations[chromosome]['c'])
        ret.finalize()
        # ret.average_template_length = ret.length / ret.total
        return ret

    @cython.ccall
    def pileup_a_chromosome(self,
                            chrom: bytes,
                            scale_factor: cython.float = 1.0,
                            baseline_value: cython.float = 0.0) -> list:
        """Compute a coverage pileup for a single chromosome.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        scale_factor : float, optional
            Value used to scale the resulting coverage.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` with numpy arrays describing
            the pileup breakpoints and scaled coverage.
        """
        ret: list

        # fragments with one count: the same arrays from int32 sorts
        ret = pileup_LRC_as_list_equal(self.locations[chrom],
                                       scale_factor,
                                       baseline_value,
                                       self.is_sorted)
        if ret is not None:
            return ret
        return pileup_from_LRC_as_list(self.locations[chrom],
                                       scale_factor,
                                       baseline_value,
                                       left_sorted=self.is_sorted)

    @cython.ccall
    def pileup_a_chromosome_c(self,
                              chrom: bytes,
                              ds: list,
                              scale_factor_s: list,
                              baseline_value: cython.float = 0.0) -> list:
        """Project paired-end fragments into pseudo single-end pileups.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        ds : list[int]
            Fragment lengths used to build the projections.
        scale_factor_s : list[float]
            Scale factors paired with each entry in ``ds``.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` representing the merged pileup
            with the maximum value taken across projections.
        """
        prev_pileup: list
        tmp_pileup: list
        scale_factor: cython.float
        d: cython.long

        ####
        if not self.is_sorted:
            self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            # fragments with one count: the same arrays from int32 sorts
            tmp_pileup = pileup_LRC_centers_as_list_equal(self.locations[chrom],
                                                          d,
                                                          scale_factor,
                                                          baseline_value)
            if tmp_pileup is None:
                tmp_pileup = pileup_from_LRC_centers_as_list(self.locations[chrom],
                                                             d,
                                                             scale_factor,
                                                             baseline_value)

            if prev_pileup:
                prev_pileup = over_two_pv_array(prev_pileup,
                                                tmp_pileup,
                                                func="max")
            else:
                prev_pileup = tmp_pileup

        return prev_pileup

    @cython.ccall
    def pileup_bdg(self,
                   scale_factor: cython.float = 1.0,
                   baseline_value: cython.float = 0.0):
        """Build a ``bedGraphTrackI`` with pileups for every chromosome.
        
        Parameters
        ----------
        scale_factor : float, optional
            Value used to scale the coverage for each chromosome.
        baseline_value : float, optional
            Minimum value enforced on the coverage arrays.
        
        Returns
        -------
        bedGraphTrackI
            BedGraph track populated with per-chromosome pileup data.
        """
        bdg: bedGraphTrackI
        tmp_pileup: list
        chrom: bytes

        bdg = bedGraphTrackI(baseline_value=baseline_value)
        for chrom in sorted(self.get_chr_names()):
            tmp_pileup = pileup_from_LRC_as_list(self.locations[chrom],
                                                 scale_factor,
                                                 baseline_value,
                                                 left_sorted=self.is_sorted)
            bdg.add_chrom_data(chrom,
                               pyarray('i', tmp_pileup[0]),
                               pyarray('f', tmp_pileup[1]))
        return bdg

    @cython.ccall
    def pileup_bdg2(self):
        """Build a ``bedGraphTrackII`` with pileups for every chromosome.
        
        Returns
        -------
        bedGraphTrackII
            BedGraph track populated with per-chromosome pileup arrays and finalized.
        """
        bdg: bedGraphTrackII
        pv: cnp.ndarray

        bdg = bedGraphTrackII()
        for chrom in self.get_chr_names():
            pv = pileup_from_LRC(self.locations[chrom])
            bdg.add_chrom_data(chrom, pv)
        # bedGraphTrackII needs to be 'finalized'.
        bdg.finalize()
        return bdg

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def exclude(self, regions):
        """Remove fragments that overlap the provided exclusion regions.
        
        Parameters
        ----------
        regions : MACS3.Signal.Region.Regions
            Sorted region collection whose intervals should be excluded.
        
        Notes
        -----
        The operation mutates the track in place, adjusts fragment counts and lengths,
        and finishes by calling :meth:`finalize`.
        """
        k: bytes
        locs: cnp.ndarray
        locs_size: cython.ulong
        chrnames: set
        merged_regions: Regions
        regions_c: list
        selected_idx: cnp.ndarray
        regions_chrs: set
        loc_starts: cnp.ndarray(cnp.int32_t, ndim=1)
        loc_ends: cnp.ndarray(cnp.int32_t, ndim=1)
        loc_counts: cnp.ndarray(cnp.uint16_t, ndim=1)
        region_starts: cnp.ndarray(cnp.int32_t, ndim=1)
        region_ends: cnp.ndarray(cnp.int32_t, ndim=1)
        region_idx: cnp.ndarray(cnp.int32_t, ndim=1)
        valid_mask: cnp.ndarray
        overlap_mask: cnp.ndarray
        removed_sizes: cnp.ndarray
        kept_counts: cnp.ndarray
        i: cython.int
        n_regions_c: cython.int
        total_counts: cython.ulonglong
        chrom_total: cython.ulonglong

        if not self.is_sorted:
            self.sort()

        assert isinstance(regions, Regions)
        merged_regions = Regions()
        for k in sorted(regions.regions.keys()):
            merged_regions.regions[k] = regions.regions[k][:]
            merged_regions.total += len(merged_regions.regions[k])
        merged_regions.merge_overlap()
        regions_chrs = set(merged_regions.regions.keys())

        chrnames = self.get_chr_names()
        total_counts = 0

        for k in chrnames:      # for each chromosome
            locs = self.locations[k]
            locs_size = self.size[k]
            # let's check if k is in regions_chr
            if k not in regions_chrs:
                total_counts += cython.cast(cython.ulonglong, np.sum(locs['c']))
                continue

            # discard overlapping reads and make a new locations[k]
            # initialize boolean array as all TRUE, or all being kept
            regions_c = merged_regions.regions[k]
            loc_starts = locs['l']
            loc_ends = locs['r']
            loc_counts = locs['c']
            n_regions_c = len(regions_c)
            region_starts = np.empty(n_regions_c, dtype=np.int32)
            region_ends = np.empty(n_regions_c, dtype=np.int32)
            for i in range(n_regions_c):
                region_starts[i] = regions_c[i][0]
                region_ends[i] = regions_c[i][1]

            region_idx = np.searchsorted(region_ends, loc_starts, side='right').astype(np.int32, copy=False)
            valid_mask = region_idx < n_regions_c
            overlap_mask = np.zeros(locs_size, dtype=bool)
            overlap_mask[valid_mask] = region_starts[region_idx[valid_mask]] < loc_ends[valid_mask]
            selected_idx = np.logical_not(overlap_mask)

            if np.any(overlap_mask):
                removed_sizes = (loc_ends[overlap_mask] - loc_starts[overlap_mask]).astype(np.uint64, copy=False)
                self.length -= cython.cast(cython.ulonglong,
                                           np.sum(removed_sizes * loc_counts[overlap_mask].astype(np.uint64, copy=False)))

            if np.all(selected_idx):
                total_counts += cython.cast(cython.ulonglong, np.sum(loc_counts))
                continue

            self.locations[k] = locs[selected_idx]
            self.barcodes[k] = self.barcodes[k][selected_idx]
            self.size[k] = self.locations[k].shape[0]

            if self.size[k] == 0:
                del self.size[k]
                del self.locations[k]
                del self.barcodes[k]
                if k in self.rlengths:
                    del self.rlengths[k]
                continue

            kept_counts = self.locations[k]['c']
            chrom_total = cython.cast(cython.ulonglong, np.sum(kept_counts))
            total_counts += chrom_total

        self.total = total_counts
        assert self.total > 0, "Error: no fragments in PETrackII"
        self.is_sorted = True
        self.average_template_length = cython.cast(cython.float,
                                                   self.length) / self.total
        return
    
    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.cfunc
    def _two_pointer_sweep(self, regions):
        peak_idx: cython.int
        n_regions_c: cython.int
        n_cells: cython.int
        peak_counter: cython.int
        peak_base: cython.int
        start: cython.int
        end: cython.int
        chrom: bytes
        chrom_str: str
        barcode_items: list
        regions_c: list
        barcode_ids: cnp.ndarray(cnp.int32_t, ndim=1)
        barcode_id_to_row: cnp.ndarray(cnp.int32_t, ndim=1)
        fragment_locs: cnp.ndarray
        fragment_barcodes: cnp.ndarray(cnp.int32_t, ndim=1)
        peak_starts: cnp.ndarray(cnp.int32_t, ndim=1)
        peak_ends: cnp.ndarray(cnp.int32_t, ndim=1)
        frag_starts: cnp.ndarray(cnp.int32_t, ndim=1)
        frag_ends: cnp.ndarray(cnp.int32_t, ndim=1)
        frag_counts: cnp.ndarray(cnp.uint16_t, ndim=1)
        frag_rows: cnp.ndarray(cnp.int32_t, ndim=1)
        left_idx: cnp.ndarray(cnp.int32_t, ndim=1)
        right_idx: cnp.ndarray(cnp.int32_t, ndim=1)
        widths: cnp.ndarray(cnp.int32_t, ndim=1)
        valid_mask: cnp.ndarray
        chunk_rows: cnp.ndarray(cnp.int32_t, ndim=1)
        chunk_cols: cnp.ndarray(cnp.int32_t, ndim=1)
        chunk_data: cnp.ndarray(cnp.int32_t, ndim=1)
        valid_rows: cnp.ndarray(cnp.int32_t, ndim=1)
        valid_counts: cnp.ndarray(cnp.int32_t, ndim=1)
        valid_left: cnp.ndarray(cnp.int32_t, ndim=1)
        valid_widths: cnp.ndarray(cnp.int32_t, ndim=1)
        chunk_offsets: cnp.ndarray(cnp.int32_t, ndim=1)
        repeated_offsets: cnp.ndarray(cnp.int32_t, ndim=1)
        repeated_left: cnp.ndarray(cnp.int32_t, ndim=1)
        intra_offsets: cnp.ndarray(cnp.int32_t, ndim=1)
        rows_arr: cnp.ndarray(cnp.int32_t, ndim=1)
        cols_arr: cnp.ndarray(cnp.int32_t, ndim=1)
        data_arr: cnp.ndarray(cnp.int32_t, ndim=1)
        chunk_nnz: cython.int
        row_chunks: list
        col_chunks: list
        data_chunks: list
        peak_names: list
        peak_data: list
        peak_names_append: object
        peak_data_append: object

        import pandas as pd
        from scipy import sparse
        import anndata as ad

        peak_names = []
        peak_data = []
        row_chunks = []
        col_chunks = []
        data_chunks = []
        peak_names_append = peak_names.append
        peak_data_append = peak_data.append

        barcode_items = sorted(self.barcode_dict.items(), key=itemgetter(1))
        barcodes = [b.decode() if isinstance(b, (bytes, bytearray)) else str(b) for b, _ in barcode_items]
        n_cells = len(barcodes)
        if n_cells:
            barcode_ids = np.fromiter((barcode_id for _, barcode_id in barcode_items),
                                      dtype=np.int32,
                                      count=n_cells)
            barcode_id_to_row = np.full(int(barcode_ids[-1]) + 1, -1, dtype=np.int32)
            barcode_id_to_row[barcode_ids] = np.arange(n_cells, dtype=np.int32)
        else:
            barcode_id_to_row = np.zeros(0, dtype=np.int32)

        regions.sort()
        peak_counter = 0

        for chrom in sorted(regions.regions.keys()):
            if chrom not in self.locations:
                continue

            regions_c = regions.regions[chrom]
            if not regions_c:
                continue

            n_regions_c = len(regions_c)
            peak_starts = np.empty(n_regions_c, dtype=np.int32)
            peak_ends = np.empty(n_regions_c, dtype=np.int32)
            peak_base = peak_counter
            chrom_str = chrom.decode() if isinstance(chrom, (bytes, bytearray)) else str(chrom)
            for peak_idx, (start, end) in enumerate(regions_c):
                peak_counter += 1
                peak_names_append(f"peak_{peak_counter}")
                peak_data_append((chrom_str, start, end))
                peak_starts[peak_idx] = start
                peak_ends[peak_idx] = end

            fragment_locs = self.locations[chrom]
            if len(fragment_locs) == 0:
                continue

            fragment_barcodes = self.barcodes[chrom]
            frag_starts = fragment_locs['l']
            frag_ends = fragment_locs['r']
            frag_counts = fragment_locs['c']
            frag_rows = barcode_id_to_row[fragment_barcodes]
            left_idx = np.searchsorted(peak_ends, frag_starts, side='right').astype(np.int32, copy=False)
            right_idx = np.searchsorted(peak_starts, frag_ends, side='left').astype(np.int32, copy=False)
            widths = right_idx - left_idx
            valid_mask = np.logical_and(frag_rows >= 0, widths > 0)
            chunk_nnz = int(widths[valid_mask].sum())
            if chunk_nnz:
                valid_rows = frag_rows[valid_mask]
                valid_counts = frag_counts[valid_mask].astype(np.int32, copy=False)
                valid_left = left_idx[valid_mask]
                valid_widths = widths[valid_mask]
                chunk_rows = np.repeat(valid_rows, valid_widths)
                chunk_data = np.repeat(valid_counts, valid_widths)
                chunk_offsets = np.empty(valid_widths.shape[0] + 1, dtype=np.int32)
                chunk_offsets[0] = 0
                np.cumsum(valid_widths, out=chunk_offsets[1:])
                repeated_offsets = np.repeat(chunk_offsets[:-1], valid_widths)
                repeated_left = np.repeat(valid_left, valid_widths)
                intra_offsets = np.arange(chunk_nnz, dtype=np.int32) - repeated_offsets
                chunk_cols = peak_base + repeated_left + intra_offsets
                row_chunks.append(chunk_rows)
                col_chunks.append(chunk_cols)
                data_chunks.append(chunk_data)

        n_peaks = peak_counter
        obs = pd.DataFrame(index=barcodes)
        var = pd.DataFrame(peak_data, columns=['chrom', 'start', 'end'], index=peak_names)
        if row_chunks:
            rows_arr = np.concatenate(row_chunks)
            cols_arr = np.concatenate(col_chunks)
            data_arr = np.concatenate(data_chunks)
            x = sparse.csr_matrix((data_arr, (rows_arr, cols_arr)), shape=(n_cells, n_peaks), dtype=np.int32)
        else:
            x = sparse.csr_matrix((n_cells, n_peaks), dtype=np.int32)
        adata_peaks_loop = ad.AnnData(X=x, obs=obs, var=var)
        return adata_peaks_loop

    def return_anndata(self, regions):
        """
        Build barcode × peak AnnData.

        Parameters
        ----------
        regions : MACS3.Signal.Region.Regions
            Sorted region collection whose intervals should be excluded.
            A merged copy is used so the sweep operates on
            non-overlapping or adjacent-collapsed intervals.
        """
        merged_regions: Regions
        chrom: bytes

        merged_regions = Regions()
        for chrom in sorted(regions.regions.keys()):
            merged_regions.regions[chrom] = regions.regions[chrom][:]
            merged_regions.total += len(merged_regions.regions[chrom])
        merged_regions.merge_overlap()
        return self._two_pointer_sweep(merged_regions)
        
    @cython.ccall
    def sample_percent(self,
                       percent: cython.float,
                       seed: cython.int = -1):
        """Down-sample fragments in place so counts reflect a given percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of total counts to keep per chromosome between 0 and 1 (inclusive).
        seed : int, optional
            Deterministic seed for the RNG; a negative value uses NumPy's global state.
        
        Notes
        -----
        Fragments are sampled proportionally to their counts by expanding to an index
        vector, shuffling, and collapsing counts for the retained entries. Aggregate
        statistics are recomputed and the result is resorted.
        """
        k: bytes
        loc: cnp.ndarray
        bar: cnp.ndarray
        counts: cnp.ndarray
        n: cython.uint
        n_sample: cython.uint
        idx_flat: cnp.ndarray
        unique_idx: cnp.ndarray
        new_counts: cnp.ndarray
        new_locs: cnp.ndarray
        new_bars: cnp.ndarray

        assert 0.0 <= percent <= 1.0, "percent must be in [0, 1]"
        chrnames = sorted(self.get_chr_names())

        # Setup shuffling logic like PETrackI
        if seed >= 0:
            info(f"#   A random seed {seed} has been used")
            rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
            rs_shuffle = rs.shuffle
        else:
            rs_shuffle = np.random.shuffle

        self.length = 0
        self.total = 0
        self.average_template_length = 0.0

        for k in chrnames:
            loc = self.locations[k]
            bar = self.barcodes[k]
            counts = loc['c']
            n = int(counts.sum())
            n_sample = int(round(n * percent))
            if n == 0 or n_sample == 0:
                self.locations[k] = loc[:0]
                self.barcodes[k] = bar[:0]
                self.size[k] = 0
                continue

            # Flatten: build an array of indices into loc, repeated by count
            idx_flat = np.repeat(np.arange(len(loc)), counts)
            rs_shuffle(idx_flat)
            idx_flat = idx_flat[:n_sample]

            # Recount: count how many times each index is chosen
            unique_idx, new_counts = np.unique(idx_flat, return_counts=True)
            # Compose new arrays
            new_locs = loc[unique_idx].copy()
            new_locs['c'] = new_counts
            new_bars = bar[unique_idx].copy()
            self.locations[k] = new_locs
            self.barcodes[k] = new_bars
            self.size[k] = len(new_locs)
            self.length += np.sum((new_locs['r'] - new_locs['l']) * new_locs['c'])
            self.total += np.sum(new_locs['c'])

        if self.total > 0:
            self.average_template_length = float(self.length) / self.total
        else:
            self.average_template_length = 0.0
        self.sort()
        return

    @cython.ccall
    def sample_percent_copy(self,
                            percent: cython.float,
                            seed: cython.int = -1):
        """Return a down-sampled copy whose counts reflect a given percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of total counts to keep per chromosome between 0 and 1 (inclusive).
        seed : int, optional
            Deterministic seed for the RNG; a negative value uses NumPy's global state.
        
        Returns
        -------
        PETrackII
            New track containing the sampled fragments with metadata preserved.
        
        Notes
        -----
        Fragments are sampled proportionally to their counts and the returned track is
        sorted with reference lengths copied from the source track.
        """
        k: bytes
        loc: cnp.ndarray
        bar: cnp.ndarray
        counts: cnp.ndarray
        n: cython.uint
        n_sample: cython.uint
        idx_flat: cnp.ndarray
        unique_idx: cnp.ndarray
        new_counts: cnp.ndarray
        new_locs: cnp.ndarray
        new_bars: cnp.ndarray

        assert 0.0 <= percent <= 1.0, "percent must be in [0, 1]"
        chrnames = sorted(self.get_chr_names())

        # Setup shuffling logic like PETrackI
        if seed >= 0:
            info(f"#   A random seed {seed} has been used")
            rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
            rs_shuffle = rs.shuffle
        else:
            rs_shuffle = np.random.shuffle

        ret = PETrackII(anno=self.annotation, buffer_size=self.buffer_size)
        ret.barcode_dict = dict(self.barcode_dict)
        ret.barcode_last_n = self.barcode_last_n

        ret.length = 0
        ret.total = 0

        for k in chrnames:
            loc = self.locations[k]
            bar = self.barcodes[k]
            counts = loc['c']
            n = int(counts.sum())
            n_sample = int(round(n * percent))
            if n == 0 or n_sample == 0:
                ret.locations[k] = loc[:0]
                ret.barcodes[k] = bar[:0]
                ret.size[k] = 0
                ret.buf_size[k] = 0
                continue

            idx_flat = np.repeat(np.arange(len(loc)), counts)
            rs_shuffle(idx_flat)
            idx_flat = idx_flat[:n_sample]
            unique_idx, new_counts = np.unique(idx_flat, return_counts=True)
            new_locs = loc[unique_idx].copy()
            new_locs['c'] = new_counts
            new_bars = bar[unique_idx].copy()
            ret.locations[k] = new_locs
            ret.barcodes[k] = new_bars
            ret.size[k] = len(new_locs)
            ret.buf_size[k] = len(new_locs)
            ret.length += np.sum((new_locs['r'] - new_locs['l']) * new_locs['c'])
            ret.total += np.sum(new_locs['c'])

        if ret.total > 0:
            ret.average_template_length = float(ret.length) / ret.total
        else:
            ret.average_template_length = 0.0
        ret.set_rlengths(self.get_rlengths())
        ret.sort()
        return ret

    @cython.ccall
    def sample_num(self,
                   samplesize: cython.ulong,
                   seed: cython.int = -1):
        """Down-sample fragments in place so total counts approximate ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target total count across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent`.
        
        Notes
        -----
        The method converts ``samplesize`` into a sampling fraction using the current
        total count and reuses :meth:`sample_percent`.
        """
        chrnames: set
        n_total: cython.uint
        chr_totals: dict
        k: bytes
        percent: cython.float

        chrnames = self.get_chr_names()
        n_total = 0
        chr_totals = {}
        for k in chrnames:
            chr_totals[k] = self.locations[k]['c'].sum()
            n_total += chr_totals[k]
        percent = 0.0 if n_total == 0 else min(samplesize / n_total, 1.0)
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def sample_num_copy(self,
                        samplesize: cython.ulong,
                        seed: cython.int = -1):
        """Return a down-sampled copy whose total counts approximate ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target total count across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent_copy`.
        
        Returns
        -------
        PETrackII
            New track containing the sampled fragments.
        """
        chrnames: set
        n_total: cython.uint
        chr_totals: dict
        k: bytes
        percent: cython.float

        chrnames = self.get_chr_names()
        n_total = 0
        chr_totals = {}
        for k in chrnames:
            chr_totals[k] = self.locations[k]['c'].sum()
            n_total += chr_totals[k]
        percent = 0.0 if n_total == 0 else min(samplesize / n_total, 1.0)
        return self.sample_percent_copy(percent, seed)

    @cython.ccall
    def pileup_bdg_hmmr(self,
                        mapping: list,
                        baseline_value: cython.float = 0.0) -> list:
        """Generate HMMRATAC-style pileups for every chromosome.
        
        Parameters
        ----------
        mapping : list
            Weight mapping produced by HMMRATAC EM training describing the short,
            mono-, di-, and tri-nucleosomal signals.
        baseline_value : float, optional
            Reserved parameter for API compatibility; not currently applied.
        
        Returns
        -------
        list
            List of dictionaries mirroring ``mapping`` where each dictionary maps
            chromosome names to pileup arrays returned by
            :func:`pileup_from_LR_hmmratac`.
        """
        ret_pileup: list
        i: cython.uint
        chroms: set
        chrom: bytes
        locs: cnp.ndarray
        counts: cnp.ndarray
        LR_expanded: cnp.ndarray
        idx: cnp.ndarray

        ret_pileup = []
        for i in range(len(mapping)):
            ret_pileup.append({})

        chroms = self.get_chr_names()
        for i in range(len(mapping)):
            for chrom in sorted(chroms):
                locs = self.locations[chrom]
                counts = locs['c']
                # Efficient numpy "explode"
                if locs.shape[0] == 0 or counts.sum() == 0:
                    LR_expanded = np.zeros((0,),
                                           dtype=[('l', 'i4'), ('r', 'i4')])
                else:
                    idx = np.repeat(np.arange(locs.shape[0]), counts)
                    LR_expanded = np.empty((len(idx),),
                                           dtype=[('l', 'i4'), ('r', 'i4')])
                    LR_expanded['l'] = locs['l'][idx]
                    LR_expanded['r'] = locs['r'][idx]
                ret_pileup[i][chrom] = pileup_from_LR_hmmratac(LR_expanded, mapping[i])
        return ret_pileup
