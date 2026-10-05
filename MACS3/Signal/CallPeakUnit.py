# cython: language_level=3
# Time-stamp: <2025-10-16 17:09:16 Tao Liu>

"""Compute MACS3 peak-calling scores and helper statistics.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import struct
from tempfile import mkstemp
import os
from stat import S_ISREG
import threading
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
import cython
import cython.cimports.numpy as cnp
# from numpy cimport int32_t, int64_t, float32_t, float64_t
from cython.cimports.cpython import bool

# ------------------------------------
# C lib
# ------------------------------------
from cython.cimports.libc.stdio import FILE, fopen, fprintf, fclose
from cython.cimports.libc.math import exp, log10, log1p, erf, sqrt
from cython.cimports.libc.stdlib import calloc, free, realloc
from cython.cimports.libc.string import memcpy
from cython.cimports.posix.unistd import close
from cython.cimports.MACS3.Signal.cmalloc import (macs3_release_free_heap,
                                                  macs3_advise_huge_pages)
from cython.cimports.MACS3.Signal.cprefetch import macs3_prefetch_w

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.SignalProcessing import maxima, enforce_peakyness
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
from MACS3.Signal.Prob import (poisson_cdf, poisson_cdf_Q_log10_many,
                                poisson_cdf_Q_log10_prepare,
                                poisson_cdf_Q_log10_range)
from MACS3.Utilities.Logger import logging

logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info
# --------------------------------------------
# cached pscore function and LR_asym functions
# --------------------------------------------
logLR_dict = {}

# One slot of the p-score cache. The key packs the observed count
# (high 32 bits) with the float32 bits of lambda (low 32 bits), so a
# lookup needs no Python object. ``first`` and ``length`` hold the
# p-score histogram of one pass of __cal_pvalue_qvalue_table: the
# first segment that used the slot, as (chromosome index << 32) +
# segment index within the chromosome (-1: not used in this pass), and
# the summed length of the segments that used it. ``used`` is 0 for an
# empty slot and 1 for a slot holding a key; the p-score of a slot
# index_unscored added is not computed.
PScoreSlot = cython.struct(key=cython.ulonglong,
                           pscore=cython.float,
                           used=cython.int,
                           first=cython.longlong,
                           length=cython.longlong)

# A key index_deferred or merge_histogram added without its p-score,
# and the slot it was put in, which a later grow may have moved.
PScorePending = cython.struct(key=cython.ulonglong,
                              slot=cython.ulonglong)

PSCORE_CACHE_INIT_SIZE: cython.Py_ssize_t = 1 << 16

# score_pending hands out the pending keys to its threads this many at
# a time.
PENDING_CHUNK = cython.declare(cython.Py_ssize_t, 8192)


@cython.cfunc
@cython.inline
@cython.profile(False)
@cython.linetrace(False)
@cython.exceptval(check=False)
def pscore_key(observed: cython.int,
               lam: cython.float) -> cython.ulonglong:
    """Pack ``observed`` and the bits of float32 ``lam`` in 64 bits."""
    bits: cython.uint

    memcpy(cython.address(bits), cython.address(lam), 4)
    return (cython.cast(cython.ulonglong,
                        cython.cast(cython.uint, observed)) << 32) | bits


@cython.cfunc
@cython.inline
@cython.nogil
@cython.profile(False)
@cython.linetrace(False)
@cython.exceptval(check=False)
def pscore_hash(key: cython.ulonglong) -> cython.ulonglong:
    """Mix the bits of ``key``; the caller masks the result."""
    h: cython.ulonglong

    h = key * ((cython.cast(cython.ulonglong, 0x9E3779B9) << 32) |
               cython.cast(cython.ulonglong, 0x7F4A7C15))
    return h ^ (h >> 32)


@cython.final
@cython.cclass
class PScoreCache:
    """Open-addressing hash table of ``-log10`` Poisson upper-tail
    p-scores keyed by (observed, lambda), with linear probing.

    Each p-score is computed once by ``poisson_cdf`` and stored as
    float32, exactly as the dict of packed (int, float32 bits) keys it
    replaces does.
    """
    slots: cython.pointer(PScoreSlot)
    mask: cython.ulonglong      # capacity - 1; capacity is a power of 2
    n: cython.Py_ssize_t        # number of used slots
    pending: cython.pointer(PScorePending)  # p-scores left to compute
    n_pending: cython.Py_ssize_t
    pending_capacity: cython.Py_ssize_t
    pending_max_observed: cython.uint  # largest observed count pending

    def __cinit__(self):
        self.slots = cython.cast(cython.pointer(PScoreSlot),
                                 calloc(PSCORE_CACHE_INIT_SIZE,
                                        cython.sizeof(PScoreSlot)))
        if self.slots == cython.NULL:
            raise MemoryError()
        self.mask = PSCORE_CACHE_INIT_SIZE - 1
        self.n = 0
        self.pending = cython.NULL
        self.n_pending = 0
        self.pending_capacity = 0
        self.pending_max_observed = 0

    def __dealloc__(self):
        free(self.slots)
        free(self.pending)

    @cython.cfunc
    @cython.inline
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def index(self,
              observed: cython.int,
              lam: cython.float) -> cython.Py_ssize_t:
        """Return the slot holding the p-score of (observed, lam),
        computing it on a miss."""
        key: cython.ulonglong
        i: cython.ulonglong

        key = pscore_key(observed, lam)
        i = pscore_hash(key) & self.mask
        while self.slots[i].used:
            if self.slots[i].key == key:
                return i
            i = (i + 1) & self.mask
        return self.insert(i, key, observed, lam)

    @cython.cfunc
    @cython.inline
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def index_deferred(self,
                       observed: cython.int,
                       lam: cython.float) -> cython.Py_ssize_t:
        """Return the slot of (observed, lam), as index does, but on
        a miss leave its p-score to compute_pending, which must run
        before the slot's p-score is read."""
        key: cython.ulonglong
        i: cython.ulonglong

        key = pscore_key(observed, lam)
        i = pscore_hash(key) & self.mask
        while self.slots[i].used:
            if self.slots[i].key == key:
                return i
            i = (i + 1) & self.mask
        if observed < 0 or not lam > 0:
            # poisson_cdf refuses these; let it raise now, as it did
            return self.insert(i, key, observed, lam)
        return self.insert_deferred(i, key)

    @cython.cfunc
    @cython.inline
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def index_unscored(self,
                       observed: cython.int,
                       lam: cython.float) -> cython.Py_ssize_t:
        """Return the slot of (observed, lam), as index does, but on
        a miss add the key without its p-score, which is then never
        computed in this cache. For a worker process whose
        histogram's keys merge_histogram scores in another."""
        key: cython.ulonglong
        i: cython.ulonglong

        key = pscore_key(observed, lam)
        i = pscore_hash(key) & self.mask
        while self.slots[i].used:
            if self.slots[i].key == key:
                return i
            i = (i + 1) & self.mask
        if observed < 0 or not lam > 0:
            # poisson_cdf refuses these; let it raise now, as it did
            return self.insert(i, key, observed, lam)
        if (self.n + 1) * 2 > cython.cast(cython.Py_ssize_t, self.mask + 1):
            self.grow()
            i = pscore_hash(key) & self.mask
            while self.slots[i].used:
                i = (i + 1) & self.mask
        self.slots[i].key = key
        self.slots[i].pscore = 0
        self.slots[i].used = 1
        self.slots[i].first = -1
        self.slots[i].length = 0
        self.n += 1
        return i

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def insert_deferred(self,
                        i: cython.ulonglong,
                        key: cython.ulonglong) -> cython.Py_ssize_t:
        """Store ``key`` in the empty slot ``i``, growing the table
        first if needed, as insert does, and add it to the keys whose
        p-score compute_pending computes."""
        if self.n_pending == self.pending_capacity:
            self.reserve_pending(1)
        if (self.n + 1) * 2 > cython.cast(cython.Py_ssize_t, self.mask + 1):
            self.grow()
            i = pscore_hash(key) & self.mask
            while self.slots[i].used:
                i = (i + 1) & self.mask
        self.slots[i].key = key
        self.slots[i].pscore = 0
        self.slots[i].used = 1
        self.slots[i].first = -1
        self.slots[i].length = 0
        self.n += 1
        self.pending[self.n_pending].key = key
        self.pending[self.n_pending].slot = i
        self.n_pending += 1
        if cython.cast(cython.uint, key >> 32) > self.pending_max_observed:
            self.pending_max_observed = cython.cast(cython.uint, key >> 32)
        return i

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def reserve_pending(self, extra: cython.Py_ssize_t) -> cython.int:
        """Make room for ``extra`` more pending keys, doubling the
        buffer, from at least 4096 entries, until they fit."""
        p: cython.pointer(PScorePending)
        capacity: cython.Py_ssize_t

        if self.n_pending + extra <= self.pending_capacity:
            return 0
        capacity = max(2 * self.pending_capacity, 4096)
        while capacity < self.n_pending + extra:
            capacity *= 2
        p = cython.cast(cython.pointer(PScorePending),
                        realloc(self.pending,
                                capacity * cython.sizeof(PScorePending)))
        if p == cython.NULL:
            raise MemoryError()
        self.pending = p
        self.pending_capacity = capacity
        return 0

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    def compute_pending(self):
        """Compute the p-score of every key index_deferred added
        since the last call, all in one poisson_cdf_Q_log10_many call,
        which gives the values poisson_cdf gives one by one, and store
        each in its slot as insert would."""
        observed: cnp.ndarray
        lams: cnp.ndarray
        values: cnp.ndarray
        o_ptr: cython.pointer(cython.uint)
        l_ptr: cython.pointer(cython.double)
        v_ptr: cython.pointer(cython.double)
        n: cython.Py_ssize_t
        e: cython.Py_ssize_t
        key: cython.ulonglong
        i: cython.ulonglong
        bits: cython.uint
        lam: cython.float

        n = self.n_pending
        if n == 0:
            return
        self.n_pending = 0
        self.pending_max_observed = 0
        observed = np.empty(n, dtype="u4")
        lams = np.empty(n, dtype="f8")
        values = np.empty(n, dtype="f8")
        o_ptr = cython.cast(cython.pointer(cython.uint), observed.data)
        l_ptr = cython.cast(cython.pointer(cython.double), lams.data)
        v_ptr = cython.cast(cython.pointer(cython.double), values.data)
        for e in range(n):
            key = self.pending[e].key
            o_ptr[e] = cython.cast(cython.uint, key >> 32)
            bits = cython.cast(cython.uint, key & 0xFFFFFFFF)
            memcpy(cython.address(lam), cython.address(bits), 4)
            l_ptr[e] = lam
        poisson_cdf_Q_log10_many(observed, lams, values)
        for e in range(n):
            key = self.pending[e].key
            i = self.pending[e].slot
            if not (self.slots[i].used and self.slots[i].key == key):
                # the table grew after the key was added
                i = pscore_hash(key) & self.mask
                while not (self.slots[i].used and self.slots[i].key == key):
                    i = (i + 1) & self.mask
            self.slots[i].pscore = -1.0 * v_ptr[e]
        return

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def insert(self,
               i: cython.ulonglong,
               key: cython.ulonglong,
               observed: cython.int,
               lam: cython.float) -> cython.Py_ssize_t:
        """Compute the p-score of (observed, lam) and store it in the
        empty slot ``i``, growing the table first if needed."""
        val: cython.float

        # calculate and cache
        val = -1.0 * poisson_cdf(observed, lam, False, True)
        if (self.n + 1) * 2 > cython.cast(cython.Py_ssize_t, self.mask + 1):
            self.grow()
            i = pscore_hash(key) & self.mask
            while self.slots[i].used:
                i = (i + 1) & self.mask
        self.slots[i].key = key
        self.slots[i].pscore = val
        self.slots[i].used = 1
        self.slots[i].first = -1
        self.slots[i].length = 0
        self.n += 1
        return i

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def grow(self) -> cython.int:
        """Double the capacity and re-insert every used slot. A table
        of 4 MB or more is advised to use huge pages: the re-insertion
        writes every page of it, and the lookups that follow land at
        random in it."""
        old: cython.pointer(PScoreSlot)
        new: cython.pointer(PScoreSlot)
        old_capacity: cython.ulonglong
        new_mask: cython.ulonglong
        i: cython.ulonglong
        j: cython.ulonglong

        old = self.slots
        old_capacity = self.mask + 1
        new = cython.cast(cython.pointer(PScoreSlot),
                          calloc(old_capacity * 2, cython.sizeof(PScoreSlot)))
        if new == cython.NULL:
            raise MemoryError()
        macs3_advise_huge_pages(new, old_capacity * 2 * cython.sizeof(PScoreSlot))
        new_mask = old_capacity * 2 - 1
        for j in range(old_capacity):
            if old[j].used:
                i = pscore_hash(old[j].key) & new_mask
                while new[i].used:
                    i = (i + 1) & new_mask
                new[i] = old[j]
        self.slots = new
        self.mask = new_mask
        free(old)
        return 0

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    def reset_histogram(self):
        """Mark every slot as unused by the coming histogram pass."""
        j: cython.ulonglong

        for j in range(self.mask + 1):
            self.slots[j].first = -1
            self.slots[j].length = 0

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    def take_histogram(self, keys: cython.bint = False) -> tuple:
        """Return the p-score, first segment and summed length of
        every slot used since the last reset, as three arrays, and
        reset those slots.

        With ``keys``, the first array holds each slot's key (uint64)
        instead of its p-score, for merge_histogram in another
        process.
        """
        slot: cython.pointer(PScoreSlot)
        s: cython.Py_ssize_t
        n_used: cython.Py_ssize_t
        values: cnp.ndarray
        firsts: cnp.ndarray
        lengths: cnp.ndarray
        ps_ptr: cython.pointer(cython.float)
        ky_ptr: cython.pointer(cython.ulonglong)
        fi_ptr: cython.pointer(cython.longlong)
        ln_ptr: cython.pointer(cython.longlong)

        values = np.empty(self.n, dtype="u8" if keys else "f4")
        firsts = np.empty(self.n, dtype="i8")
        lengths = np.empty(self.n, dtype="i8")
        ps_ptr = cython.cast(cython.pointer(cython.float), values.data)
        ky_ptr = cython.cast(cython.pointer(cython.ulonglong), values.data)
        fi_ptr = cython.cast(cython.pointer(cython.longlong), firsts.data)
        ln_ptr = cython.cast(cython.pointer(cython.longlong), lengths.data)
        n_used = 0
        for s in range(cython.cast(cython.Py_ssize_t, self.mask + 1)):
            slot = cython.address(self.slots[s])
            if slot.used and slot.first >= 0:
                if keys:
                    ky_ptr[n_used] = slot.key
                else:
                    ps_ptr[n_used] = slot.pscore
                fi_ptr[n_used] = slot.first
                ln_ptr[n_used] = slot.length
                n_used += 1
                slot.first = -1
                slot.length = 0
        return (values[:n_used], firsts[:n_used], lengths[:n_used])

    @cython.cfunc
    @cython.nogil
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(check=False)
    def merge_histogram_some(self,
                             keys: cython.pointer(cython.ulonglong),
                             firsts: cython.pointer(cython.longlong),
                             lengths: cython.pointer(cython.longlong),
                             start: cython.Py_ssize_t,
                             end: cython.Py_ssize_t) -> cython.Py_ssize_t:
        """Add entries start to end - 1 of a histogram to the slots of
        their keys, until a new key would make the table grow. Return
        the index of the first entry not added. The pending buffer
        must have room for every new key."""
        e: cython.Py_ssize_t
        key: cython.ulonglong
        i: cython.ulonglong
        slot: cython.pointer(PScoreSlot)

        for e in range(start, end):
            key = keys[e]
            i = pscore_hash(key) & self.mask
            while self.slots[i].used and self.slots[i].key != key:
                i = (i + 1) & self.mask
            slot = cython.address(self.slots[i])
            if slot.used:
                if slot.first < 0 or firsts[e] < slot.first:
                    slot.first = firsts[e]
                slot.length += lengths[e]
                continue
            if (self.n + 1) * 2 > cython.cast(cython.Py_ssize_t, self.mask + 1):
                return e
            slot.key = key
            slot.pscore = 0
            slot.used = 1
            slot.first = firsts[e]
            slot.length = lengths[e]
            self.n += 1
            self.pending[self.n_pending].key = key
            self.pending[self.n_pending].slot = i
            self.n_pending += 1
            if cython.cast(cython.uint, key >> 32) > self.pending_max_observed:
                self.pending_max_observed = cython.cast(cython.uint, key >> 32)
        return end

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    def merge_histogram(self, keys: cnp.ndarray, firsts: cnp.ndarray,
                        lengths: cnp.ndarray):
        """Add a histogram that take_histogram(True) returned in
        another process, keys (uint64), first segments and summed
        lengths (int64), to this table's: a slot keeps the smaller
        first segment and sums the lengths. A key this table does not
        hold yet is added without its p-score, which score_pending
        computes. The GIL is released while the entries are added,
        except to grow the table."""
        k_ptr: cython.pointer(cython.ulonglong)
        f_ptr: cython.pointer(cython.longlong)
        l_ptr: cython.pointer(cython.longlong)
        e: cython.Py_ssize_t
        n: cython.Py_ssize_t

        assert keys.dtype == np.uint64 and firsts.dtype == np.int64
        assert lengths.dtype == np.int64
        assert (keys.flags.c_contiguous and firsts.flags.c_contiguous and
                lengths.flags.c_contiguous)
        n = keys.shape[0]
        assert firsts.shape[0] == n and lengths.shape[0] == n
        self.reserve_pending(n)
        k_ptr = cython.cast(cython.pointer(cython.ulonglong), keys.data)
        f_ptr = cython.cast(cython.pointer(cython.longlong), firsts.data)
        l_ptr = cython.cast(cython.pointer(cython.longlong), lengths.data)
        e = 0
        while True:
            with cython.nogil:
                e = self.merge_histogram_some(k_ptr, f_ptr, l_ptr, e, n)
            if e == n:
                break
            self.grow()
        return

    @cython.cfunc
    @cython.nogil
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(check=False)
    def decode_pending_some(self,
                            observed: cython.pointer(cython.uint),
                            lams: cython.pointer(cython.double),
                            start: cython.Py_ssize_t,
                            end: cython.Py_ssize_t) -> cython.int:
        """Set observed[e] and lams[e] to the (observed, lambda) of
        pending key e, for start <= e < end."""
        e: cython.Py_ssize_t
        key: cython.ulonglong
        bits: cython.uint
        lam: cython.float

        for e in range(start, end):
            key = self.pending[e].key
            observed[e] = cython.cast(cython.uint, key >> 32)
            bits = cython.cast(cython.uint, key)    # the low 32 bits
            memcpy(cython.address(lam), cython.address(bits), 4)
            lams[e] = lam
        return 0

    @cython.cfunc
    @cython.nogil
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(check=False)
    def store_pending_some(self,
                           values: cython.pointer(cython.double),
                           pscores: cython.pointer(cython.float),
                           firsts: cython.pointer(cython.longlong),
                           lengths: cython.pointer(cython.longlong),
                           start: cython.Py_ssize_t,
                           end: cython.Py_ssize_t) -> cython.int:
        """Store -values[e] as the p-score of pending key e, for start
        <= e < end, as compute_pending does. With ``pscores`` not
        NULL, also set entry e of the three arrays to the slot's
        p-score, first segment and summed length, and reset the
        slot's, as take_histogram does. Threads may store disjoint
        ranges at once: each key has its own slot, and the table does
        not change shape meanwhile."""
        e: cython.Py_ssize_t
        key: cython.ulonglong
        i: cython.ulonglong
        slot: cython.pointer(PScoreSlot)

        for e in range(start, end):
            key = self.pending[e].key
            # a grow after the key was added moves its slot
            i = pscore_hash(key) & self.mask
            while self.slots[i].key != key or not self.slots[i].used:
                i = (i + 1) & self.mask
            slot = cython.address(self.slots[i])
            slot.pscore = -1.0 * values[e]
            if pscores != cython.NULL:
                pscores[e] = slot.pscore
                firsts[e] = slot.first
                lengths[e] = slot.length
                slot.first = -1
                slot.length = 0
        return 0

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    def score_pending(self, n_threads: cython.int, histogram: cython.bint):
        """Compute the p-score of every pending key, as compute_pending
        does, on up to ``n_threads`` threads, the calling one included.

        With ``histogram``, every slot used since the last reset must
        be pending; return their p-scores, first segments and summed
        lengths, in the order of the pending keys, and reset them, as
        take_histogram does in its order.
        """
        observed: cnp.ndarray
        lams: cnp.ndarray
        values: cnp.ndarray
        pscores: cnp.ndarray = None
        firsts: cnp.ndarray = None
        lengths: cnp.ndarray = None
        n: cython.Py_ssize_t
        chunk: cython.Py_ssize_t
        threads: list
        state: list

        n = self.n_pending
        observed = np.empty(n, dtype="u4")
        lams = np.empty(n, dtype="f8")
        values = np.empty(n, dtype="f8")
        if histogram:
            pscores = np.empty(n, dtype="f4")
            firsts = np.empty(n, dtype="i8")
            lengths = np.empty(n, dtype="i8")
        # every pending key passed index_deferred's or a worker's
        # index_unscored's check, so poisson_cdf takes its lambda
        chunk = PENDING_CHUNK
        if (n_threads < 2 or n <= chunk or
                not poisson_cdf_Q_log10_prepare(self.pending_max_observed)):
            n_threads = 1
        # the next pending key not handed out, a lock, and the first
        # exception a thread raised
        state = [0, threading.Lock(), None]
        threads = [threading.Thread(target=score_pending_chunks,
                                    args=(self, observed, lams, values,
                                          pscores, firsts, lengths,
                                          chunk, state))
                   for _ in range(min(n_threads, (n + chunk - 1) // chunk) - 1)]
        for t in threads:
            t.start()
        try:
            score_pending_chunks(self, observed, lams, values,
                                 pscores, firsts, lengths, chunk, state)
        finally:
            for t in threads:
                t.join()
        if state[2] is not None:
            raise state[2]
        self.n_pending = 0
        self.pending_max_observed = 0
        if histogram:
            return (pscores, firsts, lengths)
        return None


def score_pending_chunks(cache: PScoreCache,
                         observed: cnp.ndarray,
                         lams: cnp.ndarray,
                         values: cnp.ndarray,
                         pscores,
                         firsts,
                         lengths,
                         chunk: cython.Py_ssize_t,
                         state: list):
    """Thread side of PScoreCache.score_pending: take the next
    ``chunk`` pending keys until none is left, and decode them, compute
    their p-scores and store them, without the GIL. ``pscores``,
    ``firsts`` and ``lengths`` are the histogram arrays, or all None.
    The first exception raised stops every thread at its next chunk
    and is kept in state[2]."""
    n: cython.Py_ssize_t = observed.shape[0]
    start: cython.Py_ssize_t
    end: cython.Py_ssize_t
    o_ptr: cython.pointer(cython.uint)
    l_ptr: cython.pointer(cython.double)
    v_ptr: cython.pointer(cython.double)
    ps_ptr: cython.pointer(cython.float) = cython.NULL
    fi_ptr: cython.pointer(cython.longlong) = cython.NULL
    ln_ptr: cython.pointer(cython.longlong) = cython.NULL
    histogram: cnp.ndarray

    o_ptr = cython.cast(cython.pointer(cython.uint), observed.data)
    l_ptr = cython.cast(cython.pointer(cython.double), lams.data)
    v_ptr = cython.cast(cython.pointer(cython.double), values.data)
    if pscores is not None:
        histogram = pscores
        ps_ptr = cython.cast(cython.pointer(cython.float), histogram.data)
        histogram = firsts
        fi_ptr = cython.cast(cython.pointer(cython.longlong), histogram.data)
        histogram = lengths
        ln_ptr = cython.cast(cython.pointer(cython.longlong), histogram.data)
    lock = state[1]
    try:
        while True:
            with lock:
                if state[2] is not None:
                    return
                start = state[0]
                state[0] = start + chunk
            if start >= n:
                return
            end = min(start + chunk, n)
            with cython.nogil:
                cache.decode_pending_some(o_ptr, l_ptr, start, end)
            poisson_cdf_Q_log10_range(observed, lams, values, start, end)
            with cython.nogil:
                cache.store_pending_some(v_ptr, ps_ptr, fi_ptr, ln_ptr,
                                         start, end)
    except BaseException as e:
        with lock:
            if state[2] is None:
                state[2] = e
        return

pscore_cache = cython.declare(PScoreCache, PScoreCache())


@cython.cfunc
@cython.inline
@cython.profile(False)
@cython.linetrace(False)
def get_pscore(observed: cython.int,
               expectation: cython.float) -> cython.float:
    """Return cached ``-log10`` Poisson tail probability for
    (``observed``, ``expectation``)."""
    i: cython.Py_ssize_t

    # index() may grow the table, so read slots only after it returns
    i = pscore_cache.index(observed, expectation)
    return pscore_cache.slots[i].pscore


# Sort keys of pscore_order: the high 32 bits order the p-scores from
# the largest down, with 0.0 and -0.0 on PSCORE_ORDER_ZERO and every
# NaN, last, on PSCORE_ORDER_NAN; the low 32 bits are the slot index.
PSCORE_ORDER_ZERO = cython.declare(cython.ulonglong, 0x7FFFFFFF)
PSCORE_ORDER_NAN = cython.declare(cython.ulonglong, 0xFFFFFFFF)
PSCORE_ORDER_INDEX = cython.declare(cython.ulonglong, 0xFFFFFFFF)


@cython.cfunc
def pscore_order(pscores: cnp.ndarray, firsts: cnp.ndarray) -> cnp.ndarray:
    """Return the order in which __cal_pvalue_qvalue_table merges the
    histogram slots (float32 ``pscores``, int64 ``firsts``): by
    descending p-score, as np.lexsort((firsts, -pscores)) orders
    them, except that slots with an identical p-score stay in slot
    order instead of the order of ``first``.

    The merge sums the lengths of equal p-scores and keeps the
    p-score of the first slot of each run, so the order within a run
    matters only where equal p-scores differ in bits, 0.0 and -0.0,
    and for NaN, which equals nothing and so is never merged. The
    0.0/-0.0 run starts with its slot of smallest ``first`` and the
    NaN slots, last, are in order of ``first``, so the merged
    p-scores and lengths are those lexsort's order gives. lexsort
    itself is used when a slot index does not fit in 32 bits.

    One sort of distinct 64-bit integer keys (the high 32 bits from
    the p-score, the low 32 the slot index) replaces lexsort's two
    stable sorts.
    """
    n: cython.Py_ssize_t
    s: cython.Py_ssize_t
    z_start: cython.Py_ssize_t
    z_best: cython.Py_ssize_t
    n_nan: cython.Py_ssize_t
    v: cython.float
    bits: cython.uint
    hi: cython.ulonglong
    keys: cnp.ndarray
    order: cnp.ndarray
    tail: cnp.ndarray
    k_ptr: cython.pointer(cython.ulonglong)
    ps_ptr: cython.pointer(cython.float)
    fi_ptr: cython.pointer(cython.longlong)

    assert pscores.dtype == np.float32 and firsts.dtype == np.int64
    assert pscores.flags.c_contiguous and firsts.flags.c_contiguous
    n = pscores.shape[0]
    if cython.cast(cython.ulonglong, n) > PSCORE_ORDER_INDEX:
        return np.lexsort((firsts, -pscores)).astype(np.intp, copy=False)

    keys = np.empty(n, dtype="u8")
    k_ptr = cython.cast(cython.pointer(cython.ulonglong), keys.data)
    ps_ptr = cython.cast(cython.pointer(cython.float), pscores.data)
    fi_ptr = cython.cast(cython.pointer(cython.longlong), firsts.data)
    n_nan = 0
    for s in range(n):
        v = ps_ptr[s]
        if v != v:
            hi = PSCORE_ORDER_NAN
            n_nan += 1
        elif v == 0:
            hi = PSCORE_ORDER_ZERO
        else:
            memcpy(cython.address(bits), cython.address(v), 4)
            # a positive float's bits grow with it, a negative one's
            # with its magnitude: flip the positive ones below
            # PSCORE_ORDER_ZERO, keep the negative ones above it
            if (bits >> 31) == 0:
                bits ^= 0x7FFFFFFF
            hi = bits
        k_ptr[s] = (hi << 32) | cython.cast(cython.ulonglong, s)

    # the keys are distinct, so any sort gives this order
    keys.sort()

    # keys -> slot indices, finding the zero slot with the smallest
    # first on the way
    z_start = -1
    z_best = -1
    for s in range(n):
        hi = k_ptr[s] >> 32
        k_ptr[s] &= PSCORE_ORDER_INDEX
        if hi == PSCORE_ORDER_ZERO:
            if z_start < 0:
                z_start = s
                z_best = s
            elif fi_ptr[k_ptr[s]] < fi_ptr[k_ptr[z_best]]:
                z_best = s
    if z_best != z_start:
        (k_ptr[z_start], k_ptr[z_best]) = (k_ptr[z_best], k_ptr[z_start])
    order = keys.view(np.intp)
    if n_nan > 1:
        tail = order[n - n_nan:]
        order[n - n_nan:] = tail[np.argsort(firsts[tail], kind="stable")]
    return order


# One slot of the p-score -> q-score table. ``key`` is qscore_key()
# of the p-score, or QSCORE_EMPTY when the slot is free.
QScoreSlot = cython.struct(key=cython.uint,
                           qscore=cython.float)

QSCORE_EMPTY = cython.declare(cython.uint, 0xFFFFFFFF)


@cython.cfunc
@cython.inline
@cython.profile(False)
@cython.linetrace(False)
@cython.exceptval(check=False)
def qscore_key(pscore: cython.float) -> cython.uint:
    """The bits of float32 ``pscore``, with 0.0 and -0.0 sent to 0 and
    every NaN to 0x7FC00000, so that two non-NaN keys are equal exactly
    when a dict keyed by their Python floats treats them as one key
    (``a == b``). Never QSCORE_EMPTY, which is a NaN. (A dict keeps
    NaN keys apart; this table holds them as one. callpeak's p-scores
    are never NaN.)"""
    bits: cython.uint

    if pscore == 0:
        return 0
    if pscore != pscore:
        return 0x7FC00000
    memcpy(cython.address(bits), cython.address(pscore), 4)
    return bits


@cython.final
@cython.cclass
class QScoreTable:
    """The p-score to q-score table of ``CallerFromAlignments``: an
    open-addressing hash table, with linear probing, from float32
    p-score to float32 q-score.

    It holds what the dict ``pqtable`` given the same puts would hold:
    a later put of an equal key overwrites the value, so ``get``
    returns the float32 value the dict returns, and raises the same
    KeyError for a p-score the dict does not hold, without a Python
    object per lookup.
    """
    slots: cython.pointer(QScoreSlot)
    mask: cython.ulonglong      # capacity - 1; capacity is a power of 2
    n: cython.Py_ssize_t        # number of used slots

    def __cinit__(self, n_hint: cython.Py_ssize_t = 0):
        capacity: cython.ulonglong = 16

        while capacity < cython.cast(cython.ulonglong, n_hint) * 2:
            capacity *= 2
        self.slots = cython.NULL
        self.alloc(capacity)
        self.n = 0

    def __dealloc__(self):
        free(self.slots)

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def alloc(self, capacity: cython.ulonglong) -> cython.int:
        """Point ``slots`` at ``capacity`` free slots."""
        j: cython.ulonglong
        slots: cython.pointer(QScoreSlot)

        slots = cython.cast(cython.pointer(QScoreSlot),
                            calloc(capacity, cython.sizeof(QScoreSlot)))
        if slots == cython.NULL:
            raise MemoryError()
        for j in range(capacity):
            slots[j].key = QSCORE_EMPTY
        self.slots = slots
        self.mask = capacity - 1
        return 0

    @cython.cfunc
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1)
    def put(self,
            pscore: cython.float,
            qscore: cython.float) -> cython.int:
        """Map ``pscore`` to ``qscore``, replacing an equal key's."""
        key: cython.uint
        i: cython.ulonglong
        j: cython.ulonglong
        old: cython.pointer(QScoreSlot)
        old_capacity: cython.ulonglong

        if (self.n + 1) * 2 > cython.cast(cython.Py_ssize_t, self.mask + 1):
            old = self.slots
            old_capacity = self.mask + 1
            self.alloc(old_capacity * 2)
            for j in range(old_capacity):
                if old[j].key != QSCORE_EMPTY:
                    i = pscore_hash(old[j].key) & self.mask
                    while self.slots[i].key != QSCORE_EMPTY:
                        i = (i + 1) & self.mask
                    self.slots[i] = old[j]
            free(old)
        key = qscore_key(pscore)
        i = pscore_hash(key) & self.mask
        while self.slots[i].key != QSCORE_EMPTY:
            if self.slots[i].key == key:
                self.slots[i].qscore = qscore
                return 0
            i = (i + 1) & self.mask
        self.slots[i].key = key
        self.slots[i].qscore = qscore
        self.n += 1
        return 0

    @cython.cfunc
    @cython.inline
    @cython.profile(False)
    @cython.linetrace(False)
    @cython.exceptval(-1.0, check=True)
    def get(self, pscore: cython.float) -> cython.float:
        """The q-score of ``pscore``; KeyError(pscore) if absent, as
        a dict lookup raises. (-1 only makes the caller
        check for an exception; the table never holds a negative
        q-score.)"""
        key: cython.uint
        i: cython.ulonglong

        key = qscore_key(pscore)
        i = pscore_hash(key) & self.mask
        while self.slots[i].key != QSCORE_EMPTY:
            if self.slots[i].key == key:
                return self.slots[i].qscore
            i = (i + 1) & self.mask
        raise KeyError(pscore)


@cython.cfunc
def get_logLR_asym(t: tuple) -> cython.float:
    """Return asymmetric ``log10`` likelihood ratio for parameters ``t``."""
    val: cython.float
    x: cython.float
    y: cython.float

    if t in logLR_dict:
        return logLR_dict[t]
    else:
        x = t[0]
        y = t[1]
        # calculate and cache
        if x > y:
            val = (x*(log10(x)-log10(y))+y-x)
        elif x < y:
            val = (x*(-log10(x)+log10(y))-y+x)
        else:
            val = 0
        logLR_dict[t] = val
        return val

# ------------------------------------
# constants
# ------------------------------------


LOG10_E: cython.float = 0.43429448190325176

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.cfunc
def clean_up_ndarray(x: cnp.ndarray):
    """Resize ``x`` to zero length, releasing its underlying buffer."""
    # clean numpy ndarray in two steps
    i: cython.long

    i = x.shape[0] // 2
    x.resize(100000 if i > 100000 else i, refcheck=False)
    x.resize(0, refcheck=False)
    return


@cython.cfunc
@cython.inline
@cython.nogil
@cython.profile(False)
@cython.linetrace(False)
@cython.exceptval(check=False)
def pair_pv_i4_f4(t_p: cython.p_char, t_ps: cython.Py_ssize_t,
                  t_v: cython.p_char, t_vs: cython.Py_ssize_t,
                  lt: cython.Py_ssize_t,
                  c_p: cython.p_char, c_ps: cython.Py_ssize_t,
                  c_v: cython.p_char, c_vs: cython.Py_ssize_t,
                  lc: cython.Py_ssize_t,
                  ret_p: cython.pointer(cython.int),
                  ret_t: cython.pointer(cython.float),
                  ret_c: cython.pointer(cython.float)) -> cython.Py_ssize_t:
    """Merge the int32/float32 pileups [t_p, t_v] (lt entries) and
    [c_p, c_v] (lc entries), each array given by its data pointer and
    its stride in bytes, into ret_p/ret_t/ret_c, and return the number
    of entries written.

    Each step writes the same values as the loop in
    ``CallerFromAlignments.__chrom_pair_treat_ctrl``: the smaller of the
    two positions, t_v[it] and c_v[ic]; then it advances whichever side
    had the smaller position, or both when they are equal. The sides
    advance by the comparisons' results rather than by branches: on
    real pileups the side that advances changes about one step in
    three, which a branch mispredicts. Inlined at a call with constant
    strides, the stride multiplications compile away. The caller
    guarantees that each ret array holds lt + lc entries.
    """
    ir: cython.Py_ssize_t = 0
    it: cython.Py_ssize_t = 0
    ic: cython.Py_ssize_t = 0
    tp: cython.int
    cp: cython.int

    while it < lt and ic < lc:
        tp = cython.cast(cython.pointer(cython.int), t_p + it * t_ps)[0]
        cp = cython.cast(cython.pointer(cython.int), c_p + ic * c_ps)[0]
        ret_p[ir] = tp if tp <= cp else cp
        ret_t[ir] = cython.cast(cython.pointer(cython.float),
                                t_v + it * t_vs)[0]
        ret_c[ir] = cython.cast(cython.pointer(cython.float),
                                c_v + ic * c_vs)[0]
        ir += 1
        it += tp <= cp
        ic += cp <= tp
    return ir


@cython.cfunc
@cython.inline
def chi2_k1_cdf(x: cython.float) -> cython.float:
    """Return CDF of chi-square(1) evaluated at ``x``."""
    return erf(sqrt(x/2))


@cython.cfunc
@cython.inline
def log10_chi2_k1_cdf(x: cython.float) -> cython.float:
    """Return log10 CDF of chi-square(1) evaluated at ``x``."""
    return log10(erf(sqrt(x/2)))


@cython.cfunc
@cython.inline
def chi2_k2_cdf(x: cython.float) -> cython.float:
    """Return CDF of chi-square(2) evaluated at ``x``."""
    return 1 - exp(-x/2)


@cython.cfunc
@cython.inline
def log10_chi2_k2_cdf(x: cython.float) -> cython.float:
    """Return log10 CDF of chi-square(2) evaluated at ``x``."""
    return log1p(- exp(-x/2)) * LOG10_E


@cython.cfunc
@cython.inline
def chi2_k4_cdf(x: cython.float) -> cython.float:
    """Return CDF of chi-square(4) evaluated at ``x``."""
    return 1 - exp(-x/2) * (1 + x/2)


@cython.cfunc
@cython.inline
def log10_chi2_k4_CDF(x: cython.float) -> cython.float:
    """Return log10 CDF of chi-square(4) evaluated at ``x``."""
    return log1p(- exp(-x/2) * (1 + x/2)) * LOG10_E


@cython.cfunc
@cython.inline
def apply_multiple_cutoffs(multiple_score_arrays: list,
                           multiple_cutoffs: list) -> cnp.ndarray:
    """Count how many scores exceed their corresponding cutoffs."""
    i: cython.int
    ret: cnp.ndarray

    ret = multiple_score_arrays[0] > multiple_cutoffs[0]

    for i in range(1, len(multiple_score_arrays)):
        ret += multiple_score_arrays[i] > multiple_cutoffs[i]

    return ret


@cython.cfunc
@cython.inline
def get_from_multiple_scores(multiple_score_arrays: list,
                             index: cython.int) -> list:
    """Return scores at ``index`` from each array in ``multiple_score_arrays``."""
    ret: list = []
    i: cython.int

    for i in range(len(multiple_score_arrays)):
        ret.append(multiple_score_arrays[i][index])
    return ret


@cython.cfunc
@cython.inline
def get_logFE(x: cython.float,
              y: cython.float) -> cython.float:
    """Return ``log10`` fold enrichment for ``x`` over ``y``."""
    return log10(x/y)


@cython.cfunc
@cython.inline
def get_subtraction(x: cython.float,
                    y: cython.float) -> cython.float:
    """Return the difference ``x - y``."""
    return x - y


@cython.cfunc
@cython.inline
def getitem_then_subtract(peakset: list,
                          start: cython.int) -> list:
    """Return peak starts relative to ``start`` for axillary metrics."""
    a: list

    a = [x["start"] for x in peakset]
    for i in range(len(a)):
        a[i] = a[i] - start
    return a


@cython.cfunc
@cython.inline
def left_sum(data, pos: cython.int,
             width: cython.int) -> cython.int:
    """Return cumulative sum within ``width`` bases to the left of ``pos``."""
    return sum([data[x] for x in data if x <= pos and x >= pos - width])


@cython.cfunc
@cython.inline
def right_sum(data,
              pos: cython.int,
              width: cython.int) -> cython.int:
    """Return cumulative sum within ``width`` bases to the right of ``pos``."""
    return sum([data[x] for x in data if x >= pos and x <= pos + width])


@cython.cfunc
@cython.inline
def left_forward(data,
                 pos: cython.int,
                 window_size: cython.int) -> cython.int:
    """Return incremental count entering ``pos`` from the left window."""
    return data.get(pos, 0) - data.get(pos-window_size, 0)


@cython.cfunc
@cython.inline
def right_forward(data,
                  pos: cython.int,
                  window_size: cython.int) -> cython.int:
    """Return incremental count exiting ``pos`` toward the right window."""
    return data.get(pos + window_size, 0) - data.get(pos, 0)


@cython.cfunc
def median_from_value_length(value: cnp.ndarray,
                             length: list) -> cython.float:
    """
    """
    tmp: list
    c: cython.int
    tmp_l: cython.int
    tmp_v: cython.float
    mid_l: cython.float

    c = 0
    tmp = sorted(list(zip(value, length)))
    mid_l = sum(length)/2
    for (tmp_v, tmp_l) in tmp:
        c += tmp_l
        if c > mid_l:
            return tmp_v


@cython.cfunc
def mean_from_value_length(value: cnp.ndarray,
                           length: list) -> cython.float:
    """take of: list values and of: list corresponding lengths,
    calculate the mean.  An important function for bedGraph type of
    data.

    """
    i: cython.int
    tmp_l: cython.int
    ln: cython.int
    tmp_v: cython.double
    sum_v: cython.double
    tmp_sum: cython.double
    ret: cython.float

    sum_v = 0
    ln = 0

    for i in range(len(length)):
        tmp_l = length[i]
        tmp_v = cython.cast(cython.double, value[i])
        tmp_sum = tmp_v * tmp_l
        sum_v = tmp_sum + sum_v
        ln += tmp_l

    ret = cython.cast(cython.float, (sum_v/ln))

    return ret


@cython.cfunc
def find_optimal_cutoff(x: list, y: list) -> tuple:
    """Return the best cutoff x and y.

    We assume that total peak length increase exponentially while
    decreasing cutoff value. But while cutoff decreases to a point
    that background noises are captured, total length increases much
    faster. So we fit a linear model by taking the first 10 points,
    then look for the largest cutoff that


    *Currently, it is coded as a useless function.
    """
    npx: cnp.ndarray
    npy: cnp.ndarray
    npA: cnp.ndarray
    ln: cython.long
    i: cython.long
    m: cython.float
    c: cython.float             # slop and intercept
    sst: cython.float           # sum of squared total
    sse: cython.float           # sum of squared error
    rsq: cython.float           # R-squared

    ln = len(x)
    assert ln == len(y)
    npx = np.array(x)
    npy = np.log10(np.array(y))
    npA = np.vstack([npx, np.ones(len(npx))]).T

    for i in range(10, ln):
        # at least the largest 10 points
        m, c = np.linalg.lstsq(npA[:i], npy[:i], rcond=None)[0]
        sst = sum((npy[:i] - np.mean(npy[:i])) ** 2)
        sse = sum((npy[:i] - m*npx[:i] - c) ** 2)
        rsq = 1 - sse/sst
        # print i, x[i], y[i], m, c, rsq
    return (1.0, 1.0)


# ------------------------------------
# Worker processes for the p-score table
# ------------------------------------

# the most worker processes one call uses
MAX_WORKERS: cython.int = 8
# the caller whose p-score table is being built, read by forked workers
_caller_for_workers = None
# the threads shutting down pools whose results are all in
_pool_shutdowns = []


def shut_down_later(pool):
    """Shut ``pool`` down, waiting for its workers to exit, on a thread
    of its own, so that this process goes on meanwhile; a worker takes
    milliseconds to exit. join_pool_shutdowns waits for the thread."""
    try:
        thread = threading.Thread(target=pool.shutdown)
        thread.start()
    except RuntimeError:
        # no thread to be had: wait here
        pool.shutdown()
        return
    _pool_shutdowns.append(thread)


def join_pool_shutdowns():
    """Wait until every pool given to shut_down_later is shut down.

    Each pool step calls this before it forks, so that no thread of a
    finished pool runs when this process forks, and the p-score table
    step calls it when it is done. The call-peaks step does not: its
    workers exit while this process writes the peaks (the threads that
    shut a pool down run only when this process releases the GIL, which
    adding the peaks does not do), and the interpreter waits for the
    thread, not a daemon, at its exit."""
    global _pool_shutdowns
    threads = _pool_shutdowns
    _pool_shutdowns = []
    for thread in threads:
        thread.join()


def n_available_workers() -> int:
    """Return min(MAX_WORKERS, cores this process may run on), or 1
    when worker processes cannot be forked: on a platform without
    fork, or inside a daemonic process, which may not have children."""
    if "fork" not in multiprocessing.get_all_start_methods():
        return 1
    if multiprocessing.current_process().daemon:
        return 1
    if hasattr(os, "sched_getaffinity"):
        return min(MAX_WORKERS, len(os.sched_getaffinity(0)))
    return min(MAX_WORKERS, os.cpu_count() or 1)


def n_locations(track, chrom: bytes) -> int:
    """Return the number of reads or fragments of a track on a chromosome."""
    locs = track.get_locations_by_chr(chrom)
    if isinstance(locs, (list, tuple)):
        return sum([x.shape[0] for x in locs])
    return locs.shape[0]


# header of a pileup file: the dtype and length of each of the three
# arrays [pos, treat, ctrl], which follow it as raw bytes
PILEUP_FILE_HEADER = struct.Struct("<8sq8sq8sq")


def write_pileup_file(fd: int, arrays: list):
    """Write a chromosome's [pos, treat, ctrl] arrays to the open file
    descriptor ``fd`` and close it."""
    header: list = []

    with open(fd, "wb") as f:
        for a in arrays:
            if a.ndim != 1 or a.dtype.hasobject:
                raise TypeError("pileup arrays must be 1-D and numeric")
            header += [a.dtype.str.encode(), a.shape[0]]
        f.write(PILEUP_FILE_HEADER.pack(*header))
        for a in arrays:
            f.write(a)
    return


def read_pileup_file(filename) -> list:
    """Read back the [pos, treat, ctrl] arrays written by
    write_pileup_file."""
    arrays: list = []
    i: cython.int

    with open(filename, "rb") as f:
        header = PILEUP_FILE_HEADER.unpack(f.read(PILEUP_FILE_HEADER.size))
        for i in range(0, 6, 2):
            a = np.empty(header[i + 1], dtype=header[i].rstrip(b"\0").decode())
            if f.readinto(a) != a.nbytes:
                raise EOFError("pileup file %s is truncated" % filename)
            arrays.append(a)
    return arrays


@cython.final
@cython.cclass
class SlotQScores:
    """The q-score of the p-score in each slot of ``pscore_cache``,
    looked up in ``qtable`` the first time this process meets the
    slot, so that the q-score of a position takes one hash lookup
    instead of two.

    An entry holds the complemented bits of the q-score, or 0 for a
    slot not looked up yet; a q-score whose bits are 0xFFFFFFFF (a
    NaN) is looked up again each time. The entries are for the slots
    of a cache of capacity ``mask + 1``; a cache that grows has moved
    its slots, and needs new entries. A slot's p-score never changes
    once it is used, and a QScoreTable is filled before any q-score is
    read from it.
    """
    entries: cython.pointer(cython.uint)
    mask: cython.ulonglong
    qtable: QScoreTable

    def __cinit__(self, qtable: QScoreTable, mask: cython.ulonglong):
        self.entries = cython.cast(cython.pointer(cython.uint),
                                   calloc(mask + 1, cython.sizeof(cython.uint)))
        if self.entries == cython.NULL:
            raise MemoryError()
        self.mask = mask
        self.qtable = qtable

    def __dealloc__(self):
        free(self.entries)


slot_qscores = cython.declare(SlotQScores, None)


@cython.cfunc
def slot_qscores_for(qtable: QScoreTable) -> SlotQScores:
    """The SlotQScores of ``qtable`` for the slots ``pscore_cache``
    has now, kept from one call to the next while both stay the
    same."""
    global slot_qscores

    if (slot_qscores is None or slot_qscores.qtable is not qtable or
            slot_qscores.mask != pscore_cache.mask):
        slot_qscores = SlotQScores(qtable, pscore_cache.mask)
    return slot_qscores


@cython.final
@cython.cclass
class SlotAboveCutoff:
    """Whether the q-score of the p-score in each slot of
    ``pscore_cache`` is above ``threshold``, looked up in ``qtable``
    the first time this process meets the slot: 0 for a slot not
    looked up yet, 1 for a q-score at or below the threshold (or a
    NaN), 2 for one above it. Like SlotQScores, the flags are for the
    slots of a cache of capacity ``mask + 1``.
    """
    flags: cython.pointer(cython.uchar)
    mask: cython.ulonglong
    qtable: QScoreTable
    threshold: cython.float

    def __cinit__(self, qtable: QScoreTable, mask: cython.ulonglong,
                  threshold: cython.float):
        self.flags = cython.cast(cython.pointer(cython.uchar),
                                 calloc(mask + 1, 1))
        if self.flags == cython.NULL:
            raise MemoryError()
        self.mask = mask
        self.qtable = qtable
        self.threshold = threshold

    def __dealloc__(self):
        free(self.flags)


slot_above = cython.declare(SlotAboveCutoff, None)


@cython.cfunc
def slot_above_for(qtable: QScoreTable,
                   threshold: cython.float) -> SlotAboveCutoff:
    """The SlotAboveCutoff of ``qtable`` and ``threshold`` for the
    slots ``pscore_cache`` has now, kept from one call to the next
    while all three stay the same."""
    global slot_above

    if (slot_above is None or slot_above.qtable is not qtable or
            slot_above.mask != pscore_cache.mask or
            slot_above.threshold != threshold):
        slot_above = SlotAboveCutoff(qtable, pscore_cache.mask, threshold)
    return slot_above


def float32_threshold(cutoff):
    """The float32 t for which numpy's ``a > cutoff``, ``a`` a float32
    array, is ``a > t`` in C float arithmetic for every element, or
    None when there is none to be found here (a NaN or infinite
    cutoff, or one beyond float32's range).

    Whether numpy compares in float32 or in float64 depends on the
    version and on the type of ``cutoff``, so numpy itself is asked:
    the comparison is monotone in the element, and its threshold lies
    between t and the next float32 up exactly when numpy finds t not
    above ``cutoff`` and that next float32 above it.
    """
    if not (-3.0e38 < cutoff < 3.0e38):
        return None
    t = np.float32(cutoff)
    if (np.array([t], dtype="f4") > cutoff)[0]:
        t = np.nextafter(t, np.float32(-np.inf))
    check = np.array([t, np.nextafter(t, np.float32(np.inf))],
                     dtype="f4") > cutoff
    if check[0] or not check[1]:
        return None
    return float(t)


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def qscore_above_indices(qtable: QScoreTable,
                         array1: cnp.ndarray,
                         array2: cnp.ndarray,
                         threshold: cython.float) -> cnp.ndarray:
    """The indices, in increasing order and as int32, of the positions
    whose q-score is above ``threshold``: what ``mask_indices`` returns
    for the q-score array of ``array1`` (treatment) and ``array2``
    (control) compared with a cutoff that ``float32_threshold`` maps
    to ``threshold``, without the q-score array or the mask.

    The positions are looked up in the order and with the calls
    __cal_qscore makes, so a p-score missing from ``qtable`` raises
    the same KeyError at the same position.
    """
    a1_ptr: cython.pointer(cython.float)
    a2_ptr: cython.pointer(cython.float)
    o_ptr: cython.pointer(cython.int)
    n: cython.Py_ssize_t
    i: cython.Py_ssize_t
    count: cython.Py_ssize_t
    slot: cython.Py_ssize_t
    f: cython.uchar
    q: cython.float
    sa: SlotAboveCutoff
    out: cnp.ndarray

    assert array1.shape[0] == array2.shape[0]
    n = array1.shape[0]
    # pages past the last index written are never touched
    out = np.empty(n, dtype="i4")
    a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
    a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
    o_ptr = cython.cast(cython.pointer(cython.int), out.data)

    sa = slot_above_for(qtable, threshold)
    count = 0
    for i in range(n):
        slot = pscore_cache.index(cython.cast(cython.int, a1_ptr[i]),
                                  a2_ptr[i])
        if pscore_cache.mask != sa.mask:
            # the cache grew on a miss, which moved its slots
            sa = slot_above_for(qtable, threshold)
        f = sa.flags[slot]
        if f == 0:
            q = qtable.get(pscore_cache.slots[slot].pscore)
            f = 2 if q > threshold else 1
            sa.flags[slot] = f
        if f == 2:
            o_ptr[count] = cython.cast(cython.int, i)
            count += 1
    return out[:count]


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def mask_indices(mask: cnp.ndarray) -> cnp.ndarray:
    """The indices of the nonzero entries of a 1-D array, in
    increasing order, as int32, which is what
    ``np.arange(n, dtype="i4")[np.nonzero(mask)[0]]`` returns."""
    m_ptr: cython.pointer(cython.char)
    o_ptr: cython.pointer(cython.int)
    n: cython.Py_ssize_t
    j: cython.Py_ssize_t
    count: cython.Py_ssize_t
    out: cnp.ndarray

    if mask.dtype != np.bool_ or not mask.flags.c_contiguous:
        return np.arange(mask.shape[0], dtype="i4")[np.nonzero(mask)[0]]
    n = mask.shape[0]
    m_ptr = cython.cast(cython.pointer(cython.char), mask.data)
    count = 0
    for j in range(n):
        count += m_ptr[j] != 0
    out = np.empty(count, dtype="i4")
    o_ptr = cython.cast(cython.pointer(cython.int), out.data)
    count = 0
    for j in range(n):
        if m_ptr[j] != 0:
            o_ptr[count] = cython.cast(cython.int, j)
            count += 1
    return out


# CallerFromAlignments.destroy keeps a descriptor on at most this many
# unlinked pileup files, those of at least this many bytes, for a
# thread to close
CLOSE_LATER_MAX_FILES = cython.declare(cython.Py_ssize_t, 256)
CLOSE_LATER_MIN_BYTES = cython.declare(cython.longlong, 1 << 20)


def close_unlinked_files(fds: list):
    """Close ``fds``, the last references to pileup files that are
    already unlinked; the kernel frees each file's pages in this
    close. The closes run without the GIL, so that a thread running
    this one waits for the GIL once, not once per file."""
    n: cython.Py_ssize_t = len(fds)
    i: cython.Py_ssize_t
    c_fds: cython.pointer(cython.int)

    c_fds = cython.cast(cython.pointer(cython.int),
                        calloc(n + 1, cython.sizeof(cython.int)))
    if c_fds == cython.NULL:
        raise MemoryError()
    for i in range(n):
        c_fds[i] = fds[i]
    with cython.nogil:
        for i in range(n):
            close(c_fds[i])
    free(c_fds)
    return


def close_unlinked_files_later(fds: list, sizes: list):
    """Close ``fds``, descriptors on unlinked files of ``sizes`` bytes,
    on up to n_available_workers() threads, so that the kernel frees
    the files' pages, tens of milliseconds per GB on one thread, while
    the caller goes on. Each thread gets about the same number of
    bytes: the largest file goes to the thread with the fewest bytes
    so far. The threads are not daemons, so the interpreter waits for
    them before it exits. A thread that cannot be started closes its
    files here."""
    n_threads: cython.Py_ssize_t = min(n_available_workers(), len(fds))
    groups: list = [[] for _ in range(n_threads)]
    loads: list = [0] * n_threads
    j: cython.Py_ssize_t
    k: cython.Py_ssize_t

    for j in sorted(range(len(fds)), key=sizes.__getitem__, reverse=True):
        k = loads.index(min(loads))
        groups[k].append(fds[j])
        loads[k] += sizes[j]
    for k in range(n_threads):
        try:
            threading.Thread(target=close_unlinked_files,
                             args=(groups[k],)).start()
        except RuntimeError:
            # no thread to be had: free the pages here
            close_unlinked_files(groups[k])
    return


# how many puts ahead the q-score loops prefetch a slot of the table
QTABLE_PREFETCH_AHEAD = cython.declare(cython.Py_ssize_t, 16)


@cython.cfunc
@cython.inline
@cython.profile(False)
@cython.linetrace(False)
@cython.exceptval(check=False)
def prefetch_qscore_slot(table: QScoreTable,
                         pscore: cython.float) -> cython.void:
    """Start loading the slot at which ``table.put(pscore, ...)`` will
    begin its probe; a hint only, which changes nothing in the table."""
    macs3_prefetch_w(cython.address(
        table.slots[pscore_hash(qscore_key(pscore)) & table.mask]))


# ------------------------------------
# Classes
# ------------------------------------
@cython.cclass
class CallerFromAlignments:
    """Compute pileups, scores, and peaks from FWTrack/PETrack alignments.

    Attributes:
        treat: Treatment track (FWTrack or PETrackI/PETrackII).
        ctrl: Control track (FWTrack or PETrackI/PETrackII).
        d: Fragment extension size for treatment (unused in PE mode).
        ctrl_d_s: Fragment extension size for control.
        treat_scaling_factor: Scaling factor for treatment.
        ctrl_scaling_factor_s: Scaling factors for control.
        lambda_bg: Minimum local bias to fill missing values.
        chromosomes: Common chromosome names in treat/control.
        pseudocount: Pseudocount used for logLR/FE/logFE calculations.
        bedGraph_filename_prefix: Prefix for pileup/lambda bedGraph outputs.
        end_shift: Shift applied to read ends before extension.
        trackline: Whether to emit UCSC track lines in bedGraph outputs.
        save_bedGraph: Whether to save pileup/lambda bedGraph files.
        save_SPMR: Whether to save pileup normalized per million reads.
        no_lambda_flag: Whether to ignore local lambda (use global only).
        PE_mode: Whether treatment is paired-end.
    """
    treat: object            # FWTrack or PETrackI/II object for ChIP
    ctrl: object             # FWTrack or PETrackI/II object for Control

    d: cython.int                           # extension size for ChIP
    # extension sizes for Control. Can be multiple values
    ctrl_d_s: list
    treat_scaling_factor: cython.float       # scaling factor for ChIP
    # scaling factor for Control, corresponding to each extension size.
    ctrl_scaling_factor_s: list
    # minimum local bias to fill missing values
    lambda_bg: cython.float
    # name of common chromosomes in ChIP and Control data
    chromosomes: list
    # the pseudocount used to calcuate logLR, FE or logFE
    pseudocount: cython.double
    # prefix will be added to _pileup.bdg for treatment and
    # _lambda.bdg for control
    bedGraph_filename_prefix: bytes
    # shift of cutting ends before extension
    end_shift: cython.int
    # whether trackline should be saved in bedGraph
    trackline: bool
    # whether to save pileup and local bias in bedGraph files
    save_bedGraph: bool
    # whether to save pileup normalized by sequencing depth in million reads
    save_SPMR: bool
    # whether ignore local bias, and to use global bias instead
    no_lambda_flag: bool
    # whether it's in PE mode, will be detected during initiation
    PE_mode: bool

    # temporary data buffer
    # temporary [position, treat_pileup, ctrl_pileup] for a given chromosome
    chr_pos_treat_ctrl: list
    bedGraph_treat_filename: bytes
    bedGraph_control_filename: bytes
    bedGraph_treat_f: cython.pointer(FILE)
    bedGraph_ctrl_f: cython.pointer(FILE)

    # data needed to be pre-computed before peak calling
    # remember pvalue->qvalue convertion, for lookups from C; empty
    # until the table is computed
    qtable: QScoreTable
    # the same mapping in a dict, filled only by __pre_computes
    # (--cutoff-analysis), which reads it back
    pqtable: dict
    # whether the pvalue of whole genome is all calculated. If yes,
    # it's OK to calculate q-value.
    pvalue_all_done: bool
    # record for each pvalue cutoff, how many peaks can be called
    pvalue_npeaks: dict
    # record for each pvalue cutoff, the total length of called peaks
    pvalue_length: dict
    # automatically decide the p-value cutoff (can be translated into
    # qvalue cutoff) based on p-value to total peak length analysis.
    optimal_p_cutoff: cython.float
    # file to save the pvalue-npeaks-totallength table
    cutoff_analysis_filename: bytes
    # Record the names of temporary files for storing pileup values of
    # each chromosome
    pileup_data_files: dict

    def __init__(self,
                 treat,
                 ctrl,
                 d: cython.int = 200,
                 ctrl_d_s: list = [200, 1000, 10000],
                 treat_scaling_factor: cython.float = 1.0,
                 ctrl_scaling_factor_s: list = [1.0, 0.2, 0.02],
                 stderr_on: bool = False,
                 pseudocount: cython.float = 1,
                 end_shift: cython.int = 0,
                 lambda_bg: cython.float = 0,
                 save_bedGraph: bool = False,
                 bedGraph_filename_prefix: str = "PREFIX",
                 bedGraph_treat_filename: str = "TREAT.bdg",
                 bedGraph_control_filename: str = "CTRL.bdg",
                 cutoff_analysis_filename: str = "TMP.txt",
                 save_SPMR: bool = False):
        """Initialize.

        A calculator is unique to each comparison of treat and
        control. Treat_depth and ctrl_depth should not be changed
        during calculation.

        treat and ctrl are either FWTrack or PETrackI objects.

        treat_depth and ctrl_depth are effective depth in million:
                                    sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        d, sregion, lregion: d is the fragment size, sregion is the
                             small region size, lregion is the large
                             region size

        pseudocount: a pseudocount used to calculate logLR, FE or
                     logFE. Please note this value will not be changed
                     with normalization method. So if you really want
                     to pseudocount: set 1 per million reads, it: set
                     after you normalize treat and control by million
                     reads by `change_normalizetion_method(ord('M'))`.
        
        Examples:
        .. code-block:: python

            from MACS3.Signal.CallPeakUnit import CallerFromAlignments
            caller = CallerFromAlignments(treat, ctrl, d=200)
        """
        chr1: set
        chr2: set
        p: cython.float

        # decide PE mode
        if isinstance(treat, FWTrack):
            self.PE_mode = False
        elif isinstance(treat, PETrackI):
            self.PE_mode = True
        elif isinstance(treat, PETrackII):
            self.PE_mode = True
        else:
            raise Exception("Should be FWTrack or PETrackI/II object!")
        # decide if there is control
        self.treat = treat
        if ctrl:
            self.ctrl = ctrl
        else:                   # while there is no control
            self.ctrl = treat
        self.trackline = False
        self.d = d              # note, self.d doesn't make sense in PE mode
        self.ctrl_d_s = ctrl_d_s  # note, self.d doesn't make sense in PE mode
        self.treat_scaling_factor = treat_scaling_factor
        self.ctrl_scaling_factor_s = ctrl_scaling_factor_s
        self.end_shift = end_shift
        self.lambda_bg = lambda_bg
        self.pqtable = {}
        self.qtable = QScoreTable()
        self.save_bedGraph = save_bedGraph
        self.save_SPMR = save_SPMR
        self.bedGraph_filename_prefix = bedGraph_filename_prefix.encode()
        self.bedGraph_treat_filename = bedGraph_treat_filename.encode()
        self.bedGraph_control_filename = bedGraph_control_filename.encode()
        if not self.ctrl_d_s or not self.ctrl_scaling_factor_s:
            self.no_lambda_flag = True
        else:
            self.no_lambda_flag = False
        self.pseudocount = pseudocount
        # get the common chromosome names from both treatment and control
        chr1 = set(self.treat.get_chr_names())
        chr2 = set(self.ctrl.get_chr_names())
        self.chromosomes = sorted(list(chr1.intersection(chr2)))

        self.pileup_data_files = {}
        self.pvalue_length = {}
        self.pvalue_npeaks = {}
        # step for optimal cutoff is 0.3 in -log10pvalue, we try from
        # pvalue 1E-10 (-10logp=10) to 0.5 (-10logp=0.3)
        for p in np.arange(0.3, 10, 0.3):
            self.pvalue_length[p] = 0
            self.pvalue_npeaks[p] = 0
        self.optimal_p_cutoff = 0
        self.cutoff_analysis_filename = cutoff_analysis_filename.encode()

    @cython.ccall
    def destroy(self):
        """Remove temporary pileup files created during peak calling.

        Each file is unlinked here, so the temporary directory is left
        as before. Freeing a file's pages (2.75 GB for a large control)
        is the slow part, and the kernel does it when the last
        reference to the file goes: a descriptor opened before the
        unlink holds that reference, and threads close the
        descriptors while the caller goes on to write the peaks. The
        threads are not daemons, so the interpreter waits for them
        before it exits.
        """
        f: bytes
        fds: list = []
        sizes: list = []

        try:
            for f in self.pileup_data_files.values():
                # os.path.isfile(f), keeping the size
                try:
                    st = os.stat(f)
                except (OSError, ValueError):
                    continue
                if not S_ISREG(st.st_mode):
                    continue
                # a bounded number of descriptors, on files large
                # enough to be worth it
                if (st.st_size >= CLOSE_LATER_MIN_BYTES and
                        len(fds) < CLOSE_LATER_MAX_FILES):
                    try:
                        fds.append(os.open(f, os.O_RDONLY | os.O_CLOEXEC))
                        sizes.append(st.st_size)
                    except OSError:
                        pass
                os.unlink(f)
        except BaseException:
            close_unlinked_files(fds)
            raise
        if fds:
            close_unlinked_files_later(fds, sizes)
        return

    @cython.ccall
    def set_pseudocount(self, pseudocount: cython.float):
        """Update the pseudocount used in scoring."""
        self.pseudocount = pseudocount

    @cython.ccall
    def enable_trackline(self):
        """Enable UCSC track line output when writing bedGraphs."""
        self.trackline = True

    @cython.cfunc
    def pileup_treat_ctrl_a_chromosome(self, chrom: bytes):
        """After this function is called, self.chr_pos_treat_ctrl will
        be reand: set assigned to the pileup values of the given
        chromosome.

        """
        treat_pv: list
        ctrl_pv: list
        temp_filename: str

        assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom

        # check backup file of pileup values. If not exists, create
        # it. Otherwise, load them instead of calculating new pileup
        # values.
        if chrom in self.pileup_data_files:
            try:
                self.chr_pos_treat_ctrl = read_pileup_file(self.pileup_data_files[chrom])
                return
            except Exception:
                pass

        # reor: set clean existing self.chr_pos_treat_ctrl
        if self.chr_pos_treat_ctrl:     # not a beautiful way to clean
            clean_up_ndarray(self.chr_pos_treat_ctrl[0])
            clean_up_ndarray(self.chr_pos_treat_ctrl[1])
            clean_up_ndarray(self.chr_pos_treat_ctrl[2])

        if self.PE_mode:
            treat_pv = self.treat.pileup_a_chromosome(chrom,
                                                      self.treat_scaling_factor,
                                                      baseline_value=0.0)
        else:
            treat_pv = self.treat.pileup_a_chromosome(chrom,
                                                      self.d,
                                                      self.treat_scaling_factor,
                                                      baseline_value=0.0,
                                                      directional=True,
                                                      end_shift=self.end_shift)

        if not self.no_lambda_flag:
            if self.PE_mode:
                # note, we pileup up PE control as SE control because
                # we assume the bias only can be captured at the
                # surrounding regions of cutting sites from control experiments.
                ctrl_pv = self.ctrl.pileup_a_chromosome_c(chrom,
                                                          self.ctrl_d_s,
                                                          self.ctrl_scaling_factor_s,
                                                          baseline_value=self.lambda_bg)
            else:
                ctrl_pv = self.ctrl.pileup_a_chromosome_c(chrom,
                                                          self.ctrl_d_s,
                                                          self.ctrl_scaling_factor_s,
                                                          baseline_value=self.lambda_bg,
                                                          directional=False)
        else:
            # a: set global lambda
            ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg,],
                                                  dtype="f4")]

        self.chr_pos_treat_ctrl = self.__chrom_pair_treat_ctrl(treat_pv,
                                                               ctrl_pv)

        # clean treat_pv and ctrl_pv
        treat_pv = []
        ctrl_pv = []

        # save data to temporary file, through the descriptor mkstemp
        # opened: reopening the name with truncation ("wb") makes ext4
        # flush the file to disk when it is closed
        temp_fd, temp_filename = mkstemp()
        self.pileup_data_files[chrom] = temp_filename.encode()
        try:
            write_pileup_file(temp_fd, self.chr_pos_treat_ctrl)
        except Exception:
            # fail to write then remove the key in pileup_data_files
            self.pileup_data_files.pop(chrom)
        return

    @cython.cfunc
    def __chrom_pair_treat_ctrl(self, treat_pv, ctrl_pv) -> list:
        """*private* Pair treat and ctrl pileup for each region.

        treat_pv and ctrl_pv are [np.ndarray, np.ndarray].

        return [p, t, c] list, each element is a numpy array.
        """
        ir: cython.long         # index of ret_p/t/c
        it: cython.long         # index of t_p/v
        ic: cython.long         # index of c_p/v
        lt: cython.long
        lc: cython.long
        t_p: cnp.ndarray
        c_p: cnp.ndarray
        ret_p: cnp.ndarray
        t_v: cnp.ndarray
        c_v: cnp.ndarray
        ret_t: cnp.ndarray
        ret_c: cnp.ndarray
        t_ps: cython.Py_ssize_t
        t_vs: cython.Py_ssize_t
        c_ps: cython.Py_ssize_t
        c_vs: cython.Py_ssize_t

        [t_p, t_v] = treat_pv
        [c_p, c_v] = ctrl_pv

        lt = t_p.shape[0]
        lc = c_p.shape[0]

        chrom_max_len = lt + lc

        if (t_p.ndim == 1 and t_v.ndim == 1 and c_p.ndim == 1 and
                c_v.ndim == 1 and t_p.dtype == np.int32 and
                c_p.dtype == np.int32 and t_v.dtype == np.float32 and
                c_v.dtype == np.float32 and t_v.shape[0] >= lt and
                c_v.shape[0] >= lc):
            # every pileup function returns int32 positions and float32
            # values: merge them in C, into arrays it fills up to the
            # length they are cut to. Anything else takes the loop below.
            ret_p = np.empty(chrom_max_len, dtype="i4")
            ret_t = np.empty(chrom_max_len, dtype="f4")
            ret_c = np.empty(chrom_max_len, dtype="f4")
            t_ps = t_p.strides[0]
            t_vs = t_v.strides[0]
            c_ps = c_p.strides[0]
            c_vs = c_v.strides[0]
            if t_ps == 4 and t_vs == 4 and c_ps == 4 and c_vs == 4:
                # contiguous arrays
                ir = pair_pv_i4_f4(t_p.data, 4, t_v.data, 4, lt,
                                   c_p.data, 4, c_v.data, 4, lc,
                                   cython.cast(cython.pointer(cython.int),
                                               ret_p.data),
                                   cython.cast(cython.pointer(cython.float),
                                               ret_t.data),
                                   cython.cast(cython.pointer(cython.float),
                                               ret_c.data))
            else:
                # e.g. the fields of a structured array
                ir = pair_pv_i4_f4(t_p.data, t_ps, t_v.data, t_vs, lt,
                                   c_p.data, c_ps, c_v.data, c_vs, lc,
                                   cython.cast(cython.pointer(cython.int),
                                               ret_p.data),
                                   cython.cast(cython.pointer(cython.float),
                                               ret_t.data),
                                   cython.cast(cython.pointer(cython.float),
                                               ret_c.data))
            ret_p.resize(ir, refcheck=False)
            ret_t.resize(ir, refcheck=False)
            ret_c.resize(ir, refcheck=False)
            return [ret_p, ret_t, ret_c]

        ret_p = np.zeros(chrom_max_len, dtype="i4")  # position
        ret_t = np.zeros(chrom_max_len, dtype="f4")  # value from treatment
        ret_c = np.zeros(chrom_max_len, dtype="f4")  # value from control

        ir = 0
        it = 0
        ic = 0

        while it < lt and ic < lc:
            if t_p[it] < c_p[ic]:
                # clip a region from pre_p to p1, then pre_p: set as p1.
                ret_p[ir] = t_p[it]
                ret_t[ir] = t_v[it]
                ret_c[ir] = c_v[ic]
                ir += 1
                # call for the next p1 and v1
                it += 1
            elif t_p[it] > c_p[ic]:
                # clip a region from pre_p to p2, then pre_p: set as p2.
                ret_p[ir] = c_p[ic]
                ret_t[ir] = t_v[it]
                ret_c[ir] = c_v[ic]
                ir += 1
                # call for the next p2 and v2
                ic += 1
            else:
                # from pre_p to p1 or p2, then pre_p: set as p1 or p2.
                ret_p[ir] = t_p[it]
                ret_t[ir] = t_v[it]
                ret_c[ir] = c_v[ic]
                ir += 1
                # call for the next p1, v1, p2, v2.
                it += 1
                ic += 1

        ret_p.resize(ir, refcheck=False)
        ret_t.resize(ir, refcheck=False)
        ret_c.resize(ir, refcheck=False)
        return [ret_p, ret_t, ret_c]

    @cython.cfunc
    def __cal_score(self,
                    array1: cnp.ndarray(cython.float, ndim=1),
                    array2: cnp.ndarray(cython.float, ndim=1),
                    cal_func) -> cnp.ndarray:
        """Apply ``cal_func`` element-wise across two equally sized arrays."""
        i: cython.long
        s: cnp.ndarray(cython.float, ndim=1)

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")
        for i in range(array1.shape[0]):
            s[i] = cal_func(array1[i], array2[i])
        return s

    @cython.cfunc
    def _pscore_stat_a_chromosome(self, i: cython.long,
                                  score: cython.bint):
        """Pile up chromosome ``i``, which also saves its pileup to a
        temporary file, and add the length of each of its segments to
        the histogram held in the p-score cache, in the slot of the
        segment's (observed, lambda).

        With ``score``, the p-score of each new slot is computed;
        without, it is not, as in a worker process whose histogram
        another process scores.

        A slot's ``first`` is set by the first segment that uses it in
        the pass, so the chromosomes of one pass are added in
        increasing ``i``.
        """
        pos_array: cnp.ndarray
        treat_array: cnp.ndarray
        ctrl_array: cnp.ndarray
        pre_p: cython.long
        j: cython.long
        this_l: cython.long
        pos_view: cython.pointer(cython.int)
        treat_value_view: cython.pointer(cython.float)
        ctrl_value_view: cython.pointer(cython.float)
        cache: PScoreCache
        slot: cython.pointer(PScoreSlot)
        s: cython.Py_ssize_t
        n_seg: cython.longlong

        cache = pscore_cache
        pre_p = 0
        # segments are numbered in the serial order, by chromosome and
        # then by position; positions are 32-bit, so a chromosome has
        # fewer than 2**32 segments
        n_seg = cython.cast(cython.longlong, i) << 32

        self.pileup_treat_ctrl_a_chromosome(self.chromosomes[i])
        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        pos_view = cython.cast(cython.pointer(cython.int),
                               pos_array.data)
        treat_value_view = cython.cast(cython.pointer(cython.float),
                                       treat_array.data)
        ctrl_value_view = cython.cast(cython.pointer(cython.float),
                                      ctrl_array.data)

        # a new slot's p-score is computed after the loop, with the
        # chromosome's other new ones; the histogram needs only the slot
        for j in range(pos_array.shape[0]):
            if score:
                s = cache.index_deferred(cython.cast(cython.int,
                                                     treat_value_view[0]),
                                         ctrl_value_view[0])
            else:
                s = cache.index_unscored(cython.cast(cython.int,
                                                     treat_value_view[0]),
                                         ctrl_value_view[0])
            slot = cython.address(cache.slots[s])
            this_l = pos_view[0] - pre_p
            if slot.first < 0:
                slot.first = n_seg
                slot.length = this_l
            else:
                slot.length += this_l
            n_seg += 1
            pre_p = pos_view[0]
            pos_view += 1
            treat_value_view += 1
            ctrl_value_view += 1
        if score:
            cache.compute_pending()
        return

    @cython.cfunc
    def _chromosome_batches(self, n_workers: cython.int) -> list:
        """Group the chromosome indices into tasks for the worker
        processes, largest first.

        A chromosome's cost is taken as its number of reads or
        fragments in treatment and control. A chromosome at least
        1/(32 n_workers) of the total is a task of its own; smaller
        ones are grouped until a task reaches that size, so that 50,000
        contigs do not become 50,000 tasks.
        """
        i: cython.long
        w: cython.long
        total: cython.long
        target: cython.long
        batch_weight: cython.long
        order: list
        batches: list
        batch: list

        order = []
        total = 0
        for i in range(len(self.chromosomes)):
            w = n_locations(self.treat, self.chromosomes[i])
            if self.ctrl is not self.treat:
                w += n_locations(self.ctrl, self.chromosomes[i])
            order.append((-w, i))
            total += w
        order.sort()

        target = total // (32 * n_workers)
        batches = []
        batch = []
        batch_weight = 0
        for (w, i) in order:
            batch.append(i)
            batch_weight -= w
            if batch_weight >= target:
                batches.append(batch)
                batch = []
                batch_weight = 0
        if batch:
            batches.append(batch)
        return batches

    @cython.cfunc
    def _pscore_stat_in_workers(self, n_workers: cython.int) -> tuple:
        """Build the p-score histogram of all chromosomes in forked
        worker processes, which share the tracks copy-on-write.

        Each worker fills its own copy of the p-score cache with keys
        but no p-scores. For each task it returns the names of the
        files holding its chromosomes' pileups and, from
        PScoreCache.take_histogram(True), the key, first segment and
        summed length of every slot its segments used. As each task
        finishes, while the workers run the others, its histogram is
        merged into this process's cache: a slot keeps the first
        segment that used it in any task and the summed length.
        ``first`` numbers the segments in the serial loop's order
        across all chromosomes, so each slot's is the serial pass's.

        Once every task is in, the p-score of each distinct new key is
        computed once, on up to n_workers threads, and stored in
        the cache, so that the call-peaks step's workers, forked
        later, start with them. Returns the p-score, first segment and
        summed length of every slot used, as take_histogram does; they
        sort by p-score and then ``first`` exactly as the serial
        pass's do.
        """
        global _caller_for_workers
        batches: list
        files: dict
        all_new: cython.bint
        chrom: bytes

        batches = self._chromosome_batches(n_workers)
        files = {}
        # every worker's copy of the cache starts with no slot used
        pscore_cache.reset_histogram()
        # with no key held before the pass, every slot it uses is new
        all_new = pscore_cache.n == 0
        # the workers need not inherit the heap's free memory
        macs3_release_free_heap()
        join_pool_shutdowns()
        _caller_for_workers = self
        try:
            pool = ProcessPoolExecutor(max_workers=min(n_workers, len(batches)),
                                       mp_context=multiprocessing.get_context("fork"))
            try:
                futures = {pool.submit(pscore_stat_worker, batches[k]): k
                           for k in range(len(batches))}
                try:
                    for future in as_completed(futures):
                        (chrom_files, keys, firsts, lengths) = future.result()
                        for (chrom, temp_filename) in chrom_files:
                            files[chrom] = temp_filename
                        pscore_cache.merge_histogram(keys, firsts, lengths)
                except BaseException:
                    pool.shutdown(cancel_futures=True)
                    raise
            except BaseException:
                pool.shutdown()
                raise
        finally:
            _caller_for_workers = None
        # every result is in: the workers exit while this process
        # merges the histograms; __cal_pvalue_qvalue_table waits for
        # them when it is done
        shut_down_later(pool)

        for chrom in self.chromosomes:
            temp_filename = files[chrom]
            if temp_filename is not None:
                self.pileup_data_files[chrom] = temp_filename
        if all_new:
            return pscore_cache.score_pending(n_workers, True)
        pscore_cache.score_pending(n_workers, False)
        return pscore_cache.take_histogram()

    @cython.cfunc
    def __cal_pvalue_qvalue_table(self):
        """Populate ``self.qtable`` by scanning all chromosomes."""
        # pre_l: cython.long
        l: cython.long
        i: cython.long
        j: cython.long
        # pre_v: cython.float
        v: cython.float
        q: cython.float
        pre_q: cython.float
        N: cython.long
        k: cython.long
        f: cython.float
        n_workers: cython.int
        s: cython.Py_ssize_t
        n_used: cython.Py_ssize_t
        n_unique: cython.Py_ssize_t
        pscores: cnp.ndarray
        firsts: cnp.ndarray
        lengths: cnp.ndarray
        order: cnp.ndarray
        unique_values: cnp.ndarray
        unique_lengths: cnp.ndarray
        ps_ptr: cython.pointer(cython.float)
        ln_ptr: cython.pointer(cython.longlong)
        order_ptr: cython.pointer(cython.Py_ssize_t)
        uv_ptr: cython.pointer(cython.float)
        ul_ptr: cython.pointer(cython.longlong)
        qtable: QScoreTable

        debug("Start to calculate pvalue stat...")

        # pscore_stat, the total length of the segments at each
        # p-score, is first summed per (observed, lambda) slot of the
        # p-score cache, then merged by p-score.
        n_workers = n_available_workers()
        if n_workers > 1 and len(self.chromosomes) > 1:
            (pscores, firsts, lengths) = self._pscore_stat_in_workers(n_workers)
            # where the serial loop leaves i, which the last loop
            # below starts from when there are no p-scores at all
            i = len(self.chromosomes) - 1
        else:
            pscore_cache.reset_histogram()
            for i in range(len(self.chromosomes)):
                self._pscore_stat_a_chromosome(i, True)
            (pscores, firsts, lengths) = pscore_cache.take_histogram()

        # the slots used in this pass: p-score, first segment, length
        n_used = pscores.shape[0]
        ps_ptr = cython.cast(cython.pointer(cython.float), pscores.data)
        ln_ptr = cython.cast(cython.pointer(cython.longlong), lengths.data)

        # Distinct p-scores in descending order, with their summed
        # lengths, as a dict keyed by p-score would hold them: equal
        # p-scores merge, and the value kept is the one met first,
        # which matters only for 0.0 and -0.0.
        order = pscore_order(pscores, firsts)
        order_ptr = cython.cast(cython.pointer(cython.Py_ssize_t), order.data)
        unique_values = np.empty(n_used, dtype="f4")
        unique_lengths = np.empty(n_used, dtype="i8")
        uv_ptr = cython.cast(cython.pointer(cython.float), unique_values.data)
        ul_ptr = cython.cast(cython.pointer(cython.longlong),
                             unique_lengths.data)
        n_unique = 0
        for s in range(n_used):
            if n_unique > 0 and ps_ptr[order_ptr[s]] == uv_ptr[n_unique - 1]:
                ul_ptr[n_unique - 1] += ln_ptr[order_ptr[s]]
            else:
                uv_ptr[n_unique] = ps_ptr[order_ptr[s]]
                ul_ptr[n_unique] = ln_ptr[order_ptr[s]]
                n_unique += 1

        N = 0
        for s in range(n_unique):
            N += ul_ptr[s]     # total length
        k = 1                          # rank
        f = -log10(N)
        # pre_v = -2147483647
        # pre_l = 0
        pre_q = 2147483647      # save the previous q-value

        # Only the C table is filled: pqtable, a Python-level map, is
        # read only by --cutoff-analysis, whose __pre_computes fills
        # its own. Each put lands on a random slot of a table larger
        # than the cache, so the slot of a put some puts ahead is
        # prefetched.
        qtable = QScoreTable(n_unique)
        self.qtable = qtable
        for i in range(n_unique):
            if i + QTABLE_PREFETCH_AHEAD < n_unique:
                prefetch_qscore_slot(qtable, uv_ptr[i + QTABLE_PREFETCH_AHEAD])
            v = uv_ptr[i]
            l = ul_ptr[i]
            q = v + (log10(k) + f)
            if q > pre_q:
                q = pre_q
            if q <= 0:
                q = 0
                break
            # q = max(0,min(pre_q,q))           # make q-score monotonic
            qtable.put(v, q)
            pre_q = q
            k += l
        # bottom rank pscores all have qscores 0
        for j in range(i, n_unique):
            if j + QTABLE_PREFETCH_AHEAD < n_unique:
                prefetch_qscore_slot(qtable, uv_ptr[j + QTABLE_PREFETCH_AHEAD])
            v = uv_ptr[j]
            qtable.put(v, 0)
        # wait for the pool's workers, which exit while this step runs
        join_pool_shutdowns()
        return

    @cython.cfunc
    def __pre_computes(self,
                       max_gap: cython.int = 50,
                       min_length: cython.int = 200):
        """After this function is called, self.pqtable and
        self.pvalue_length is built. All chromosomes will be
        iterated. So it will take some time.

        """
        chrom: bytes
        pos_array: cnp.ndarray
        treat_array: cnp.ndarray
        ctrl_array: cnp.ndarray
        score_array: cnp.ndarray
        pscore_stat: dict
        n: cython.long
        pre_p: cython.long
        this_p: cython.long
        j: cython.long
        l: cython.long
        i: cython.long
        q: cython.float
        pre_q: cython.float
        this_v: cython.float
        v: cython.float
        cutoff: cython.float
        N: cython.long
        k: cython.long
        this_l: cython.long
        f: cython.float
        unique_values: list
        above_cutoff: cnp.ndarray
        above_cutoff_endpos: cnp.ndarray
        above_cutoff_startpos: cnp.ndarray
        peak_content: list
        peak_length: cython.long
        total_l: cython.long
        total_p: cython.long
        tmplist: list

        # above cutoff start position pointer
        acs_ptr: cython.pointer(cython.int)
        # above cutoff end position pointer
        ace_ptr: cython.pointer(cython.int)
        # position array pointer
        pos_array_ptr: cython.pointer(cython.int)
        # score array pointer
        score_array_ptr: cython.pointer(cython.float)

        debug("Start to calculate pvalue stat...")

        # tmpcontains: list a of: list log pvalue cutoffs from 0.3 to 10
        tmplist = [round(x, 5)
                   for x in sorted(list(np.arange(0.3, 10.0, 0.3)),
                                   reverse=True)]

        pscore_stat = {}      # dict()
        # print (list(pscore_stat.keys()))
        # print (list(self.pvalue_length.keys()))
        # print (list(self.pvalue_npeaks.keys()))
        for i in range(len(self.chromosomes)):
            chrom = self.chromosomes[i]
            self.pileup_treat_ctrl_a_chromosome(chrom)
            [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

            score_array = self.__cal_pscore(treat_array, ctrl_array)

            for n in range(len(tmplist)):
                cutoff = tmplist[n]
                total_l = 0           # total length in potential peak
                total_p = 0

                # get the regions with scores above cutoffs this is
                # not an optimized method. It would be better to store
                # score array in a 2-D ndarray?
                above_cutoff = np.nonzero(score_array > cutoff)[0]
                # end positions of regions where score is above cutoff
                above_cutoff_endpos = pos_array[above_cutoff]
                # start positions of regions where score is above cutoff
                above_cutoff_startpos = pos_array[above_cutoff-1]

                if above_cutoff_endpos.size == 0:
                    continue

                # first bit of region above cutoff
                acs_ptr = cython.cast(cython.pointer(cython.int),
                                      above_cutoff_startpos.data)
                ace_ptr = cython.cast(cython.pointer(cython.int),
                                      above_cutoff_endpos.data)

                peak_content = [(acs_ptr[0], ace_ptr[0]),]
                lastp = ace_ptr[0]
                acs_ptr += 1
                ace_ptr += 1

                for i in range(1, above_cutoff_startpos.size):
                    tl = acs_ptr[0] - lastp
                    if tl <= max_gap:
                        peak_content.append((acs_ptr[0], ace_ptr[0]))
                    else:
                        peak_length = peak_content[-1][1] - peak_content[0][0]
                        # if the peak is too small, reject it
                        if peak_length >= min_length:
                            total_l += peak_length
                            total_p += 1
                        peak_content = [(acs_ptr[0], ace_ptr[0]),]
                    lastp = ace_ptr[0]
                    acs_ptr += 1
                    ace_ptr += 1

                if peak_content:
                    peak_length = peak_content[-1][1] - peak_content[0][0]
                    # if the peak is too small, reject it
                    if peak_length >= min_length:
                        total_l += peak_length
                        total_p += 1
                self.pvalue_length[cutoff] = self.pvalue_length.get(cutoff, 0) + total_l
                self.pvalue_npeaks[cutoff] = self.pvalue_npeaks.get(cutoff, 0) + total_p

            pos_array_ptr = cython.cast(cython.pointer(cython.int),
                                        pos_array.data)
            score_array_ptr = cython.cast(cython.pointer(cython.float),
                                          score_array.data)

            pre_p = 0
            for i in range(pos_array.shape[0]):
                this_p = pos_array_ptr[0]
                this_l = this_p - pre_p
                this_v = score_array_ptr[0]
                if this_v in pscore_stat:
                    pscore_stat[this_v] += this_l
                else:
                    pscore_stat[this_v] = this_l
                pre_p = this_p  # pos_array[i]
                pos_array_ptr += 1
                score_array_ptr += 1

        # debug ("make pscore_stat cost %.5f seconds" % t)

        # add all pvalue cutoffs from cutoff-analysis part. So that we
        # can get the corresponding qvalues for them.
        for cutoff in tmplist:
            if cutoff not in pscore_stat:
                pscore_stat[cutoff] = 0

        N = sum(pscore_stat.values())  # total length
        k = 1                           # rank
        f = -log10(N)
        pre_q = 2147483647              # save the previous q-value

        self.pqtable = {}
        # sorted(unique_values,reverse=True)
        unique_values = sorted(list(pscore_stat.keys()), reverse=True)
        self.qtable = QScoreTable(len(unique_values))
        for i in range(len(unique_values)):
            v = unique_values[i]
            l = pscore_stat[v]
            q = v + (log10(k) + f)
            if q > pre_q:
                q = pre_q
            if q <= 0:
                q = 0
                break
            # q = max(0,min(pre_q,q))           # make q-score monotonic
            self.pqtable[v] = q
            self.qtable.put(v, q)
            pre_q = q
            k += l
        for j in range(i, len(unique_values)):
            v = unique_values[j]
            self.pqtable[v] = 0
            self.qtable.put(v, 0)

        # write pvalue and total length of predicted peaks
        # this is the output from cutoff-analysis
        fhd = open(self.cutoff_analysis_filename, "w")
        fhd.write("pscore\tqscore\tnpeaks\tlpeaks\tavelpeak\n")
        x = []
        y = []
        for cutoff in tmplist:
            if self.pvalue_npeaks[cutoff] > 0:
                fhd.write("%.2f\t%.2f\t%d\t%d\t%.2f\n" %
                          (cutoff, self.pqtable[cutoff],
                           self.pvalue_npeaks[cutoff],
                           self.pvalue_length[cutoff],
                           self.pvalue_length[cutoff]/self.pvalue_npeaks[cutoff]))
                x.append(cutoff)
                y.append(self.pvalue_length[cutoff])
        fhd.close()
        info("#3 Analysis of cutoff vs num of peaks or total length has been saved in %s" % self.cutoff_analysis_filename)
        # info("#3 Suggest a cutoff...")
        # optimal_cutoff, optimal_length = find_optimal_cutoff(x, y)
        # info("#3 -10log10pvalue cutoff %.2f will call approximately %.0f bps regions as significant regions" % (optimal_cutoff, optimal_length))
        # print (list(pqtable.keys()))
        # print (list(self.pvalue_length.keys()))
        # print (list(self.pvalue_npeaks.keys()))
        return

    @cython.ccall
    def call_peaks(self,
                   scoring_function_symbols: list,
                   score_cutoff_s: list,
                   min_length: cython.int = 200,
                   max_gap: cython.int = 50,
                   call_summits: bool = False,
                   cutoff_analysis: bool = False):
        """Call narrow peaks for all chromosomes.

        Args:
            scoring_function_symbols: Symbols for score functions.
                Use ``'p'`` (pscore), ``'q'`` (qscore), ``'f'`` (fold change),
                or ``'s'`` (subtraction). Example: ``['p', 'q']``.
            score_cutoff_s: Cutoff values corresponding to ``scoring_function_symbols``.
            min_length: Minimum peak length.
            max_gap: Maximum gap of insignificant regions within a peak.
            call_summits: Whether to call sub-peaks (summits).
            cutoff_analysis: Whether to compute cutoff-vs-peak metrics.

        Returns:
            PeakIO: Collection of called peaks.

        Examples:
            .. code-block:: python

                peaks = caller.call_peaks(['p'], [5.0], min_length=200)
        """
        chrom: bytes
        tmp_bytes: bytes
        n_workers: cython.int

        peaks = PeakIO()

        # prepare p-q table
        if self.qtable.n == 0:
            info("#3 Pre-compute pvalue-qvalue table...")
            if cutoff_analysis:
                info("#3 Cutoff vs peaks called will be analyzed!")
                self.__pre_computes(max_gap=max_gap, min_length=min_length)
            else:
                self.__cal_pvalue_qvalue_table()

        # prepare bedGraph file
        if self.save_bedGraph:
            self.bedGraph_treat_f = fopen(self.bedGraph_treat_filename, "w")
            self.bedGraph_ctrl_f = fopen(self.bedGraph_control_filename, "w")

            info("#3 In the peak calling step, the following will be performed simultaneously:")
            info("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... %s" %
                 self.bedGraph_filename_prefix.decode() + "_treat_pileup.bdg")
            info("#3   Write bedGraph files for control lambda (after scaling if necessary)... %s" %
                 self.bedGraph_filename_prefix.decode() + "_control_lambda.bdg")

            if self.save_SPMR:
                info("#3   --SPMR is requested, so pileup will be normalized by sequencing depth in million reads.")
            elif self.treat_scaling_factor == 1:
                info("#3   Pileup will be based on sequencing depth in treatment.")
            else:
                info("#3   Pileup will be based on sequencing depth in control.")

            if self.trackline:
                # this line is REQUIRED by the wiggle format for UCSC browser
                tmp_bytes = ("track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf(self.bedGraph_treat_f, tmp_bytes)
                tmp_bytes = ("track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf(self.bedGraph_ctrl_f, tmp_bytes)

        info("#3 Call peaks for each chromosome...")
        n_workers = n_available_workers()
        if n_workers > 1 and len(self.chromosomes) > 1 and not self.save_bedGraph:
            # the bedGraph files are written chromosome by chromosome
            # through one file handle, so -B keeps the serial loop
            self._call_peaks_in_workers(peaks,
                                        n_workers,
                                        scoring_function_symbols,
                                        score_cutoff_s,
                                        min_length,
                                        max_gap,
                                        call_summits)
        else:
            for chrom in self.chromosomes:
                # treat/control bedGraph will be saved if requested by user.
                self.__chrom_call_peak_using_certain_criteria(peaks,
                                                              chrom,
                                                              scoring_function_symbols,
                                                              score_cutoff_s,
                                                              min_length,
                                                              max_gap,
                                                              call_summits,
                                                              self.save_bedGraph)

        # close bedGraph file
        if self.save_bedGraph:
            fclose(self.bedGraph_treat_f)
            fclose(self.bedGraph_ctrl_f)
            self.save_bedGraph = False

        return peaks

    @cython.cfunc
    def _call_peaks_a_chromosome(self,
                                 chrom: bytes,
                                 scoring_function_symbols: list,
                                 score_cutoff_s: list,
                                 min_length: cython.int,
                                 max_gap: cython.int,
                                 call_summits: bool) -> list:
        """Call the peaks of one chromosome, without bedGraph output,
        and return them in the order they were added, each as the
        tuple (start, end, summit, peak_score, pileup, pscore,
        fold_change, qscore) of the values the PeakContent holds."""
        chrom_peaks: object
        p: object
        state: tuple
        out: list

        chrom_peaks = PeakIO()
        self.__chrom_call_peak_using_certain_criteria(chrom_peaks,
                                                      chrom,
                                                      scoring_function_symbols,
                                                      score_cutoff_s,
                                                      min_length,
                                                      max_gap,
                                                      call_summits,
                                                      False)
        out = []
        for p in chrom_peaks.get_data_from_chrom(chrom):
            # (chrom, start, end, length, summit, score, pileup,
            #  pscore, fc, qscore, name)
            state = p.__getstate__()
            out.append((state[1], state[2], state[4], state[5],
                        state[6], state[7], state[8], state[9]))
        return out

    @cython.cfunc
    def _call_peaks_in_workers(self,
                               peaks,
                               n_workers: cython.int,
                               scoring_function_symbols: list,
                               score_cutoff_s: list,
                               min_length: cython.int,
                               max_gap: cython.int,
                               call_summits: bool):
        """Call the peaks of all chromosomes in forked worker
        processes, which share the tracks and the p/q-value table
        copy-on-write, and add them to ``peaks`` in chromosome order,
        as the serial loop does.

        Each worker returns, for each of its chromosomes, the name of
        the file holding the pileup and the chromosome's peaks.
        """
        global _caller_for_workers
        batches: list
        tasks: list
        results: dict
        part: list
        chrom_peaks: list
        chrom: bytes
        t: tuple

        batches = self._chromosome_batches(n_workers)
        tasks = [(batch, scoring_function_symbols, score_cutoff_s,
                  min_length, max_gap, call_summits) for batch in batches]
        results = {}
        # the workers need not inherit the heap's free memory
        macs3_release_free_heap()
        join_pool_shutdowns()
        _caller_for_workers = self
        try:
            pool = ProcessPoolExecutor(max_workers=min(n_workers, len(batches)),
                                       mp_context=multiprocessing.get_context("fork"))
            try:
                for part in pool.map(call_peaks_worker, tasks):
                    for (chrom, temp_filename, chrom_peaks) in part:
                        results[chrom] = (temp_filename, chrom_peaks)
            except BaseException:
                pool.shutdown()
                raise
        finally:
            _caller_for_workers = None
        # every result is in: the workers exit while this process goes
        # on to write the peaks; the next pool's fork, or the
        # interpreter's exit, waits for them
        shut_down_later(pool)

        for chrom in self.chromosomes:
            (temp_filename, chrom_peaks) = results[chrom]
            if temp_filename is not None:
                self.pileup_data_files[chrom] = temp_filename
            for t in chrom_peaks:
                peaks.add(chrom, t[0], t[1],
                          summit=t[2],
                          peak_score=t[3],
                          pileup=t[4],
                          pscore=t[5],
                          fold_change=t[6],
                          qscore=t[7])
        return

    @cython.cfunc
    def __chrom_call_peak_using_certain_criteria(self,
                                                 peaks,
                                                 chrom: bytes,
                                                 scoring_function_s: list,
                                                 score_cutoff_s: list,
                                                 min_length: cython.int,
                                                 max_gap: cython.int,
                                                 call_summits: bool,
                                                 save_bedGraph: bool):
        """ Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object, the return value of this function
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """
        i: cython.int
        s: str
        pos_array: cnp.ndarray
        above_cutoff_index_array: cnp.ndarray
        treat_array: cnp.ndarray
        ctrl_array: cnp.ndarray
        score_array_s: list  # to: list keep different types of scores
        n_above: cython.long
        g0: cython.long
        k: cython.long
        lastp: cython.long
        ts: cython.long
        ti: cython.long
        pos_ptr: cython.pointer(cython.int)
        acia_ptr: cython.pointer(cython.int)
        treat_array_ptr: cython.pointer(cython.float)
        ctrl_array_ptr: cython.pointer(cython.float)

        assert len(scoring_function_s) == len(score_cutoff_s), "number of functions and cutoffs should be the same!"

        # first, build pileup, self.chr_pos_treat_ctrl
        # this step will be speeped up if pqtable is pre-computed.
        self.pileup_treat_ctrl_a_chromosome(chrom)
        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        # while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
        if save_bedGraph:
            self.__write_bedGraph_for_a_chromosome(chrom)

        # keep all types of scores needed
        # t0 = ttime()
        score_array_s = []
        # A q-score cutoff alone, without summits, needs no score
        # array: the segments above it are found from the q-score of
        # each slot of the p-score cache, and the q-score of a peak's
        # summit is looked up when the peak is closed. score_array_s
        # stays empty.
        threshold = None
        if (not call_summits and len(scoring_function_s) == 1 and
                scoring_function_s[0] == 'q'):
            threshold = float32_threshold(score_cutoff_s[0])
        if threshold is not None:
            above_cutoff_index_array = qscore_above_indices(self.qtable,
                                                            treat_array,
                                                            ctrl_array,
                                                            threshold)
        else:
            for i in range(len(scoring_function_s)):
                s = scoring_function_s[i]
                if s == 'p':
                    score_array_s.append(self.__cal_pscore(treat_array,
                                                           ctrl_array))
                elif s == 'q':
                    score_array_s.append(self.__cal_qscore(treat_array,
                                                           ctrl_array))
                elif s == 'f':
                    score_array_s.append(self.__cal_FE(treat_array,
                                                       ctrl_array))
                elif s == 's':
                    score_array_s.append(self.__cal_subtraction(treat_array,
                                                                ctrl_array))

            # get the regions with scores above cutoffs: the indices
            # of the segments above them, in increasing order
            above_cutoff_index_array = mask_indices(
                apply_multiple_cutoffs(score_array_s, score_cutoff_s))
        n_above = above_cutoff_index_array.shape[0]

        if n_above == 0:
            # nothing above cutoff
            return

        # Segment ti spans [pos[ti - 1], pos[ti]), and segment 0
        # [0, pos[0]). A segment above cutoff joins the region of the
        # one before it when the gap between them is at most max_gap.
        # Each region, segments [g0, g1) of above_cutoff_index_array,
        # is closed when the next segment is too far or at the end.
        acia_ptr = cython.cast(cython.pointer(cython.int),
                               above_cutoff_index_array.data)
        pos_ptr = cython.cast(cython.pointer(cython.int), pos_array.data)
        treat_array_ptr = cython.cast(cython.pointer(cython.float),
                                      treat_array.data)
        ctrl_array_ptr = cython.cast(cython.pointer(cython.float),
                                     ctrl_array.data)

        g0 = 0
        lastp = pos_ptr[acia_ptr[0]]
        for k in range(1, n_above):
            ti = acia_ptr[k]
            ts = pos_ptr[ti - 1]
            if ts - lastp > max_gap:
                # smooth length is min_length, i.e. fragment size 'd'
                self.__close_peak_region(peaks, chrom, g0, k, acia_ptr,
                                         pos_ptr, treat_array_ptr,
                                         ctrl_array_ptr, min_length,
                                         call_summits, score_array_s,
                                         score_cutoff_s)
                g0 = k
            lastp = pos_ptr[ti]
        # save the last peak
        self.__close_peak_region(peaks, chrom, g0, n_above, acia_ptr,
                                 pos_ptr, treat_array_ptr, ctrl_array_ptr,
                                 min_length, call_summits, score_array_s,
                                 score_cutoff_s)
        return

    @cython.cfunc
    def __close_peak_region(self,
                            peaks,
                            chrom: bytes,
                            g0: cython.long,
                            g1: cython.long,
                            acia_ptr: cython.pointer(cython.int),
                            pos_ptr: cython.pointer(cython.int),
                            treat_ptr: cython.pointer(cython.float),
                            ctrl_ptr: cython.pointer(cython.float),
                            min_length: cython.int,
                            call_summits: bool,
                            score_array_s: list,
                            score_cutoff_s: list):
        """Close the region of the segments above cutoff
        acia_ptr[g0:g1], with smoothing length min_length.

        With call_summits, the region goes to __close_peak_with_subpeaks
        as the list of its segments' (start, end, treat, ctrl, index)
        tuples. Otherwise it is closed here as __close_peak_wo_subpeaks
        closes that list, without building it.
        """
        k: cython.long
        kk: cython.long
        k_reset: cython.long
        n_tie: cython.long
        midindex: cython.long
        n_seen: cython.long
        start: cython.long
        end: cython.long
        ts: cython.long
        te: cython.long
        tii: cython.long
        tp: cython.float
        cp: cython.float
        ti: cython.int
        tstart: cython.int
        tend: cython.int
        summit_pos: cython.int
        i: cython.int
        tscore: cython.double
        summit_value: cython.double
        summit_treat: cython.double
        summit_ctrl: cython.double
        summit_p_score: cython.double
        summit_q_score: cython.double
        peak_content: list

        if call_summits:
            peak_content = []
            for k in range(g0, g1):
                tii = acia_ptr[k]
                ts = pos_ptr[tii - 1] if tii > 0 else 0
                te = pos_ptr[tii]
                tp = treat_ptr[tii]
                cp = ctrl_ptr[tii]
                peak_content.append((ts, te, tp, cp, tii))
            self.__close_peak_with_subpeaks(peak_content,
                                            peaks,
                                            min_length,
                                            chrom,
                                            min_length,
                                            score_array_s,
                                            score_cutoff_s=score_cutoff_s)
            return

        ti = acia_ptr[g0]
        start = pos_ptr[ti - 1] if ti > 0 else 0
        end = pos_ptr[acia_ptr[g1 - 1]]
        if end - start < min_length:
            return  # if the peak is too small, reject it

        # The summit is the middle one of the segments, in order, that
        # hold the highest treatment pileup: k_reset is the first of
        # them and n_tie their number. A pileup of 0 so far gives way
        # to any next segment.
        summit_value = 0
        k_reset = g0
        n_tie = 0
        for k in range(g0, g1):
            tscore = treat_ptr[acia_ptr[k]]
            if summit_value == 0 or summit_value < tscore:
                k_reset = k
                n_tie = 1
                summit_value = tscore
            elif summit_value == tscore:
                n_tie += 1
        midindex = (n_tie + 1) // 2 - 1
        kk = k_reset
        n_seen = 0
        k = k_reset + 1
        while n_seen < midindex:
            tscore = treat_ptr[acia_ptr[k]]
            if tscore == summit_value:
                n_seen += 1
                kk = k
            k += 1

        ti = acia_ptr[kk]
        tstart = pos_ptr[ti - 1] if ti > 0 else 0
        tend = pos_ptr[ti]
        summit_pos = (tend + tstart) // 2
        summit_treat = treat_ptr[ti]
        summit_ctrl = ctrl_ptr[ti]

        # this is a double-check to see if the summit can pass cutoff values.
        if score_array_s:
            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > score_array_s[i][ti]:
                    return  # not passed, then disgard this peak.

        summit_p_score = get_pscore(cython.cast(cython.int,
                                                summit_treat),
                                    summit_ctrl)
        summit_q_score = self.qtable.get(summit_p_score)
        if not score_array_s:
            # no score array (a q-score cutoff alone): segment ti's
            # entry in it would be summit_q_score, as a numpy float32
            if score_cutoff_s[0] > np.float32(summit_q_score):
                return

        peaks.add(chrom,           # chromosome
                  start,           # start
                  end,             # end
                  summit=summit_pos,     # summit position
                  peak_score=summit_q_score,  # score at summit
                  pileup=summit_treat,    # pileup
                  pscore=summit_p_score,  # pvalue
                  fold_change=(summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount),  # fold change
                  qscore=summit_q_score  # qvalue
                  )
        return

    @cython.cfunc
    def __close_peak_wo_subpeaks(self,
                                 peak_content: list,
                                 peaks,
                                 min_length: cython.int,
                                 chrom: bytes,
                                 smoothlen: cython.int,
                                 score_array_s: list,
                                 score_cutoff_s: list = []) -> bool:
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peak_content contains [start, end, treat_p, ctrl_p, index_in_score_array]

        peaks: a PeakIO object

        """
        summit_pos: cython.int
        tstart: cython.int
        tend: cython.int
        summit_index: cython.int
        i: cython.int
        midindex: cython.int
        ttreat_p: cython.double
        tctrl_p: cython.double
        tscore: cython.double
        summit_treat: cython.double
        summit_ctrl: cython.double
        summit_p_score: cython.double
        summit_q_score: cython.double
        tlist_scores_p: cython.int

        peak_length = peak_content[-1][1] - peak_content[0][0]
        if peak_length >= min_length:  # if the peak is too small, reject it
            tsummit = []
            summit_pos = 0
            summit_value = 0
            for i in range(len(peak_content)):
                (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
                tscore = ttreat_p  # use pscore as general score to find summit
                if not summit_value or summit_value < tscore:
                    tsummit = [(tend + tstart) // 2,]
                    tsummit_index = [i,]
                    summit_value = tscore
                elif summit_value == tscore:
                    # remember continuous summit values
                    tsummit.append((tend + tstart) // 2)
                    tsummit_index.append(i)
            # the middle of all highest points in peak region is defined as summit
            midindex = (len(tsummit) + 1) // 2 - 1
            summit_pos = tsummit[midindex]
            summit_index = tsummit_index[midindex]

            summit_treat = peak_content[summit_index][2]
            summit_ctrl = peak_content[summit_index][3]

            # this is a double-check to see if the summit can pass cutoff values.
            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > score_array_s[i][peak_content[summit_index][4]]:
                    return False  # not passed, then disgard this peak.

            summit_p_score = get_pscore(cython.cast(cython.int,
                                                    summit_treat),
                                        summit_ctrl)
            summit_q_score = self.qtable.get(summit_p_score)

            peaks.add(chrom,           # chromosome
                      peak_content[0][0],           # start
                      peak_content[-1][1],          # end
                      summit=summit_pos,     # summit position
                      peak_score=summit_q_score,  # score at summit
                      pileup=summit_treat,    # pileup
                      pscore=summit_p_score,  # pvalue
                      fold_change=(summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount),  # fold change
                      qscore=summit_q_score  # qvalue
                      )
            # start a new peak
            return True

    @cython.cfunc
    def __close_peak_with_subpeaks(self,
                                   peak_content: list,
                                   peaks,
                                   min_length: cython.int,
                                   chrom: bytes,
                                   smoothlen: cython.int,
                                   score_array_s: list,
                                   score_cutoff_s: list = [],
                                   min_valley: cython.float = 0.9) -> bool:
        """Algorithm implemented by Ben, to profile the pileup signals
        within a peak region then find subpeak summits. This method is
        highly recommended for TFBS or DNAase I sites.

        """
        tstart: cython.int
        tend: cython.int
        summit_index: cython.int
        summit_offset: cython.int
        peak_start: cython.int
        peak_end: cython.int
        start: cython.int
        end: cython.int
        i: cython.int
        start_boundary: cython.int
        m: cython.int
        n: cython.int
        ttreat_p: cython.double
        tctrl_p: cython.double
        tscore: cython.double
        summit_treat: cython.double
        summit_ctrl: cython.double
        summit_p_score: cython.double
        summit_q_score: cython.double
        peakdata: cnp.ndarray(cython.float, ndim=1)
        peakindices: cnp.ndarray(cython.int, ndim=1)
        summit_offsets: cnp.ndarray(cython.int, ndim=1)
        mapped_summits: cnp.ndarray
        tlist_scores_p: cython.int

        peak_start = peak_content[0][0]
        peak_end = peak_content[-1][1]
        peak_length = peak_end - peak_start

        if peak_length < min_length:
            return  # if the region is too small, reject it

        # Add 10 bp padding to peak region so that we can get true minima
        start = max(peak_start - 10, 0)
        end = peak_end + 10
        start_boundary = peak_start - start

        # save the scores (qscore) for each position in this region
        peakdata = np.zeros(end - start, dtype='f4')
        # save the indices for each position in this region
        # Use -1 for positions that do not belong to an above-cutoff chunk.
        # Smoothing may place a maximum inside a gap between such chunks, and
        # zero would incorrectly associate that maximum with peak_content[0].
        peakindices = np.full(end - start, -1, dtype='i4')
        for i in range(len(peak_content)):
            (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
            tscore = ttreat_p  # use pileup as general score to find summit
            m = tstart - start
            n = tend - start
            peakdata[m:n] = tscore
            peakindices[m:n] = i

        # offsets are the indices for summits in peakdata/peakindices array.
        summit_offsets = maxima(peakdata, smoothlen)

        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak_wo_subpeaks(peak_content,
                                                 peaks,
                                                 min_length,
                                                 chrom,
                                                 smoothlen,
                                                 score_array_s,
                                                 score_cutoff_s)
        else:
            # remove maxima that occurred in padding
            m = np.searchsorted(summit_offsets,
                                start_boundary)
            n = np.searchsorted(summit_offsets,
                                peak_length + start_boundary)
            summit_offsets = summit_offsets[m:n]

        summit_offsets = enforce_peakyness(peakdata, summit_offsets)

        # print "enforced:",summit_offsets
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak_wo_subpeaks(peak_content,
                                                 peaks,
                                                 min_length,
                                                 chrom,
                                                 smoothlen,
                                                 score_array_s,
                                                 score_cutoff_s)

        # indices are those point to peak_content
        summit_indices = peakindices[summit_offsets]
        mapped_summits = summit_indices >= 0
        summit_offsets = summit_offsets[mapped_summits]
        summit_indices = summit_indices[mapped_summits]

        if summit_offsets.shape[0] == 0:
            # Smoothed maxima can occur in below-cutoff gaps. Fall back to a
            # summit selected directly from the supported peak chunks.
            return self.__close_peak_wo_subpeaks(peak_content,
                                                 peaks,
                                                 min_length,
                                                 chrom,
                                                 smoothlen,
                                                 score_array_s,
                                                 score_cutoff_s)

        for summit_offset, summit_index in list(zip(summit_offsets,
                                                    summit_indices)):

            summit_treat = peak_content[summit_index][2]
            summit_ctrl = peak_content[summit_index][3]

            summit_p_score = get_pscore(cython.cast(cython.int,
                                                    summit_treat),
                                        summit_ctrl)
            summit_q_score = self.qtable.get(summit_p_score)

            for i in range(len(score_cutoff_s)):
                if score_cutoff_s[i] > score_array_s[i][peak_content[summit_index][4]]:
                    return False  # not passed, then disgard this summit.

            peaks.add(chrom,
                      peak_content[0][0],
                      peak_content[-1][1],
                      summit=start + summit_offset,
                      peak_score=summit_q_score,
                      pileup=summit_treat,
                      pscore=summit_p_score,
                      fold_change=(summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount),  # fold change
                      qscore=summit_q_score
                      )
        # start a new peak
        return True

    @cython.cfunc
    def __cal_pscore(self,
                     array1: cnp.ndarray,
                     array2: cnp.ndarray) -> cnp.ndarray:
        """Compute ``-log10`` Poisson p-scores element-wise for two arrays."""

        i: cython.long
        array1_size: cython.long
        s: cnp.ndarray
        a1_ptr: cython.pointer(cython.float)
        a2_ptr: cython.pointer(cython.float)
        s_ptr: cython.pointer(cython.float)

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")

        a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
        a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
        s_ptr = cython.cast(cython.pointer(cython.float), s.data)

        array1_size = array1.shape[0]

        for i in range(array1_size):
            s_ptr[0] = get_pscore(cython.cast(cython.int,
                                              a1_ptr[0]),
                                  a2_ptr[0])
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    @cython.cfunc
    def __cal_qscore(self,
                     array1: cnp.ndarray,
                     array2: cnp.ndarray) -> cnp.ndarray:
        """Map p-scores to q-scores using the precomputed ``pqtable``.

        The q-score of each slot of the p-score cache is kept in
        ``slot_qscores``, so a position takes one hash lookup."""
        i: cython.long
        s: cnp.ndarray
        a1_ptr: cython.pointer(cython.float)
        a2_ptr: cython.pointer(cython.float)
        s_ptr: cython.pointer(cython.float)
        qtable: QScoreTable = self.qtable
        sq: SlotQScores
        slot: cython.Py_ssize_t
        e: cython.uint
        q: cython.float

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")

        a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
        a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
        s_ptr = cython.cast(cython.pointer(cython.float), s.data)

        sq = slot_qscores_for(qtable)
        for i in range(array1.shape[0]):
            slot = pscore_cache.index(cython.cast(cython.int, a1_ptr[0]),
                                      a2_ptr[0])
            if pscore_cache.mask != sq.mask:
                # the cache grew on a miss, which moved its slots
                sq = slot_qscores_for(qtable)
            e = sq.entries[slot]
            if e == 0:
                q = qtable.get(pscore_cache.slots[slot].pscore)
                memcpy(cython.address(e), cython.address(q), 4)
                sq.entries[slot] = ~e
                s_ptr[0] = q
            else:
                e = ~e
                memcpy(s_ptr, cython.address(e), 4)
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    @cython.cfunc
    def __cal_logLR(self,
                    array1: cnp.ndarray,
                    array2: cnp.ndarray) -> cnp.ndarray:
        """Compute asymmetric log-likelihood ratios element-wise."""
        i: cython.long
        s: cnp.ndarray
        a1_ptr: cython.pointer(cython.float)
        a2_ptr: cython.pointer(cython.float)
        s_ptr: cython.pointer(cython.float)

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")

        a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
        a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
        s_ptr = cython.cast(cython.pointer(cython.float), s.data)

        for i in range(array1.shape[0]):
            s_ptr[0] = get_logLR_asym((a1_ptr[0] + self.pseudocount,
                                       a2_ptr[0] + self.pseudocount))
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    @cython.cfunc
    def __cal_logFE(self,
                    array1: cnp.ndarray,
                    array2: cnp.ndarray) -> cnp.ndarray:
        """Compute log fold enrichment with the configured pseudocount."""
        i: cython.long
        s: cnp.ndarray
        a1_ptr: cython.pointer(cython.float)
        a2_ptr: cython.pointer(cython.float)
        s_ptr: cython.pointer(cython.float)

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")

        a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
        a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
        s_ptr = cython.cast(cython.pointer(cython.float), s.data)

        for i in range(array1.shape[0]):
            s_ptr[0] = get_logFE(a1_ptr[0] + self.pseudocount,
                                 a2_ptr[0] + self.pseudocount)
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    @cython.cfunc
    def __cal_FE(self,
                 array1: cnp.ndarray,
                 array2: cnp.ndarray) -> cnp.ndarray:
        """Compute linear fold enrichment with the configured pseudocount."""
        i: cython.long
        s: cnp.ndarray
        a1_ptr: cython.pointer(cython.float)
        a2_ptr: cython.pointer(cython.float)
        s_ptr: cython.pointer(cython.float)

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")

        a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
        a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
        s_ptr = cython.cast(cython.pointer(cython.float), s.data)

        for i in range(array1.shape[0]):
            s_ptr[0] = (a1_ptr[0] + self.pseudocount) / (a2_ptr[0] + self.pseudocount)
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    @cython.cfunc
    def __cal_subtraction(self,
                          array1: cnp.ndarray,
                          array2: cnp.ndarray) -> cnp.ndarray:
        """Compute treatment-control subtraction element-wise."""
        i: cython.long
        s: cnp.ndarray
        a1_ptr: cython.pointer(cython.float)
        a2_ptr: cython.pointer(cython.float)
        s_ptr: cython.pointer(cython.float)

        assert array1.shape[0] == array2.shape[0]
        s = np.zeros(array1.shape[0], dtype="f4")

        a1_ptr = cython.cast(cython.pointer(cython.float), array1.data)
        a2_ptr = cython.cast(cython.pointer(cython.float), array2.data)
        s_ptr = cython.cast(cython.pointer(cython.float), s.data)

        for i in range(array1.shape[0]):
            s_ptr[0] = a1_ptr[0] - a2_ptr[0]
            s_ptr += 1
            a1_ptr += 1
            a2_ptr += 1
        return s

    @cython.cfunc
    def __write_bedGraph_for_a_chromosome(self, chrom: bytes) -> bool:
        """Write treat/control values for a certain chromosome into a
        specified file handler.

        """
        pos_array: cnp.ndarray
        treat_array: cnp.ndarray
        ctrl_array: cnp.ndarray
        pos_array_ptr: cython.pointer(cython.int)
        treat_array_ptr: cython.pointer(cython.float)
        ctrl_array_ptr: cython.pointer(cython.float)
        l: cython.int
        i: cython.int
        p: cython.int
        pre_p_t: cython.int
        # current position, previous position for treat, previous position for control
        pre_p_c: cython.int
        pre_v_t: cython.float
        pre_v_c: cython.float
        v_t: cython.float
        # previous value for treat, for control, current value for treat, for control
        v_c: cython.float
        # 1 if save_SPMR is false, or depth in million if save_SPMR is
        # true. Note, while piling up and calling peaks, treatment and
        # control have been scaled to the same depth, so we need to
        # find what this 'depth' is.
        denominator: cython.float
        ft: cython.pointer(FILE)
        fc: cython.pointer(FILE)

        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl
        pos_array_ptr = cython.cast(cython.pointer(cython.int),
                                    pos_array.data)
        treat_array_ptr = cython.cast(cython.pointer(cython.float),
                                      treat_array.data)
        ctrl_array_ptr = cython.cast(cython.pointer(cython.float),
                                     ctrl_array.data)

        if self.save_SPMR:
            if self.treat_scaling_factor == 1:
                # in this case, control has been asked to be scaled to depth of treatment
                denominator = self.treat.total/1e6
            else:
                # in this case, treatment has been asked to be scaled to depth of control
                denominator = self.ctrl.total/1e6
        else:
            denominator = 1.0

        l = pos_array.shape[0]

        if l == 0:              # if there is no data, return
            return False

        ft = self.bedGraph_treat_f
        fc = self.bedGraph_ctrl_f
        # t_write_func = self.bedGraph_treat.write
        # c_write_func = self.bedGraph_ctrl.write

        pre_p_t = 0
        pre_p_c = 0
        pre_v_t = treat_array_ptr[0]/denominator
        pre_v_c = ctrl_array_ptr[0]/denominator
        treat_array_ptr += 1
        ctrl_array_ptr += 1

        for i in range(1, l):
            v_t = treat_array_ptr[0]/denominator
            v_c = ctrl_array_ptr[0]/denominator
            p = pos_array_ptr[0]
            pos_array_ptr += 1
            treat_array_ptr += 1
            ctrl_array_ptr += 1

            if abs(pre_v_t - v_t) > 1e-5:  # precision is 5 digits
                fprintf(ft, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_t, p, pre_v_t)
                pre_v_t = v_t
                pre_p_t = p

            if abs(pre_v_c - v_c) > 1e-5:  # precision is 5 digits
                fprintf(fc, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_c, p, pre_v_c)
                pre_v_c = v_c
                pre_p_c = p

        p = pos_array_ptr[0]
        # last one
        fprintf(ft, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_t, p, pre_v_t)
        fprintf(fc, b"%s\t%d\t%d\t%.5f\n", chrom, pre_p_c, p, pre_v_c)

        return True

    @cython.ccall
    def call_broadpeaks(self,
                        scoring_function_symbols: list,
                        lvl1_cutoff_s: list,
                        lvl2_cutoff_s: list,
                        min_length: cython.int = 200,
                        lvl1_max_gap: cython.int = 50,
                        lvl2_max_gap: cython.int = 400,
                        cutoff_analysis: bool = False):
        """This function try to find enriched regions within which,
        scores are continuously higher than a given cutoff for level
        1, and link them using the gap above level 2 cutoff with a
        maximum length of lvl2_max_gap.

        Args:
            scoring_function_symbols: Symbols for score functions.
                Use ``'p'`` (pscore), ``'q'`` (qscore), ``'f'`` (fold change),
                or ``'s'`` (subtraction). Example: ``['p', 'q']``.
            lvl1_cutoff_s: Cutoffs for highly enriched regions.
            lvl2_cutoff_s: Cutoffs for linkage regions.
            min_length: Minimum peak length.
            lvl1_max_gap: Maximum gap to merge nearby peaks.
            lvl2_max_gap: Maximum length of linkage regions.
            cutoff_analysis: Whether to compute cutoff-vs-peak metrics.

        Returns:
            tuple: ``(PeakIO, BroadPeakIO)`` for level-1 peaks and broad regions.

        Examples:
            .. code-block:: python

                lvl1, broad = caller.call_broadpeaks(['p'], [5.0], [2.0])
        """
        i: cython.int
        j: cython.int
        chrom: bytes
        lvl1peaks: object
        lvl1peakschrom: object
        lvl1: object
        lvl2peaks: object
        lvl2peakschrom: object
        lvl2: object
        broadpeaks: object
        chrs: set
        tmppeakset: list

        lvl1peaks = PeakIO()
        lvl2peaks = PeakIO()

        # prepare p-q table
        if self.qtable.n == 0:
            info("#3 Pre-compute pvalue-qvalue table...")
            if cutoff_analysis:
                info("#3 Cutoff value vs broad region calls will be analyzed!")
                self.__pre_computes(max_gap=lvl2_max_gap, min_length=min_length)
            else:
                self.__cal_pvalue_qvalue_table()

        # prepare bedGraph file
        if self.save_bedGraph:

            self.bedGraph_treat_f = fopen(self.bedGraph_treat_filename, "w")
            self.bedGraph_ctrl_f = fopen(self.bedGraph_control_filename, "w")
            info("#3 In the peak calling step, the following will be performed simultaneously:")
            info("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... %s" % self.bedGraph_filename_prefix.decode() + "_treat_pileup.bdg")
            info("#3   Write bedGraph files for control lambda (after scaling if necessary)... %s" % self.bedGraph_filename_prefix.decode() + "_control_lambda.bdg")

            if self.trackline:
                # this line is REQUIRED by the wiggle format for UCSC browser
                tmp_bytes = ("track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf(self.bedGraph_treat_f, tmp_bytes)
                tmp_bytes = ("track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'%s\'\"\n" % self.bedGraph_filename_prefix).encode()
                fprintf(self.bedGraph_ctrl_f, tmp_bytes)

        info("#3 Call peaks for each chromosome...")
        for chrom in self.chromosomes:
            self.__chrom_call_broadpeak_using_certain_criteria(lvl1peaks,
                                                               lvl2peaks,
                                                               chrom,
                                                               scoring_function_symbols,
                                                               lvl1_cutoff_s,
                                                               lvl2_cutoff_s,
                                                               min_length,
                                                               lvl1_max_gap,
                                                               lvl2_max_gap,
                                                               self.save_bedGraph)

        # close bedGraph file
        if self.save_bedGraph:
            fclose(self.bedGraph_treat_f)
            fclose(self.bedGraph_ctrl_f)
            # self.bedGraph_ctrl.close()
            self.save_bedGraph = False

        # now combine lvl1 and lvl2 peaks
        chrs = lvl1peaks.get_chr_names()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in sorted(chrs):
            lvl1peakschrom = lvl1peaks.get_data_from_chrom(chrom)
            lvl2peakschrom = lvl2peaks.get_data_from_chrom(chrom)
            lvl1peakschrom_next = iter(lvl1peakschrom).__next__
            tmppeakset = []             # to temporarily store lvl1 region inside a lvl2 region
            # our assumption is lvl1 regions should be included in lvl2 regions
            try:
                lvl1 = lvl1peakschrom_next()
                for i in range(len(lvl2peakschrom)):
                    # for each lvl2 peak, find all lvl1 peaks inside
                    # I assume lvl1 peaks can be ALL covered by lvl2 peaks.
                    lvl2 = lvl2peakschrom[i]

                    while True:
                        if lvl2["start"] <= lvl1["start"] and lvl1["end"] <= lvl2["end"]:
                            tmppeakset.append(lvl1)
                            lvl1 = lvl1peakschrom_next()
                        else:
                            # make a hierarchical broad peak
                            # print lvl2["start"], lvl2["end"], lvl2["score"]
                            self.__add_broadpeak(broadpeaks,
                                                 chrom,
                                                 lvl2,
                                                 tmppeakset)
                            tmppeakset = []
                            break
            except StopIteration:
                # no more strong (aka lvl1) peaks left
                self.__add_broadpeak(broadpeaks,
                                     chrom,
                                     lvl2,
                                     tmppeakset)
                tmppeakset = []
                # add the rest lvl2 peaks
                for j in range(i+1, len(lvl2peakschrom)):
                    self.__add_broadpeak(broadpeaks,
                                         chrom,
                                         lvl2peakschrom[j],
                                         tmppeakset)

        return broadpeaks

    @cython.cfunc
    def __chrom_call_broadpeak_using_certain_criteria(self,
                                                      lvl1peaks,
                                                      lvl2peaks,
                                                      chrom: bytes,
                                                      scoring_function_s: list,
                                                      lvl1_cutoff_s: list,
                                                      lvl2_cutoff_s: list,
                                                      min_length: cython.int,
                                                      lvl1_max_gap: cython.int,
                                                      lvl2_max_gap: cython.int,
                                                      save_bedGraph: bool):
        """Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object

        scoring_function_s: symbols of functions to calculate score as
        score=f(x, y) where x is treatment pileup, and y is control
        pileup

        save_bedGraph : whether or not to save pileup and control into
        a bedGraph file

        """
        i: cython.int
        s: str
        above_cutoff: cnp.ndarray
        above_cutoff_endpos: cnp.ndarray
        above_cutoff_startpos: cnp.ndarray
        pos_array: cnp.ndarray
        treat_array: cnp.ndarray
        ctrl_array: cnp.ndarray
        above_cutoff_index_array: cnp.ndarray
        score_array_s: list          # to: list keep different types of scores
        peak_content: list
        acs_ptr: cython.pointer(cython.int)
        ace_ptr: cython.pointer(cython.int)
        acia_ptr: cython.pointer(cython.int)
        treat_array_ptr: cython.pointer(cython.float)
        ctrl_array_ptr: cython.pointer(cython.float)

        assert len(scoring_function_s) == len(lvl1_cutoff_s), "number of functions and cutoffs should be the same!"
        assert len(scoring_function_s) == len(lvl2_cutoff_s), "number of functions and cutoffs should be the same!"

        # first, build pileup, self.chr_pos_treat_ctrl
        self.pileup_treat_ctrl_a_chromosome(chrom)
        [pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl

        # while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
        if save_bedGraph:
            self.__write_bedGraph_for_a_chromosome(chrom)

        # keep all types of scores needed
        score_array_s = []
        for i in range(len(scoring_function_s)):
            s = scoring_function_s[i]
            if s == 'p':
                score_array_s.append(self.__cal_pscore(treat_array,
                                                       ctrl_array))
            elif s == 'q':
                score_array_s.append(self.__cal_qscore(treat_array,
                                                       ctrl_array))
            elif s == 'f':
                score_array_s.append(self.__cal_FE(treat_array,
                                                   ctrl_array))
            elif s == 's':
                score_array_s.append(self.__cal_subtraction(treat_array,
                                                            ctrl_array))

        # lvl1 : strong peaks
        peak_content = []           # to store points above cutoff

        # get the regions with scores above cutoffs
        above_cutoff = np.nonzero(apply_multiple_cutoffs(score_array_s,
                                                         lvl1_cutoff_s))[0]  # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff_index_array = np.arange(pos_array.shape[0],
                                             dtype="int32")[above_cutoff]  # indices
        above_cutoff_endpos = pos_array[above_cutoff]  # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1]  # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        acs_ptr = cython.cast(cython.pointer(cython.int),
                              above_cutoff_startpos.data)
        ace_ptr = cython.cast(cython.pointer(cython.int),
                              above_cutoff_endpos.data)
        acia_ptr = cython.cast(cython.pointer(cython.int),
                               above_cutoff_index_array.data)
        treat_array_ptr = cython.cast(cython.pointer(cython.float),
                                      treat_array.data)
        ctrl_array_ptr = cython.cast(cython.pointer(cython.float),
                                     ctrl_array.data)

        ts = acs_ptr[0]
        te = ace_ptr[0]
        ti = acia_ptr[0]
        tp = treat_array_ptr[ti]
        cp = ctrl_array_ptr[ti]

        peak_content.append((ts, te, tp, cp, ti))
        acs_ptr += 1            # move ptr
        ace_ptr += 1
        acia_ptr += 1
        lastp = te

        # peak_content.append((above_cutoff_startpos[0], above_cutoff_endpos[0], treat_array[above_cutoff_index_array[0]], ctrl_array[above_cutoff_index_array[0]], score_array_s, above_cutoff_index_array[0]))
        for i in range(1, above_cutoff_startpos.size):
            ts = acs_ptr[0]
            te = ace_ptr[0]
            ti = acia_ptr[0]
            acs_ptr += 1
            ace_ptr += 1
            acia_ptr += 1
            tp = treat_array_ptr[ti]
            cp = ctrl_array_ptr[ti]
            tl = ts - lastp
            if tl <= lvl1_max_gap:
                # append
                peak_content.append((ts, te, tp, cp, ti))
                lastp = te
            else:
                # close
                self.__close_peak_for_broad_region(peak_content,
                                                   lvl1peaks,
                                                   min_length,
                                                   chrom,
                                                   lvl1_max_gap//2,
                                                   score_array_s)
                peak_content = [(ts, te, tp, cp, ti),]
                lastp = te      # above_cutoff_endpos[i]

        # save the last peak
        if peak_content:
            self.__close_peak_for_broad_region(peak_content,
                                               lvl1peaks,
                                               min_length,
                                               chrom,
                                               lvl1_max_gap//2,
                                               score_array_s)

        # lvl2 : weak peaks
        peak_content = []           # to store points above cutoff

        # get the regions with scores above cutoffs

        # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
        above_cutoff = np.nonzero(apply_multiple_cutoffs(score_array_s,
                                                         lvl2_cutoff_s))[0]

        above_cutoff_index_array = np.arange(pos_array.shape[0],
                                             dtype="i4")[above_cutoff] # indices
        above_cutoff_endpos = pos_array[above_cutoff]  # end positions of regions where score is above cutoff
        above_cutoff_startpos = pos_array[above_cutoff-1]  # start positions of regions where score is above cutoff

        if above_cutoff.size == 0:
            # nothing above cutoff
            return

        if above_cutoff[0] == 0:
            # first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
            above_cutoff_startpos[0] = 0

        # first bit of region above cutoff
        acs_ptr = cython.cast(cython.pointer(cython.int),
                              above_cutoff_startpos.data)
        ace_ptr = cython.cast(cython.pointer(cython.int),
                              above_cutoff_endpos.data)
        acia_ptr = cython.cast(cython.pointer(cython.int),
                               above_cutoff_index_array.data)
        treat_array_ptr = cython.cast(cython.pointer(cython.float),
                                      treat_array.data)
        ctrl_array_ptr = cython.cast(cython.pointer(cython.float),
                                     ctrl_array.data)

        ts = acs_ptr[0]
        te = ace_ptr[0]
        ti = acia_ptr[0]
        tp = treat_array_ptr[ti]
        cp = ctrl_array_ptr[ti]
        peak_content.append((ts, te, tp, cp, ti))
        acs_ptr += 1            # move ptr
        ace_ptr += 1
        acia_ptr += 1

        lastp = te
        for i in range(1, above_cutoff_startpos.size):
            # for everything above cutoff
            ts = acs_ptr[0]     # get the start
            te = ace_ptr[0]     # get the end
            ti = acia_ptr[0]    # get the index

            acs_ptr += 1        # move ptr
            ace_ptr += 1
            acia_ptr += 1
            tp = treat_array_ptr[ti]  # get the treatment pileup
            cp = ctrl_array_ptr[ti]  # get the control pileup
            tl = ts - lastp  # get the distance from the current point to last position of existing peak_content

            if tl <= lvl2_max_gap:
                # append
                peak_content.append((ts, te, tp, cp, ti))
                lastp = te
            else:
                # close
                self.__close_peak_for_broad_region(peak_content,
                                                   lvl2peaks,
                                                   min_length,
                                                   chrom,
                                                   lvl2_max_gap//2,
                                                   score_array_s)

                peak_content = [(ts, te, tp, cp, ti),]
                lastp = te

        # save the last peak
        if peak_content:
            self.__close_peak_for_broad_region(peak_content,
                                               lvl2peaks,
                                               min_length,
                                               chrom,
                                               lvl2_max_gap//2,
                                               score_array_s)

        return

    @cython.cfunc
    def __close_peak_for_broad_region(self,
                                      peak_content: list,
                                      peaks,
                                      min_length: cython.int,
                                      chrom: bytes,
                                      smoothlen: cython.int,
                                      score_array_s: list,
                                      score_cutoff_s: list = []) -> bool:
        """Close the broad peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        peak_content contains [start, end, treat_p, ctrl_p, list_scores]

        peaks: a BroadPeakIO object

        """
        tstart: cython.int
        tend: cython.int
        i: cython.int
        ttreat_p: cython.double
        tctrl_p: cython.double
        tlist_pileup: list
        tlist_control: list
        tlist_length: list
        tlist_scores_p: cython.int
        tarray_pileup: cnp.ndarray
        tarray_control: cnp.ndarray
        tarray_pscore: cnp.ndarray
        tarray_qscore: cnp.ndarray
        tarray_fc: cnp.ndarray

        peak_length = peak_content[-1][1] - peak_content[0][0]
        if peak_length >= min_length:  # if the peak is too small, reject it
            tlist_pileup = []
            tlist_control = []
            tlist_length = []
            for i in range(len(peak_content)):  # each position in broad peak
                (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
                tlist_pileup.append(ttreat_p)
                tlist_control.append(tctrl_p)
                tlist_length.append(tend - tstart)

            tarray_pileup = np.array(tlist_pileup, dtype="f4")
            tarray_control = np.array(tlist_control, dtype="f4")
            tarray_pscore = self.__cal_pscore(tarray_pileup, tarray_control)
            tarray_qscore = self.__cal_qscore(tarray_pileup, tarray_control)
            tarray_fc = self.__cal_FE(tarray_pileup, tarray_control)

            peaks.add(chrom,           # chromosome
                      peak_content[0][0],  # start
                      peak_content[-1][1],  # end
                      summit=0,
                      peak_score=mean_from_value_length(tarray_qscore, tlist_length),
                      pileup=mean_from_value_length(tarray_pileup, tlist_length),
                      pscore=mean_from_value_length(tarray_pscore, tlist_length),
                      fold_change=mean_from_value_length(tarray_fc, tlist_length),
                      qscore=mean_from_value_length(tarray_qscore, tlist_length),
                      )
            # if chrom == "chr1" and  peak_content[0][0] == 237643 and peak_content[-1][1] == 237935:
            #    print tarray_qscore, tlist_length
            # start a new peak
            return True

    @cython.cfunc
    def __add_broadpeak(self,
                        bpeaks,
                        chrom: bytes,
                        lvl2peak: object,
                        lvl1peakset: list):
        """Internal function to create broad peak.

        *Note* lvl1peakset/strong_regions might be empty
        """

        blockNum: cython.int
        start: cython.int
        end: cython.int
        blockSizes: bytes
        blockStarts: bytes
        thickStart: bytes
        thickEnd: bytes

        start = lvl2peak["start"]
        end = lvl2peak["end"]

        if not lvl1peakset:
            # will complement by adding 1bps start and end to this region
            # may change in the future if gappedPeak format was improved.
            bpeaks.add(chrom, start, end,
                       score=lvl2peak["score"],
                       thickStart=(b"%d" % start),
                       thickEnd=(b"%d" % end),
                       blockNum=2,
                       blockSizes=b"1,1",
                       blockStarts=(b"0,%d" % (end-start-1)),
                       pileup=lvl2peak["pileup"],
                       pscore=lvl2peak["pscore"],
                       fold_change=lvl2peak["fc"],
                       qscore=lvl2peak["qscore"])
            return bpeaks

        thickStart = b"%d" % (lvl1peakset[0]["start"])
        thickEnd = b"%d" % (lvl1peakset[-1]["end"])
        blockNum = len(lvl1peakset)
        blockSizes = b",".join([b"%d" % y for y in [x["length"] for x in lvl1peakset]])
        blockStarts = b",".join([b"%d" % x for x in getitem_then_subtract(lvl1peakset, start)])

        # add 1bp left and/or right block if necessary
        if int(thickStart) != start:
            # add 1bp left block
            thickStart = b"%d" % start
            blockNum += 1
            blockSizes = b"1,"+blockSizes
            blockStarts = b"0,"+blockStarts
        if int(thickEnd) != end:
            # add 1bp right block
            thickEnd = b"%d" % end
            blockNum += 1
            blockSizes = blockSizes + b",1"
            blockStarts = blockStarts + b"," + (b"%d" % (end-start-1))

        bpeaks.add(chrom, start, end,
                   score=lvl2peak["score"],
                   thickStart=thickStart,
                   thickEnd=thickEnd,
                   blockNum=blockNum,
                   blockSizes=blockSizes,
                   blockStarts=blockStarts,
                   pileup=lvl2peak["pileup"],
                   pscore=lvl2peak["pscore"],
                   fold_change=lvl2peak["fc"],
                   qscore=lvl2peak["qscore"])
        return bpeaks


def pscore_stat_worker(chroms: list) -> tuple:
    """Worker side of CallerFromAlignments._pscore_stat_in_workers.

    Runs in a forked process. Adds the chromosomes whose indices are
    in ``chroms`` to this process's p-score histogram, in increasing
    index order, without computing p-scores, and returns [(name, file
    holding its pileup or None)] for them with the histogram's key,
    first-segment and length arrays. Taking the histogram resets it
    for the worker's next task.
    """
    caller: CallerFromAlignments = _caller_for_workers
    i: cython.long
    chrom: bytes
    files: list = []

    for i in sorted(chroms):
        caller._pscore_stat_a_chromosome(i, False)
        chrom = caller.chromosomes[i]
        files.append((chrom, caller.pileup_data_files.get(chrom)))
    (keys, firsts, lengths) = pscore_cache.take_histogram(True)
    return (files, keys, firsts, lengths)


def call_peaks_worker(task: tuple) -> list:
    """Worker side of CallerFromAlignments._call_peaks_in_workers.

    Runs in a forked process. ``task`` is (chromosome indices, scoring
    function symbols, cutoffs, min_length, max_gap, call_summits). For
    each chromosome, return (name, file holding its pileup or None,
    its peaks).
    """
    caller: CallerFromAlignments = _caller_for_workers
    i: cython.long
    chrom: bytes
    out: list = []

    (chroms, scoring_function_symbols, score_cutoff_s,
     min_length, max_gap, call_summits) = task
    for i in chroms:
        chrom = caller.chromosomes[i]
        chrom_peaks = caller._call_peaks_a_chromosome(chrom,
                                                      scoring_function_symbols,
                                                      score_cutoff_s,
                                                      min_length,
                                                      max_gap,
                                                      call_summits)
        out.append((chrom, caller.pileup_data_files.get(chrom), chrom_peaks))
    return out
