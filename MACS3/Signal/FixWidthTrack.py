# cython: language_level=3
# Time-stamp: <2025-11-10 15:24:27 Tao Liu>

"""Module for FWTrack classes.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import io

# ------------------------------------
# MACS3 modules
# ------------------------------------

from MACS3.IO.PeakIO import PeakIO
from MACS3.Signal.PileupV2 import pileup_from_PN_shifted, over_two_pv_array
from MACS3.Signal.Pileup import se_all_in_one_pileup_max3
from MACS3.Utilities.Threads import MIN_THREADED_SIZE, map_chromosomes

# ------------------------------------
# Other modules
# ------------------------------------
import cython
import numpy as np
from cython.cimports.cpython import bool
import cython.cimports.numpy as cnp
from cython.cimports.libc.stdint import INT32_MAX as INT_MAX
from cython.cimports.MACS3.Signal.crefcount import Py_REFCNT

# ------------------------------------
# constants
# ------------------------------------
# finalize and filter_dup hand a chromosome of at least this many
# positions to map_chromosomes, and run any other on the spot
_MIN_THREADED_SIZE = cython.declare(cython.long, MIN_THREADED_SIZE)

# ------------------------------------
# Misc functions
# ------------------------------------

# The C fast path of FWTrack.filter_dup. FWTrack._filter_dup_strand runs
# it through filter_dup_plain on a plain array (is_plain_i4_array), in
# place when the array owns its data (owns_buffer) and nothing else
# holds it, and otherwise runs the original loop, so any other array
# gives the original result or exception.

@cython.cfunc
@cython.inline
def is_plain_i4_array(a: cnp.ndarray) -> cython.bint:
    """Whether `a` is a one-dimensional, C-contiguous, aligned array of
    native int32, the layout filter_dup_i4 reads through a pointer."""
    return (cnp.PyArray_TYPE(a) == cnp.NPY_INT32 and
            cnp.PyArray_NDIM(a) == 1 and cnp.PyArray_ISCARRAY_RO(a))


@cython.cfunc
@cython.nogil
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def filter_dup_i4(src: cython.pointer(cython.int), size: cython.Py_ssize_t,
                  dst: cython.pointer(cython.int),
                  maxnum: cython.int) -> cython.Py_ssize_t:
    """Copy src[0:size] to dst, keeping at most `maxnum` positions of
    each run of equal values, and return how many were written.

    The loop of FWTrack.filter_dup, in C: the first position is always
    kept, even at `maxnum` 0, and the same comparisons are made on the
    same int32 values. `size` must be at least 1 and `dst` must have
    room for `size` values.
    """
    i_old: cython.Py_ssize_t
    i_new: cython.Py_ssize_t = 1
    n: cython.int = 1
    p: cython.int
    current_loc: cython.int = src[0]

    dst[0] = current_loc
    for i_old in range(1, size):
        p = src[i_old]
        if p == current_loc:
            n += 1
        else:
            current_loc = p
            n = 1
        if n <= maxnum:
            dst[i_new] = p
            i_new += 1
    return i_new


@cython.cfunc
def filter_dup_plain(src, dst, maxnum: cython.int) -> cython.Py_ssize_t:
    """Run filter_dup_i4 from array `src` into array `dst`, without
    the GIL, and return how many positions it wrote. Returns -1,
    touching nothing, unless both are plain int32 arrays
    (is_plain_i4_array), `dst` is writable and `src` is not empty and
    not longer than `dst`; the caller then runs its own loop."""
    a: cnp.ndarray
    b: cnp.ndarray
    src_p: cython.pointer(cython.int)
    dst_p: cython.pointer(cython.int)
    size: cython.Py_ssize_t
    kept: cython.Py_ssize_t

    if not isinstance(src, cnp.ndarray) or not isinstance(dst, cnp.ndarray):
        return -1
    a = src
    b = dst
    if not (is_plain_i4_array(a) and is_plain_i4_array(b) and
            cnp.PyArray_ISWRITEABLE(b)):
        return -1
    if a.shape[0] < 1 or a.shape[0] > b.shape[0]:
        return -1
    src_p = cython.cast(cython.pointer(cython.int), a.data)
    dst_p = cython.cast(cython.pointer(cython.int), b.data)
    size = a.shape[0]
    with cython.nogil:
        kept = filter_dup_i4(src_p, size, dst_p, maxnum)
    return kept


@cython.cfunc
def owns_buffer(a) -> cython.bint:
    """Whether `a` is an array that owns its data and is no view of
    another object, as `a.resize(..., refcheck=False)` requires."""
    b: cnp.ndarray

    if not isinstance(a, cnp.ndarray):
        return False
    b = a
    return (cnp.PyArray_CHKFLAGS(b, cnp.NPY_ARRAY_OWNDATA) and
            cnp.PyArray_BASE(b) == cython.NULL)

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class FWTrack:
    """Fixed-width fragment track grouped by chromosome.
    
    Stores plus- and minus-strand 5' cut positions in numpy arrays and exposes
    utilities for sorting, filtering, sampling, and pileup generation.
    """
    locations: dict
    pointer: dict
    buf_size: dict
    rlengths: dict
    is_sorted: bool
    is_destroyed: bool
    total = cython.declare(cython.ulong, visibility="public")
    annotation = cython.declare(str, visibility="public")
    buffer_size = cython.declare(cython.long, visibility="public")
    length = cython.declare(cython.ulonglong, visibility="public")
    fw = cython.declare(cython.int, visibility="public")

    def __init__(self,
                 fw: cython.int = 0,
                 anno: str = "",
                 buffer_size: cython.long = 100000):
        """Initialize an empty fixed-width track.
        
        Parameters
        ----------
        fw : int, optional
            Fixed fragment width (bp) used when estimating coverage and region length.
        anno : str, optional
            Annotation label retained with the track metadata.
        buffer_size : int, optional
            Number of positions allocated per growth chunk for each strand array.
        """
        self.fw = fw
        self.locations = {}    # location pairs: two strands
        self.pointer = {}      # location pairs
        self.buf_size = {}     # location pairs
        self.is_sorted = False
        self.total = 0           # total tags
        self.annotation = anno   # need to be figured out
        # lengths of reference sequences, e.g. each chromosome in a genome
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.is_destroyed = False

    @cython.ccall
    def destroy(self):
        """Release numpy buffers held by the track.
        
        All per-chromosome arrays are resized to zero so the memory footprint returns
        to the allocator, and the track is marked as destroyed.
        """
        chrs: set
        chromosome: bytes

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.locations:
                self.locations[chromosome][0].resize(self.buffer_size,
                                                     refcheck=False)
                self.locations[chromosome][0].resize(0,
                                                     refcheck=False)
                self.locations[chromosome][1].resize(self.buffer_size,
                                                     refcheck=False)
                self.locations[chromosome][1].resize(0,
                                                     refcheck=False)
                self.locations[chromosome] = [None, None]
                self.locations.pop(chromosome)
        self.is_destroyed = True
        return

    @cython.ccall
    def add_loc(self,
                chromosome: bytes,
                fiveendpos: cython.int,
                strand: cython.int):
        """Append a 5' cut position to the track.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) that owns the cut.
        fiveendpos : int
            Zero-based 5' coordinate of the cut site.
        strand : int
            Strand flag where ``0`` denotes plus and ``1`` denotes minus.
        
        Notes
        -----
        Positions are stored in strand-specific numpy arrays keyed by chromosome, and
        the strand pointer is advanced as new positions are appended.
        """
        i: cython.int
        b: cython.int
        arr: cnp.ndarray

        if chromosome not in self.locations:
            self.buf_size[chromosome] = [self.buffer_size, self.buffer_size]
            self.locations[chromosome] = [np.zeros(self.buffer_size, dtype='i4'),
                                          np.zeros(self.buffer_size, dtype='i4')]
            self.pointer[chromosome] = [0, 0]
            self.locations[chromosome][strand][0] = fiveendpos
            self.pointer[chromosome][strand] = 1
        else:
            i = self.pointer[chromosome][strand]
            b = self.buf_size[chromosome][strand]
            arr = self.locations[chromosome][strand]
            if b == i:
                b += self.buffer_size
                arr.resize(b, refcheck=False)
                self.buf_size[chromosome][strand] = b
            arr[i] = fiveendpos
            self.pointer[chromosome][strand] += 1
        return

    @cython.ccall
    def add_loc_arrays(self, chromosome: bytes, plus, minus):
        """Append many 5' cut positions of one chromosome to the track.

        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) that owns the cuts.
        plus : numpy.ndarray
            int32 positions of plus-strand cuts, in the order to append.
        minus : numpy.ndarray
            int32 positions of minus-strand cuts, in the order to append.

        Notes
        -----
        Leaves the track exactly as calling ``add_loc`` for each element
        of ``plus`` with strand 0 and each element of ``minus`` with
        strand 1 would: the same arrays, pointers and buffer sizes.
        """
        strand: cython.int
        i: cython.long
        n: cython.long
        b: cython.long
        arr: cnp.ndarray

        if len(plus) == 0 and len(minus) == 0:
            return
        if self.buffer_size <= 0:
            # add_loc's own behaviour, errors included
            for x in plus:
                self.add_loc(chromosome, x, 0)
            for x in minus:
                self.add_loc(chromosome, x, 1)
            return

        if chromosome not in self.locations:
            self.buf_size[chromosome] = [self.buffer_size, self.buffer_size]
            self.locations[chromosome] = [np.zeros(self.buffer_size, dtype='i4'),
                                          np.zeros(self.buffer_size, dtype='i4')]
            self.pointer[chromosome] = [0, 0]
        for strand in range(2):
            positions = minus if strand else plus
            n = len(positions)
            if n == 0:
                continue
            i = self.pointer[chromosome][strand]
            b = self.buf_size[chromosome][strand]
            arr = self.locations[chromosome][strand]
            if i + n > b:
                # grow in steps of buffer_size, as add_loc does
                while b < i + n:
                    b += self.buffer_size
                arr.resize(b, refcheck=False)
                self.buf_size[chromosome][strand] = b
            arr[i:i + n] = positions
            self.pointer[chromosome][strand] = i + n
        return

    @cython.ccall
    def finalize(self):
        """Shrink arrays and sort per-strand coordinates in place.
        
        Each chromosome's plus- and minus-strand arrays are resized to the observed
        counts, sorted ascending, and aggregate counters such as ``total`` and
        ``length`` are refreshed. Call this after loading data.
        """
        c: bytes
        size: cython.long
        big: list
        sizes: list
        kept: cython.ulong

        self.total = 0

        # each chromosome on its own: the large ones on threads
        # (map_chromosomes), any other here and now
        big = []
        sizes = []
        for c in self.get_chr_names():
            size = self.pointer[c][0] + self.pointer[c][1]
            if size >= _MIN_THREADED_SIZE:
                big.append(c)
                sizes.append(size)
            else:
                self.total += self._finalize_chrom(c)
        for kept in map_chromosomes(self._finalize_chrom, big, sizes):
            self.total += kept

        self.is_sorted = True
        self.length = self.fw * self.total
        return

    @cython.ccall
    def _finalize_chrom(self, c: bytes) -> cython.ulong:
        """finalize's work on chromosome `c`: shrink each strand's
        array to its count and sort it, and return the two sizes'
        sum for the caller to add to total. Touches nothing shared, so
        chromosomes can run on separate threads."""
        self.locations[c][0].resize(self.pointer[c][0], refcheck=False)
        self.locations[c][0].sort()
        self.locations[c][1].resize(self.pointer[c][1], refcheck=False)
        self.locations[c][1].sort()
        return self.locations[c][0].size + self.locations[c][1].size

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
    def get_locations_by_chr(self, chromosome: bytes):
        """Return the strand-specific arrays for a chromosome.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name, provided as bytes.
        
        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Pair of numpy arrays ``(plus, minus)`` containing 5' positions.
        
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
        """Return a sorted set of chromosome names stored in the track.
        
        Returns
        -------
        set
            Sorted chromosome names (bytes) that currently have positions.
        """
        return set(sorted(self.locations.keys()))

    @cython.ccall
    def sort(self):
        """Sort per-strand coordinate arrays for every chromosome.
        
        Positions are ordered ascending on each strand and the ``is_sorted`` flag is
        set to ``True`` once sorting completes.
        """
        c: bytes
        chrnames: set

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c][0].sort()
            self.locations[c][1].sort()

        self.is_sorted = True
        return

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.cfunc
    def _filter_dup_strand(self, k: bytes, strand: cython.int,
                           maxnum: cython.int,
                           kept_total: cython.pointer(cython.ulong)):
        """filter_dup for one strand of chromosome `k`.

        Returns locations[k][strand] itself when it holds at most one
        position. Otherwise returns an array of the positions kept,
        sets pointer[k][strand] to their number and adds it to
        kept_total[0]. A plain int32 array (filter_dup_plain) that owns
        its data, is held by nothing but locations[k], and holds at
        most pointer[k][strand] + 1 positions, the room the original's
        new array had, is filtered in place and shrunk, and returned;
        any other gets a new array and the old one is released.
        """
        p: cython.int
        n: cython.int
        current_loc: cython.int
        # index for old array, and index for new one
        i_old: cython.ulong
        i_new: cython.ulong
        size: cython.ulong
        kept: cython.Py_ssize_t
        locs: cnp.ndarray(cython.int, ndim=1)
        new_locs: cnp.ndarray(cython.int, ndim=1)

        i_new = 0
        locs = self.locations[k][strand]
        size = locs.shape[0]
        if len(locs) <= 1:
            return locs         # do nothing
        # in place: filter_dup_i4 allows src == dst, since each write
        # trails its read. Only when the references to locs are
        # locations[k]'s and this function's, so that no view or other
        # holder of the array sees it change
        if (size <= self.pointer[k][strand] + 1 and
                Py_REFCNT(locs) == 2 and owns_buffer(locs)):
            kept = filter_dup_plain(locs, locs, maxnum)
            if kept >= 0:
                locs.resize(kept, refcheck=False)
                kept_total[0] += kept
                self.pointer[k][strand] = kept
                return locs
        new_locs = np.zeros(self.pointer[k][strand] + 1, dtype='i4')
        # the loop below, in C, when both arrays are plain int32
        kept = filter_dup_plain(locs, new_locs, maxnum)
        if kept >= 0:
            i_new = kept
        else:
            new_locs[i_new] = locs[i_new]  # first item
            i_new += 1
            # the number of tags in the current location
            n = 1
            current_loc = locs[0]
            for i_old in range(1, size):
                p = locs[i_old]
                if p == current_loc:
                    n += 1
                else:
                    current_loc = p
                    n = 1
                if n <= maxnum:
                    new_locs[i_new] = p
                    i_new += 1
        new_locs.resize(i_new, refcheck=False)
        kept_total[0] += i_new
        self.pointer[k][strand] = i_new
        # free memory?
        # I know I should shrink it to 0 size directly,
        # however, on Mac OSX, it seems directly assigning 0
        # doesn't do a thing.
        locs.resize(self.buffer_size, refcheck=False)
        locs.resize(0, refcheck=False)
        # hope there would be no mem leak...
        return new_locs

    @cython.ccall
    def filter_dup(self, maxnum: cython.int = -1) -> cython.ulong:
        """Limit duplicate 5' positions to a maximum count per strand.
        
        Parameters
        ----------
        maxnum : int, optional
            Maximum number of occurrences allowed per coordinate. A negative value
            disables duplicate filtering.
        
        Returns
        -------
        int
            Total number of retained positions across both strands after filtering.
        
        Notes
        -----
        The track must be sorted before filtering. Coordinates exceeding ``maxnum``
        are discarded, pointers are updated, and ``total``/``length`` are recomputed.
        """
        k: bytes
        size: cython.long
        big: list
        sizes: list
        kept: cython.ulong

        if maxnum < 0:
            return self.total         # do nothing

        if not self.is_sorted:
            self.sort()

        self.total = 0
        self.length = 0

        # each chromosome on its own: the large ones on threads
        # (map_chromosomes), any other here and now
        big = []
        sizes = []
        for k in self.get_chr_names():
            size = len(self.locations[k][0]) + len(self.locations[k][1])
            if size >= _MIN_THREADED_SIZE:
                big.append(k)
                sizes.append(size)
            else:
                self.total += self._filter_dup_chrom(k, maxnum)
        for kept in map_chromosomes(self._filter_dup_chrom, big, sizes,
                                    maxnum):
            self.total += kept

        self.length = self.fw * self.total
        return self.total

    @cython.ccall
    def _filter_dup_chrom(self, k: bytes, maxnum: cython.int) -> cython.ulong:
        """filter_dup's work on chromosome `k`, + strand and then -
        strand: returns what it adds to total. Touches nothing shared
        but the existing key k of locations, so chromosomes can run on
        separate threads."""
        kept: cython.ulong = 0
        new_plus: cnp.ndarray(cython.int, ndim=1)
        new_minus: cnp.ndarray(cython.int, ndim=1)

        new_plus = self._filter_dup_strand(k, 0, maxnum, cython.address(kept))
        new_minus = self._filter_dup_strand(k, 1, maxnum, cython.address(kept))
        self.locations[k] = [new_plus, new_minus]
        return kept

    @cython.ccall
    def sample_percent(self, percent: cython.float, seed: cython.int = -1):
        """Down-sample positions in place by a fixed percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of positions to keep per strand between 0 and 1 (inclusive).
        seed : int, optional
            Seed forwarded to NumPy's RNG; a negative value uses global state.
        
        Notes
        -----
        Sampling is performed independently for plus and minus strands by shuffling
        each array, resizing to the requested fraction, and restoring sort order.
        Aggregate counters ``total`` and ``length`` are refreshed.
        """
        num: cython.int  # num: number of reads allowed on a certain chromosome
        k: bytes
        chrnames: set

        self.total = 0
        self.length = 0

        chrnames = self.get_chr_names()

        if seed >= 0:
            np.random.seed(seed)

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            num = cython.cast(cython.int,
                              round(self.locations[k][0].shape[0] * percent, 5))
            np.random.shuffle(self.locations[k][0])
            self.locations[k][0].resize(num, refcheck=False)
            self.locations[k][0].sort()
            self.pointer[k][0] = self.locations[k][0].shape[0]

            num = cython.cast(cython.int,
                              round(self.locations[k][1].shape[0] * percent, 5))
            np.random.shuffle(self.locations[k][1])
            self.locations[k][1].resize(num, refcheck=False)
            self.locations[k][1].sort()
            self.pointer[k][1] = self.locations[k][1].shape[0]

            self.total += self.pointer[k][0] + self.pointer[k][1]

        self.length = self.fw * self.total
        return

    @cython.ccall
    def sample_num(self, samplesize: cython.ulong, seed: cython.int = -1):
        """Down-sample positions in place so the total approximates ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target number of positions across both strands.
        seed : int, optional
            Seed forwarded to :meth:`sample_percent`.
        
        Notes
        -----
        The method converts ``samplesize`` into a sampling fraction using the current
        ``total`` and reuses :meth:`sample_percent`.
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize) / self.total
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def print_to_bed(self, fhd=None):
        """Stream the track as BED records.
        
        Parameters
        ----------
        fhd : io.IOBase, optional
            Writable file-like object. Defaults to ``sys.stdout``.
        
        Notes
        -----
        Emits one record per stored position with fixed-width intervals derived from
        ``fw`` and strand-specific orientation.
        """
        i: cython.int
        p: cython.int
        k: bytes
        chrnames: set

        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, io.IOBase)
        assert self.fw > 0, "FWTrack object .fw should be set larger than 0!"

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            plus = self.locations[k][0]

            for i in range(plus.shape[0]):
                p = plus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k.decode(),
                                                      p,
                                                      p + self.fw,
                                                      "+"))

            minus = self.locations[k][1]

            for i in range(minus.shape[0]):
                p = minus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k.decode(),
                                                      p-self.fw,
                                                      p,
                                                      "-"))
        return

    @cython.ccall
    def extract_region_tags(self, chromosome: bytes,
                            startpos: cython.int, endpos: cython.int) -> tuple:
        """Collect positions within a genomic window for both strands.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome identifier to query.
        startpos : int
            Inclusive start coordinate of the window.
        endpos : int
            Inclusive end coordinate of the window.
        
        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Pair of numpy arrays ``(plus, minus)`` containing positions inside the
            requested window.
        
        Notes
        -----
        The track is sorted on demand before performing the windowed lookup.
        """
        i: cython.int
        pos: cython.int
        rt_plus: np.ndarray(cython.int, ndim=1)
        rt_minus: np.ndarray(cython.int, ndim=1)
        temp: list
        chrnames: set

        if not self.is_sorted:
            self.sort()

        chrnames = self.get_chr_names()
        assert chromosome in chrnames, "chromosome %s can't be found in the FWTrack object." % chromosome

        (plus, minus) = self.locations[chromosome]

        temp = []
        for i in range(plus.shape[0]):
            pos = plus[i]
            if pos < startpos:
                continue
            elif pos > endpos:
                break
            else:
                temp.append(pos)
        rt_plus = np.array(temp)

        temp = []
        for i in range(minus.shape[0]):
            pos = minus[i]
            if pos < startpos:
                continue
            elif pos > endpos:
                break
            else:
                temp.append(pos)
        rt_minus = np.array(temp)
        return (rt_plus, rt_minus)

    @cython.ccall
    def compute_region_tags_from_peaks(self, peaks: PeakIO,
                                       func,
                                       window_size: cython.int = 100,
                                       cutoff: cython.float = 5.0) -> list:
        """Apply a summary function to tags collected around peak regions.
        
        Parameters
        ----------
        peaks : MACS3.IO.PeakIO.PeakIO
            Peak container providing genomic intervals and metadata.
        func : callable
            Callback invoked as ``func(chrom, plus, minus, startpos, endpos, ...)``
            for each peak. The callable must accept ``window_size`` and ``cutoff``
            keyword arguments.
        window_size : int, optional
            Half-window size added on each side of every peak when collecting tags.
        cutoff : float, optional
            Additional threshold passed to ``func``.
        
        Returns
        -------
        list
            Results returned by ``func`` for each processed peak.
        
        Notes
        -----
        Both the track and the ``peaks`` object are sorted before iteration, and
        per-chromosome state is reused to avoid rescanning arrays.
        """
        m: cython.int
        i: cython.int
        j: cython.int
        pos: cython.int
        startpos: cython.int
        endpos: cython.int

        plus: cnp.ndarray(cython.int, ndim=1)
        minus: cnp.ndarray(cython.int, ndim=1)
        rt_plus: cnp.ndarray(cython.int, ndim=1)
        rt_minus: cnp.ndarray(cython.int, ndim=1)

        chrom: bytes
        name: bytes

        temp: list
        retval: list
        pchrnames: set
        chrnames: set

        pchrnames = peaks.get_chr_names()
        retval = []

        # this object should be sorted
        if not self.is_sorted:
            self.sort()
        # PeakIO object should be sorted
        peaks.sort()

        chrnames = self.get_chr_names()

        for chrom in sorted(pchrnames):
            assert chrom in chrnames, "chromosome %s can't be found in the FWTrack object." % chrom
            (plus, minus) = self.locations[chrom]
            cpeaks = peaks.get_data_from_chrom(chrom)
            prev_i = 0
            prev_j = 0
            for m in range(len(cpeaks)):
                startpos = cpeaks[m]["start"] - window_size
                endpos = cpeaks[m]["end"] + window_size
                name = cpeaks[m]["name"]

                temp = []
                for i in range(prev_i, plus.shape[0]):
                    pos = plus[i]
                    if pos < startpos:
                        continue
                    elif pos > endpos:
                        prev_i = i
                        break
                    else:
                        temp.append(pos)
                rt_plus = np.array(temp, dtype="i4")

                temp = []
                for j in range(prev_j, minus.shape[0]):
                    pos = minus[j]
                    if pos < startpos:
                        continue
                    elif pos > endpos:
                        prev_j = j
                        break
                    else:
                        temp.append(pos)
                rt_minus = np.array(temp, dtype="i4")

                retval.append(func(chrom, rt_plus, rt_minus, startpos, endpos,
                                   name=name,
                                   window_size=window_size,
                                   cutoff=cutoff))
                # rewind window_size
                for i in range(prev_i, 0, -1):
                    if plus[prev_i] - plus[i] >= window_size:
                        break
                prev_i = i

                for j in range(prev_j, 0, -1):
                    if minus[prev_j] - minus[j] >= window_size:
                        break
                prev_j = j
                # end of a loop

        return retval

    @cython.ccall
    def pileup_a_chromosome(self, chrom: bytes,
                            d: cython.long,
                            scale_factor: cython.float = 1.0,
                            baseline_value: cython.float = 0.0,
                            directional: bool = True,
                            end_shift: cython.int = 0) -> list:
        """Compute a coverage pileup for a single chromosome.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        d : int
            Extension length applied in the 3' direction unless ``directional`` is
            ``False``.
        scale_factor : float, optional
            Value used to scale the resulting coverage.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        directional : bool, optional
            If ``False``, extend cuts symmetrically to both sides by ``d / 2``.
        end_shift : int, optional
            Shift applied to the 5' cuts before extension; positive values move
            toward the 3' direction.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` with numpy arrays describing
            the pileup breakpoints and scaled coverage.
        """
        five_shift: cython.long
        # adjustment to 5' end and 3' end positions to make a fragment
        three_shift: cython.long
        rlength: cython.long
        chrlengths: dict
        tmp_pileup: list

        chrlengths = self.get_rlengths()
        rlength = chrlengths[chrom]

        # adjust extension length according to 'directional' and
        # 'halfextension' setting.
        if directional:
            # only extend to 3' side
            five_shift = - end_shift
            three_shift = end_shift + d
        else:
            # both sides
            five_shift = d//2 - end_shift
            three_shift = end_shift + d - d//2

        tmp_pileup = pileup_from_PN_shifted(self.locations[chrom][0],
                                            self.locations[chrom][1],
                                            five_shift,
                                            three_shift,
                                            rlength,
                                            scale_factor,
                                            baseline_value)
        return tmp_pileup

    @cython.ccall
    def pileup_a_chromosome_c(self, chrom: bytes, ds: list,
                              scale_factor_s: list,
                              baseline_value: cython.float = 0.0,
                              directional: bool = True,
                              end_shift: cython.int = 0) -> list:
        """Compute a control pileup using multiple extension lengths.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        ds : list[int]
            Extension lengths used to build individual pileups.
        scale_factor_s : list[float]
            Scale factors paired with each entry in ``ds``.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        directional : bool, optional
            If ``False``, extend cuts symmetrically to both sides by ``d / 2``.
        end_shift : int, optional
            Shift applied to the 5' cuts before extension; positive values move
            toward the 3' direction.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` representing the merged pileup
            where the maximum value is taken across the supplied extensions.
        """
        d: cython.long
        five_shift: cython.long
        # adjustment to 5' end and 3' end positions to make a fragment
        three_shift: cython.long
        rlength: cython.long
        chrlengths: dict
        five_shift_s: list = []
        three_shift_s: list = []
        tmp_pileup: list
        prev_pileup: list

        chrlengths = self.get_rlengths()
        rlength = chrlengths[chrom]
        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        # adjust extension length according to 'directional' and
        # 'halfextension' setting.
        for d in ds:
            if directional:
                # only extend to 3' side
                five_shift_s.append(- end_shift)
                three_shift_s.append(end_shift + d)
            else:
                # both sides
                five_shift_s.append(d//2 - end_shift)
                three_shift_s.append(end_shift + d - d//2)

        # three windows (d, slocal and llocal): their pileups and the
        # maximum in one sweep, the same arrays as the loop below
        if len(ds) == 3:
            prev_pileup = se_all_in_one_pileup_max3(self.locations[chrom][0],
                                                    self.locations[chrom][1],
                                                    five_shift_s,
                                                    three_shift_s,
                                                    rlength,
                                                    scale_factor_s,
                                                    baseline_value)
            if prev_pileup is not None:
                return prev_pileup

        prev_pileup = None

        for i in range(len(ds)):
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]
            tmp_pileup = pileup_from_PN_shifted(self.locations[chrom][0],
                                                self.locations[chrom][1],
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


@cython.inline
@cython.cfunc
def left_sum(data,
             pos: cython.int,
             width: cython.int) -> cython.int:
    return sum([data[x] for x in data if x <= pos and x >= pos - width])


@cython.inline
@cython.cfunc
def right_sum(data,
              pos: cython.int,
              width: cython.int) -> cython.int:
    return sum([data[x] for x in data if x >= pos and x <= pos + width])


@cython.inline
@cython.cfunc
def left_forward(data,
                 pos: cython.int,
                 window_size: cython.int) -> cython.int:
    return data.get(pos, 0) - data.get(pos-window_size, 0)


@cython.inline
@cython.cfunc
def right_forward(data,
                  pos: cython.int,
                  window_size: cython.int) -> cython.int:
    return data.get(pos + window_size, 0) - data.get(pos, 0)
