# cython: language_level=3
# Time-stamp: <2024-10-15 10:20:32 Tao Liu>
"""Module Description: Build shifting model

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# Python modules
# ------------------------------------
import os
import threading

# ------------------------------------
# MACS3 modules
# ------------------------------------
# from MACS3.Utilities.Constants import *
from MACS3.Signal.PileupV2 import naive_quick_pileup, naive_call_peaks

# ------------------------------------
# Other modules
# ------------------------------------
import cython
from cython.cimports.cpython import bool
import numpy as np
import cython.cimports.numpy as cnp

# ------------------------------------
# C lib
# ------------------------------------
from cython.cimports.libc.stdlib import malloc, realloc, free
from cython.cimports.libc.string import memcpy

# ------------------------------------
# Peak finding on one strand in C
# ------------------------------------
# naive_call_peaks(naive_quick_pileup(tags, extension), min_v, max_v)
# for a C-contiguous int32 tag array, with the [p, v] points fed to the
# scan as they are made instead of being stored, and the peak content
# kept as C values instead of (pre_p, p, v) tuples. Every step and every
# C type is the one those two functions and __close_peak use, so the
# peaks are the same; it runs without the GIL, so strands can run on
# threads.

# at most this many threads find the peaks, and no more than the cores
# this process may run on
_MODEL_MAX_THREADS = cython.declare(cython.int, 8)

_PeakScan = cython.struct(
    min_v=cython.float,
    max_v=cython.float,
    max_gap=cython.int,
    min_length=cython.int,
    found=cython.bint,          # a point above min_v has been seen
    pre_p=cython.int,           # the previous point's p (0 before the first)
    first_start=cython.int,     # peak_content[0][0]
    last_end=cython.int,        # peak_content[-1][1]
    summit_value=cython.float,  # __close_peak's summit_value so far
    mids=cython.p_int,          # __close_peak's tsummit so far
    n_mids=cython.Py_ssize_t,
    cap_mids=cython.Py_ssize_t,
    summits=cython.p_int,       # the returned (summit, height) pairs
    heights=cython.p_float,
    n_peaks=cython.Py_ssize_t,
    cap_peaks=cython.Py_ssize_t,
    failed=cython.bint)         # an allocation failed


@cython.cfunc
@cython.inline
@cython.nogil
@cython.exceptval(check=False)
def _scan_push_mid(st: cython.pointer(_PeakScan), mid: cython.int) -> cython.void:
    new_cap: cython.Py_ssize_t
    new_mids: cython.p_int

    if st.n_mids == st.cap_mids:
        new_cap = 2 * st.cap_mids if st.cap_mids > 0 else 64
        new_mids = cython.cast(cython.p_int,
                               realloc(st.mids, new_cap * cython.sizeof(cython.int)))
        if new_mids == cython.NULL:
            st.failed = True
            return
        st.mids = new_mids
        st.cap_mids = new_cap
    st.mids[st.n_mids] = mid
    st.n_mids += 1


@cython.cfunc
@cython.inline
@cython.nogil
@cython.exceptval(check=False)
def _scan_add(st: cython.pointer(_PeakScan), tstart: cython.int,
              tend: cython.int, tvalue: cython.float) -> cython.void:
    """One pass of __close_peak's summit loop, for the entry
    (tstart, tend, tvalue) appended to the peak content."""
    # int((tend+tstart)/2): C int sum, true division as a double,
    # truncated toward zero
    mid: cython.int = cython.cast(cython.int,
                                  cython.cast(cython.double, tend + tstart) / 2.0)

    if not st.summit_value or st.summit_value < tvalue:
        st.n_mids = 0
        _scan_push_mid(st, mid)
        st.summit_value = tvalue
    elif st.summit_value == tvalue:
        _scan_push_mid(st, mid)


@cython.cfunc
@cython.inline
@cython.nogil
@cython.exceptval(check=False)
def _scan_start(st: cython.pointer(_PeakScan), tstart: cython.int,
                tend: cython.int, tvalue: cython.float) -> cython.void:
    """peak_content = [(tstart, tend, tvalue),]"""
    st.first_start = tstart
    st.last_end = tend
    st.summit_value = 0
    st.n_mids = 0
    _scan_add(st, tstart, tend, tvalue)


@cython.cfunc
@cython.nogil
@cython.exceptval(check=False)
def _scan_close(st: cython.pointer(_PeakScan)) -> cython.void:
    """The end of __close_peak: pick the summit and keep the peak if
    its height is below max_v."""
    summit: cython.int
    new_cap: cython.Py_ssize_t
    new_summits: cython.p_int
    new_heights: cython.p_float

    # tsummit[int((len(tsummit)+1)/2)-1], len(tsummit) >= 1
    summit = st.mids[(st.n_mids + 1) // 2 - 1]
    if st.summit_value < st.max_v:
        if st.n_peaks == st.cap_peaks:
            new_cap = 2 * st.cap_peaks if st.cap_peaks > 0 else 256
            new_summits = cython.cast(cython.p_int,
                                      realloc(st.summits, new_cap * cython.sizeof(cython.int)))
            if new_summits == cython.NULL:
                st.failed = True
                return
            st.summits = new_summits
            new_heights = cython.cast(cython.p_float,
                                      realloc(st.heights, new_cap * cython.sizeof(cython.float)))
            if new_heights == cython.NULL:
                st.failed = True
                return
            st.heights = new_heights
            st.cap_peaks = new_cap
        st.summits[st.n_peaks] = summit
        st.heights[st.n_peaks] = st.summit_value
        st.n_peaks += 1


@cython.cfunc
@cython.inline
@cython.nogil
@cython.exceptval(check=False)
def _scan_point(st: cython.pointer(_PeakScan), p: cython.int,
                vf: cython.float) -> cython.void:
    """naive_call_peaks' scan, for the next [p, v] point."""
    # naive_call_peaks reads the float32 value into a double
    v: cython.double = vf

    if not st.found:
        # looking for the first region above min_v
        if v > st.min_v:
            st.found = True
            _scan_start(st, st.pre_p, p, vf)
        st.pre_p = p
        return
    if v <= st.min_v:           # not be detected as 'peak'
        st.pre_p = p
        return
    # the gap and length tests were on Python ints, so no overflow
    if cython.cast(cython.longlong, st.pre_p) - st.last_end <= st.max_gap:
        st.last_end = p
        _scan_add(st, st.pre_p, p, vf)
    else:
        if cython.cast(cython.longlong, st.last_end) - st.first_start >= st.min_length:
            _scan_close(st)
        _scan_start(st, st.pre_p, p, vf)
    st.pre_p = p


@cython.cfunc
@cython.nogil
@cython.exceptval(check=False)
def _pileup_scan(tags: cython.p_int, l: cython.Py_ssize_t,
                 ext: cython.int, st: cython.pointer(_PeakScan)) -> cython.void:
    """naive_quick_pileup's merge of the start (tag - ext, clipped at 0)
    and end (tag + ext) positions, each [p, v] point passed to
    _scan_point, then the last peak closed. The int32 arithmetic wraps
    as numpy's does."""
    i_s: cython.Py_ssize_t = 0
    i_e: cython.Py_ssize_t = 0
    i: cython.Py_ssize_t
    pileup: cython.int = 0
    p: cython.int
    pre_p: cython.int
    sp: cython.int
    ep: cython.int
    uext: cython.uint = cython.cast(cython.uint, ext)

    sp = cython.cast(cython.int, cython.cast(cython.uint, tags[0]) - uext)
    if sp < 0:
        sp = 0
    ep = cython.cast(cython.int, cython.cast(cython.uint, tags[0]) + uext)
    pre_p = min(sp, ep)

    if pre_p != 0:
        # the first chunk of 0
        _scan_point(st, pre_p, 0)

    while i_s < l and i_e < l:
        sp = cython.cast(cython.int, cython.cast(cython.uint, tags[i_s]) - uext)
        if sp < 0:
            sp = 0
        ep = cython.cast(cython.int, cython.cast(cython.uint, tags[i_e]) + uext)
        if sp < ep:
            p = sp
            if p != pre_p:
                _scan_point(st, p, cython.cast(cython.float, pileup))
                pre_p = p
            pileup += 1
            i_s += 1
        elif sp > ep:
            p = ep
            if p != pre_p:
                _scan_point(st, p, cython.cast(cython.float, pileup))
                pre_p = p
            pileup -= 1
            i_e += 1
        else:
            i_s += 1
            i_e += 1

    # add rest of end positions
    for i in range(i_e, l):
        p = cython.cast(cython.int, cython.cast(cython.uint, tags[i]) + uext)
        if p != pre_p:
            _scan_point(st, p, cython.cast(cython.float, pileup))
            pre_p = p
        pileup -= 1

    # save the last peak
    if st.found:
        if cython.cast(cython.longlong, st.last_end) - st.first_start >= st.min_length:
            _scan_close(st)


@cython.cfunc
def _c_scan_takes(taglist, extension) -> cython.bint:
    """True when _naive_find_peaks_c gives this tag list's peaks: a 1-d
    C-contiguous native int32 array of at least 2 tags, and an extension
    numpy's int32 arithmetic takes."""
    if not isinstance(taglist, np.ndarray):
        return False
    return (taglist.ndim == 1 and taglist.dtype == np.int32 and
            taglist.flags.c_contiguous and taglist.shape[0] >= 2 and
            -2147483648 <= extension <= 2147483647)


@cython.cfunc
@cython.exceptval(check=False)
def _scan_init(st: cython.pointer(_PeakScan), min_v: cython.float,
               max_v: cython.float) -> cython.void:
    """The scan state before the first point, with the default max_gap
    and min_length."""
    st.min_v = min_v
    st.max_v = max_v
    st.max_gap = 50
    st.min_length = 200
    st.found = False
    st.pre_p = 0
    st.first_start = 0
    st.last_end = 0
    st.summit_value = 0
    st.mids = cython.NULL
    st.n_mids = 0
    st.cap_mids = 0
    st.summits = cython.NULL
    st.heights = cython.NULL
    st.n_peaks = 0
    st.cap_peaks = 0
    st.failed = False


@cython.ccall
def _naive_find_peaks_c(taglist: cnp.ndarray, extension: cython.int,
                        min_v: cython.float, max_v: cython.float) -> list:
    """naive_call_peaks(naive_quick_pileup(taglist, extension), min_v,
    max_v) with the default max_gap and min_length, for a tag list
    _c_scan_takes accepts. The scan runs without the GIL."""
    st: _PeakScan
    tags: cython.p_int = cython.cast(cython.p_int, taglist.data)
    l: cython.Py_ssize_t = taglist.shape[0]
    k: cython.Py_ssize_t
    ret: list

    _scan_init(cython.address(st), min_v, max_v)
    try:
        with cython.nogil:
            _pileup_scan(tags, l, extension, cython.address(st))
        if st.failed:
            raise MemoryError()
        ret = []
        for k in range(st.n_peaks):
            ret.append((st.summits[k], st.heights[k]))
    finally:
        free(st.mids)
        free(st.summits)
        free(st.heights)
    return ret


@cython.ccall
def _find_peaks_np(taglist: cnp.ndarray, extension: cython.int,
                   min_v: cython.float, max_v: cython.float) -> tuple:
    """The peaks _naive_find_peaks_c gives, as (summits, heights): an
    int32 and a float32 array holding the same C values its tuples
    are made of, so that no Python object is made per peak."""
    st: _PeakScan
    tags: cython.p_int = cython.cast(cython.p_int, taglist.data)
    l: cython.Py_ssize_t = taglist.shape[0]
    summits: cnp.ndarray
    heights: cnp.ndarray

    _scan_init(cython.address(st), min_v, max_v)
    try:
        with cython.nogil:
            _pileup_scan(tags, l, extension, cython.address(st))
        if st.failed:
            raise MemoryError()
        summits = np.empty(st.n_peaks, dtype=np.int32)
        heights = np.empty(st.n_peaks, dtype=np.float32)
        if st.n_peaks > 0:
            memcpy(summits.data, st.summits,
                   st.n_peaks * cython.sizeof(cython.int))
            memcpy(heights.data, st.heights,
                   st.n_peaks * cython.sizeof(cython.float))
    finally:
        free(st.mids)
        free(st.summits)
        free(st.heights)
    return (summits, heights)


def _peaks_as_list(peaks) -> list:
    """A strand's peaks as the [(summit, height)] list
    __naive_find_peaks returns: peaks itself when it is that list, or
    the list _naive_find_peaks_c would have made from the arrays of
    _find_peaks_np: Python ints from the C ints, Python floats from
    the C floats."""
    if isinstance(peaks, list):
        return peaks
    (summits, heights) = peaks
    return list(zip(summits.tolist(), heights.tolist()))


def _n_peaks(peaks) -> int:
    """The number of peaks in a list or an arrays pair."""
    if isinstance(peaks, list):
        return len(peaks)
    return peaks[0].shape[0]


def _peak_at(peaks, k: int) -> tuple:
    """Peak k as the (summit, height) tuple of the list form."""
    if isinstance(peaks, list):
        return peaks[k]
    return (int(peaks[0][k]), float(peaks[1][k]))


@cython.cfunc
@cython.nogil
@cython.exceptval(check=False)
def _pair_centers_c(ps: cython.p_int, ph: cython.p_float,
                    ip_max: cython.long,
                    ms: cython.p_int, mh: cython.p_float,
                    im_max: cython.long,
                    peaksize: cython.int,
                    out: cython.pointer(cython.p_int),
                    cap: cython.pointer(cython.Py_ssize_t)) -> cython.Py_ssize_t:
    """PeakModel.__find_pair_center on C arrays: the same loop, the
    same C types and the same expressions, the centers written to
    *out (grown with realloc, *cap its capacity). Returns the number of
    centers, -1 when a minus peak of height 0 meets the ratio test
    (where __find_pair_center raises ZeroDivisionError), or -2 when an
    allocation fails."""
    ip: cython.long = 0
    im: cython.long = 0
    im_prev: cython.long = 0
    flag_find_overlap: cython.bint = False
    pp: cython.int
    mp: cython.int
    pn: cython.float
    mn: cython.float
    n: cython.Py_ssize_t = 0
    new_cap: cython.Py_ssize_t
    new_out: cython.p_int

    while ip < ip_max and im < im_max:
        pp = ps[ip]
        pn = ph[ip]
        mp = ms[im]
        mn = mh[im]
        if pp-peaksize > mp:
            # move minus
            im += 1
        elif pp+peaksize < mp:
            # move plus
            ip += 1
            im = im_prev    # search minus peaks from previous index
            flag_find_overlap = False
        else:               # overlap!
            if not flag_find_overlap:
                flag_find_overlap = True
                # only the first index is recorded
                im_prev = im
            if mn == 0:
                return -1
            # number tags in plus and minus peak region are comparable...
            if pn/mn < 2 and pn/mn > 0.5:
                if pp < mp:
                    if n == cap[0]:
                        new_cap = 2 * cap[0] if cap[0] > 0 else 1024
                        new_out = cython.cast(cython.p_int,
                                              realloc(out[0], new_cap * cython.sizeof(cython.int)))
                        if new_out == cython.NULL:
                            return -2
                        out[0] = new_out
                        cap[0] = new_cap
                    out[0][n] = (pp+mp)//2
                    n += 1
            im += 1
    return n


@cython.cfunc
def _pair_centers_np(plus: tuple, minus: tuple, peaksize: cython.int):
    """The paired centers __find_pair_center gives for the peaks of two
    strands in the arrays form of _find_peaks_np, as an int32 array
    (the values its Python ints have), or None where it would raise."""
    ps_a: cnp.ndarray = plus[0]
    ph_a: cnp.ndarray = plus[1]
    ms_a: cnp.ndarray = minus[0]
    mh_a: cnp.ndarray = minus[1]
    ps: cython.p_int = cython.cast(cython.p_int, ps_a.data)
    ph: cython.p_float = cython.cast(cython.p_float, ph_a.data)
    ms: cython.p_int = cython.cast(cython.p_int, ms_a.data)
    mh: cython.p_float = cython.cast(cython.p_float, mh_a.data)
    ip_max: cython.long = ps_a.shape[0]
    im_max: cython.long = ms_a.shape[0]
    out: cython.p_int = cython.NULL
    cap: cython.Py_ssize_t = 0
    n: cython.Py_ssize_t
    centers: cnp.ndarray

    try:
        with cython.nogil:
            n = _pair_centers_c(ps, ph, ip_max, ms, mh, im_max, peaksize,
                                cython.address(out), cython.address(cap))
        if n == -2:
            raise MemoryError()
        if n == -1:
            return None
        centers = np.empty(n, dtype=np.int32)
        if n > 0:
            memcpy(centers.data, out, n * cython.sizeof(cython.int))
    finally:
        free(out)
    return centers


@cython.cfunc
@cython.nogil
@cython.exceptval(check=False)
def _add_line_c(pos1_ptr: cython.p_int, i1_max: cython.int,
                pos2_ptr: cython.p_int, i2_max: cython.int,
                start_ptr: cython.p_int, end_ptr: cython.p_int,
                max_index: cython.int, psize_adjusted1: cython.int,
                half_expansion: cython.double) -> cython.void:
    """The model's line projection, as PeakModel.__model_add_line
    did it: project each tag in pos2 within psize_adjusted1 of a center
    in pos1 onto start and end."""
    i1: cython.int = 0          # index for pos1
    i2: cython.int = 0          # index for pos2
    # index for pos2 in previous pos1
    # [pos1-self.peaksize,pos1+self.peaksize] region
    i2_prev: cython.int = 0
    p1: cython.int
    p2: cython.int
    s: cython.int
    e: cython.int
    flag_find_overlap: cython.bint = False

    while i1 < i1_max and i2 < i2_max:
        p1 = pos1_ptr[i1]
        p2 = pos2_ptr[i2]

        if p1-psize_adjusted1 > p2:
            # move pos2
            i2 += 1
        elif p1+psize_adjusted1 < p2:
            # move pos1
            i1 += 1
            i2 = i2_prev    # search minus peaks from previous index
            flag_find_overlap = False
        else:               # overlap!
            if not flag_find_overlap:
                flag_find_overlap = True
                # only the first index is recorded
                i2_prev = i2
            # project; p1-psize_adjusted1 <= p2 <= p1+psize_adjusted1
            # keeps 0 <= s and e <= max_index after the clamps
            s = cython.cast(cython.int, p2-half_expansion-p1+psize_adjusted1)
            if s < 0:
                s = 0
            start_ptr[s] += 1
            e = cython.cast(cython.int, p2+half_expansion-p1+psize_adjusted1)
            if e > max_index:
                e = max_index
            end_ptr[e] -= 1
            i2 += 1


@cython.cfunc
def _add_line_arrays(pos1_a: cnp.ndarray, pos2_a: cnp.ndarray,
                     start_a: cnp.ndarray, end_a: cnp.ndarray,
                     psize_adjusted1: cython.int,
                     half_expansion: cython.double):
    """_add_line_c on int32 C-contiguous arrays, without the GIL."""
    pos1_ptr: cython.p_int = cython.cast(cython.p_int, pos1_a.data)
    pos2_ptr: cython.p_int = cython.cast(cython.p_int, pos2_a.data)
    start_ptr: cython.p_int = cython.cast(cython.p_int, start_a.data)
    end_ptr: cython.p_int = cython.cast(cython.p_int, end_a.data)
    i1_max: cython.int = pos1_a.shape[0]
    i2_max: cython.int = pos2_a.shape[0]
    max_index: cython.int = start_a.shape[0] - 1

    with cython.nogil:
        _add_line_c(pos1_ptr, i1_max, pos2_ptr, i2_max, start_ptr, end_ptr,
                    max_index, psize_adjusted1, half_expansion)


def _model_threads() -> int:
    """min(_MODEL_MAX_THREADS, the cores this process may run on)."""
    if hasattr(os, "sched_getaffinity"):
        return min(_MODEL_MAX_THREADS, len(os.sched_getaffinity(0)))
    return min(_MODEL_MAX_THREADS, os.cpu_count() or 1)


def _run_on_threads(order: list, func, make_state=None) -> list:
    """func(task, state) for each task in order, taken in that order by
    up to _model_threads() threads (the calling one included), all
    joined before returning. Each thread's state is make_state() (None
    without it); returns the states of the threads that ran. func is
    expected to release the GIL for its work."""
    n_threads = min(_model_threads(), len(order))
    lock = threading.Lock()
    next_task = iter(order)
    errors = []
    states = []

    def work():
        try:
            state = make_state() if make_state is not None else None
        except BaseException as e:
            with lock:
                errors.append(e)
            return
        with lock:
            states.append(state)
        while True:
            with lock:
                if errors:
                    return
                task = next(next_task, None)
            if task is None:
                return
            try:
                func(task, state)
            except BaseException as e:
                with lock:
                    errors.append(e)
                return

    threads = [threading.Thread(target=work) for _ in range(n_threads - 1)]
    for t in threads:
        t.start()
    try:
        work()
    finally:
        for t in threads:
            t.join()
    if errors:
        raise errors[0]
    return states


def _find_peaks_on_threads(tasks: list, found: list, extension,
                           min_v, max_v):
    """found[k] = _find_peaks_np(tags, ...) for each (k, tags) in
    tasks, the largest first, on up to _model_threads() threads."""
    def scan(task, state):
        found[task[0]] = _find_peaks_np(task[1], extension, min_v, max_v)

    _run_on_threads(sorted(tasks, key=lambda t: -t[1].shape[0]), scan)


def _add_lines_on_threads(tasks: list, lines: list, psize_adjusted1,
                          half_expansion):
    """For each (centers, tags, k) in tasks, _add_line_c's
    projection of the tags around the centers onto the (start, end)
    arrays lines[k], the largest tag arrays first, on up to
    _model_threads() threads. Each thread projects onto zeroed arrays
    of its own, which are added to lines once all are joined: every
    projection adds 1 to a start and -1 to an end entry, so the sums
    are those of the serial loop, in any order."""
    def make_state():
        return [(np.zeros_like(start), np.zeros_like(end))
                for (start, end) in lines]

    def project(task, state):
        (start, end) = state[task[2]]
        _add_line_arrays(task[0], task[1], start, end, psize_adjusted1,
                         half_expansion)

    states = _run_on_threads(sorted(tasks, key=lambda t: -t[1].shape[0]),
                             project, make_state)
    for state in states:
        for k in range(len(lines)):
            np.add(lines[k][0], state[k][0], out=lines[k][0])
            np.add(lines[k][1], state[k][1], out=lines[k][1])


class NotEnoughPairsException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


@cython.cclass
class PeakModel:
    """Peak Model class.
    """
    # this can be PETrackI or FWTrack
    treatment: object
    # genome size
    gz: cython.double
    max_pairnum: cython.int
    umfold: cython.int
    lmfold: cython.int
    bw: cython.int
    d_min: cython.int
    tag_expansion_size: cython.int

    info: object
    debug: object
    warn: object
    error: object

    summary: str
    max_tags: cython.int
    peaksize: cython.int

    plus_line = cython.declare(cnp.ndarray, visibility="public")
    minus_line = cython.declare(cnp.ndarray, visibility="public")
    shifted_line = cython.declare(cnp.ndarray, visibility="public")
    xcorr = cython.declare(cnp.ndarray, visibility="public")
    ycorr = cython.declare(cnp.ndarray, visibility="public")

    d = cython.declare(cython.int, visibility="public")
    scan_window = cython.declare(cython.int, visibility="public")
    min_tags = cython.declare(cython.int, visibility="public")
    alternative_d = cython.declare(list, visibility="public")

    def __init__(self, opt, treatment, max_pairnum: cython.int = 500):
        # , double gz = 0, int umfold=30, int lmfold=10, int bw=200,
        # int ts = 25, int bg=0, bool quiet=False):
        self.treatment = treatment
        self.gz = opt.gsize
        self.umfold = opt.umfold
        self.lmfold = opt.lmfold
        # opt.tsize| test 10bps. The reason is that we want the best
        # 'lag' between left & right cutting sides. A tag will be
        # expanded to 10bps centered at cutting point.
        self.tag_expansion_size = 10
        # discard any predicted fragment sizes < d_min
        self.d_min = opt.d_min
        self.bw = opt.bw
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn
        self.error = opt.warn
        self.max_pairnum = max_pairnum

    @cython.ccall
    def build(self):
        """Build the model. Main function of PeakModel class.

        1. prepare self.d, self.scan_window, self.plus_line,
        self.minus_line and self.shifted_line.

        2. find paired + and - strand peaks

        3. find the best d using x-correlation
        """
        paired_peakpos: dict
        num_paired_peakpos: cython.long
        c: bytes                # chromosome

        self.peaksize = 2*self.bw
        # mininum unique hits on single strand, decided by lmfold
        self.min_tags = int(round(float(self.treatment.total) *
                                  self.lmfold *
                                  self.peaksize / self.gz / 2))
        # maximum unique hits on single strand, decided by umfold
        self.max_tags = int(round(float(self.treatment.total) *
                                  self.umfold *
                                  self.peaksize / self.gz / 2))
        self.debug(f"#2 min_tags: {self.min_tags}; max_tags:{self.max_tags}; ")
        self.info("#2 looking for paired plus/minus strand peaks...")
        # find paired + and - strand peaks
        paired_peakpos = self.__find_paired_peaks()

        num_paired_peakpos = 0
        for c in list(paired_peakpos.keys()):
            num_paired_peakpos += len(paired_peakpos[c])

        self.info("#2 Total number of paired peaks: %d" % (num_paired_peakpos))

        if num_paired_peakpos < 100:
            self.error(f"#2 MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only find {num_paired_peakpos}! Please make your MFOLD range broader and try again. If MACS3 still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead.")
            self.error("#2 Process for pairing-model is terminated!")
            raise NotEnoughPairsException("No enough pairs to build model")

        # build model, find the best d using cross-correlation
        self.__paired_peak_model(paired_peakpos)

    def __str__(self):
        """For debug...

        """
        return """
Summary of Peak Model:
  Baseline: %d
  Upperline: %d
  Fragment size: %d
  Scan window size: %d
""" % (self.min_tags, self.max_tags, self.d, self.scan_window)

    @cython.cfunc
    def __find_paired_peaks(self) -> dict:
        """Call paired peaks from fwtrackI object.

        Return paired peaks center positions.
        """
        i: cython.int
        chrs: list
        chrom: bytes
        plus_tags: cnp.ndarray(cython.int, ndim=1)
        minus_tags: cnp.ndarray(cython.int, ndim=1)
        plus_peaksinfo: object  # a list, or the arrays of _find_peaks_np
        minus_peaksinfo: object
        n_plus: cython.long
        n_minus: cython.long
        centers: object
        paired_peaks_pos: dict  # return
        found: list

        chrs = list(self.treatment.get_chr_names())
        chrs.sort()
        # the peaks of every strand the C scan takes, found up front
        found = self.__find_peaks_c(chrs)
        paired_peaks_pos = {}
        for i in range(len(chrs)):
            chrom = chrs[i]
            self.debug(f"Chromosome: {chrom}")
            # extract tag positions
            [plus_tags, minus_tags] = self.treatment.get_locations_by_chr(chrom)
            # look for + strand peaks
            plus_peaksinfo = found[2*i]
            if plus_peaksinfo is None:
                plus_peaksinfo = self.__naive_find_peaks(plus_tags)
            n_plus = _n_peaks(plus_peaksinfo)
            self.debug("Number of unique tags on + strand: %d" % (plus_tags.shape[0]))
            self.debug("Number of peaks in + strand: %d" % (n_plus))
            if n_plus:
                self.debug(f"plus peaks: first - {_peak_at(plus_peaksinfo, 0)} ... last - {_peak_at(plus_peaksinfo, n_plus-1)}")
            # look for - strand peaks
            minus_peaksinfo = found[2*i+1]
            if minus_peaksinfo is None:
                minus_peaksinfo = self.__naive_find_peaks(minus_tags)
            n_minus = _n_peaks(minus_peaksinfo)
            self.debug("Number of unique tags on - strand: %d" % (minus_tags.shape[0]))
            self.debug("Number of peaks in - strand: %d" % (n_minus))
            if n_minus:
                self.debug(f"minus peaks: first - {_peak_at(minus_peaksinfo, 0)} ... last - {_peak_at(minus_peaksinfo, n_minus-1)}")
            if not n_plus or not n_minus:
                self.debug("Chrom %s is discarded!" % (chrom))
                continue
            else:
                # both strands from the C scan: pair them in C
                centers = None
                if not isinstance(plus_peaksinfo, list) and not isinstance(minus_peaksinfo, list):
                    centers = _pair_centers_np(plus_peaksinfo, minus_peaksinfo, self.peaksize)
                if centers is None:
                    paired_peaks_pos[chrom] = self.__find_pair_center(_peaks_as_list(plus_peaksinfo),
                                                                      _peaks_as_list(minus_peaksinfo))
                else:
                    # __find_pair_center's messages
                    self.debug(f"ip_max: {n_plus}; im_max: {n_minus}")
                    if centers.shape[0]:
                        self.debug(f"Paired centers: first - {int(centers[0])} ... second - {int(centers[-1])} ")
                    paired_peaks_pos[chrom] = centers
                self.debug("Number of paired peaks in this chromosome: %d" % (len(paired_peaks_pos[chrom])))
        return paired_peaks_pos

    @cython.cfunc
    def __find_peaks_c(self, chrs: list) -> list:
        """__naive_find_peaks of both strands of every chromosome in
        chrs, at [2*i] (+) and [2*i+1] (-) for chrs[i], for the tag
        lists the C scan takes, on up to min(8, cores) threads, as the
        arrays of _find_peaks_np; None for the others.
        """
        i: cython.int
        found: list
        tasks: list

        extension = int(self.peaksize/2)
        found = [None] * (2 * len(chrs))
        tasks = []
        for i in range(len(chrs)):
            [plus_tags, minus_tags] = self.treatment.get_locations_by_chr(chrs[i])
            if _c_scan_takes(plus_tags, extension):
                tasks.append((2*i, plus_tags))
            if _c_scan_takes(minus_tags, extension):
                tasks.append((2*i+1, minus_tags))
        if tasks:
            _find_peaks_on_threads(tasks, found, extension,
                                   self.min_tags, self.max_tags)
        return found

    @cython.cfunc
    def __naive_find_peaks(self,
                           taglist: cnp.ndarray(cython.int, ndim=1)) -> list:
        """Naively call peaks based on tags counting.

        Return peak positions and the tag number in peak region by a tuple list[(pos,num)].
        """
        peak_info: list
        pileup_array: list

        # store peak pos in every peak region and unique tag number in
        # every peak region
        peak_info = []

        # less than 2 tags, no need to call peaks, return []
        if taglist.shape[0] < 2:
            return peak_info

        # build pileup by extending both side to half peak size
        pileup_array = naive_quick_pileup(taglist, int(self.peaksize/2))
        peak_info = naive_call_peaks(pileup_array,
                                     self.min_tags,
                                     self.max_tags)

        return peak_info

    @cython.cfunc
    def __paired_peak_model(self, paired_peakpos: dict):
        """Use paired peak positions and treatment tag positions to
        build the model.

        Modify self.(d, model_shift size and scan_window size. and
        extra, plus_line, minus_line and shifted_line for plotting).

        """
        window_size: cython.int
        i: cython.int
        chroms: list
        tasks: list
        centers: cnp.ndarray

        tags_plus: cnp.ndarray(cython.int, ndim=1)
        tags_minus: cnp.ndarray(cython.int, ndim=1)
        plus_start: cnp.ndarray(cython.int, ndim=1)
        plus_end: cnp.ndarray(cython.int, ndim=1)
        minus_start: cnp.ndarray(cython.int, ndim=1)
        minus_end: cnp.ndarray(cython.int, ndim=1)
        plus_line: cnp.ndarray(cython.int, ndim=1)
        minus_line: cnp.ndarray(cython.int, ndim=1)

        plus_data: cnp.ndarray
        minus_data: cnp.ndarray
        xcorr: cnp.ndarray
        ycorr: cnp.ndarray
        i_l_max: cnp.ndarray

        window_size = 1+2*self.peaksize+self.tag_expansion_size
        # for plus strand pileup
        self.plus_line = np.zeros(window_size, dtype="i4")
        # for minus strand pileup
        self.minus_line = np.zeros(window_size, dtype="i4")
        # for fast pileup
        plus_start = np.zeros(window_size, dtype="i4")
        # for fast pileup
        plus_end = np.zeros(window_size, dtype="i4")
        # for fast pileup
        minus_start = np.zeros(window_size, dtype="i4")
        # for fast pileup
        minus_end = np.zeros(window_size, dtype="i4")
        self.debug("start model_add_line...")
        chroms = list(paired_peakpos.keys())

        # every paired peak has plus line and minus line
        tasks = []
        for i in range(len(chroms)):
            # paired centers come from C ints, so int32 holds them exactly
            centers = np.array(paired_peakpos[chroms[i]], dtype="i4")
            (tags_plus, tags_minus) = self.treatment.get_locations_by_chr(chroms[i])
            tasks.append((centers, np.ascontiguousarray(tags_plus, dtype="i4"), 0))
            tasks.append((centers, np.ascontiguousarray(tags_minus, dtype="i4"), 1))
        # the half window, and half the expansion as a double
        _add_lines_on_threads(tasks,
                              [(plus_start, plus_end), (minus_start, minus_end)],
                              self.peaksize + self.tag_expansion_size // 2,
                              self.tag_expansion_size / 2)

        self.__count(plus_start, plus_end, self.plus_line)
        self.__count(minus_start, minus_end, self.minus_line)

        self.debug("start X-correlation...")
        # Now I use cross-correlation to find the best d
        plus_line = self.plus_line
        minus_line = self.minus_line

        # normalize first
        minus_data = (minus_line - minus_line.mean())/(minus_line.std()*len(minus_line))
        plus_data = (plus_line - plus_line.mean())/(plus_line.std()*len(plus_line))

        # cross-correlation
        ycorr = np.correlate(minus_data, plus_data, mode="full")[window_size-self.peaksize:window_size+self.peaksize]
        xcorr = np.linspace(len(ycorr)//2*-1, len(ycorr)//2, num=len(ycorr))

        # smooth correlation values to get rid of local maximums from small fluctuations.
        # window size is by default 11.
        ycorr = smooth(ycorr, window="flat")

        # all local maximums could be alternative ds.
        i_l_max = np.r_[False, ycorr[1:] > ycorr[:-1]] & np.r_[ycorr[:-1] > ycorr[1:], False]
        i_l_max = np.where(i_l_max)[0]
        i_l_max = i_l_max[xcorr[i_l_max] > self.d_min]
        i_l_max = i_l_max[np.argsort(ycorr[i_l_max])[::-1]]

        self.alternative_d = sorted([int(x) for x in xcorr[i_l_max]])
        assert len(self.alternative_d) > 0, "No proper d can be found! Tweak --mfold?"

        self.d = xcorr[i_l_max[0]]

        self.ycorr = ycorr
        self.xcorr = xcorr

        self.scan_window = max(self.d, self.tag_expansion_size)*2

        self.info("#2 Model building with cross-correlation: Done")

        return True

    @cython.cfunc
    def __count(self,
                start: cnp.ndarray(cython.int, ndim=1),
                end: cnp.ndarray(cython.int, ndim=1),
                line: cnp.ndarray(cython.int, ndim=1)):
        """
        """
        i: cython.int
        pileup: cython.long

        pileup = 0
        for i in range(line.shape[0]):
            pileup += start[i] + end[i]
            line[i] = pileup
        return

    @cython.cfunc
    def __find_pair_center(self,
                           pluspeaks: list,
                           minuspeaks: list):
        # index for plus peaks
        ip: cython.long = 0
        # index for minus peaks
        im: cython.long = 0
        # index for minus peaks in previous plus peak
        im_prev: cython.long = 0
        pair_centers: list
        ip_max: cython.long
        im_max: cython.long
        flag_find_overlap: bool
        pp: cython.int
        mp: cython.int
        pn: cython.float
        mn: cython.float

        pair_centers = []
        ip_max = len(pluspeaks)
        im_max = len(minuspeaks)
        self.debug(f"ip_max: {ip_max}; im_max: {im_max}")
        flag_find_overlap = False
        while ip < ip_max and im < im_max:
            # for (peakposition, tagnumber in peak)
            (pp, pn) = pluspeaks[ip]
            (mp, mn) = minuspeaks[im]
            if pp-self.peaksize > mp:
                # move minus
                im += 1
            elif pp+self.peaksize < mp:
                # move plus
                ip += 1
                im = im_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    # only the first index is recorded
                    im_prev = im
                # number tags in plus and minus peak region are comparable...
                if pn/mn < 2 and pn/mn > 0.5:
                    if pp < mp:
                        pair_centers.append((pp+mp)//2)
                im += 1
        if pair_centers:
            self.debug(f"Paired centers: first - {pair_centers[0]} ... second - {pair_centers[-1]} ")
        return pair_centers


# smooth function from SciPy cookbook:
# http://www.scipy.org/Cookbook/SignalSmooth
@cython.ccall
def smooth(x,
           window_len: cython.int = 11,
           window: str = 'hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be
                    an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming',
                'bartlett', 'blackman' flat window will produce a
                moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array
          instead of a string

    NOTE: length(output) != length(input), to correct this: return
          y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'flat':        # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[(window_len//2):-(window_len//2)]
