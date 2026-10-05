# cython: language_level=3
# Time-stamp: <2025-02-12 18:14:41 Tao Liu>

"""Module Description: For pileup functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
# from MACS3.Utilities.Constants import *
import cython
from cython.cimports.MACS3.Signal.cPosValCalculation import single_end_pileup as c_single_end_pileup
from cython.cimports.MACS3.Signal.cPosValCalculation import write_pv_array_to_bedGraph as c_write_pv_array_to_bedGraph
from cython.cimports.MACS3.Signal.cPosValCalculation import PosVal
from cython.cimports.MACS3.Signal.cPosValCalculation import quick_pileup as c_quick_pileup

# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
import cython.cimports.numpy as cnp
from cython.cimports.cpython import bool

# ------------------------------------
# C lib
# ------------------------------------
from cython.cimports.libc.stdlib import free
from cython.cimports.libc.string import memcpy
from cython.cimports.libc.stdint import INT32_MAX, INT32_MIN
from cython.cimports.libc.math import INFINITY

# ------------------------------------
# utility internal functions
# ------------------------------------


@cython.cfunc
@cython.inline
def mean(a: float, b: float) -> float:
    return (a + b) / 2


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
def quiet_nan(x: cython.float) -> cython.float:
    """Return x as a trip through a Python float would: unchanged,
    except that a signalling NaN comes back quiet, as every IEEE 754
    float -> double conversion makes it. The C compiler folds a plain
    (float)(double)x to x, so the quiet bit is set by hand.
    """
    u: cython.uint
    memcpy(cython.address(u), cython.address(x), 4)
    # a NaN's bits, sign aside, are above those of inf
    u |= cython.cast(cython.uint, (u & 0x7fffffff) > 0x7f800000) << 22
    memcpy(cython.address(x), cython.address(u), 4)
    return x


@cython.cfunc
def clean_up_ndarray(x: cnp.ndarray):
    """ Clean up numpy array in two steps
    """
    i: cython.long

    i = x.shape[0] // 2
    x.resize(100000 if i > 100000 else i, refcheck=False)
    x.resize(0, refcheck=False)
    return


@cython.cfunc
def fix_coordinates(poss: cnp.ndarray, rlength: cython.int) -> cnp.ndarray:
    """Fix the coordinates.
    """
    i: cython.long
    ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int),
                                                  poss.data)  # pointer

    # fix those negative coordinates
    for i in range(poss.shape[0]):
        if ptr[i] < 0:
            ptr[i] = 0
        else:
            break

    # fix those over-boundary coordinates
    for i in range(poss.shape[0]-1, -1, -1):
        if ptr[i] > rlength:
            ptr[i] = rlength
        else:
            break
    return poss

# ------------------------------------
# functions
# ------------------------------------

# ------------------------------------
# functions for pileup_cmd only
# ------------------------------------

# This function uses pure C code for pileup


@cython.ccall
def pileup_and_write_se(trackI,
                        output_filename: bytes,
                        d: cython.int,
                        scale_factor: cython.float,
                        baseline_value: float = 0.0,
                        directional: bool = True,
                        halfextension: bool = True):
    """Pile up a single-end track and write a bedGraph file.

    This is a thin Cython wrapper that computes pileup using the
    C-accelerated routines in ``cPosValCalculation``.

    Parameters
    ----------
    trackI
        Single-end track object (FWTrackI).
    output_filename
        Output bedGraph path.
    d
        Fragment length estimate.
    scale_factor
        Scalar applied to pileup values.
    baseline_value
        Minimum output value per bin. Default is ``0.0``.
    directional
        If ``True``, extend reads only to 3' direction. If ``False``,
        extend to both sides. Default is ``True``.
    halfextension
        If ``True``, compute shift values from ``d`` using the half-
        extension scheme. Default is ``True``.

    Returns
    -------
    None

    Notes
    -----
    This function is currently only used by the ``macs3 pileup`` command.
    """
    five_shift: cython.long
    three_shift: cython.long
    i: cython.long
    rlength: cython.long
    l_data: cython.long = 0
    chroms: list
    n_chroms: cython.int
    chrom: bytes
    plus_tags: cnp.ndarray
    minus_tags: cnp.ndarray
    chrlengths: dict = trackI.get_rlengths()
    plus_tags_pos: cython.pointer(cython.int)
    minus_tags_pos: cython.pointer(cython.int)
    py_bytes: bytes
    chrom_char: cython.pointer(cython.char)
    _data: cython.pointer(PosVal)

    # This block should be reused to determine the actual shift values
    if directional:
        # only extend to 3' side
        if halfextension:
            five_shift = d//-4  # five shift is used to move cursor towards 5' direction to find the start of fragment
            three_shift = d*3//4  # three shift is used to move cursor towards 3' direction to find the end of fragment
        else:
            five_shift = 0
            three_shift = d
    else:
        # both sides
        if halfextension:
            five_shift = d//4
            three_shift = five_shift
        else:
            five_shift = d//2
            three_shift = d - five_shift
    # end of the block

    chroms = list(chrlengths.keys())
    n_chroms = len(chroms)

    fh = open(output_filename, "w")
    fh.write("")
    fh.close()

    for i in range(n_chroms):
        chrom = chroms[i]
        (plus_tags, minus_tags) = trackI.get_locations_by_chr(chrom)
        rlength = cython.cast(cython.long, chrlengths[chrom])
        plus_tags_pos = cython.cast(cython.p_int, plus_tags.data)
        minus_tags_pos = cython.cast(cython.p_int, minus_tags.data)

        _data = c_single_end_pileup(plus_tags_pos,
                                    plus_tags.shape[0],
                                    minus_tags_pos,
                                    minus_tags.shape[0],
                                    five_shift,
                                    three_shift,
                                    0,
                                    rlength,
                                    scale_factor,
                                    baseline_value,
                                    cython.address(l_data))

        # write
        py_bytes = chrom
        chrom_char = py_bytes
        c_write_pv_array_to_bedGraph(_data,
                                     l_data,
                                     chrom_char,
                                     output_filename,
                                     1)

        # clean
        free(_data)
    return

# function to pileup BAMPE/BEDPE stored in PETrackI object and write to a BEDGraph file
# this function uses c function


@cython.ccall
def pileup_and_write_pe(petrackI,
                        output_filename: bytes,
                        scale_factor: float = 1,
                        baseline_value: float = 0.0):
    """Pile up a paired-end track and write a bedGraph file.

    This is a thin Cython wrapper that computes pileup using the
    C-accelerated routines in ``cPosValCalculation``.

    Parameters
    ----------
    petrackI
        Paired-end track object (PETrackI). Must provide
        ``get_rlengths()`` and ``get_locations_by_chr(chrom)``.
    output_filename
        Output bedGraph path as ``bytes``.
    scale_factor
        Scalar applied to pileup values. Default is ``1``.
    baseline_value
        Minimum output value per bin. Default is ``0.0``.

    Returns
    -------
    None

    Notes
    -----
    This function is currently only used by the ``macs3 pileup`` command.

    Examples
    --------
    .. code-block:: python

        pileup_and_write_pe(
            petrackI,
            b"out.bedGraph",
            scale_factor=1.0,
            baseline_value=0.0,
        )
    """
    chrlengths: dict = petrackI.get_rlengths()
    chroms: list
    n_chroms: cython.int
    i: cython.long
    chrom: bytes
    locs: cnp.ndarray
    locs0: cnp.ndarray
    locs1: cnp.ndarray
    start_pos: cython.pointer(cython.int)
    end_pos: cython.pointer(cython.int)
    py_bytes: bytes
    chrom_char: cython.pointer(cython.char)
    _data: cython.pointer(PosVal)
    l_data: cython.long = 0

    chroms = list(chrlengths.keys())
    n_chroms = len(chroms)

    fh = open(output_filename, "w")
    fh.write("")
    fh.close()

    for i in range(n_chroms):
        chrom = chroms[i]
        locs = petrackI.get_locations_by_chr(chrom)

        locs0 = np.sort(locs['l'])
        locs1 = np.sort(locs['r'])
        start_pos = cython.cast(cython.p_int, locs0.data)  # <int *> locs0.data
        end_pos = cython.cast(cython.p_int, locs1.data)  # <int *> locs1.data

        _data = c_quick_pileup(start_pos,
                               end_pos,
                               locs0.shape[0],
                               scale_factor,
                               baseline_value,
                               cython.address(l_data))

        # write
        py_bytes = chrom
        chrom_char = py_bytes
        c_write_pv_array_to_bedGraph(_data,
                                     l_data,
                                     chrom_char,
                                     output_filename,
                                     1)

        # clean
        free(_data)
    return

# ------------------------------------
# functions for other codes
# ------------------------------------

# general pileup function implemented in cython


@cython.ccall
def se_all_in_one_pileup(plus_tags: cnp.ndarray,
                         minus_tags: cnp.ndarray,
                         five_shift: cython.long,
                         three_shift: cython.long,
                         rlength: cython.int,
                         scale_factor: cython.float,
                         baseline_value: cython.float) -> list:
    """Return pileup given 5' end of fragment at plus or minus strand
    separately, and given shift at both direction to recover a
    fragment. This function is for single end sequencing library
    only. Please directly use 'quick_pileup' function for Pair-end
    library.

    It contains a super-fast and simple algorithm proposed by Jie
    Wang. It will take sorted start positions and end positions, then
    compute pileup values.

    It will return a pileup result in similar structure as
    bedGraph. There are two python arrays:

    [end positions, values] or '[p,v] array' in other description for
    functions within MACS3.

    Two arrays have the same length and can be matched by index. End
    position at index x (p[x]) record continuous value of v[x] from
    p[x-1] to p[x].
    Parameters
    ----------
    plus_tags
        Sorted 5' end positions on the plus strand.
    minus_tags
        Sorted 5' end positions on the minus strand.
    five_shift
        Shift applied toward the 5' direction.
    three_shift
        Shift applied toward the 3' direction.
    rlength
        Chromosome length; coordinates are clipped to ``[0, rlength]``.
    scale_factor
        Scalar applied to pileup values.
    baseline_value
        Minimum output value per bin.

    Returns
    -------
    list
        ``[p, v]`` where ``p`` is an ``int32`` numpy array of end
        positions and ``v`` is a ``float32`` numpy array of values.

    Examples
    --------
    .. code-block:: python

        plus = np.array([10, 50, 100], dtype="i4")
        minus = np.array([30, 80], dtype="i4")
        p, v = se_all_in_one_pileup(
            plus,
        minus,
        five_shift=-25,
        three_shift=75,
        rlength=1000,
        scale_factor=1.0,
        baseline_value=0.0,
        )

    """
    p: cython.int
    pre_p: cython.int
    pileup: cython.int = 0

    i_s: cython.long = 0        # index of start_poss
    i_e: cython.long = 0        # index of end_poss
    i: cython.long
    I: cython.long = 0
    lx: cython.long
    w: cython.long

    start_poss: cnp.ndarray
    end_poss: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray

    # pointers are used for numpy arrays
    start_poss_ptr: cython.pointer(cython.int)
    end_poss_ptr: cython.pointer(cython.int)
    ret_p_ptr: cython.pointer(cython.int)
    ret_v_ptr: cython.pointer(cython.float)

    start_poss = np.concatenate((plus_tags-five_shift, minus_tags-three_shift))
    start_poss.sort()

    # A tag's end is its start plus w = five_shift + three_shift, so
    # the sorted ends are the sorted starts plus w: the same int32
    # values, in order unless one of them overflows.
    w = five_shift + three_shift
    lx = start_poss.shape[0]
    if (start_poss.dtype == np.int32 and
            -2147483648 <= w <= 2147483647 and
            (lx == 0 or
             (cython.cast(cython.long, start_poss[0]) + w >= -2147483648 and
              cython.cast(cython.long, start_poss[lx - 1]) + w <= 2147483647))):
        end_poss = start_poss + w
    else:
        end_poss = np.concatenate((plus_tags+three_shift, minus_tags+five_shift))
        end_poss.sort()

    # fix negative coordinations and those extends over end of chromosomes
    start_poss = fix_coordinates(start_poss, rlength)
    end_poss = fix_coordinates(end_poss, rlength)

    lx = start_poss.shape[0]

    start_poss_ptr = cython.cast(cython.pointer(cython.int),
                                 start_poss.data)  # <int32_t *> start_poss.data
    end_poss_ptr = cython.cast(cython.pointer(cython.int),
                               end_poss.data)  # <int32_t *> end_poss.data

    ret_p = np.zeros(2 * lx, dtype="i4")
    ret_v = np.zeros(2 * lx, dtype="f4")

    ret_p_ptr = cython.cast(cython.pointer(cython.int), ret_p.data)
    ret_v_ptr = cython.cast(cython.pointer(cython.float), ret_v.data)

    tmp = [ret_p, ret_v]        # for (endpos,value)

    if start_poss.shape[0] == 0:
        return tmp
    pre_p = min(start_poss_ptr[0], end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = max(0, baseline_value)
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1

    # pre_v = pileup

    assert start_poss.shape[0] == end_poss.shape[0]
    lx = start_poss.shape[0]

    while i_s < lx and i_e < lx:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup += 1
            i_s += 1
            start_poss_ptr += 1
        elif start_poss_ptr[0] > end_poss_ptr[0]:
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            i_e += 1
            end_poss_ptr += 1
        else:
            i_s += 1
            i_e += 1
            start_poss_ptr += 1
            end_poss_ptr += 1

    if i_e < lx:
        # add rest of end positions
        for i in range(i_e, lx):
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    # clean mem
    clean_up_ndarray(start_poss)
    clean_up_ndarray(end_poss)

    # resize
    ret_p.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)

    return tmp


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
def lower_bound_i4(t: cython.pointer(cython.int), n: cython.long,
                   x: cython.long) -> cython.long:
    """The first index of the sorted t[0:n] whose value is >= x."""
    lo: cython.long = 0
    hi: cython.long = n
    mid: cython.long

    while lo < hi:
        mid = (lo + hi) >> 1
        if t[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@cython.cfunc
@cython.inline
@cython.exceptval(check=False)
def clip_position(x: cython.long, rlength: cython.int) -> cython.long:
    """x clipped to [0, rlength], as fix_coordinates clips."""
    return min(max(x, 0), rlength)


@cython.cfunc
@cython.exceptval(check=False)
def last_breakpoint(t: cython.pointer(cython.int), n: cython.long,
                    five_shift: cython.long, three_shift: cython.long,
                    rlength: cython.int) -> cython.long:
    """The last position of the [p, v] array that se_all_in_one_pileup
    returns for starts t[0:n] - five_shift and ends t[0:n] +
    three_shift (t sorted, n > 0), each clipped to [0, rlength]; -1
    when that array is empty.

    It is the last position where the number of starts differs from
    the number of ends, or the first start position when no later
    one does and it is not 0.
    """
    i_s: cython.long = n - 1
    i_e: cython.long = n - 1
    first: cython.long = clip_position(t[0] - five_shift, rlength)
    p: cython.long
    s: cython.long
    e: cython.long

    while True:
        # the first start is the smallest position, and the last
        # one walked
        p = -1
        if i_e >= 0:
            p = clip_position(t[i_e] + three_shift, rlength)
        if i_s >= 0:
            p = max(p, clip_position(t[i_s] - five_shift, rlength))
        s = 0
        e = 0
        while i_s >= 0 and clip_position(t[i_s] - five_shift, rlength) == p:
            s += 1
            i_s -= 1
        while i_e >= 0 and clip_position(t[i_e] + three_shift, rlength) == p:
            e += 1
            i_e -= 1
        if p == first:
            return p if first != 0 else -1
        if s != e:
            return p


@cython.ccall
def se_all_in_one_pileup_max3(plus_tags: cnp.ndarray,
                              minus_tags: cnp.ndarray,
                              five_shift_s: list,
                              three_shift_s: list,
                              rlength: cython.int,
                              scale_factor_s: list,
                              baseline_value: cython.float) -> list:
    """The pileups of three windows and their maximum in one sweep.

    Returns the [p, v] array that this computes, bit for bit:

        prev = None
        for k in range(3):
            tmp = se_all_in_one_pileup(plus_tags, minus_tags,
                                       five_shift_s[k], three_shift_s[k],
                                       rlength, scale_factor_s[k],
                                       baseline_value)
            prev = over_two_pv_array(prev, tmp, func="max") if prev else tmp

    or None when the inputs are outside what the sweep covers (not
    three windows, not int32 tags, no tags, rlength < 1, a window of
    width five_shift + three_shift <= 0, a scale factor that is not
    finite and > 0, shifts or positions near the int32 limits); the
    caller then runs the loop above.

    A window's start positions are its tags' 5' ends minus
    five_shift on the plus strand and minus three_shift on the minus
    strand, so they are T - five_shift, where T is the sorted plus
    tags together with the minus tags minus (three_shift -
    five_shift); its ends are T + three_shift. Windows whose
    three_shift - five_shift agree share one T: one sort instead of
    two per window. The sweep walks the six position streams (starts
    and ends of three windows) together, clipped to [0, rlength] as
    fix_coordinates clips them; a window's pileup is the number of
    its starts walked minus the number of its ends.

    A window's own pileup has a breakpoint wherever its pileup
    changes (a start and an end at one position cancel), and at its
    first start position F when F != 0 even if it does not change
    there; its value on the segment that ends at p is
    max(pileup * scale_factor, baseline_value) with the pileup
    counted before p, which is max(0, baseline_value) for p <= F.
    Inside (0, rlength) the pileup does change at F, as no end can
    be there, and max(0 * scale_factor, baseline_value) has the bits
    of max(0, baseline_value) for a scale factor > 0. The maximum of
    the three has a breakpoint wherever one of them has, up to the
    first of their last breakpoints, and the value over_two_pv_array
    gives: the second window's value against the first's, then the
    third's against that (no value is NaN, so no quieting is due).
    """
    n_p: cython.long
    n_m: cython.long
    n: cython.long
    k: cython.int
    bound: cython.long = 0
    lim: cython.long = 1 << 28
    ts: list = [None, None, None]
    bases: dict = {}
    t_arr: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray
    out_p: cython.pointer(cython.int)
    out_v: cython.pointer(cython.float)

    # per window: base T, shifts, scale, first start position
    t0: cython.pointer(cython.int)
    t1: cython.pointer(cython.int)
    t2: cython.pointer(cython.int)
    f0: cython.long
    f1: cython.long
    f2: cython.long
    h0: cython.long
    h1: cython.long
    h2: cython.long
    sc0: cython.float
    sc1: cython.float
    sc2: cython.float
    first0: cython.long
    first1: cython.long
    first2: cython.long
    first_v: cython.float
    # per stream (a: starts, b: ends): the next position not walked
    pa0: cython.pointer(cython.int)
    pb0: cython.pointer(cython.int)
    pa1: cython.pointer(cython.int)
    pb1: cython.pointer(cython.int)
    pa2: cython.pointer(cython.int)
    pb2: cython.pointer(cython.int)
    m01: cython.long
    m23: cython.long
    m45: cython.long
    # pileups before the current position, and after it
    pu0: cython.int
    pu1: cython.int
    pu2: cython.int
    nu0: cython.int
    nu1: cython.int
    nu2: cython.int
    na0: cython.int
    nb0: cython.int
    na1: cython.int
    nb1: cython.int
    na2: cython.int
    nb2: cython.int
    w0: cython.float
    w1: cython.float
    w2: cython.float
    m: cython.float
    bp0: cython.int
    bp1: cython.int
    bp2: cython.int
    last: cython.long
    p: cython.long
    q: cython.long
    stop: cython.long
    I: cython.long = 0

    if (len(five_shift_s) != 3 or len(three_shift_s) != 3 or
            len(scale_factor_s) != 3):
        return None
    if (plus_tags is None or minus_tags is None or
            plus_tags.dtype != np.int32 or minus_tags.dtype != np.int32):
        return None
    n_p = plus_tags.shape[0]
    n_m = minus_tags.shape[0]
    n = n_p + n_m
    if n == 0 or rlength < 1:
        return None
    for k in range(3):
        if not (-lim <= five_shift_s[k] <= lim and
                -lim <= three_shift_s[k] <= lim and
                five_shift_s[k] + three_shift_s[k] > 0):
            return None
        bound = max(bound, abs(five_shift_s[k]), abs(three_shift_s[k]),
                    abs(three_shift_s[k] - five_shift_s[k]))
    sc0 = scale_factor_s[0]
    sc1 = scale_factor_s[1]
    sc2 = scale_factor_s[2]
    # finite and > 0 (a NaN fails every comparison)
    if not (0 < sc0 < INFINITY and 0 < sc1 < INFINITY and 0 < sc2 < INFINITY):
        return None

    # T for each distinct three_shift - five_shift, with a sentinel
    # after its n positions; every value T +/- a shift is exact in
    # int32 and below the sentinel's
    for k in range(3):
        delta = three_shift_s[k] - five_shift_s[k]
        if delta not in bases:
            t_arr = np.empty(n + 1, dtype=np.int32)
            t_arr[:n_p] = plus_tags
            np.subtract(minus_tags, delta, out=t_arr[n_p:n])
            t_arr[:n].sort()
            if (t_arr[0] < INT32_MIN + bound or
                    t_arr[n - 1] > INT32_MAX - 2 * bound - 1):
                return None
            t_arr[n] = INT32_MAX
            bases[delta] = t_arr
        ts[k] = bases[delta]

    t_arr = ts[0]
    t0 = cython.cast(cython.pointer(cython.int), t_arr.data)
    t_arr = ts[1]
    t1 = cython.cast(cython.pointer(cython.int), t_arr.data)
    t_arr = ts[2]
    t2 = cython.cast(cython.pointer(cython.int), t_arr.data)
    f0 = five_shift_s[0]
    f1 = five_shift_s[1]
    f2 = five_shift_s[2]
    h0 = three_shift_s[0]
    h1 = three_shift_s[1]
    h2 = three_shift_s[2]
    first_v = max(0, baseline_value)
    first0 = clip_position(t0[0] - f0, rlength)
    first1 = clip_position(t1[0] - f1, rlength)
    first2 = clip_position(t2[0] - f2, rlength)

    ret_p = np.empty(6 * n + 1, dtype="i4")
    ret_v = np.empty(6 * n + 1, dtype="f4")
    out_p = cython.cast(cython.pointer(cython.int), ret_p.data)
    out_v = cython.cast(cython.pointer(cython.float), ret_v.data)

    # position 0: every stream value <= 0 is clipped to it; no
    # window has a breakpoint at 0
    pa0 = t0 + lower_bound_i4(t0, n, f0 + 1)
    pb0 = t0 + lower_bound_i4(t0, n, 1 - h0)
    pa1 = t1 + lower_bound_i4(t1, n, f1 + 1)
    pb1 = t1 + lower_bound_i4(t1, n, 1 - h1)
    pa2 = t2 + lower_bound_i4(t2, n, f2 + 1)
    pb2 = t2 + lower_bound_i4(t2, n, 1 - h2)
    pu0 = cython.cast(cython.int, pa0 - pb0)
    pu1 = cython.cast(cython.int, pa1 - pb1)
    pu2 = cython.cast(cython.int, pa2 - pb2)

    # positions inside (0, rlength): no stream value there is
    # clipped, and the sentinel's values are >= stop
    stop = min(cython.cast(cython.long, rlength), INT32_MAX - bound)
    p = min(min(min(pa0[0] - f0, pb0[0] + h0), min(pa1[0] - f1, pb1[0] + h1)),
            min(pa2[0] - f2, pb2[0] + h2))
    while p < stop:
        # walk every stream's values at p, one per stream per round
        while True:
            pa0 += pa0[0] - f0 == p
            pb0 += pb0[0] + h0 == p
            pa1 += pa1[0] - f1 == p
            pb1 += pb1[0] + h1 == p
            pa2 += pa2[0] - f2 == p
            pb2 += pb2[0] + h2 == p
            # (the unsigned level keeps the C compiler from chaining
            # the five comparisons one after another; every value
            # here is > 0)
            m01 = min(pa0[0] - f0, pb0[0] + h0)
            m23 = min(pa1[0] - f1, pb1[0] + h1)
            m45 = min(pa2[0] - f2, pb2[0] + h2)
            q = cython.cast(cython.long,
                            min(min(cython.cast(cython.ulong, m01),
                                    cython.cast(cython.ulong, m23)),
                                cython.cast(cython.ulong, m45)))
            if q != p:
                break
        nu0 = cython.cast(cython.int, pa0 - pb0)
        nu1 = cython.cast(cython.int, pa1 - pb1)
        nu2 = cython.cast(cython.int, pa2 - pb2)
        w0 = max(pu0 * sc0, baseline_value)
        w1 = max(pu1 * sc1, baseline_value)
        w2 = max(pu2 * sc2, baseline_value)
        m = w1 if w1 > w0 else w0
        out_p[I] = cython.cast(cython.int, p)
        out_v[I] = w2 if w2 > m else m
        I += (nu0 != pu0) | (nu1 != pu1) | (nu2 != pu2)
        pu0 = nu0
        pu1 = nu1
        pu2 = nu2
        p = q

    # position rlength, to which every stream value at or above it is
    # clipped
    if rlength <= INT32_MAX - bound:
        p = rlength
        na0 = cython.cast(cython.int, n - lower_bound_i4(t0, n, rlength + f0))
        nb0 = cython.cast(cython.int, n - lower_bound_i4(t0, n, rlength - h0))
        na1 = cython.cast(cython.int, n - lower_bound_i4(t1, n, rlength + f1))
        nb1 = cython.cast(cython.int, n - lower_bound_i4(t1, n, rlength - h1))
        na2 = cython.cast(cython.int, n - lower_bound_i4(t2, n, rlength + f2))
        nb2 = cython.cast(cython.int, n - lower_bound_i4(t2, n, rlength - h2))
        if na0 or nb0 or na1 or nb1 or na2 or nb2:
            w0 = first_v if p <= first0 else max(pu0 * sc0, baseline_value)
            w1 = first_v if p <= first1 else max(pu1 * sc1, baseline_value)
            w2 = first_v if p <= first2 else max(pu2 * sc2, baseline_value)
            bp0 = (na0 != nb0) if p != first0 else (first0 != 0)
            bp1 = (na1 != nb1) if p != first1 else (first1 != 0)
            bp2 = (na2 != nb2) if p != first2 else (first2 != 0)
            m = w1 if w1 > w0 else w0
            out_p[I] = cython.cast(cython.int, p)
            out_v[I] = w2 if w2 > m else m
            I += bp0 | bp1 | bp2

    # the maximum ends at the first of the windows' last breakpoints
    last = min(min(last_breakpoint(t0, n, f0, h0, rlength),
                   last_breakpoint(t1, n, f1, h1, rlength)),
               last_breakpoint(t2, n, f2, h2, rlength))
    while I > 0 and out_p[I - 1] > last:
        I -= 1

    ret_p.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)
    return [ret_p, ret_v]

# quick pileup implemented in cython


@cython.ccall
def quick_pileup(start_poss: cnp.ndarray,
                 end_poss: cnp.ndarray,
                 scale_factor: cython.float,
                 baseline_value: cython.float) -> list:
    """Compute pileup from fragment start/end positions.

    The algorithm is a fast sweep over sorted start and end positions
    (Jie Wang). It returns a [p, v] array compatible with bedGraph
    semantics.

    Parameters
    ----------
    start_poss
        Sorted fragment start positions.
    end_poss
        Sorted fragment end positions.
    scale_factor
        Scalar applied to pileup values.
    baseline_value
        Minimum output value per bin.

    Returns
    -------
    list
        ``[p, v]`` where ``p`` is an ``int32`` numpy array of end
        positions and ``v`` is a ``float32`` numpy array of values.

    Examples
    --------
    .. code-block:: python
    
        starts = np.array([10, 50, 100], dtype="i4")
        ends = np.array([40, 90, 140], dtype="i4")
        p, v = quick_pileup(starts, ends, scale_factor=1.0, baseline_value=0.0)
    """
    p: cython.int
    pre_p: cython.int
    pileup: cython.int = 0

    i_s: cython.long = 0        # index of plus_tags
    i_e: cython.long = 0        # index of minus_tags
    i: cython.long
    I: cython.long = 0
    ls: cython.long = start_poss.shape[0]
    le: cython.long = end_poss.shape[0]
    l: cython.long = ls + le

    start_poss: cnp.ndarray
    end_poss: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray

    tmp: list

    # pointers are used for numpy arrays
    start_poss_ptr: cython.pointer(cython.int)
    end_poss_ptr: cython.pointer(cython.int)
    ret_p_ptr: cython.pointer(cython.int)
    ret_v_ptr: cython.pointer(cython.float)

    start_poss_ptr = cython.cast(cython.pointer(cython.int),
                                 start_poss.data)  # <int32_t *> start_poss.data
    end_poss_ptr = cython.cast(cython.pointer(cython.int),
                               end_poss.data)  # <int32_t *> end_poss.data

    ret_p = np.zeros(l, dtype="i4")
    ret_v = np.zeros(l, dtype="f4")

    ret_p_ptr = cython.cast(cython.pointer(cython.int), ret_p.data)
    ret_v_ptr = cython.cast(cython.pointer(cython.float), ret_v.data)

    tmp = [ret_p, ret_v]        # for (endpos,value)

    if ls == 0:
        return tmp
    pre_p = min(start_poss_ptr[0], end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = max(0, baseline_value)
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1

    # pre_v = pileup

    while i_s < ls and i_e < le:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup += 1
            # if pileup > max_pileup:
            #    max_pileup = pileup
            i_s += 1
            start_poss_ptr += 1
        elif start_poss_ptr[0] > end_poss_ptr[0]:
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            i_e += 1
            end_poss_ptr += 1
        else:
            i_s += 1
            i_e += 1
            start_poss_ptr += 1
            end_poss_ptr += 1

    if i_e < le:
        # add rest of end positions
        for i in range(i_e, le):
            p = end_poss_ptr[0]
            # for p in minus_tags[i_e:]:
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    ret_p.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)

    return tmp

# quick pileup implemented in cython


@cython.ccall
def naive_quick_pileup(sorted_poss: cnp.ndarray, extension: int) -> list:
    """Simple pileup by extending each tag symmetrically.

    Each input position is extended left and right by ``extension``.
    The input must be sorted; no sorting or validation is performed.

    Parameters
    ----------
    sorted_poss
        Sorted tag positions.
    extension
        Extension size applied to both sides.

    Returns
    -------
    list
        ``[p, v]`` where ``p`` is an ``int32`` numpy array of end
        positions and ``v`` is a ``float32`` numpy array of values.
    """
    p: cython.int
    pre_p: cython.int
    pileup: cython.int = 0

    i_s: cython.long = 0  # index of plus_tags
    i_e: cython.long = 0  # index of minus_tags
    i: cython.long
    I: cython.long = 0
    l: cython.long = sorted_poss.shape[0]

    start_poss: cnp.ndarray
    end_poss: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray

    # pointers are used for numpy arrays
    start_poss_ptr: cython.pointer(cython.int)
    end_poss_ptr: cython.pointer(cython.int)
    ret_p_ptr: cython.pointer(cython.int)
    ret_v_ptr: cython.pointer(cython.float)

    start_poss = sorted_poss - extension
    start_poss[start_poss < 0] = 0
    end_poss = sorted_poss + extension

    start_poss_ptr = cython.cast(cython.pointer(cython.int),
                                 start_poss.data)  # <int32_t *> start_poss.data
    end_poss_ptr = cython.cast(cython.pointer(cython.int),
                               end_poss.data)  # <int32_t *> end_poss.data

    ret_p = np.zeros(2*l, dtype="i4")
    ret_v = np.zeros(2*l, dtype="f4")

    ret_p_ptr = cython.cast(cython.pointer(cython.int), ret_p.data)
    ret_v_ptr = cython.cast(cython.pointer(cython.float), ret_v.data)

    if l == 0:
        raise Exception("length is 0")

    pre_p = min(start_poss_ptr[0], end_poss_ptr[0])

    if pre_p != 0:
        # the first chunk of 0
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = 0
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1

    # pre_v = pileup

    while i_s < l and i_e < l:
        if start_poss_ptr[0] < end_poss_ptr[0]:
            p = start_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = pileup
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup += 1
            i_s += 1
            start_poss_ptr += 1
        elif start_poss_ptr[0] > end_poss_ptr[0]:
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = pileup
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            i_e += 1
            end_poss_ptr += 1
        else:
            i_s += 1
            i_e += 1
            start_poss_ptr += 1
            end_poss_ptr += 1

    # add rest of end positions
    if i_e < l:
        for i in range(i_e, l):
            p = end_poss_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = pileup
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_poss_ptr += 1

    ret_p.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)

    return [ret_p, ret_v]

# general function to compare two pv arrays in cython.


@cython.ccall
def over_two_pv_array(pv_array1: list,
                      pv_array2: list,
                      func: str = "max") -> list:
    """Merge two [p, v] arrays with a pointwise reducer.

    Parameters
    ----------
    pv_array1
        First ``[p, v]`` array, same as output from quick_pileup function.
    pv_array2
        Second ``[p, v]`` array, same as output from quick_pileup function.
    func
        Reducer for overlapping regions. One of ``"max"``, ``"min"``,
        or ``"mean"``. Default is ``"max"``.

    Returns
    -------
    list
        Merged ``[p, v]`` array with the reducer applied to overlap
        regions.

    """
    # pre_p: cython.int

    l1: cython.long
    l2: cython.long
    i1: cython.long = 0
    i2: cython.long = 0
    I: cython.long = 0

    a1_pos: cnp.ndarray
    a2_pos: cnp.ndarray
    ret_pos: cnp.ndarray
    a1_v: cnp.ndarray
    a2_v: cnp.ndarray
    ret_v: cnp.ndarray

    # pointers are used for numpy arrays
    a1_pos_ptr: cython.pointer(cython.int)
    a2_pos_ptr: cython.pointer(cython.int)
    ret_pos_ptr: cython.pointer(cython.int)
    a1_v_ptr: cython.pointer(cython.float)
    a2_v_ptr: cython.pointer(cython.float)
    ret_v_ptr: cython.pointer(cython.float)

    # the reducer, applied in C: 0 is Python's max(v1, v2), 1 its
    # min(v1, v2), 2 mean(v1, v2)
    op: cython.int
    v1: cython.float
    v2: cython.float
    p1: cython.int
    p2: cython.int

    if func == "max":
        op = 0
    elif func == "min":
        op = 1
    elif func == "mean":
        op = 2
    else:
        raise Exception("Invalid function")

    [a1_pos, a1_v] = pv_array1
    [a2_pos, a2_v] = pv_array2
    ret_pos = np.zeros(a1_pos.shape[0] + a2_pos.shape[0], dtype="i4")
    ret_v = np.zeros(a1_pos.shape[0] + a2_pos.shape[0], dtype="f4")

    a1_pos_ptr = cython.cast(cython.pointer(cython.int), a1_pos.data)
    a1_v_ptr = cython.cast(cython.pointer(cython.float), a1_v.data)
    a2_pos_ptr = cython.cast(cython.pointer(cython.int), a2_pos.data)
    a2_v_ptr = cython.cast(cython.pointer(cython.float), a2_v.data)
    ret_pos_ptr = cython.cast(cython.pointer(cython.int), ret_pos.data)
    ret_v_ptr = cython.cast(cython.pointer(cython.float), ret_v.data)

    l1 = a1_pos.shape[0]
    l2 = a2_pos.shape[0]

    # pre_p = 0
    # remember the previous position in the new bedGraphTrackI object ret

    while i1 < l1 and i2 < l2:
        v1 = a1_v_ptr[i1]
        v2 = a2_v_ptr[i2]
        p1 = a1_pos_ptr[i1]
        p2 = a2_pos_ptr[i2]
        # Python's max(v1, v2) and min(v1, v2) return v2 only when it
        # compares greater (less) than v1, so v1 when the two are equal
        # (0.0 and -0.0) or either is NaN. These calls once took the
        # values through Python floats and back; quiet_nan keeps that.
        if op == 0:
            ret_v_ptr[I] = quiet_nan(v2 if v2 > v1 else v1)
        elif op == 1:
            ret_v_ptr[I] = quiet_nan(v2 if v2 < v1 else v1)
        else:
            ret_v_ptr[I] = mean(v1, v2)
        # the smaller position ends this segment; an array whose
        # position it is moves on, both when they are equal
        ret_pos_ptr[I] = p1 if p1 <= p2 else p2
        i1 += p1 <= p2
        i2 += p2 <= p1
        I += 1

    ret_pos.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)
    return [ret_pos, ret_v]


@cython.ccall
def naive_call_peaks(pv_array: list, min_v: cython.float,
                     max_v: cython.float = 1e30,
                     max_gap: cython.int = 50,
                     min_length: cython.int = 200):
    """Identify peak summits from a [p, v] array.

    Parameters
    ----------
    pv_array
        ``[p, v]`` array as produced by pileup functions.
    min_v
        Minimum value to be considered part of a peak.
    max_v
        Maximum allowed summit height.
        Default is ``1e30``.
    max_gap
        Maximum gap (in bp) allowed between adjacent peak segments to
        be merged. Default is ``50``.
    min_length
        Minimum peak length (in bp). Default is ``200``.

    Returns
    -------
    list
        List of ``(summit, height)`` tuples.

    """

    pre_p: cython.int
    p: cython.int
    i: cython.int
    x: cython.long           # index used for searching the first peak
    v: cython.double
    peak_content: list     # (pre_p, p, v) for each region in the peak
    ret_peaks: list = []   # returned peak summit and height

    peak_content = []
    (ps, vs) = pv_array
    if (isinstance(ps, np.ndarray) and isinstance(vs, np.ndarray) and
            ps.ndim == 1 and vs.ndim == 1 and
            ps.shape[0] == vs.shape[0] and
            ps.dtype == np.int32 and vs.dtype == np.float32 and
            ps.flags.c_contiguous and vs.flags.c_contiguous):
        # the arrays naive_quick_pileup returns
        return __naive_call_peaks_i4f4(ps, vs, min_v, max_v, max_gap,
                                       min_length)
    psn = iter(ps).__next__  # assign the next function to a viable to speed up
    vsn = iter(vs).__next__
    x = 0
    pre_p = 0                   # remember previous position
    while True:
        # find the first region above min_v
        try:         # try to read the first data range for this chrom
            p = psn()
            v = vsn()
        except Exception:
            break
        x += 1                  # index for the next point
        if v > min_v:
            peak_content = [(pre_p, p, v),]
            pre_p = p
            break               # found the first range above min_v
        else:
            pre_p = p

    for i in range(x, len(ps)):
        # continue scan the rest regions
        p = psn()
        v = vsn()
        if v <= min_v:          # not be detected as 'peak'
            pre_p = p
            continue
        # for points above min_v
        # if the gap is allowed
        # gap = pre_p - peak_content[-1][1] or the dist between pre_p and the last p
        if pre_p - peak_content[-1][1] <= max_gap:
            peak_content.append((pre_p, p, v))
        else:
            # when the gap is not allowed, close this peak IF length is larger than min_length
            if peak_content[-1][1] - peak_content[0][0] >= min_length:
                __close_peak(peak_content, ret_peaks, max_v, min_length)
            # reset and start a new peak
            peak_content = [(pre_p, p, v),]
        pre_p = p

    # save the last peak
    if peak_content:
        if peak_content[-1][1] - peak_content[0][0] >= min_length:
            __close_peak(peak_content, ret_peaks, max_v, min_length)
    return ret_peaks


@cython.cfunc
def __naive_call_peaks_i4f4(ps: cnp.ndarray, vs: cnp.ndarray,
                            min_v: cython.float,
                            max_v: cython.float,
                            max_gap: cython.int,
                            min_length: cython.int) -> list:
    """naive_call_peaks for C-contiguous int32 positions and float32
    values of equal length.

    The same scan, reading the arrays through pointers instead of
    one numpy scalar per point. Points above min_v go through the
    same Python-level code as in naive_call_peaks, so the peaks are
    the same.
    """
    n: cython.long = ps.shape[0]
    ps_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int),
                                                     ps.data)
    vs_ptr: cython.pointer(cython.float) = cython.cast(cython.pointer(cython.float),
                                                       vs.data)
    i: cython.long
    x: cython.long = 0
    pre_p: cython.int = 0
    p: cython.int
    v: cython.double
    peak_content: list = []
    ret_peaks: list = []

    # find the first region above min_v
    while x < n:
        p = ps_ptr[x]
        v = vs_ptr[x]
        x += 1
        if v > min_v:
            peak_content = [(pre_p, p, v),]
            pre_p = p
            break
        else:
            pre_p = p

    for i in range(x, n):
        # continue scan the rest regions
        p = ps_ptr[i]
        v = vs_ptr[i]
        if v <= min_v:          # not be detected as 'peak'
            pre_p = p
            continue
        if pre_p - peak_content[-1][1] <= max_gap:
            peak_content.append((pre_p, p, v))
        else:
            if peak_content[-1][1] - peak_content[0][0] >= min_length:
                __close_peak(peak_content, ret_peaks, max_v, min_length)
            peak_content = [(pre_p, p, v),]
        pre_p = p

    # save the last peak
    if peak_content:
        if peak_content[-1][1] - peak_content[0][0] >= min_length:
            __close_peak(peak_content, ret_peaks, max_v, min_length)
    return ret_peaks


@cython.cfunc
def __close_peak(peak_content,
                 peaks,
                 max_v: cython.float,
                 min_length: cython.int):
    """Internal function to find the summit and height

    If the height is larger than max_v, skip
    """
    tsummit: list = []
    summit: cython.int = 0
    summit_value: cython.float = 0
    tstart: cython.int
    tend: cython.int
    tvalue: cython.float

    for (tstart, tend, tvalue) in peak_content:
        if not summit_value or summit_value < tvalue:
            tsummit = [int((tend+tstart)/2),]
            summit_value = tvalue
        elif summit_value == tvalue:
            tsummit.append(int((tend+tstart)/2))
    summit = tsummit[int((len(tsummit)+1)/2)-1]
    if summit_value < max_v:
        peaks.append((summit, summit_value))
    return
