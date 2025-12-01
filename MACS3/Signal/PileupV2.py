# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-04-11 13:26:28 Tao Liu>

"""Module Description:

New pileup algorithm based on p-v array idea. It's modified from
original algorithm in MACS2 proposed by Jie Wang. Now we allow
different weights for different genomic ranges, and simplify the
approach. The basic idea is to remember for each position, how we
modify the 'pileup' value while walking on the chromosome.

For a genomic range r_i covering genomic positions from `s_i` to `e_i`
to be piled up, we assign the start position a value `w_i`, and the
end -`w_i`, so we will use a tuple of position and weight to remember
these operations: (`s_i`, `w_i`) and (`e_i`, `-w_i`). Then all N
ranges will be made into an array (2D) of position and weights as:

PV = [ (s_0, w_0), (e_0, -w_0), (s_1, w_1), (e_1, -w_1), ... (s_i, w_i),
       (e_i, -w_i), ..., (s_N, w_N), (e_N, -w_N) ]

Then the array PV will be sorted by the first dimension, aka the
position, no matter the position is from start or end positions from
ranges.

PV_sorted = [ (p_0, v_0), (p_1, v_1), ... , (p_i, v_i), ..., (p_{2N}, v_{2N}) ]

The pileup algorithm to produce a bedGraph style pileup (another p-v
array as in Pileup.py) can be simply described as:

set the initial pileup z as 0 or a given value, and a start position s
as 0, and an end position e as not-defined.

for i from 0 to 2N in PV_sorted:
    1: z = z + v_i
    2: e = p_i
    3: save the pileup from position s to e is z,
       in bedGraph style is to only save (e, z)
    4: s = e

The pileup result is in bedGraph style as in Pileup.py, where p is the
end of a region, and v is the pileup value. The result is stored in a NumPy array with dtype=[('p','i4'),('v','f4')]. Note that this is different from the Pilup.py where the pileup is stored in a list as [np.ndarray(dtype="i4"), np.ndarray(dtype="f4")].

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
from turtle import mode
from weakref import ref
import numpy as np
import cython
import cython.cimports.numpy as cnp

# ------------------------------------
# from cython.cimports.libc.stdlib import malloc, free, qsort

# ------------------------------------
# utility internal functions
# ------------------------------------

@cython.cfunc
@cython.inline
def mean(a: float, b: float) -> float:
    return (a + b) / 2


@cython.ccall
def mapping_function_always_1(L: cython.int, R: cython.int) -> cython.float:
    # always return 1, useful while the weight is already 1, or in
    # case of simply piling up fragments for coverage.
    return 1.0


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
def make_PV_from_LR(LR_array: cnp.ndarray,
                    mapping_func=mapping_function_always_1) -> cnp.ndarray:
    """Make sorted PV array from a LR array for certain chromosome in a
    PETrackI object. The V/weight will be assigned as
    `mapping_func( L, R )` or simply 1 if mapping_func is the default.

    LR array is an np.ndarray as with dtype
    [('l','i4'),('r','i4')] with length of N

    PV array is an np.ndarray with
    dtype=[('p','i4'),('v','f4')] with length of 2N
    """
    l_LR: cython.ulong
    l_PV: cython.ulong
    i: cython.ulong
    L: cython.int
    R: cython.int
    PV: cnp.ndarray

    l_LR = LR_array.shape[0]
    l_PV = 2 * l_LR
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    for i in range(l_LR):
        (L, R) = LR_array[i]
        PV[i*2] = (L, mapping_func(L, R))
        PV[i*2 + 1] = (R, -1.0 * mapping_func(L, R))
    PV.sort(order='p')
    return PV


@cython.cfunc
def make_PV_from_LRC(LRC_array: cnp.ndarray,
                     mapping_func=mapping_function_always_1) -> cnp.ndarray:
    """Make sorted PV array from a LR array for certain chromosome in a
    PETrackII object. The V/weight will be assigned as
    `mapping_func( L, R )` or simply 1 if mapping_func is the default.

    LRC array is an np.ndarray as with dtype
    [('l','i4'),('r','i4'),('c','u2')] with length of N

    PV array is an np.ndarray with
    dtype=[('p','i4'),('v','f4')] with length of 2N
    """
    l_LRC: cython.ulong
    l_PV: cython.ulong
    i: cython.ulong
    L: cython.int
    R: cython.int
    C: cython.ushort
    PV: cnp.ndarray

    l_LRC = LRC_array.shape[0]
    l_PV = 2 * l_LRC
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    for i in range(l_LRC):
        (L, R, C) = LRC_array[i]
        PV[i*2] = (L, C*mapping_func(L, R))
        PV[i*2 + 1] = (R, -1.0 * C * mapping_func(L, R))
    PV.sort(order='p')
    return PV


@cython.cfunc
def make_PV_from_PN(P_array: cnp.ndarray, N_array: cnp.ndarray,
                    extsize: cython.int) -> cnp.ndarray:
    """Make sorted PV array from two arrays for certain chromosome in
    a FWTrack object. P_array is for the 5' end positions in plus
    strand, and N_array is for minus strand. We don't support weight
    in this case since all positions should be extended with a fixed
    'extsize'.

    P_array or N_array is an np.ndarray with dtype='i4'

    PV array is an np.ndarray with
    dtype=[('p','i4'),('v','f4')] with length of 2N
    """
    l_PN: cython.ulong
    l_PV: cython.ulong
    i: cython.ulong
    L: cython.int
    R: cython.int
    PV: cnp.ndarray

    l_PN = P_array.shape[0]
    assert l_PN == N_array.shape[0]
    l_PV = 4 * l_PN
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    for i in range(l_PN):
        L = P_array[i]
        R = L + extsize
        PV[i*2] = (L, 1)
        PV[i*2 + 1] = (R, -1)
    for i in range(l_PN):
        R = N_array[i]
        L = R - extsize
        PV[(l_PN + i)*2] = (L, 1)
        PV[(l_PN + i)*2 + 1] = (R, -1)
    PV.sort(order='p')
    return PV


@cython.cfunc
def pileup_PV(PV_array: cnp.ndarray) -> cnp.ndarray:
    """The pileup algorithm to produce a bedGraph style pileup (another
    p-v array as in Pileup.py) can be simply described as:

    set the initial pileup z as 0 or a given value, and a start
    position s as 0, and an end position e as not-defined.

    for i from 0 to 2N in PV_sorted:
        z = z + v_i
        e = p_i
        save the pileup from position s to e is z --  in bedGraph style is to only save (e, z)
        s = e
    """
    z: cython.float
    e: cython.int
    s: cython.int
    c: cython.int
    v: cython.float
    i: cython.ulong
    pre_z: cython.float

    # this is in bedGraph style as in Pileup.pyx, p is the end of a
    # region, and v is the pileup value.
    pileup_PV: cnp.ndarray
    z = 0
    pre_z = -10000
    s = 0
    pileup_PV = np.zeros(shape=PV_array.shape[0], dtype=[('p', 'i4'),
                                                         ('v', 'f4')])
    c = 0
    for i in range(PV_array.shape[0]):
        e = PV_array[i]['p']
        v = PV_array[i]['v']
        # make sure only to record the final value for the same position
        if e != s:
            # merge the p-v pair with the previous pair if the same v is found
            if z == pre_z:
                pileup_PV[c-1]['p'] = e
            else:
                pileup_PV[c] = (e, z)
                c += 1
                pre_z = z
        z += v
        s = e
    pileup_PV.resize(c, refcheck=False)
    # assert z == 0
    return pileup_PV

# ------------------------------------
# public python functions
# ------------------------------------


@cython.ccall
def pileup_from_LR_hmmratac(LR_array: cnp.ndarray,
                            mapping_dict: dict) -> cnp.ndarray:
    # this function is specifically designed for piling up fragments
    # for `hmmratac`.
    #
    # As for `hmmratac`, the weight depends on the length of the
    # fragment, aka, value of R-L. Therefore, we need a mapping_dict
    # for mapping length to weight.
    l_LR: cython.ulong
    l_PV: cython.ulong
    i: cython.ulong
    L: cython.int
    R: cython.int
    PV: cnp.ndarray
    pileup: cnp.ndarray

    l_LR = LR_array.shape[0]
    l_PV = 2 * l_LR
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    for i in range(l_LR):
        (L, R) = LR_array[i]
        PV[i*2] = (L, mapping_dict[R - L])
        PV[i*2 + 1] = (R, -1 * mapping_dict[R - L])
    PV.sort(order='p')
    pileup = pileup_PV(PV)
    clean_up_ndarray(PV)
    return pileup


@cython.ccall
def pileup_from_LR(LR_array: cnp.ndarray,
                   mapping_func=mapping_function_always_1) -> cnp.ndarray:
    """This function will pile up the ndarray containing left and
    right positions, which is typically from PETrackI object. It's
    useful when generating the pileup of a single chromosome is
    needed.

    User needs to provide a numpy array of left and right positions,
    with dtype=[('l','i4'),('r','i4')]. User also needs to
    provide a mapping function to map the left and right position to
    certain weight.

    """
    PV_array: cnp.ndarray
    pileup: cnp.ndarray

    PV_array = make_PV_from_LR(LR_array, mapping_func=mapping_func)
    pileup = pileup_PV(PV_array)
    clean_up_ndarray(PV_array)
    return pileup


@cython.ccall
def pileup_from_LRC(LRC_array: cnp.ndarray,
                    mapping_func=mapping_function_always_1) -> cnp.ndarray:
    """This function will pile up the ndarray containing left and
    right positions and the counts, which is typically from PETrackII
    object which . It's useful when generating the pileup of a single
    chromosome is needed.

    User needs to provide a numpy array of left and right positions
    and the counts, with
    dtype=[('l','i4'),('r','i4'),('c','u2')]. User also needs to
    provide a mapping function to map the left and right position to
    certain weight.

    """
    PV_array: cnp.ndarray
    pileup: cnp.ndarray

    PV_array = make_PV_from_LRC(LRC_array, mapping_func=mapping_func)
    pileup = pileup_PV(PV_array)
    clean_up_ndarray(PV_array)
    return pileup


@cython.ccall
def pileup_from_PN(P_array: cnp.ndarray, N_array: cnp.ndarray,
                   extsize: cython.int) -> cnp.ndarray:
    """This function will pile up the ndarray containing plus strand
    (positive) and minus (negative) strand positions of all reads, which is
    typically from FWTrackI object. It's useful when generating the
    pileup of a single chromosome is needed.

    """
    PV_array: cnp.ndarray
    pileup: cnp.ndarray

    PV_array = make_PV_from_PN(P_array, N_array, extsize)
    pileup = pileup_PV(PV_array)
    clean_up_ndarray(PV_array)
    return pileup


@cython.cfunc
def _merge_sorted_arrays(arr1: cnp.ndarray, arr2: cnp.ndarray) -> cnp.ndarray:
    """Merge two individually sorted int32 arrays into one sorted array."""
    l1: cython.long = arr1.shape[0]
    l2: cython.long = arr2.shape[0]
    l: cython.long = l1 + l2
    i1: cython.long = 0
    i2: cython.long = 0
    I: cython.long = 0
    ret: cnp.ndarray = np.empty(l, dtype="i4")

    arr1_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), arr1.data)
    arr2_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), arr2.data)
    ret_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), ret.data)

    while i1 < l1 and i2 < l2:
        if arr1_ptr[0] <= arr2_ptr[0]:
            ret_ptr[0] = arr1_ptr[0]
            arr1_ptr += 1
            i1 += 1
        else:
            ret_ptr[0] = arr2_ptr[0]
            arr2_ptr += 1
            i2 += 1
        ret_ptr += 1
        I += 1

    while i1 < l1:
        ret_ptr[0] = arr1_ptr[0]
        ret_ptr += 1
        arr1_ptr += 1
        i1 += 1
        I += 1

    while i2 < l2:
        ret_ptr[0] = arr2_ptr[0]
        ret_ptr += 1
        arr2_ptr += 1
        i2 += 1
        I += 1

    ret.resize(I, refcheck=False)
    return ret


@cython.cfunc
def _pileup_sorted_bounds(start_poss: cnp.ndarray,
                          end_poss: cnp.ndarray,
                          scale_factor: cython.float,
                          baseline_value: cython.float) -> cnp.ndarray:
    """Pileup using already-sorted start/end bounds (two-pointer walk)."""
    p: cython.int
    pre_p: cython.int
    pileup: cython.int = 0
    i_s: cython.long = 0
    i_e: cython.long = 0
    I: cython.long = 0
    ls: cython.long = start_poss.shape[0]
    le: cython.long = end_poss.shape[0]
    l: cython.long = ls + le

    if ls == 0:
        return np.zeros(shape=0, dtype=[('p', 'i4'), ('v', 'f4')])

    start_view: cnp.ndarray = start_poss
    end_view: cnp.ndarray = end_poss
    start_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), start_view.data)
    end_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), end_view.data)

    ret_p = np.zeros(l, dtype="i4")
    ret_v = np.zeros(l, dtype="f4")
    ret_p_view: cnp.ndarray = ret_p
    ret_v_view: cnp.ndarray = ret_v
    ret_p_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), ret_p_view.data)
    ret_v_ptr: cython.pointer(cython.float) = cython.cast(cython.pointer(cython.float), ret_v_view.data)

    pre_p = start_ptr[0] if start_ptr[0] < end_ptr[0] else end_ptr[0]
    if pre_p != 0:
        ret_p_ptr[0] = pre_p
        ret_v_ptr[0] = max(0, baseline_value)
        ret_p_ptr += 1
        ret_v_ptr += 1
        I += 1

    while i_s < ls and i_e < le:
        if start_ptr[0] < end_ptr[0]:
            p = start_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup += 1
            i_s += 1
            start_ptr += 1
        elif start_ptr[0] > end_ptr[0]:
            p = end_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            i_e += 1
            end_ptr += 1
        else:
            i_s += 1
            i_e += 1
            start_ptr += 1
            end_ptr += 1

    if i_e < le:
        for _ in range(i_e, le):
            p = end_ptr[0]
            if p != pre_p:
                ret_p_ptr[0] = p
                ret_v_ptr[0] = max(pileup * scale_factor, baseline_value)
                ret_p_ptr += 1
                ret_v_ptr += 1
                I += 1
                pre_p = p
            pileup -= 1
            end_ptr += 1

    ret_p.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)
    ret = np.empty(I, dtype=[('p', 'i4'), ('v', 'f4')])
    ret['p'] = ret_p
    ret['v'] = ret_v
    return ret


@cython.ccall
def quick_pileup(start_poss: cnp.ndarray,
                 end_poss: cnp.ndarray,
                 scale_factor: cython.float,
                 baseline_value: cython.float) -> cnp.ndarray:
    """Return pileup given plus strand and minus strand positions of fragments.

    The output is a NumPy array with dtype=[('p','i4'),('v','f4')].
    """
    return _pileup_sorted_bounds(start_poss, end_poss, scale_factor, baseline_value)


@cython.ccall
def se_all_in_one_pileup(plus_tags: cnp.ndarray,
                         minus_tags: cnp.ndarray,
                         five_shift: cython.long,
                         three_shift: cython.long,
                         rlength: cython.int,
                         scale_factor: cython.float,
                         baseline_value: cython.float) -> cnp.ndarray:
    """Pileup single-end tags into a PV array with dtype=[('p','i4'),('v','f4')]."""
    start_poss: cnp.ndarray
    end_poss: cnp.ndarray
    pileup: cnp.ndarray
    start_plus = plus_tags - five_shift
    start_minus = minus_tags - three_shift
    np.clip(start_plus, 0, rlength, out=start_plus)
    np.clip(start_minus, 0, rlength, out=start_minus)
    end_plus = plus_tags + three_shift
    end_minus = minus_tags + five_shift
    np.clip(end_plus, 0, rlength, out=end_plus)
    np.clip(end_minus, 0, rlength, out=end_minus)

    start_poss = _merge_sorted_arrays(start_plus, start_minus)
    end_poss = _merge_sorted_arrays(end_plus, end_minus)

    pileup = _pileup_sorted_bounds(start_poss,
                                   end_poss,
                                   scale_factor,
                                   baseline_value)
    clean_up_ndarray(start_poss)
    clean_up_ndarray(end_poss)

    if pileup.shape[0] == 0:
        return pileup
    if scale_factor != 1:
        pileup['v'] *= scale_factor
    if baseline_value != 0:
        pileup['v'] = np.maximum(pileup['v'], baseline_value)
    return pileup


@cython.ccall
def naive_quick_pileup(sorted_poss: cnp.ndarray, extension: int) -> cnp.ndarray:
    """Simple pileup, every tag will be extended left and right with
    length `extension`.

    Note: Assumption is that `sorted_poss` has to be pre-sorted!
    """
    start_poss: cnp.ndarray
    end_poss: cnp.ndarray

    if sorted_poss.shape[0] == 0:
        raise Exception("length is 0")

    start_poss = sorted_poss - extension
    np.clip(start_poss, 0, None, out=start_poss)
    end_poss = sorted_poss + extension
    # arrays remain sorted after constant shift
    return _pileup_sorted_bounds(start_poss, end_poss, 1.0, 0.0)


@cython.ccall
def over_two_pv_array(pv_array1: cnp.ndarray,
                      pv_array2: cnp.ndarray,
                      func: str = "max") -> cnp.ndarray:
    """Merge two position-value arrays, keeping dtype=[('p','i4'),('v','f4')]."""
    l1: cython.long
    l2: cython.long
    i1: cython.long = 0
    i2: cython.long = 0
    I: cython.long = 0
    ret: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray

    if func == "max":
        call_func = max
    elif func == "min":
        call_func = min
    elif func == "mean":
        call_func = mean
    else:
        raise Exception("Invalid function")

    l1 = pv_array1.shape[0]
    l2 = pv_array2.shape[0]
    if l1 == 0:
        return pv_array2.copy()
    if l2 == 0:
        return pv_array1.copy()

    a1_p: cnp.ndarray = np.asarray(pv_array1['p'], dtype="i4", order="C")
    a2_p: cnp.ndarray = np.asarray(pv_array2['p'], dtype="i4", order="C")
    a1_v: cnp.ndarray = np.asarray(pv_array1['v'], dtype="f4", order="C")
    a2_v: cnp.ndarray = np.asarray(pv_array2['v'], dtype="f4", order="C")

    ret_p = np.empty(l1 + l2, dtype="i4")
    ret_v = np.empty(l1 + l2, dtype="f4")

    a1_p_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), a1_p.data)
    a2_p_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), a2_p.data)
    a1_v_ptr: cython.pointer(cython.float) = cython.cast(cython.pointer(cython.float), a1_v.data)
    a2_v_ptr: cython.pointer(cython.float) = cython.cast(cython.pointer(cython.float), a2_v.data)
    ret_p_ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int), ret_p.data)
    ret_v_ptr: cython.pointer(cython.float) = cython.cast(cython.pointer(cython.float), ret_v.data)

    while i1 < l1 and i2 < l2:
        p1 = a1_p_ptr[0]
        p2 = a2_p_ptr[0]
        ret_v_ptr[0] = call_func(a1_v_ptr[0], a2_v_ptr[0])
        if p1 < p2:
            ret_p_ptr[0] = p1
            i1 += 1
            a1_p_ptr += 1
            a1_v_ptr += 1
        elif p1 > p2:
            ret_p_ptr[0] = p2
            i2 += 1
            a2_p_ptr += 1
            a2_v_ptr += 1
        else:
            ret_p_ptr[0] = p1
            i1 += 1
            a1_p_ptr += 1
            a1_v_ptr += 1
            i2 += 1
            a2_p_ptr += 1
            a2_v_ptr += 1
        I += 1
        ret_p_ptr += 1
        ret_v_ptr += 1

    # while i1 < l1:
    #     ret_p[I] = a1_p[i1]
    #     ret_v[I] = a1_v[i1]
    #     i1 += 1
    #     I += 1
    #     a1_p_ptr += 1
    #     a1_v_ptr += 1
    #     ret_p_ptr += 1
    #     ret_v_ptr += 1

    # while i2 < l2:
    #     ret_p[I] = a2_p[i2]
    #     ret_v[I] = a2_v[i2]
    #     i2 += 1
    #     I += 1
    #     a2_p_ptr += 1
    #     a2_v_ptr += 1
    #     ret_p_ptr += 1
    #     ret_v_ptr += 1

    ret_p.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)
    ret = np.empty(I, dtype=[('p', 'i4'), ('v', 'f4')])
    ret['p'] = ret_p
    ret['v'] = ret_v
    return ret


@cython.cfunc
def __close_peak(peak_content,
                 peaks,
                 max_v: cython.float,
                 min_length: cython.int):
    """Internal function to find the summit and height."""
    tsummit: list = []
    summit: cython.int = 0
    summit_value: cython.float = 0
    tstart: cython.int
    tend: cython.int
    tvalue: cython.float

    for (tstart, tend, tvalue) in peak_content:
        if not summit_value or summit_value < tvalue:
            tsummit = [int((tend + tstart) / 2), ]
            summit_value = tvalue
        elif summit_value == tvalue:
            tsummit.append(int((tend + tstart) / 2))
    summit = tsummit[int((len(tsummit) + 1) / 2) - 1]
    if summit_value < max_v:
        peaks.append((summit, summit_value))
    return


@cython.ccall
def naive_call_peaks(pv_array: cnp.ndarray, min_v: cython.float,
                     max_v: cython.float = 1e30,
                     max_gap: cython.int = 50,
                     min_length: cython.int = 200):
    pre_p: cython.int
    p: cython.int
    i: cython.int
    x: cython.long
    v: cython.double
    peak_content: list
    ret_peaks: list = []

    peak_content = []
    psn = iter(pv_array['p']).__next__
    vsn = iter(pv_array['v']).__next__
    x = 0
    pre_p = 0
    while True:
        try:
            p = psn()
            v = vsn()
        except Exception:
            break
        x += 1
        if v > min_v:
            peak_content = [(pre_p, p, v), ]
            pre_p = p
            break
        else:
            pre_p = p

    for i in range(x, len(pv_array)):
        p = psn()
        v = vsn()
        if v <= min_v:
            pre_p = p
            continue
        if pre_p - peak_content[-1][1] <= max_gap:
            peak_content.append((pre_p, p, v))
        else:
            if peak_content[-1][1] - peak_content[0][0] >= min_length:
                __close_peak(peak_content, ret_peaks, max_v, min_length)
            peak_content = [(pre_p, p, v), ]
        pre_p = p

    if peak_content:
        if peak_content[-1][1] - peak_content[0][0] >= min_length:
            __close_peak(peak_content, ret_peaks, max_v, min_length)
    return ret_peaks


@cython.cfunc
def _write_pv_array_to_bedgraph(pv_array: cnp.ndarray,
                                chrom,
                                output_filename: bytes):
    chrom_str: str
    pre_p: cython.int = 0
    p: cython.int
    v: cython.float

    chrom_str = chrom.decode() if isinstance(chrom, (bytes, bytearray)) else str(chrom)
    with open(output_filename, "a") as fh:
        for (p, v) in pv_array:
            if p > pre_p:
                fh.write(f"{chrom_str}\t{pre_p}\t{int(p)}\t{float(v)}\n")
            pre_p = p
    return


@cython.ccall
def pileup_and_write_se(trackI,
                        output_filename: bytes,
                        d: cython.int,
                        scale_factor: cython.float,
                        baseline_value: float = 0.0,
                        directional: bool = True,
                        halfextension: bool = True):
    """Pileup a FWTrackI object and write the pileup result into a bedGraph file."""
    five_shift: cython.long
    three_shift: cython.long
    i: cython.long
    rlength: cython.long
    chroms: list
    n_chroms: cython.int
    chrom: bytes
    plus_tags: cnp.ndarray
    minus_tags: cnp.ndarray
    chrlengths: dict = trackI.get_rlengths()

    if directional:
        if halfextension:
            five_shift = d//-4
            three_shift = d*3//4
        else:
            five_shift = 0
            three_shift = d
    else:
        if halfextension:
            five_shift = d//4
            three_shift = five_shift
        else:
            five_shift = d//2
            three_shift = d - five_shift

    chroms = list(chrlengths.keys())
    n_chroms = len(chroms)

    fh = open(output_filename, "w")
    fh.write("")
    fh.close()

    for i in range(n_chroms):
        chrom = chroms[i]
        (plus_tags, minus_tags) = trackI.get_locations_by_chr(chrom)
        rlength = cython.cast(cython.long, chrlengths[chrom])
        pileup = se_all_in_one_pileup(plus_tags,
                                      minus_tags,
                                      five_shift,
                                      three_shift,
                                      rlength,
                                      scale_factor,
                                      baseline_value)
        _write_pv_array_to_bedgraph(pileup, chrom, output_filename)
    return


@cython.ccall
def pileup_and_write_pe(petrackI,
                        output_filename: bytes,
                        scale_factor: float = 1,
                        baseline_value: float = 0.0):
    """Pileup a PETrackI object and write the pileup result into a bedGraph file."""
    chrlengths: dict = petrackI.get_rlengths()
    chroms: list
    n_chroms: cython.int
    i: cython.long
    chrom: bytes
    locs: cnp.ndarray
    pileup: cnp.ndarray

    chroms = list(chrlengths.keys())
    n_chroms = len(chroms)

    fh = open(output_filename, "w")
    fh.write("")
    fh.close()

    for i in range(n_chroms):
        chrom = chroms[i]
        locs = petrackI.get_locations_by_chr(chrom)
        pileup = pileup_from_LR(locs)
        if pileup.shape[0] > 0:
            if scale_factor != 1:
                pileup['v'] *= scale_factor
            if baseline_value != 0:
                pileup['v'] = np.maximum(pileup['v'], baseline_value)
        _write_pv_array_to_bedgraph(pileup, chrom, output_filename)
    return
