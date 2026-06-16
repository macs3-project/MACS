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

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
import numpy as np
import cython
import cython.cimports.numpy as cnp
from cython.cimports.libc.stdio import FILE, fopen, fclose, fprintf

# ------------------------------------
# from cython.cimports.libc.stdlib import malloc, free, qsort

# ------------------------------------
# utility internal functions
# ------------------------------------


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
def fix_coordinates(poss: cnp.ndarray, rlength: cython.int) -> cnp.ndarray:
    """Clip a sorted coordinate array to [0, rlength]."""
    i: cython.long
    ptr: cython.pointer(cython.int) = cython.cast(cython.pointer(cython.int),
                                                  poss.data)

    for i in range(poss.shape[0]):
        if ptr[i] < 0:
            ptr[i] = 0
        else:
            break

    for i in range(poss.shape[0]-1, -1, -1):
        if ptr[i] > rlength:
            ptr[i] = rlength
        else:
            break
    return poss


@cython.cfunc
def _pileup_sorted_unit(start_poss: cnp.ndarray,
                        end_poss: cnp.ndarray) -> cnp.ndarray:
    """Fast sweep for unit-weight sorted start/end positions."""
    p: cython.int
    pre_p: cython.int = 0
    z: cython.float = 0
    pre_z: cython.float = -10000
    i_s: cython.long = 0
    i_e: cython.long = 0
    c: cython.long = 0
    ls: cython.long = start_poss.shape[0]
    le: cython.long = end_poss.shape[0]
    ret: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray
    start_ptr: cython.pointer(cython.int)
    end_ptr: cython.pointer(cython.int)
    ret_p_ptr: cython.pointer(cython.int)
    ret_v_ptr: cython.pointer(cython.float)

    if ls + le == 0:
        return np.zeros(shape=0, dtype=[('p', 'i4'), ('v', 'f4')])

    ret_p = np.zeros(shape=ls + le, dtype="i4")
    ret_v = np.zeros(shape=ls + le, dtype="f4")
    start_ptr = cython.cast(cython.pointer(cython.int), start_poss.data)
    end_ptr = cython.cast(cython.pointer(cython.int), end_poss.data)
    ret_p_ptr = cython.cast(cython.pointer(cython.int), ret_p.data)
    ret_v_ptr = cython.cast(cython.pointer(cython.float), ret_v.data)

    while i_s < ls or i_e < le:
        if i_s < ls and (i_e >= le or start_ptr[i_s] < end_ptr[i_e]):
            p = start_ptr[i_s]
        elif i_e < le and (i_s >= ls or end_ptr[i_e] < start_ptr[i_s]):
            p = end_ptr[i_e]
        else:
            p = start_ptr[i_s]

        if p != pre_p:
            if z == pre_z:
                ret_p_ptr[c-1] = p
            else:
                ret_p_ptr[c] = p
                ret_v_ptr[c] = z
                c += 1
                pre_z = z
            pre_p = p

        while i_s < ls and start_ptr[i_s] == p:
            z += 1
            i_s += 1
        while i_e < le and end_ptr[i_e] == p:
            z -= 1
            i_e += 1

    ret_p.resize(c, refcheck=False)
    ret_v.resize(c, refcheck=False)
    ret = np.zeros(shape=c, dtype=[('p', 'i4'), ('v', 'f4')])
    ret['p'] = ret_p
    ret['v'] = ret_v
    clean_up_ndarray(ret_p)
    clean_up_ndarray(ret_v)
    return ret


@cython.cfunc
def _pileup_sorted_weighted(start_poss: cnp.ndarray,
                            start_weights: cnp.ndarray,
                            end_poss: cnp.ndarray,
                            end_weights: cnp.ndarray) -> cnp.ndarray:
    """Fast sweep for weighted sorted start/end positions."""
    p: cython.int
    pre_p: cython.int = 0
    z: cython.float = 0
    pre_z: cython.float = -10000
    i_s: cython.long = 0
    i_e: cython.long = 0
    c: cython.long = 0
    ls: cython.long = start_poss.shape[0]
    le: cython.long = end_poss.shape[0]
    ret: cnp.ndarray
    ret_p: cnp.ndarray
    ret_v: cnp.ndarray
    start_ptr: cython.pointer(cython.int)
    end_ptr: cython.pointer(cython.int)
    start_w_ptr: cython.pointer(cython.float)
    end_w_ptr: cython.pointer(cython.float)
    ret_p_ptr: cython.pointer(cython.int)
    ret_v_ptr: cython.pointer(cython.float)

    if ls + le == 0:
        return np.zeros(shape=0, dtype=[('p', 'i4'), ('v', 'f4')])

    ret_p = np.zeros(shape=ls + le, dtype="i4")
    ret_v = np.zeros(shape=ls + le, dtype="f4")
    start_ptr = cython.cast(cython.pointer(cython.int), start_poss.data)
    end_ptr = cython.cast(cython.pointer(cython.int), end_poss.data)
    start_w_ptr = cython.cast(cython.pointer(cython.float),
                              start_weights.data)
    end_w_ptr = cython.cast(cython.pointer(cython.float), end_weights.data)
    ret_p_ptr = cython.cast(cython.pointer(cython.int), ret_p.data)
    ret_v_ptr = cython.cast(cython.pointer(cython.float), ret_v.data)

    while i_s < ls or i_e < le:
        if i_s < ls and (i_e >= le or start_ptr[i_s] < end_ptr[i_e]):
            p = start_ptr[i_s]
        elif i_e < le and (i_s >= ls or end_ptr[i_e] < start_ptr[i_s]):
            p = end_ptr[i_e]
        else:
            p = start_ptr[i_s]

        if p != pre_p:
            if z == pre_z:
                ret_p_ptr[c-1] = p
            else:
                ret_p_ptr[c] = p
                ret_v_ptr[c] = z
                c += 1
                pre_z = z
            pre_p = p

        while i_s < ls and start_ptr[i_s] == p:
            z += start_w_ptr[i_s]
            i_s += 1
        while i_e < le and end_ptr[i_e] == p:
            z -= end_w_ptr[i_e]
            i_e += 1

    ret_p.resize(c, refcheck=False)
    ret_v.resize(c, refcheck=False)
    ret = np.zeros(shape=c, dtype=[('p', 'i4'), ('v', 'f4')])
    ret['p'] = ret_p
    ret['v'] = ret_v
    clean_up_ndarray(ret_p)
    clean_up_ndarray(ret_v)
    return ret


@cython.cfunc
def _pileup_from_PN_shifted(P_array: cnp.ndarray,
                            N_array: cnp.ndarray,
                            five_shift: cython.long,
                            three_shift: cython.long,
                            rlength: cython.int) -> cnp.ndarray:
    """Pile up plus/minus read ends using V1-compatible shifts."""
    start_poss: cnp.ndarray
    end_poss: cnp.ndarray

    start_poss = np.concatenate((P_array-five_shift, N_array-three_shift))
    end_poss = np.concatenate((P_array+three_shift, N_array+five_shift))
    start_poss.sort()
    end_poss.sort()
    start_poss = fix_coordinates(start_poss, rlength)
    end_poss = fix_coordinates(end_poss, rlength)
    return _pileup_sorted_unit(start_poss, end_poss)


@cython.cfunc
def _write_pv_to_bedGraph(pv: cnp.ndarray,
                          chrom: bytes,
                          output_filename: bytes,
                          scale_factor: cython.float,
                          baseline_value: cython.float):
    """Append one chromosome PV array to a bedGraph file."""
    fh: cython.pointer(FILE)
    p: cnp.ndarray
    v: cnp.ndarray
    p_ptr: cython.pointer(cython.int)
    v_ptr: cython.pointer(cython.float)
    i: cython.long
    l_pv: cython.long = pv.shape[0]
    pre: cython.int = 0
    pos: cython.int
    value: cython.float
    py_bytes: bytes = chrom
    chrom_char: cython.pointer(cython.char) = py_bytes

    fh = fopen(output_filename, b"a")
    if fh == cython.NULL:
        raise OSError("Unable to open bedGraph output file")

    p = np.ascontiguousarray(pv['p'])
    v = np.ascontiguousarray(pv['v'])
    p_ptr = cython.cast(cython.pointer(cython.int), p.data)
    v_ptr = cython.cast(cython.pointer(cython.float), v.data)

    for i in range(l_pv):
        pos = p_ptr[i]
        value = v_ptr[i] * scale_factor
        if value < baseline_value:
            value = baseline_value
        fprintf(fh, b"%s\t%d\t%d\t%.5f\n", chrom_char, pre, pos, value)
        pre = pos
    fclose(fh)
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
    weight: cython.float
    PV: cnp.ndarray
    LR_l: cnp.ndarray
    LR_r: cnp.ndarray
    PV_p: cnp.ndarray
    PV_v: cnp.ndarray

    l_LR = LR_array.shape[0]
    l_PV = 2 * l_LR
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    LR_l = LR_array['l']
    LR_r = LR_array['r']
    PV_p = PV['p']
    PV_v = PV['v']
    for i in range(l_LR):
        L = LR_l[i]
        R = LR_r[i]
        weight = mapping_func(L, R)
        PV_p[i*2] = L
        PV_v[i*2] = weight
        PV_p[i*2 + 1] = R
        PV_v[i*2 + 1] = -1.0 * weight
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
    weight: cython.float
    PV: cnp.ndarray
    LRC_l: cnp.ndarray
    LRC_r: cnp.ndarray
    LRC_c: cnp.ndarray
    PV_p: cnp.ndarray
    PV_v: cnp.ndarray

    l_LRC = LRC_array.shape[0]
    l_PV = 2 * l_LRC
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    LRC_l = LRC_array['l']
    LRC_r = LRC_array['r']
    LRC_c = LRC_array['c']
    PV_p = PV['p']
    PV_v = PV['v']
    for i in range(l_LRC):
        L = LRC_l[i]
        R = LRC_r[i]
        C = LRC_c[i]
        weight = C * mapping_func(L, R)
        PV_p[i*2] = L
        PV_v[i*2] = weight
        PV_p[i*2 + 1] = R
        PV_v[i*2 + 1] = -1.0 * weight
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
    PV_p: cnp.ndarray
    PV_v: cnp.ndarray

    l_PN = P_array.shape[0]
    assert l_PN == N_array.shape[0]
    l_PV = 4 * l_PN
    PV = np.zeros(shape=l_PV, dtype=[('p', 'i4'), ('v', 'f4')])
    PV_p = PV['p']
    PV_v = PV['v']
    for i in range(l_PN):
        L = P_array[i]
        R = L + extsize
        PV_p[i*2] = L
        PV_v[i*2] = 1
        PV_p[i*2 + 1] = R
        PV_v[i*2 + 1] = -1
    for i in range(l_PN):
        R = N_array[i]
        L = R - extsize
        PV_p[(l_PN + i)*2] = L
        PV_v[(l_PN + i)*2] = 1
        PV_p[(l_PN + i)*2 + 1] = R
        PV_v[(l_PN + i)*2 + 1] = -1
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
    """Pile up LR fragments using HMMRATAC
    """
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
    start_poss: cnp.ndarray
    end_poss: cnp.ndarray

    if mapping_func is mapping_function_always_1:
        start_poss = np.sort(LR_array['l'])
        end_poss = np.sort(LR_array['r'])
        pileup = _pileup_sorted_unit(start_poss, end_poss)
        clean_up_ndarray(start_poss)
        clean_up_ndarray(end_poss)
        return pileup

    PV_array = make_PV_from_LR(LR_array, mapping_func=mapping_func)
    pileup = pileup_PV(PV_array)
    clean_up_ndarray(PV_array)
    return pileup


@cython.ccall
def pileup_from_LRC(LRC_array: cnp.ndarray,
                    mapping_func=mapping_function_always_1) -> cnp.ndarray:
    """This function will pile up the ndarray containing left and
    right positions and the counts, which is typically from PETrackII
    object. It's useful when generating the pileup of a single
    chromosome is needed.

    User needs to provide a numpy array of left and right positions
    and the counts, with
    dtype=[('l','i4'),('r','i4'),('c','u2')]. User also needs to
    provide a mapping function to map the left and right position to
    certain weight.

    """
    PV_array: cnp.ndarray
    pileup: cnp.ndarray
    start_poss: cnp.ndarray
    end_poss: cnp.ndarray
    start_weights: cnp.ndarray
    end_weights: cnp.ndarray
    indices: cnp.ndarray
    start_indices: cnp.ndarray

    if mapping_func is mapping_function_always_1:
        start_indices = np.argsort(LRC_array['l'])
        start_poss = LRC_array['l'][start_indices]
        start_weights = LRC_array['c'][start_indices].astype("f4", copy=False)
        indices = np.argsort(LRC_array['r'])
        end_poss = LRC_array['r'][indices]
        end_weights = LRC_array['c'][indices].astype("f4", copy=False)
        pileup = _pileup_sorted_weighted(start_poss,
                                         start_weights,
                                         end_poss,
                                         end_weights)
        clean_up_ndarray(start_poss)
        clean_up_ndarray(start_weights)
        clean_up_ndarray(end_poss)
        clean_up_ndarray(end_weights)
        return pileup

    PV_array = make_PV_from_LRC(LRC_array, mapping_func=mapping_func)
    pileup = pileup_PV(PV_array)
    clean_up_ndarray(PV_array)
    return pileup


@cython.ccall
def pileup_from_PN(P_array: cnp.ndarray, N_array: cnp.ndarray,
                   extsize: cython.int) -> cnp.ndarray:
    """This function will pile up the ndarray containing plus
    (positive) and minus (negative) positions of all reads, which is
    typically from FWTrackI object. It's useful when generating the
    pileup of a single chromosome is needed.

    """
    pileup: cnp.ndarray
    start_poss: cnp.ndarray
    end_poss: cnp.ndarray

    start_poss = np.concatenate((P_array, N_array-extsize))
    end_poss = np.concatenate((P_array+extsize, N_array))
    start_poss.sort()
    end_poss.sort()
    pileup = _pileup_sorted_unit(start_poss, end_poss)
    clean_up_ndarray(start_poss)
    clean_up_ndarray(end_poss)
    return pileup


@cython.ccall
def pileup_and_write_se(trackI,
                        output_filename: bytes,
                        d: cython.int,
                        scale_factor: cython.float,
                        baseline_value: cython.float = 0.0,
                        directional: cython.bint = True,
                        halfextension: cython.bint = True):
    """Pile up a single-end track and write a bedGraph file."""
    five_shift: cython.long
    three_shift: cython.long
    rlength: cython.long
    chrom: bytes
    plus_tags: cnp.ndarray
    minus_tags: cnp.ndarray
    chrlengths: dict = trackI.get_rlengths()
    pv: cnp.ndarray
    fh: cython.pointer(FILE)

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

    fh = fopen(output_filename, b"w")
    if fh == cython.NULL:
        raise OSError("Unable to open bedGraph output file")
    fclose(fh)

    for chrom in list(chrlengths.keys()):
        (plus_tags, minus_tags) = trackI.get_locations_by_chr(chrom)
        rlength = cython.cast(cython.long, chrlengths[chrom])
        pv = _pileup_from_PN_shifted(plus_tags,
                                     minus_tags,
                                     five_shift,
                                     three_shift,
                                     rlength)
        _write_pv_to_bedGraph(pv,
                              chrom,
                              output_filename,
                              scale_factor,
                              baseline_value)
        clean_up_ndarray(pv)
    return


@cython.ccall
def pileup_and_write_pe(petrackI,
                        output_filename: bytes,
                        scale_factor: cython.float = 1.0,
                        baseline_value: cython.float = 0.0):
    """Pile up a paired-end track and write a bedGraph file."""
    chrom: bytes
    locs: cnp.ndarray
    pv: cnp.ndarray
    chrlengths: dict = petrackI.get_rlengths()
    fh: cython.pointer(FILE)

    fh = fopen(output_filename, b"w")
    if fh == cython.NULL:
        raise OSError("Unable to open bedGraph output file")
    fclose(fh)

    for chrom in sorted(chrlengths.keys()):
        locs = petrackI.get_locations_by_chr(chrom)
        if locs.dtype.names is not None and "c" in locs.dtype.names:
            pv = pileup_from_LRC(locs)
        else:
            pv = pileup_from_LR(locs)
        _write_pv_to_bedGraph(pv,
                              chrom,
                              output_filename,
                              scale_factor,
                              baseline_value)
        clean_up_ndarray(pv)
    return
