# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-02-05 12:40:02 Tao Liu>

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

# ------------------------------------
# utility internal functions
# ------------------------------------


@cython.cfunc
@cython.inline
def mean(a: float, b: float) -> float:
    return (a + b) / 2


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
    """ Pileup a FWTrackI object and write the pileup result into a
    bedGraph file.

    This function calls pure C functions from cPosValCalculation.

    This function is currently only used `macs3 pileup` cmd.
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
    """ Pileup a PETrackI object and write the pileup result into a
    bedGraph file.

    This function calls pure C functions from cPosValCalculation.

    This function is currently only used `macs3 pileup` cmd.
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

    """
    p: cython.int
    pre_p: cython.int
    pileup: cython.int = 0

    i_s: cython.long = 0        # index of start_poss
    i_e: cython.long = 0        # index of end_poss
    i: cython.long
    I: cython.long = 0
    lx: cython.long

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
    end_poss = np.concatenate((plus_tags+three_shift, minus_tags+five_shift))

    # sort
    start_poss.sort()
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

# quick pileup implemented in cython


@cython.ccall
def quick_pileup(start_poss: cnp.ndarray,
                 end_poss: cnp.ndarray,
                 scale_factor: cython.float,
                 baseline_value: cython.float) -> list:
    """Return pileup given plus strand and minus strand positions of fragments.

    A super-fast and simple algorithm proposed by Jie Wang. It will
    take sorted start positions and end positions, then compute pileup
    values.

    It will return a pileup result in similar structure as
    bedGraph. There are two python arrays:

    [end positions, values] or [p,v]

    Two arrays have the same length and can be matched by index. End
    position at index x (p[x]) record continuous value of v[x] from
    p[x-1] to p[x].

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
    """Simple pileup, every tag will be extended left and right with
    length `extension`.

    Note: Assumption is that `poss` has to be pre-sorted! There is no
    check on whether it's sorted.

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
    """Merge two position-value arrays. For intersection regions, take
    the maximum value within region.

    pv_array1 and pv_array2 are [p,v] type lists, same as the output
    from quick_pileup function. 'p' and 'v' are numpy arrays of int32
    and float32.

    available operations are 'max', 'min', and 'mean'
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

    if func == "max":
        f = max
    elif func == "min":
        f = min
    elif func == "mean":
        f = mean
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
        ret_v_ptr[0] = f(a1_v_ptr[0], a2_v_ptr[0])
        I += 1
        if a1_pos_ptr[0] < a2_pos_ptr[0]:
            # clip a region from pre_p to p1, then set pre_p as p1.
            ret_pos_ptr[0] = a1_pos_ptr[0]
            ret_pos_ptr += 1
            ret_v_ptr += 1
            # pre_p = a1_pos_ptr[0]
            # call for the next p1 and v1
            a1_pos_ptr += 1
            a1_v_ptr += 1
            i1 += 1
        elif a1_pos_ptr[0] > a2_pos_ptr[0]:
            # clip a region from pre_p to p2, then set pre_p as p2.
            ret_pos_ptr[0] = a2_pos_ptr[0]
            ret_pos_ptr += 1
            ret_v_ptr += 1
            # pre_p = a2_pos_ptr[0]
            # call for the next p1 and v1
            a2_pos_ptr += 1
            a2_v_ptr += 1
            i2 += 1
        else:
            # from pre_p to p1 or p2, then set pre_p as p1 or p2.
            ret_pos_ptr[0] = a1_pos_ptr[0]
            ret_pos_ptr += 1
            ret_v_ptr += 1
            # pre_p = a1_pos_ptr[0]
            # call for the next p1, v1, p2, v2.
            a1_pos_ptr += 1
            a1_v_ptr += 1
            i1 += 1
            a2_pos_ptr += 1
            a2_v_ptr += 1
            i2 += 1

    ret_pos.resize(I, refcheck=False)
    ret_v.resize(I, refcheck=False)
    return [ret_pos, ret_v]


@cython.ccall
def naive_call_peaks(pv_array: list, min_v: cython.float,
                     max_v: cython.float = 1e30,
                     max_gap: cython.int = 50,
                     min_length: cython.int = 200):

    pre_p: cython.int
    p: cython.int
    i: cython.int
    x: cython.long           # index used for searching the first peak
    v: cython.double
    peak_content: list     # (pre_p, p, v) for each region in the peak
    ret_peaks: list = []   # returned peak summit and height

    peak_content = []
    (ps, vs) = pv_array
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
