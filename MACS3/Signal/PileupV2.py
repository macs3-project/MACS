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
    PV_array: cnp.ndarray
    pileup: cnp.ndarray

    PV_array = make_PV_from_PN(P_array, N_array, extsize)
    pileup = pileup_PV(PV_array)
    clean_up_ndarray(PV_array)
    return pileup
