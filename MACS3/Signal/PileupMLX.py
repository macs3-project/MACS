"""Pileup functions mirroring :mod:`MACS3.Signal.PileupV2` using NumPy."""

from __future__ import annotations

from typing import Callable

import numpy as np


def mapping_function_always_1(L: int, R: int) -> float:
    return 1.0


def _pileup_from_pv_numpy(positions_np: np.ndarray, deltas_np: np.ndarray) -> np.ndarray:
    """Pure NumPy implementation for speed and correctness."""
    if positions_np.size == 0:
        return np.zeros(shape=0, dtype=[("p", "i4"), ("v", "f4")])

    order = np.argsort(positions_np, kind="quicksort")
    p_sorted = positions_np[order]
    d_sorted = deltas_np[order].astype(np.float64, copy=False)

    uniq_pos, idx, counts = np.unique(p_sorted, return_index=True, return_counts=True)
    summed = np.add.reduceat(d_sorted, idx)  # float64 accumulation
    if summed.size:
        exclusive = np.empty_like(summed, dtype=np.float64)
        exclusive[0] = 0.0
        np.cumsum(summed[:-1], dtype=np.float64, out=exclusive[1:])
    else:
        exclusive = summed

    pos_np = uniq_pos.astype(np.int32, copy=False)
    val_np = exclusive

    if pos_np.size > 1:
        # keep the LAST coordinate for each run of identical values
        changes = np.empty_like(val_np, dtype=bool)
        changes[0] = False
        changes[1:] = val_np[1:] != val_np[:-1]
        run_ends = np.nonzero(changes)[0]
        run_ends = np.concatenate((run_ends, np.array([val_np.size], dtype=run_ends.dtype)))
        last_idx = run_ends - 1
        pos_np = pos_np[last_idx]
        val_np = val_np[last_idx]

    if pos_np.size and pos_np[0] == 0 and val_np[0] == 0:
        pos_np = pos_np[1:]
        val_np = val_np[1:]

    ret = np.empty(pos_np.shape[0], dtype=[("p", "i4"), ("v", "f4")])
    ret["p"] = pos_np
    ret["v"] = val_np.astype(np.float32, copy=False)
    return ret


def _pileup_from_pv(positions_mx, deltas_mx) -> np.ndarray:
    """Core pileup routine; run purely in NumPy."""
    pos_np = np.asarray(positions_mx, dtype=np.int32, order="C")
    deltas_np = np.asarray(deltas_mx, dtype=np.float32, order="C")
    return _pileup_from_pv_numpy(pos_np, deltas_np)


def _pileup_from_bounds(start_poss, end_poss) -> np.ndarray:
    """Pileup using start/end bounds (arrays of equal length)."""
    start_np = np.asarray(start_poss, dtype=np.int32)
    end_np = np.asarray(end_poss, dtype=np.int32)
    positions = np.concatenate([start_np, end_np]).astype(np.int32, copy=False)
    deltas = np.concatenate(
        [
            np.ones_like(start_np, dtype=np.float32),
            -np.ones_like(end_np, dtype=np.float32),
        ]
    )
    return _pileup_from_pv_numpy(positions, deltas)


def _weights_from_mapping(L: np.ndarray,
                          R: np.ndarray,
                          mapping_func: Callable[[int, int], float]) -> np.ndarray:
    if mapping_func is mapping_function_always_1:
        return np.ones(L.shape[0], dtype=np.float32)
    return np.fromiter(
        (mapping_func(int(l), int(r)) for l, r in zip(L, R)),
        dtype=np.float32,
        count=L.shape[0],
    )


def pileup_from_LR_hmmratac(LR_array: np.ndarray,
                            mapping_dict: dict) -> np.ndarray:
    L = np.asarray(LR_array["l"], dtype=np.int32)
    R = np.asarray(LR_array["r"], dtype=np.int32)
    length = R - L
    weights = np.asarray([mapping_dict[int(x)] for x in length], dtype=np.float32)

    positions = np.empty(L.shape[0] * 2, dtype=np.int32)
    deltas = np.empty(L.shape[0] * 2, dtype=np.float32)
    positions[0::2] = L
    positions[1::2] = R
    deltas[0::2] = weights
    deltas[1::2] = -weights
    return _pileup_from_pv(positions, deltas)


def pileup_from_LR(LR_array: np.ndarray,
                   mapping_func: Callable[[int, int], float] = mapping_function_always_1) -> np.ndarray:
    L = np.asarray(LR_array["l"], dtype=np.int32)
    R = np.asarray(LR_array["r"], dtype=np.int32)
    weights = _weights_from_mapping(L, R, mapping_func)

    positions = np.empty(L.shape[0] * 2, dtype=np.int32)
    deltas = np.empty(L.shape[0] * 2, dtype=np.float32)
    positions[0::2] = L
    positions[1::2] = R
    deltas[0::2] = weights
    deltas[1::2] = -weights
    return _pileup_from_pv(positions, deltas)


def pileup_from_LRC(LRC_array: np.ndarray,
                    mapping_func: Callable[[int, int], float] = mapping_function_always_1) -> np.ndarray:
    L = np.asarray(LRC_array["l"], dtype=np.int32)
    R = np.asarray(LRC_array["r"], dtype=np.int32)
    C = np.asarray(LRC_array["c"], dtype=np.float32)
    weights = _weights_from_mapping(L, R, mapping_func) * C

    positions = np.empty(L.shape[0] * 2, dtype=np.int32)
    deltas = np.empty(L.shape[0] * 2, dtype=np.float32)
    positions[0::2] = L
    positions[1::2] = R
    deltas[0::2] = weights
    deltas[1::2] = -weights
    return _pileup_from_pv(positions, deltas)


def pileup_from_PN(P_array: np.ndarray, N_array: np.ndarray,
                   extsize: int) -> np.ndarray:
    P = np.asarray(P_array, dtype=np.int32)
    N = np.asarray(N_array, dtype=np.int32)
    if P.shape[0] != N.shape[0]:
        raise ValueError("P_array and N_array must share the same length")

    plus_start = P
    plus_end = P + extsize
    minus_end = N
    minus_start = N - extsize

    positions = np.concatenate([plus_start, plus_end, minus_start, minus_end]).astype(np.int32, copy=False)
    deltas = np.concatenate([
        np.ones_like(plus_start, dtype=np.float32),
        -np.ones_like(plus_end, dtype=np.float32),
        np.ones_like(minus_start, dtype=np.float32),
        -np.ones_like(minus_end, dtype=np.float32),
    ])
    return _pileup_from_pv(positions, deltas)


def quick_pileup(start_poss: np.ndarray,
                 end_poss: np.ndarray,
                 scale_factor: float,
                 baseline_value: float) -> np.ndarray:
    pileup = _pileup_from_bounds(start_poss, end_poss)
    if pileup.shape[0] == 0:
        return pileup
    if scale_factor != 1:
        pileup["v"] *= scale_factor
    if baseline_value != 0:
        pileup["v"] = np.maximum(pileup["v"], baseline_value)
    return pileup


def se_all_in_one_pileup(plus_tags: np.ndarray,
                         minus_tags: np.ndarray,
                         five_shift: int,
                         three_shift: int,
                         rlength: int,
                         scale_factor: float,
                         baseline_value: float) -> np.ndarray:
    plus = np.asarray(plus_tags, dtype=np.int32)
    minus = np.asarray(minus_tags, dtype=np.int32)

    start_plus = plus - five_shift
    start_minus = minus - three_shift
    end_plus = plus + three_shift
    end_minus = minus + five_shift

    np.clip(start_plus, 0, rlength, out=start_plus)
    np.clip(start_minus, 0, rlength, out=start_minus)
    np.clip(end_plus, 0, rlength, out=end_plus)
    np.clip(end_minus, 0, rlength, out=end_minus)

    start_poss = np.concatenate([start_plus, start_minus]).astype(np.int32, copy=False)
    end_poss = np.concatenate([end_plus, end_minus]).astype(np.int32, copy=False)

    positions = np.concatenate([start_poss, end_poss]).astype(np.int32, copy=False)
    deltas = np.concatenate(
        [
            np.ones_like(start_poss, dtype=np.float32),
            -np.ones_like(end_poss, dtype=np.float32),
        ]
    )

    pileup = _pileup_from_pv_numpy(positions, deltas)

    if pileup.shape[0] == 0:
        return pileup
    if scale_factor != 1:
        pileup["v"] *= scale_factor
    if baseline_value != 0:
        pileup["v"] = np.maximum(pileup["v"], baseline_value)
    return pileup


def naive_quick_pileup(sorted_poss: np.ndarray, extension: int) -> np.ndarray:
    poss = np.asarray(sorted_poss, dtype=np.int32)
    if poss.shape[0] == 0:
        raise ValueError("length is 0")
    start_poss = poss - extension
    start_poss[start_poss < 0] = 0
    end_poss = poss + extension
    return quick_pileup(start_poss, end_poss, 1.0, 0.0)


def over_two_pv_array(pv_array1,
                      pv_array2,
                      func: str = "max") -> np.ndarray:
    """Merge two PV arrays with element-wise reduction over intervals."""
    if func == "max":
        reducer = np.maximum
    elif func == "min":
        reducer = np.minimum
    elif func == "mean":
        reducer = lambda a, b: (a + b) / 2.0
    else:
        raise ValueError("Invalid function")

    if pv_array1.shape[0] == 0:
        return pv_array2.copy()
    if pv_array2.shape[0] == 0:
        return pv_array1.copy()

    p1 = np.asarray(pv_array1["p"], dtype=np.int32)
    v1 = np.asarray(pv_array1["v"], dtype=np.float32)
    p2 = np.asarray(pv_array2["p"], dtype=np.int32)
    v2 = np.asarray(pv_array2["v"], dtype=np.float32)

    i1 = i2 = 0
    l1 = p1.shape[0]
    l2 = p2.shape[0]

    ret_p: list[int] = []
    ret_v: list[float] = []

    while i1 < l1 and i2 < l2:
        p1_cur = p1[i1]
        p2_cur = p2[i2]
        ret_v.append(float(reducer(v1[i1], v2[i2])))
        if p1_cur < p2_cur:
            ret_p.append(int(p1_cur))
            i1 += 1
        elif p1_cur > p2_cur:
            ret_p.append(int(p2_cur))
            i2 += 1
        else:
            ret_p.append(int(p1_cur))
            i1 += 1
            i2 += 1

    ret = np.empty(len(ret_p), dtype=[("p", "i4"), ("v", "f4")])
    ret["p"] = np.array(ret_p, dtype=np.int32)
    ret["v"] = np.array(ret_v, dtype=np.float32)
    return ret


def naive_call_peaks(pv_array: np.ndarray,
                     min_v: float,
                     max_v: float = 1e30,
                     max_gap: int = 50,
                     min_length: int = 200):
    ret_peaks = []
    peak_content = []
    psn = iter(pv_array["p"]).__next__
    vsn = iter(pv_array["v"]).__next__
    pre_p = 0
    x = 0
    while True:
        try:
            p = psn()
            v = vsn()
        except Exception:
            break
        x += 1
        if v > min_v:
            peak_content = [(pre_p, p, v)]
            pre_p = p
            break
        else:
            pre_p = p

    for _ in range(x, len(pv_array)):
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
            peak_content = [(pre_p, p, v)]
        pre_p = p

    if peak_content and peak_content[-1][1] - peak_content[0][0] >= min_length:
        __close_peak(peak_content, ret_peaks, max_v, min_length)
    return ret_peaks


def __close_peak(peak_content,
                 peaks,
                 max_v: float,
                 min_length: int):
    tsummit = []
    summit_value = 0

    for (tstart, tend, tvalue) in peak_content:
        if summit_value < tvalue:
            tsummit = [int((tend + tstart) / 2)]
            summit_value = tvalue
        elif summit_value == tvalue:
            tsummit.append(int((tend + tstart) / 2))
    summit = tsummit[int((len(tsummit) + 1) / 2) - 1]
    if summit_value < max_v and peak_content[-1][1] - peak_content[0][0] >= min_length:
        peaks.append((summit, summit_value))
    return


def _write_pv_array_to_bedgraph(pv_array: np.ndarray,
                                chrom,
                                output_filename: bytes):
    chrom_str = chrom.decode() if isinstance(chrom, (bytes, bytearray)) else str(chrom)
    pre_p = 0
    with open(output_filename, "a") as fh:
        for (p, v) in pv_array:
            if p > pre_p:
                fh.write(f"{chrom_str}\t{pre_p}\t{int(p)}\t{float(v)}\n")
            pre_p = p
    return


def pileup_and_write_se(trackI,
                        output_filename: bytes,
                        d: int,
                        scale_factor: float,
                        baseline_value: float = 0.0,
                        directional: bool = True,
                        halfextension: bool = True):
    if directional:
        if halfextension:
            five_shift = d // -4
            three_shift = d * 3 // 4
        else:
            five_shift = 0
            three_shift = d
    else:
        if halfextension:
            five_shift = d // 4
            three_shift = five_shift
        else:
            five_shift = d // 2
            three_shift = d - five_shift

    chrlengths = trackI.get_rlengths()
    chroms = list(chrlengths.keys())

    with open(output_filename, "w"):
        pass

    for chrom in chroms:
        plus_tags, minus_tags = trackI.get_locations_by_chr(chrom)
        rlength = int(chrlengths[chrom])
        pileup = se_all_in_one_pileup(plus_tags,
                                      minus_tags,
                                      five_shift,
                                      three_shift,
                                      rlength,
                                      scale_factor,
                                      baseline_value)
        _write_pv_array_to_bedgraph(pileup, chrom, output_filename)
    return


def pileup_and_write_pe(petrackI,
                        output_filename: bytes,
                        scale_factor: float = 1.0,
                        baseline_value: float = 0.0):
    chrlengths = petrackI.get_rlengths()
    chroms = list(chrlengths.keys())

    with open(output_filename, "w"):
        pass

    for chrom in chroms:
        locs = petrackI.get_locations_by_chr(chrom)
        pileup = pileup_from_LR(locs)
        if pileup.shape[0] > 0:
            if scale_factor != 1:
                pileup["v"] *= scale_factor
            if baseline_value != 0:
                pileup["v"] = np.maximum(pileup["v"], baseline_value)
        _write_pv_array_to_bedgraph(pileup, chrom, output_filename)
    return
