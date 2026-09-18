"""Benchmark and consistency checks across Pileup, PileupV2, and PileupMLX.

Example:
    python scripts/benchmark_pileup_mlx.py --reads 200000 --chrom-size 1000000 --repeats 3
"""

from __future__ import annotations

import argparse
import time
import tracemalloc
import resource
from typing import Callable, Tuple

import numpy as np

from MACS3.Signal import Pileup as P0
from MACS3.Signal import PileupV2 as P2
from MACS3.Signal import PileupMLX as P3


Result = Tuple[object, float, float, float]  # result, seconds, tracemalloc peak (MiB), max RSS (MiB)


def _compare_structured(a: np.ndarray, b: np.ndarray, rtol=1e-5, atol=1e-8) -> bool:
    return np.array_equal(a["p"], b["p"]) and np.allclose(a["v"], b["v"], rtol=rtol, atol=atol)


def _first_diff_structured(a: np.ndarray, b: np.ndarray, rtol=1e-5, atol=1e-8) -> str:
    n = min(a.shape[0], b.shape[0])
    for i in range(n):
        if a["p"][i] != b["p"][i] or not np.isclose(a["v"][i], b["v"][i], rtol=rtol, atol=atol):
            return f"idx={i}, a=({a['p'][i]}, {a['v'][i]}), b=({b['p'][i]}, {b['v'][i]})"
    if a.shape[0] != b.shape[0]:
        if a.shape[0] > b.shape[0]:
            return f"len diff: len(a)={a.shape[0]}, len(b)={b.shape[0]}, extra a idx={n}, a=({a['p'][n]}, {a['v'][n]})"
        return f"len diff: len(a)={a.shape[0]}, len(b)={b.shape[0]}, extra b idx={n}, b=({b['p'][n]}, {b['v'][n]})"
    return "no diff"


def _list_to_struct(pv_list) -> np.ndarray:
    """Convert [p, v] list into structured array dtype=[('p','i4'),('v','f4')]."""
    p, v = pv_list
    ret = np.empty(p.shape[0], dtype=[("p", "i4"), ("v", "f4")])
    ret["p"] = p
    ret["v"] = v
    return ret


def _measure(func: Callable, *args, repeats: int = 3, **kwargs) -> Result:
    """Run func multiple times and report median time and peak memory usage."""
    times = []
    peaks = []
    rss_peaks = []
    result = None
    for _ in range(repeats):
        tracemalloc.start()
        rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        t0 = time.perf_counter()
        result = func(*args, **kwargs)
        dt = time.perf_counter() - t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        times.append(dt)
        peaks.append(peak)
        rss_peaks.append(max(rss_before, rss_after))
    peak_mem = max(peaks) / (1024 * 1024)
    peak_rss = max(rss_peaks) / 1024  # ru_maxrss is in KB on macOS/Linux
    return result, float(np.median(times)), float(peak_mem), float(peak_rss)


def build_datasets(reads: int, chrom_size: int, seed: int = 1):
    rng = np.random.default_rng(seed)

    # PE-style bounds for quick_pileup and merging
    start_poss = rng.integers(0, chrom_size - 300, size=reads, dtype=np.int32)
    frag_len = rng.integers(50, 300, size=reads, dtype=np.int32)
    end_poss = start_poss + frag_len
    start_poss.sort()
    end_poss.sort()

    start_poss_b = rng.integers(0, chrom_size - 300, size=reads, dtype=np.int32)
    frag_len_b = rng.integers(50, 300, size=reads, dtype=np.int32)
    end_poss_b = start_poss_b + frag_len_b
    start_poss_b.sort()
    end_poss_b.sort()

    # SE tags
    plus_tags = rng.integers(0, chrom_size, size=reads, dtype=np.int32)
    minus_tags = rng.integers(0, chrom_size, size=reads, dtype=np.int32)
    plus_tags.sort()
    minus_tags.sort()
    naive_pos = np.sort(np.concatenate([plus_tags[: reads // 2], minus_tags[: reads // 2]]))

    # LR arrays
    l_bounds = rng.integers(0, chrom_size - 300, size=reads, dtype=np.int32)
    l_length = rng.integers(50, 300, size=reads, dtype=np.int32)
    r_bounds = l_bounds + l_length
    lr_array = np.empty(reads, dtype=[("l", "i4"), ("r", "i4")])
    lr_array["l"] = l_bounds
    lr_array["r"] = r_bounds

    pn_array = naive_pos  # reuse for pileup_from_PN

    mapping_dict = {int(k): float(k) * 0.01 for k in np.unique(l_length)}

    return {
        "start_poss": start_poss,
        "end_poss": end_poss,
        "start_poss_b": start_poss_b,
        "end_poss_b": end_poss_b,
        "plus_tags": plus_tags,
        "minus_tags": minus_tags,
        "naive_pos": naive_pos,
        "rlength": chrom_size,
        "lr_array": lr_array,
        "pn_array": pn_array,
        "mapping_dict": mapping_dict,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", type=int, default=100_000, help="Number of fragments/reads to simulate.")
    ap.add_argument("--chrom-size", type=int, default=1_000_000, help="Chromosome size for random coordinates.")
    ap.add_argument("--extension", type=int, default=150, help="Extension length for naive pileup.")
    ap.add_argument("--repeats", type=int, default=3, help="Number of repetitions per benchmark.")
    ap.add_argument("--seed", type=int, default=1, help="Random seed.")
    args = ap.parse_args()

    data = build_datasets(args.reads, args.chrom_size, seed=args.seed)

    cases = [
        {
            "name": "quick_pileup",
            "runner_p0": lambda: _list_to_struct(P0.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0)),
            "runner_p2": lambda: P2.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
            "runner_p3": lambda: P3.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
        {
            "name": "se_all_in_one_pileup",
            "runner_p0": lambda: _list_to_struct(
                P0.se_all_in_one_pileup(
                    data["plus_tags"], data["minus_tags"], -50, 100, data["rlength"], 1.0, 0.0
                )
            ),
            "runner_p2": lambda: P2.se_all_in_one_pileup(
                data["plus_tags"], data["minus_tags"], -50, 100, data["rlength"], 1.0, 0.0
            ),
            "runner_p3": lambda: P3.se_all_in_one_pileup(
                data["plus_tags"], data["minus_tags"], -50, 100, data["rlength"], 1.0, 0.0
            ),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
        {
            "name": "naive_quick_pileup",
            "runner_p0": lambda: _list_to_struct(P0.naive_quick_pileup(data["naive_pos"], args.extension)),
            "runner_p2": lambda: P2.naive_quick_pileup(data["naive_pos"], args.extension),
            "runner_p3": lambda: P3.naive_quick_pileup(data["naive_pos"], args.extension),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
        {
            "name": "over_two_pv_array",
            "runner_p0": lambda: _list_to_struct(
                P0.over_two_pv_array(
                    P0.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
                    P0.quick_pileup(data["start_poss_b"], data["end_poss_b"], 1.0, 0.0),
                    func="max",
                )
            ),
            "runner_p2": lambda: P2.over_two_pv_array(
                P2.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
                P2.quick_pileup(data["start_poss_b"], data["end_poss_b"], 1.0, 0.0),
                func="max",
            ),
            "runner_p3": lambda: P3.over_two_pv_array(
                P3.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
                P3.quick_pileup(data["start_poss_b"], data["end_poss_b"], 1.0, 0.0),
                func="max",
            ),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
        {
            "name": "naive_call_peaks",
            "runner_p0": lambda: np.array(
                P0.naive_call_peaks(
                    P0.naive_quick_pileup(data["naive_pos"], args.extension), min_v=1.0, max_gap=75, min_length=100
                ),
                dtype=[("p", "i4"), ("v", "f4")],
            ),
            "runner_p2": lambda: np.array(
                P2.naive_call_peaks(
                    P2.naive_quick_pileup(data["naive_pos"], args.extension), min_v=1.0, max_gap=75, min_length=100
                ),
                dtype=[("p", "i4"), ("v", "f4")],
            ),
            "runner_p3": lambda: np.array(
                P3.naive_call_peaks(
                    P3.naive_quick_pileup(data["naive_pos"], args.extension), min_v=1.0, max_gap=75, min_length=100
                ),
                dtype=[("p", "i4"), ("v", "f4")],
            ),
            "comparator": lambda a, b: np.array_equal(a, b),
            "diff": lambda a, b: f"len a={len(a)}, len b={len(b)}",
        },
        {
            "name": "pileup_from_LR",
            "runner_p0": None,  # not available in Pileup.py
            "runner_p2": lambda: P2.pileup_from_LR(data["lr_array"]),
            "runner_p3": lambda: P3.pileup_from_LR(data["lr_array"]),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
        {
            "name": "pileup_from_LR_hmmratac",
            "runner_p0": None,
            "runner_p2": lambda: P2.pileup_from_LR_hmmratac(data["lr_array"], data["mapping_dict"]),
            "runner_p3": lambda: P3.pileup_from_LR_hmmratac(data["lr_array"], data["mapping_dict"]),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
        {
            "name": "pileup_from_PN",
            "runner_p0": None,
            "runner_p2": lambda: P2.pileup_from_PN(data["pn_array"], data["pn_array"], args.extension),
            "runner_p3": lambda: P3.pileup_from_PN(data["pn_array"], data["pn_array"], args.extension),
            "comparator": _compare_structured,
            "diff": _first_diff_structured,
        },
    ]

    print(f"Simulated reads: {args.reads:,}, chrom size: {args.chrom_size:,}, repeats: {args.repeats}")
    for case in cases:
        res_p2, t_p2, mem_p2, rss_p2 = _measure(case["runner_p2"], repeats=args.repeats)
        res_p3, t_p3, mem_p3, rss_p3 = _measure(case["runner_p3"], repeats=args.repeats)
        match_p0 = diff_p0 = None
        t_p0 = mem_p0 = rss_p0 = None
        if case["runner_p0"] is not None:
            res_p0, t_p0, mem_p0, rss_p0 = _measure(case["runner_p0"], repeats=args.repeats)
            match_p0 = case["comparator"](res_p0, res_p2)
            if not match_p0 and case["diff"] is not None:
                diff_p0 = case["diff"](res_p0, res_p2)

        match_p3 = case["comparator"](res_p3, res_p2)
        diff_p3 = ""
        if not match_p3 and case["diff"] is not None:
            diff_p3 = case["diff"](res_p3, res_p2)

        parts = [
            f"{case['name']:24s}",
            f"P0 match: {match_p0}" if case["runner_p0"] is not None else "P0 match: N/A",
            f"P2 time: {t_p2*1000:8.2f} ms, py {mem_p2:6.2f} MiB, rss {rss_p2:6.2f} MiB",
            f"P3 time: {t_p3*1000:8.2f} ms, py {mem_p3:6.2f} MiB, rss {rss_p3:6.2f} MiB",
            f"P0 time: {t_p0*1000:8.2f} ms, py {mem_p0:6.2f} MiB, rss {rss_p0:6.2f} MiB" if case["runner_p0"] is not None else "P0 time: N/A",
        ]
        if diff_p0:
            parts.append(f"P0 diff: {diff_p0}")
        if diff_p3:
            parts.append(f"P3 diff: {diff_p3}")
        print(" | ".join(parts))


if __name__ == "__main__":
    main()
