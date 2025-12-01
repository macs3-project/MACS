"""Benchmark and consistency checks for Pileup.py vs PileupV2.py.

This script compares speed and memory usage while ensuring equivalent
results across the two implementations. It focuses on:
    - quick_pileup
    - se_all_in_one_pileup
    - naive_quick_pileup
    - over_two_pv_array
    - naive_call_peaks

Example:
    python scripts/benchmark_pileup_performance.py --reads 200000 --chrom-size 1000000 --repeats 3
"""

from __future__ import annotations

import argparse
import time
import tracemalloc
import resource
from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np

from MACS3.Signal import Pileup as P1
from MACS3.Signal import PileupV2 as P2


Result = Tuple[np.ndarray, float, float, float]  # result, seconds, tracemalloc peak (MiB), max RSS (MiB)


def _list_to_struct(pv_list) -> np.ndarray:
    """Convert [p, v] list into structured array dtype=[('p','i4'),('v','f4')]."""
    p, v = pv_list
    ret = np.empty(p.shape[0], dtype=[('p', 'i4'), ('v', 'f4')])
    ret['p'] = p
    ret['v'] = v
    return ret


def _compare_structured(a: np.ndarray, b: np.ndarray, rtol=1e-5, atol=1e-8) -> bool:
    return np.array_equal(a['p'], b['p']) and np.allclose(a['v'], b['v'], rtol=rtol, atol=atol)


def _first_diff_structured(a: np.ndarray, b: np.ndarray, rtol=1e-5, atol=1e-8) -> str:
    n = min(a.shape[0], b.shape[0])
    for i in range(n):
        if a['p'][i] != b['p'][i] or not np.isclose(a['v'][i], b['v'][i], rtol=rtol, atol=atol):
            return f"idx={i}, a=({a['p'][i]}, {a['v'][i]}), b=({b['p'][i]}, {b['v'][i]})"
    if a.shape[0] != b.shape[0]:
        if a.shape[0] > b.shape[0]:
            return f"len diff: len(a)={a.shape[0]}, len(b)={b.shape[0]}, extra a idx={n}, a=({a['p'][n]}, {a['v'][n]})"
        else:
            return f"len diff: len(a)={a.shape[0]}, len(b)={b.shape[0]}, extra b idx={n}, b=({b['p'][n]}, {b['v'][n]})"
    return "no diff"


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
    # convert bytes to MiB
    peak_mem = max(peaks) / (1024 * 1024)
    peak_rss = max(rss_peaks) / 1024  # ru_maxrss is in KB on Linux/macOS
    return result, float(np.median(times)), float(peak_mem), float(peak_rss)


def build_datasets(reads: int, chrom_size: int, seed: int = 1):
    rng = np.random.default_rng(seed)

    # PE data for quick_pileup and merging
    start_poss = rng.integers(0, chrom_size - 300, size=reads, dtype=np.int32)
    frag_len = rng.integers(50, 300, size=reads, dtype=np.int32)
    end_poss = start_poss + frag_len
    start_poss.sort()
    end_poss.sort()

    # second set for over_two_pv_array
    start_poss_b = rng.integers(0, chrom_size - 300, size=reads, dtype=np.int32)
    frag_len_b = rng.integers(50, 300, size=reads, dtype=np.int32)
    end_poss_b = start_poss_b + frag_len_b
    start_poss_b.sort()
    end_poss_b.sort()

    # SE data for se_all_in_one_pileup and naive_quick_pileup
    plus_tags = rng.integers(0, chrom_size, size=reads, dtype=np.int32)
    minus_tags = rng.integers(0, chrom_size, size=reads, dtype=np.int32)
    plus_tags.sort()
    minus_tags.sort()

    naive_pos = np.sort(np.concatenate([plus_tags[: reads // 2], minus_tags[: reads // 2]]))

    return {
        "start_poss": start_poss,
        "end_poss": end_poss,
        "start_poss_b": start_poss_b,
        "end_poss_b": end_poss_b,
        "plus_tags": plus_tags,
        "minus_tags": minus_tags,
        "naive_pos": naive_pos,
        "rlength": chrom_size,
    }


@dataclass
class Case:
    name: str
    runner_v1: Callable[[], np.ndarray]
    runner_v2: Callable[[], np.ndarray]
    comparator: Callable[[np.ndarray, np.ndarray], bool]


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
        Case(
            name="quick_pileup",
            runner_v1=lambda: _list_to_struct(
                P1.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0)
            ),
            runner_v2=lambda: P2.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
            comparator=_compare_structured,
        ),
        Case(
            name="se_all_in_one_pileup",
            runner_v1=lambda: _list_to_struct(
                P1.se_all_in_one_pileup(
                    data["plus_tags"], data["minus_tags"], -50, 100, data["rlength"], 1.0, 0.0
                )
            ),
            runner_v2=lambda: P2.se_all_in_one_pileup(
                data["plus_tags"], data["minus_tags"], -50, 100, data["rlength"], 1.0, 0.0
            ),
            comparator=_compare_structured,
        ),
        Case(
            name="naive_quick_pileup",
            runner_v1=lambda: _list_to_struct(P1.naive_quick_pileup(data["naive_pos"], args.extension)),
            runner_v2=lambda: P2.naive_quick_pileup(data["naive_pos"], args.extension),
            comparator=_compare_structured,
        ),
        Case(
            name="over_two_pv_array",
            runner_v1=lambda: _list_to_struct(
                P1.over_two_pv_array(
                    P1.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
                    P1.quick_pileup(data["start_poss_b"], data["end_poss_b"], 1.0, 0.0),
                    func="max",
                )
            ),
            runner_v2=lambda: P2.over_two_pv_array(
                P2.quick_pileup(data["start_poss"], data["end_poss"], 1.0, 0.0),
                P2.quick_pileup(data["start_poss_b"], data["end_poss_b"], 1.0, 0.0),
                func="max",
            ),
            comparator=_compare_structured,
        ),
        Case(
            name="naive_call_peaks",
            runner_v1=lambda: np.array(
                P1.naive_call_peaks(
                    P1.naive_quick_pileup(data["naive_pos"], args.extension), min_v=1.0, max_gap=75, min_length=100
                ),
                dtype=[('p', 'i4'), ('v', 'f4')],
            ),
            runner_v2=lambda: np.array(
                P2.naive_call_peaks(
                    P2.naive_quick_pileup(data["naive_pos"], args.extension), min_v=1.0, max_gap=75, min_length=100
                ),
                dtype=[('p', 'i4'), ('v', 'f4')],
            ),
            comparator=lambda a, b: np.array_equal(a, b),
        ),
    ]

    print(f"Simulated reads: {args.reads:,}, chrom size: {args.chrom_size:,}, repeats: {args.repeats}")
    for case in cases:
        res1, t1, mem1, rss1 = _measure(case.runner_v1, repeats=args.repeats)
        res2, t2, mem2, rss2 = _measure(case.runner_v2, repeats=args.repeats)
        consistent = case.comparator(res1, res2)
        diff_info = ""
        if not consistent:
            if isinstance(res1, np.ndarray) and res1.dtype.names == ('p', 'v'):
                diff_info = _first_diff_structured(res1, res2)
            else:
                diff_info = f"len(res1)={len(res1)}, len(res2)={len(res2)}"
        print(
            f"{case.name:24s}"
            f" | match: {consistent}"
            f" | Pileup.py: {t1*1000:8.2f} ms, peak_py {mem1:6.2f} MiB, rss {rss1:6.2f} MiB"
            f" | PileupV2: {t2*1000:8.2f} ms, peak_py {mem2:6.2f} MiB, rss {rss2:6.2f} MiB"
            f"{' | diff: ' + diff_info if diff_info else ''}"
        )


if __name__ == "__main__":
    main()
