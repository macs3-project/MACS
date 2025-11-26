"""Benchmark BedGraph overlay into ScoreTrackII vs ScoreTrackMLX."""

from __future__ import annotations

import argparse
import time
from typing import Dict, List, Tuple

import numpy as np

from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.BedGraphMLX import BedGraphTrackMLX
from MACS3.Signal.ScoreTrack import ScoreTrackMLX


def build_synthetic_bedgraph(num_bins: int, bin_size: int):
    pos = np.arange(bin_size, (num_bins + 1) * bin_size, bin_size, dtype=np.int64)
    rng = np.random.default_rng(42)
    chip = rng.gamma(shape=5.0, scale=2.0, size=num_bins).astype(np.float32)
    control = rng.gamma(shape=3.0, scale=2.0, size=num_bins).astype(np.float32)
    return pos, chip, control


def fill_bedgraph(track, pos, values, chrom=b"chr1"):
    start = 0
    for end, val in zip(pos, values):
        track.add_loc(chrom, start, int(end), float(val))
        start = int(end)
    return track


def peaks_to_intervals(peaks) -> Dict[bytes, List[Tuple[int, int]]]:
    intervals: Dict[bytes, List[Tuple[int, int]]] = {}
    for chrom in peaks.get_chr_names():
        spans = sorted((int(p['start']), int(p['end'])) for p in peaks.get_data_from_chrom(chrom))
        intervals[chrom] = spans
    return intervals


def peak_summary(peaks):
    intervals = peaks_to_intervals(peaks)
    total = sum(len(v) for v in intervals.values())
    bp = sum(end - start for spans in intervals.values() for start, end in spans)
    return total, bp


def benchmark(num_bins: int, bin_size: int, cutoff: float):
    chrom = b"chr1"
    pos, chip, control = build_synthetic_bedgraph(num_bins, bin_size)

    # CPU path via bedGraphTrackI -> ScoreTrackII
    bdg_treat = fill_bedgraph(bedGraphTrackI(), pos, chip, chrom)
    bdg_ctrl = fill_bedgraph(bedGraphTrackI(), pos, control, chrom)
    start = time.perf_counter()
    cpu_track = bdg_treat.make_ScoreTrackII_for_macs(bdg_ctrl, depth1=30.0, depth2=30.0)
    build_cpu_time = time.perf_counter() - start
    start = time.perf_counter()
    cpu_track.change_score_method(ord("p"))
    cpu_p_time = time.perf_counter() - start
    start = time.perf_counter()
    cpu_track.change_score_method(ord("q"))
    cpu_q_time = time.perf_counter() - start
    start = time.perf_counter()
    cpu_peaks = cpu_track.call_peaks(cutoff=cutoff)
    cpu_peak_time = time.perf_counter() - start

    # MLX path via BedGraphTrackMLX -> ScoreTrackMLX
    bdg_treat_m = fill_bedgraph(BedGraphTrackMLX(), pos, chip, chrom)
    bdg_ctrl_m = fill_bedgraph(BedGraphTrackMLX(), pos, control, chrom)
    start = time.perf_counter()
    gpu_track = bdg_treat_m.make_ScoreTrackMLX_for_macs(bdg_ctrl_m, depth1=30.0, depth2=30.0)
    build_gpu_time = time.perf_counter() - start
    start = time.perf_counter()
    gpu_track.compute_pvalue()
    gpu_p_time = time.perf_counter() - start
    start = time.perf_counter()
    gpu_track.compute_qvalue()
    gpu_q_time = time.perf_counter() - start
    start = time.perf_counter()
    gpu_peaks = gpu_track.call_peaks(cutoff=cutoff)
    gpu_peak_time = time.perf_counter() - start

    cpu_total, cpu_bp = peak_summary(cpu_peaks)
    gpu_total, gpu_bp = peak_summary(gpu_peaks)

    return {
        "cpu_build": build_cpu_time,
        "cpu_p": cpu_p_time,
        "cpu_q": cpu_q_time,
        "cpu_call": cpu_peak_time,
        "gpu_build": build_gpu_time,
        "gpu_p": gpu_p_time,
        "gpu_q": gpu_q_time,
        "gpu_call": gpu_peak_time,
        "cpu_total": cpu_total,
        "cpu_bp": cpu_bp,
        "gpu_total": gpu_total,
        "gpu_bp": gpu_bp,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", type=int, default=200000, help="Number of genomic bins to simulate")
    parser.add_argument("--bin-size", type=int, default=50, help="Bin size in bp")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Score cutoff for peak calling")
    args = parser.parse_args()

    results = benchmark(args.bins, args.bin_size, args.cutoff)
    print(f"CPU build: {results['cpu_build']:.3f}s | p-value: {results['cpu_p']:.3f}s | q-value: {results['cpu_q']:.3f}s | call_peaks: {results['cpu_call']:.3f}s")
    print(f"MLX build: {results['gpu_build']:.3f}s | p-value: {results['gpu_p']:.3f}s | q-value: {results['gpu_q']:.3f}s | call_peaks: {results['gpu_call']:.3f}s")
    cpu_total = results["cpu_build"] + results["cpu_p"] + results["cpu_q"] + results["cpu_call"]
    gpu_total = results["gpu_build"] + results["gpu_p"] + results["gpu_q"] + results["gpu_call"]
    print(f"Total CPU: {cpu_total:.3f}s | Total MLX: {gpu_total:.3f}s")
    if gpu_total > 0:
        print(f"Overall speedup (CPU/MLX): {cpu_total / gpu_total:.2f}x")
    print(f"Peaks: CPU {results['cpu_total']} ({results['cpu_bp']} bp) | MLX {results['gpu_total']} ({results['gpu_bp']} bp)")


if __name__ == "__main__":
    main()
