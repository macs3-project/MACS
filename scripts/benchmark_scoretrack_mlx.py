"""Benchmark utility comparing CPU ScoreTrackII vs GPU ScoreTrackMLX."""

from __future__ import annotations

import argparse
import time

import numpy as np

from MACS3.Signal.ScoreTrack import ScoreTrackII, ScoreTrackMLX

try:  # optional statistics about the MLX runtime
    import mlx.core as mx  # type: ignore
except Exception:  # pragma: no cover - mlx might be missing
    mx = None


def build_synthetic_data(num_bins: int, bin_size: int):
    pos = np.arange(bin_size, (num_bins + 1) * bin_size, bin_size, dtype=np.int64)
    rng = np.random.default_rng(42)
    chip = rng.gamma(shape=5.0, scale=2.0, size=num_bins).astype(np.float32)
    control = rng.gamma(shape=3.0, scale=2.0, size=num_bins).astype(np.float32)
    return pos, chip, control


def populate_scoretrack(track, pos, chip, control, chrom=b"chr1"):
    if isinstance(track, ScoreTrackII):
        track.add_chromosome(chrom, len(pos))
    else:
        track.add_chromosome(chrom)
    for endpos, treatment, ctrl in zip(pos, chip, control):
        track.add(chrom, int(endpos), float(treatment), float(ctrl))
    track.finalize()
    return track


def benchmark(num_bins: int, bin_size: int, cutoff: float, use_mlx: bool = True):
    chrom = b"chr1"
    pos, chip, control = build_synthetic_data(num_bins, bin_size)

    # CPU path
    cpu_track = populate_scoretrack(ScoreTrackII(30.0, 30.0), pos, chip, control, chrom)
    start = time.perf_counter()
    cpu_track.change_score_method(ord('p'))
    cpu_p_time = time.perf_counter() - start
    start = time.perf_counter()
    cpu_track.change_score_method(ord('q'))
    cpu_q_time = time.perf_counter() - start
    start = time.perf_counter()
    cpu_track.call_peaks(cutoff=cutoff)
    cpu_peak_time = time.perf_counter() - start

    # MLX path
    gpu_track = populate_scoretrack(ScoreTrackMLX(30.0, 30.0, approx_pvalue=True), pos, chip, control, chrom)
    start = time.perf_counter()
    gpu_track.compute_pvalue()
    gpu_p_time = time.perf_counter() - start
    start = time.perf_counter()
    gpu_track.compute_qvalue()
    gpu_q_time = time.perf_counter() - start
    start = time.perf_counter()
    gpu_track.call_peaks(cutoff=cutoff)
    gpu_peak_time = time.perf_counter() - start

    return {
        "cpu_p": cpu_p_time,
        "cpu_q": cpu_q_time,
        "cpu_call": cpu_peak_time,
        "gpu_p": gpu_p_time,
        "gpu_q": gpu_q_time,
        "gpu_call": gpu_peak_time,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", type=int, default=50000, help="Number of genomic bins to simulate")
    parser.add_argument("--bin-size", type=int, default=50, help="Bin size in bp")
    parser.add_argument("--cutoff", type=float, default=1.0, help="Score cutoff for peak calling")
    parser.add_argument(
        "--mlx",
        choices=["auto", "cpu"],
        default="auto",
        help="Use MLX tensors (auto) or force the MLX implementation to run on NumPy/CPU.",
    )
    args = parser.parse_args()

    print(f"MLX default device: {mx.default_device() if mx is not None else 'Unavailable'}")

    use_mlx = args.mlx != "cpu"
    if ScoreTrackMLX is None:
        raise SystemExit("ScoreTrackMLX is unavailable; ensure MACS3.Signal.ScoreTrackMLX can be imported.")

    results = benchmark(args.bins, args.bin_size, args.cutoff, use_mlx=use_mlx)
    device_label = "MLX tensor" if use_mlx and mx is not None else "NumPy fallback"
    print(f"CPU p-value: {results['cpu_p']:.3f}s | q-value: {results['cpu_q']:.3f}s | call_peaks: {results['cpu_call']:.3f}s")
    print(f"{device_label} p-value: {results['gpu_p']:.3f}s | q-value: {results['gpu_q']:.3f}s | call_peaks: {results['gpu_call']:.3f}s")
    total_cpu = results['cpu_p'] + results['cpu_q'] + results['cpu_call']
    total_gpu = results['gpu_p'] + results['gpu_q'] + results['gpu_call']
    print(f"Total CPU: {total_cpu:.3f}s | Total {device_label}: {total_gpu:.3f}s")
    if total_gpu > 0:
        print(f"Overall speedup (CPU/GPU): {total_cpu / total_gpu:.2f}x")


if __name__ == "__main__":
    main()
