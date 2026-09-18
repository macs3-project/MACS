#!/usr/bin/env python
"""Quick benchmark for ScoreTrackTorch.

Generates synthetic pileups and runs p/q-value computation and peak calling
on the requested torch backend. Designed to compare CPU vs GPU (mps/cuda)
timings on developer machines.
"""

import argparse
import time

import torch

from MACS3.Signal.ScoreTrackTorch import ScoreTrackTorch


def _sync(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elif device.type == "mps":
        torch.mps.synchronize()


def make_track(n: int, backend: str) -> ScoreTrackTorch:
    track = ScoreTrackTorch(treat_depth=1.0, ctrl_depth=1.0, backend=backend)
    track.add_chromosome(b"chr1", n)
    # generate on CPU to avoid MPS per-element transfer overhead
    positions = torch.arange(1, n + 1, device="cpu", dtype=torch.int64)
    treatment = torch.abs(torch.randn(n, device="cpu")) * 20 + 5
    control = torch.abs(torch.randn(n, device="cpu")) * 10 + 5

    # bulk transfer once to target device; avoid Python loops that fragment MPS memory
    chrom_data = track.data[b"chr1"]
    chrom_data.pos[:n] = positions.to(track.device)
    chrom_data.treatment[:n] = treatment.to(track.device)
    chrom_data.control[:n] = control.to(track.device)
    chrom_data.length = n
    return track


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", default="mps", choices=["cpu", "cuda", "mps"], help="torch backend to test")
    parser.add_argument(
        "--points",
        type=int,
        default=50_000,
        help="number of intervals to generate (reduce for slower GPU backends)",
    )
    parser.add_argument("--cutoff", type=float, default=5.0, help="score cutoff for peak calling")
    args = parser.parse_args()

    track = make_track(args.points, args.backend)
    device = track.device
    print(f"Using device: {device}")

    timings = {}

    start = time.perf_counter()
    track.compute_pvalue()
    _sync(device)
    timings["pvalue"] = time.perf_counter() - start

    start = time.perf_counter()
    track.compute_qvalue()
    _sync(device)
    timings["qvalue"] = time.perf_counter() - start

    start = time.perf_counter()
    peaks = track.call_peaks(cutoff=args.cutoff, min_length=50, max_gap=20)
    _sync(device)
    timings["call_peaks"] = time.perf_counter() - start

    print(f"Intervals: {args.points}, peaks found: {len(peaks.peaks.get(b'chr1', []))}")
    for name, dt in timings.items():
        print(f"{name:12s}: {dt:0.4f} sec")

    if device.type == "mps":
        try:
            mem = torch.mps.current_allocated_memory()
            print(f"MPS memory used: {mem/1024/1024:0.2f} MB")
        except Exception:
            pass
    elif device.type == "cuda":
        try:
            mem = torch.cuda.memory_allocated(device)
            print(f"CUDA memory used: {mem/1024/1024:0.2f} MB")
        except Exception:
            pass


if __name__ == "__main__":
    main()
