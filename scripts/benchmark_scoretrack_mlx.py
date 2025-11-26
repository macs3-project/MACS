"""Benchmark utility comparing CPU ScoreTrackII vs GPU ScoreTrackMLX."""

from __future__ import annotations

import argparse
import time
from typing import Dict, List, Optional, Tuple

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


def _peak_field(peak, field: str) -> int:
    """Return integer field from PeakContent regardless of backing type."""
    try:
        return int(getattr(peak, field))
    except AttributeError:
        return int(peak[field])


def peaks_to_intervals(peaks) -> Dict[bytes, List[Tuple[int, int]]]:
    """Return sorted intervals per chromosome for a PeakIO-like object."""
    intervals: Dict[bytes, List[Tuple[int, int]]] = {}
    for chrom in peaks.get_chr_names():
        chrom_peaks = peaks.get_data_from_chrom(chrom)
        spans = sorted((_peak_field(p, "start"), _peak_field(p, "end")) for p in chrom_peaks)
        intervals[chrom] = spans
    return intervals


def _to_numpy_local(x):
    if mx is not None:
        try:
            return np.asarray(mx.asarray(x))
        except Exception:
            pass
    try:
        return np.asarray(x)
    except Exception:
        return np.array(x)


def snapshot_scores(track):
    """Return {chrom: (pos, score)} as NumPy copies to avoid later mutation."""
    snap: Dict[bytes, Tuple[np.ndarray, np.ndarray]] = {}
    for chrom in track.get_chr_names():
        pos, _, _, score = track.get_data_by_chr(chrom)
        snap[chrom] = (_to_numpy_local(pos).copy(), _to_numpy_local(score).copy())
    return snap


def _interval_overlap(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    overlap = 0
    i = j = 0
    while i < len(a) and j < len(b):
        a_start, a_end = a[i]
        b_start, b_end = b[j]
        start = max(a_start, b_start)
        end = min(a_end, b_end)
        if start < end:
            overlap += end - start
        if a_end <= b_end:
            i += 1
        else:
            j += 1
    return overlap


def _count_overlapping(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    """Count how many intervals in ``a`` touch any interval in ``b``."""
    count = 0
    j = 0
    for a_start, a_end in a:
        while j < len(b) and b[j][1] <= a_start:
            j += 1
        k = j
        while k < len(b) and b[k][0] < a_end:
            if a_start < b[k][1] and a_end > b[k][0]:
                count += 1
                break
            k += 1
    return count


def _first_nonoverlapping(source: Dict[bytes, List[Tuple[int, int]]],
                          target: Dict[bytes, List[Tuple[int, int]]],
                          limit: int = 3) -> List[Tuple[str, int, int]]:
    examples: List[Tuple[str, int, int]] = []
    for chrom, intervals in source.items():
        other = target.get(chrom, [])
        j = 0
        for start, end in intervals:
            while j < len(other) and other[j][1] <= start:
                j += 1
            k = j
            found = False
            while k < len(other) and other[k][0] < end:
                if start < other[k][1] and end > other[k][0]:
                    found = True
                    break
                k += 1
            if not found:
                chrom_str = chrom.decode() if isinstance(chrom, (bytes, bytearray)) else str(chrom)
                examples.append((chrom_str, start, end))
                if len(examples) >= limit:
                    return examples
    return examples


def summarize_peak_consistency(cpu_peaks, mlx_peaks):
    cpu_intervals = peaks_to_intervals(cpu_peaks)
    mlx_intervals = peaks_to_intervals(mlx_peaks)

    cpu_bp = sum(end - start for spans in cpu_intervals.values() for start, end in spans)
    mlx_bp = sum(end - start for spans in mlx_intervals.values() for start, end in spans)
    cpu_total = sum(len(spans) for spans in cpu_intervals.values())
    mlx_total = sum(len(spans) for spans in mlx_intervals.values())

    overlap_bp = 0
    cpu_hits = 0
    mlx_hits = 0
    for chrom in set(cpu_intervals) | set(mlx_intervals):
        cpu_spans = cpu_intervals.get(chrom, [])
        mlx_spans = mlx_intervals.get(chrom, [])
        overlap_bp += _interval_overlap(cpu_spans, mlx_spans)
        cpu_hits += _count_overlapping(cpu_spans, mlx_spans)
        mlx_hits += _count_overlapping(mlx_spans, cpu_spans)

    union_bp = cpu_bp + mlx_bp - overlap_bp
    jaccard = overlap_bp / union_bp if union_bp else 0.0

    return {
        "cpu_total": cpu_total,
        "mlx_total": mlx_total,
        "cpu_bp": cpu_bp,
        "mlx_bp": mlx_bp,
        "overlap_bp": overlap_bp,
        "jaccard": jaccard,
        "cpu_hits": cpu_hits,
        "mlx_hits": mlx_hits,
        "cpu_examples": _first_nonoverlapping(cpu_intervals, mlx_intervals),
        "mlx_examples": _first_nonoverlapping(mlx_intervals, cpu_intervals),
    }


from typing import Optional

def summarize_score_differences(cpu_scores, mlx_scores, label: str, cutoff: Optional[float] = None):
    """Compare per-bin scores (p or q) between CPU and MLX snapshots."""
    diffs: List[np.ndarray] = []
    coverage_cpu = coverage_mlx = overlap = 0
    examples: List[Tuple[float, str, int, float, float]] = []
    mismatched: List[str] = []

    for chrom in sorted(set(cpu_scores) | set(mlx_scores)):
        if chrom not in cpu_scores or chrom not in mlx_scores:
            missing = "CPU" if chrom not in cpu_scores else "MLX"
            mismatched.append(f"{chrom!r} missing on {missing}")
            continue
        cpos, cscore = cpu_scores[chrom]
        mpos, mscore = mlx_scores[chrom]
        if cpos.shape != mpos.shape or not np.array_equal(cpos, mpos):
            mismatched.append(f"{chrom!r} positions differ (CPU {cpos.shape}, MLX {mpos.shape})")
            continue
        diff = np.abs(cscore - mscore)
        diffs.append(diff)
        if diff.size:
            top_idx = np.argpartition(diff, -3)[-3:]
            for idx in top_idx:
                delta = float(diff[idx])
                examples.append((delta,
                                 chrom.decode() if isinstance(chrom, (bytes, bytearray)) else str(chrom),
                                 int(cpos[idx]),
                                 float(cscore[idx]),
                                 float(mscore[idx])))
        if cutoff is not None:
            cpu_mask = cscore >= cutoff
            mlx_mask = mscore >= cutoff
            coverage_cpu += int(cpu_mask.sum())
            coverage_mlx += int(mlx_mask.sum())
            overlap += int(np.logical_and(cpu_mask, mlx_mask).sum())

    if diffs:
        all_diff = np.concatenate(diffs)
        stats = {
            "mean": float(np.mean(all_diff)),
            "median": float(np.median(all_diff)),
            "p99": float(np.percentile(all_diff, 99)),
            "max": float(np.max(all_diff)),
            "n": int(all_diff.size),
        }
    else:
        stats = {"mean": 0.0, "median": 0.0, "p99": 0.0, "max": 0.0, "n": 0}

    examples.sort(key=lambda x: x[0], reverse=True)
    examples = examples[:3]
    return {
        "label": label,
        "stats": stats,
        "examples": examples,
        "mismatched": mismatched,
        "coverage_cpu": coverage_cpu,
        "coverage_mlx": coverage_mlx,
        "overlap": overlap,
    }


def benchmark(num_bins: int, bin_size: int, cutoff: float, approx_pvalue: bool):
    chrom = b"chr1"
    pos, chip, control = build_synthetic_data(num_bins, bin_size)
    # CPU path
    cpu_track = populate_scoretrack(ScoreTrackII(30.0, 30.0), pos, chip, control, chrom)
    start = time.perf_counter()
    cpu_track.change_score_method(ord('p'))
    cpu_p_time = time.perf_counter() - start
    cpu_p_scores = snapshot_scores(cpu_track)
    start = time.perf_counter()
    cpu_track.change_score_method(ord('q'))
    cpu_q_time = time.perf_counter() - start
    cpu_q_scores = snapshot_scores(cpu_track)
    start = time.perf_counter()
    cpu_peaks = cpu_track.call_peaks(cutoff=cutoff)
    cpu_peak_time = time.perf_counter() - start

    # MLX path
    gpu_track = populate_scoretrack(ScoreTrackMLX(30.0, 30.0, approx_pvalue=approx_pvalue), pos, chip, control, chrom)
    start = time.perf_counter()
    gpu_track.compute_pvalue()
    gpu_p_time = time.perf_counter() - start
    gpu_p_scores = snapshot_scores(gpu_track)
    start = time.perf_counter()
    gpu_track.compute_qvalue()
    gpu_q_time = time.perf_counter() - start
    gpu_q_scores = snapshot_scores(gpu_track)
    start = time.perf_counter()
    gpu_peaks = gpu_track.call_peaks(cutoff=cutoff)
    gpu_peak_time = time.perf_counter() - start

    return {
        "cpu_p": cpu_p_time,
        "cpu_q": cpu_q_time,
        "cpu_call": cpu_peak_time,
        "gpu_p": gpu_p_time,
        "gpu_q": gpu_q_time,
        "gpu_call": gpu_peak_time,
        "cpu_peaks": cpu_peaks,
        "gpu_peaks": gpu_peaks,
        "cpu_p_scores": cpu_p_scores,
        "cpu_q_scores": cpu_q_scores,
        "gpu_p_scores": gpu_p_scores,
        "gpu_q_scores": gpu_q_scores,
        "pos": pos,
        "chip": chip,
        "control": control,
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
    parser.add_argument(
        "--mlx-pvalue",
        choices=["approx", "exact"],
        default="approx",
        help="Use approximate normal-tail p-values on MLX (fast, default) or exact Poisson on CPU fallback.",
    )
    parser.add_argument(
        "--top-diffs",
        type=int,
        default=3,
        help="Show up to this many largest score differences with pileup context.",
    )
    args = parser.parse_args()

    if mx is not None and args.mlx == "cpu":
        mx.set_default_device(mx.cpu)
    print(f"MLX default device: {mx.default_device() if mx is not None else 'Unavailable'}")

    if ScoreTrackMLX is None:
        raise SystemExit("ScoreTrackMLX is unavailable; ensure MACS3.Signal.ScoreTrackMLX can be imported.")

    approx_pvalue = args.mlx_pvalue == "approx"
    results = benchmark(args.bins, args.bin_size, args.cutoff, approx_pvalue=approx_pvalue)
    device_label = "MLX tensor" if mx is not None else "NumPy fallback"
    print(f"CPU p-value: {results['cpu_p']:.3f}s | q-value: {results['cpu_q']:.3f}s | call_peaks: {results['cpu_call']:.3f}s")
    print(f"{device_label} p-value: {results['gpu_p']:.3f}s | q-value: {results['gpu_q']:.3f}s | call_peaks: {results['gpu_call']:.3f}s")
    total_cpu = results['cpu_p'] + results['cpu_q'] + results['cpu_call']
    total_gpu = results['gpu_p'] + results['gpu_q'] + results['gpu_call']
    print(f"Total CPU: {total_cpu:.3f}s | Total {device_label}: {total_gpu:.3f}s")
    if total_gpu > 0:
        print(f"Overall speedup (CPU/GPU): {total_cpu / total_gpu:.2f}x")
    consistency = summarize_peak_consistency(results["cpu_peaks"], results["gpu_peaks"])
    cpu_overlap_pct = (consistency["cpu_hits"] / consistency["cpu_total"] * 100) if consistency["cpu_total"] else 0.0
    mlx_overlap_pct = (consistency["mlx_hits"] / consistency["mlx_total"] * 100) if consistency["mlx_total"] else 0.0
    print("Peak consistency (CPU vs MLX):")
    print(f"- Peaks called: CPU {consistency['cpu_total']} | MLX {consistency['mlx_total']}")
    print(f"- Basepair coverage: CPU {consistency['cpu_bp']} bp | MLX {consistency['mlx_bp']} bp | Overlap {consistency['overlap_bp']} bp | Jaccard {consistency['jaccard']:.3f}")
    print(f"- Overlap counts: CPU {consistency['cpu_hits']}/{consistency['cpu_total']} ({cpu_overlap_pct:.1f}%) | MLX {consistency['mlx_hits']}/{consistency['mlx_total']} ({mlx_overlap_pct:.1f}%)")
    if consistency["cpu_examples"]:
        preview = ", ".join(f"{c}:{s}-{e}" for c, s, e in consistency["cpu_examples"])
        print(f"  CPU-only examples (first {len(consistency['cpu_examples'])}): {preview}")
    if consistency["mlx_examples"]:
        preview = ", ".join(f"{c}:{s}-{e}" for c, s, e in consistency["mlx_examples"])
        print(f"  MLX-only examples (first {len(consistency['mlx_examples'])}): {preview}")
    # Score-level deltas to help diagnose approximate p-values
    q_diff = summarize_score_differences(results["cpu_q_scores"], results["gpu_q_scores"], label="q", cutoff=args.cutoff)
    p_diff = summarize_score_differences(results["cpu_p_scores"], results["gpu_p_scores"], label="p")
    print(f"Score deltas (-log10 {q_diff['label']}): mean {q_diff['stats']['mean']:.4f}, median {q_diff['stats']['median']:.4f}, "
          f"p99 {q_diff['stats']['p99']:.4f}, max {q_diff['stats']['max']:.4f} across {q_diff['stats']['n']} bins")
    if q_diff["coverage_cpu"] or q_diff["coverage_mlx"]:
        print(f"  Bins >= cutoff {args.cutoff}: CPU {q_diff['coverage_cpu']} | MLX {q_diff['coverage_mlx']} | overlap {q_diff['overlap']}")
    if q_diff["examples"]:
        preview = ", ".join(f"{c}:{pos} CPU {cs:.2f} vs MLX {gs:.2f} (Δ{d:.2f})"
                            for d, c, pos, cs, gs in q_diff["examples"])
        print(f"  Largest q-score diffs: {preview}")
    if p_diff["examples"]:
        preview = ", ".join(f"{c}:{pos} CPU {cs:.2f} vs MLX {gs:.2f} (Δ{d:.2f})"
                            for d, c, pos, cs, gs in p_diff["examples"])
        print(f"  Largest p-score diffs: {preview}")
    if p_diff["mismatched"] or q_diff["mismatched"]:
        print("  Position mismatches found:", "; ".join(p_diff["mismatched"] + q_diff["mismatched"]))
    if args.top_diffs > 0 and args.mlx_pvalue == "approx":
        pos = results["pos"]
        chip = results["chip"]
        control = results["control"]
        pos_to_idx = {int(p): i for i, p in enumerate(pos)}
        def show_examples(examples, label: str):
            rows = []
            for d, chrom, position, cpu_score, mlx_score in examples[:args.top_diffs]:
                idx = pos_to_idx.get(int(position))
                if idx is None:
                    continue
                rows.append(
                    f"{chrom}:{position} chip {chip[idx]:.2f} ctrl {control[idx]:.2f} "
                    f"CPU {label} {cpu_score:.2f} vs MLX {label} {mlx_score:.2f} (Δ{d:.2f})"
                )
            if rows:
                print(f"  {label.upper()} mismatch context:")
                for line in rows:
                    print(f"    - {line}")
        show_examples(q_diff["examples"], "q")
        show_examples(p_diff["examples"], "p")


if __name__ == "__main__":
    main()
