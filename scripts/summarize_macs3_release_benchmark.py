#!/usr/bin/env python3
"""Summarize a current-versus-baseline MACS3 release benchmark."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


PEAK_SUFFIXES = (
    "_peaks.xls",
    "_peaks.narrowPeak",
    "_peaks.broadPeak",
    "_peaks.gappedPeak",
)


def read_runs(path: Path) -> tuple[dict[str, list[tuple[float, int]]], dict[str, str]]:
    groups: dict[str, list[tuple[float, int]]] = defaultdict(list)
    refs: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            revision = row["revision"]
            refs[revision] = row["ref"]
            if row["status"] != "0":
                continue
            if not row["seconds_wall"] or not row["max_rss_bytes"]:
                continue
            groups[revision].append(
                (float(row["seconds_wall"]), int(row["max_rss_bytes"]))
            )
    return groups, refs


def calculate_summary(
    groups: dict[str, list[tuple[float, int]]], refs: dict[str, str]
) -> dict[str, dict[str, float | int | str]]:
    summary: dict[str, dict[str, float | int | str]] = {}
    for revision in ("baseline", "current"):
        values = groups.get(revision, [])
        if not values:
            continue
        wall = [value[0] for value in values]
        rss_mb = [value[1] / 1024 / 1024 for value in values]
        summary[revision] = {
            "ref": refs.get(revision, ""),
            "n": len(values),
            "median_wall_s": statistics.median(wall),
            "min_wall_s": min(wall),
            "max_wall_s": max(wall),
            "median_peak_rss_mb": statistics.median(rss_mb),
        }
    return summary


def write_summary(summary: dict[str, dict[str, float | int | str]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "revision",
                "ref",
                "n",
                "median_wall_s",
                "min_wall_s",
                "max_wall_s",
                "median_peak_rss_mb",
            ]
        )
        for revision in ("baseline", "current"):
            row = summary.get(revision)
            if row is None:
                continue
            writer.writerow(
                [
                    revision,
                    row["ref"],
                    row["n"],
                    f"{row['median_wall_s']:.3f}",
                    f"{row['min_wall_s']:.3f}",
                    f"{row['max_wall_s']:.3f}",
                    f"{row['median_peak_rss_mb']:.1f}",
                ]
            )


def comparison_values(
    summary: dict[str, dict[str, float | int | str]]
) -> dict[str, float] | None:
    if "baseline" not in summary or "current" not in summary:
        return None
    baseline = summary["baseline"]
    current = summary["current"]
    baseline_wall = float(baseline["median_wall_s"])
    current_wall = float(current["median_wall_s"])
    baseline_rss = float(baseline["median_peak_rss_mb"])
    current_rss = float(current["median_peak_rss_mb"])
    return {
        "wall_ratio": current_wall / baseline_wall,
        "speedup": baseline_wall / current_wall,
        "wall_change_pct": (current_wall / baseline_wall - 1) * 100,
        "rss_ratio": current_rss / baseline_rss,
        "rss_change_pct": (current_rss / baseline_rss - 1) * 100,
    }


def write_comparison_tsv(values: dict[str, float] | None, path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "current_vs_baseline_wall_ratio",
                "baseline_over_current_speedup",
                "wall_change_pct",
                "current_vs_baseline_rss_ratio",
                "rss_change_pct",
            ]
        )
        if values is not None:
            writer.writerow(
                [
                    f"{values['wall_ratio']:.3f}",
                    f"{values['speedup']:.3f}",
                    f"{values['wall_change_pct']:+.1f}",
                    f"{values['rss_ratio']:.3f}",
                    f"{values['rss_change_pct']:+.1f}",
                ]
            )


def write_markdown(
    summary: dict[str, dict[str, float | int | str]],
    values: dict[str, float] | None,
    path: Path,
) -> None:
    lines = [
        "## MACS3 release benchmark",
        "",
        "| Revision | Ref | Runs | Median wall (s) | Min–max (s) | Median RSS (MB) |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for revision in ("baseline", "current"):
        row = summary.get(revision)
        if row is None:
            continue
        lines.append(
            f"| {revision} | `{row['ref']}` | {row['n']} | "
            f"{row['median_wall_s']:.3f} | {row['min_wall_s']:.3f}–"
            f"{row['max_wall_s']:.3f} | {row['median_peak_rss_mb']:.1f} |"
        )
    lines.append("")
    if values is None:
        lines.append("A comparison could not be calculated because one revision had no successful runs.")
    else:
        lines.extend(
            [
                f"Current/baseline wall-time ratio: **{values['wall_ratio']:.3f}** "
                f"({values['wall_change_pct']:+.1f}%; speedup "
                f"**{values['speedup']:.3f}×**).",
                "",
                f"Current/baseline peak-RSS ratio: **{values['rss_ratio']:.3f}** "
                f"({values['rss_change_pct']:+.1f}%).",
            ]
        )
    path.write_text("\n".join(lines) + "\n")


def count_data_rows(path: Path) -> int:
    with path.open(errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def write_peak_counts(results_dir: Path, path: Path) -> None:
    peak_files = sorted(
        candidate
        for candidate in results_dir.rglob("*")
        if candidate.is_file() and candidate.name.endswith(PEAK_SUFFIXES)
    )
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "rows"])
        for peak_file in peak_files:
            writer.writerow(
                [str(peak_file.relative_to(results_dir.parent)), count_data_rows(peak_file)]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        nargs="?",
        default=str(Path.home() / "benchmarks" / "macs3-release-benchmark"),
        help="Benchmark root containing results/benchmark_runs.tsv.",
    )
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    results_dir = root / "results"
    runs_tsv = results_dir / "benchmark_runs.tsv"
    if not runs_tsv.exists():
        raise SystemExit(f"Missing benchmark run table: {runs_tsv}")

    groups, refs = read_runs(runs_tsv)
    summary = calculate_summary(groups, refs)
    values = comparison_values(summary)
    write_summary(summary, results_dir / "summary.tsv")
    write_comparison_tsv(values, results_dir / "comparison.tsv")
    write_markdown(summary, values, results_dir / "comparison.md")
    write_peak_counts(results_dir, results_dir / "output_peak_counts.tsv")


if __name__ == "__main__":
    main()

