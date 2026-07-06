#!/usr/bin/env python3
"""Summarize MACS version speed/memory survey outputs."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


PEAK_SUFFIXES = ("_peaks.xls", "_peaks.narrowPeak", "_peaks.broadPeak", "_peaks.gappedPeak")


def read_successful_runs(path: Path) -> dict[str, list[tuple[float, int]]]:
    groups: dict[str, list[tuple[float, int]]] = defaultdict(list)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("status") != "0":
                continue
            wall = row.get("seconds_wall", "")
            rss = row.get("max_rss_bytes", "")
            if not wall or not rss:
                continue
            groups[row["version"]].append((float(wall), int(rss)))
    return groups


def write_summary(groups: dict[str, list[tuple[float, int]]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "version",
                "n",
                "median_wall_s",
                "min_wall_s",
                "max_wall_s",
                "median_peak_rss_mb",
            ]
        )
        for version in sorted(groups):
            vals = groups[version]
            wall = [item[0] for item in vals]
            rss_mb = [item[1] / 1024 / 1024 for item in vals]
            writer.writerow(
                [
                    version,
                    len(vals),
                    f"{statistics.median(wall):.3f}",
                    f"{min(wall):.3f}",
                    f"{max(wall):.3f}",
                    f"{statistics.median(rss_mb):.1f}",
                ]
            )


def count_data_rows(path: Path) -> int:
    count = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            count += 1
    return count


def peak_files(results_dir: Path) -> list[Path]:
    files = []
    for path in results_dir.rglob("*"):
        if path.is_file() and path.name.endswith(PEAK_SUFFIXES):
            files.append(path)
    return sorted(files)


def write_peak_counts(results_dir: Path, path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "rows"])
        for peak_file in peak_files(results_dir):
            writer.writerow([str(peak_file.relative_to(results_dir.parent)), count_data_rows(peak_file)])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        nargs="?",
        default=str(Path.home() / "benchmarks" / "macs-version-survey"),
        help="Benchmark root containing results/benchmark_runs.tsv.",
    )
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    results_dir = root / "results"
    runs_tsv = results_dir / "benchmark_runs.tsv"
    if not runs_tsv.exists():
        raise SystemExit(f"Missing benchmark run table: {runs_tsv}")

    groups = read_successful_runs(runs_tsv)
    write_summary(groups, results_dir / "summary.tsv")
    write_peak_counts(results_dir, results_dir / "output_peak_counts.tsv")


if __name__ == "__main__":
    main()
