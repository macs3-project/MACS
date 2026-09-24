# MACS performance benchmarks

This page reports measured `callpeak` performance on one reference dataset.
These results describe this workload on GitHub-hosted Linux runners; they are
not a guarantee of the same speed or memory use on other datasets or machines.

## MACS 3.0.5 versus 3.0.4

The [September 24, 2026 release benchmark](https://github.com/macs3-project/MACS/actions/runs/35957533788)
compared the `v3.0.4` tag (`e8ff040`) with the 3.0.5 release candidate on
`main` (`e0810f7`). Both revisions report their expected version strings.

| Revision | Median wall time | Range | Median peak memory | Repeats |
| --- | ---: | ---: | ---: | ---: |
| MACS 3.0.4 (`v3.0.4`) | 24.59 s | 24.20–24.73 s | 288.9 MB | 3 |
| MACS 3.0.5 (`e0810f7`) | 17.75 s | 17.62–17.95 s | 292.8 MB | 3 |

For this dataset, MACS 3.0.5 took **27.8% less wall time**
(**1.385× speedup**) and used similar peak memory (+1.4%). All
six runs produced 36,411 `narrowPeak` rows; matching row counts do not imply
identical peak calls. The [run artifact](https://github.com/macs3-project/MACS/actions/runs/35957533788)
contains the individual timings, output counts, input checksums, revision IDs,
and command logs.

The benchmark runs `macs3 callpeak` on the 5M-read CTCF treatment and input
BED files in `test/`, with `-f BED -g hs -q 0.01`, on an Ubuntu 22.04 x86-64
runner. Each revision is installed in a separate, cloned Conda environment
with Python 3.12 and matched dependencies. Execution order alternates between
repeats; wall time and maximum resident memory come from GNU `time`.

Because GitHub-hosted runner performance can vary, the speedup compares the
two revisions **within this run**; timings from separate runs should not be
compared directly.

## MACS v1, MACS2, and MACS3 survey

The [September 24, 2026 major-version benchmark](https://github.com/macs3-project/MACS/actions/runs/35957562638)
ran MACS v1, MACS2, and the 3.0.5 release candidate on the same CTCF and
input files. Each version ran three times on one Ubuntu 22.04 x86-64 runner.

| Version and source commit | Median wall time | Range | Median peak memory | Repeats |
| --- | ---: | ---: | ---: | ---: |
| MACS v1 (`a662072`) | 67.34 s | 65.77–68.06 s | 367.8 MB | 3 |
| MACS2 (`b18703b`) | 31.87 s | 31.74–31.99 s | 298.1 MB | 3 |
| MACS3 3.0.5 (`e0810f7`) | 32.12 s | 32.02–32.37 s | 298.2 MB | 3 |

MACS v1 used Python 2.7, MACS2 used Python 3.9, and MACS3 used Python 3.12.
The major versions also differ in dependencies, defaults, and peak-calling
behavior. They produced different peak-output row counts in this run, so the
timings are a software-performance survey rather than a
like-for-like algorithmic speedup. In particular, MACS2 and MACS3 had similar
wall time and peak memory on this runner; no major-version speedup is claimed.

This survey used a different runner from the release comparison above; its
32.12 s MACS3 time should not be compared with the 17.75 s figure from that
separate run. The [survey artifact](https://github.com/macs3-project/MACS/actions/runs/35957562638)
contains the exact commands, revision IDs, raw timing logs, and output counts.

## Reproducing the reports

The two manual GitHub Actions workflows are [MACS3 Release Benchmark](https://github.com/macs3-project/MACS/actions/workflows/macs3-release-benchmark.yml)
and [MACS Version Benchmark](https://github.com/macs3-project/MACS/actions/workflows/macs-version-benchmark.yml).
Run them from the revision of interest with three measured repeats. The release
comparison defaults to `v3.0.4` as its baseline. Both workflows publish their
raw results as downloadable artifacts; the release workflow also displays its
comparison in the run summary.
