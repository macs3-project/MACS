# MACS Version Speed and Memory Survey

This developer benchmark compares MACS v1, MACS2, and MACS3 on the same
5M-read CTCF BED dataset. It is a pragmatic software-performance survey, not a
formal methods benchmark and not a claim that biological outputs should be
identical across major versions.

Use the manual GitHub Actions workflow for the recommended run. It benchmarks
all versions on one GitHub-hosted x86_64 Ubuntu runner, avoiding Apple Silicon
Rosetta overhead and avoiding Docker.

The benchmark uses these refs by default:

- MACS v1: `origin/macs_v1`
- MACS2: `origin/macs_v2`
- MACS3: the current workflow checkout, usually the current branch or PR head

The primary task is:

- treatment: `test/CTCF_12878_5M.bed.gz`
- control: `test/Input_12878_5M.bed.gz`
- format: `BED`
- genome: `hs`
- repeats: `3` by default

## GitHub Actions Run

Open the repository's **Actions** tab, choose **MACS Version Benchmark**, and
run the workflow manually. For a first smoke test, use `repeats=1`. For the
reported survey, use `repeats=3` or higher.

The workflow:

- runs on `ubuntu-22.04` x86_64
- installs GNU `time`
- creates isolated conda environments for MACS v1, MACS2, and MACS3
- installs each version from its matching worktree, including required
  submodules for MACS3's bundled fermi-lite/SIMDe code
- runs `scripts/run_macs_version_survey.sh` with `TIME_MODE=linux`
- uploads a `macs-version-benchmark-results` artifact

The workflow uses:

```bash
MACS1_CMD="conda run -n macs-v1-py2 macs"
MACS2_CMD="conda run -n macs-v2-survey macs2"
MACS3_CMD="conda run -n macs-v3-survey macs3"
```

MACS v1 installs the historical command as `macs`, not `macs14`, on the
`origin/macs_v1` branch.

## Outputs

The uploaded artifact contains:

- `results/benchmark_runs.tsv`: one row per version and repeat
- `results/summary.tsv`: median/min/max wall time and median peak RSS
- `results/output_peak_counts.tsv`: row counts for generated peak files
- `logs/git_revisions.txt`: exact commit and describe output for each ref
- `logs/input_sha256.txt`: input checksums
- `logs/version_checks.txt`: command/version checks
- `logs/*.stdout.txt` and `logs/*.time_stderr.txt`: per-run command logs

On Linux, peak resident memory is parsed from `/usr/bin/time -v` as `Maximum
resident set size (kbytes)` and converted to bytes. On macOS, the same harness
can parse `/usr/bin/time -l` as `maximum resident set size`, already in bytes.

To regenerate summaries from an artifact or local result directory:

```bash
python scripts/summarize_macs_version_survey.py /path/to/macs-version-survey
```

## Local Caveat

Do not use an Apple Silicon local run for the main comparison if MACS v1 must
run under Rosetta. Native `osx-arm64` conda does not provide Python 2.7 from the
usual channels, while Rosetta x86_64 would add translation overhead and make the
comparison less fair. The GitHub Actions x86_64 workflow avoids that problem.

Local macOS runs remain possible for harness debugging:

```bash
ROOT=/path/to/benchmark/root \
REPEATS=1 \
TIME_MODE=macos \
MACS1_CMD="macs" \
MACS2_CMD="macs2" \
MACS3_CMD="macs3" \
bash scripts/run_macs_version_survey.sh
```

## Reporting

Use a table like this when reporting results:

| Version | Source ref | Python | Command | Repeats | Median wall time (s) | Median peak RSS (MB) | Notes |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| MACS v1 | `origin/macs_v1` commit | Python 2.7 | `macs ...` | 3 | TBD | TBD | GitHub Actions x86_64 |
| MACS2 | `origin/macs_v2` commit | Python 3.9 | `macs2 callpeak ...` | 3 | TBD | TBD | GitHub Actions x86_64 |
| MACS3 | current branch/PR head | Python 3.11 | `macs3 callpeak ...` | 3 | TBD | TBD | GitHub Actions x86_64 |

Do not require identical peak counts across MACS v1, MACS2, and MACS3. Defaults
and algorithms differ across major versions; this survey targets computational
behavior on the same input files.
