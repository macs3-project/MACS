# Benchmarking performance defaults

This is the evidence protocol for [#773](https://github.com/macs3-project/MACS/issues/773).
The inventory in [performance defaults](performance_defaults.md) describes
policy, not optimal values. This change does not retune buffers, task thresholds,
cache reserves, huge-page policy, prefetch distances or vendor constants.

## Revisions and matched builds

Compare three compiled revisions, on the same host and dependencies:

| Label | Revision | Purpose |
| --- | --- | --- |
| baseline | `170556de25475fc59b1e219e898c21b7eecd5dd8` | Original merged PR #756 implementation |
| extraction | `324f779a270910c8befda97443912b1175704823` | #767 default-preserving extraction, including Cython 3.0 layout correction |
| current | `c1287de2d9b955267aea426cb7112e7dbd2b224f` | Implementation after #772; benchmark-only changes follow this commit |

The historical builds do not have `--threads` or `--process-backend`.
The harness records their **legacy default** runs, without pretending that an
environment variable or a patched constant supplies equivalent CPU controls.
On Linux, apply the same CPU affinity/address-space profile to all revisions.
Legacy fixed pool policies may still oversubscribe that allocation. Current
budgets include the caller, and process workers use at most budget minus one.
Backend diagnostics identify small-workload and host safety fallbacks. On
macOS the current automatic peak-process backend is intentionally serial;
native parser/track threads can still run. See [execution policy](execution.md).

The manual `.github/workflows/performance-defaults.yml` workflow builds
baseline, extraction and current wheels from source archives and the exact
SIMDe submodule commit, then installs each into its own environment. Dependencies
are resolved once and shared, avoiding different NumPy/Cython versions between
revisions. It runs Linux x86_64 and macOS x86_64/arm64, alternates revision
order between repetitions, and uploads raw outputs, commands, environment
metadata and build logs. Failed/mismatched results fail the workflow.

## Reproduce locally

Activate `macs3dev`; use `uv run` for tooling. `psutil` is a benchmark runner
dependency, not a MACS3 runtime dependency. Build each revision with the same
Python, NumPy, Cython, compiler and dependencies, using the workflow's archive
and wheel steps. Include the submodule: a bare `git archive` omits its contents.
Preserve build logs and wheel checksums alongside the results.

```sh
conda activate macs3dev
uv run --with build python -m build --wheel --no-isolation
uv run python scripts/benchmark_performance.py prepare --out /tmp/macs3-perf/data

# Installed wheel: choose its interpreter; no checkout PYTHONPATH.
uv run python scripts/benchmark_performance.py run \
  --manifest /tmp/macs3-perf/data/manifest.json \
  --python /path/to/current-env/bin/python \
  --label current --revision c1287de2d9b955267aea426cb7112e7dbd2b224f \
  --budgets 1,2,4,8 --repeats 3 --out /tmp/macs3-perf/current

# For an in-place build, select both package and matching CLI explicitly.
uv run python scripts/benchmark_performance.py run \
  --manifest /tmp/macs3-perf/data/manifest.json \
  --package-root "$PWD" --cli "$PWD/bin/macs3" \
  --label current-inplace --revision c1287de2d9b955267aea426cb7112e7dbd2b224f \
  --budgets 1,4 --repeats 1 --out /tmp/macs3-perf/current-inplace

# Repeat run with baseline and extraction interpreters and matching SHAs.
uv run python scripts/benchmark_performance.py compare \
  /tmp/macs3-perf/baseline /tmp/macs3-perf/extraction /tmp/macs3-perf/current \
  --out /tmp/macs3-perf/comparison.md
```

Set `PYTHONPATH` only when explicitly choosing an in-place installation.
Target subprocesses execute from their result directories, using the selected
interpreter. Metadata refuses pure-Python substitutes for native modules and
records the paths and SHA256 of the actual compiled extensions. Check the
recorded package path and CLI before interpreting a comparison. Each result
directory is new; existing evidence is never overwritten by `run`.

The two historical revisions have a pre-existing `PETrackII.destroy()` bug
(`dict` has no `resize`). The isolated allocator diagnostic records that exact
exception in `allocator_cleanup_notes` and then releases memory at process
exit. It does not change the historical source or suppress CLI failures.

For Linux constraints, add `--profile limited --affinity 0,1
--address-space-mib 2048`, replacing CPU numbers with two CPUs from
`os.sched_getaffinity(0)`. This is a per-process **virtual address-space** limit,
not a cap on RSS or the process tree. For a real shared memory ceiling, run the
same command inside a VM/container/cgroup with a recorded memory limit; retain
`memory.max`, `memory.peak`, `memory.events`, CPU quota and affinity information.
Do not describe an RLIMIT_AS run as a cgroup memory benchmark. macOS requires a
limited VM for this profile. Resource failures are evidence, not successful
timings. The default per-target timeout is 180 seconds, with process-group kill
and child reaping on timeout. Increase it explicitly for biological datasets.

## Datasets and runtime

`prepare` creates BGZF BAM and BGZF FRAG files from seed 773. Each BAM contains
complete 36 bp read pairs for 150 bp fragments; FRAG contains the corresponding
intervals, count one and 8,192 rotating barcode keys. A seeded dispersed
background and sixteen enriched regions per contig produce nonempty peaks.
Treatment and control have independent positions. The inputs are sorted;
each compressed file has a SHA256 recorded and verified before measurement.

| Cases | Fragments per treatment/control | Contigs | Purpose |
| --- | ---: | ---: | --- |
| small-bam / small-bampe / small-frag | 6,000 | 3 | Startup and serial fallback |
| large-bam / large-bampe / large-frag | 500,000 | 3 | Cross parser/track/process thresholds |
| many-bam / many-bampe / many-frag | 500,000 | 256 | Many small contigs and scheduling/storage overhead |

Use `--small`, `--large`, `--contigs` to scale generation. BAM/BAMPE share the
same input bytes but exercise different APIs. All revisions receive the same
`callpeak -t ... -c ... -f ... -g CONTIGS*2000000 --keep-dup all --nomodel
--extsize 150 -B -n bench` options. Model fitting is intentionally excluded;
this is parser/track/scoring/caller evidence, not a model-fitting benchmark.
Default allocation buffers remain active. External biological inputs can be
specified with the same manifest schema (`name`, `format`, `treat`, `control`,
`genome_size`, `sha256`); record origin, accession, license, preprocessing and
runtime in the report. Generated data does not establish biological workload
representativeness. Never commit large generated inputs.

The small integration test uses 900/3,000 fragments and twelve contigs, runs
nine CLI jobs plus the full cache probe, and takes about 8 seconds on the local
Intel Mac. Full-matrix expected runtime and individual observed timings are
recorded in the [initial report](performance_benchmark_report.md). The workflow
allows 90 minutes including three native builds; individual run timeouts remain
bounded. Synthetic scale and runtime are part of the evidence, not universal
performance expectations.

## Measurements and equality gates

Every run preserves its specification, command, merged CLI log, exit status,
timeout status, sampling errors, timings, outputs and scientific output hashes.
`metadata.json` includes Python/dependency versions, CPU/platform, compiler
configuration, relevant environment settings, resource profile and the expected
run matrix. Compiler configuration describes Python's build settings; the
actual extension flags must also be taken from the retained wheel build log.
Third-party OpenMP/BLAS pools are held at one thread and hash seed at zero.
There are no unreported warm-up exclusions. Use at least three repetitions,
alternate revision order, avoid simultaneous builds/tests, and report medians
and ranges from one host/run. Cache and allocator diagnostic probes are outside
the measured CLI jobs.

* Wall time includes Python/CLI startup, output writing, cleanup and monitoring
  overhead. The default 20 ms polling interval adds up to one polling interval
  of exit-detection delay. Keep it identical across compared revisions.
* Sampled process-tree RSS sums resident bytes across live descendants. Fork
  shared pages are counted repeatedly: this is not PSS or unique physical
  memory. Separately report OS high-water RSS for the caller and the largest
  waited-for child; neither is a simultaneous tree peak.
* Temp storage is sampled logical file length under the private `TMPDIR`.
  Linux also follows `/proc/PID/fd` for unlinked open files, deduplicating
  device/inode pairs across processes. Short-lived peaks between samples can
  be missed. On macOS deleted descriptors are not enumerated by psutil, so the
  value is explicitly a **partial lower bound**. It cannot close the full temp
  storage acceptance criterion. Directory scans alone are not an exact peak.
* After CLI exit, retained temp bytes or sampling permission failures on live
  processes fail the run. Processes that have exited or entered the Linux
  kernel exit path can lose procfs access before reaping; these expected exit
  races are ignored. Native module failure is not converted into a skip or
  pure-Python run.

`compare` rejects missing/duplicate planned runs, input-definition differences,
failed invocations, cache retention regressions, missing/empty peak artifacts,
and scientific byte differences. It compares peak tables, summit BED, treatment
and control bedGraphs; XLS comments contain volatile invocation/version metadata
and are the only excluded bytes. Row-count equality alone is insufficient.
Current budget/backend outputs must match. Historical differences caused by
intervening scientific fixes must be investigated with a reproducer and exact
affected hashes, and documented separately; do not relax the gate or silently
accept a mismatch. These generated cases cover narrow peaks; broad/summit-mode
equivalence remains covered by `test_native_cli.py`.

## Log cache and repeated jobs

In a fresh target interpreter, `performance_probe.py` calls the public
`poisson_cdf_Q_log10_prepare` at 500, 65,536 and 4,194,302 counts, then repeats
the maximum and returns to 500. Weak references to NumPy allocations made
during preparation measure live table bytes without accessing private Cython
symbols. Only preparation is instrumented. Each phase also executes three
vector/scalar scoring jobs and checks exact float bit patterns. Record growth
time, array count, retained bytes, allocation count and process high-water RSS.
The final two float64 tables occupy exactly 64 MiB, remain available to later
jobs, and must not allocate again on the warm/smaller requests. Retired arrays
must be gone once no reader retains them. Temporary construction/overlapping
readers can exceed 64 MiB; this is a current-generation budget, not a process
memory ceiling. Existing concurrency regressions cover active-reader overlap.

The original `170556d` uses native `realloc` buffers, so its table payload is
not visible to this NumPy observer. Its payload/count fields are explicitly
null and only RSS/growth timing and repeated scoring equivalence are measured.
`compare` records this observation limit for that exact historical revision;
the 64 MiB payload/lifetime gate applies to extraction and current builds.

Current builds also create and destroy three independent public serial peak
callers in that warm-cache interpreter, checking output hashes and retained
table bytes after each cleanup. Historical builds lack the current explicit
one-CPU execution policy, so their cache probe reports repeated probability
jobs only. The extraction already contains peak-process safety controls;
the probe does not synthesize an all-phase CPU budget for it. Releasing
a caller does not clear the module's cache. A fresh CLI process releases it on
exit. CPU pools, allocator retention and RSS need separate interpretation from
the measured NumPy table payload.

## Profiling and vendor ownership

Normal wheels compile with `-O3 -g0` in `setup.py`. Python-only profiling cannot
attribute native/nogil time. For native profiling, make a separate development
build with debug symbols (remove the trailing `-g0`, use `-g`), retain its exact
source diff/flags and Cython-generated C, and use platform native profiling
tools. Add Cython line tracing only in a separate instrumented build when
needed; its timing is not comparable to a release wheel. Match optimization,
CPU architecture and dependencies for release timing. Do not merge diagnostic
flags as a performance-default change.

MACS maintainers own the integration and provenance of bundled code:

| Component | Provenance owner and update boundary |
| --- | --- |
| `MACS3/IO/libdeflate` | IO maintainers; vendored libdeflate v1.26 as recorded by `setup.py`, license `COPYING`. Record upstream tag/commit and local integration diff for any update. |
| `MACS3/fermi-lite` | Variant-calling maintainers; retain upstream/local-change history and `LICENSE.txt`. |
| `MACS3/fermi-lite/lib` | SIMDe submodule pinned to `8cd136a43bae7ab9b82316179b9cef8887726778` for these revisions; retain its `COPYING` and exact submodule SHA. |
| NumPy-derived allocator code | Signal maintainers; retain `MACS3/Signal/NUMPY_LICENSE.txt` and provenance comments when changing the allocator. |

Keep upstream algorithm constants in their upstream owner. Benchmark evidence
does not authorize editing vendor tuning constants. Vendor updates require a
separate review of provenance, licenses, architectures, scientific equivalence
and native distribution checks. The repository cannot currently substantiate
an upstream commit for every historical vendored file; document that gap rather
than inventing a provenance claim.

## Release checklist

Run `uv run pytest test` and the focused native tests with `--run-slow`, build
wheels, and run the new subprocess tests against both an in-place build and an
installed wheel from a clean checkout without extensions. The native regression
runner and compatibility workflow include these tests with psutil installed.
Preserve Linux/macOS results and compare all revisions before considering the
evidence complete. Uncollected platforms/profiles remain open in the report.
Any proposed retuning must be a separate commit/PR with workload evidence,
tradeoffs, numerical/byte equality and release validation; a benchmark report
alone establishes no optimal chunk size, cache reserve or startup threshold.
