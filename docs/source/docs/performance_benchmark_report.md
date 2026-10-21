# Performance-defaults benchmark report (#773)

Reviewed on October 20, 2026. This report combines the initial local measurements
with a completed Linux and macOS workflow. No performance defaults or vendor
constants were changed.

## Findings in plain language

**Recommendation: close #773 and its parent tracker #758 once this report update
is merged. Keep the current performance defaults.** The requested benchmarks,
output checks and documentation are complete. Closing these tickets means the
planned follow-up work is complete; it does not mean every setting is optimal
or that the current implementation is always faster.

All **972 benchmark invocations** in the hosted run succeeded. The peak results
and signal files matched across all tested revisions, CPU budgets and backends
within each platform. The cache checks also passed: the current logarithm
tables stay within their 64 MiB budget and are reused by later jobs.
The file comparison excludes command/version comments in XLS headers; the
scientific data must still match exactly.

Linux timings were close to the original merged PR #756 baseline. Mac timings
were more variable, and some current medians were slower. For example, the
large FRAG case on the Intel Mac rose from 1.953 to 2.520 seconds, while the
large BAM case on the ARM Mac rose from 0.761 to 0.971 seconds. These slower
results remain part of the report. Three repetitions on hosted runners are
not enough to separate all timing variation from the effects of safety fixes
and execution-policy changes.

The results do not show a consistent benefit from raising the CPU budget or
identify a better set of defaults. Changing buffer sizes, cache budgets or
worker thresholds would need a separate proposal and measurements supporting
that particular change.

## Completed cross-platform run

[GitHub Actions run 39374432219](https://github.com/macs3-project/MACS/actions/runs/39374432219)
tested commit `85f6d08fe97e5249d586c86d9e0bbc179e2774f2`, including the workflow
fix subsequently merged in PR #795. All three runners built and installed the
original baseline (`170556d`), the default-preserving extraction (`324f779`),
and the current implementation. **The baseline is merged PR #756, not the
`v3.0.5` release tag.**

Each profile used three repetitions of nine generated cases: small and large
BAM, BAMPE and FRAG, plus many-contig versions of each format. Current runs
requested budgets of 1, 2, 4 and 8, with serial and automatic backends. A requested
budget is an upper limit; available CPUs and the safety policy can reduce the
actual worker count. The workflow alternated revision order between rounds.

| Platform and profile | CLI invocations | Output/cache checks | Runtime change from baseline | Runtime change from extraction |
| --- | ---: | --- | --- | --- |
| Linux, host | 243 | Passed | 0.3% faster to 0.7% slower | 0.5% to 2.8% slower |
| Linux, constrained | 243 | Passed | 0.7% faster to 5.7% slower | 0.8% faster to 5.1% slower |
| macOS Intel | 243 | Passed | 14.1% faster to 29.0% slower | 36.3% faster to 23.4% slower |
| macOS ARM | 243 | Passed | 3.3% faster to 27.6% slower | 16.7% faster to 7.0% slower |

The runtime ranges above compare medians for the six large/many-contig cases,
using the **requested 8/auto** policy for the current build and legacy defaults
for both historical revisions. They are ranges across workloads, not confidence
intervals or matched-policy speedups. The Intel Mac's baseline many-BAMPE runs,
for example, ranged from 3.138 to 9.791 seconds. This variation limits how much
we can infer from a single percentage.

The [Linux summary](https://github.com/macs3-project/MACS/actions/runs/39374432219#summary-118469270581),
[Intel Mac summary](https://github.com/macs3-project/MACS/actions/runs/39374432219#summary-118469270500)
and [ARM Mac summary](https://github.com/macs3-project/MACS/actions/runs/39374432219#summary-118469270316)
contain the full median/minimum/maximum tables, including the extraction and
every tested CPU policy. All three report no output mismatches or failed runs.
The workflow took 35 minutes 38 seconds overall, including builds and artifact
uploads. Each summary group contains three measured repetitions.

Linux's constrained profile limited CPU affinity to at most two CPUs and set
a 2 GiB address-space limit for each process. This checks restricted execution;
it does not test a shared 2 GiB physical-memory limit for the whole process tree.
The environments used Python 3.14.8 on Linux and 3.14.7 on the Macs, with NumPy
2.3.5. Each artifact records the resolved dependency versions, build logs,
revision identifiers, input hashes, per-run measurements and outputs.

The extraction and current builds passed the cache checks on every platform:
two arrays totaling 64 MiB remained after preparation, warm requests did not
allocate new tables, and repeated current callers retained the same cache
budget after cleanup. The original baseline's native cache allocation is not
visible to this observer; its limitation is reported rather than counted as a
measured 64 MiB payload.

Memory and temporary-storage measurements need care. Summed process RSS counts
shared pages more than once, so it cannot prove unique physical-memory savings.
Mac temp measurements omit deleted files that are still open. All successful
runs finished without retained temp files or sampling errors, but a zero sample
does not prove that no temporary storage was used.

The [protocol](performance_benchmark.md) records reproducible commands,
measurement definitions, profiling-build instructions and vendor ownership.
The full hosted artifacts are available from the run page; they expire on
November 19, 2026. The figures above preserve the main findings in this report.

## Initial local measurements

The following detailed tables preserve the original October 19 local run.
They describe one Intel Mac and should not be read as the hosted-run results.
This is a synthetic-workload report; it does not establish optimal buffer
sizes, task thresholds, cache reserves or prefetch distances.

The [protocol](performance_benchmark.md) specifies the revisions, generated
datasets, commands, measurement semantics and gates. Local validation completed
the full suite (4,899 passed, 365 skipped; 693.60 s), built all three wheels,
and checked the new subprocess launchers against both the in-place build and
an installed wheel from a clean checkout without in-place extensions.
The installed-wheel run including related native CLI/probability regressions
passed 24 tests with one Linux-only skip. After the cache-observer correction,
the focused launcher tests passed in both modes (3 passed, 1 Linux-only skip).

The local host is an Intel Core i5-8259U at 2.30 GHz (4 physical/8 logical CPUs),
16 GiB RAM, macOS 15.7.7 x86_64. The matched environments use Python 3.14.6,
NumPy 2.3.5, SciPy 1.16.3, Cython 3.3.0 and setuptools 84.0.0. Wheels were
built with Apple clang 17.0.0 and the normal `-O3 -g0` extension flags.
All tested tracks reported NumPy's `default_allocator` on this host.

### Local run coverage

All **243 CLI invocations succeeded**. Each of the nine cases has 27 runs:
three baseline defaults, three extraction defaults, and three repetitions of
seven current policies (1/serial plus 2, 4 and 8 with serial and auto).
All five scientific artifacts matched byte for byte within each case after
excluding XLS comment metadata. Cache payload/reuse/lifetime gates passed for
the extraction and current builds; the original native-cache payload remains
unobserved by the NumPy tracker.

The complete driver took **598.18 s (9.97 minutes)**, including metadata/cache
probes and checksum verification, excluding builds and tests. Each outer round
used `run --repeats 1 --budgets 1,2,4,8` with the same generated manifest and
an explicitly selected installed wheel package root, CLI and Python interpreter.
Revision order was baseline/extraction/current, current/extraction/baseline,
then baseline/extraction/current. There were no simultaneous builds or tests.

The 585 KiB {download}`audit data <benchmarks/773-macos.json>` contains every
run, exact CLI arguments, input/extension/wheel hashes, environment versions,
cache phases and all policy variants. Local path prefixes are normalized to
`<checkout>` and `<artifact-root>`; the scientific hashes are unchanged.

### Local wall time

Median [minimum, maximum] seconds, three repetitions per cell. These are
within-host implementation comparisons with different recorded historical
policies; they are not matched serial-budget speedups.

| Case | Baseline legacy | Extraction legacy | Current 1/serial | Current 8/auto |
| --- | ---: | ---: | ---: | ---: |
| small-bam | 0.525 [0.511, 0.568] | 0.521 [0.511, 1.445] | 0.530 [0.512, 0.536] | 0.492 [0.490, 0.511] |
| small-bampe | 0.538 [0.536, 0.545] | 0.485 [0.478, 0.655] | 0.500 [0.492, 0.540] | 0.528 [0.494, 0.631] |
| small-frag | 0.546 [0.532, 0.550] | 0.485 [0.481, 0.490] | 0.499 [0.494, 0.522] | 0.498 [0.492, 0.499] |
| large-bam | 1.769 [1.762, 1.795] | 1.875 [1.800, 1.897] | 1.913 [1.887, 1.962] | 1.867 [1.806, 1.914] |
| large-bampe | 2.151 [2.146, 2.156] | 2.239 [2.185, 2.268] | 2.319 [2.223, 2.327] | 2.269 [2.224, 2.377] |
| large-frag | 2.311 [2.224, 2.385] | 2.443 [2.408, 2.520] | 2.654 [2.554, 2.852] | 2.452 [2.418, 2.499] |
| many-bam | 3.082 [3.079, 3.113] | 3.322 [3.161, 3.410] | 3.484 [3.429, 3.537] | 3.311 [3.206, 3.362] |
| many-bampe | 3.636 [3.591, 3.645] | 3.811 [3.807, 3.873] | 4.011 [3.915, 4.248] | 3.833 [3.749, 4.363] |
| many-frag | 3.796 [3.781, 3.872] | 4.040 [3.956, 4.501] | 4.170 [4.112, 4.249] | 4.021 [4.004, 4.247] |

Across all policies, the longest invocation took 4.625 s. Most full-size runs
take about 2–4 s on this host; small cases take about half a second.
Startup/cache effects remain visible rather than being discarded as warm-up.

For the six large/many-contig cases, current 8/auto medians were 5.4–7.4% slower
than the original baseline, and between 0.5% faster and 1.3% slower than the
extraction. The comparisons include intervening safety fixes and execution
policy changes; they do not isolate the cost of a specific fallback or
extraction. The table does not justify retuning a default.

### Local selected CPU budgets

Current medians in seconds. Automatic peak processes fall back to serial
on this Mac; native parser/track work can still use the selected budget.

| Case | 1/serial | 2/auto | 4/auto | 8/auto |
| --- | ---: | ---: | ---: | ---: |
| small-bam | 0.530 | 0.492 | 0.497 | 0.492 |
| small-bampe | 0.500 | 0.494 | 0.493 | 0.528 |
| small-frag | 0.499 | 0.523 | 0.502 | 0.498 |
| large-bam | 1.913 | 1.901 | 1.798 | 1.867 |
| large-bampe | 2.319 | 2.280 | 2.184 | 2.269 |
| large-frag | 2.654 | 2.462 | 2.492 | 2.452 |
| many-bam | 3.484 | 3.316 | 3.315 | 3.311 |
| many-bampe | 4.011 | 3.903 | 3.750 | 3.833 |
| many-frag | 4.170 | 4.458 | 3.995 | 4.021 |

These three-repetition synthetic results do not support a universal optimal
CPU budget. Explicit serial-backend runs at 2/4/8 are also preserved in the
audit data and passed the same output gate.

### Local RSS and temporary storage

Medians of sampled per-run peaks in MiB. Tree RSS double-counts shared pages
across forked processes; lower values cannot be read as unique physical-memory
savings. Temp values exclude deleted descriptors on macOS and can miss brief
peaks between 20 ms samples. **Zero means no bytes observed, not no temp use.**
All runs finished with zero retained temp bytes and no sampling errors.

| Case | Baseline tree RSS | Current 1 RSS | Current 8/auto RSS | Baseline temp | Current 8/auto temp |
| --- | ---: | ---: | ---: | ---: | ---: |
| small-bam | 46.36 | 36.71 | 46.97 | 0.52 | 0.52 |
| small-bampe | 47.93 | 40.00 | 47.29 | 0.00 | 0.79 |
| small-frag | 40.09 | 39.05 | 39.76 | 0.00 | 0.62 |
| large-bam | 141.53 | 63.58 | 87.34 | 27.17 | 27.17 |
| large-bampe | 189.47 | 78.71 | 100.63 | 36.63 | 36.63 |
| large-frag | 204.07 | 90.75 | 93.90 | 33.63 | 33.63 |
| many-bam | 128.34 | 53.33 | 73.55 | 43.71 | 43.71 |
| many-bampe | 138.80 | 56.62 | 93.79 | 66.24 | 66.24 |
| many-frag | 158.36 | 71.44 | 90.05 | 56.19 | 56.19 |

### Local cache growth and lifetime

The observer sees 16,384 bytes at the minimum prepared prefix, 1,048,608
bytes after preparation at 65,536, and 67,108,864 bytes at the cap in both
NumPy-backed implementations. The final generation has exactly two arrays.
Warm cap and smaller requests retain the same generation without new table
allocations. Three independent current callers per interpreter produce equal
outputs and retain the same 64 MiB after cleanup.

| Revision | Median preparation at cap (ms) | Median post-cap process high-water RSS (MiB) | Live table payload |
| --- | ---: | ---: | --- |
| baseline | 95.289 | 97.90 | unobserved native realloc |
| extraction | 81.842 | 99.46 | 64 MiB / 2 arrays |
| current | 80.160 | 99.48 | 64 MiB / 2 arrays |

Preparation timings include the temporary allocation observer. Warm timings
in the JSON include Python instrumentation overhead and should not be treated
as native API latency. OS high-water RSS includes imports, allocator retention
and temporary construction; it is not the table payload alone.

## Evidence status

| Coverage | Status |
| --- | --- |
| macOS x86_64, generated small/large BAM/BAMPE/FRAG and 256 contigs | Measured; all output gates passed |
| Current requested budgets 1, 2, 4, 8; explicit serial and automatic policies | Measured; all output gates passed |
| Current 64 MiB log-cache payload, warm reuse and repeated caller cleanup | Public-API probe implemented and validated |
| Linux host and constrained affinity/address space | Measured in the hosted run; all output/cache checks passed |
| Shared cgroup/VM memory ceiling | Not measured; a separate physical-memory limit would inform future tuning |
| macOS arm64 | Measured in the hosted run; all output/cache checks passed |
| Complete unlinked-file temp accounting on macOS | Unavailable to this sampler; reported as a partial lower bound |
| Independent biological large-dataset replication | Not measured; useful for future tuning on real datasets |

The requested cross-platform profiles have now been measured and reviewed.
The remaining measurement limits stay documented; they do not establish a
reason to change a default or add a new requirement for closing this benchmark
ticket. Future tuning must retain the output-equality checks.

## Historical diagnostic differences

Both historical revisions expose an unrelated `PETrackII.destroy()` bug:
the method calls `resize` on its barcode dictionary. The allocator-only probe
records that exact exception and exits its isolated process. CLI timing and
output checks are not bypassed. The historical source is unchanged.

The original `170556d` log cache uses native `realloc`, so NumPy allocation
weak references cannot measure its payload. Its table-payload fields are null;
RSS growth, preparation time and repeated scoring outputs remain measured.
The extraction and current cache payload/lifetime gates use actual retained
NumPy arrays. Historical builds lack the current explicit one-CPU execution
policy, so their repeated jobs exercise probability APIs instead. The
extraction already includes peak-process safety controls, but the probe does
not synthesize an all-phase CPU budget for it.

## Closing #773 and the parent tracker #758

**Both tickets are ready to close once this report update is merged.** #773
asked for reproducible benchmark evidence before retuning. It did not require
a speedup or a new set of defaults. The completed work records the comparisons,
the Mac safety fallback, CPU policies, cache lifetime, resource measurements
and their limits. No tuning change is proposed.

The parent [tracker #758](https://github.com/macs3-project/MACS/issues/758)
contains 15 ordered implementation tickets. As checked on October 20, the other
14 tickets (#759 through #772) are all closed as completed. The
[native regression and packaging run](https://github.com/macs3-project/MACS/actions/runs/39364431609)
passed all nine jobs, covering Linux x64/ARM and macOS Intel/ARM, in-place
compiled extensions, source and wheel packaging, and installed-wheel tests
from clean checkouts. The defaults, backend choices and benchmark results are
documented. Completing #773 therefore completes the remaining task in #758.

Closing the tracker records completion of that plan. It does not certify all
possible datasets, memory limits or future dependency versions. Those broader
questions can be addressed by specific future changes with their own evidence.

## Guidance for future tuning

macOS serial peak-process execution is an intentional safety fallback from
#761. Native parser/track threads can still use the selected budget. The
historical legacy-default measurements are not matched serial-budget controls.
Report them as implementation comparisons with different recorded policies.

Keep the current values until Linux process startup/batching, real memory
constraints, allocator residency and representative biological inputs provide
enough evidence for a specific change. In particular, the 262,144-location
process threshold is still a conservative heuristic, not a measured universal
break-even point. Any retuning belongs in a separate commit or PR with the
output-equivalence gate and release validation repeated.
