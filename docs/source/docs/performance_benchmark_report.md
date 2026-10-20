# Initial performance-defaults evidence (#773)

Measured on October 19, 2026. No performance defaults or vendor constants were
changed. This is a synthetic-workload report; it does not establish optimal
buffer sizes, task thresholds, cache reserves or prefetch distances.

The [protocol](performance_benchmark.md) specifies the revisions, generated
datasets, commands, measurement semantics and gates. Local validation completed
the full suite (4,899 passed, 365 skipped; 693.60 s), built all three wheels,
and checked the new subprocess launchers against both the in-place build and
an installed wheel from a clean checkout without in-place extensions.
The installed-wheel run including related native CLI/probability regressions
passed 24 tests with one Linux-only skip. After the cache-observer correction,
the focused launcher tests passed in both modes (3 passed, 1 Linux-only skip).

The host is an Intel Core i5-8259U at 2.30 GHz (4 physical/8 logical CPUs),
16 GiB RAM, macOS 15.7.7 x86_64. The matched environments use Python 3.14.6,
NumPy 2.3.5, SciPy 1.16.3, Cython 3.3.0 and setuptools 84.0.0. Wheels were
built with Apple clang 17.0.0 and the normal `-O3 -g0` extension flags.
All tested tracks reported NumPy's `default_allocator` on this host.

## Measured results

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

### Wall time

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

### Selected CPU budgets

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

### RSS and temporary storage

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

### Cache growth and lifetime

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
| Linux host and constrained affinity/address space | Manual workflow added; measurements not collected locally |
| Shared cgroup/VM memory ceiling | Pending; RLIMIT_AS is not a shared RSS cap |
| macOS arm64 | Manual workflow added; measurements not collected locally |
| Complete unlinked-file temp accounting on macOS | Unavailable to this sampler; reported as a partial lower bound |
| Independent biological large-dataset replication | Pending |

Ticket #773's cross-platform evidence remains open until the pending profiles
are measured and reviewed. A configured workflow is not evidence of a passing
run. Native failure/output mismatches must remain visible in the artifacts.

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

## Interpretation and next evidence

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
