# CPU budgets and process backends

Use `--threads N` to set one CPU budget for a MACS3 invocation:

```sh
macs3 callpeak -t treatment.bam -c control.bam -g hs --threads 4
macs3 callpeak -t treatment.bam -g hs --threads 1
```

The budget includes the calling thread. A budget of **one** runs all MACS3
parsing, track processing, model building, scoring, peak calling and cleanup
serially, with no MACS3 worker threads or child processes. The option also
applies to `predictd`, `pileup`, `filterdup`, `randsample`, `refinepeak`,
`hmmratac` and `callvar`. Serial bedGraph commands do not need this option.

The default request is **eight**, clamped to detected available CPUs. Eight
is not a hard maximum: a larger explicit request can use more CPUs if they
are available. Zero, negative, noninteger and excessively large values are
errors. Setting a larger budget does not guarantee a speedup; small phases
can run serially or use fewer participants than the budget allows.

## Configuration precedence

An explicit API execution policy takes precedence; conflicting explicit
overrides are rejected. Otherwise the order is:

1. Explicit `--threads N` or API `threads=N`.
2. The `MACS3_NUM_THREADS` environment variable, a positive decimal integer.
3. The default request of eight.

An explicit count ignores an unused environment value, even if that value is
invalid. MACS3 snapshots the chosen settings once per invocation and passes
them through its phases. Changing the environment later does not change
existing objects. Independent jobs have independent budgets; there is no
global scheduler dividing a single budget between jobs.

`callvar` retains its historical default of **one** when neither an explicit
count nor `MACS3_NUM_THREADS` is provided. Its existing `-m` /
`--multiple-processing` option is a compatibility spelling of the budget;
different explicit values for it and `--threads` conflict. Nonpositive legacy
values now raise an error rather than being silently converted to one.

## Available resources

Detection combines positive host CPU counts, process CPU counts when the
Python/OS supports them, and process affinity. Unsupported queries and
affinity errors fall back to the other signals. Unknown counts and empty
affinity never produce a zero budget; the final fallback is one.

On Linux, MACS3 also reads cgroup v2 `cpu.max` and cgroup v1 CPU quota/period
files, using process membership and mount metadata. It considers the most
restrictive quota among visible ancestors. A quota of 250,000 microseconds
per 100,000-microsecond period permits two whole CPU participants; a positive
quota below one CPU permits one. Unlimited, missing, inaccessible or malformed
quota metadata does not abort the job. Other detection signals remain usable.
Debug logging records the detected limits and unavailable metadata.

This is a snapshot. Hidden cgroup ancestors, scheduler allocations not exposed
through OS limits and later changes to resource assignments cannot be inferred
reliably. Supply a smaller explicit budget when necessary.

## Threads, processes and cleanup

For a budget `B`, threaded phases use the caller and at most `B-1` helpers.
BGZF decompression and FRAG line parsing share the same native thread pool.
Peak process phases use at most `B-1` children, reserving one slot for the
parent coordinator. Thus the default total budget eight permits at most
seven peak workers; budget two permits one child. Workers do not start nested
MACS3 pools. Parent p-score threads start after process workers have joined
and can use the full thread budget.

Cleanup also obeys the budget. Unlinked pileup files are closed on the caller
at budget one; additional cleanup threads are joined before `destroy()`
returns. Cleanup no longer continues in the background during a later phase.

This supplies the configurable chromosome-process support requested in
[issue #665](https://github.com/macs3-project/MACS/issues/665).

`callpeak` and `callvar` also accept `--process-backend {auto,serial}`:

- `auto` permits fork only on Linux, from the main Python thread, outside a
  daemon process, with no other detected Python threads and with fork support.
  Eligibility is checked again immediately before pool creation. Other hosts
  use serial process execution.
- `serial` disables process pools. Other phases can still use budgeted threads.
  Use `--threads 1` for completely serial MACS3 execution.

Native library threads and concurrent thread creation cannot be detected
reliably. Embedded hosts that cannot establish fork safety should select
`serial`. Logs report the requested/effective budget, its source, process
selection and safety fallbacks. Logged maxima describe allowed capacity;
actual phase worker counts can be smaller.

## Python API

```python
from MACS3.Utilities.Execution import available_cpus, resolve_execution_policy
from MACS3.IO.Parser import BAMParser
from MACS3.Signal.CallPeakUnit import CallerFromAlignments

policy = resolve_execution_policy(threads=4, process_backend="serial")
parser = BAMParser("treatment.bam", execution_policy=policy)
track = parser.build_fwtrack()
track.finalize()
caller = CallerFromAlignments(track, None, execution_policy=policy)
try:
    peaks = caller.call_peaks(["q"], [2.0])
finally:
    caller.destroy()
```

Parser, FWTrack, PETrackI/II, PeakModel and CallerFromAlignments constructors
accept optional `execution_policy=` and `threads=` keywords. Each exposes a
read-only `execution_policy`; the policy itself is immutable. With no supplied
policy, each constructor resolves its own invocation settings. `guess_parser`
accepts the same settings and resolves once before probing candidate formats.
High-level commands resolve once for all their objects.

Returned tracks and their samples, conversions and barcode subsets inherit
the original policy. Appending uses the parser's budget for parsing and retains
the supplied track's budget for later track processing. A caller's policy
governs scoring/calling independently of its already-built input tracks.
Pass the same policy explicitly to objects belonging to one invocation.

Legacy `Threads.MAX_THREADS` and `CallPeakUnit.MAX_WORKERS` names remain
informational aliases of the internal default. Rebinding them does not
configure jobs. Use the resolver and constructors instead.

MACS3's budget controls MACS3 pools. BLAS, OpenMP, scikit-learn and HMM libraries
have their own resource settings, such as `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS`. MACS3 does not rewrite those
environment variables or promise to control third-party pool sizes.

## Chromosome work and process startup

Peak calling uses stored treatment locations plus stored control locations
when the control is a distinct track. A paired-end fragment is one location;
a FRAG record is one location regardless of its count. Only chromosomes common
to treatment and control participate, as before. Chromosomes empty in both
inputs contribute no work, histogram segments, peaks or bedGraph rows. An
empty comparison returns an empty peak collection; requested bedGraph files
are empty (apart from requested track lines), and empty cutoff analysis writes
its column header. A chromosome empty in only one input remains in the
comparison: treatment coverage is zero or control uses the background baseline.

Both the p-score histogram phase and narrow/summit process phase require at
least **262,144 stored locations** across at least two nonempty chromosomes.
Smaller workloads run serially even when the invocation permits processes.
This internal minimum is a conservative startup heuristic, not a measured
Linux break-even point or a promise of speedup above the threshold. Local
macOS observations found three-chromosome serial call/cleanup workloads up to
196,608 stored locations took roughly 3–7 ms. Linux startup costs, local-control
windows, data distribution, hardware and filesystem costs vary. Budgets and
host-safety checks remain upper limits; they do not force pool creation.

For `W` eligible workers, the scheduler targets at most `32 * W` tasks. It
uses a positive ceiling of total locations divided by that task limit, sends
large chromosomes individually, and groups smaller chromosomes until they
reach the target. Only the final batch may fall short. Empty chromosomes
never become tasks. Chromosomes are not split, so a dominant chromosome can
still limit parallelism. Task submission is largest first; histogram ties and
returned peaks retain serial chromosome/position order regardless of completion
order. Each phase joins its workers before the next phase begins.
