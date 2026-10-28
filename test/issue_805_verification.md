# Ticket #805: singleton-strand duplicate-filter totals

- Ticket: https://github.com/macs3-project/MACS/issues/805
- Aggregate: https://github.com/macs3-project/MACS/issues/757 (item 8)
- Branch: `fix/macs/805-singleton-strand-totals`
- Baseline: `90b61de74af2ddd836283753d1c18aacee9f4b11`, verified as current
  upstream `main` before implementation.
- Environment: macOS x86_64, macs3dev Python 3.14.6, Cython 3.3.0,
  NumPy 2.5.1.

## Verification before the fix

Built the unchanged revision with `python setup.py build_ext --inplace -j 4`.
Confirmed that `MACS3.Signal.FixWidthTrack` loaded the compiled `.so`.
The original two-line BED reproduction exited 1 with default filtering:
2 input tags, 0 retained tags, and `AssertionError: Lambda must > 0`.
The same command with `--keep-dup all` exited 0.

The four-distinct-read example (three plus-strand tags on chr1 and one
minus-strand tag on chr2) wrote four BED records but reported three retained
and redundancy 0.25. Public `FWTrack` probes confirmed the defect on both
strands, on separate chromosomes, and alongside duplicate groups. Empty
tracks remained empty and disabling filtering preserved totals.

CLI probes ran in a temporary working directory, with the parent's Python
interpreter and MACS3 package root, and a 30-second timeout. Empty `filterdup`
input raised `ValueError: No alignments found in input files`; a malformed
one-column BED raised `IndexError`. Additional regression verification found
that empty `callpeak` input raises `ZeroDivisionError` when calculating
redundancy. These existing invalid-input errors are outside this count fix;
all terminate without a native crash or hang.

Before rebuilding the fix, the new regressions produced **18 failed,
15 passed in 6.35s**. Failures covered singleton API totals, both CLI strands,
automatic duplicate limits, and mixed-input CLI totals.

## Change and coverage

`FWTrack._filter_dup_strand` adds its zero-or-one-element array size to the
per-chromosome retained count before returning. The shared helper handles
both strands and both serial and threaded chromosome execution. Coordinates,
duplicate thresholds, and allocation behavior are unchanged.

`test/test_filter_dup_singletons.py` adds 33 public-API/CLI cases:

- Empty, plus-only, minus-only, both-strand, and multiple-chromosome tracks.
- Singleton strands alongside distinct reads or duplicate groups.
- Limits 1 and 2 and disabled filtering; exact coordinates, return values,
  totals and lengths, including a repeated pass after duplicate collapse.
- Two-read `callpeak` on each strand, with fixed, automatic, or disabled
  duplicate filtering; successful output creation and retained counts.
- Mixed-input `filterdup` with and without an actual duplicate; exact output
  coordinates, retained totals, and redundancy rates.
- Bounded empty/malformed CLI probes for both commands.

CLI tests use `run_macs3` from `test/conftest.py`, which derives `PYTHONPATH`
from the parent's imported MACS3 package and uses `sys.executable`, temporary
working directories, disabled core dumps, and a 30-second timeout.

## Fixed reproduction

The original default-filtering `callpeak` command now exits 0 with two retained
tags and zero redundancy; the `--keep-dup all` control still succeeds.
The mixed-input `filterdup` command reports four retained tags, redundancy
0.00, and writes exactly four records. Empty/malformed errors are unchanged.

Local scripts and logs: `/private/tmp/macs805-verification/`, including
`verify.py`, `baseline-results.log`, `fixed-results.log`, and
`pytest-before-fix.log`.

## Validation

Rebuilt the in-place extensions and ran:

```sh
python -m pytest test/test_filter_dup_singletons.py \
  test/test_FixWidthTrack.py test/test_filterdup_cmd.py \
  test/test_callpeak_cmd.py -q
```

Result: **509 passed in 146.08s**, including all 33 new cases.

Built `macs3-3.1.0-cp314-cp314-macosx_11_0_x86_64.whl` using
`python -m build --wheel --no-isolation` from a fresh `git archive HEAD`
with the changed source and new tests overlaid. Populated the pinned
fermi-lite SIMDe submodule from the existing local copy. Installed the wheel
with `pip install --no-deps --no-index --target` into an isolated temporary
directory without changing the conda environment.

Ran the same four-module pytest selection from the clean source checkout
using `python -P -m pytest`, disabled plugin autoload, and `PYTHONPATH` with
the installed wheel and temporary test dependencies. Confirmed before testing
that MACS3, FixWidthTrack and Parser came from the wheel installation and
that the checkout's `MACS3` directory contained no `.so` files. CLI children
selected the same installed package through `macs3_env()`.

Result: **509 passed in 146.85s**, including all changed subprocess tests.
Neither selection skipped tests. This was targeted validation; the entire
repository suite was not run. `git diff --check` passed.

Logs in the local verification directory: `build-fixed.log`,
`pytest-inplace.log`, `build-wheel.log`, `install-wheel.log`, `pytest-wheel.log`.
Related issues: #805 and #757.

## macOS CI follow-up

The first PR run failed on macOS/Python 3.13 in the pre-existing
`test_native_stream_failure_releases_workers_before_next_load[build-bam]`.
The native cleanup probe observed `stream retained native workers (1 -> 2)`
after a successful load. The run otherwise had 5,075 passing tests, including
all 33 singleton regressions, and 468 skipped tests. The x64, non-x64 and
native compatibility workflows succeeded. Run:
https://github.com/macs3-project/MACS/actions/runs/40770352244

The pool already calls `pthread_join` for every worker before freeing its
storage. The macOS test counter used the raw `task_threads` count, which can
include a joined pthread's Mach thread during teardown. An independent C
probe with 10,000 three-worker create/join cycles reproduced extra raw Mach
threads after completed joins in 1,345 cycles; mapping the same snapshots
through `pthread_from_mach_thread_np` reported zero extra pthreads. A separate
held-worker probe counted all three live waiting workers. The initial
six-case MACS3 stream cleanup selection passed locally before the counter
change, consistent with an intermittent observation race.

The macOS branch of the shared test counter now counts Mach ports that still
map to pthreads using public `pthread_from_mach_thread_np`, with its pointer
return type declared explicitly. Joined remnants are excluded; live pthread
workers, including blocked ones, still count. Snapshot ports and memory are
released in `finally`, including on observation failure. Linux observation,
production parsing/pool code, timeouts, and strict cleanup bounds are unchanged.

Apple's implementation searches the pthread list under its lock and removes
exiting/joined threads from that list:
https://github.com/apple-oss-distributions/libpthread/blob/main/src/pthread.c

Added deterministic macOS-counter regressions for joined remnants versus a
still-live worker, pointer return type, enumeration failure, and snapshot
cleanup on lookup failure. Before the fix these produced 2 failures and
1 passing control. Logs and independent probes are in the existing temporary
verification directory (`thread-counter-before.log`, `mach-join-probe.c`,
`mach-join-results.log`, `mach-map-probe.py`).

Follow-up validation ran `test/test_Execution_paths.py`,
`test/test_native_streams.py`, and `test/test_filter_dup_singletons.py`:

- In-place compiled build: **132 passed, 43 skipped in 92.62s**.
- Installed wheel, from the clean source checkout without in-place extensions:
  **132 passed, 43 skipped in 87.07s**. Only changed test/helper files were
  overlaid; the production wheel is unchanged by this test-only follow-up.
- Both selections include all six build/append stream-failure cleanup cases,
  actual native worker-budget observation, the three new counter regressions,
  and all 33 singleton cases. Slow/host-specific skips use existing markers.
- `git diff --check` passed.

Final follow-up logs: `thread-counter-inplace.log` and
`thread-counter-wheel.log` in `/private/tmp/macs805-verification/`.
