# Disjoint chromosome-name verification

Tickets: [#802](https://github.com/macs3-project/MACS/issues/802),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

## Baseline

Verified compiled current-main revision
`fb28bda92295511908e8de40b36fd3e268b91ae1` on macOS with Python 3.14.6,
Cython 3.3.0, NumPy 2.5.1, and pysam 0.24.1. GitHub's main-branch API
confirmed this revision before verification. Built in place using
`python setup.py build_ext --inplace -j 4` with the pinned SIMDe submodule.

The original seed-0 reproduction, with 2,000 BED reads each on `chr1`
and `chr2`, exits 1 with `TypeError: sequence item 0: expected str
instance, bytes found`, before listing either chromosome set. An otherwise
identical run using the treatment as its own control succeeds.

Public calls to both commands' `check_names()` using compiled `FWTrack`
objects reproduce the TypeError for disjoint names, one empty track, and
invalid UTF-8 chromosome bytes. Both empty tracks instead raise
`SystemExit(None)`. The actual CLI command with two empty BED files prints
the no-common-chromosome error but exits **0**. A BED containing a
non-numeric start position also reaches the failing name diagnostic when
its chromosome differs from the control's.

The `callvar` helper is unused by its CLI. Disjoint, indexed BAM headers
instead produce an `AssertionError` from the first-chromosome membership
check, without listing the names. That assertion can also reject headers
which share a later chromosome, and disappears under Python optimization.

Empty BEDPE and FRAG inputs were checked as well: they fail with
`ZeroDivisionError` in their loaders before chromosome validation. An empty
`PETrackI.finalize()` also divides by zero. These existing empty-input
defects belong to [#808](https://github.com/macs3-project/MACS/issues/808)
and are unchanged here. Tests of empty paired-end track names use newly
constructed tracks without finalization.

## Change and coverage

Both commands use a shared chromosome-set check. On disjoint sets it
prints the assembly guidance and sorted treatment/control name lists,
then raises `SystemExit(1)`. Bytes are decoded as UTF-8 with invalid bytes
escaped, so malformed names cannot replace the useful diagnostic with a
decoding traceback. String-valued callers remain supported. Matching names
retain the existing silent return behavior; comparisons still use the
original names, not their display strings.

`callvar` applies this check to all BAM header names before opening the
output VCF, replacing the first-name assertion. Indexed BAMs with no
reference names are handled without indexing an empty list.

`test_chromosome_names.py` contains 50 regressions covering compiled
single-end and paired-end tracks, multiple sorted names, empty sets,
Unicode and invalid UTF-8, BED/BEDPE/FRAG CLI diagnostics, the original
reproduction and successful control, and indexed BAM header checks.
Header tests include empty sets, overlap only beyond the first name, and
execution with Python assertions disabled. Existing command tests are
updated for the failure status and consistent diagnostic.

Every CLI subprocess has a 30-second timeout and uses `run_macs3`, which
launches the parent's interpreter with `macs3_env` pointing at the parent's
MACS3 installation. Compiled-track tests assert that native extensions are
loaded. The wheel validation additionally checks that both the extension
and new utility import from the isolated wheel installation.

Repeating the baseline probes after the fix confirms exit 1 for all
disjoint/empty-set checks and the original, empty BED, malformed BED, and
disjoint-header CLI commands. Overlapping-track calls and the valid BED
CLI control succeed. Failure diagnostics contain the chromosome lists
without a traceback.

## Validation

- Compiled in-place build succeeded.
- New regression module: **50 passed** in 19.00 seconds.
- Broader in-place run of `test_chromosome_names.py`,
  `test_callpeak_cmd.py`, and `test_callvar_cmd.py`: **383 passed,
  1 failed, 1 skipped** in 202.59 seconds. This run collected the earlier
  41-case version of the new module; the subsequent standalone run above
  covers its final 50 cases.
- The failure is the existing
  `test_callvar_fermi_overlap_longer_than_reads_aborts`: it expects
  `SIGABRT` but this build raises the empty-unitig `TypeError` instead.
  Repeated that exact test on an unmodified archive of baseline `main`,
  using the same compiled extensions: it fails identically (5.75 seconds).
  This is the fermi/empty-unitig behavior tracked in
  [#813](https://github.com/macs3-project/MACS/issues/813), not a name-check
  regression.
- Clean-source wheel build succeeded with
  `python -m build --wheel --no-isolation`, including the pinned SIMDe
  submodule. Installed using `pip --no-deps --no-index --target` into an
  isolated temporary directory. The clean source tree has no in-place
  extensions; validation uses `python -P -m pytest` with the installed
  directory on `PYTHONPATH`.
- Installed-wheel run of the same three modules, including all 50 final
  regressions: **392 passed, 1 failed, 1 skipped** in 180.18 seconds.
  The only failure is the same baseline-confirmed fermi-overlap test;
  all chromosome-name tests pass. The skip is the suite's opt-in slow test.
- Verified that the wheel's two command modules and new utility exactly
  match the working-tree sources. `git diff --check` passed.

The test interpreter is `/private/tmp/macs800-verification/test-env/bin/python`,
which inherits the `macs3dev` scientific dependencies. Pytest plugin
autoload is disabled. Temporary probe inputs and logs are under
`/private/tmp/macs802-verification`. The full repository suite was not run.
