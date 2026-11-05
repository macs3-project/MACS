# BEDPE/FRAG coordinate validation verification

Issue: [#815](https://github.com/macs3-project/MACS/issues/815),
part of [#757](https://github.com/macs3-project/MACS/issues/757), item 18.
Branch: `fix/macs/815-nonnumeric-fragment-coordinates`.
Baseline: current upstream `main` at
`b7929196a98032219afbef43b00ce3e0996aacb7`.

## Compiled baseline

Rebuilt the unchanged source using the `macs3dev` Python 3.14 interpreter
with `python setup.py build_ext --inplace` (Python 3.14.6, Cython 3.3.0,
NumPy 2.5.1, pytest 9.2.0). The imported parser was the
in-place `MACS3/IO/Parser.cpython-314-darwin.so`.
`uv run` could not run in the sandbox (cache permissions, then a macOS
system-configuration panic with a temporary cache), so verification uses
Python directly from the required conda environment.

Ran 192 public-API probes, each in a subprocess with a 20-second timeout,
the parent's interpreter and MACS3 installation, and a temporary working
directory. The matrix covers BEDPE and FRAG, build and append, normal
Python and `-O`, and records with and without a final newline.

| Input | Baseline result in both formats |
| --- | --- |
| `chrom start end barcode count` (tab-separated) | Misleading interval-order `AssertionError`; under `-O`, accepted as a zero-length fragment |
| Start `start` or an empty start, end `20` | Accepted with start zero, including without `-O` |
| Start `10`, end `end` | Interval-order `AssertionError`; under `-O`, accepted with negative fragment length |
| Start `10junk` or end `20junk` | Numeric prefix accepted |
| Equal/reversed coordinates | Interval-order `AssertionError`; under `-O`, accepted |
| Empty/comment-only file | `ValueError: No valid fragments found in BEDPE/FRAG input` |
| Leading `#` header followed by a valid row | Correctly skipped; valid fragment retained |
| Valid row | Correct coordinates, count, and mean fragment length |

These reproduce both the inaccurate diagnostic and optimization-dependent
validation on current `main`, rather than relying on the original report's
older revision. Temporary baseline evidence is in `/tmp/macs815_verify.py`,
`/tmp/macs815-results.json`, `/tmp/macs815-summary.txt`, and
`/tmp/macs-815-build.log`.

## Changes and coverage

BEDPE and FRAG now check complete ASCII decimal coordinate fields before
converting them to signed 32-bit track positions. Existing optional signs,
surrounding ASCII whitespace, and leading zeros remain accepted. Negative
starts and empty chromosomes are still skipped when their coordinates are
valid. Invalid or out-of-range coordinates raise `ValueError` identifying
the start/end field, its value, and the offending row. Interval-order checks
also raise explicit `ValueError`, including under `python -O`.

FRAG's usual digit-only C paths remain in place. Its general C walk now
accepts only complete digit fields; signs, longer fields, and malformed
coordinates fall back to the checked line parser. Fragment count parsing
and capping are unchanged.

`test_Parser_fragment_coordinates.py` adds 40 public-API subprocess tests.
All validation assertions run in the parent, so optimized child Python
cannot remove them. Each child uses `sys.executable`, `macs3_env`, a
`tmp_path` script and working directory, and a 60-second timeout. The
parent also checks that the child imported the exact same compiled parser.

Coverage includes both formats and build/append operations, normal and
optimized Python, plain/gzip/BGZF input, first/later rows, complete lines
and EOF tails, empty and partial coordinates, column headers, decimal and
exponent notation, signs without digits, underscores, NUL bytes, int32
overflow/underflow, very long fields, and equal/reversed intervals. Valid
controls verify coordinates, counts, parser statistics, barcode assignment,
existing append records, signed/padded coordinates, int32 limits, skipped
records, empty/comment-only input, and incomplete rows. A >64 KiB BGZF
control exercises serial and two-thread FRAG parsing, with malformed
records in the middle and at EOF. Existing parser/EOF/singleton regressions
expect the new explicit exception type for invalid intervals. Two callpeak
regression modules now require the accurate non-numeric-coordinate
diagnostic instead of the old interval-order message.

## Validation

- In-place build: passed.
- Initial parser/EOF run: **314 passed, 99 skipped** in 79.68 seconds.
- Final in-place coordinate and singleton regressions: **46 passed** in
  32.32 seconds.
- Installed-wheel run of every parser module plus singleton FRAG regressions:
  **444 passed, 99 skipped** in 147.65 seconds.
- Both affected callpeak CLI modules: **75 passed** in-place (62.30 seconds)
  and **75 passed** against the installed wheel (59.47 seconds).
- Final parser source and all changed tests match the clean source snapshot.
- `git diff --check`: passed.
- Full `python -m pytest -q test` run: **5,476 passed, 491 skipped,
  4 failed** in 1,073.85 seconds. All four failures were stale diagnostic
  or exception expectations in tests already imported before their updates:
  two callpeak `--llocal 0` malformed-input cases, the callpeak local-window
  malformed-input case, and the singleton FRAG invalid-interval case.
  The final versions of their complete modules pass in-place and against
  the installed wheel, as recorded above. The full run reported no other
  failures.
- `python -m pytest -q test --lf`: **4 passed, 434 deselected** in 6.96
  seconds (pytest limits collection to the previously failing modules).
  All four full-suite failures are resolved; the entire suite was not
  rerun a second time.
- Existing skips were retained; no new skip was added.

The wheel is built from a temporary clean source snapshot containing the
tracked baseline and these changes, with the bundled SIMDe headers supplied
from the workspace. Installation goes into a separate temporary target.
The wheel runner starts with `python -P`, imports MACS3 from that target,
checks the compiled parser path and source equality, and verifies that the
clean source's `MACS3/` tree contains no in-place extensions before invoking
pytest on the clean source tests. This exercises pytest's subprocess
launcher against the installed wheel, rather than probing an extracted
wheel alone.

Final evidence: `/tmp/macs815-fix-build.log`,
`/tmp/macs815-targeted.log`, `/tmp/macs815-final-inplace.log`,
`/tmp/macs815-singleton.log`, `/tmp/macs815-cli-tests.log`,
`/tmp/macs815-full-suite.log`, `/tmp/macs815-full-suite-rerun.log`, and
`/tmp/macs815-release/` (`build-wheel.log`, `install-wheel.log`,
`wheel-tests.log`, `wheel-cli-tests.log`, the clean source snapshot,
installed package, and built wheel).
