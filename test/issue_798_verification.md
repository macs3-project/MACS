# Singleton FRAG chromosome verification

Tickets: [#798](https://github.com/macs3-project/MACS/issues/798),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

Verified against `88c82d62d239d47b9c8ce9e545eeaaa6b7f695de` on macOS
with Python 3.14 in `macs3dev`. The reported ownership defect was already
fixed by `90c31505d55eedc7e05219579b16fc705fd31956`: the sorted-start
branch of `pileup_from_LRC_as_list` copies its input field before resizing
temporary arrays. This change adds regression coverage and corrects an
outdated fast-path comment; it does not change the pileup algorithm.

## Original reproduction

Generated the ticket's 200 chr1 records and varied chr2 between zero,
one, and two overlapping records. Ran both original commands in temporary
directories with the parent's compiled MACS3 installation, with a 60-second
timeout per command.

| Command | No chr2 records | One chr2 record | Two chr2 records |
| --- | --- | --- | --- |
| `callpeak -f FRAG -g 1e6 -t frag.tsv -n x` | exit 0 | exit 0 | exit 0 |
| `hmmratac -f FRAG -i frag.tsv -n h` | insufficient training regions | insufficient training regions | insufficient training regions |

All three HMMRATAC runs pass the affected pileup stage, then raise
`Not enough training regions!`. The tiny synthetic input cannot establish
an end-to-end HMM training success. Its affected public API,
`PETrackII.pileup_bdg`, is tested directly with exact expected coverage.

## Regression coverage

`test_singleton_frag.py` checks:

- Empty, singleton, and two-record arrays with sorted and unsorted starts,
  counts 1, 7, and 65535, scaling, and a nonzero baseline.
- Read-only and reversed views, and unchanged inputs after cleanup.
- FRAG parsing followed by both `pileup_a_chromosome` and `pileup_bdg`,
  checking exact coverage and preservation across repeated calls.
- Empty files, missing columns, equal coordinates, and reversed coordinates.
- The original callpeak reproduction with zero, one, and two chr2 records,
  including exact bedGraph coverage.

Native probes use temporary script files, timeouts, `sys.executable`, and
`macs3_env`. They verify that imported extension paths belong to the
parent's installation. CLI probes use the existing `run_macs3` fixture.

## Validation

- Fresh in-place build: `python setup.py build_ext --inplace` succeeded.
- Targeted run of `test_singleton_frag.py`, `test_PileupV2.py`,
  `test_PairedEndTrack.py`, and `test_Parser.py`: 375 passed, 89 skipped.
- Full suite, `python -m pytest test -q`: 4911 passed, 367 skipped,
  302 warnings in 757.58 seconds. Default slow-test exclusions apply.
- Wheel build: `python -m build --wheel --no-isolation` succeeded from
  a clean source snapshot including the pinned SIMDe submodule.
- Installed-wheel run: all 6 tests in `test_singleton_frag.py` passed.
  Installed the wheel into a temporary package directory, put that
  directory first in `PYTHONPATH`, and invoked `python -P -m pytest`
  from the clean source tree, whose `MACS3/` contained no in-place
  extensions. Confirmed the loaded extension came from the installed
  wheel. This exercises pytest's subprocess launcher as well as the
  native probes, rather than running a probe directly from a wheel.

Test tooling was supplied by `uv run --no-project --with pytest --with build
--with psutil` inside `macs3dev`; the conda environment itself was not modified.
