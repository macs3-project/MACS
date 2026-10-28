# Reverse-strand SAM CIGAR verification

Tickets: [#800](https://github.com/macs3-project/MACS/issues/800),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

## Baseline

Verified against compiled revision
`7b739c56af5b7fa73866808ab3447547de9090b5` on macOS with Python 3.14.6,
Cython 3.3.0, and NumPy 2.5.1 in `macs3dev`.

The original reproduction, with 2,000 seeded forward alignments and one
flag-16 alignment, fails in `filterdup -f SAM -g 1e6` with
`TypeError: an integer is required`. The forward-only control exits 0.
Public `SAMParser.build_fwtrack` and `append_fwtrack` calls also fail for
reverse alignments, including flag 83 (reverse first mate of a proper pair)
and CIGARs combining reference-consuming operations with clipping,
insertions, and padding. Flag 147 (second mate) remains filtered.

Empty and blank-only inputs produce empty tracks. A header-only SAM file
raises `IndexError` because the header skipper rewinds to the last header.
Truncated records raise `IndexError`. Missing (`*`) and unrecognized CIGARs
currently contribute zero reference span; this permissive behavior is
preserved rather than introducing new validation rules.

The updated regression tests run against the baseline extension report
**20 failed, 256 passed**. Every failure concerns reverse SAM parsing or
header-only handling.

## Change and coverage

Convert matched CIGAR length bytes with `int(x)` before summing reference
operations M/D/N/=/X. Rewind the header reader only when an alignment line
has actually been found, leaving empty and header-only streams at EOF.

`test_Parser.py` now enables the previously excluded reverse SAM flag and
CIGAR cases, includes reverse reads in format detection, and checks exact
SAM/BAM coordinate agreement for mixed CIGAR operations and paired flags.

`test_sam_reverse.py` exercises append loading for single reads and reverse
first mates, plain and gzip input, empty/header-only/blank inputs, truncated
records, and missing or invalid CIGARs. The CLI regression repeats the
original seeded reproduction, checks the reverse BED interval
`chr1:4999-5049`, and compares all forward output to the forward-only control.

Subprocess probes use 30-second timeouts, `sys.executable`, the existing
`macs3_env` helper, and temporary working directories. Public API probes
verify that their compiled parser belongs to the parent's installation.

Exploratory input without a final newline exposed the shared text loader's
EOF loop. This is the separate defect tracked by
[#804](https://github.com/macs3-project/MACS/issues/804); it affects SAM as
well as the formats listed there. This change keeps that loader untouched.

## Validation

- In-place build: `python setup.py build_ext --inplace -j 4` succeeded.
- Focused in-place run: `python -m pytest test/test_sam_reverse.py
  test/test_Parser.py -q` — **276 passed**.
- Full suite: **5,205 passed, 95 skipped, 20 failed** in 863.38 seconds.
  Nineteen failures were existing worker-lifecycle probes whose subprocesses
  replace `PYTHONPATH` and therefore could not import pytest supplied only
  through the parent's `PYTHONPATH`. The remaining failure is the baseline
  `callvar` failure described below.
- Corrected-environment follow-up: `python -m pytest
  test/test_thread_lifecycle.py test/test_sam_reverse.py test/test_Parser.py
  -q` — **295 passed** in 36.18 seconds. All nineteen worker-lifecycle
  failures are resolved by the temporary test environment; no repository
  changes were needed. The full suite was not repeated after this
  environment-only correction.
- Clean-source wheel build: `python -m build --wheel --no-isolation`
  succeeded, including the pinned SIMDe submodule.
- Installed-wheel run: **276 passed** for the same two test modules.
  Installed the wheel into a temporary package directory, put it first in
  `PYTHONPATH`, and invoked `python -P -m pytest` from the clean source tree,
  which contained no in-place extensions. Confirmed the parser loaded from
  the installed wheel. This exercises the pytest subprocess launchers as
  well as the public APIs and CLI, rather than running a probe from an
  extracted wheel alone.
- The unrelated existing test
  `test_callvar_cmd.py::test_callvar_fermi_overlap_longer_than_reads_aborts`
  fails in isolation on both this build and a pre-fix installed wheel whose
  SAM parser still contains the original cast. It expects SIGABRT for an
  overlap longer than the reads, but receives exit 1 with
  `TypeError: 'int' object is not subscriptable` in the assembly path.
  The relevant `RACollection.py` sources are identical between the builds.

Initial test tooling is supplied from an existing temporary uv environment
via `PYTHONPATH`, using the `macs3dev` interpreter and its scientific
dependencies. For the worker-lifecycle follow-up, a temporary virtual
environment uses `macs3dev` system packages and contains copies of the cached
test tooling, making pytest available even when a child resets `PYTHONPATH`.
The conda environment is unchanged. Plugin autoload is disabled for pytest;
the repository's `conftest.py` fixtures and default slow-test exclusions apply.
