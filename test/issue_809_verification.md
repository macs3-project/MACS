# Ticket #809: disabled large local window without a control

- Ticket: https://github.com/macs3-project/MACS/issues/809
- Aggregate: https://github.com/macs3-project/MACS/issues/757 (item 12)
- Branch: `fix/macs/809-llocal-zero-without-control`
- Verified baseline: `3a3d85d7c638367512842036fc14a8cd891f90fd`, fetched
  from `origin/main` before verification.
- Environment: macOS x86_64, Python 3.14.6 in `macs3dev`.

## Baseline verification

Built the unchanged baseline using `python setup.py build_ext --inplace`.
Confirmed that the parent and child processes used the checkout's compiled
`PeakDetect.cpython-314-darwin.so`. Each CLI probe used the parent's
interpreter and package root, a temporary working directory, and a
30-second timeout.

Re-ran the original seed-1 BED reproduction with 2,000 reads,
`-g 1e5 --nomodel --extsize 200 --llocal 0`, and equivalent BEDPE/FRAG
inputs. The baseline script, generated inputs, individual logs, and
`results.json` are in `/private/tmp/macs809-verification/`.

| Variant | Baseline | After fix |
| --- | --- | --- |
| Original BED reproduction | Exit 1, `ZeroDivisionError` | Exit 0 |
| BEDPE or FRAG, no control, `--llocal 0` | Same error | Exit 0 |
| BED/BEDPE/FRAG, no control, both windows zero | Same error | Exit 0 |
| Same formats with `--nolambda` | Exit 0 | Exit 0 |
| Same formats with default windows | Exit 0 | Exit 0 |
| Same formats with control and `--llocal 0` | Exit 0 | Exit 0 |
| Empty BED | Exit 1, no usable alignments | Same diagnostic |
| Empty BEDPE | Exit 1, no valid fragments | Same diagnostic |
| BEDPE with nonnumeric coordinates | Exit 1, coordinate assertion | Same diagnostic |
| BED with nonnumeric coordinates | Exit 1, local-window division | Exit 0 |

The BED parser converts nonnumeric coordinates using `atoi`, accepting the
last case as a zero-position read. The fix removes its downstream
local-window division error; coordinate parsing is unchanged. Separate
regressions verify that invalid BED strands and invalid BEDPE/FRAG
coordinates still terminate with diagnostics. No probe crashed natively
or timed out. After-fix probe artifacts are in `/private/tmp/macs809-after/`.
All 18 peak/summit files from the default-window, with-control/zero-window,
and `--nolambda` probes match their baseline files byte for byte.

## Change and regression coverage

The no-control peak detector now calculates a large-window background only
when local lambda is enabled and the large window is nonzero. Otherwise it
passes empty local-window lists to the caller, retaining the existing
global background. Positive-window and with-control calculations are
unchanged. CLI help and callpeak documentation explain the no-control
`--llocal 0` behavior.

Replaced the existing test that expected `ZeroDivisionError` with single-end
and paired-end assertions on the caller's background configuration. Added
real public-API comparisons of peak fields and bedGraphs for both track
types and both p-value and q-value calling. CLI comparisons cover BED,
BEDPE, BAMPE, and FRAG; narrow and broad peaks; and enabled/disabled small
windows. They compare nonempty peak files and bedGraphs against
`--nolambda`, and independently check the global background value.
Invalid-input subprocesses have 30-second timeouts. All CLI tests reuse
`run_macs3`, which selects the parent's MACS3 installation and runs from
`tmp_path`.

Before rebuilding the fix, four public-API cases and one representative
BEDPE CLI case failed with the reported division by zero. After rebuilding
with `python setup.py build_ext --inplace --parallel 4`, the focused suite
passed: **165 passed in 65.58 s**.

## Validation commands

Commands run in `macs3dev`. `uv run --no-project` panicked during startup in
this environment, so validation used the environment's Python directly.

```sh
python setup.py build_ext --inplace --parallel 4
python -m pytest test/test_callpeak_llocal_zero.py test/test_PeakDetect.py \
  test/test_callpeak_local_windows.py -q
python -m pytest test -q
python -m build --wheel --no-isolation \
  --outdir /private/tmp/macs809-validation/wheels \
  /private/tmp/macs809-validation/clean-source
```

The clean source copy was extracted with `git archive HEAD` and overlaid
with the changed files, including the new regressions. The initial build
failed because `git archive` omits submodule contents. The SIMD header
submodule was then archived at its pinned revision,
`8cd136a43bae7ab9b82316179b9cef8887726778`, into the clean source copy.
The wheel is built
from that copy, installed into a separate venv, and tested from the same
source copy without in-place extensions. The wheel-test launcher asserts
that MACS3, PeakDetect, and Parser load from the venv, then runs pytest
against the real test files. The launcher and build artifacts are in
`/private/tmp/macs809-validation/`.

The wheel build succeeded:
`macs3-3.1.0-cp314-cp314-macosx_11_0_x86_64.whl`.
It was installed with `pip install --no-deps --ignore-installed` into a
venv using `macs3dev`'s dependencies via `--system-site-packages`.
From the clean source directory, with `PYTHONPATH` empty:

```sh
/private/tmp/macs809-validation/venv/bin/python \
  /private/tmp/macs809-validation/run_wheel_tests.py
```

This ran the three focused test modules plus `test_callpeak_cmd.py`:
**407 passed, 29 skipped in 199.01 s**. The package and compiled modules
were confirmed under the venv's `site-packages`; no in-place extensions
were present in the clean source package. Build, installation, and pytest
logs are `/private/tmp/macs809-wheel-build-retry.log`,
`/private/tmp/macs809-wheel-install.log`, and
`/private/tmp/macs809-wheel-tests.log`.

The full in-place suite completed successfully: **5,226 passed, 468 skipped,
302 warnings in 1,046.84 s (17:26)**. The warnings come from hmmlearn's use
of NumPy's deprecated array-shape assignment. The full log is
`/private/tmp/macs809-full-suite.log`. `git diff --check` also passed.
