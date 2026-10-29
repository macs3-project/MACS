# Ticket #807: fragments longer than local background windows

- Ticket: https://github.com/macs3-project/MACS/issues/807
- Aggregate: https://github.com/macs3-project/MACS/issues/757 (item 10)
- Branch: `fix/macs/807-long-fragment-local-windows`
- Verified baseline: `f37383588b435d934516810acace8f98828ef775`, matching
  the local `origin/main` reference at verification time.
- Environment: macOS x86_64, Python 3.14.6, Cython 3.3.0, NumPy 2.5.1.

## Baseline verification

Built the unchanged baseline with `python setup.py build_ext --inplace
--parallel 4`. Confirmed `MACS3.Signal.PeakDetect` loaded the in-place
`PeakDetect.cpython-314-darwin.so` extension.

Re-ran the ticket's seed-1 reproduction with 1,000 treatment fragments of
1,500 bp and 1,000 control fragments of 200 bp, using `-f BEDPE -g 2e6`.
Every CLI probe used the same interpreter and package root as the parent,
a temporary working directory, and a 30-second timeout.

| Variant | Exit | Result |
| --- | --- | --- |
| Original reproduction | 1 | `AssertionError: 1000 can't be smaller than 1500.0!` |
| Original plus `--nolambda` | 1 | Same assertion |
| Explicit `--slocal 1500 --llocal 10000` | 0 | Completed |
| `--slocal 0` | 0 | Completed |
| `--slocal 0 --llocal 1000` | 1 | Large-window assertion |
| `--slocal 0 --llocal 0` | 0 | Completed |
| No control | 0 | Completed |
| BED, `--nomodel --extsize 1500`, with control | 1 | Small-window assertion |
| Same BED run plus `--nolambda` | 1 | Same assertion |
| Empty BEDPE treatment | 1 | `No valid fragments found in BEDPE input` |
| Malformed BEDPE treatment | 1 | Parser assertion about fragment coordinates |

No probe crashed natively or timed out. Logs, inputs, results and the probe
script are under `/private/tmp/macs807-verification/`.

The control path in `PeakDetect` asserts that `d <= slocal <= llocal`
(for enabled windows) before checking `nolambda`. Paired-end `d` is the
mean treatment fragment length; single-end `d` is the modeled or explicit
extension length. The underlying control pileup uses these windows and
their matching `d / window` scaling factors.

## Implementation

1. Replaced the fragment/window assertions with runtime handling in the
   peak detector, where the final fragment length is available to both CLI
   and public-API callers. An enabled window shorter than `d` is enlarged
   to at least `ceil(d)`, with a warning naming the option, old size and
   new size. Disabled (`0`) windows are preserved. Negative or inverted
   enabled window settings raise an explicit `ValueError`.
2. Bypassed local-window calculations and their size/order checks when
   `--nolambda` is active. Control scaling factors use the effective window
   sizes. Ordinary valid-input calculations are preserved.
3. Added public-API/CLI regressions for paired-end and single-end long
   fragments, fractional paired-end means, lengths exceeding both windows,
   exact boundaries, disabled windows, `--nolambda`, and malformed window
   settings. Background bedGraphs and peak records are compared against
   runs with explicitly enlarged windows, in both scaling directions.
   Tests cover BED, BEDPE, BAMPE and FRAG inputs, plus bounded
   empty/malformed input probes. Updated existing assertion-based tests
   for the new behavior.
4. Updated CLI help and callpeak documentation to explain automatic window
   enlargement and warnings. When a window is adjusted, the peaks XLS
   header records both effective local window sizes in addition to the
   original command and arguments.

## Fixed reproduction and compatibility

All nine valid-input probes in the baseline table now exit 0, including
the original reproduction, `--nolambda`, the smaller large-window case,
and both single-end cases. Empty and malformed BEDPE inputs retain their
previous errors, without native crashes or timeouts.

For the four baseline configurations that already succeeded (explicit
windows, disabled small window, both windows disabled, and no control),
peak and summit files are byte-for-byte identical before and after the
fix. Fixed probe logs and outputs are under
`/private/tmp/macs807-verification/fixed/`.

## Validation

The focused new regression and peak-detector suite passed **138 tests**.
The final expanded in-place suite passed **423 tests in 209.56s**, with no
skips:

```sh
python -m pytest test/test_callpeak_local_windows.py test/test_PeakDetect.py \
  test/test_callpeak_cmd.py test/test_Parser_bampe_validation.py --run-slow -q
```

Optional BAM coverage ran with pysam 0.24.1 installed in a temporary
dependency directory. The entire repository suite was not run.

Built `macs3-3.1.0-cp314-cp314-macosx_11_0_x86_64.whl` from a fresh
`git archive HEAD` with changed sources and tests overlaid, and the pinned
SIMDe submodule copied from the local build. Installed it with
`pip install --no-deps --no-index --target` into a dedicated temporary
directory. Verified that the wheel sources and regression tests match
the working tree, and that the clean source checkout contains no `.so`
files. The wheel-test runner checks that MACS3, Parser and PeakDetect load
from the installed wheel, with Parser and PeakDetect compiled.

Ran the same four-module selection with `--run-slow -q` from that clean
checkout using `python -P`, the installed wheel first on `PYTHONPATH`,
and pytest plugin autoload disabled. Result: **423 passed in 192.85s**,
with no skips. This includes every changed subprocess test. `git diff
--check` also passed.

CLI subprocesses use `run_macs3`, which selects `sys.executable` and the
parent's imported MACS3 package root, caps thread pools, disables core
dumps, and enforces a timeout. New CLI probes use a 30-second timeout.

The Sphinx HTML documentation build succeeded with 233 repository
warnings; none reference the modified `callpeak.md` page. Logs, scripts,
wheel and build outputs are under `/private/tmp/macs807-verification/`.
Any eventual fix PR should link both #807 and #757.
