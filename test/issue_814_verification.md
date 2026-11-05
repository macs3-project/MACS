# Degenerate HMM model verification

Issue: [#814](https://github.com/macs3-project/MACS/issues/814),
part of [#757](https://github.com/macs3-project/MACS/issues/757).
Branch: `fix/macs/814-degenerate-hmm-models`.
Baseline: current `origin/main` at
`ebdc4417ab304b2cd9e00486800e6f349d06ef9c`.

## Compiled baseline verification

Built the unchanged HMMR_HMM and Prob extensions from this revision with
the macs3dev interpreter: Python 3.14.6, Cython 3.3.0, NumPy 2.5.1,
hmmlearn 0.3.3.
The full in-place build initially lacked SIMDe headers; supplied the
locally available headers from earlier verification work and retried.
The full `python setup.py build_ext --inplace -j 4` retry succeeded.
Each API probe runs in a separate subprocess with a 30-second timeout,
the parent's interpreter and MACS3 package, and an assertion that the
child imports the exact same compiled HMMR_HMM extension as the parent.

| Probe | Baseline result |
| --- | --- |
| Valid full-covariance model | Saved; finite posterior probabilities |
| Deliberately indefinite covariance, diagonal entries including `1e12` and `-3e-5` | Saved; decoding raises `ValueError: 'covars' must be symmetric, positive-definite` |
| Singular covariance | Saved; decoding accepts it through hmmlearn's covariance jitter |
| Asymmetric covariance | Saved; decoding accepts it |
| NaN or infinite covariance | Saved; decoding raises `ValueError: array must not contain infs or NaNs` |
| Zero transition row | Saved; decoding raises a transition normalization error |
| Empty prediction/training data | Rejected with `ValueError` |
| Prediction lengths that do not sum to sample count | Rejected with `ValueError` |
| 200 identical four-feature training observations | Training returns a model with `converged == True` and rank-deficient covariances; decoding produces finite posteriors |
| One varied training sequence of 240 observations | Training and decoding succeed with positive-definite covariances |

The intentionally invalid covariance probes modify the public `covars_`
array of an otherwise valid GaussianHMM, then call MACS3's public
`hmm_model_save` and `hmm_predict`. These verify missing validation;
they do **not** reproduce the original training failure. The reporter's
small BAM is unavailable, and the constant-data example does not reproduce
the reported decoding exception either.

Baseline code returned the fitted model without validating covariance
matrices. The CLI then logged `monitor_.converged`, assigned states, and
saved the model before decoding. In installed hmmlearn, the convergence
property can be true after a negative likelihood change or reaching its
iteration limit. It is not a covariance validity check.

Existing `test/test_HMMR_HMM.py`: **51 passed, 1 skipped**, with 196
hmmlearn/NumPy deprecation warnings, in 42.55 seconds.
The legacy prediction test is explicitly skipped for sklearn/hmmlearn
version variation.

Temporary evidence: `/private/tmp/macs814_probe.py`,
`/private/tmp/macs814-baseline-results.json`,
`/private/tmp/macs814-baseline-build.log`, and
`/private/tmp/macs814-baseline-hmm-tests.log`.

## Implementation

Added shared validation after training and at the prediction, save, and
load boundaries. It checks the number of states/features, emission
dimensions and finite values, nonnegative normalized probabilities, and
nonnegative Poisson rates. Gaussian full/shared covariance matrices must
be symmetric and pass Cholesky without jitter; diagonal/spherical
variances must be positive. The stored covariance representation is
checked because hmmlearn's public getter expands several covariance
types and can conceal invalid dimensions. Validation does not modify
parameters or depend on the convergence monitor.

Invalid models raise `ValueError` identifying the parameter (and the
state for a full covariance failure), with guidance on representative
training regions, `-l/--lower`, `-u/--upper`, and `-t/--training`.
The existing CLI error handler logs validation errors and exits with
status 1 before model saving, convergence reporting, or decoding.
Public save validation also preserves an existing destination file.
Documentation explains why
`-c/--prescan-cutoff` affects decoding candidates rather than training.

`test_HMMR_validation.py` covers invalid full/diagonal covariances,
nonfinite and malformed emission parameters, invalid state probability
arrays, invalid Poisson rates, all four public boundaries, constant-data
training, empty/malformed observations, and CLI failure before output
including `--modelonly`. A controlled final-fit fixture verifies that
the CLI rejects the reported covariance pattern before saving. This
fixture is explicitly distinct from reproducing the unavailable BAM's
training failure. The constant-data fit is an actual training example
whose rank-deficient result is now rejected.

Valid controls include a single training sequence for every supported
training covariance type (full, diagonal, tied, spherical), tiny positive
variances, zero Poisson rates, zero probability entries, and a valid model
that has not converged. Manually configured models can still infer
their feature count, and a single-state Poisson model can use flat rates.
Existing tests retain exact JSON output, trained parameter expectations,
model round trips, forward/backward
posterior comparisons, and ordinary CLI output checks. The existing
two-state saved-model test now expects rejection during load rather
than delayed failure during prediction.

Against the baseline compiled HMM extension, the initial regression set reported
**75 failed, 15 passed** in 46.07 seconds. Failures include rejection
timing and actionable-diagnostic assertions for models that hmmlearn
already rejected with a less useful error. The HMM extension was still
the unmodified baseline build. Each new training/CLI probe uses the parent's
interpreter and installation, with a 30-second timeout. The training
probe also verifies the exact compiled extension path.

## Validation

- In-place Cython build succeeded with `python setup.py build_ext
  --inplace -j 4`.
- Broader HMM suite (`test_HMMR_validation`, `test_HMMR_HMM`,
  `test_hmmratac_cmd`, `test_HMMR_EM`, `test_HMMR_Signal_Processing`,
  `test_HMMR_poisson`): **423 passed, 7 skipped** in 193.06 seconds.
- In-place CLI rejection checks and all five slow real-data controls:
  **11 passed** in 274.39 seconds. Gaussian BAMPE/BEDPE, Poisson, and
  fragment/barcode runs retain the expected peak regions or existing
  output checks. The standard-result comparisons require Jaccard > 0.99.
- Built a wheel from a Git archive plus the changed sources/tests and
  local SIMDe headers, with no preexisting extensions or build artifacts.
  Installed with `pip --no-deps --no-index --no-compile --target` into a
  separate temporary directory. The runner uses `python -P`, asserts
  that MACS3 and the compiled HMM module come from that target, and
  asserts that the clean source package has no in-place `.so` files.
- Broader installed-wheel suite with `--run-slow`: **428 passed,
  2 skipped** in 440.86 seconds.
- After the final integer-covariance JSON preservation adjustment,
  rebuilt both packages and repeated the two affected API modules,
  including every new subprocess regression: **145 passed, 1 skipped**
  in-place (61.49 seconds) and **145 passed, 1 skipped** against the
  rebuilt installed wheel (63.59 seconds). The integer covariance model
  retains byte-for-byte JSON on load/save. Checked that the installed
  HMM and CLI source files, and the clean source regression/probe files,
  match the workspace's final versions.
- `git diff --check` passes.

The broader wheel run's two skips are the explicitly disabled legacy
prediction test and an existing synthetic BAM generation test requiring
unavailable optional pysam. The final API runs skip only that legacy
prediction test. Warnings are hmmlearn/NumPy deprecations and the expected
duplicate-cluster warning from a malformed-input training control.

Further evidence: `/private/tmp/macs814-final-inplace.log`,
`/private/tmp/macs814-realdata-inplace.log`,
`/private/tmp/macs814-release-inplace.log`, and
`/private/tmp/macs814-verification/` (`final-build-wheel.log`,
`final-install-wheel.log`, `final-wheel-tests.log`,
`release-wheel-tests.log`, clean source, and the installed wheel).
