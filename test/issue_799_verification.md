# Empty chromosomes after down-sampling

Tickets: [#799](https://github.com/macs3-project/MACS/issues/799),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

Baseline: `8db88cb31b1b49b4f0ea79db2c4150f6c3e66dd3` (current `main`
when work began), macOS x86_64, Python 3.14.6 in `macs3dev`, with a fresh
`python setup.py build_ext --inplace` build.

## Baseline findings

The original BED reproduction exits 0 with and without `--down-sample`,
including with treatment and control swapped. Equivalent BEDPE and BAMPE
inputs also exit 0 in all four combinations: 12 successful CLI runs.
The BED data generator and options reproduce the ticket exactly. Paired-end
fixtures use 1000 versus 200 chr1 fragments and one chr2 fragment per file,
with 200 bp fragments and 50 bp paired BAM reads.

Peak calling already handles empty chromosomes: `CallPeakUnit` checks
`n_locations` before requesting treatment and control pileups. These guards
were introduced by `90c31505d55eedc7e05219579b16fc705fd31956`.

The public pileup APIs remain affected. Both empty
`pileup_from_PN_shifted` input and empty `pileup_from_LR_as_list` input
raise `Exception: length is 0`. A plus-only or minus-only input succeeds;
the failure requires the combined interval set to be empty.

After `sample_percent(0.2, seed=7)` leaves chr2 empty while retaining two
chr1 reads/fragments, both `FWTrack` and `PETrackI` still raise the same
exception in `pileup_a_chromosome` and `pileup_a_chromosome_c`. Control
projections with one and three windows both reproduce the failure.

Running the new regression module against the original compiled extension
produced **2 failed, 13 passed**: the legacy pileup and sampled-track probes
failed with `length is 0`; malformed-field and CLI tests passed.

## Change and regression coverage

The shared legacy unit-weight sweep now returns empty int32 positions and
float32 values before accessing coordinate pointers. Sampling, nonempty
coverage, scaling, and baseline behavior retain their existing algorithms.

`test_downsample_empty.py` checks:

- Empty and singleton intervals, empty combined strands, plus-only and
  minus-only input, exact coverage, scaling, and nonzero baselines.
- Repeated calls, unchanged inputs, and read-only fragment arrays.
- Sampled-empty chromosomes in both track types, treatment pileups,
  one/two/three control windows, and retained coverage on another chromosome.
- Missing interval fields and mismatched control windows remain rejected.
- The original BED CLI reproduction and equivalent BEDPE/BAMPE cases,
  both treatment/control directions, with and without down-sampling.

Native probes use temporary scripts, `sys.executable`, `macs3_env`, a
30-second timeout, and disabled core dumps. They assert that imported
extension paths belong to the parent's installation. CLI tests use
`run_macs3` with 60-second timeouts. BAMPE fixtures require pysam.

These tests keep reads on another chromosome; the separate all-reads-removed
sampling behavior is tracked by #808.

## Validation

- Fresh in-place build: `python setup.py build_ext --inplace` succeeded.
- Focused in-place run: **585 passed**, no skips, in 147.33 seconds:
  `python -m pytest test/test_downsample_empty.py test/test_PileupV2.py
  test/test_FixWidthTrack.py test/test_PairedEndTrack.py
  test/test_callpeak_cmd.py test/test_singleton_frag.py -q`.
- Wheel build: `python -m build --wheel --no-isolation` succeeded from a
  clean source snapshot containing the changed source/test and pinned SIMDe
  revision `8cd136a43bae7ab9b82316179b9cef8887726778`.
- Installed-wheel run: **21 passed**, no skips, in 16.02 seconds:
  `python -P -m pytest test/test_downsample_empty.py
  test/test_singleton_frag.py -q`. Installed the wheel with
  `pip install --no-deps --no-index --target <temporary-directory>`, put
  that directory first in `PYTHONPATH`, and ran pytest from the clean
  source tree, whose `MACS3/` contained no in-place extensions. Confirmed
  the parent's extension came from the installed wheel; native probes also
  assert their extension paths belong to the parent's package directory.
- `git diff --check` passed.

Test tooling (pytest, build, psutil, and pysam) was installed in a temporary
uv environment. Its site-packages directory was added to `PYTHONPATH` when
using the `macs3dev` interpreter; the conda environment was not modified.
