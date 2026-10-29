# Empty tracks and sampled fragment sets: #808

Ticket: https://github.com/macs3-project/MACS/issues/808

Parent report: https://github.com/macs3-project/MACS/issues/757 (item 11)

## Baseline verification

Verified `f37383588b435d934516810acace8f98828ef775` on macOS x86_64,
using Python 3.14 in the `macs3dev` conda environment. This revision matched
local `main` and `origin/main`. Rebuilt the checkout with
`python setup.py build_ext --inplace` before rerunning the reproductions.
Parser, PairedEndTrack and HMMR_EM loaded their compiled `.so` extensions.

Public-API and CLI probes used temporary working directories, the parent
Python interpreter and MACS3 package, and 30-second subprocess timeouts.
Small synthetic BGZF BAM fixtures were written directly, without pysam.

| Reported case | Result on the rebuilt baseline |
| --- | --- |
| Empty BED treatment, default duplicate filtering | `callpeak` exits 1 with `ZeroDivisionError` |
| Empty BED control | Exits 1 with a misleading chromosome mismatch |
| Empty BED treatment, `--keep-dup all --nomodel` | Exits 1 with a downstream background-lambda `ValueError`; the old report's success behavior no longer occurs |
| Empty BED, `filterdup` default | Already raises `ValueError: No alignments found in input files` |
| Empty BED, `filterdup --keep-dup all` | Fails the positive tag-width assertion |
| BAM with improper pairs, R2-only records, or unmapped mates, read as BAM | `callpeak` exits 1 with `ZeroDivisionError` |
| Same BAMs read as BAMPE, build and append | Already reject input with `ValueError: No valid fragments found in BAMPE input` |
| One pair per chromosome, `randsample -f BAMPE -p 50` | Exits 1 with `ZeroDivisionError` |
| Both `PETrackI.sample_percent` methods round all chromosomes to zero | `ZeroDivisionError` |
| Sparse BAMPE, default `hmmratac` EM sampling | `ZeroDivisionError` |
| 200 BAMPE fragments, all 2 kbp, default `hmmratac` EM | `ZeroDivisionError` after length selection |
| Empty `PETrackI.finalize()` | Additional related `ZeroDivisionError` |

The baseline parser-validation selection and existing EM empty-range test
passed: **7 passed, 42 deselected**. The latter explicitly expected the
division-by-zero defect and has been updated.

## Changes and regression coverage

Empty paired-end tracks now have zero total, length and mean fragment
length. Percentage sampling keeps the existing per-chromosome rounding;
sampling to zero is valid. Number sampling, finalization, duplicate
filtering and fragment-length extraction also handle empty tracks.
`PETrackII` can return an empty barcode subset or blacklist result;
subset reference metadata is copied so dropping empty chromosomes does
not remove the original track's reference lengths.

Commands requiring alignments reject empty tracks with a logged diagnostic
and exit status 1. The checks run before chromosome-name comparisons,
redundancy calculations or assigning a zero tag width. BAM diagnostics
explain the proper-pair flag requirement. Existing BAMPE parser rejection
is retained. An intentionally empty `randsample` result succeeds and writes
an empty BED file.

EM checks for a zero sample before extracting lengths, and checks again
after length selection. The errors distinguish these causes and suggest
changing sampling/range parameters or using `hmmratac --no-fragem` with
suitable fragment parameters.

Coverage includes:

- `test_empty_tracks.py`: bounded public-API lifecycle/sampling probes,
  empty treatment/control, duplicate modes, empty BEDPE/FRAG files,
  BAM flag filtering, intentionally empty random samples, both EM failure
  cases, unmatched barcodes and blacklist removal of all fragments.
- `test_Parser_bampe_validation.py`: build and append on improper pairs,
  R2-only and mate-unmapped input, alongside existing empty, malformed,
  mixed-invalid and valid BAM coverage in raw, gzip and BGZF formats.
- `test_HMMR_EM.py`: empty/sparse/zero-percentage samples and lengths
  below and above the accepted range, alongside existing valid fits.
- `test_PairedEndTrack.py`: empty FRAG finalization and barcode subsets,
  including preservation of the original track and reference lengths.
- `test_filter_dup_singletons.py`: empty CLI cases now require the
  descriptive error; the malformed BED case is unchanged.
- `test_chromosome_names.py` and `test_hmmratac_cmd.py`: older expectations
  for empty-track chromosome mismatches and unmatched-barcode assertions
  now require graceful empty-track diagnostics. Nonempty disjoint-name
  cases retain their chromosome-specific errors.

All subprocess tests use `run_macs3` or `macs3_env` from `conftest.py`,
selecting the parent's installed MACS3 package and `sys.executable`.
Native/API probes run scripts in `tmp_path`, with a 30-second timeout.

## Validation

- Rebuilt the final in-place extensions with
  `python setup.py build_ext --inplace` in `macs3dev`.
- Final full in-place suite: **5150 passed, 468 skipped** in 765.83 seconds
  with `python -m pytest test -q`. The run reported 302 existing NumPy
  deprecation warnings from hmmlearn.
- Wheel build passed with `python -m build --wheel --no-isolation`.
  `uv` crashed during startup in this sandbox, so validation used the
  conda environment's Python directly.
- Final installed-wheel regressions: **482 passed, 18 skipped** in
  185.41 seconds. The selection included `test_empty_tracks.py`,
  `test_Parser_bampe_validation.py`, `test_HMMR_EM.py`,
  `test_PairedEndTrack.py`, `test_filter_dup_singletons.py`,
  `test_chromosome_names.py` and `test_hmmratac_cmd.py`.
- The wheel was installed with pip into a temporary virtual environment.
  Tests ran from a clean source snapshot containing the updated tests and
  no in-place `.so` extensions. Before collection, the runner verified
  that MACS3, Parser, PairedEndTrack and HMMR_EM loaded from that virtual
  environment's `site-packages`. Subprocesses selected the same package
  through `macs3_env`.
- The first full run exposed four old assertions for the intentionally
  replaced errors: three empty-BED chromosome mismatch cases and one
  unmatched-barcode assertion. All four corrected cases passed in both
  build layouts. The final full-suite run uses the corrected tests.

Local build and test logs are in `/private/tmp/macs808-*.log`; the clean
source snapshot, virtual environment and wheel are under
`/private/tmp/macs808-validation/`.
