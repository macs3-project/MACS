# Missing-SEQ BAM tag-size verification

Issue: [#810](https://github.com/macs3-project/MACS/issues/810),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

Branch: `fix/macs/810-bam-missing-seq-tag-size`.
Baseline: current `origin/main`,
`6c1611b831bee9e7c3fa8151278e1d10cc00687c`.

## Baseline verification

Built an unmodified Git archive with `python setup.py build_ext --inplace
-j 4`, using the macs3dev-based test interpreter (Python 3.14.6,
Cython 3.3.0, NumPy 2.5.1, pysam 0.24.1) and the locally available SIMDe
headers. Confirmed the probe imports the compiled parser extension.

The original 2,000-record alternating-strand `50M`/missing-SEQ reproduction
still returns `BAMParser.tsize() == 0`:

| Command | Baseline result |
| --- | --- |
| `filterdup -f BAM` | Exit 1; `FWTrack` positive-width assertion |
| `randsample -f BAM -p 50` | Exit 1; same assertion |
| Either command with `-s 50` | Exit 0; size 50 |
| `filterdup -f AUTO` | Exit 1; `File is not of a valid BAM format! 0` |
| `pileup -f BAM` | Exit 0; reports size 0 |
| `callpeak -f BAM --nomodel` | Exit 0; reports size 0 |

Additional missing-SEQ probes with a complex CIGAR, absent CIGAR, and
`0M` produce the same zero-size behavior across these four commands.
Ordinary `50M` reads with SEQ report size 50 and all four commands succeed.
Header-only BAM raises `ValueError: No alignments found in BAM file.`
through `tsize()`; all four commands reject it. These results verify the
reported failure and zero-size propagation, not the downstream scientific
effect on called peaks. Every CLI invocation has a 30-second timeout and
uses the same interpreter and compiled package as the parent probe.

Existing `test/test_bam_short.py`: **54 passed, 1 failed** (20.15 seconds).
Malformed framing and short-file tests pass. The existing failure is a
message assertion for header-only BAM with explicit `-s`: the test expects
`No alignments`, but current main returns `No usable alignments remain`.
It predates this work and does not indicate an accepted empty input.

Local evidence: `/private/tmp/macs810-verification/` contains the build log,
`verify.py`, `results.json`, individual CLI logs, and `short-bam-tests.log`.
These temporary artifacts are not committed fixtures.

## Change and coverage

`BAMParser.tsize()` retains positive `l_seq` values. When `l_seq` is zero,
it infers query length from `M`, `I`, `S`, `=`, and `X` CIGAR operations.
Deletions, reference skips, hard clips, and padding do not consume query
bases and do not contribute. This estimates the same quantity as `l_seq`;
using reference span would inflate tag size for spliced reads. The
first-ten-record sample, integer mean, cache, and successful rewind remain.

Inference rejects negative `l_seq`, missing read names, out-of-bounds CIGAR
arrays, unknown or zero-length operations, and lengths exceeding `INT_MAX`.
If a sampled missing-SEQ record has no query bases, it raises `ValueError`
instead of including zero in the average, even in a mixed sample. The
message describes how to supply usable lengths and suggests `-s/--tsize`
conditionally, since `pileup` has no such option. Explicit `-f BAM -s 50`
continues to bypass inference. AUTO detection now recognizes inferable
missing-SEQ BAM files. No command options or track-loading code changed.

`test_bam_missing_seq.py` uses small synthetic BAM records without optional
pysam dependencies. It covers raw/gzip/BGZF input, both public parser
factories, mixed SEQ presence, complex CIGARs, both strands and exact 5'
coordinates, short inputs, the sampling boundary, integer means, caching,
stream reuse, header-only input, and malformed/unusable lengths. CLI tests
repeat the original 2,000-record reproduction for `filterdup`, `randsample`,
`pileup`, and `callpeak --nomodel`, in BAM and AUTO modes. They compare
outputs with explicit-size controls, check exact BED coordinates and sample
counts, and verify actionable errors and explicit-size recovery. Every
subprocess has a 30-second timeout and uses `sys.executable` and `macs3_env`;
API probes additionally assert the exact parent's compiled extension path.

Against unmodified compiled main, the new module reports **18 failed,
3 passed**; the passes are the explicit-size bypass controls.
The existing empty-input test assertion now matches current main's
`No usable alignments` diagnostic when explicit BAM size bypasses inference.
The older parser test expecting missing-SEQ BAM detection to fail now
checks successful detection, the inferred size, and preservation of reads.

SAM missing-SEQ inference, new CLI options, and BAMPE fragment loading are
outside this change. Existing empty/truncated-input tests remain in the
validation suite.

## Validation

- In-place build: `python setup.py build_ext --inplace -j 4` succeeded.
- Final changed test modules (`test_bam_missing_seq.py`, `test_bam_short.py`,
  `test_Parser.py`): **335 passed** in 38.84 seconds in-place.
- Clean-source wheel: `python -m build --wheel --no-isolation` succeeded.
  Installed with `pip --no-deps --no-index --target` into a fresh temporary
  directory. The clean source's `MACS3/` tree contains no in-place extensions.
  The runner uses `python -P`, asserts the parent's parser is the installed
  compiled extension, and launches pytest against that source checkout.
- Final same three modules against the installed wheel: **335 passed** in
  37.49 seconds. This validates every changed subprocess test with the
  parent's installation selected by `macs3_env`.
- Broader wheel coverage (`test_bam_missing_seq.py`, `test_bam_short.py`,
  `test_Parser.py`, `test_Parser_bam_truncation.py`,
  `test_Parser_bampe_validation.py`, `test_Parser_imports.py`,
  `test_filterdup_cmd.py`, `test_randsample_cmd.py`, `test_pileup_cmd.py`,
  `test_callpeak_cmd.py`): **993 passed, 1 failed** in 250.05 seconds.
  This run collected the older parser assertion expecting detection to fail
  before that assertion was updated. The corrected full parser module then
  passed all **259** tests against the wheel; the final 335-test run above
  also includes that correction.
- Full in-place suite (`python -m pytest test -q`): **5618 passed,
  95 skipped, 2 failed** in 941.51 seconds. One failure was the old parser
  assertion collected before it was corrected; the final 335-test runs in
  both environments above pass its replacement. The other failure is the
  pre-existing callvar case confirmed separately below. No other failures
  occurred. The complete suite was not repeated after the test-only
  assertion correction.
- Separate baseline confirmation: the full suite's
  `test_callvar_fermi_overlap_longer_than_reads_aborts` failure reproduces on
  unmodified compiled revision `6c1611b831bee9e7c3fa8151278e1d10cc00687c`.
  Both baseline and fixed builds exit 1 with
  `TypeError: 'int' object is not subscriptable`, whereas the test expects
  `SIGABRT`. This pre-existing callvar failure is outside #810. Evidence is
  in `callvar-baseline.log` and `callvar-inplace.log`.
- `git diff --check` passed.

All runs use the existing `macs800-verification/test-env` interpreter,
which inherits the `macs3dev` scientific dependencies. Pytest plugin
autoload is disabled. Logs, the isolated baseline build, the clean wheel
source, wheel installation, and validation runner are retained under
`/private/tmp/macs810-verification/`.
