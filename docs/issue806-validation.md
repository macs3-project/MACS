# BAMPE validation: #806

Issue: https://github.com/macs3-project/MACS/issues/806

Parent report: https://github.com/macs3-project/MACS/issues/757

## Verification before the fix

Verified revision `90b61de74af2ddd836283753d1c18aacee9f4b11` on macOS
with Python 3.14 in `macs3dev`. Rebuilt the in-place Cython extensions using
`conda run -n macs3dev uv run --no-project python setup.py build_ext --inplace`.
The loaded parser was `MACS3/IO/Parser.cpython-314-darwin.so`.

The equivalent reproduction writes small BGZF BAM files directly, without
requiring pysam. Each file contains 20 proper-pair first-mate records with
50M CIGARs, forward position 100, mate position 250, and the indicated BAM
core fields. The normal TLEN is 200; the header declares chr1 and chr2.
Public parser and CLI probes run in temporary directories with timeouts
and the parent's MACS3 installation.

| Input | Result before the fix |
| --- | --- |
| Valid same-reference TLEN 200 | 20 fragments, mean 200, total length 4000 |
| TLEN 0 | 20 fragments, mean 0, total length 0 |
| Mate on chr2 | Accepted as 20 chr1 fragments |
| Mate reference -1 or outside header | Accepted as 20 chr1 fragments |
| Header only | `ZeroDivisionError: float division` |
| TLEN -2147483648 | Negative mean; unsigned total length wraps |
| TLEN 2147483647 at position 100 | Fragment endpoint overflows |
| Negative read or mate position | Accepted |

Running `callpeak -f BAMPE -t zero_tlen.bam -g 1e6 -n x` exits 1
with `AssertionError: Lambda must > 0, however we got 0`, reproducing
the original report on current main.

## Fix and regression coverage

The BAMPE loader skips invalid references, missing positions, zero TLEN,
and template lengths whose absolute value or inferred endpoint cannot fit
signed 32-bit coordinates. Valid negative TLEN remains supported. Both
track construction and append reject files with no accepted fragments,
using `ValueError` with the input filename, consistent with BEDPE's
existing empty-input convention. `callpeak` logs that error and exits 1.
Peak detection also rejects unusable paired-end tracks and a nonpositive
background rate before invoking the Poisson scorer.

`test/test_Parser_bampe_validation.py` exercises public build and append
APIs through timeout-bounded subprocesses for raw, gzip and BGZF input.
It covers invalid-only and mixed files, empty and flag-filtered input,
positive/negative TLEN, the maximum representable endpoint, second-mate
filtering, and preservation of an existing track on a rejected append.
CLI tests cover empty/invalid treatment and control and compare peak,
summit and bedGraph output from valid and mixed input. Tests in
`test/test_PeakDetect.py` exercise the public peak-calling guards.

## Validation

- Full in-place suite: **5067 passed, 468 skipped** with
  `conda run -n macs3dev uv run --no-project python -m pytest test -q`.
  The run took 874 seconds and reported 302 NumPy deprecation warnings
  from hmmlearn.
- Clean-source wheel build: passed with `python -m build --wheel
  --no-isolation` in `macs3dev`.
- Installed-wheel regression run: **551 passed, 128 skipped**. The source
  snapshot had no in-place extensions. The parent pytest process loaded
  `MACS3/IO/Parser.cpython-314-darwin.so` from the temporary environment's
  `site-packages`; CLI and native probes used `macs3_env` and ran from
  temporary directories. This run included `test_Parser_bampe_validation.py`,
  `test_Parser_bam_truncation.py`, `test_Parser.py`, `test_PeakDetect.py` and
  `test_callpeak_cmd.py`.
