# Short BAM verification

Tickets: [#801](https://github.com/macs3-project/MACS/issues/801),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

## Baseline

Verified against compiled current-main revision
`1fb5318e57a80081a33c927d9a6e3fda2a6f1485` on macOS with Python 3.14.6,
Cython 3.3.0, NumPy 2.5.1, and pysam 0.24.1. The baseline was built from
an unmodified Git archive with the pinned SIMDe submodule populated.

The original pysam reproduction was repeated for every record count from
0 through 11. Public `BAMParser.tsize()` and `filterdup -f BAM` / `-f AUTO`
fail with `struct.error: unpack requires a buffer of 4 bytes` for 0–9
records. Counts 10 and 11 succeed. Explicit BAM with `-s 50` succeeds for
1–11 records; header-only BAM instead reaches `ZeroDivisionError` in
filterdup's redundant-rate calculation.

Empty streams, truncated headers, partial record-size fields, and partial
record bodies were also probed. All were rejected, but tag-size inference
reported low-level unpack errors for incomplete records. The track loader
already reports `EOFError` for partial record sizes and bodies when `-s 50`
bypasses inference. No probes crashed or timed out.

The final new regression module against the baseline reports
**39 failed, 16 passed**. These failures cover short-file inference and AUTO
detection, header-only diagnostics, truncated-record diagnostics, and invalid
record sizes.

## Change and coverage

`BAMParser.tsize()` averages the available records up to the existing limit
of ten, stopping only at a clean record boundary. Header-only BAM raises
`ValueError` with a "No alignments" message. Partial size fields and bodies
raise `EOFError`; record sizes smaller than the 32-byte BAM core are
rejected before reading the body. Successful inference still rewinds the
stream and caches the integer-truncated mean.

When duplicate filtering receives an empty track, `filterdup` now raises a
clear "No alignments" error before its zero-denominator calculation. This
covers header-only BAM with explicit tag size. The `--keep-dup all` path
retains its existing behavior.

`test_bam_short.py` covers raw and BGZF files with 1–12 records, integer
means and the first-ten limit, format detection, stream rewinding, caching
after track loading closes the input, exact track and BED coordinates,
explicit tag size, header-only files, and raw/gzip malformed inputs.
The original nine-record CLI command is compared with its explicit-size
control. Automatic duplicate filtering's existing small-sample cutoff is
preserved; coordinate tests use `--keep-dup 1` to retain distinct reads.

Every subprocess has a 30-second timeout and uses the parent's interpreter
and MACS3 installation through `macs3_env`. API probes additionally assert
that the child imports exactly the parent's compiled parser extension.

## Validation

- In-place build: `python setup.py build_ext --inplace -j 4` succeeded.
- Focused regression module: **55 passed**.
- In-place parser/filterdup suite: **483 passed** in 85.13 seconds:
  `python -m pytest test/test_bam_short.py test/test_Parser.py
  test/test_Parser_bam_truncation.py test/test_filterdup_cmd.py
  test/test_Parser_imports.py -q`.
- Clean-source wheel: `python -m build --wheel --no-isolation` succeeded,
  including the pinned SIMDe submodule.
- Installed-wheel run: **483 passed** in 66.27 seconds for the same five
  modules. Installed with `pip --no-deps --no-index --target` into an isolated
  temporary directory. Ran `python -P -m pytest` from the clean source tree,
  which contained no in-place extensions, with the installed wheel directory
  on `PYTHONPATH`. Confirmed the imported parser belongs to that installation;
  child API probes assert the exact same compiled extension path.
- `git diff --check` passed. The full repository suite was not run.

Tests use the existing temporary `macs800-verification/test-env` interpreter,
which inherits the `macs3dev` scientific dependencies and includes pytest,
pysam, and build tooling. Plugin autoload is disabled. Baseline extensions
and the clean wheel source are isolated in temporary directories; the
in-place run imports extensions from this worktree.
