# callvar assembly and spliced-read verification

Issue: [#813](https://github.com/macs3-project/MACS/issues/813),
part of [#757](https://github.com/macs3-project/MACS/issues/757).
Branch: `fix/macs/813-callvar-assembly-edge-cases`.
Baseline: `origin/main` at `b8dc7d694ad299d2bfff5521e0d57c433fe27c04`.

## Compiled baseline

Built unmodified source with `python setup.py build_ext --inplace -j 4`,
Python 3.14.6, Cython 3.3.0, NumPy 2.5.1, and local SIMDe headers.
Pysam 0.24.1 generated sorted/indexed BAMs with MD tags against a
deterministic 20 kbp reference. Probes recorded the compiled extension
paths, used the parent's interpreter and MACS3 package, disabled core
dumps, and bounded each subprocess to 10 or 30 seconds.

| Input / option | Baseline result |
| --- | --- |
| 25 bp reads, `-F on`, overlap 30 | Empty assembly followed by `TypeError: 'int' object is not subscriptable` |
| 25 bp reads with deletion, `-F auto`, overlap 30 | Empty assembly; first-pass fallback succeeds |
| Same short inputs, overlap 20 | Assembly succeeds |
| One `25M2000N25M` read, `-F off` or `auto` | `IndexError` indexing a shortened reference |
| Same spliced read, `-F on` | Empty-assembly `TypeError` |
| Four overlapping 30 bp reads, overlap 30, `-F on` | Empty-assembly `TypeError` |
| Overlap -1 | SIGSEGV (-11) |
| Overlap 0 | Accepted despite the documented minimum of 1 |
| Empty indexed BAM or empty peak file | Exit 0, header-only VCF |
| Malformed peak coordinate | Exit 1, `ValueError` |

The default build defines `NDEBUG`, so native C assertions are disabled.
Recompiled the unchanged `mrope.c` with `-UNDEBUG` and relinked
RACollection in a separate package copy, retaining the other baseline
objects. The 25 bp forced-assembly and auto-with-deletion probes both
then reproduce **SIGABRT (-6)** at `mrope.c:230`, `mr_insert_multi`,
`len > 0 && s[len-1] == 0`. Overlap 20 succeeds in that build as well.

Existing focused tests on baseline: **230 passed, 1 skipped, 1 failed**
in 42.28 seconds. The failure is the test expecting SIGABRT; the default
build reaches the downstream empty-unitig TypeError instead.

## Changes

- Validate overlap in both CLI option handling and the assembly API,
  before input loading/native work. The supported range is **1–63**:
  fermi-lite's hash implementation asserts a maximum k-mer length of 63.
  At small valid overlaps, cap the hash prefix length so native hash
  shifts remain nonnegative. Keep a halved retry overlap at least 1.
- In fermi-lite, return no index when filtering leaves no sequence to
  insert, and return zero unitigs rather than passing that missing index
  into graph construction. The sequence buffers and rope are freed on
  this path; the caller still owns and frees its options/count buffers.
  Avoid division by zero when no k-mers meet the coverage threshold.
  No overlap is silently clamped to the read length, and mixed inputs
  can still assemble using their longer reads.
- Skip CIGAR `N` alignments in treatment and control before constructing
  callvar's reference collection, reporting the sample and count per
  peak. Keep ordinary alignments in mixed inputs. All-spliced treatment
  peaks follow the existing no-reads skip behavior. This is the
  skip-spliced-reads policy explicitly permitted by #813.
- Check reconstructed sequence length against its reference span before
  assigning it into the reference buffer. Direct API calls with a
  mismatched reconstruction raise `ValueError` instead of resizing the
  buffer silently.
- Always finish the peak after an empty assembly. Auto mode writes any
  first-pass variants; forced assembly produces no calls for that peak
  and continues to later peaks.

CLI help and the callvar documentation describe these behaviors.
The existing test requiring a native abort now requires safe completion
of both affected peaks.

## Regression coverage and validation

`test_callvar_assembly.py` and `callvar_assembly_probe.py` cover short
reads, equal overlap, reads eliminated by filtering, ambiguous bases,
mixed read lengths in treatment/control, minimum and maximum overlap,
invalid overlap including integer overflow, spliced and malformed
reference spans, all three assembly modes, mixed/all-spliced samples,
empty inputs, first-pass fallback, smaller-overlap recovery, and a
successful deletion call in a later peak following an empty assembly.
CLI results are compared to ordinary-read controls and nonassembly
first-pass calls. Existing valid-input VCF expectations are retained.

Each new native API probe runs in a subprocess with a 30-second timeout
and asserts that its compiled RACollection path exactly matches the
parent's. CLI probes use `run_macs3`/`macs3_env`, which pins the package
to the parent's imported installation, with a 30-second timeout.

- In-place focused suite (`test_callvar_assembly`, `test_ReadAlignment`,
  `test_RACollection`, `test_callvar_cmd`): **280 passed, 1 skipped**
  in 71.95 seconds. The skip is the existing slow full-real-data test;
  the ordinary two-peak upstream VCF comparison passes.
- Assertion-enabled package: all **14 native API regressions pass**;
  the original three short-read CLI probes also exit 0, including the
  auto-mode first-pass fallback and overlap-20 recovery.
- Built `macs3-3.1.0-cp314-cp314-macosx_11_0_x86_64.whl` using
  `python -m build --wheel --no-isolation` from a fresh source copy
  containing tracked sources plus the new tests and SIMDe headers, with
  no copied extension binaries or build artifacts. Installed it into a
  separate target with `pip install --no-deps --no-compile --target ...`.
  Repeated the same focused suite from the clean source checkout:
  **280 passed, 1 skipped** in 81.49 seconds. Before pytest, the driver
  asserted that MACS3, RACollection, and ReadAlignment were loaded from
  the installation target, and that the clean source package contained
  no `.so` files. New subprocess API probes checked the exact parent's
  extension path again.
- Full in-place suite (`python -m pytest -q test`): **5,730 passed,
  101 skipped** in 1,025.90 seconds (17:05). The 302 warnings are existing
  hmmlearn/NumPy deprecation warnings. Slow tests remain opt-in under the
  repository's default test configuration.

Local evidence is under `/private/tmp/macs813-verification/`: baseline
`build.log`, `results.json`, `assert-results.json`, `edge-results.json`,
and `baseline-tests.log`; fixed `build-fixed.log`, `targeted-inplace.log`,
and `assert-api-tests.log`; wheel `build-wheel.log`, `install-wheel.log`,
`targeted-wheel.log`, and `validate_installation.py`; full-suite
`full-suite.log`; plus `after/`
containing the repeated CLI probes and assertion-enabled build.
Temporary BAMs/logs are not committed fixtures.
