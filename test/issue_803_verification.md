# Pileup AUTO format: verification and fix

Tickets: [#803](https://github.com/macs3-project/MACS/issues/803),
part of [#757](https://github.com/macs3-project/MACS/issues/757).

## Verified baseline

Fetched current `main` and created `fix/macs/803-pileup-auto-format` at
`d53d85a7f97e0a546c562c0e718612a391d28094`.
Built the revision in place with Python 3.14.6, Cython 3.3.0, NumPy 2.3.5:

```sh
python setup.py build_ext --inplace -j 4
```

The probe imported `MACS3/IO/Parser.cpython-314-darwin.so`, used that same
package and interpreter in CLI subprocesses, ran from a temporary directory,
disabled core dumps, and bounded each invocation with a 20-second timeout.

Original input:

```text
chr1\t100\t150\tr\t0\t+
chr1\t200\t250\tr\t0\t-
```

Here `\t` denotes a literal tab. Commands:

```sh
macs3 pileup -i reads.bed -o default.bdg
macs3 pileup -i reads.bed -f AUTO -o auto.bdg
macs3 pileup -i reads.bed -f BED -o bed.bdg
```

| Input | Default format | Explicit AUTO | Explicit BED |
| --- | --- | --- | --- |
| Original two-read BED | Exit 1, rejects AUTO | Exit 1, rejects AUTO | Exit 0, expected coverage |
| Gzipped original BED | Exit 1, rejects AUTO | Exit 1, rejects AUTO | Exit 0, identical coverage |
| Empty BED | Exit 1, rejects AUTO | Exit 1, rejects AUTO | Exit 0, empty bedGraph |
| `not an alignment` followed by newline | Exit 1, rejects AUTO | Exit 1, rejects AUTO | Exit 1, IndexError in BED parsing |

No probe crashed or timed out. Successful original-input coverage:

```text
chr1\t0\t50\t0.00000
chr1\t50\t100\t1.00000
chr1\t100\t250\t2.00000
chr1\t250\t300\t1.00000
```

Local probe script, logs, and full results are under
`/private/tmp/macs803-verification/`.

## Cause and scope

`opt_validate_pileup()` has no AUTO branch, despite the CLI default and
documentation promising detection. Other command validators already assign
`guess_parser` for AUTO. The pileup loaders pass the buffer size and execution
policy to the selected parser, so they can use the existing detector.

The actual CLI dispatches to `pileup_v2_cmd.run()`, which shares the validator
and loaders with `pileup_cmd.run()`. Existing pileup tests check the AUTO
argument default but do not exercise successful AUTO detection.

## Implementation plan

1. Add the missing AUTO-to-`guess_parser` branch in `opt_validate_pileup()`,
   following neighboring validators. Preserve explicit format handling and
   keep BAMPE, BEDPE, and FRAG explicitly selected, as documented.
2. Add CLI regression tests using `run_macs3` and its parent-installation
   environment. Cover omitted `-f` and explicit AUTO, compare bedGraph bytes
   with explicit BED and hand-derived expected coverage, and include gzip,
   multiple input files, extension size, and bidirectional pileup.
3. Exercise other detectable single-end formats with small existing fixture
   writers and compare AUTO output with their explicit-format output. Include
   SAM and short BAM files. Exercise the public direct command entry point as
   well as the actual CLI dispatch.
4. Cover empty and malformed AUTO inputs in timeout-bounded subprocesses.
   Require a nonzero exit without a native crash or hang when detection is
   impossible. Preserve explicit BED's existing empty-file behavior. Avoid
   changing shared detector heuristics or unrelated parser behavior unless
   verification demonstrates a defect necessary to resolve this ticket.
5. Run the focused regressions and existing pileup/parser tests on the
   in-place compiled build. Build and install a wheel, then run changed
   subprocess tests from a clean source checkout with no in-place extensions,
   verifying that parent and children use the installed package.
6. Extend this evidence with fixed results and test commands. Any fix PR
   should link both #803 and #757.

## Implemented fix

Added the missing AUTO branch in `opt_validate_pileup()`. It assigns the
existing public `guess_parser`, matching the other validators. No shared
parser or pileup algorithm changes are included.

Added 24 public command/CLI regression cases to `test/test_pileup_cmd.py`:
default and explicit AUTO for the original BED (plain and gzip), five working
detectable single-end formats, one- and two-record BAMs, mixed BED/SAM inputs,
custom extension size, bidirectional coverage, empty/malformed inputs, and
the direct public command entry point. CLI subprocesses use `run_macs3`,
which selects the parent's MACS3 installation and interpreter, runs from
`tmp_path`, disables core dumps, and applies a 20-second timeout.

Re-running all 12 original probes after the fix gives:

| Input | Default format | Explicit AUTO | Explicit BED |
| --- | --- | --- | --- |
| Original two-read BED | Exit 0, expected coverage | Exit 0, expected coverage | Exit 0, same coverage |
| Gzipped original BED | Exit 0, expected coverage | Exit 0, expected coverage | Exit 0, same coverage |
| Empty BED | Exit 1, cannot detect format | Exit 1, cannot detect format | Exit 0, empty bedGraph |
| Malformed text | Exit 1, parser exception | Exit 1, parser exception | Exit 1, parser exception |

No probe crashed or timed out. Full baseline and fixed probe results are in
`baseline-results.json` and `fixed-results.json` under the local verification
directory.

## Additional existing parser limitations

The exploratory seven-format run exposed issues outside the missing pileup
validator branch. These are reproduced independently through the unchanged
`guess_parser()`/parser APIs, without using the pileup validator:

- ELANDMULTI input generated by the existing `elandmulti_lines()` fixture is
  misdetected as ELAND by the shared detector and produces an empty track.
  Explicit `-f ELANDMULTI` instead raises `TypeError: an integer is required`
  when converting byte-valued hit counts. It cannot serve as a successful
  explicit-output reference for this ticket.
- BOWTIE input succeeds with explicit `-f BOWTIE`, but shared AUTO detection
  raises `IndexError` while the SAM candidate probes its shorter rows.

Both formats were investigated with bounded subprocesses; evidence is in
`legacy-formats.log`. They are excluded from the successful-output test
matrix rather than marked as skipped or expected failures. This fix connects
pileup to the existing detector; it does not claim to repair those detector
and legacy parser limitations.

## Validation

The initial existing pileup suite passed **87 tests, with 10 skips** before
the fix. Installing pysam 0.24.1 into a temporary dependency directory enabled
the SAM/BAM fixtures that were otherwise skipped.

The final in-place command/validator/pileup run passed **669 tests**:

```sh
PYTHONPATH=/private/tmp/macs803-verification/deps python -m pytest \
  test/test_pileup_cmd.py test/test_OptValidator.py \
  test/test_Pileup.py test/test_PileupV2.py -q
```

All commands use the development `.venv/bin/python` backed by the `macs3dev`
interpreter. The native baseline extensions remain unchanged because the
production fix is in the Python validator.

The parser, import-order, EOF, and BAM truncation regression run passed
**359 tests**:

```sh
PYTHONPATH=/private/tmp/macs803-verification/deps python -m pytest \
  test/test_Parser.py test/test_Parser_imports.py \
  test/test_Parser_frag_eof.py test/test_Parser_bam_truncation.py -q
```

The wheel build succeeded:
`macs3-3.1.0-cp314-cp314-macosx_11_0_x86_64.whl`. It was built with
`python -m build --wheel --no-isolation` from a clean `git archive` of the
baseline plus the changed validator. The pinned SIMDe submodule was populated
from commit `8cd136a43bae7ab9b82316179b9cef8887726778`.

Installed that wheel with `pip install --no-deps --no-index --target` into
`/private/tmp/macs803-verification/installed-wheel`. Ran the same four-module
command/validator/pileup pytest selection from the clean source directory,
with `python -P -m pytest`, `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`, and
`PYTHONPATH` containing only the installed-wheel root and temporary pysam
dependencies. The checkout had no in-place `.so` files. Before pytest, import
assertions confirmed both `Parser` and `OptValidator` came from the wheel
installation. All **669 tests passed**, including all changed subprocess
tests, with no skips. The parent-derived `macs3_env` selected that installed
package for CLI children.

Final logs: `pytest-inplace-final.log` (669 passed),
`pytest-parser-final.log` (359 passed), and `pytest-wheel.log` (669 passed),
all under `/private/tmp/macs803-verification/`. `git diff --check` also passed.
These are targeted validation runs; the entire repository suite was not run.
