# Native compatibility regression checks

Run builds and tests with the `macs3dev` environment active. The native
compatibility runner checks which compiled modules pytest imports before it
collects the public API and CLI regressions:

```sh
uv run python setup.py build_ext --inplace
uv run python scripts/run_native_regressions.py --run-slow
```

The default suite includes exact PETrackII ordering against the installed
NumPy's structured `argsort`, barcode association/subsets, repeated cleanup,
compressed framing errors, mixed-member fallbacks, fresh parser/stream import
order, repeated stream failure/recovery, and BAM/BAMPE/FRAG CLI
artifacts with multiple inputs, CPU budgets and backend policies. Existing
execution, model, parser-boundary, scheduling and lifecycle regressions are
included by the runner. Tests call public APIs and run native failure probes
with timeouts, process-group cleanup, temporary scripts and temporary working
directories. Subprocess imports follow the parent pytest installation.

`--run-slow` additionally exercises corruption beyond the Python header
reader's 10 MB lookahead, CLI inputs above the process threshold, and allocator
growth. Most individual large cases take about 1–5 seconds and use under
150 MiB; read-ahead stream cases additionally generate 5 MiB of seeded random
data and records spanning multiple buffers. Fixtures are generated in pytest's
temporary directory. Linux tests
select available/disabled/missing/invalid huge-page discovery using a small
`LD_PRELOAD` shim that redirects only two sysfs reads. Real track allocation,
mapping, growth, sorting and release remain intact. NumPy's handler diagnostics
verify which branch ran. This checks allocator behavior, not physical huge-page
residency. macOS verifies default-allocator fallback. PETrackII uses NumPy's
default allocator on both platforms; the custom handler applies to PETrackI.

The partial process-start regression fails the second `Process.start` after
one real child starts. It checks child joining and temporary-file cleanup.
The host-eligibility case runs on Linux; a separate controlled eligibility
case can exercise real fork cleanup on other POSIX hosts. Production macOS
safety/fallback assertions remain in the backend suite.

For distribution validation, build an sdist, rebuild a wheel from the extracted
sdist, then verify both artifacts from the matching checkout:

```sh
python scripts/check_native_distribution.py /path/to/macs3.tar.gz /path/to/macs3.whl
```

The checker verifies native build inputs and licenses against source bytes,
compiled wheel modules, absence of test fixtures/in-place binaries from the
sdist, and third-party notices in the wheel's license directory. Install the
wheel and run pytest from a separate clean checkout with no in-place extensions:

```sh
python scripts/run_native_regressions.py --installed --run-slow
```

The runner rejects source imports in installed mode. Running a probe only from
an extracted wheel does not validate pytest's subprocess launch behavior.

`.github/workflows/native-compatibility.yml` performs both installation modes
with Python 3.12/3.13/3.14, NumPy 1.26 where supported, selected older NumPy 2.x
versions, and current NumPy 2.x. Linux and macOS x86_64/arm64 jobs report the
resolved versions. Build-isolation constraints pin the same NumPy/Cython used
at runtime. The existing full-suite workflows continue to run separately.
Other architecture workflows retain their manual status. Configured matrix
jobs establish compatibility only once their builds and tests pass.

Version selection follows NumPy's [Python 3.12 support in 1.26](https://numpy.org/doc/2.1/release/1.26.0-notes.html),
[Python 3.13 support in 2.1](https://numpy.org/doc/2.3/release/2.1.0-notes.html),
and [3.14 wheel availability in 2.3.3](https://pypi.org/project/numpy/2.3.3/).
