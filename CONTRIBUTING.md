# Contributing to MACS3

Thank you for helping improve MACS3. Bug reports, documentation improvements,
tests, code, and ideas are all welcome. These guidelines explain how to get
started and help maintainers review your contribution.

## Community standards

Please follow the project's [Code of Conduct](https://github.com/macs3-project/MACS/blob/main/CODE_OF_CONDUCT.md)
in issues, discussions, pull requests, and other project spaces. To report
unacceptable behavior, use the private contact listed in that document rather
than posting the report publicly.

## Ask a question or propose an idea

Check the [MACS3 documentation](https://macs3-project.github.io/MACS/) and
search existing [discussions](https://github.com/macs3-project/MACS/discussions)
and [issues](https://github.com/macs3-project/MACS/issues) first. Use
Discussions for usage questions, analysis advice, and early-stage ideas. Open
an issue for a reproducible bug or a concrete feature request. The repository's
issue chooser provides templates for both.

When reporting a bug, include:

- The exact command and input format, plus what you expected and what happened.
- The output of `macs3 --version`, your Python version, operating system, and
  CPU architecture.
- Relevant error messages or logs, and the smallest shareable input that
  reproduces the problem. Remove private or sensitive data before posting.

For a feature request, describe the use case, the behavior you would like,
and any alternatives you considered. If the design is still open, start a
Discussion before writing code.

## Set up a development environment

MACS3 requires Python 3.12 or later. On Linux or macOS, clone the repository
with its submodules and install the development requirements in an isolated
environment:

```bash
git clone --recurse-submodules https://github.com/macs3-project/MACS.git
cd MACS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .
```

If you already use Conda, activate your development environment before the
installation commands instead of creating `.venv`. The sources are organized
under `MACS3/`, tests under `test/`, and Sphinx documentation under
`docs/source/`.

## Make and check a change

Keep changes focused and follow the surrounding code's conventions. Use PEP 8
for Python, keep command-line help and documentation in sync with CLI changes,
and add or update tests for changed behavior. Keep new test fixtures small;
describe any data-heavy tests and their expected runtime.

Run the relevant tests before opening a pull request:

```bash
python -m pytest test
# Or target a module while developing:
python -m pytest test/test_PeakModel.py -k peak
```

For packaging changes, also check that a wheel builds:

```bash
python -m pip install build
python -m build --wheel
```

For documentation changes, install the documentation tools and build the site:

```bash
python -m pip install sphinx myst-parser sphinx-rtd-theme sphinx-autodoc-typehints
make -C docs html
```

## Open a pull request

Create a branch from `main` and open a focused pull request against `main`.
Explain the problem and your solution, link any related issue, and list the
tests and build checks you ran. Call out new CLI options or compatibility
changes. For user-facing behavior, include a short before-and-after command
example or screenshot when it helps reviewers see the effect.

Check the CI results after opening the pull request. If a check fails for a
reason you believe is unrelated to your change, explain that in a comment;
maintainers can then investigate or rerun it. Reviewers may request further
changes before merging.
