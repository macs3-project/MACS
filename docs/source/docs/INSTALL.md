# Install MACS3

MACS3 3.0.5 requires **Python 3.12 or later**. The project tests Python
3.12, 3.13, and 3.14 on Linux and macOS. Windows is not currently supported.

We recommend installing MACS3 in a dedicated virtual environment. The
following commands use `python3`; check that it refers to a supported Python
version on your system:

```bash
python3 --version
python3 -m venv macs3env
source macs3env/bin/activate
```

If you use Conda, activate a Python 3.12-or-later environment instead of
creating a `venv`.

## Install from PyPI

```bash
python -m pip install --upgrade pip
python -m pip install macs3
macs3 --version
```

To upgrade an existing installation, run
`python -m pip install --upgrade macs3` in the environment where MACS3 is
installed.

The PyPI release is a source distribution. Installing it may compile MACS3's
Cython extensions, so you need a working C compiler and Python development
headers. On macOS, install the Xcode Command Line Tools. On Linux, install the
compiler toolchain and Python headers provided by your distribution.

## Install from Bioconda

In a Conda environment configured with the `conda-forge` and `bioconda`
channels, install MACS3 with:

```bash
conda install macs3
macs3 --version
```

Bioconda packages may become available after the PyPI release. Check the
[Bioconda MACS3 package](https://bioconda.github.io/recipes/macs3/README.html)
to confirm that version 3.0.5 is available before using this method for the
new release; use PyPI if you need 3.0.5 sooner.

## Install from source

Clone the repository with its submodules, activate a supported Python
environment, and install from the checkout:

```bash
git clone --recurse-submodules https://github.com/macs3-project/MACS.git
cd MACS
python -m pip install .
macs3 --version
```

If you already have a checkout, run `git submodule update --init --recursive`
before installing. The same `python -m pip install .` command works after
unpacking a release source archive.

To build a wheel locally instead, install the build frontend and run:

```bash
python -m pip install build
python -m build --wheel
```

The wheel is written to `dist/`; install it with `python -m pip install`
followed by its filename.

## Dependencies and build requirements

The package installer resolves the runtime dependencies declared in
`pyproject.toml`: NumPy >=1.25, SciPy >=1.12, pandas >=2.0, AnnData >=0.10,
hmmlearn >=0.3.2, and scikit-learn >=1.3. Source builds use Cython >=3.0,
NumPy >=1.25, and setuptools >=80.0. MACS3's Cython extension sources use
the `.py` suffix and are compiled during the build. The `cykhash` package is
**not** a dependency of MACS3 3.0.5.

The `callvar` extension also uses bundled fermi-lite code and the SIMDe
submodule. If you install from a Git checkout, initialize that submodule as
shown above.
