# INSTALL Guide For MACS3
Time-stamp: <2024-02-15 12:07:37 Tao Liu>

Please check the following instructions to complete your installation.

## Prerequisites

Here we list some prerequisites for installing and running MACS3. But
if you are using conda or pip to install, the installer will check the
dependencies and install them if necessary. Please note that, we
haven't tested installation on any Windows OS, so currently only Linux
and Mac OS systems are supported.

### Python3

MACS v3 requires Python3. We have tested MACS in Python3.9 to 3.12. 

### NumPy, hmmlearn

MACS requires [NumPy](https://numpy.org/)>=1.19 (>=1.24 recommended)
and [hmmlearn](https://hmmlearn.readthedocs.io/)>=0.3 during
installation. Note that hmmlearn further requires
[SciPy](https://scipy.org/) and
[scikit-learn](https://scikit-learn.org/).

### Cython

[Cython](http://cython.org/) is required to translate .pyx codes to .c
code. The version of Cython has to be >=0.29. We recommend Cython
version >= 3.

### cykhash

[cykhash](https://github.com/realead/cykhash) is a fast and efficient
hash implementation in Cython. It is used to replace python dictionary
in MACS3 codes. Since it requires Cython, make sure you install Cython
first, then install cykhash. 

### fermi-lite and simde

A newly added `callvar` subcommand in MACS3 uses
[fermi-lite](https://github.com/lh3/fermi-lite) to assemble the DNA
sequence in a peak region while necessary. A modified fermi-lite has
been included in MACS3 package. Since fermi-lite was implemented using
intel SSE2 intrinsics for x86 CPUs, we added
[simde](https://github.com/simd-everywhere/simde) as submodule to
solve the compatibility issues on non-x86 architectures. Note that, we
may remove this submodule and add simde in *dependencies* of MACS3
later.

### GCC and Python-dev 

GCC is required to compile `.c` codes in MACS v3 package, and python 
header files are needed. If you are using Mac OSX, we recommend you 
install Xcode; if you are using Linux, you need to make sure 
`python-dev` package is installed -- the actual package name depends 
on the Linux OS distribution, you are using. 

## Prepare a virtual Python environment 

We strongly recommend installing your MACS program in a virtual
environment, so that you have full control of your installation and
won't mess up with your system libraries. To learn about virtual
environment, read [this
article](https://docs.python.org/3/library/venv.html). A simple way to
create a virtual environment of Python3 is

`$ python3 -m venv MACS3env/`

Then activate it by

`$ source MACS3env/bin/activate`

You can also let the virtual enviroment access system-wide libraries,
the Python libraries globally installed in your operating system, so
that you may not need to install certain dependencies later. 

`$ python3 -m venv --system-site-packages`

If you use 'conda', it will also provide virtual environment. Please
read:
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)
or [miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Install through PyPI

The easiest way to install MACS is through PyPI system. Get `pip` if
it's not available in your system. If you create a virtual environment
as described before, your `pip` command will install everything under
the folder you specified previously through `python3 -m env` command,
or to your active conda environment. 

Then under the command line, type `pip install macs3`. PyPI will
install dependencies automatically if it is absent. By default, `pip`
will install the newest version of dependencies that satisfy the
requirements of MACS3. When you run the command without virtual
environment, you may need to be the root user or system administrator
so as to complete the installation. Please contact the system
administrator if you want their help.

To upgrade MACS3, type `pip install --upgrade macs3`. It will check
currently installed MACS3, compare the version with the one on PyPI
repository, download and install a newer version while necessary.

## Install from source through pip

MACS uses `pip` for source code installations. To install a source
distribution of MACS, unpack the distribution tarball, or clone Git
repository with `git clone --recurse-submodules
git@github.com:macs3-project/MACS.git`.  Go to the directory where you
cloned MACS from github, and simply run the install command:

 `$ pip install .`

The command will treat the current working directory as the 'package'
to be installed. The behavior will be the same as `pip install macs3`
as described in the previous section "Install through PyPI".

