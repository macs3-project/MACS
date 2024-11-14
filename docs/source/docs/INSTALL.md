# INSTALL Guide For MACS3

Please check the following instructions to complete your installation.

## Prerequisites

Here we list some prerequisites for installing and running MACS3. But
if you are using conda or pip to install, the installer will check the
dependencies and install them if necessary. Therefore, this section is
for reference purpose, and if you are just looking for steps to
install MACS3, please go to the next section. Please note that, we
haven't tested installation on any Windows OS, so currently only Linux
and Mac OS systems are supported.

### Python3

MACS v3 requires Python3. We have tested MACS in Python3.9 to 3.12. 

### NumPy, hmmlearn, Scipy, scikit-learn

MACS calls functions from [NumPy](https://numpy.org/) and
[hmmlearn](https://hmmlearn.readthedocs.io/). Since hmmlearn further
requires [scikit-learn](https://scikit-learn.org/) which requires
[SciPy](https://scipy.org/), and these libraries are crucial for
reproducing your results, we also add them into the requirement list
with specific version numbers. So here is the list of the required
python libraries that will impact the numerical calculation in MACS3:

 - numpy>=1.25 <2.0.0
 - hmmlearn>=0.3.2
 - scikit-learn>=1.3
 - scipy>=1.12

### Cython

[Cython](http://cython.org/) is required to translate .pyx codes to .c
code. The version of Cython has to be >=3.0. 

### cykhash

[cykhash](https://github.com/realead/cykhash) is a fast and efficient
hash implementation in Cython. It can be seen as the cython
implementation of
[khash](https://github.com/attractivechaos/klib/blob/master/khash.h). It
is used to replace python dictionary in MACS3 codes. We require
cykhash version 2.

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

### GCC, Python-dev, meson ... 

GCC is required to compile `.c` codes in MACS v3 package, and python
header files are needed. If you are using Mac OSX, we recommend you
install Xcode; if you are using Linux, you need to make sure
`python-dev` package is installed -- the actual package name depends
on the Linux OS distribution, you are using. Also, since the most
recent Numpy/Scipy use [meson](https://mesonbuild.com/) to build the
libraries, make sure they have been installed.

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

If you use 'conda' through 'miniforge' project, it will also provide
virtual environment. Please read:
[miniforge](https://github.com/conda-forge/miniforge). For example,
after installing 'conda', you can use `conda create -n MACS3` to
create a new environment called 'MACS3' then switch to this
environment with `conda activate MACS3`.

There is another solution, [pipx](https://pipx.pypa.io/), to make a
clean isolated python environment to run python tools such as
MACS3. We won't explore it here but if you are insterested in it,
please click the link above and read the tutorial.

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

## Install through conda

To install MACS3 using 'conda', simply execute `conda install -c
bioconda macs3` in your conda environment. This command installs MACS3
along with its dependencies from the Bioconda channel. Please ensure
conda is installed and a dedicated conda environment has been created
and activated beforehand for a smooth installation process.

## Install from source through pip

MACS uses `pip` for source code installations. To install a source
distribution of MACS, unpack the distribution tarball, or clone Git
repository with `git clone --recurse-submodules 
https://github.com/macs3-project/MACS.git`.  Go to the directory where you
cloned MACS from github or unpacked from tarball, and simply run the
install command:

 `$ pip install .`

The command will treat the current working directory as the 'package'
to be installed. The behavior will be the same as `pip install macs3`
as described in the previous section "Install through PyPI".

You can also install MACS3 from source code in a "modern" two-steps
way. First, use the build system to make a wheel (in this case, you
need to install `build` first by `pip install build`):

`$ python -m build --wheel`

This will make a '.whl' file under 'dist' directory. Then you can
install the wheel through `pip`:

`$ pip install dist/MACS3-3.x.x-x-x-x.whl`

Please use the real filename in the 'dist' directory.

