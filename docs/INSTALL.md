# INSTALL Guide For MACS3
Time-stamp: <2023-06-16 11:46:32 Tao Liu>

Please check the following instructions to complete your installation.

## Prerequisites

Here we list some prerequisites for installing and running MACS3. But
if you are using conda or pip to install, the installer will check the
dependencies and install them if necessary.

### Python3

MACS v3.x.x requires Python3. We have tested MACS in Python3.8 to 3.11. 

### NumPy, hmmlearn

MACS requires NumPy>=1.19 (>=1.24 recommended) and hmmlearn>=0.3 during installation. Note that hmmlearn further requires SciPy and sklearn (aka scikit-learn). 

### Cython

[Cython](http://cython.org/) is required to translate .pyx codes to .c
code. The version of Cython has to be >=0.29.

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
header files are needed. If you are using Mac OSX, I recommend you 
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

`$ python3 -m venv MyPythonEnv/`

Then activate it by

`$ source MyPythonEnv/bin/activate`

If you use 'conda', it will also provide virtual environment. Please 
read: [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Install through PyPI

The easiest way to install MACS is through PyPI system. Get `pip` if
it's not available in your system. If you create a virtual environment
as described before, your `pip` command will install everything under
the folder you specified previously through `python3 -m env` command,
or to your active conda environment. 

Then under the command line, type `pip install macs3`. PyPI will
install dependencies automatically if it is absent.

To upgrade MACS3, type `pip install --upgrade macs3`. It will check
currently installed MACS3, compare the version with the one on PyPI
repository, download and install a newer version while necessary.

## Install from source

MACS uses `pip` for source code installations. To install a source 
distribution of MACS, unpack the distribution tarball, or clone Git 
repository with `git clone --recurse-submodules git@github.com:taoliu/MACS.git`. 
Go to the directory where you cloned MACS from github, and simply
run the install command:

 `$ pip install .`

By default, the script will install python library and executable
codes according to the environment. When you run the command under
virtualenv or conda environment, the script will install to the virtual
environment instead. When you run the command without virtual environment, 
you may need to be root or administrator of the machine so as to 
complete the installation. Please contact the system administrator
if you want their help. 

## Configure environment variables

*Note*, if you are using a virtual environment, you should skip this
section since all the corresponding environment variables have been
correctly set while you `activate` the environment.

After running the setup script, you might need to add the install
location to your `PYTHONPATH` and `PATH` environment variables. The
process for doing this varies on each platform, but the general
concept is the same across platforms.

### PYTHONPATH

To set up your `PYTHONPATH` environment variable, you'll need to add
the value `PREFIX/lib/pythonX.Y/site-packages` to your existing
`PYTHONPATH`. In this value, X.Y stands for the major–minor version of
Python you are using (such as 3.7; you can find this with
`sys.version[:3]` from a Python command line). `PREFIX` is the install
prefix where you installed MACS. If you did not specify a prefix on
the command line, MACS will be installed using Python's sys.prefix
value.

On Linux, using bash, I include the new value in my `PYTHONPATH` by
adding this line to my `~/.bashrc`::

 `$ export
 PYTHONPATH=/home/taoliu/lib/python3.7/site-packages:$PYTHONPATH`

Using Windows, you need to open up the system properties dialog and
locate the tab labeled Environment. Add your value to the `PYTHONPATH`
variable, or create a new `PYTHONPATH` variable if there isn't one
already.

### PATH

Just like your `PYTHONPATH`, you'll also need to add a new value to
your PATH environment variable so that you can use the MACS command
line directly. Unlike the `PYTHONPATH` value, however, this time
you'll need to add `PREFIX/bin` to your PATH environment variable. The
process for updating this is the same as described above for the
`PYTHONPATH` variable::

 `$ export PATH=/home/myaccount/bin:$PATH`

