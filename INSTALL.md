# INSTALL Guide For MACS
Time-stamp: <2019-09-30 14:47:45 taoliu>

Please check the following instructions to complete your installation.

## Prerequisites

Python3

[Numpy](http://www.scipy.org/Download) (>=1.15)

GCC is required to compile `.c` codes in MACS v2 package, and python
header files are needed. If you are using Mac OSX, I recommend you
install Xcode; if you are using Linux, you need to make sure
`python-dev` is installed.

[Cython](http://cython.org/) (>=0.29) is required to generate `.c`
files from `.pyx` files using `setup.py` script.

## Prepare a virtual Python environment 

I strongly recommend to install your MACS2 program in a virtual
environment, so that you have full control of your installation and
won't mess up with your system libraries. To learn about virtual
environment, read
[this article](https://docs.python.org/3/library/venv.html). A simple
way to create a virtual environment of Python3 is

`$ python3 -m venv MyPythonEnv/`

Then active it by

`$ source MyPythonEnv/bin/active`

## Install through PyPI

The easiest way to install MACS2 is through PyPI system. Get `pip` if
it's not available in your system. If you create a virtual environment
as described before, your `pip` command will install everything under
the folder you specified previously through `python3 -m env` command.

Then under the command line, type `pip install MACS2`. PyPI will
install Numpy automatically if it is absent. 

To upgrade MACS2, type `pip install --upgrade MACS2`. It will check
currently installed MACS2, compare the version with the one on PyPI
repository, download and install newer version while necessary.

## Install from source

MACS uses Python's distutils tools for source installations. To
install a source distribution of MACS, unpack the distribution tarball
and open up a command terminal. Go to the directory where you unpacked
MACS, and simply run the install script:

 `$ python setup.py install`

By default, the script will install python library and executable
codes globally, which means you need to be root or administrator of
the machine so as to complete the installation. Please contact the
administrator of that machine if you want their help. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script. Use the –help option to see a brief list of available options:

 `$ python setup.py --help`

For example, if I want to install everything under my own HOME
directory, use this command:

 `$ python setup.py install --prefix /home/taoliu/`


## Configure enviroment variables

*Note*, if you are using virtual environment, you should skip this
section since all the corresponding environment variables have been
correctly set while you `activate` the environment.

After running the setup script, you might need to add the install
location to your `PYTHONPATH` and `PATH` environment variables. The
process for doing this varies on each platform, but the general
concept is the same across platforms.

### PYTHONPATH

To set up your `PYTHONPATH` environment variable, you'll need to add the
value `PREFIX/lib/pythonX.Y/site-packages` to your existing
`PYTHONPATH`. In this value, X.Y stands for the major–minor version of
Python you are using (such as 2.7 ; you can find this with
`sys.version[:3]` from a Python command line). `PREFIX` is the install
prefix where you installed MACS. If you did not specify a prefix on
the command line, MACS will be installed using Python's sys.prefix
value.

On Linux, using bash, I include the new value in my `PYTHONPATH` by
adding this line to my `~/.bashrc`::

 `$ export PYTHONPATH=/home/taoliu/lib/python2.7/site-packages:$PYTHONPATH`

Using Windows, you need to open up the system properties dialog, and
locate the tab labeled Environment. Add your value to the `PYTHONPATH`
variable, or create a new `PYTHONPATH` variable if there isn't one
already.

### PATH

Just like your `PYTHONPATH`, you'll also need to add a new value to your
PATH environment variable so that you can use the MACS command line
directly. Unlike the `PYTHONPATH` value, however, this time you'll need
to add `PREFIX/bin` to your PATH environment variable. The process for
updating this is the same as described above for the `PYTHONPATH`
variable::

 `$ export PATH=/home/taoliu/bin:$PATH`

--
Tao Liu <vladimir.liu@gmail.com>

