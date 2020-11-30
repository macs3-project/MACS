#!/usr/bin/env python3
"""Description:

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import os
from setuptools import setup, Extension
from distutils.version import LooseVersion
import subprocess

numpy_requires = '>=1.17'
cykhash_requires = '>=1.0.2'
install_requires = [f"numpy>={numpy_requires}",f"cykhash>={cykhash_requires}"]

def main():
    if float(sys.version[:3])<3.6:
        sys.stderr.write("CRITICAL: Python version must >= 3.6!\n")
        sys.exit(1)

    cwd = os.path.abspath(os.path.dirname(__file__))

    # install required numpy
    p = subprocess.call([sys.executable, "-m", 'pip', 'install', f'numpy{numpy_requires}'],cwd=cwd)
    if p != 0:
        # Could be due to a too old pip version and build isolation, check that
        try:
            # Note, pip may not be installed or not have been used
            import pip
            if LooseVersion(pip.__version__) < LooseVersion('18.0.0'):
                raise RuntimeError("Installing requirements failed. Possibly due "
                                   "to `pip` being too old, found version {}, "
                                   "needed is >= 18.0.0.".format(pip.__version__))
            else:
                raise RuntimeError("Installing requirements failed!")
        except ImportError:
            raise RuntimeError("Installing requirement failed! `pip` has to be installed!")

    from numpy import get_include as numpy_get_include
    numpy_include_dir = [numpy_get_include()]

    # I intend to use -Ofast, however if gcc version < 4.6, this option is unavailable so...
    extra_c_args = ["-w","-O3","-ffast-math","-g0"] # for C, -Ofast implies -O3 and -ffast-math

    ext_modules = [ Extension("MACS3.Signal.Prob", ["MACS3/Signal/Prob.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.Pileup", ["MACS3/Signal/Pileup.pyx","MACS3/Signal/cPosValCalculation.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.PeakModel", ["MACS3/Signal/PeakModel.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.PeakDetect", ["MACS3/Signal/PeakDetect.pyx"], extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.SignalProcessing", ["MACS3/Signal/SignalProcessing.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.FixWidthTrack", ["MACS3/Signal/FixWidthTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.PairedEndTrack", ["MACS3/Signal/PairedEndTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.BedGraph", ["MACS3/Signal/BedGraph.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.ScoreTrack", ["MACS3/Signal/ScoreTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.CallPeakUnit", ["MACS3/Signal/CallPeakUnit.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.IO.Parser",["MACS3/IO/Parser.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.IO.PeakIO", ["MACS3/IO/PeakIO.pyx"], extra_compile_args=extra_c_args),
                    Extension("MACS3.IO.BedGraphIO", ["MACS3/IO/BedGraphIO.pyx"], extra_compile_args=extra_c_args), ]

    with open("README.md", "r") as fh:
        long_description = fh.read()

    setup(name="MACS3",
          version="3.0.0a1",
          description="Model Based Analysis for ChIP-Seq data",
          long_description = long_description,
          long_description_content_type="text/markdown",
          author='Tao Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://github.com/taoliu/MACS/',
          package_dir={'MACS3' : 'MACS3'},
          packages=['MACS3', 'MACS3.IO', 'MACS3.Signal', 'MACS3.Commands','MACS3.Utilities'],
          package_data={'MACS3':['*.pxd']},
          scripts=['bin/macs3', ],
          classifiers=[
              'Development Status :: 3 - Alpha',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3.8',
              'Programming Language :: Cython',
              ],
          install_requires=install_requires,
          setup_requires=install_requires,
          python_requires='>=3.6',
          ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()
