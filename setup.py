#!/usr/bin/env python3
# Time-stamp: <2019-12-12 13:04:02 taoliu>

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
install_requires = [f"numpy{numpy_requires}",]

def main():
    # if float(sys.version[:3])<3.6:
    #     sys.stderr.write("CRITICAL: Python version must >= 3.6!\n")
    #     sys.exit(1)

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

    ext_modules = [Extension("MACS2.Prob", ["MACS2/Prob.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                   Extension("MACS2.IO.Parser",["MACS2/IO/Parser.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("MACS2.Pileup", ["MACS2/Pileup.pyx","MACS2/cPosValCalculation.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                   Extension("MACS2.PeakModel", ["MACS2/PeakModel.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("MACS2.PeakDetect", ["MACS2/PeakDetect.pyx"], extra_compile_args=extra_c_args),
                   Extension("MACS2.Signal", ["MACS2/Signal.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("MACS2.IO.PeakIO", ["MACS2/IO/PeakIO.pyx"], extra_compile_args=extra_c_args),
                   Extension("MACS2.IO.BedGraphIO", ["MACS2/IO/BedGraphIO.pyx"], extra_compile_args=extra_c_args),                   
                   Extension("MACS2.IO.FixWidthTrack", ["MACS2/IO/FixWidthTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("MACS2.IO.PairedEndTrack", ["MACS2/IO/PairedEndTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   Extension("MACS2.IO.BedGraph", ["MACS2/IO/BedGraph.pyx"], libraries=["m"], extra_compile_args=extra_c_args),
                   Extension("MACS2.IO.ScoreTrack", ["MACS2/IO/ScoreTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                   Extension("MACS2.IO.CallPeakUnit", ["MACS2/IO/CallPeakUnit.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                   #Extension("MACS2.Statistics", ["MACS2/Statistics.pyx"], libraries=["m"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=extra_c_args),
    ]

    with open("README.md", "r") as fh:
        long_description = fh.read()
        
    setup(name="MACS2",
          version="2.2.7.2",
          description="Model Based Analysis for ChIP-Seq data",
          long_description = long_description,
          long_description_content_type="text/markdown",
          author='Tao Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://github.com/taoliu/MACS/',
          package_dir={'MACS2' : 'MACS2'},
          packages=['MACS2', 'MACS2.IO'],
          package_data={'MACS2':['*.pxd']},
          scripts=['bin/macs2', ],
          classifiers=[
              'Development Status :: 5 - Production/Stable',
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
