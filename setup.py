#!/usr/bin/env python
# Time-stamp: <2019-09-20 14:04:51 taoliu>

"""Description: 

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

Use this when you need Cython regenerate .c files.

Copyright (c) 2008-2019 Tao Liu <vladimir.liu@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@version: $Revision$
@author:  Tao Liu
@contact: vladimir.liu@gmail.com
"""

import os
import sys
from setuptools import setup, Extension

# Use build_ext from Cython if found
command_classes = {}
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
    from Cython.Build import cythonize
    has_cython = True
    # Note: if this script is run under pip, cython should be installed already due to requirements.txt setting.
except:
    has_cython = False

try: 
    from numpy import get_include as numpy_get_include 
    numpy_include_dir = [numpy_get_include()] 
except: 
    numpy_include_dir = [] 


def main():
    if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

    # I intend to use -Ofast, however if gcc version < 4.6, this option is unavailable so...
    extra_c_args = ["-w","-O3","-ffast-math"] # for C, -Ofast implies -O3 and -ffast-math

    if has_cython:
        ext_modules = [Extension("MACS2.Prob", ["MACS2/Prob.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                       Extension("MACS2.IO.Parser",["MACS2/IO/Parser.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.Pileup", ["MACS2/Pileup.pyx","MACS2/cPosValCalculation.pxd","MACS2/cPosValCalculation.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
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
                       Extension("MACS2.hashtable", ["MACS2/hashtable.pyx"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=extra_c_args),
                       Extension("MACS2.Statistics", ["MACS2/Statistics.pyx"], libraries=["m"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=extra_c_args),
                       ]
    else:
        ext_modules = [Extension("MACS2.Prob", ["MACS2/Prob.c"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                       Extension("MACS2.IO.Parser",["MACS2/IO/Parser.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.Pileup", ["MACS2/Pileup.c","MACS2/cPosValCalculation.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                       Extension("MACS2.PeakModel", ["MACS2/PeakModel.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.PeakDetect", ["MACS2/PeakDetect.c"], extra_compile_args=extra_c_args),
                       Extension("MACS2.Signal", ["MACS2/Signal.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.IO.PeakIO", ["MACS2/IO/PeakIO.c"], extra_compile_args=extra_c_args),
                       Extension("MACS2.IO.BedGraphIO", ["MACS2/IO/BedGraphIO.c"], extra_compile_args=extra_c_args),                   
                       Extension("MACS2.IO.FixWidthTrack", ["MACS2/IO/FixWidthTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.IO.PairedEndTrack", ["MACS2/IO/PairedEndTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.IO.BedGraph", ["MACS2/IO/BedGraph.c"], libraries=["m"], extra_compile_args=extra_c_args),
                       Extension("MACS2.IO.ScoreTrack", ["MACS2/IO/ScoreTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                       Extension("MACS2.IO.CallPeakUnit", ["MACS2/IO/CallPeakUnit.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                       Extension("MACS2.hashtable", ["MACS2/hashtable.c"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=extra_c_args),
                       Extension("MACS2.Statistics", ["MACS2/Statistics.c", "MACS2/cStatistics.c"], libraries=["m"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=extra_c_args),
                       ]

    with open("README.md", "r") as fh:
        long_description = fh.read()
        
    setup(name="MACS2",
          version="2.1.3.2",
          description="Model Based Analysis for ChIP-Seq data",
          long_description = long_description,
          long_description_content_type="text/markdown",
          author='Tao Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://github.com/taoliu/MACS/',
          package_dir={'MACS2' : 'MACS2'},
          packages=['MACS2', 'MACS2.IO', 'MACS2.data'],
          #package_data={'MACS2': ['data/*.dat']},          
          scripts=['bin/macs2',
                   ],
          classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',              
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python :: 2 :: Only',
              'Programming Language :: Cython',
              ],
          install_requires=[
              'numpy>=1.15',
              'cython>=0.25',
              ],
          cmdclass = command_classes,
          ext_modules = cythonize(ext_modules, language_level="2")
          )

if __name__ == '__main__':
    main()
