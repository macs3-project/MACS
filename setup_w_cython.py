#!/usr/bin/env python
# Time-stamp: <2014-10-29 02:00:22 Tao Liu>

"""Description: 

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

Use this when you need Cython regenerate .c files.

Copyright (c) 2008,2009,2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  beta
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

import os
import sys
from setuptools import setup, Extension

# Use build_ext from Cython if found
command_classes = {}
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
    has_cython = True
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

    if has_cython:
        ext_modules = [Extension("MACS2.cProb", ["MACS2/cProb.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                       Extension("MACS2.IO.cParser",["MACS2/IO/cParser.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.cPileup", ["MACS2/cPileup.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                       Extension("MACS2.cArray", ["MACS2/cArray.pyx"], extra_compile_args=["-w","-Ofast"]),                       
                       Extension("MACS2.cPeakModel", ["MACS2/cPeakModel.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.cPeakDetect", ["MACS2/cPeakDetect.pyx"], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.cSignal", ["MACS2/cSignal.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cPeakIO", ["MACS2/IO/cPeakIO.pyx"], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cBedGraphIO", ["MACS2/IO/cBedGraphIO.pyx"], extra_compile_args=["-w","-Ofast"]),                   
                       Extension("MACS2.IO.cFixWidthTrack", ["MACS2/IO/cFixWidthTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cPairedEndTrack", ["MACS2/IO/cPairedEndTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cBedGraph", ["MACS2/IO/cBedGraph.pyx"], libraries=["m"], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cScoreTrack", ["MACS2/IO/cScoreTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                       Extension("MACS2.IO.cCallPeakUnit", ["MACS2/IO/cCallPeakUnit.pyx"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.hashtable", ["MACS2/hashtable.pyx"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.Statistics", ["MACS2/Statistics.pyx", "MACS2/cStatistics.c"], libraries=["m"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=["-w","-Ofast"]),
                       ]
    else:
        ext_modules = [Extension("MACS2.cProb", ["MACS2/cProb.c"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                       Extension("MACS2.IO.cParser",["MACS2/IO/cParser.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.cPileup", ["MACS2/cPileup.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                       Extension("MACS2.cArray", ["MACS2/cArray.c"], extra_compile_args=["-w","-Ofast"]),                       
                       Extension("MACS2.cPeakModel", ["MACS2/cPeakModel.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.cPeakDetect", ["MACS2/cPeakDetect.c"], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.cSignal", ["MACS2/cSignal.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cPeakIO", ["MACS2/IO/cPeakIO.c"], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cBedGraphIO", ["MACS2/IO/cBedGraphIO.c"], extra_compile_args=["-w","-Ofast"]),                   
                       Extension("MACS2.IO.cFixWidthTrack", ["MACS2/IO/cFixWidthTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cPairedEndTrack", ["MACS2/IO/cPairedEndTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cBedGraph", ["MACS2/IO/cBedGraph.c"], libraries=["m"], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.IO.cScoreTrack", ["MACS2/IO/cScoreTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                       Extension("MACS2.IO.cCallPeakUnit", ["MACS2/IO/cCallPeakUnit.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.hashtable", ["MACS2/hashtable.c"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=["-w","-Ofast"]),
                       Extension("MACS2.Statistics", ["MACS2/Statistics.c", "MACS2/cStatistics.c"], libraries=["m"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=["-w","-Ofast"]),
                       ]

    setup(name="MACS2",
          version="2.1.0.20141030",
          description="Model Based Analysis for ChIP-Seq data",
          author='Tao Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://github.com/taoliu/MACS/',
          package_dir={'MACS2' : 'MACS2'},
          packages=['MACS2', 'MACS2.IO', 'MACS2.data'],
          package_data={'MACS2': ['data/*.dat']},          
          scripts=['bin/macs2',
                   ],
          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',              
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python',
              ],
          install_requires=[
              'numpy>=1.6',
              'cython>=0.18',
              #'scipy',
              ],
          cmdclass = command_classes,
          ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()
