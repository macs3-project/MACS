#!/usr/bin/env python
# Time-stamp: <2015-04-20 14:32:39 Tao Liu>

"""Description

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

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
    from numpy import get_include as numpy_get_include 
    numpy_include_dir = [numpy_get_include()] 
except: 
    numpy_include_dir = []
    sys.stderr.write("CRITICAL:Numpy must be installed!\n")
    sys.exit(1)

def main():
    if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

    ext_modules = [Extension("MACS2.Prob", ["MACS2/Prob.c"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                   Extension("MACS2.IO.Parser",["MACS2/IO/Parser.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.Pileup", ["MACS2/Pileup.c","MACS2/cPosValCalculation.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                   Extension("MACS2.PeakModel", ["MACS2/PeakModel.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.PeakDetect", ["MACS2/PeakDetect.c"], extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.Signal", ["MACS2/Signal.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.IO.PeakIO", ["MACS2/IO/PeakIO.c"], extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.IO.BedGraphIO", ["MACS2/IO/BedGraphIO.c"], extra_compile_args=["-w","-Ofast"]),                   
                   Extension("MACS2.IO.FixWidthTrack", ["MACS2/IO/FixWidthTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.IO.PairedEndTrack", ["MACS2/IO/PairedEndTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.IO.BedGraph", ["MACS2/IO/BedGraph.c"], libraries=["m"], extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.IO.ScoreTrack", ["MACS2/IO/ScoreTrack.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"] ),
                   Extension("MACS2.IO.CallPeakUnit", ["MACS2/IO/CallPeakUnit.c"], include_dirs=numpy_include_dir, extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.hashtable", ["MACS2/hashtable.c"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=["-w","-Ofast"]),
                   Extension("MACS2.Statistics", ["MACS2/Statistics.c", "MACS2/cStatistics.c"], libraries=["m"], include_dirs=["MACS2/",numpy_get_include()], extra_compile_args=["-w","-Ofast"]),
                   ]
    
    setup(name="MACS2",
          version="2.1.0.20150420",
          description="Model Based Analysis for ChIP-Seq data",
          author='Tao Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://github.com/taoliu/MACS/',
          package_dir={'MACS2' : 'MACS2'},
          packages=['MACS2', 'MACS2.IO'],#, 'MACS2.data'],
          #package_data={'MACS2': ['data/*.dat']},          
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
              #'scipy',
              ],
          cmdclass = command_classes,
          ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()
