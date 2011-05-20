#!/usr/bin/env python
# Time-stamp: <2011-05-19 12:02:37 Tao Liu>

"""Description

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

Copyright (c) 2008,2009,2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  beta
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

import os
import sys
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

def main():
    if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
        sys.exit(1)

    ext_modules = [Extension("MACS2.cProb", ["lib/cProb.pyx"], libraries=["m"]),
                   Extension("MACS2.cPileup", ["lib/cPileup.pyx"]),
                   Extension("MACS2.cPeakDetect", ["lib/cPeakDetect.pyx"]),                   
                   Extension("MACS2.IO.cFeatIO", ["lib/IO/cFeatIO.pyx"],),
                   ]
    

    setup(name="MACS",
          version="2.0.1",
          description="Model Based Analysis for ChIP-Seq data",
          author='Yong Zhang; Tao (Foo) Liu',
          author_email='zy@jimmy.harvard.edu; taoliu@jimmy.harvard.edu',
          url='http://liulab.dfci.harvard.edu/MACS/',
          package_dir={'MACS2' : 'lib'},
          packages=['MACS2', 'MACS2.IO'],
          scripts=['bin/macs2',
                   #'bin/elandmulti2bed','bin/elandresult2bed','bin/elandexport2bed','bin/sam2bed',
                   'bin/bdgpeakcall'],
          classifiers=[
              'Development Status :: 4 - experimental',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'License :: OSI Approved :: Artistic License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Programming Language :: Python',
              ],
          cmdclass = {'build_ext': build_ext},
          ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()
