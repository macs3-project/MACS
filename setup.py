#!/usr/bin/env python3
"""Description:

Setup script for MACS -- Model Based Analysis for ChIP-Seq data

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import os
import re
from setuptools import setup, Extension
from Cython.Build import cythonize
import subprocess
import sysconfig
import numpy

# get MACS version
exec(open("MACS3/Utilities/Constants.py").read())

# classifiers
classifiers =[\
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python :: 3.9',
              'Programming Language :: Python :: 3.10',
              'Programming Language :: Python :: 3.11',
              'Programming Language :: Cython', ]

install_requires = [ "numpy>=1.19",
                     "hmmlearn>=0.3",
                     "cykhash>=2.0",
                     "Cython>=0.29" ]

tests_requires = [ 'pytest' ]


def main():
    if sys.version_info < (3,8):
        sys.stderr.write("CRITICAL: Python version must >= 3.8!\n")
        sys.exit(1)

    # NumPy include dir
    numpy_include_dir = [ numpy.get_include() ]

    # CFLAG
    # I intend to use -Ofast, however if gcc version < 4.6, this option is unavailable so...
    # should I care about old gcc compiler?...
    extra_c_args = ["-w","-Ofast", "-g0"] # for C, -Ofast implies -O3 and -ffast-math

    # CFLAG for fermi-lite related codes
    clang = False
    icc = False
    new_gcc = False
    try:
        if os.environ['CC'] == "clang":
            clang = True
    except KeyError:
        pass

    if not clang:
        try:
            gcc_version_check = subprocess.check_output( ["gcc", "--version"], universal_newlines=True)
            if gcc_version_check.find("clang") != -1:
                clang = True
            else:
                gcc_version_check = gcc_version_check.split('\n')[0] # get the first line
                m = re.search( "\s+(\d+\.\d+)\.\d+", gcc_version_check )
                if m:
                    gcc_version = float( m[1] )
                    if gcc_version > 4.8:
                        new_gcc = True
        except subprocess.CalledProcessError:
            pass

    try:
        if os.environ['CC'] == "icc":
            icc = True
    except KeyError:
        pass

    extra_c_args_for_fermi = ["-std=gnu99","-DUSE_SIMDE", "-DSIMDE_ENABLE_NATIVE_ALIASES"]
    if icc or sysconfig.get_config_vars()['CC'] == 'icc':
        extra_c_args_for_fermi.extend(['-qopenmp-simd', '-DSIMDE_ENABLE_OPENMP'])
    elif new_gcc or clang or sysconfig.get_config_vars()['CC'] == 'clang':
        extra_c_args_for_fermi.extend(['-fopenmp-simd', '-DSIMDE_ENABLE_OPENMP'])

    # extensions, those have to be processed by Cython
    ext_modules = [ \
                    Extension("MACS3.Signal.HMMR_EM", ["MACS3/Signal/HMMR_EM.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.HMMR_Signal_Processing", ["MACS3/Signal/HMMR_Signal_Processing.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.HMMR_HMM", ["MACS3/Signal/HMMR_HMM.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.Prob", ["MACS3/Signal/Prob.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.Region", ["MACS3/Signal/Region.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.Pileup", ["MACS3/Signal/Pileup.pyx","MACS3/Signal/cPosValCalculation.c"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.PileupV2", ["MACS3/Signal/PileupV2.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),                    
                    Extension("MACS3.Signal.PeakModel", ["MACS3/Signal/PeakModel.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.PeakDetect", ["MACS3/Signal/PeakDetect.pyx"], extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.SignalProcessing", ["MACS3/Signal/SignalProcessing.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.FixWidthTrack", ["MACS3/Signal/FixWidthTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.PairedEndTrack", ["MACS3/Signal/PairedEndTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.BedGraph", ["MACS3/Signal/BedGraph.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    #Extension("MACS3.Signal.BedGraphV2", ["MACS3/Signal/BedGraphV2.pyx"], libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),                    
                    Extension("MACS3.Signal.ScoreTrack", ["MACS3/Signal/ScoreTrack.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args ),
                    Extension("MACS3.Signal.CallPeakUnit", ["MACS3/Signal/CallPeakUnit.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.VariantStat",["MACS3/Signal/VariantStat.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.ReadAlignment",["MACS3/Signal/ReadAlignment.pyx",],libraries=["m"],include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.RACollection",["MACS3/Signal/RACollection.pyx","MACS3/fermi-lite/bfc.c","MACS3/fermi-lite/bseq.c",\
                                                           "MACS3/fermi-lite/bubble.c","MACS3/fermi-lite/htab.c","MACS3/fermi-lite/ksw.c","MACS3/fermi-lite/kthread.c",\
                                                           "MACS3/fermi-lite/mag.c","MACS3/fermi-lite/misc.c","MACS3/fermi-lite/mrope.c","MACS3/fermi-lite/rld0.c",\
                                                           "MACS3/fermi-lite/rle.c","MACS3/fermi-lite/rope.c","MACS3/fermi-lite/unitig.c", "MACS3/Signal/swalign.c" ], \
                              libraries=["m","z"], include_dirs=numpy_include_dir+["./","./MACS3/fermi-lite/","./MACS3/Signal/"], extra_compile_args=extra_c_args+extra_c_args_for_fermi),
                    Extension("MACS3.Signal.UnitigRACollection",["MACS3/Signal/UnitigRACollection.pyx"],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.PosReadsInfo",["MACS3/Signal/PosReadsInfo.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.Signal.PeakVariants",["MACS3/Signal/PeakVariants.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    # IO
                    Extension("MACS3.IO.Parser",["MACS3/IO/Parser.pyx"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args),
                    Extension("MACS3.IO.PeakIO", ["MACS3/IO/PeakIO.pyx"], extra_compile_args=extra_c_args),
                    Extension("MACS3.IO.BedGraphIO", ["MACS3/IO/BedGraphIO.pyx"], extra_compile_args=extra_c_args),
                    Extension("MACS3.IO.BAM",["MACS3/IO/BAM.pyx",],libraries=["m"], include_dirs=numpy_include_dir, extra_compile_args=extra_c_args) ]


    with open("README.md", "r") as fh:
        long_description = fh.read()

    setup( name = "MACS3",
           version = MACS_VERSION,
           description = "Model Based Analysis for ChIP-Seq data",
           long_description = long_description,
           long_description_content_type = "text/markdown",
           author = 'Tao Liu',
           author_email = 'vladimir.liu@gmail.com',
           url = 'http://github.com/macs3-project/MACS/',
           package_dir = {'MACS3' : 'MACS3'},
           packages = ['MACS3', 'MACS3.IO', 'MACS3.Signal', 'MACS3.Commands','MACS3.Utilities'],
           package_data = {'MACS3':['*.pxd']},
           scripts = ['bin/macs3', ],
           classifiers = classifiers,
           install_requires = install_requires,
           setup_requires = install_requires,
           tests_require = tests_requires,
           python_requires = '>=3.9',
           ext_modules=cythonize( ext_modules ) ),
#                                  compiler_directives={'linetrace': True, 'binding': True}) )

if __name__ == '__main__':
    main()
