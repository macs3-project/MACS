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


def main():
    if sys.version_info < (3, 9):
        sys.stderr.write("CRITICAL: Python version must >= 3.9!\n")
        sys.exit(1)

    # NumPy include dir
    numpy_include_dir = [numpy.get_include()]

    # CFLAG
    extra_c_args = ["-w", "-O3", "-g0"]
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
            gcc_version_check = subprocess.check_output([os.environ.get('CC', 'gcc'), "--version"], universal_newlines=True)
            if gcc_version_check.find("clang") != -1:
                clang = True
            else:
                gcc_version_check = gcc_version_check.split('\n')[0]  # get the first line
                m = re.search(r"\s+(\d+\.\d+)\.\d+", gcc_version_check)
                if m:
                    gcc_version = float(m[1])
                    if gcc_version > 4.8:
                        new_gcc = True
        except subprocess.CalledProcessError:
            pass

    try:
        if os.environ['CC'] == "icc":
            icc = True
    except KeyError:
        pass

    extra_c_args_for_fermi = ["-std=gnu99", "-DUSE_SIMDE",
                              "-DSIMDE_ENABLE_NATIVE_ALIASES"]

    if icc or sysconfig.get_config_vars()['CC'] == 'icc':
        extra_c_args_for_fermi.extend(['-qopenmp-simd',
                                       '-DSIMDE_ENABLE_OPENMP'])
    elif new_gcc or clang or sysconfig.get_config_vars()['CC'] == 'clang':
        extra_c_args_for_fermi.extend(['-fopenmp-simd',
                                       '-DSIMDE_ENABLE_OPENMP'])

    # extensions, those have to be processed by Cython
    ext_modules = [Extension("MACS3.Signal.HMMR_EM",
                             ["MACS3/Signal/HMMR_EM.py"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.HMMR_Signal_Processing",
                             ["MACS3/Signal/HMMR_Signal_Processing.py"],
                             libraries=["m"],
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.HMMR_HMM",
                             ["MACS3/Signal/HMMR_HMM.py"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.Prob",
                             ["MACS3/Signal/Prob.py"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.Region",
                             ["MACS3/Signal/Region.py"],
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.Pileup",
                             ["MACS3/Signal/Pileup.py",
                              "MACS3/Signal/cPosValCalculation.c"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.PileupV2",
                             ["MACS3/Signal/PileupV2.py"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.PeakModel",
                             ["MACS3/Signal/PeakModel.pyx"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.PeakDetect",
                             ["MACS3/Signal/PeakDetect.pyx"],
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.SignalProcessing",
                             ["MACS3/Signal/SignalProcessing.pyx"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.FixWidthTrack",
                             ["MACS3/Signal/FixWidthTrack.pyx"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.PairedEndTrack",
                             ["MACS3/Signal/PairedEndTrack.pyx"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.BedGraph",
                             ["MACS3/Signal/BedGraph.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.ScoreTrack",
                             ["MACS3/Signal/ScoreTrack.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.CallPeakUnit",
                             ["MACS3/Signal/CallPeakUnit.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.VariantStat",
                             ["MACS3/Signal/VariantStat.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.ReadAlignment",
                             ["MACS3/Signal/ReadAlignment.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.RACollection",
                             ["MACS3/Signal/RACollection.pyx",
                              "MACS3/fermi-lite/bfc.c",
                              "MACS3/fermi-lite/bseq.c",
                              "MACS3/fermi-lite/bubble.c",
                              "MACS3/fermi-lite/htab.c",
                              "MACS3/fermi-lite/ksw.c",
                              "MACS3/fermi-lite/kthread.c",
                              "MACS3/fermi-lite/mag.c",
                              "MACS3/fermi-lite/misc.c",
                              "MACS3/fermi-lite/mrope.c",
                              "MACS3/fermi-lite/rld0.c",
                              "MACS3/fermi-lite/rle.c",
                              "MACS3/fermi-lite/rope.c",
                              "MACS3/fermi-lite/unitig.c",
                              "MACS3/Signal/swalign.c"],
                             libraries=["m", "z"],
                             include_dirs=numpy_include_dir+["./",
                                                             "./MACS3/fermi-lite/",
                                                             "./MACS3/Signal/"],
                             extra_compile_args=extra_c_args+extra_c_args_for_fermi),
                   Extension("MACS3.Signal.UnitigRACollection",
                             ["MACS3/Signal/UnitigRACollection.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.PosReadsInfo",
                             ["MACS3/Signal/PosReadsInfo.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.Signal.PeakVariants",
                             ["MACS3/Signal/PeakVariants.pyx"],
                             libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.IO.Parser",
                             ["MACS3/IO/Parser.py"],
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.IO.PeakIO",
                             ["MACS3/IO/PeakIO.pyx"],
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.IO.BedGraphIO",
                             ["MACS3/IO/BedGraphIO.pyx"],
                             extra_compile_args=extra_c_args),
                   Extension("MACS3.IO.BAM",
                             ["MACS3/IO/BAM.pyx",], libraries=["m"],
                             include_dirs=numpy_include_dir,
                             extra_compile_args=extra_c_args)]

    setup(version=MACS_VERSION,
          package_dir={'MACS3': 'MACS3'},
          packages=['MACS3', 'MACS3.IO', 'MACS3.Signal', 'MACS3.Commands', 'MACS3.Utilities'],
          package_data={'MACS3': ['*.pxd']},
          scripts=['bin/macs3', ],
          ext_modules=cythonize(ext_modules))


if __name__ == '__main__':
    main()
