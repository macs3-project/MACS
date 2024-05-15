# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-05-15 11:10:27 Tao Liu>

"""Module Description:  IO Module for bedGraph file

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import io

from MACS3.Signal.BedGraph import bedGraphTrackI
from cpython cimport bool
#import numpy as np
#cimport numpy as np
#from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# C lib
# ------------------------------------

from libc.stdio cimport *
from libc.stdlib cimport *

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class bedGraphIO:
    """File IO Class for bedGraph File.

    two publicly available member variables:

    1. bedGraph_filename: the filename for the bedGraph file

    2. data: a bedGraphTrackI class object

    There are two assumptions in the bedGraphTrackI class:

    1. Continuous: the next region should be after the previous one
    unless they are on different chromosomes;

    2. Non-overlapping: the next region should never have overlaps
    with preceding region.

    If any of the above two criteria is violated, parsering will fail.
    """
    cdef:
        public str bedGraph_filename
        public object data

    def __init__ ( self, str bedGraph_filename, object data = None):
        """f must be a filename or a file handler.

        """
        self.bedGraph_filename = bedGraph_filename
        if data:
            assert isinstance( data, bedGraphTrackI )
            self.data = data
        else:
            self.data = bedGraphTrackI()

    cpdef object read_bedGraph (self, double baseline_value=0):
        """Use this function to return a bedGraphTrackI object.

        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with
        baseline_value. Default of baseline_value is 0.
        """
        cdef bytes i

        self.data.reset_baseline( baseline_value )
        add_func = self.data.add_loc
        # python open file
        bedGraph_file = open( self.bedGraph_filename, "rb" )

        for i in bedGraph_file:
            if i.startswith(b"track"):
                continue
            elif i.startswith(b"#"):
                continue
            elif i.startswith(b"browse"):
                continue
            else:
                fs = i.split()
                add_func(fs[0],atoi(fs[1]),atoi(fs[2]),atof(fs[3]))

        bedGraph_file.close()
        return self.data

    cpdef void write_bedGraph (self, str name = "", str description = "", bool trackline = True):
        """Write all data to self.bedGraph_filename in bedGraph Format.

        name/description: the name and description in track line.
        """
        cdef:
            int pre, pos, i
            double value
            bytes chrom
            set chrs
            tuple trackcontents

        fhd = open( self.bedGraph_filename, "w" )
        if trackline:
            trackcontents = (name.replace("\"", "\\\""), description.replace("\"", "\\\""))
            fhd.write("track type=bedGraph name=\"%s\" description=\"%s\" visibility=2 alwaysZero=on\n" % trackcontents)
        chrs = self.data.get_chr_names()
        for chrom in sorted(chrs):
            (p,v) = self.data.get_data_by_chr(chrom)
            pnext = iter(p).__next__
            vnext = iter(v).__next__
            pre = 0

            for i in range(len(p)):
                pos = pnext()
                value = vnext()
                fhd.write("%s\t%d\t%d\t%.5f\n" % (chrom.decode(),pre,pos,value))
                pre = pos
        fhd.close()
        return

