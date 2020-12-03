# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-12-03 11:46:58 Tao Liu>

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
    """File Parser Class for bedGraph File.

    There are two assumptions in my bedGraphTrackI object:

    1. Continuous: the next region should be after the previous one
    unless they are on different chromosomes;

    2. Non-overlapping: the next region should never have overlaps
    with preceding region.

    If any of the above two criteria is violated, parsering will fail.
    """
    cdef:
        str bedGraph_filename

    def __init__ ( self, str bedGraph_filename ):
        """f must be a filename or a file handler.

        """
        self.bedGraph_filename = bedGraph_filename

    def build_bdgtrack (self, double baseline_value=0):
        """Use this function to return a bedGraphTrackI object.

        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with
        baseline_value. Default of baseline_value is 0.
        """
        cdef bytes i

        data = bedGraphTrackI(baseline_value=baseline_value)
        add_func = data.add_loc
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
        return data

