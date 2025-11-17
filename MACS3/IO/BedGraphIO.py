# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-02-05 12:38:24 Tao Liu>

"""Utilities for reading and writing MACS3 bedGraph files.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from array import array

from MACS3.Signal.BedGraph import bedGraphTrackI
import cython
from cython.cimports.cpython import bool

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# C lib
# ------------------------------------
from cython.cimports.libc.stdlib import atoi, atof

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class bedGraphIO:
    """Helper for loading and writing bedGraph tracks."""
    bedGraph_filename = cython.declare(str, visibility='public')
    data = cython.declare(object, visibility='public')

    def __init__(self, bedGraph_filename: str, data=None):
        """Initialise the IO wrapper for ``bedGraph_filename``.

        Args:
            bedGraph_filename: Path to the bedGraph file.
            data: Optional existing :class:`bedGraphTrackI` to populate.
        """
        self.bedGraph_filename = bedGraph_filename
        if data:
            assert isinstance(data, bedGraphTrackI)
            self.data = data
        else:
            self.data = bedGraphTrackI()

    @cython.ccall
    def read_bedGraph(self, baseline_value: cython.double = 0):
        """Load bedGraph intervals into the internal :class:`bedGraphTrackI`.

        Args:
            baseline_value: Value used to fill gaps between defined intervals.

        Returns:
            bedGraphTrackI: Populated track instance.
        """
        i: bytes

        self.data.reset_baseline(baseline_value)
        add_func = self.data.add_loc
        # python open file
        bedGraph_file = open(self.bedGraph_filename, "rb")

        for i in bedGraph_file:
            if i.startswith(b"track"):
                continue
            elif i.startswith(b"#"):
                continue
            elif i.startswith(b"browse"):
                continue
            else:
                fs = i.split()
                add_func(fs[0], atoi(fs[1]), atoi(fs[2]), atof(fs[3]))

        bedGraph_file.close()
        return self.data

    @cython.ccall
    def write_bedGraph(self, name: str = "", description: str = "",
                       trackline: bool = True):
        """Persist the current track to ``self.bedGraph_filename``.

        Args:
            name: Track name used in the optional header line.
            description: Track description used in the optional header line.
            trackline: Whether to emit a UCSC ``track`` header.
        """
        pre: cython.int
        pos: cython.int
        i: cython.int
        value: cython.double
        chrom: bytes
        chrs: set
        trackcontents: tuple
        p: array
        v: array

        fhd = open(self.bedGraph_filename, "w")
        if trackline:
            trackcontents = (name.replace("\"", "\\\""),
                             description.replace("\"", "\\\""))
            fhd.write("track type=bedGraph name=\"%s\" description=\"%s\" visibility=2 alwaysZero=on\n" % trackcontents)
        chrs = self.data.get_chr_names()
        for chrom in sorted(chrs):
            (p, v) = self.data.get_data_by_chr(chrom)
            pnext = iter(p).__next__
            vnext = iter(v).__next__
            pre = 0

            for i in range(len(p)):
                pos = pnext()
                value = vnext()
                fhd.write("%s\t%d\t%d\t%.5f\n" %
                          (chrom.decode(), pre, pos, value))
                pre = pos
        fhd.close()
        return
