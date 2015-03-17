# Time-stamp: <2015-03-13 12:50:14 Tao Liu>

"""Module Description:  IO Module for bedGraph file

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import io

from MACS2.IO.BedGraph import bedGraphTrackI,bedRegionTrackI

# ------------------------------------
# constants
# ------------------------------------
# cdef extern from "stdio.h":
#     ctypedef struct FILE
#     FILE *fopen   (const char *filename, const char  *opentype)
#     int fclose   (FILE *stream)
#     int fscanf   (FILE *stream, const char *template, ...)
#     int fprintf  (FILE *stream, const char *template, ...)
#     enum: EOF

# cdef extern from "stdlib.h":
#     ctypedef unsigned int size_t
#     size_t strlen(char *s)
#     void *malloc(size_t size)
#     void *calloc(size_t n, size_t size)
#     void free(void *ptr)
#     int strcmp(char *a, char *b)
#     char * strcpy(char *a, char *b)
#     long atol(char *str)
#     int atoi(char *str)
#     double atof(char *str)

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
        char * bedGraph_filename

    def __init__ ( self, bedGraph_filename ):
        """f must be a filename or a file handler.
        
        """
        self.bedGraph_filename = <bytes> bedGraph_filename

    def build_bdgtrack (self, double baseline_value=0):
        """Use this function to return a bedGraphTrackI object.

        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with
        baseline_value. Default of baseline_value is 0.
        """
        cdef str i

        data = bedGraphTrackI(baseline_value=baseline_value)
        add_func = data.add_loc
        # python open file
        bedGraph_file = open( self.bedGraph_filename, "r" )
        
        for i in bedGraph_file:
            if i.startswith("track"):
                continue
            elif i.startswith("#"):
                continue
            elif i.startswith("browse"):
                continue
            else:
                fs = i.split()
                add_func(fs[0],atoi(fs[1]),atoi(fs[2]),atof(fs[3]))

        bedGraph_file.close()
        return data

    def build_bdgtrack_c (self, double baseline_value=0):
        """Use this function to return a bedGraphTrackI object. This is a C version.

        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with
        baseline_value. Default of baseline_value is 0.
        """
        cdef:
            char[80] chrom      # will there be chromosome name longer than 80 characters?
            int start, end
            float value
            FILE * bedGraph_file
        
        data = bedGraphTrackI(baseline_value=baseline_value)
        add_func = data.add_loc
        # C open file
        bedGraph_file = fopen( self.bedGraph_filename, "r" )
        if bedGraph_file == NULL:
            raise Exception("File '%s' can not be opened!" % self.bedGraph_filename )
        
        while ( fscanf( bedGraph_file, "%s\t%d\t%d\t%f\n", chrom, &start, &end, &value ) != EOF ):
            add_func( chrom, start, end, value )
        fclose( bedGraph_file )
        #data.finalize()
        return data



cdef class genericBedIO:
    """File Parser Class for generic bed File with at least column #1,#2,#3,and #5.

    There are two assumptions in my bedGraphTrackI object:

    1. Continuous: the next region should be after the previous one
    unless they are on different chromosomes;
    
    2. Non-overlapping: the next region should never have overlaps
    with preceding region.

    If any of the above two criteria is violated, parsering will
    fail. You'd better use it to read peak file from MACS. Or sort BED
    by chromosome and start position.
    """
    def __init__ (self,f):
        """f must be a filename or a file handler.
        
        """
        if type(f) == str:
            self.fhd = open(f,"r")
        elif type(f) == file:
            self.fhd = f
        else:
            raise Exception("f must be a filename or a file handler.")

    def build_bedtrack (self):
        """Use this function to return a bedGraphTrackI object.

        baseline_value is the value to fill in the regions not defined
        in bedGraph. For example, if the bedGraph is like:

        chr1  100 200  1
        chr1  250 350  2

        Then the region chr1:200..250 should be filled with
        baseline_value. Default of baseline_value is 0.
        """
        cdef str i
        
        data = bedRegionTrackI() #(baseline_value=baseline_value)
        add_func = data.safe_add_loc
        chrom_itemcount = {}

        self.fhd.seek(0)
        
        for i in self.fhd:
            fs = i.split()
            add_func(fs[0],atoi(fs[1]),atoi(fs[2])) #,float(value))
        self.fhd.seek(0)
        return data



