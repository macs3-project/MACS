# cython: profile=True
# Time-stamp: <2012-04-29 18:04:37 Tao Liu>

"""Module for filter duplicate tags from paired-end data

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
from MACS2.Constants import BYTE4
from logging import debug
import numpy as np
cimport numpy as np
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

def filter_pe_dup ( treat, int maxnum=-1 ):
    """Filter the duplicated reads.

    Run it right after you add all data into this object.
    """
    cdef int32_t i_chrom, i, start, end, n, current_start, current_end
    cdef uint64_t i_old, i_new
    cdef str k
    treatF, treatR = treat
    assert treatF.total == treatR.total, "Different number of fragment starts and ends"

    if maxnum < 0: return # condition to return if not filtering
    if not treatF.__sorted: treatF.sort()
    if not treatR.__sorted: treatR.sort()
    treatF.total = 0
    treatR.total = 0
    chrnames = treatF.get_chr_names()
    for i_chrom in range(len(chrnames)): # for each chromosome
        k = chrnames [ i_chrom ]
        i_new = 0
        starts = treatF.__locations[k][0]
        ends = treatR.__locations[k][0]
        size = starts.shape[0]
        assert size == ends.shape[0], "Lists of 5' and 3' ends are not the same length"
        if size < 1:
            new_starts = starts
            new_ends = ends
            continue
                
        new_starts = np.zeros( treatF.__pointer[k][0],dtype='int32' )
        new_ends = np.zeros( treatR.__pointer[k][0],dtype='int32' )
        new_starts[ i_new ] = starts[ i_new ]
        new_ends[ i_new ] = ends[ i_new ]
        n = 1                # the number of tags in the current location
        
        current_start = starts[0]
        current_end = ends[0]
        for i_old in xrange(1, size):
            start = starts[i_old]
            end = ends[i_old]
            if start == current_start and end == current_end:
                n += 1
                if n <= maxnum:
                    new_starts[ i_new ] = start
                    new_ends[ i_new ] = end
                    i_new += 1
                else:
                    debug("Duplicate fragments found at %s:%d-%d" % (k, start, end) )
            else:
                current_start = start
                current_end = end
                new_starts[ i_new ] = start
                new_ends[ i_new ] = end
                i_new += 1
                n = 1
        new_starts.resize( i_new )
        new_ends.resize( i_new )
        treatF.total += new_starts.shape[0]
        treatR.total += new_ends.shape[0]
       # free memory?
        # I know I should shrink it to 0 size directly,
        # however, on Mac OSX, it seems directly assigning 0
        # doesn't do a thing.
        starts.resize( 100000, refcheck=False )
        starts.resize( 0, refcheck=False )
        ends.resize( 100000, refcheck=False )
        ends.resize( 0, refcheck=False )
        # hope there would be no mem leak...

        treatF.__locations[k]=[new_starts,np.array([], dtype='int32')]
        treatR.__locations[k]=[new_ends,np.array([], dtype='int32')]
    return (treatF, treatR)