# Time-stamp: <2011-02-25 22:20:05 Tao Liu>

"""Module Description: BinKeeper for Wiggle-like tracks.

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

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

import sys
import re
from bisect import insort,bisect_left,bisect_right
from array import array
# ------------------------------------
# constants
# ------------------------------------
# to determine the byte size
if array('H',[1]).itemsize == 2:
    BYTE2 = 'H'
else:
    raise Exception("BYTE2 type cannot be determined!")

if array('I',[1]).itemsize == 4:
    BYTE4 = 'I'
elif array('L',[1]).itemsize == 4:
    BYTE4 = 'L'
else:
    raise Exception("BYTE4 type cannot be determined!")

if array('f',[1]).itemsize == 4:
    FBYTE4 = 'f'
elif array('d',[1]).itemsize == 4:
    FBYTE4 = 'd'
else:
    raise Exception("BYTE4 type cannot be determined!")

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
class BinKeeperI:
    """BinKeeper keeps point data from a chromosome in a bin list.

    Example:
    >>> from taolib.CoreLib.Parser import WiggleIO
    >>> w = WiggleIO('sample.wig')
    >>> bk = w.build_binKeeper()
    >>> bk['chrI'].pp2v(1000,2000) # to extract values in chrI:1000..2000
    """
    def __init__ (self,binsize=8000,chromosomesize=1e9):
        """Initializer.

        Parameters:
        binsize : size of bin in Basepair
        chromosomesize : size of chromosome, default is 1G
        """
        self.binsize = binsize
        self.binnumber = int(chromosomesize/self.binsize)+1
        self.cage = []
        a = self.cage.append
        for i in xrange(self.binnumber):
            a([array(BYTE4,[]),array(FBYTE4,[])])

    def add ( self, p, value ):
        """Add a position into BinKeeper.

        Note: position must be sorted before adding. Otherwise, pp2v
        and pp2p will not work.
        """
        bin = p/self.binsize
        self.cage[bin][0].append(p)
        self.cage[bin][1].append(value)        

    def p2bin (self, p ):
        """Return the bin index for a position.
        
        """
        return p/self.binsize

    def p2cage (self, p):
        """Return the bin containing the position.
        
        """
        return self.cage[p/self.binsize]

    def __pp2cages (self, p1, p2):
        assert p1<=p2
        bin1 = self.p2bin(p1)
        bin2 = self.p2bin(p2)+1
        t = [array(BYTE4,[]),array(FBYTE4,[])]
        for i in xrange(bin1,bin2):
            t[0].extend(self.cage[i][0])
            t[1].extend(self.cage[i][1])            
        return t

    def pp2p (self, p1, p2):
        """Give the position list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        list of positions between p1 and p2.
        """
        (ps,vs) = self.__pp2cages(p1,p2)
        p1_in_cages = bisect_left(ps,p1)
        p2_in_cages = bisect_right(ps,p2)
        return ps[p1_in_cages:p2_in_cages]

    def pp2v (self, p1, p2):
        """Give the value list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        list of values whose positions are between p1 and p2.
        """
        (ps,vs) = self.__pp2cages(p1,p2)
        p1_in_cages = bisect_left(ps,p1)
        p2_in_cages = bisect_right(ps,p2)
        return vs[p1_in_cages:p2_in_cages]


    def pp2pv (self, p1, p2):
        """Give the (position,value) list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        list of (position,value) between p1 and p2.
        """
        (ps,vs) = self.__pp2cages(p1,p2)
        p1_in_cages = bisect_left(ps,p1)
        p2_in_cages = bisect_right(ps,p2)
        return zip(ps[p1_in_cages:p2_in_cages],vs[p1_in_cages:p2_in_cages])
