# Time-stamp: <2011-03-14 17:52:00 Tao Liu>

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
from bisect import insort,bisect_left,bisect_right,insort_right
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


class BinKeeperII:
    """BinKeeperII keeps non-overlapping interval data from a chromosome in a bin list.

    This is especially designed for bedGraph type data.

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
            a([array(BYTE4,[]),array(BYTE4,[]),array(FBYTE4,[])])

    def add ( self, startp, endp, value ):
        """Add an interval data into BinKeeper.

        Note: position must be sorted before adding. Otherwise, pp2v
        and pp2p will not work.
        """
        startbin = startp/self.binsize
        endbin = endp/self.binsize
        if startbin == endbin:
            # some intervals may only be within a bin
            j = bisect.bisect_left(self.cage[startbin][0],startp)
            self.cage[startbin][0].insert(j,startp)
            self.cage[startbin][1].insert(j,endp)
            self.cage[startbin][2].insert(j,value)
        else:
            # some intervals may cover the end of bins
            # first bin
            j = bisect.bisect_left(self.cage[startbin][0],startp)
            self.cage[startbin][0].insert(j,startp)
            self.cage[startbin][1].insert(j,(startbin+1)*self.binsize)
            self.cage[startbin][2].insert(j,value)
            # other bins fully covered
            for i in xrange(startbin+1,endbin):
                p = i*self.binsize
                j = bisect.bisect_left(self.cage[startbin][0],p)
                self.cage[startbin][0].insert(j,p)
                self.cage[startbin][1].insert(j,(i+1)*self.binsize)
                self.cage[startbin][2].insert(j,value)

                insort_right(self.cage[i][0],i*self.binsize)
                insort_right(self.cage[i][1],(i+1)*self.binsize)
                insort_right(self.cage[i][2],value)
            # last bin -- the start of this bin should be covered
            insort_right(self.cage[endbin][0],endbin*self.binsize)
            insort_right(self.cage[endbin][1],endp)
            insort_right(self.cage[endbin][2],value)

    def p2bin (self, p ):
        """Given a position, return the bin index for a position.
        
        """
        return p/self.binsize

    def p2cage (self, p):
        """Given a position, return the bin containing the position.
        
        """
        return self.cage[p/self.binsize]

    def pp2cages (self, p1, p2):
        """Given an interval, return the bins containing this interval.
        
        """
        assert p1<=p2
        bin1 = self.p2bin(p1)
        bin2 = self.p2bin(p2)
        t = [array(BYTE4,[]),array(BYTE4,[]),array(FBYTE4,[])]
        for i in xrange(bin1,bin2+1):
            t[0].extend(self.cage[i][0])
            t[1].extend(self.cage[i][1])
            t[2].extend(self.cage[i][2])                        
        return t

    def pp2intervals (self, p1, p2):
        """Given an interval, return the intervals list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        A list of intervals start and end positions (tuple) between p1 and p2.

        * Remember, I assume all intervals saved in this BinKeeperII
          are not overlapping, so if there is some overlap, this
          function will not work as expected.
        """
        (startposs,endposs,vs) = self.pp2cages(p1,p2)
        p1_in_cages = bisect_left(startposs,p1)
        p2_in_cages = bisect_right(endposs,p2)
        output_startpos_list = startposs[p1_in_cages:p2_in_cages]
        output_endpos_list = endposs[p1_in_cages:p2_in_cages]

        # check if the bin (p1_in_cages-1) covers p1
        if p1 < endposs[p1_in_cages-1]:
            # add this interval
            output_startpos_list = array(BYTE4,[p1,])+output_startpos_list
            output_endpos_list = array(BYTE4,[endposs[p1_in_cages-1],])+output_endpos_list

        # check if the bin (p2_in_cages+1) covers p2
        if p2 > startposs[p2_in_cages+1]:
            # add this interval
            output_startpos_list = array(BYTE4,[startposs[p2_in_cages+1],])+output_startpos_list
            output_endpos_list = array(BYTE4,[p2,])+output_endpos_list

        return zip(output_startpos_list,output_endpos_list)

    def pp2pvs (self, p1, p2):
        """Given an interval, return the values list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:

        A list of start, end positions, values (tuple) between p1 and
        p2. Each value represents the value in an interval. Remember
        the interval length and positions are lost in the output.

        * Remember, I assume all intervals saved in this BinKeeperII
          are not overlapping, so if there is some overlap, this
          function will not work as expected.
        """

        (startposs,endposs,vs) = self.pp2cages(p1,p2)
        p1_in_cages = bisect_left(startposs,p1)
        p2_in_cages = bisect_right(endposs,p2)
        output_startpos_list = startposs[p1_in_cages:p2_in_cages]
        output_endpos_list = endposs[p1_in_cages:p2_in_cages]
        output_value_list = vs[p1_in_cages:p2_in_cages]

        # print p1_in_cages,p2_in_cages
        # print vs
        print output_startpos_list
        print output_endpos_list
        print output_value_list

        # check if the bin (p1_in_cages-1) covers p1
        
        if p1_in_cages-1 >= 0 and p1 < self.cage[p1_in_cages-1][1]:
            # add this interval
            output_startpos_list = array(BYTE4,[p1,])+output_startpos_list
            output_endpos_list = array(BYTE4,[self.cage[p1_in_cages-1][1],])+output_endpos_list
            output_value_list = array(BYTE4,[self.cage[p1_in_cages-1][2],])+output_value_list
                

        # check if the bin (p2_in_cages+1) covers p2
        #print p2_in_cages+1,len(self.cage)
        #print p2, self.cage[p2_in_cages+1][0]
        if p2_in_cages+1 < len(self.cage) and p2 > self.cage[p2_in_cages+1][0]:
            # add this interval
            output_startpos_list = output_startpos_list+array(BYTE4,[self.cage[p2_in_cages+1][0],])
            output_endpos_list = output_endpos_list+array(BYTE4,[p2,])
            output_value_list = output_value_list+array(BYTE4,[self.cage[p2_in_cages+1][2],])

        print output_startpos_list
        print output_endpos_list
        print output_value_list

        return zip(output_startpos_list,output_endpos_list,output_value_list)

