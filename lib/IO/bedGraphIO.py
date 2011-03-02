# Time-stamp: <2011-03-02 17:28:30 Tao Liu>

"""Module Description

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
import os
import sys
import re
import shutil
from MACS14.IO.FeatIO import bedGraphTrackI
from MACS14.IO.BinKeeper import BinKeeperII

import time
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class bedGraphIO:
    """File Parser Class for bedGraph File.

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

    def build_bdgtrack (self):
        """Use this function to return a bedGraphTrackI object.

        """
        data = bedGraphTrackI()
        add_func = data.add_loc
        for i in self.fhd:
            if i.startswith("track"):
                continue
            elif i.startswith("#"):
                continue
            elif i.startswith("browse"):
                continue
            else:
                (chrom,startpos,endpos,value)=i.split()
                add_func(chrom,int(startpos),int(endpos),float(value))
        self.fhd.seek(0)
        return data

    def build_binKeeper (self,chromLenDict={},binsize=200):
        """Use this function to return a dictionary of BinKeeperII
        objects.

        chromLenDict is a dictionary for chromosome length like

        {'chr1':100000,'chr2':200000}

        bin is in bps. for detail, check BinKeeper.
        """
        data = {}

        for i in self.fhd:
            if i.startswith("track"):
                continue
            elif i.startswith("#"):
                continue
            elif i.startswith("browse"):
                continue
            else:
                (chrom,startpos,endpos,value)=i.split()

                if not data.has_key(chrom):
                    chrlength = chromLenDict.setdefault(chrom,250000000) + 10000000
                    data.setdefault(chrom,BinKeeperII(binsize=binsize,chromosomesize=chrlength))

                data[chrom].add(int(startpos),int(endpos),float(value))

        self.fhd.seek(0)
        return data

