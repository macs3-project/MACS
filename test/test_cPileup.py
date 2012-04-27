#!/usr/bin/env python
# Time-stamp: <2012-04-25 17:19:57 Tao Liu>

"""Module Description: Test functions for pileup functions.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""


import os
import sys
import unittest

from math import log10
from MACS2.cPileup import *
from MACS2.IO.cFixWidthTrack import FWTrackIII

# ------------------------------------
# Main function
# ------------------------------------

class Test_pileup(unittest.TestCase):
    """Unittest for pileup_bdg() in cPileup.pyx.

    """
    def setUp(self):
        self.maxDiff = None
        self.chrom = "chr1"
        self.plus_pos = ( 0, 1, 3, 4, 5 )
        self.minus_pos = ( 5, 6, 8, 9, 10 )
        self.d = 5
        self.scale_factor = 0.5
        self.expect = [ ( 0, 1, 1.0 ),
                        ( 1, 3, 2.0 ),
                        ( 3, 4, 3.0 ),
                        ( 4, 6, 4.0 ),
                        ( 6, 8, 3.0 ),
                        ( 8, 9, 2.0 ),                            
                        ( 9, 10, 1.0 )
                        ]
        self.expect2 = [(0, 1, 13.0),
                        (1, 3, 14.0),
                        (3, 4, 16.0),
                        (4, 6, 18.0),
                        (6, 8, 16.0),
                        (8, 9, 14.0),
                        (9, 10, 13.0)]

        self.d_s = [ 5, 10, 100 ]
        self.scale_factor_s = [ 0.5, 1, 2 ]

    def test_pileup(self):
        # build FWTrackII
        self.fwtrack2 = FWTrackIII()
        for i in self.plus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 0)
        for i in self.minus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 1)            
        self.fwtrack2.finalize()
        
        self.pileup = pileup_bdg(self.fwtrack2, self.d, halfextension=False, scale_factor = self.scale_factor)
        self.result = []
        chrs = self.pileup.get_chr_names()
        for chrom in chrs:
            (p,v) = self.pileup.get_data_by_chr(chrom)
            pnext = iter(p).next
            vnext = iter(v).next
            pre = 0
            for i in xrange(len(p)):
                pos = pnext()
                value = vnext()
                self.result.append( (pre,pos,value) )
                pre = pos
        # check result
        self.assertEqual(self.result, self.expect)

    def test_pileup_w_multiple_d_bdg ( self ):
        # build FWTrackII
        self.fwtrack2 = FWTrackIII()
        for i in self.plus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 0)
        for i in self.minus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 1)            
        self.fwtrack2.finalize()
        # pileup test
        self.pileup = pileup_w_multiple_d_bdg(self.fwtrack2, self.d_s, baseline_value=13, halfextension=False, scale_factor_s = self.scale_factor_s)
        self.result = []
        chrs = self.pileup.get_chr_names()
        for chrom in chrs:
            (p,v) = self.pileup.get_data_by_chr(chrom)
            pnext = iter(p).next
            vnext = iter(v).next
            pre = 0
            for i in xrange(len(p)):
                pos = pnext()
                value = vnext()
                self.result.append( (pre,pos,value) )
                pre = pos
        # check result
        self.assertEqual(self.result, self.expect2)
        

if __name__ == '__main__':
    unittest.main()