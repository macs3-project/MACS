#!/usr/bin/env python
# Time-stamp: <2019-09-20 15:06:44 taoliu>

"""Module Description: Test functions for pileup functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest

from math import log10
from MACS3.Pileup import *
from MACS3.IO.FixWidthTrack import FWTrack

# ------------------------------------
# Main function
# ------------------------------------

class Test_pileup(unittest.TestCase):
    """Unittest for pileup_bdg() in Pileup.pyx.

    """
    def setUp(self):
        self.maxDiff = None
        self.chrom = b"chr1"
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
        self.fwtrack2 = FWTrack()
        for i in self.plus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 0)
        for i in self.minus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 1)            
        self.fwtrack2.finalize()
        
        self.pileup = unified_pileup_bdg(self.fwtrack2, self.d, self.scale_factor, halfextension=False)
        self.result = []
        chrs = self.pileup.get_chr_names()
        for chrom in chrs:
            (p,v) = self.pileup.get_data_by_chr(chrom)
            pnext = iter(p).__next__
            vnext = iter(v).__next__
            pre = 0
            for i in range(len(p)):
                pos = pnext()
                value = vnext()
                self.result.append( (pre,pos,value) )
                pre = pos
        # check result
        self.assertEqual(self.result, self.expect)

    def test_pileup_w_multiple_d_bdg ( self ):
        # build FWTrackII
        self.fwtrack2 = FWTrack(fw=5)
        for i in self.plus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 0)
        for i in self.minus_pos:
            self.fwtrack2.add_loc(self.chrom, i, 1)            
        self.fwtrack2.finalize()
        # pileup test
        self.pileup = unified_pileup_bdg(self.fwtrack2, self.d_s, self.scale_factor_s, baseline_value=13, halfextension=False)
        self.result = []
        chrs = self.pileup.get_chr_names()
        for chrom in chrs:
            (p,v) = self.pileup.get_data_by_chr(chrom)
            pnext = iter(p).__next__
            vnext = iter(v).__next__
            pre = 0
            for i in range(len(p)):
                pos = pnext()
                value = vnext()
                self.result.append( (pre,pos,value) )
                pre = pos
        # check result
        self.assertEqual(self.result, self.expect2)
        
