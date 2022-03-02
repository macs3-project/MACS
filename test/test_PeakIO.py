#!/usr/bin/env python
# Time-stamp: <2022-03-02 18:09:48 Tao Liu>

import unittest
import sys

from MACS3.IO.PeakIO import *

class Test_Region(unittest.TestCase):
    def setUp(self):
        self.test_regions1 = [(b"chrY",0,100),
                              (b"chrY",300,500),
                              (b"chrY",700,900),
                              (b"chrY",1000,1200),
                              ]
        self.test_regions2 = [(b"chrY",100,200),
                              (b"chrY",300,400),
                              (b"chrY",600,800),
                              (b"chrY",1200,1300),
                              ]
        self.merge_result_regions = [ (b"chrY",0,200),
                                      (b"chrY",300,500),
                                      (b"chrY",600,900),
                                      (b"chrY",1000,1300),
                                      ]
        self.subpeak_n = [1,10,100,1000]



    def test_add_loc1(self):
        # make sure the shuffled sequence does not lose any elements
        self.r1 = RegionIO()
        for a in self.test_regions1:
            self.r1.add_loc(a[0],a[1],a[2])

    def test_add_loc2(self):
        # make sure the shuffled sequence does not lose any elements
        self.r2 = RegionIO()
        for a in self.test_regions2:
            self.r2.add_loc(a[0],a[1],a[2])

    def test_merge(self):
        self.mr = RegionIO()
        for a in self.test_regions1:
            self.mr.add_loc(a[0],a[1],a[2])
        for a in self.test_regions2:
            self.mr.add_loc(a[0],a[1],a[2])
        self.mr.merge_overlap()
        self.mr.write_to_bed(sys.stdout)

