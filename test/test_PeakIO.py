#!/usr/bin/env python
# Time-stamp: <2022-03-03 14:44:54 Tao Liu>

import unittest
import sys

from MACS3.IO.PeakIO import *

class Test_RegionIO(unittest.TestCase):
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

class Test_PeakIO(unittest.TestCase):
    def setUp(self):
        self.test_peaks1 = [  (b"chrY",0,100),
                              (b"chrY",300,500),
                              (b"chrY",700,900),
                              (b"chrY",1000,1200),
                              (b"chrY",1250,1450),
                              (b"chrX",1000,2000),
                              (b"chrX",3000,4000),
                              (b"chr1",100, 10000), # chr1 only has one region but overlapping with peaks2
                              (b"chr2",1000,2000), # only peaks1 has chr2
                              (b"chr4",500,800),   # chr4 only one region, and not overlapping with peaks2
                            ]
        self.test_peaks2 = [  (b"chrY",100,200),
                              (b"chrY",300,400),
                              (b"chrY",600,800),
                              (b"chrY",1100,1300),
                              (b"chrY",1700,1800),
                              (b"chrX",1100,1200),
                              (b"chrX",1300,1400),
                              (b"chr1",2000,3000),
                              (b"chr3",1000, 5000), # only peaks2 has chr3
                              (b"chr4",1000,2000),
                            ]
        self.result_exclude2from1 = [ (b"chrY",0,100),
                                      (b"chrX",3000,4000),
                                      (b"chr2",1000,2000),
                                      (b"chr4",500,800),
                                     ]
        self.exclude2from1 = PeakIO()
        for a in self.result_exclude2from1:
            self.exclude2from1.add(a[0],a[1],a[2])

    def test_exclude(self):
        r1 = PeakIO()
        for a in self.test_peaks1:
            r1.add(a[0],a[1],a[2])
        r2 = PeakIO()
        for a in self.test_peaks2:
            r2.add(a[0],a[1],a[2])
        r1.exclude(r2)
        result = str(r1)
        expected = str(self.exclude2from1)
        print( "result:\n",result )
        print( "expected:\n", expected )
        self.assertEqual( result, expected )

