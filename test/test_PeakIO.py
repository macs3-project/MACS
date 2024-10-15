#!/usr/bin/env python
# Time-stamp: <2024-10-14 21:32:21 Tao Liu>

import unittest
import sys

from MACS3.IO.PeakIO import *

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
        # print( "result:\n",result )
        # print( "expected:\n", expected )
        self.assertEqual( result, expected )

