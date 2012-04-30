#!/usr/bin/env python
# Time-stamp: <2012-04-29 17:27:36 Tao Liu>

import os
import sys
import unittest

from MACS2.IO.cPeakIO import *

class Test_Region(unittest.TestCase):

    def setUp(self):
        self.test_regions1 = [("chrY",0,100),
                              ("chrY",300,500),
                              ("chrY",700,900),
                              ("chrY",1000,1200),
                              ]
        self.test_regions2 = [("chrY",100,200),
                              ("chrY",300,400),
                              ("chrY",600,800),
                              ("chrY",1200,1300),
                              ]
        self.merge_result_regions = [ ("chrY",0,200),
                                      ("chrY",300,500),
                                      ("chrY",600,900),
                                      ("chrY",1000,1300),
                                      ]
        self.subpeak_n = [1,10,100,1000]



    def test_add_loc1(self):
        # make sure the shuffled sequence does not lose any elements
        self.r1 = Region()
        for a in self.test_regions1:
            self.r1.add_loc(a[0],a[1],a[2])

    def test_add_loc2(self):
        # make sure the shuffled sequence does not lose any elements
        self.r2 = Region()
        for a in self.test_regions2:
            self.r2.add_loc(a[0],a[1],a[2])

    def test_merge(self):
        self.mr = Region()
        for a in self.test_regions1:
            self.mr.add_loc(a[0],a[1],a[2])
        for a in self.test_regions2:
            self.mr.add_loc(a[0],a[1],a[2])            
        self.mr.merge_overlap()
        self.mr.write_to_bed(sys.stdout)

#    def test_subpeak_letters(self):
#        for i in self.subpeak_n:
#            print subpeak_letters(i)

if __name__ == '__main__':
    unittest.main()
