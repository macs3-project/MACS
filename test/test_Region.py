#!/usr/bin/env python
# Time-stamp: <2022-09-14 14:32:12 Tao Liu>

import unittest
import sys

from MACS3.Signal.Region import *

class Test_Regions(unittest.TestCase):
    def setUp( self ):
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
        self.merge_result_bedcontent = "chrY\t0\t200\nchrY\t300\t500\nchrY\t600\t900\nchrY\t1000\t1300\n"

        self.test_regions3 = [(b"chr1",0,100),
                              (b"chr1",300,500),
                              (b"chr1",700,900),
                              (b"chr1",1000,1200),
                              (b"chrY",100,200),
                              (b"chrY",300,400),
                              (b"chrY",600,800),
                              (b"chrY",1200,1300),
                              ]
        self.popped_regions = [ "chr1\t0\t100\nchr1\t300\t500\nchr1\t700\t900\n",
                                "chr1\t1000\t1200\nchrY\t100\t200\nchrY\t300\t400\n",
                                "chrY\t600\t800\nchrY\t1200\t1300\n"]

    def test_add_loc1( self ):
        # make sure the shuffled sequence does not lose any elements
        self.r1 = Regions()
        for a in self.test_regions1:
            self.r1.add_loc(a[0],a[1],a[2])

    def test_add_loc2( self ):
        # make sure the shuffled sequence does not lose any elements
        self.r2 = Regions()
        for a in self.test_regions2:
            self.r2.add_loc(a[0],a[1],a[2])

    def test_merge( self ):
        self.mr = Regions()
        for a in self.test_regions1:
            self.mr.add_loc(a[0],a[1],a[2])
        for a in self.test_regions2:
            self.mr.add_loc(a[0],a[1],a[2])
        self.mr.merge_overlap()
        self.assertEqual( str(self.mr), self.merge_result_bedcontent )

    def test_pop3( self ):
        self.r = Regions()
        for a in self.test_regions3:
            self.r.add_loc(a[0],a[1],a[2])
        # now pop 3 at a time
        ret_list_regions = []
        while self.r.total != 0:
            ret_list_regions.append( str( self.r.pop( 3 ) ) )
        self.assertEqual( ret_list_regions, self.popped_regions )
