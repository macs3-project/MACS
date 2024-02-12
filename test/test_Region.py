#!/usr/bin/env python
# Time-stamp: <2024-02-12 15:23:48 Tao Liu>

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
        # when we add test_regions1 and test_region2 into the same
        # Regions, we should get this from 'merge_overlaps'
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

        # test_regions4 and test_regions5 are used to test the intersection
        self.test_regions4 = [(b"chrY",0,100),
                              (b"chrY",300,500),
                              (b"chrY",700,900),
                              (b"chrY",1000,1200)
                              ]
        self.test_regions5 = [(b"chrY",100,200),
                              (b"chrY",300,400),
                              (b"chrY",600,800),
                              (b"chrY",1100,1150),
                              (b"chrY",1175,1300)
                              ]
        # After we add test_regions4 and test_region5 into two Regions
        # objects, we should get this from 'intersect'
        self.intersect_regions_4_vs_5 = { b"chrY": [ (300, 400),
                                                     (700, 800),
                                                     (1100,1150),
                                                     (1175,1200) ] }

        
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

    def test_intersect( self ):
        self.r4 = Regions()
        for a in self.test_regions4:
            self.r4.add_loc(a[0],a[1],a[2])        
        self.r5 = Regions()
        for a in self.test_regions5:
            self.r5.add_loc(a[0],a[1],a[2])
        self.intersect_4_vs_5 = self.r4.intersect( self.r5 )
        self.assertEqual( self.intersect_4_vs_5.regions, self.intersect_regions_4_vs_5 )
        self.intersect_4_vs_5 = self.r5.intersect( self.r4 )
        self.assertEqual( self.intersect_4_vs_5.regions, self.intersect_regions_4_vs_5 )        
