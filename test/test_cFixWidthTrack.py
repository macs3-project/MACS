#!/usr/bin/env python
# Time-stamp: <2012-04-24 18:17:48 Tao Liu>

import os
import sys
import unittest

from MACS2.IO.cFixWidthTrack import *

class Test_FWTrackIII(unittest.TestCase):

    def setUp(self):

        self.input_regions = [("chrY",0,0 ),
                              ("chrY",90,0 ),
                              ("chrY",150,0 ),
                              ("chrY",70,0 ),
                              ("chrY",80,0 ),
                              ("chrY",85,0 ),
                              ("chrY",85,0 ),
                              ("chrY",85,0 ),
                              ("chrY",85,0 ),                                    
                              ("chrY",90,1 ),
                              ("chrY",150,1 ),
                              ("chrY",70,1 ),
                              ("chrY",80,1 ),
                              ("chrY",80,1 ),
                              ("chrY",80,1 ),
                              ("chrY",85,1 ),
                              ("chrY",90,1 ),                                    
                              ]
        self.fw = 50

    def test_add_loc(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrackIII(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )         
        self.assertEqual( fw.length(), 17*self.fw )

    def test_filter_dup(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrackIII(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )      
        self.assertEqual( fw.length(), 17*self.fw )

        # filter out more than 3 tags
        fw.filter_dup( 3 )
        # one chrY:85:0 should be removed
        self.assertEqual( fw.total, 16 )

        # filter out more than 2 tags
        fw.filter_dup( 2 )
        # then, one chrY:85:0 and one chrY:80:- should be removed
        self.assertEqual( fw.total, 14 )

        # filter out more than 1 tag
        fw.filter_dup( 1 )
        # then, one chrY:85:0 and one chrY:80:1, one chrY:90:1 should be removed
        self.assertEqual( fw.total, 11 )        

    def test_sample_num(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrackIII(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )         
        self.assertEqual( fw.length(), 17*self.fw )        

        fw.sample_num( 10 )
        self.assertEqual( fw.total, 9 )
        
    def test_sample_num(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrackIII(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )         
        self.assertEqual( fw.length(), 17*self.fw )        

        fw.sample_percent( 0.5 )
        self.assertEqual( fw.total, 8 )        
        
        #fw.print_to_bed()
        #self.assertTrue( abs(result - expect) < 1e-5*result)

        #self.assertEqual(result, expect)

        #self.assertEqual(result, expect)
        #self.assertEqual_float( result, expect )


if __name__ == '__main__':
    unittest.main()
