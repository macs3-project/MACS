#!/usr/bin/env python
# Time-stamp: <2019-08-09 14:19:31 taoliu>

import unittest

from MACS3.IO.FixWidthTrack import *

class Test_FWTrack(unittest.TestCase):

    def setUp(self):

        self.input_regions = [(b"chrY",0,0 ),
                              (b"chrY",90,0 ),
                              (b"chrY",150,0 ),
                              (b"chrY",70,0 ),
                              (b"chrY",80,0 ),
                              (b"chrY",85,0 ),
                              (b"chrY",85,0 ),
                              (b"chrY",85,0 ),
                              (b"chrY",85,0 ),                                    
                              (b"chrY",90,1 ),
                              (b"chrY",150,1 ),
                              (b"chrY",70,1 ),
                              (b"chrY",80,1 ),
                              (b"chrY",80,1 ),
                              (b"chrY",80,1 ),
                              (b"chrY",85,1 ),
                              (b"chrY",90,1 ),                                    
                              ]
        self.fw = 50

    def test_add_loc(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )         
        self.assertEqual( fw.length, 17*self.fw )

    def test_filter_dup(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )      
        self.assertEqual( fw.length, 17*self.fw )

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
        fw = FWTrack(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )         
        self.assertEqual( fw.length, 17*self.fw )        

        fw.sample_num( 10 )
        self.assertEqual( fw.total, 9 )
        
    def test_sample_percent(self):
        # make sure the shuffled sequence does not lose any elements
        fw = FWTrack(fw=self.fw)
        for ( c, p, s ) in self.input_regions:
            fw.add_loc(c, p, s)
        fw.finalize()
        # roughly check the numbers...
        self.assertEqual( fw.total, 17 )         
        self.assertEqual( fw.length, 17*self.fw )        

        fw.sample_percent( 0.5 )
        self.assertEqual( fw.total, 8 )        

