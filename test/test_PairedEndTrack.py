#!/usr/bin/env python
# Time-stamp: <2019-08-09 14:19:31 taoliu>

import unittest

from MACS3.IO.PairedEndTrack import *

class Test_PETrackI(unittest.TestCase):

    def setUp(self):

        self.input_regions = [(b"chrY",0,100 ),
                              (b"chrY",70,270 ),
                              (b"chrY",70,100 ),
                              (b"chrY",80,160 ),
                              (b"chrY",80,160 ),
                              (b"chrY",80,180 ),
                              (b"chrY",80,180 ),
                              (b"chrY",85,185 ),
                              (b"chrY",85,285 ),
                              (b"chrY",85,285 ),
                              (b"chrY",85,285 ),
                              (b"chrY",85,385 ),
                              (b"chrY",90,190 ),
                              (b"chrY",90,190 ),
                              (b"chrY",90,191 ),
                              (b"chrY",150,190 ),
                              (b"chrY",150,250 ),
                              ]
        self.t = sum([ x[2]-x[1] for x in self.input_regions ])

    def test_add_loc(self):
        pe = PETrackI()
        for ( c, l, r ) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        # roughly check the numbers...
        self.assertEqual( pe.total, 17 )
        self.assertEqual( pe.length, self.t )

    def test_filter_dup(self):
        pe = PETrackI()
        for ( c, l, r ) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        # roughly check the numbers...
        self.assertEqual( pe.total, 17 )
        self.assertEqual( pe.length, self.t )

        # filter out more than 3 tags
        pe.filter_dup( 3 )
        self.assertEqual( pe.total, 17 )

        # filter out more than 2 tags
        pe.filter_dup( 2 )
        self.assertEqual( pe.total, 16 )

        # filter out more than 1 tag
        pe.filter_dup( 1 )
        self.assertEqual( pe.total, 12 )


    def test_sample_num(self):
        pe = PETrackI()
        for ( c, l, r ) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        pe.sample_num( 10 )
        self.assertEqual( pe.total, 10 )

    def test_sample_percent(self):
        pe = PETrackI()
        for ( c, l, r ) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        pe.sample_percent( 0.5 )
        self.assertEqual( pe.total, 8 )

