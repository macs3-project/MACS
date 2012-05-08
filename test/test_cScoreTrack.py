#!/usr/bin/env python
# Time-stamp: <2012-05-07 22:36:20 Tao Liu>

import os
import sys
import unittest

from MACS2.IO.cScoreTrack import *

class Test_ScoreTrackII(unittest.TestCase):

    def setUp(self):
        self.test_regions1 = [("chrY",10,100,10),
                              ("chrY",60,10,10),
                              ("chrY",110,15,20),
                              ("chrY",160,5,20),
                              ("chrY",210,20,5)]
        self.treat_edm = 10
        self.ctrl_edm = 5
        self.p_result = [508, 23, 4, 1, 110]
        self.q_result = [276, 0, 0, 0 ,0]
        self.l_result = [5720, 0, -39, -379, 436]
        self.f_result = [1000, 100, 75, 25, 400]
        self.d_result = [9000, 0, -500, -1500, 1500]
        self.m_result = [1000, 100, 150, 50, 200]
        #
        self.norm_T = [[ 10, 100,  20,   0],
                       [ 60,  10,  20,   0],
                       [110,  15,  40,   0],
                       [160,   5,  40,   0],
                       [210,  20,  10,   0]]
        self.norm_C = [[ 10,  50,  10,   0],
                       [ 60,   5,  10,   0],
                       [110,   7,  20,   0],
                       [160,   2,  20,   0],
                       [210,  10,   5,   0]]
        self.norm_M = [[ 10,  10,   2,   0],
                       [ 60,   1,   2,   0],
                       [110,   1,   4,   0],
                       [160,   0,   4,   0],
                       [210,   2,   1,   0]]
        self.norm_N = [[ 10, 100,  10,   0],  # note precision lost
                       [ 60,  10,  10,   0],
                       [110,  10,  20,   0],
                       [160,   0,  20,   0],
                       [210,  20,   5,   0]] 
        
    def assertEqual_float ( self, a, b, roundn = 5 ):
        self.assertEqual( round( a, roundn ), round( b, roundn ) )

    def test_compute_scores(self):
        s1 = scoreTrackII( self.treat_edm, self.ctrl_edm )
        s1.add_chromosome( "chrY", 5 )
        for a in self.test_regions1:
            s1.add( a[0],a[1],a[2],a[3] )

        s1.change_score_method( ord('p') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r[:,3] == self.p_result).all() )

        s1.change_score_method( ord('q') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r[:,3] == self.q_result).all() )

        s1.change_score_method( ord('l') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r[:,3] == self.l_result).all() )

        s1.change_score_method( ord('f') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r[:,3] == self.f_result).all() )        

        s1.change_score_method( ord('d') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r[:,3] == self.d_result).all() )

        s1.change_score_method( ord('m') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r[:,3] == self.m_result).all() )

    def test_normalize(self):
        s1 = scoreTrackII( self.treat_edm, self.ctrl_edm )
        s1.add_chromosome( "chrY", 5 )
        for a in self.test_regions1:
            s1.add( a[0],a[1],a[2],a[3] )

        s1.change_normalization_method( ord('T') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r == self.norm_T).all() )

        s1.change_normalization_method( ord('C') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r == self.norm_C).all() )

        s1.change_normalization_method( ord('M') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r == self.norm_M).all() )

        s1.change_normalization_method( ord('N') )
        r = s1.get_data_by_chr("chrY")
        self.assertTrue( (r == self.norm_N).all() )

if __name__ == '__main__':
    unittest.main()
