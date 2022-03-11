#!/usr/bin/env python
"""Module Description: Test functions for pileup functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest
import numpy as np
from math import log2
from MACS3.Signal.Pileup import *
from MACS3.Signal.PileupV2 import pileup_from_LR, pileup_from_PN

# ------------------------------------
# Main function
# ------------------------------------

class Test_SE_Pileup(unittest.TestCase):
    """Unittest for pileup functions in Pileup.pyx for single-end
    datasets.

    Function to test: se_all_in_one_pileup

    """
    def setUp(self):
        self.plus_pos = np.array(( 0, 1, 3 ), dtype="int32")
        self.minus_pos = np.array(( 8, 9, 10 ), dtype="int32")
        self.rlength = 100      # right end of coordinates
        # expected result from pileup_bdg_se: ( start, end, value )
        # the actual fragment length is 1+five_shift+three_shift
        self.param_1 = { "five_shift": 0,
                         "three_shift": 5,
                         "scale_factor": 0.5,
                         "baseline": 0}
        self.expect_pileup_1 = \
          [ ( 0,  1, 0.5 ),
            ( 1,  3, 1.0 ),
            ( 3,  4, 2.0 ),
            ( 4,  6, 2.5 ),
            ( 6,  8, 2.0 ),
            ( 8,  9, 1.0 ),
            ( 9, 10, 0.5 ) ]
        # expected result from pileup_w_multiple_d_bdg_se: ( start, end, value )
        self.param_2 = { "five_shift": 0,
                         "three_shift": 10,
                         "scale_factor": 2,
                         "baseline": 8}
        self.expect_pileup_2 = \
          [ ( 0,  1, 8.0 ),
            ( 1,  3, 10.0 ),
            ( 3,  8, 12.0 ),
            ( 8,  9, 10.0 ),
            ( 9,  10, 8.0 ),
            ( 10, 11, 8.0 ),
            ( 11, 13, 8.0 ) ]

    def test_pileup_1(self):
        pileup = se_all_in_one_pileup( self.plus_pos, self.minus_pos,
                                       self.param_1["five_shift"],
                                       self.param_1["three_shift"],
                                       self.rlength,
                                       self.param_1["scale_factor"],
                                       self.param_1["baseline"] )
        result = []
        (p,v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append( (pre,pos,value) )
            pre = pos
        # check result
        self.assertEqual( result, self.expect_pileup_1 )

    def test_pileup_2(self):
        pileup = se_all_in_one_pileup( self.plus_pos, self.minus_pos,
                                       self.param_2["five_shift"],
                                       self.param_2["three_shift"],
                                       self.rlength,
                                       self.param_2["scale_factor"],
                                       self.param_2["baseline"] )        
        result = []
        (p,v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append( (pre,pos,value) )
            pre = pos
        # check result
        self.assertEqual( result, self.expect_pileup_2 )

class Test_Quick_Pileup(unittest.TestCase):
    """Unittest for pileup functions in Pileup.pyx for quick-pileup.

    Function to test: quick_pileup

    """
    def setUp(self):
        self.start_pos = np.array( ( 0, 1, 3, 3, 4, 5 ), dtype="int32")
        self.end_pos = np.array( ( 5, 6, 8, 8, 9, 10 ), dtype="int32")
        # expected result from pileup_bdg_se: ( start, end, value )
        self.param_1 = { "scale_factor": 0.5,
                         "baseline": 0}
        self.expect_pileup_1 = \
          [ ( 0,  1, 0.5 ),
            ( 1,  3, 1.0 ),
            ( 3,  4, 2.0 ),
            ( 4,  6, 2.5 ),
            ( 6,  8, 2.0 ),
            ( 8,  9, 1.0 ),
            ( 9, 10, 0.5 ) ]

    def test_pileup_1(self):
        pileup = quick_pileup ( self.start_pos, self.end_pos,
                                self.param_1["scale_factor"],
                                self.param_1["baseline"] )
        result = []
        (p,v) = pileup
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
            pos = pnext()
            value = vnext()
            result.append( (pre,pos,value) )
            pre = pos
        # check result
        self.assertEqual( result, self.expect_pileup_1 )

class Test_Naive_Pileup(unittest.TestCase):
    """Unittest for pileup functions in Pileup.pyx for naive-quick-pileup.

    Function to test: naive_quick_pileup

    """
    def setUp(self):
        self.pos = np.array( ( 2, 3, 5, 5, 6, 7 ), dtype="int32")
        # expected result from pileup_bdg_se: ( start, end, value )
        self.param_1 = { "extension": 2 }
        self.expect_pileup_1 = \
          [ ( 0,  1, 1.0 ),
            ( 1,  3, 2.0 ),
            ( 3,  7, 4.0 ),
            ( 7,  8, 2.0 ),
            ( 8,  9, 1.0 ) ]

    def test_pileup_1(self):
        pileup = naive_quick_pileup ( self.pos,
                                      self.param_1["extension"] )
        result = []
        (p,v) = pileup
        print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
           pos = pnext()
           value = vnext()
           result.append( (pre,pos,value) )
           pre = pos
        #check result
        self.assertEqual( result, self.expect_pileup_1 )

class Test_Over_Two_PV_Array(unittest.TestCase):
    """Unittest for over_two_pv_array function

    Function to test: over_two_pv_array

    """
    def setUp(self):
        self.pv1 = [ np.array( ( 2, 5, 7, 8, 9, 12 ), dtype="int32" ),\
                     np.array( ( 1, 2, 3, 4, 3, 2), dtype="float32" ) ]
        self.pv2 = [ np.array( ( 1, 4, 6, 8, 10, 11 ), dtype="int32" ),\
                     np.array( ( 5, 3, 2, 1, 0, 3), dtype="float32" ) ]
        # expected result from pileup_bdg_se: ( start, end, value )
        self.expect_pv_max = \
          [ ( 0,  1, 5.0 ),
            ( 1,  2, 3.0 ),
            ( 2,  4, 3.0 ),
            ( 4,  5, 2.0 ),
            ( 5,  6, 3.0 ),
            ( 6,  7, 3.0 ),
            ( 7,  8, 4.0 ),
            ( 8,  9, 3.0 ),
            ( 9, 10, 2.0 ),
            (10, 11, 3.0 ) ]
        self.expect_pv_min = \
          [ ( 0,  1, 1.0 ),
            ( 1,  2, 1.0 ),
            ( 2,  4, 2.0 ),
            ( 4,  5, 2.0 ),
            ( 5,  6, 2.0 ),
            ( 6,  7, 1.0 ),
            ( 7,  8, 1.0 ),
            ( 8,  9, 0.0 ),
            ( 9, 10, 0.0 ),
            (10, 11, 2.0 ) ]
        self.expect_pv_mean = \
          [ ( 0,  1, 3.0 ),
            ( 1,  2, 2.0 ),
            ( 2,  4, 2.5 ),
            ( 4,  5, 2.0 ),
            ( 5,  6, 2.5 ),
            ( 6,  7, 2.0 ),
            ( 7,  8, 2.5 ),
            ( 8,  9, 1.5 ),
            ( 9, 10, 1.0 ),
            (10, 11, 2.5 ) ]

    def test_max(self):
        pileup = over_two_pv_array ( self.pv1, self.pv2, func="max" )
        result = []
        (p,v) = pileup
        print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
           pos = pnext()
           value = vnext()
           result.append( (pre,pos,value) )
           pre = pos
        #check result
        self.assertEqual( result, self.expect_pv_max )

    def test_min(self):
        pileup = over_two_pv_array ( self.pv1, self.pv2, func="min" )
        result = []
        (p,v) = pileup
        print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
           pos = pnext()
           value = vnext()
           result.append( (pre,pos,value) )
           pre = pos
        #check result
        self.assertEqual( result, self.expect_pv_min )

    def test_mean(self):
        pileup = over_two_pv_array ( self.pv1, self.pv2, func="mean" )
        result = []
        (p,v) = pileup
        print(p, v)
        pnext = iter(p).__next__
        vnext = iter(v).__next__
        pre = 0
        for i in range(len(p)):
           pos = pnext()
           value = vnext()
           result.append( (pre,pos,value) )
           pre = pos
        #check result
        self.assertEqual( result, self.expect_pv_mean )

class Test_PileupV2_PE(unittest.TestCase):
    """Unittest for pileup functions in PileupV2.pyx.

    Function to test: pileup_from_LR

    """
    def setUp(self):
        self.LR_array1 = np.array( [ (1,5), (2,6), (4,8), (4,8), (5,9), (6,10),
                                     (12,14), (13,17), (14,18), (17,19) ], dtype= [('l','int32'),('r','int32')])
        # expected result from pileup_bdg_se: ( start, end, value )
        self.expect_pileup_1 = np.array( [  ( 1, 0.0 ),
                                            ( 2, 1.0 ),
                                            ( 4, 2.0 ),
                                            ( 8, 4.0 ),
                                            ( 9, 2.0 ),
                                            ( 10, 1.0 ),
                                            ( 12, 0.0 ),
                                            ( 13, 1.0 ),
                                            ( 18, 2.0 ),
                                            ( 19, 1.0 ) ],
                                             dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )
        # with log2(length) as weight
        self.expect_pileup_2 = np.array( [  ( 1, 0.0 ),
                                            ( 2, 2.0 ),
                                            ( 4, 4.0 ),
                                            ( 8, 8.0 ),
                                            ( 9, 4.0 ),
                                            ( 10, 2.0 ),
                                            ( 12, 0.0 ),
                                            ( 13, 1.0 ),
                                            ( 14, 3.0 ),
                                            ( 17, 4.0 ),
                                            ( 18, 3.0 ),
                                            ( 19, 1.0 ) ],
                                             dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )

    def test_pileup_1(self):
        pileup = pileup_from_LR( self.LR_array1 )
        np.testing.assert_equal( pileup, self.expect_pileup_1 )

    def test_pileup_2(self):
        pileup = pileup_from_LR( self.LR_array1, lambda x,y: log2(y-x) )
        np.testing.assert_equal( pileup, self.expect_pileup_2 )

class Test_PileupV2_PE(unittest.TestCase):
    """Unittest for pileup functions in PileupV2.pyx.

    Function to test: pileup_from_PN

    """
    def setUp(self):
        self.P = np.array( ( 0, 1, 3, 3, 4, 5 ), dtype="int32")   #plus strand pos
        self.N = np.array( ( 5, 6, 8, 8, 9, 10 ), dtype="int32")  #minus strand pos
        # expected result from pileup_bdg_se: ( start, end, value )
        self.extsize = 2
        self.expect_pileup_1 = np.array( [  ( 1, 1.0 ),
                                            ( 3, 2.0 ),
                                            ( 4, 4.0 ),
                                            ( 6, 5.0 ),
                                            ( 8, 4.0 ),
                                            ( 9, 2.0 ),
                                            ( 10, 1.0 ) ],
                                            dtype=[ ( 'p', 'uint32' ), ( 'v', 'float32' ) ] )\

    def test_pileup_1(self):
        pileup = pileup_from_PN( self.P, self.N, self.extsize )
        np.testing.assert_equal( pileup, self.expect_pileup_1 )        
