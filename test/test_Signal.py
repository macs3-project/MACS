#!/usr/bin/env python
# Time-stamp: <2019-12-18 17:02:57 taoliu>

"""Module Description: Test functions for Signal.pyx

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest

from math import log10
import numpy as np
from MACS2.Signal import maxima, savitzky_golay, savitzky_golay_order2_deriv1

# ------------------------------------
# Main function
# ------------------------------------

class Test_maxima(unittest.TestCase):

    def setUp(self):
        data = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 
                 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
                 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
                 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
                 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 
                 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
                 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 
                 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
                 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 
                 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
                 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 
                 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
                 3, 3, 3, 3, 3, 3 ]
        self.signal = np.array(data, dtype=np.float32 )
        self.windowsize = 253
        self.summit = 161       # this is based on 1-deriv smoothed data
        self.smoothed162 = -2.98155597e-18

    def test_implement_smooth_here ( self ):
        signal = self.signal
        window_size = self.windowsize
        half_window = (window_size - 1) // 2
        # precompute coefficients
        b = np.array([[1, k, k**2] for k in range(-half_window, half_window+1)], dtype='int64')
        m = np.linalg.pinv(b)[1]
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
        lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
        signal = np.concatenate((firstvals, signal, lastvals))
        #print (repr(m))
        #print (repr(signal))
        ret = np.convolve( m[::-1], signal.astype("float64"), mode='valid').astype("float32")
        p = ret[162]
        print ("calculated step by step:\n", p)
        print ("expected:\n", self.smoothed162)
        self.assertAlmostEqual( p, self.smoothed162, places = 4 )
        self.assertEqual( np.sign(p), np.sign(self.smoothed162) )

    def test_implement_smooth_here2 ( self ): # try to tweak some dtypes to see if problem exists
        signal = self.signal
        window_size = self.windowsize
        half_window = (window_size - 1) // 2
        # precompute coefficients
        b = np.array([[1, k, k**2] for k in range(-half_window, half_window+1)], dtype='int32')
        m = np.linalg.pinv(b)[1]
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
        lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
        signal = np.concatenate((firstvals, signal, lastvals))
        #print (repr(m))
        #print (repr(signal))
        ret = np.convolve( m[::-1], signal.astype("float32"), mode='valid')
        p = ret[162]
        print ("calculated step by step:\n", p)
        print ("expected:\n", self.smoothed162)
        self.assertAlmostEqual( p, self.smoothed162, places = 4 )
        self.assertEqual( np.sign(p), np.sign(self.smoothed162) )        

    def test_maxima(self):
        expect = self.summit
        result = maxima( self.signal, self.windowsize )[0]
        self.assertEqual( result, expect, msg=f"Not equal: result: {result}, expected: {expect}" )

    def assertEqual_nparray1d ( self, a, b, roundn = 7 ):
        self.assertEqual( a.shape[0], b.shape[0] )
        l = a.shape[0]
        for i in range( l ):
            self.assertAlmostEqual( a[i], b[i], places = roundn, msg=f"Not equal at {i} {a[i]} {b[i]}" )
            
