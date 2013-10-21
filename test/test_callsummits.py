#!/usr/bin/env python
# Time-stamp: <2013-10-18 16:27:40 Tao Liu>

import os
import sys
import unittest

from MACS2.IO.cCallPeakUnit import *
import numpy as np
from MACS2.IO.cPeakIO import PeakIO
from math import factorial
from random import normalvariate

class Test_CallSummits ( unittest.TestCase ):
    
    def setUp( self ):
        self.range             = [   0, 2000 ]
        self.binding_sites     = [ 300, 500, 700 ]
        self.binding_strength  = [ 60,  45,  55 ] # approximate binding affility 
        self.binding_width     = [ 150, 150, 150 ]# binding width, left and right sides are cutting sites
        self.cutting_variation = 50              # variation at the correct cutting sites
        self.tag_size          = 50
        self.test_tags_file    = "random_test.bed"
        self.genome_size       = 10000
        
        self.plus_tags         = [ ]
        self.minus_tags        = [ ]
        
        for i in range( len(self.binding_sites) ):
            j = 0
            while j <= self.binding_strength[ i ]:
                x = int( normalvariate( self.binding_sites[ i ] - self.binding_width[ i ]/2,
                                        self.cutting_variation ) )
                if x > self.range[ 0 ] and x + self.tag_size < self.range[ 1 ]:
                    self.plus_tags.append( x )
                    j += 1
            
            j = 0
            while j <= self.binding_strength[ i ]:
                x = int( normalvariate( self.binding_sites[ i ] + self.binding_width[ i ]/2,
                                        self.cutting_variation ) )
                if x - self.tag_size > self.range[ 0 ] and x < self.range[ 1 ]:
                    self.minus_tags.append( x )
                    j += 1

        self.plus_tags = sorted(self.plus_tags)
        self.minus_tags = sorted(self.minus_tags)

        #print self.plus_tags
        #print self.minus_tags

        self.result_peak = PeakIO()

        # write reads in bed files
        fhd = open( self.test_tags_file, "w" )
        for x in self.plus_tags:
            fhd.write( "chr1\t%d\t%d\t.\t0\t+\n" % ( x, x + self.tag_size ) )
        for x in self.minus_tags:
            fhd.write( "chr1\t%d\t%d\t.\t0\t-\n" % ( x - self.tag_size, x ) )

    def test_pileup ( self ):
        pass


    # def test_wo_subpeak ( self ):
    #     peak_content = self.test_peak_content
    #     tsummit = []
    #     summit_pos   = 0
    #     summit_value = 0
    #     for i in range(len(peak_content)):
    #         (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
    #         tscore = ttreat_p #self.pqtable[ get_pscore(int(ttreat_p), tctrl_p) ] # use qscore as general score to find summit
    #         if not summit_value or summit_value < tscore:
    #             tsummit = [(tend + tstart) / 2, ]
    #             tsummit_index = [ i, ]
    #             summit_value = tscore
    #         elif summit_value == tscore:
    #             # remember continuous summit values
    #             tsummit.append(int((tend + tstart) / 2))
    #             tsummit_index.append( i )
    #     # the middle of all highest points in peak region is defined as summit
    #     print "wo, all:",tsummit
    #     midindex = int((len(tsummit) + 1) / 2) - 1
    #     summit_pos    = tsummit[ midindex ]
    #     summit_index  = tsummit_index[ midindex ]
    #     print "wo:",summit_pos

    # def test_w_subpeak ( self ):
    #     peak_content = self.test_peak_content

    #     smoothlen = 20

    #     peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            
    #     # Add 10 bp padding to peak region so that we can get true minima
    #     end = peak_content[ -1 ][ 1 ] + 10
    #     start = peak_content[ 0 ][ 0 ] - 10
    #     if start < 0:
    #         start_boundary = 5 + start # this is the offset of original peak boundary in peakdata list.
    #         start = 0
    #     else:
    #         start_boundary = 5 # this is the offset of original peak boundary in peakdata list.

    #     peakdata = np.zeros(end - start, dtype='float32') # save the scores (qscore) for each position in this region
    #     peakindices = np.zeros(end - start, dtype='int32') # save the indices for each position in this region
    #     for i in range(len(peak_content)):
    #         (tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i]
    #         #tscore = self.pqtable[ get_pscore(int(ttreat_p), tctrl_p) ] # use qscore as general score to find summit
    #         tscore = ttreat_p # use pileup as general score to find summit
    #         m = tstart - start + start_boundary
    #         n = tend - start + start_boundary
    #         peakdata[m:n] = tscore
    #         peakindices[m:n] = i

    #     np.set_printoptions(precision=1, suppress=True)
    #     #print "before smoothed data:", len(peakdata), peakdata
    #     print "maximum points:", np.where(peakdata == peakdata.max())
    #     summit_offsets = maxima(peakdata, smoothlen) # offsets are the indices for summits in peakdata/peakindices array.
    #     print "summit_offsets:", summit_offsets

    #     m = np.searchsorted(summit_offsets, start_boundary)
    #     n = np.searchsorted(summit_offsets, peak_length + start_boundary, 'right')
    #     summit_offsets = summit_offsets[m:n]
    #     print "summit_offsets adjusted:", summit_offsets        
        
    #     summit_offsets = enforce_peakyness(peakdata, summit_offsets)
    #     print "summit_offsets enforced:", summit_offsets        
        
    #     summit_indices = peakindices[summit_offsets] # indices are those point to peak_content
    #     summit_offsets -= start_boundary
    #     print "summit_offsets final:", summit_offsets        

    #     for summit_offset, summit_index in zip(summit_offsets, summit_indices):
    #         print "w:",start+summit_offset

def maxima ( signal, window_size=51 ):
    """return the local maxima in a signal after applying a 2nd order
    Savitsky-Golay (polynomial) filter using window_size specified  
    """
    #data1 = savitzky_golay(signal, window_size, order=2, deriv=1)
    #data2 = savitzky_golay(signal, window_size, order=2, deriv=2)
    #m = np.where(np.diff(np.sign( data1 )) <= -1)[0].astype('int32')

    #data1 = savitzky_golay_order2(signal, window_size, deriv=1)
    data1 = savitzky_golay(signal, window_size, order=2, deriv=1)
    m = np.where( np.diff( np.sign( data1 ) ) <= -1)[0].astype('int32')
    #m = np.where( np.logical_and( data2 < 0 , abs(data1) <= 1e-10) ) [0].astype('int32')
    return m


def savitzky_golay_order2(signal, window_size, deriv=0):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    if window_size % 2 != 1: window_size += 1
    half_window = (window_size - 1) / 2
    # precompute coefficients
    b = np.mat([[1, k, k**2] for k in range(-half_window, half_window+1)],
               dtype='int64')
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
    lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
    signal = np.concatenate((firstvals, signal, lastvals))
    ret = np.convolve( m, signal.astype('float64'), mode='valid').astype('float32')
    return ret

def savitzky_golay(y, window_size, order=2, deriv=0, rate=1):

    if window_size % 2 != 1: window_size += 1

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

if __name__ == '__main__':
    unittest.main()
