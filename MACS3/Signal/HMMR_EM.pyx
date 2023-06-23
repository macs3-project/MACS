# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-04-15 11:24:58 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
from math import sqrt
import logging
import MACS3.Utilities.Logger

logger = logging.getLogger(__name__)
debug   = logger.debug
info    = logger.info
# ------------------------------------
# Other modules
# ------------------------------------

import numpy as np
cimport numpy as np
from cpython cimport bool

# ------------------------------------
# MACS3 modules
# ------------------------------------
#from scipy.stats import norm    # there is another implemented function MACS3.Signal.Prob.pnorm
from MACS3.Signal.Prob import pnorm2

# ------------------------------------
# Misc functions
# ------------------------------------
cdef tuple online_update( float x, long c, float m, float s):
    cdef:
        float delta
    c += 1
    if c == 1:
        return( 1, x, float(0.0) )
    delta = x - m
    m += delta / c
    s += delta * (x - m)
    return (c, m, s)

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
cdef inline float get_weighted_density( int x, float m, float v, w ):
    """Description:
    
    parameters:
      1. x: the observed value
      2. m: the mean of gaussian
      3. v: the variance of the gaussian
      4. w: the weight
    return value:
    """
    return w * pnorm2( float(x), m, v )

cdef int return_greater( np.ndarray data ):
    """
    Return the index of the largest value in an array of doubles
    @param data: an Array of doubles 
    @return an integer representing thre index of the largest value in the inputted array
    """
    cdef:
        int largest_index
        np.ndarray largest_inds
    largest_inds = np.argwhere(data == np.amax(data)) # find the indices of the most likely category(s) (max values)
    if len(largest_inds) > 1:    # if there are multiple "most-likely" categories, ignore data point
        largest_index = -1
    else:
        largest_index = largest_inds[0][0]
    return largest_index
    
# ------------------------------------
# Classes
# ------------------------------------

cdef class HMMR_EM:
    """ Main HMMR EM class.

    This EM trainer will find the optimal mean and stddev of three of
    the four modes -- mono-nucleosomal, di-nucloeosomal, and
    three-nucleosomal fragments. Please note that the mean and stddev
    of the short fragment won't be optimized in this approach and only
    rely on the user's input.
    
    """
    cdef:
        public np.ndarray fragMeans    # fragment length mean for each of the three modes
        public np.ndarray fragStddevs  # fragment length standard deviation for each of the three modes
        np.ndarray fragVars            # we need variances to calculate PDF of normal distribution
        int min_fraglen
        int max_fraglen
        float epsilon          # maximum difference to call the value is converged
        int maxIter             # maximum iternation
        float jump
        int seed                # random seed for downsampling
        bool converged          # wheter the EM is converged
        float sample_percentage
        object __petrack          # PETrackI object
        np.ndarray __data         # data for fragment lengths
        np.ndarray __weights
        
    def __init__ ( self, object petrack, list init_means, list init_stddevs , int min_fraglen = 100, int max_fraglen = 1000, float sample_percentage  = 10, float epsilon = 0.05, int maxIter = 20, float jump = 1.5, int seed = 12345):
        """Initialize HMMR_EM object. The first three parameters are required.

        parameters:
            1. petrack: a MACS3.Signal.PairedEndTrack.PETrackI object
            2. init_means: list of initial means of fragments, for mono, di, and tri signals
            3. init_stddevs: list of initial stddevs of fragments, for mono, di, and tri signals
            4. min_fraglen
            5. max_fraglen
            6. sample_percentage: downsample the original data to get the lengths distribution, default 10
            7. epsilon
            8. maxIter
            9. jump: to emplify the difference
            9. seed
        """
        cdef:
            float cutoff1, cutoff2
            long sum1, sum2, sum3, counter1, counter2, counter3
        # initial values
        self.min_fraglen = min_fraglen
        self.max_fraglen = max_fraglen
        self.epsilon = epsilon
        self.maxIter = maxIter
        self.jump = jump
        self.seed = seed
        self.converged = False
        self.fragMeans = np.array(init_means, dtype=float)
        self.fragStddevs = np.array(init_stddevs, dtype=float)
        self.fragVars = (np.array(init_stddevs, dtype=float)) ** 2
        self.sample_percentage = sample_percentage

        # first, let's prepare the lengths data
        # sample down
        self.__petrack = petrack.sample_percent_copy( self.sample_percentage/100.0, seed = self.seed ) # may need to provide seed option for init function
        self.__data = self.__petrack.fraglengths()
        #fhd = open("allfrag.txt","w")
        #for d in self.__data:
        #    fhd.write(f"{d}\n")
        #fhd.close()
        # then we only keep those with fragment lengths within certain range
        self.__data = self.__data[ np.logical_and( self.__data >= self.min_fraglen, self.__data <= self.max_fraglen ) ]

        info(f"# Downsampled {self.__data.shape[0]} fragments will be used for EM training...")
        # next, we will calculate the weights -- ie the proportion of fragments in each length category
        cutoff1 = (init_means[ 1 ] - init_means[ 0 ])/2 + init_means[ 0 ]
        cutoff2 = (init_means[ 2 ] - init_means[ 1 ])/2 + init_means[ 1 ]

        sum3 = len( self.__data )
        sum2 = sum( self.__data < cutoff2 )
        sum1 = sum( self.__data < cutoff1 )

        counter3 = sum3 - sum2
        counter2 = sum2 - sum1
        counter1 = sum1

        self.__weights = np.array([ counter1/sum3, counter2/sum3, counter3/sum3])
        debug( f"initial: means: {self.fragMeans}, stddevs: {self.fragStddevs}, weights: {self.__weights}" )
        self.__learn()
        return

    cdef bool __learn(self):
        """Description: When we train the mean and stddev for 

        parameters:
        return value:
        """
        cdef:
            int itr = 0         # number of iterations
            int i
            int counter         # number of modes that has been converged
            np.ndarray old_means, old_stddevs, old_weights

        old_means = np.zeros( 3, dtype=float )
        old_stddevs = np.zeros( 3, dtype=float )
        old_weights = np.zeros( 3, dtype=float )
        
        self.converged = False
        while self.converged == False:
            for i in range( 3 ):                
                old_means[i] = self.fragMeans[i]
                old_stddevs[i] = self.fragStddevs[i]
                old_weights[i] = self.__weights[i]

            self.__iterate()
            itr += 1
            debug( f"after iteration {itr}: means: {self.fragMeans}, stddevs: {self.fragStddevs}, weights: {self.__weights}" )
            counter = 0
            for i in range( 3 ):
                if abs(old_means[i] - self.fragMeans[i]) < self.epsilon and abs(old_weights[i] - self.__weights[i]) < self.epsilon and abs(old_stddevs[i] - self.fragStddevs[i]) < self.epsilon:
                    counter += 1
            if counter == 3:
                self.converged = True
                info( f"# Reached convergence after {itr} iterations" )
            if itr >= self.maxIter:
                info( f"# Reached maximum number ({self.maxIter}) of iterations" )
                break
        return self.converged

    cdef void __iterate(self):
        """Description: This is a private function only used by HMMR_EM class

        parameters:
        return value:
        """
        cdef:
            np.ndarray temp, counter, means, variances, s, c
            long total
            int i, j, index

        temp = np.zeros(3, dtype=float) # for each category, the likelihood
        counter = np.zeros(3, dtype=int) # for each category, the number of data points/fragment
        total = 0 # total number of data points/fragments assigned to
                  # three categories
        means = np.zeros(3, dtype=float) # for each category, the new mean
        variances = np.zeros(3, dtype=float) # for each category, the new variances
        s = np.zeros(3, dtype=float)         # for each category, __s
                                             # is for storing
                                             # intermediate values for
                                             # the online algorithm
        c = np.zeros(3, dtype=long)          # for each category, __c
                                             # is for storing
                                             # intermediate values for
                                             # the online algorithm
        for i in range( len( self.__data ) ):
            for j in range( 3 ):
                # for each category: mono, di, tri- (3 in total), we get the likelihoods
                #debug( f"j:{j} {self.__data[i]}, {self.fragMeans[j]}, {self.fragVars[j]}, {self.__weights[j]}" )
                temp[j] = get_weighted_density( self.__data[i], self.fragMeans[j], self.fragVars[j], self.__weights[j] )
            # now look for the most likely category, as `index`
            index = return_greater( temp )

            # then we will update __means and __stds
            if index != -1: # If we can find a mostly likely category
                ##---- update with online algorithm --
                (c[ index ], means[ index ], s[ index ]) = online_update( self.__data[ i ], c[ index ], means[ index ], s[ index ] )
                variances[index] = s[ index ]/c[ index ]
                total += 1
            #debug( f"index: {index} mean: {means[index]}" )

        # Set fragMeans, fragStddevs, and fragVars
        #debug( f"counts: {c[0]} {c[1]} {c[2]}" )
        #debug( f"means: {means[0]} {means[1]} {means[2]}" )
        #debug( f"vars: {variances[0]} {variances[1]} {variances[2]}")
        try:
            for j in range( 3 ):
                if c[ j ] == 0:
                    # if there is no fragment size in this category, ignore
                    continue
                self.fragMeans[ j ] = self.fragMeans[ j ] + self.jump*(means[ j ] - self.fragMeans[ j ])
                self.fragVars[ j ] = self.fragVars[ j ] + self.jump*(variances[ j ] - self.fragVars[ j ])
                self.fragStddevs[ j ] = sqrt( self.fragVars[ j ] )
                self.__weights[ j ] = c[ j ] / total
        except ValueError:
            print(' ValueError:  Adjust --means and --stddev options and re-run command')
        debug( f"After this iteration, {total} fragments have been assigned with either of the three modes" )
        return
