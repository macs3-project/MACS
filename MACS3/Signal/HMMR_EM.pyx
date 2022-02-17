# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-02-17 11:21:42 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
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
from MACS3.Signal.Prob import pnorm

# ------------------------------------
# Misc functions
# ------------------------------------

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
cdef inline float get_weighted_density( x, mean, stddev, weight ):
    """Description:
    
    parameters:
    return value:
    """
    return weight * pnorm( x, mean, stddev )
    # return np.multiply(weight, normal_dist_density) # make sure both are np.array types

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
    if len(largest_inds) > 1: # if there are multiple "most-likely" categories, ignore data point
        largest_index = -1
    else:
        largest_index = largest_inds[0][0]
    return largest_index
    
 # cdef int return_greater( list data ):
 #    """
 #    Return the index of the largest value in an array of doubles
 #    @param data: an Array of doubles 
 #    @return an integer representing thre index of the largest value in the inputted array
 #    """
 #    cdef:
 #        int largest_index
 #        float greatest_value = -1.0
 #        int i
 #    # TL: there must be a easy function in numpy/scipy or even Python itself for this
 #    # also, this function doesn't use any 'self', so could be a general function outside of this class
 #    largest_index = -1
 #    greatest_value = -1.0
 #    for i in range(0, len(data)):
 #        if data[i] > greatest_value:
 #            greatest_value = data[i]
 #            largest_index = i
 #    for i in range(0, len(data)): # if there are other mostly likely category, ignore this data point
 #        if i != largest_index:
 #            if data[i] == greatest_value:
 #                largest_index = -1
 #    return largest_index

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
        # define the values that should be returned
        np.ndarray fragMeans    # fragment length mean for each of the three modes
        np.ndarray fragStddevs  # fragment length standard deviation for each of the three modes
        int min_fraglen
        int max_fraglen
        float episilon          # maximum difference to call the value is converged
        int maxIter             # maximum iternation
        float jump              # amplify the difference
        int seed                # random seed for downsampling
        bool converged          # wheter the EM is converged

        # private variables
        object __petrack          # PETrackI object
        np.ndarray __data         # data for fragment lengths
        np.ndarray __weight
        
    def __init__ ( self, object petrack, list init_means, list init_stddevs , int min_fraglen = 100, int max_fraglen = 1000, float sample_percentage  = 10, float epsilon = 0.0005, int maxIter = 20, float jump = 1.5, int seed = 12345):
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
            9. jump
            10. seed
        """
        cdef:
            float cutoff1, cutoff2
            long sum1, sum2, sum3, counter1, counter2, counter3
        # initial values
        self.petrack = petrack # we may need to use a deepcopy
        self.min_fraglen = min_fraglen
        self.max_fraglen = max_fraglen
        self.epsilon = epsilon
        self.maxIter = maxIter
        self.jump = jump
        self.seed = seed
        self.converged = False
        self.fragMeans = np.array(init_means)
        self.fragStddevs = np.array(init_stddevs)
        self.sample_percentage = sample_percentage

        # first, let's prepare the lengths data
        # sample down
        self.petrack.sample_percent( self.sample_percentage, seed = self.seed ) # may need to provide seed option for init function
        self.__data = self.petrack.fraglengths()
        # then we only keep those with fragment lengths within certain range
        self.__data = self.__data[ np.logical_and( self.__data >= self.min_fraglen, self.__data <= self.max_fraglen ) ]

        # next, we will calculate the weights -- ie the proportion of fragments in each length category
        cutoff1 = init_means[ 1 ] - init_means[ 0 ]/2 + init_means[ 0 ]
        cutoff2 = init_means[ 2 ] - init_means[ 1 ]/2 + init_means[ 1 ]

        sum3 = len( self.__data )
        sum2 = sum( self.__data < cutoff2 )
        sum1 = sum( self.__data < cutoff1 )

        counter3 = sum3 - sum2
        counter2 = sum2 - sum1
        counter1 = sum1

        self.__weight = np.array([ counter1/sum3, counter2/sum3, counter3/sum3])
        self.learn()
        return

    cdef bool learn(self):
        """Description: When we train the mean and stddev for 

        parameters:
        return value:
        """
        cdef:
            int itr = 0         # number of iterations
            int i
            int counter         # number of modes that has been converged
            np.ndarray old_means, old_stddevs, old_weights
            
        old_means = np.array( self.fragMeans )
        old_stddevs = np.array( self.fragStddevs )
        old_weights = np.array( self.__weights )
        
        self.converged = False
        while self.converged == False:
            for i in range( 0, 3 ):                
                old_means[i] = self.fragMeans[i]
                old_stddevs[i] = self.fragStddevs[i]
                old_weights[i] = self.__weights[i]

            self.__iterate()

            counter = 0
            for i in range( 0, 3 ):
                if abs(old_means[i] - self.fragMeans[i]) < self.epsilon and abs(old_weights[i] - self.__weights[i]) < self.epsilon and abs(old_stddevs[i] - self.fragStddevs[i]) < self.epsilon:
                    counter += 1
            if counter == 3:
                self.converged = True
            itr += 1
            if itr >= self.maxIter:
                break
        return self.converged

    cdef void __iterate(self):
        """Description: This is a private function only used by HMMR_EM class

        parameters:
        return value:
        """
        cdef:
            np.ndarray temp, counter, means, stds
            long total
            int i, j, index

        temp = np.zeros(3, dtype=float) # for each category, the likelihood
        counter = np.zeros(3, dtype=int) # for each category, the number of data points/fragment
        total = 0 # total number of data points/fragments assigned to
                  # three categories
        means = np.zeros(3, dtype=float)              # for each category, the new mean
        stds = np.zeros(3, dtype=float)               # for each category, the new stddev
        for i in range( 0, len( self.__data ) ):
            for j in range( 0, 3 ):
                # for each category: mono, di, tri- (4 in total)
                temp[j] = get_weighted_density( self.__data[i], self.fragMeans[j], self.fragStddevs[j], self.__weights[j] )
            # now look for the most likely category, as `index`
            index = return_greater( temp )

            # is this too large of a file to save means,stds for each iteration?
            if index != -1: # If we can find a mostly likely category
                ##----
                means[index] = np.mean( self.__data[0:i] ) #check - do we want means[index] or means[i] or means.append()
                stds[index] = np.std(self.__data[0:i])
                ##---- This block is wrong. We need to implement incremental way to calculate mean and stddev
                ## https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
                ## Check this
                counter[index] += 1.0 # check on this - do we want it to be a list with increasing values?
                total += 1
                
        for j in range( 0, 3 ): # what are the lengths of means, stds, data, mu?
            # we will emplify the difference between new and old means/stds, and update weights.
            self.fragMeans[ j ] = self.fragMeans[ j ] + (self.jump * (means[ j ] - self.fragMeans[ j ]))
            self.fragStddevs[ j ] = self.fragStddevs[ j ] + (self.jump * (stds[ j ] - self.fragStddevs[ j ]))
            self.weights[ j ] = counter[ j ] / total
