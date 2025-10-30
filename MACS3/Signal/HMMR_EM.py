# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-07-24 15:46:59 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
from math import sqrt
import cython
import numpy as np
from cython.cimports import numpy as cnp
from MACS3.Signal.Prob import pnorm2
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.Utilities.Logger import logging

logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.cfunc
def online_update(x: cython.float, c: cython.long,
                  m: cython.float, s: cython.float) -> cython.tuple:
    delta: cython.float

    c += 1
    if c == 1:
        return (1, x, float(0.0))
    delta = x - m
    m += delta / c
    s += delta * (x - m)
    return (c, m, s)


# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html


@cython.inline
@cython.cfunc
def get_weighted_density(x: cython.int, m: cython.float, v:
                         cython.float, w: cython.float) -> cython.float:
    """Description:

    parameters:
      1. x: the observed value
      2. m: the mean of gaussian
      3. v: the variance of the gaussian
      4. w: the weight
    return value:
    """
    return w * pnorm2(float(x), m, v)


@cython.cfunc
def return_greater(data: cnp.ndarray) -> cython.int:
    """
    Return the index of the largest value in an array of doubles
    @param data: an Array of doubles
    @return an integer representing thre index of the largest value in the inputted array
    """

    largest_index: cython.int
    largest_inds: cnp.ndarray
    largest_inds = np.argwhere(data == np.amax(data))  # find the indices of the most likely category(s) (max values)
    if len(largest_inds) > 1:    # if there are multiple "most-likely" categories, ignore data point
        largest_index = -1
    else:
        largest_index = largest_inds[0][0]
    return largest_index

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class HMMR_EM:
    """ Main HMMR EM class.

    This EM trainer will find the optimal mean and stddev of three of
    the four modes -- mono-nucleosomal, di-nucloeosomal, and
    three-nucleosomal fragments. Please note that the mean and stddev
    of the short fragment won't be optimized in this approach and only
    rely on the user's input.

    """
    # fragment length mean for each of the three modes
    fragMeans = cython.declare(cnp.ndarray, visibility='public')
    # fragment length standard deviation for each of the three modes
    fragStddevs = cython.declare(cnp.ndarray, visibility='public')
    fragVars: cnp.ndarray  # we need variances to calculate PDF of normal distribution
    min_fraglen: cython.int
    max_fraglen: cython.int
    epsilon: cython.float  # maximum difference to call the value is converged
    maxIter: cython.int  # maximum iteration
    jump: cython.float
    seed: cython.int            # random seed for downsampling
    converged: cython.bint      # wheter the EM is converged
    sample_percentage: cython.float
    __petrack: PETrackI         # PETrackI object
    __data: cnp.ndarray          # data for fragment lengths
    __weights: cnp.ndarray

    def __init__(self, petrack: PETrackI, init_means: cython.list,
                 init_stddevs: cython.list,
                 min_fraglen: cython.int = 100,
                 max_fraglen: cython.int = 1000,
                 sample_percentage: cython.float = 10,
                 epsilon: cython.float = 0.05,
                 maxIter: cython.int = 20,
                 jump: cython.float = 1.5, seed: cython.int = 12345):
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
        cutoff1: cython.float
        cutoff2: cython.float
        sum1: cython.long
        sum2: cython.long
        sum3: cython.long
        counter1: cython.long
        counter2: cython.long
        counter3: cython.long

        # initial values
        self.min_fraglen = min_fraglen
        self.max_fraglen = max_fraglen
        self.epsilon = epsilon
        self.maxIter = maxIter
        self.jump = jump
        self.seed = seed
        self.converged = False
        self.fragMeans = np.array(init_means, dtype='f4')
        self.fragStddevs = np.array(init_stddevs, dtype='f4')
        self.fragVars = (np.array(init_stddevs, dtype='f4')) ** 2
        self.sample_percentage = sample_percentage

        # first, let's prepare the lengths data
        # sample down
        self.__petrack = petrack.sample_percent_copy(self.sample_percentage/100.0, seed=self.seed)  # may need to provide seed option for init function
        self.__data = self.__petrack.fraglengths()
        # fhd = open("allfrag.txt","w")
        # for d in self.__data:
        #    fhd.write(f"{d}\n")
        # fhd.close()
        # then we only keep those with fragment lengths within certain range
        self.__data = self.__data[np.logical_and(self.__data >= self.min_fraglen, self.__data <= self.max_fraglen)]

        info(f"# Downsampled {self.__data.shape[0]} fragments will be used for EM training...")
        # next, we will calculate the weights -- ie the proportion of fragments in each length category
        cutoff1 = (init_means[1] - init_means[0])/2 + init_means[0]
        cutoff2 = (init_means[2] - init_means[1])/2 + init_means[1]

        sum3 = len(self.__data)
        sum2 = sum(self.__data < cutoff2)
        sum1 = sum(self.__data < cutoff1)

        counter3 = sum3 - sum2
        counter2 = sum2 - sum1
        counter1 = sum1

        self.__weights = np.array([counter1/sum3, counter2/sum3, counter3/sum3])
        debug(f"initial: means: {self.fragMeans}, stddevs: {self.fragStddevs}, weights: {self.__weights}")
        self.__learn()
        return

    @cython.cfunc
    def __learn(self) -> cython.bint:
        """Description: When we train the mean and stddev

        parameters:
        return value:
        """
        itr: cython.int = 0
        i: cython.int
        counter: cython.int
        old_means: cnp.ndarray = np.zeros(3, dtype='f4')
        old_stddevs: cnp.ndarray = np.zeros(3, dtype='f4')
        old_weights: cnp.ndarray = np.zeros(3, dtype='f4')

        self.converged = False
        while self.converged is False:
            for i in range(3):       
                old_means[i] = self.fragMeans[i]
                old_stddevs[i] = self.fragStddevs[i]
                old_weights[i] = self.__weights[i]

            self.__iterate()
            itr += 1
            debug(f"after iteration {itr}: means: {self.fragMeans}, stddevs: {self.fragStddevs}, weights: {self.__weights}")
            counter = 0
            for i in range(3):
                if abs(old_means[i] - self.fragMeans[i]) < self.epsilon and abs(old_weights[i] - self.__weights[i]) < self.epsilon and abs(old_stddevs[i] - self.fragStddevs[i]) < self.epsilon:
                    counter += 1
            if counter == 3:
                self.converged = True
                info(f"# Reached convergence after {itr} iterations")
            if itr >= self.maxIter:
                info(f"# Reached maximum number ({self.maxIter}) of iterations")
                break
        return self.converged

    @cython.cfunc
    def __iterate(self):
        """Description: This is a private function only used by HMMR_EM class

        parameters:
        return value:
        """
        temp: cnp.ndarray
        means: cnp.ndarray
        variances: cnp.ndarray
        s: cnp.ndarray
        c: cnp.ndarray
        total: cython.long
        i: cython.int
        j: cython.int
        index: cython.int

        temp = np.zeros(3, dtype='f4')  # for each category, the likelihood
        total = 0  # total number of data points/fragments assigned to three categories
        means = np.zeros(3, dtype='f4')  # for each category, the new mean
        variances = np.zeros(3, dtype='f4')  # for each category, the new variances
        s = np.zeros(3, dtype='f4')          # for each category, __s
                                            # is for storing
                                            # intermediate values
                                            # for the online
                                            # algorithm
        c = np.zeros(3, dtype='i8')          # for each category, __c
                                            # is for storing
                                            # intermediate values for
                                            # the online algorithm
        for i in range(len(self.__data)):
            for j in range(3):
                # for each category: mono, di, tri- (3 in total), we get the likelihoods
                # debug( f"j:{j} {self.__data[i]}, {self.fragMeans[j]}, {self.fragVars[j]}, {self.__weights[j]}" )
                temp[j] = get_weighted_density(self.__data[i], self.fragMeans[j], self.fragVars[j], self.__weights[j])
            # now look for the most likely category, as `index`
            index = return_greater(temp)

            # then we will update __means and __stds
            if index != -1:  # If we can find a mostly likely category
                # #---- update with online algorithm --
                (c[index], means[index], s[index]) = online_update(self.__data[i], c[index], means[index], s[index])
                variances[index] = s[index]/c[index]
                total += 1
            # debug( f"index: {index} mean: {means[index]}" )

        # Set fragMeans, fragStddevs, and fragVars
        # debug( f"counts: {c[0]} {c[1]} {c[2]}" )
        # debug( f"means: {means[0]} {means[1]} {means[2]}" )
        # debug( f"vars: {variances[0]} {variances[1]} {variances[2]}")
        try:
            for j in range(3):
                if c[j] == 0:
                    # if there is no fragment size in this category, ignore
                    continue
                self.fragMeans[j] = self.fragMeans[j] + self.jump*(means[j] - self.fragMeans[j])
                self.fragVars[j] = self.fragVars[j] + self.jump*(variances[j] - self.fragVars[j])
                self.fragStddevs[j] = sqrt(self.fragVars[j])
                self.__weights[j] = c[j] / total
        except ValueError:
            print(' ValueError:  Adjust --means and --stddevs options and re-run command')
        debug(f"After this iteration, {total} fragments have been assigned with either of the three modes")
        return
