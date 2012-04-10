#!/usr/bin/env python
# Time-stamp: <2012-02-26 19:33:56 Tao Liu>

import os
import sys
import unittest

from math import log10
from MACS2.cStat import *
from pymc.Matplot import mean
import numpy.random as rand


class Test_MCMCPoissonPosteriorRatio(unittest.TestCase):
    def setUp(self):
        self.p1 = [10,5]
        self.p2 = [10,5]
        self.p3 = [10,5]
        self.p4 = [10,5]        

    def test_func(self):
        # call func
        rand.seed([10])
        P = MCMCPoissonPosteriorRatio(15000,5000,self.p1[0],self.p1[1])
        P = sorted(P)
        print self.p1[0],self.p1[1],P[100],mean(P),P[-100]

        rand.seed([10])
        P = MCMCPoissonPosteriorRatio(15000,5000,self.p2[0],self.p2[1])
        P = sorted(P)
        print self.p2[0],self.p2[1],P[100],mean(P),P[-100]

        rand.seed([10])
        P = MCMCPoissonPosteriorRatio(15000,5000,self.p3[0],self.p3[1])
        P = sorted(P)
        print self.p3[0],self.p3[1],P[100],mean(P),P[-100]

        rand.seed([10])
        P = MCMCPoissonPosteriorRatio(15000,5000,self.p4[0],self.p4[1])
        P = sorted(P)
        print self.p4[0],self.p4[1],P[100],mean(P),P[-100]        

        #self.assertTrue( abs(result - expect) < 1e-5*result)

if __name__ == '__main__':
    unittest.main()
