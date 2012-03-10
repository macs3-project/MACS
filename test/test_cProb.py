#!/usr/bin/env python
# Time-stamp: <2012-03-09 09:50:12 Tao Liu>

"""Module Description: Test functions to calculate probabilities.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""


import os
import sys
import unittest

from math import log10
from MACS2.cProb import *

# ------------------------------------
# Main function
# ------------------------------------

class Test_factorial(unittest.TestCase):

    def setUp(self):
        self.n1 = 100
        self.n2 = 10
        self.n3 = 1

    def test_factorial_big_n1(self):
        expect = 9.332622e+157
        result = factorial(self.n1)
        self.assertTrue( abs(result - expect) < 1e-5*result)

    def test_factorial_median_n2(self):
        expect = 3628800
        result = factorial(self.n2)
        self.assertEqual(result, expect)

    def test_factorial_small_n3(self):
        expect = 1
        result = factorial(self.n3)
        self.assertEqual(result, expect)

class Test_poisson_cdf(unittest.TestCase):

    def setUp(self):
        # n, lam
        self.n1 = (80,100)
        self.n2 = (200,100)
        self.n3 = (100,1000)
        self.n4 = (1500,1000)

    def test_poisson_cdf_n1(self):
        expect = (round(0.9773508,5),round(0.02264918,5))
        result = (round(poisson_cdf(self.n1[0],self.n1[1],False),5),
                  round(poisson_cdf(self.n1[0],self.n1[1],True),5))
        self.assertEqual( result, expect )

    def test_poisson_cdf_n2(self):
        expect = (round(log10(4.626179e-19),4),
                  round(log10(1),4))
        result = (round(log10(poisson_cdf(self.n2[0],self.n2[1],False)),4),
                  round(log10(poisson_cdf(self.n2[0],self.n2[1],True)),4))
        self.assertEqual( result, expect )

    def test_poisson_cdf_n3(self):
        expect = (round(log10(1),4),
                  round(log10(6.042525e-293),4))
        result = (round(log10(poisson_cdf(self.n3[0],self.n3[1],False)),4),
                  round(log10(poisson_cdf(self.n3[0],self.n3[1],True)),4))
        self.assertEqual( result, expect )

    def test_poisson_cdf_n4(self):
        expect = (round(log10(2.097225e-49),4),
                  round(log10(1),4))
        result = (round(log10(poisson_cdf(self.n4[0],self.n4[1],False)),4),
                  round(log10(poisson_cdf(self.n4[0],self.n4[1],True)),4))
        self.assertEqual( result, expect )

class Test_binomial_cdf(unittest.TestCase):

    def setUp(self):
        # x, a, b
        self.n1 = (20,1000,0.01)
        self.n2 = (200,1000,0.01)

    def test_binomial_cdf_n1(self):
        expect = (round(0.001496482,5),round(0.9985035,5))
        result = (round(binomial_cdf(self.n1[0],self.n1[1],self.n1[2],False),5),
                  round(binomial_cdf(self.n1[0],self.n1[1],self.n1[2],True),5))
        self.assertEqual( result, expect )

    def test_binomial_cdf_n2(self):
        expect = (round(log10(8.928717e-190),4),
                  round(log10(1),4))
        result = (round(log10(binomial_cdf(self.n2[0],self.n2[1],self.n2[2],False)),4),
                  round(log10(binomial_cdf(self.n2[0],self.n2[1],self.n2[2],True)),4))
        self.assertEqual( result, expect )

class Test_binomial_cdf_inv(unittest.TestCase):

    def setUp(self):
        # x, a, b
        self.n1 = (0.1,1000,0.01)
        self.n2 = (0.01,1000,0.01)

    def test_binomial_cdf_inv_n1(self):
        expect = 6
        result = binomial_cdf_inv(self.n1[0],self.n1[1],self.n1[2])
        self.assertEqual( result, expect )

    def test_poisson_cdf_inv_n2(self):
        expect = 3
        result = binomial_cdf_inv(self.n2[0],self.n2[1],self.n2[2])
        self.assertEqual( result, expect )

class Test_binomial_pdf(unittest.TestCase):

    def setUp(self):
        # x, a, b
        self.n1 = (20,1000,0.01)
        self.n2 = (200,1000,0.01)

    def test_binomial_cdf_inv_n1(self):
        expect = round(0.001791878,5)
        result = round(binomial_pdf(self.n1[0],self.n1[1],self.n1[2]),5)
        self.assertEqual( result, expect )

    def test_poisson_cdf_inv_n2(self):
        expect = round(log10(2.132196e-188),4)
        result = binomial_pdf(self.n2[0],self.n2[1],self.n2[2])
        result = round(log10(result),4)
        self.assertEqual( result, expect )

if __name__ == '__main__':
    unittest.main()
