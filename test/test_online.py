#!/usr/bin/env python
"""Module Description: Test functions for online algorithm functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest
from math import sqrt
# ------------------------------------
# Main function
# ------------------------------------


def update_m_s(x, c, m, s):
    c += 1
    delta = x - m
    m += delta / c
    s += delta * (x - m)
    return (c, m, s)


class Test_Mean(unittest.TestCase):
    def setUp(self):
        # np.random.randint(1,1000,size=10)
        self.data = [410, 710, 891, 312, 453, 460, 237, 622, 268]
        # [np.mean(self.data[0:x]) for x in range(0,len(self.data)+1)]
        self.correct_mean = [410.0, 560.0, 670.3333333333334, 580.75, 555.2, 539.3333333333334, 496.14285714285717, 511.875, 484.77777777777777]
        # [np.std(self.data[0:x]) for x in range(0,len(self.data)+1)]
        self.correct_stddev = [0.0, 150.0, 198.36050234078579, 231.48582569997671, 213.25984150795946, 197.8852080261573, 211.55845431425504, 202.2247743848414, 205.48737608036913]

    def test_mean(self):
        c = 1
        m = self.data[0]
        s = 0
        for i in range(1, len(self.data)):
            (c, m, s) = update_m_s(self.data[i], c, m, s)
            mean = m
            std = sqrt(s/c)
            self.assertAlmostEqual(self.correct_mean[i], mean, 2)
            self.assertAlmostEqual(self.correct_stddev[i], std, 2)
