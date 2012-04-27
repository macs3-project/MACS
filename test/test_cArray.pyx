#!/usr/bin/env python
# Time-stamp: <2012-04-22 23:13:41 Tao Liu>

"""Module Description: Test functions for pileup functions.

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
from random import randint

from MACS2.cArray import IntArray

# ------------------------------------
# Main function
# ------------------------------------

class Test_cArray(unittest.TestCase):
    """Unittest for pileup_bdg() in cPileup.pyx.

    """
    def setUp(self):
        self.l = 1000000
        self.d_s = [ 5, 10, 100 ]
        self.scale_factor_s = [ 0.5, 1, 2 ]

    def test_IntArray(self):
        self.a = IntArray(self.l)
        for i in range(1000000):
            self.a.put(randint(0,1000000))
        self.a.sort()
        ai = iter(self.a)
        ai.next()

if __name__ == '__main__':
    unittest.main()
