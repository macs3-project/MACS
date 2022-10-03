#!/usr/bin/env python
# Time-stamp: <2022-03-01 00:16:44 Tao Liu>

import unittest

from MACS3.Signal.BedGraph import *

class Test_bedGraphTrackI_add_loc(unittest.TestCase):

    def setUp(self):
        self.test_regions1 = [(b"chrY",0,10593155,0.0),
                              (b"chrY",10593155,10597655,0.0066254580149)]

    def test_add_loc1(self):
        # make sure the shuffled sequence does not lose any elements
        bdg = bedGraphTrackI()
        for a  in self.test_regions1:
            bdg.add_loc(a[0],a[1],a[2],a[3])

        #self.assertTrue( abs(result - expect) < 1e-5*result)

        #self.assertEqual(result, expect)

        #self.assertEqual(result, expect)

class Test_bedGraphTrackI_overlie(unittest.TestCase):

    def setUp(self):
        self.test_regions1 = [(b"chrY",0,70,0.00),
                              (b"chrY",70,80,7.00),
                              (b"chrY",80,150,0.00)]
        self.test_regions2 = [(b"chrY",0,85,0.00),
                              (b"chrY",85,90,75.00),
                              (b"chrY",90,155,0.00)]
        self.test_regions3 = [(b"chrY",0,75,20.0),
                              (b"chrY",75,90,35.0),
                              (b"chrY",90,150,10.0)]
        self.test_overlie_max = [(b"chrY",  0, 75, 20.0), # max
                                 (b"chrY", 75, 85, 35.0),
                                 (b"chrY", 85, 90, 75.0),
                                 (b"chrY", 90,150, 10.0),
                                 (b"chrY",150,155, 0.0)]
        self.test_overlie_mean = [(b"chrY",  0, 70, 6.66667), # mean
                                  (b"chrY", 70, 75, 9),
                                  (b"chrY", 75, 80, 14),
                                  (b"chrY", 80, 85, 11.66667),
                                  (b"chrY", 85, 90, 36.66667),
                                  (b"chrY", 90,150, 3.33333),
                                  (b"chrY",150,155, 0.0)]
        self.test_overlie_fisher = [(b"chrY",  0, 70, (0,0,20), 92.10340371976183, 1.1074313239555578e-17, 16.9557 ), # fisher
                                    (b"chrY", 70, 75, (7,0,20), 124.33959502167846, 1.9957116587802055e-24, 23.6999 ),
                                    (b"chrY", 75, 80, (7,0,35), 193.41714781149986, 4.773982707347631e-39, 38.3211 ),
                                    (b"chrY", 80, 85, (0,0,35), 161.1809565095832, 3.329003070922764e-32, 31.4777),
                                    (b"chrY", 85, 90, (0,75,35), 506.56872045869005, 3.233076792862357e-106, 105.4904),
                                    (b"chrY", 90,150, (0,0,10), 46.051701859880914, 2.8912075645386016e-08, 7.5389),
                                    (b"chrY",150,155, (0,0,0), 0, 1.0, 0.0)]
        self.test_sum = []
        self.test_product = []
        self.test_subtract = []
        self.test_divide = []

        self.bdg1 = bedGraphTrackI()
        self.bdg2 = bedGraphTrackI()
        self.bdg3 = bedGraphTrackI()
        for a in self.test_regions1:
            self.bdg1.add_loc(a[0],a[1],a[2],a[3])

        for a in self.test_regions2:
            self.bdg2.add_loc(a[0],a[1],a[2],a[3])

        for a in self.test_regions3:
            self.bdg3.add_loc(a[0],a[1],a[2],a[3])

    def assertEqual_float ( self, a, b, roundn = 4 ):
        self.assertEqual( round( a, roundn ), round( b, roundn ) )

    def test_overlie_max(self):
        bdgb = self.bdg1.overlie([self.bdg2,self.bdg3], func="max")

        chrom = b"chrY"
        (p,v) = bdgb.get_data_by_chr(chrom)
        pre = 0
        for i in range(len(p)):
            pos = p[i]
            value = v[i]
            self.assertEqual_float( self.test_overlie_max[i][1], pre )
            self.assertEqual_float( self.test_overlie_max[i][2], pos )
            self.assertEqual_float( self.test_overlie_max[i][3], value )
            pre = pos

    def test_overlie_mean(self):
        bdgb = self.bdg1.overlie([self.bdg2,self.bdg3], func="mean")

        chrom = b"chrY"
        (p,v) = bdgb.get_data_by_chr(chrom)
        pre = 0
        for i in range(len(p)):
            pos = p[i]
            value = v[i]
            self.assertEqual_float( self.test_overlie_mean[i][1], pre )
            self.assertEqual_float( self.test_overlie_mean[i][2], pos )
            self.assertEqual_float( self.test_overlie_mean[i][3], value )
            pre = pos

    def test_overlie_fisher(self):
        bdgb = self.bdg1.overlie([self.bdg2,self.bdg3], func="fisher")

        chrom = b"chrY"
        (p,v) = bdgb.get_data_by_chr(chrom)
        pre = 0
        for i in range(len(p)):
            pos = p[i]
            value = v[i]
            self.assertEqual_float( self.test_overlie_fisher[i][1], pre )
            self.assertEqual_float( self.test_overlie_fisher[i][2], pos )
            self.assertEqual_float( self.test_overlie_fisher[i][6], value )
            pre = pos

