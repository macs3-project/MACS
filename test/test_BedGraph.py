#!/usr/bin/env python
# Time-stamp: <2019-08-09 12:47:51 taoliu>

import unittest

from MACS2.IO.BedGraph import *

class Test_bedGraphTrackI_add_loc(unittest.TestCase):

    def setUp(self):
        self.test_regions1 = [("chrY",0,10593155,0.0),
                              ("chrY",10593155,10597655,0.0066254580149)]

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
        self.test_cslregions1 = [("chrY",0,70,0.00),
                              ("chrY",70,80,0.07),
                              ("chrY",80,150,0.00),
                              ("chrY",150,160,0.07),
                              ("chrY",160,190,0.00)]
        self.test_cdregions2 = [("chrY",0,85,0.00),
                             ("chrY",85,90,0.75),
                             ("chrY",90,155,0.00),
                             ("chrY",155,165,0.75),
                             ("chrY",165,200,0.00)]
        self.test_overlie_result = [("chrY",0,70,0.0),
                                    ("chrY",70,80,0.07),
                                    ("chrY",80,85,0.0),
                                    ("chrY",85,90,0.75),
                                    ("chrY",90,150,0.0),
                                    ("chrY",150,155,0.07),
                                    ("chrY",155,165,0.75),
                                    ("chrY",165,190,0.0)]

    def assertEqual_float ( self, a, b, roundn = 5 ):
        self.assertEqual( round( a, roundn ), round( b, roundn ) )

    def test_overlie(self):
        bdg1 = bedGraphTrackI()
        bdg2 = bedGraphTrackI()
        for a in self.test_cslregions1:
            bdg1.safe_add_loc(a[0],a[1],a[2],a[3])

        for a in self.test_cdregions2:
            bdg2.safe_add_loc(a[0],a[1],a[2],a[3])

        bdgb = bdg1.overlie(bdg2)

        chrom = "chrY"
        (p,v) = bdgb.get_data_by_chr(chrom)
        pre = 0
        for i in xrange(len(p)):
            pos = p[i]
            value = v[i]
            self.assertEqual_float( self.test_overlie_result[i][1], pre )
            self.assertEqual_float( self.test_overlie_result[i][2], pos )
            self.assertEqual_float( self.test_overlie_result[i][3], value )            
            pre = pos        


if __name__ == '__main__':
    unittest.main()
