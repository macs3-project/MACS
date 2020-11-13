#!/usr/bin/env python
# Time-stamp: <2019-12-12 14:42:28 taoliu>

import unittest

from MACS3.IO.Parser import *

class Test_auto_guess ( unittest.TestCase ):

    def setUp ( self ):
        self.bedfile = "test/tiny.bed.gz"
        self.bedpefile = "test/tiny.bedpe.gz"
        self.samfile = "test/tiny.sam.gz"
        self.bamfile = "test/tiny.bam"

    def test_guess_parser_bed ( self ):
        p = guess_parser( self.bedfile )
        self.assertTrue( p.is_gzipped() )
        self.assertTrue( isinstance(p, BEDParser) )

    def test_guess_parser_sam ( self ):
        p = guess_parser( self.samfile )
        self.assertTrue( p.is_gzipped() )
        self.assertTrue( isinstance(p, SAMParser) )

    def test_guess_parser_bam ( self ):
        p = guess_parser( self.bamfile )
        self.assertTrue( p.is_gzipped() )
        self.assertTrue( isinstance(p, BAMParser) )



