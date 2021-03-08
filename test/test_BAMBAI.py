#!/usr/bin/env python
# Time-stamp: <2021-03-05 17:17:12 Tao Liu>

import unittest

from MACS3.IO.BAM import *

class Test_bai ( unittest.TestCase ):

    def setUp ( self ):
        self.bamfile = "test/tiny.bam"
        self.baifile = "test/tiny.bam.bai"
        self.bai = BAIFile( self.baifile )
        
    def test_load ( self ):
        meta = self.bai.get_metadata_by_refseq( 0 )
        print( meta )

    def test_get_chunks_bin ( self ):
        chunks = self.bai.get_chunks_by_bin( 0, 4728 )
        print( chunks )

    def test_get_chunks_list_bins ( self ):
        chunks = self.bai.get_chunks_by_list_bins( 0, [ 591, 4694, 4728] )
        print( chunks )

    def test_get_chunks_region ( self ):
        chunks = self.bai.get_chunks_by_region( 0, 557000, 851600 )
        print( chunks )

    def test_get_chunks_list_regions ( self ):
        chunks = self.bai.get_chunks_by_list_regions( 0, [ (500000, 600000), (800000, 850000) ] )
        print( chunks )


