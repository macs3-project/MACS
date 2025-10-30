#!/usr/bin/env python
# Time-stamp: <2025-09-29 14:18:42 Tao Liu>

import unittest

from MACS3.IO.BAM import BAIFile, BAMaccessor


class Test_BAIFile(unittest.TestCase):

    def setUp(self):
        self.bamfile = "test/tiny.bam"
        self.baifile = "test/tiny.bam.bai"
        self.bai = BAIFile(self.baifile)

    def test_load(self):
        expected = {'ref_beg': 50724864, 'ref_end': 2709258240,
                    'n_mapped': 971, 'n_unmapped': 0}
        meta = self.bai.get_metadata_by_refseq(0)
        self.assertDictEqual(expected, meta)

    def test_get_chunks_bin(self):
        expected = [(2178362419, 2178363422)]
        chunks = self.bai.get_chunks_by_bin(0, 4728)
        self.assertListEqual(chunks, expected)

    def test_get_chunks_list_bins(self):
        expected = [(50741975, 50751615),
                    (2178362419, 2178363422),
                    (2178363422, 2709258240)]
        chunks = self.bai.get_chunks_by_list_of_bins(0, [591, 4694, 4728])
        self.assertListEqual(chunks, expected)

    def test_get_chunks_region(self):
        expected = [(1130775424, 1130820164), (1130820164, 1130822629),
                    (1130822629, 2178352280), (2178352280, 2178353973),
                    (2178353973, 2178356807), (2178356807, 2178358139),
                    (2178358139, 2178359954), (2178359954, 2178362419),
                    (2178362419, 2178363422), (2178363422, 2709258240)]
        chunks = self.bai.get_chunks_by_region(0, 557000, 851600)
        self.assertListEqual(chunks, expected)

    def test_get_chunks_list_regions(self):
        expected = [(1130770109, 1130772126), (1130772126, 1130775424),
                    (1130775424, 1130820164), (2178363422, 2709258240)]
        chunks = self.bai.get_chunks_by_list_of_regions(0, [(500000, 600000),
                                                            (800000, 850000)])
        self.assertListEqual(chunks, expected)


class Test_BAM_w_BAI(unittest.TestCase):

    def setUp(self):
        self.bamfile = "test/tiny.bam"
        self.baifile = "test/tiny.bam.bai"
        self.bam = BAMaccessor(self.bamfile)

    def test_get_reads_large_dedup(self):
        c = b"chr10"
        s = 100
        e = 900000
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=1)
        self.assertEqual(len(a), 941)

    def test_get_reads_large(self):
        c = b"chr10"
        s = 100
        e = 900000
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=100)
        self.assertEqual(len(a), 971)

    def test_get_reads_mod(self):
        c = b"chr10"
        s = 100000
        e = 550000
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=100)
        self.assertEqual(len(a), 520)

    def test_get_reads_few(self):
        c = b"chr10"
        s = 142999
        e = 163000
        expected = """chr10	145930	145966	MAGNUM:8:93:1052:1153#0	36	+
chr10	148133	148169	MAGNUM:8:108:19082:3782#0	36	-
chr10	148930	148966	MAGNUM:8:34:18198:13544#0	36	+
chr10	152575	152611	ROCKFORD:1:114:13238:15292#0	36	+
chr10	152927	152963	ROCKFORD:1:21:4055:6592#0	36	-
chr10	153995	154031	MAGNUM:8:44:12716:12267#0	36	+
chr10	156153	156189	MAGNUM:8:74:16381:18650#0	36	-
chr10	159259	159295	ROCKFORD:1:112:5095:4754#0	36	+"""
        a = self.bam.get_reads_in_region(c, s, e, maxDuplicate=100)
        result = "\n".join(map(str, a))
        self.assertEqual(len(a), 8)
        self.assertEqual(result, expected)

    def test_get_empty1(self):
        # region outside of bams
        c = b"chr10"
        s = 961000
        e = 963000
        a = self.bam.get_reads_in_region(c, s, e)
        self.assertEqual(len(a), 0)

    def test_get_empty2(self):
        # small region missed in bam
        c = b"chr10"
        s = 70000
        e = 90000
        a = self.bam.get_reads_in_region(c, s, e)
        self.assertEqual(len(a), 0)
