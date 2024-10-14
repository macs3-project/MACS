#!/usr/bin/env python
# Time-stamp: <2024-10-11 16:20:17 Tao Liu>

import unittest

from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII


class Test_PETrackI(unittest.TestCase):
    def setUp(self):

        self.input_regions = [(b"chrY", 0, 100),
                              (b"chrY", 70, 270),
                              (b"chrY", 70, 100),
                              (b"chrY", 80, 160),
                              (b"chrY", 80, 160),
                              (b"chrY", 80, 180),
                              (b"chrY", 80, 180),
                              (b"chrY", 85, 185),
                              (b"chrY", 85, 285),
                              (b"chrY", 85, 285),
                              (b"chrY", 85, 285),
                              (b"chrY", 85, 385),
                              (b"chrY", 90, 190),
                              (b"chrY", 90, 190),
                              (b"chrY", 90, 191),
                              (b"chrY", 150, 190),
                              (b"chrY", 150, 250),
                              ]
        self.t = sum([x[2]-x[1] for x in self.input_regions])

    def test_add_loc(self):
        pe = PETrackI()
        for (c, l, r) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        # roughly check the numbers...
        self.assertEqual(pe.total, 17)
        self.assertEqual(pe.length, self.t)

    def test_filter_dup(self):
        pe = PETrackI()
        for (c, l, r) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        # roughly check the numbers...
        self.assertEqual(pe.total, 17)
        self.assertEqual(pe.length, self.t)

        # filter out more than 3 tags
        pe.filter_dup(3)
        self.assertEqual(pe.total, 17)

        # filter out more than 2 tags
        pe.filter_dup(2)
        self.assertEqual(pe.total, 16)

        # filter out more than 1 tag
        pe.filter_dup(1)
        self.assertEqual(pe.total, 12)

    def test_sample_num(self):
        pe = PETrackI()
        for (c, l, r) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        pe.sample_num(10)
        self.assertEqual(pe.total, 10)

    def test_sample_percent(self):
        pe = PETrackI()
        for (c, l, r) in self.input_regions:
            pe.add_loc(c, l, r)
        pe.finalize()
        pe.sample_percent(0.5)
        self.assertEqual(pe.total, 8)


class Test_PETrackII(unittest.TestCase):
    def setUp(self):
        self.input_regions = [(b"chrY", 0, 100, b"0w#AAACGAAAGACTCGGA", 2),
                              (b"chrY", 70, 170, b"0w#AAACGAAAGACTCGGA", 1),
                              (b"chrY", 80, 190, b"0w#AAACGAAAGACTCGGA", 1),
                              (b"chrY", 85, 180, b"0w#AAACGAAAGACTCGGA", 3),
                              (b"chrY", 100, 190, b"0w#AAACGAAAGACTCGGA", 1),
                              (b"chrY", 0, 100, b"0w#AAACGAACAAGTAACA", 1),
                              (b"chrY", 70, 170, b"0w#AAACGAACAAGTAACA", 2),
                              (b"chrY", 80, 190, b"0w#AAACGAACAAGTAACA", 1),
                              (b"chrY", 85, 180, b"0w#AAACGAACAAGTAACA", 1),
                              (b"chrY", 100, 190, b"0w#AAACGAACAAGTAACA", 3),
                              (b"chrY", 10, 110, b"0w#AAACGAACAAGTAAGA", 1),
                              (b"chrY", 50, 160, b"0w#AAACGAACAAGTAAGA", 2),
                              (b"chrY", 100, 170, b"0w#AAACGAACAAGTAAGA", 3)
                              ]
        self.t = sum([(x[2]-x[1]) * x[4] for x in self.input_regions])

    def test_add_frag(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_frag(c, l, r, b, C)
        pe.finalize()
        # roughly check the numbers...
        self.assertEqual(pe.total, 22)
        self.assertEqual(pe.length, self.t)

    def test_subset(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_frag(c, l, r, b, C)
        pe.finalize()
        pe_subset = pe.subset(set([b'0w#AAACGAACAAGTAACA', b"0w#AAACGAACAAGTAAGA"]))
        # roughly check the numbers...
        self.assertEqual(pe_subset.total, 14)
        self.assertEqual(pe_subset.length, 1305)
