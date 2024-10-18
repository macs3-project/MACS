#!/usr/bin/env python
# Time-stamp: <2024-10-15 16:07:27 Tao Liu>

import unittest
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
import numpy as np


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
        self.pileup_p = np.array([10, 50, 70, 80, 85, 100, 110, 160, 170, 180, 190], dtype="i4")
        self.pileup_v = np.array([3.0, 4.0, 6.0, 9.0, 11.0, 15.0, 19.0, 18.0, 16.0, 10.0, 6.0], dtype="f4")
        self.peak_str = "chrom:chrY	start:80	end:180	name:peak_1	score:19	summit:105\n"
        self.subset_barcodes = {b'0w#AAACGAACAAGTAACA', b"0w#AAACGAACAAGTAAGA"}
        self.subset_pileup_p = np.array([10, 50, 70, 80, 85, 100, 110, 160, 170, 180, 190], dtype="i4")
        self.subset_pileup_v = np.array([1.0, 2.0, 4.0, 6.0, 7.0, 8.0, 13.0, 12.0, 10.0, 5.0, 4.0], dtype="f4")
        self.subset_peak_str = "chrom:chrY	start:100	end:170	name:peak_1	score:13	summit:105\n"

        self.t = sum([(x[2]-x[1]) * x[4] for x in self.input_regions])

    def test_add_frag(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_loc(c, l, r, b, C)
        pe.finalize()
        # roughly check the numbers...
        self.assertEqual(pe.total, 22)
        self.assertEqual(pe.length, self.t)

        # subset
        pe_subset = pe.subset(self.subset_barcodes)
        # roughly check the numbers...
        self.assertEqual(pe_subset.total, 14)
        self.assertEqual(pe_subset.length, 1305)

    def test_pileup(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_loc(c, l, r, b, C)
        pe.finalize()
        bdg = pe.pileup_bdg()
        d = bdg.get_data_by_chr(b'chrY')  # (p, v) of ndarray
        np.testing.assert_array_equal(d[0], self.pileup_p)
        np.testing.assert_array_equal(d[1], self.pileup_v)

        pe_subset = pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg()
        d = bdg.get_data_by_chr(b'chrY')  # (p, v) of ndarray
        np.testing.assert_array_equal(d[0], self.subset_pileup_p)
        np.testing.assert_array_equal(d[1], self.subset_pileup_v)

    def test_pileup2(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_loc(c, l, r, b, C)
        pe.finalize()
        bdg = pe.pileup_bdg2()
        d = bdg.get_data_by_chr(b'chrY')  # (p, v) of ndarray
        np.testing.assert_array_equal(d['p'], self.pileup_p)
        np.testing.assert_array_equal(d['v'], self.pileup_v)

        pe_subset = pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg2()
        d = bdg.get_data_by_chr(b'chrY')  # (p, v) of ndarray
        np.testing.assert_array_equal(d['p'], self.subset_pileup_p)
        np.testing.assert_array_equal(d['v'], self.subset_pileup_v)

    def test_callpeak(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_loc(c, l, r, b, C)
        pe.finalize()
        bdg = pe.pileup_bdg()  # bedGraphTrackI object
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.peak_str)

        pe_subset = pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg()
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.subset_peak_str)

    def test_callpeak2(self):
        pe = PETrackII()
        for (c, l, r, b, C) in self.input_regions:
            pe.add_loc(c, l, r, b, C)
        pe.finalize()
        bdg = pe.pileup_bdg2()  # bedGraphTrackII object
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.peak_str)

        pe_subset = pe.subset(self.subset_barcodes)
        bdg = pe_subset.pileup_bdg2()
        peaks = bdg.call_peaks(cutoff=10, min_length=20, max_gap=10)
        self.assertEqual(str(peaks), self.subset_peak_str)
