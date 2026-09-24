#!/usr/bin/env python
# Time-stamp: <2025-09-29 15:22:11 Tao Liu>

"""Module Description: Test functions for Signal.pyx

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import unittest
import pytest

import numpy as np
from MACS3.Signal.SignalProcessing import enforce_peakyness, maxima

# ------------------------------------
# Main function
# ------------------------------------


class Test_maxima(unittest.TestCase):
    def setUp(self):
        self.signal = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4,
                                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                                6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8,
                                8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                                6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4,
                                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3], dtype="float32")
        self.windowsize = 253
        self.summit = 161       # this is based on 1-deriv smoothed data
        self.smoothed162 = -2.98155597e-18  # the value is based on python3.7+numpy1.17.4; from smoothing func with np.convolve

    # this test uses only numpy functions. we aim to capture strange behavior under specific py+numpy
    # @pytest.mark.skip(reason="it fails under some combinations of py+np with unknown reason.")
    def test_implement_smooth_here(self):
        signal = self.signal
        window_size = self.windowsize
        half_window = (window_size - 1) // 2
        # precompute coefficients
        b = np.array([[1, k, k**2] for k in range(-half_window, half_window+1)],
                     dtype='int64')
        m = np.linalg.pinv(b)[1]
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
        lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
        signal = np.concatenate((firstvals, signal, lastvals))
        ret = np.convolve(m[::-1], signal.astype("float64"), mode='valid').astype("float32")  # convolve function seems to have odd ret with spec py+numpy
        p = ret[162]
        print("calculated step by step:\n", p)
        print("expected:\n", self.smoothed162)
        self.assertAlmostEqual(p, self.smoothed162, places=16)
        self.assertEqual(np.sign(p), np.sign(self.smoothed162))

    def test_maxima(self):
        expect = self.summit
        result = maxima(self.signal, self.windowsize)[0]
        self.assertEqual(result, expect, msg=f"Not equal: result: {result}, expected: {expect}")

    def assertEqual_nparray1d(self, a, b, places=7):
        self.assertEqual(a.shape[0], b.shape[0])
        l = a.shape[0]
        for i in range(l):
            self.assertAlmostEqual(a[i], b[i], places=places, msg=f"Not equal at {i} {a[i]} {b[i]}")


def test_enforce_peakyness_clips_both_sides():
    """Regression test for the inactive left scan in issue #748."""
    signal = np.zeros(400, dtype="f4")
    signal[:100] = np.linspace(30, 1, 100)
    signal[100:200] = 1
    signal[200:220] = 1 + np.arange(20) * 1.5
    signal[220:240] = 1 + np.arange(20)[::-1] * 1.5
    signal[240:300] = 1
    signal[300:] = np.linspace(1, 30, 100)
    candidates = np.array([0, 219, 399], dtype="i4")

    np.testing.assert_array_equal(
        enforce_peakyness(signal, candidates),
        np.array([0, 399], dtype="i4"),
    )


def _make_width_boundary_signal(width):
    """Return three candidates with center support of exactly ``width``."""
    signal = np.ones(300, dtype="f4")
    signal[:40] = np.linspace(20, 3, 40)
    signal[260:] = np.linspace(3, 20, 40)
    start = 150 - width // 2
    support = 2 + np.minimum(np.arange(width), np.arange(width)[::-1])
    signal[start:start + width] = support
    center = start + int(np.argmax(support))
    return signal, np.array([0, center, 299], dtype="i4"), center


@pytest.mark.parametrize(("width", "expected"), [(49, False), (50, True)])
def test_enforce_peakyness_minimum_width(width, expected):
    """Fifty nonnegative bases, including zero boundaries, are required."""
    signal, candidates, center = _make_width_boundary_signal(width)

    retained = enforce_peakyness(signal, candidates)

    assert (center in retained) is expected


def test_enforce_peakyness_is_reverse_symmetric():
    """Candidate filtering should not depend on signal orientation."""
    signal = np.zeros(300, dtype="f4")
    signal[:120] = 6 + np.arange(120) * 0.02
    signal[120:180] = 5
    signal[180:] = 6 + np.arange(120) * 0.02
    candidates = np.array([119, 299], dtype="i4")
    reverse_candidates = (
        len(signal) - 1 - candidates[::-1]
    ).astype("i4")

    forward = enforce_peakyness(signal, candidates)
    reverse = enforce_peakyness(signal[::-1].copy(), reverse_candidates)
    reverse_mapped = len(signal) - 1 - reverse[::-1]

    np.testing.assert_array_equal(forward, candidates)
    np.testing.assert_array_equal(forward, reverse_mapped)
