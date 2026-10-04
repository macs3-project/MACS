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

from scipy.signal import savgol_filter

from MACS3.Signal.SignalProcessing import (enforce_peakyness,
                                           enforce_valleys,
                                           savitzky_golay_order2_deriv1)

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


# ------------------------------------
# Comprehensive tests (added below the original ones)
# ------------------------------------
# C-only helpers are covered through their callers:
#   internal_minima, sqrt, is_valid_peak, too_flat, hard_clip
#                                              -> enforce_peakyness

def f4(values):
    return np.asarray(values, dtype="f4")


def i4(values):
    return np.asarray(values, dtype="i4")


def macs_pad(y, half_window, dtype="f4"):
    """The padding both Savitzky-Golay functions apply, in the input's
    precision (float32 for MACS3's callers): the left end becomes
    y0 - |y_k - y0| and the right end y_-1 + |y_-1-k - y_-1| for
    k = half_window..1."""
    y = np.asarray(y, dtype=dtype)
    first = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    last = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    return np.concatenate((first, y, last)).astype("f8")


def sg_deriv1_reference(y, window, dtype="f4"):
    """Slope of the least-squares quadratic over a centred window. For a
    symmetric window the linear coefficient of the quadratic fit equals
    that of the straight-line fit, c_k = k / sum(k**2)."""
    if window % 2 != 1:
        window += 1
    h = (window - 1) // 2
    k = np.arange(-h, h + 1, dtype="f8")
    c = k / np.sum(k * k)
    padded = macs_pad(y, h, dtype)
    return np.array([np.dot(c, padded[i:i + window]) for i in range(len(y))])


def maxima_reference(signal, window, dtype="f4"):
    window = window // 2 * 2 + 1
    d = np.round(sg_deriv1_reference(signal, window, dtype), 16)
    return np.where(np.diff(np.sign(d)) <= -1)[0]


def two_gaussians(n=200, c1=60.3, c2=140.6, sd=8.0, h2=0.7):
    i = np.arange(n, dtype="f8")
    y = (np.exp(-(i - c1)**2 / (2 * sd**2))
         + h2 * np.exp(-(i - c2)**2 / (2 * sd**2)))
    return f4(y)


def random_signal(n=200, seed=0):
    return f4(np.random.default_rng(seed).random(n) * 100)


# ------------------------------------
# savitzky_golay_order2_deriv1
# ------------------------------------

def test_sg_deriv1_quadratic_by_hand():
    # y = k**2, window 3: c = (-1/2, 0, 1/2), interior slope (y[n+1] -
    # y[n-1]) / 2 = 2n. Left pad y0 - |y1 - y0| = -1 gives (1 + 1)/2 = 1
    # at n = 0; right pad 36 + |25 - 36| = 47 gives (47 - 25)/2 = 11.
    y = f4([0, 1, 4, 9, 16, 25, 36])
    result = savitzky_golay_order2_deriv1(y, 3)
    assert result == pytest.approx([1, 2, 4, 6, 8, 10, 11], abs=1e-12)


def test_sg_deriv1_increasing_line_is_constant_slope():
    # rising ends are padded by point reflection, so the line continues
    y = f4(3 * np.arange(20) + 2)
    assert savitzky_golay_order2_deriv1(y, 5) == pytest.approx(
        np.full(20, 3.0), abs=1e-12)


def test_sg_deriv1_decreasing_line_mirrored_ends():
    # y = 50 - 2k: a falling start pads as y0 - |y_k - y0| = y_k (mirror),
    # so the slope at n = 0 is 0; at n = 1 the window is (48, 50, 48, 46,
    # 44) and sum(k*y)/10 = -1.2. The end is symmetric.
    y = f4(50 - 2 * np.arange(20))
    expected = np.full(20, -2.0)
    expected[[0, -1]] = 0.0
    expected[[1, -2]] = -1.2
    assert savitzky_golay_order2_deriv1(y, 5) == pytest.approx(expected,
                                                               abs=1e-12)


@pytest.mark.parametrize("window", [3, 5, 11, 51, 10])
def test_sg_deriv1_matches_reference_convolution(window):
    y = random_signal()
    result = savitzky_golay_order2_deriv1(y, window)
    assert result.dtype == np.float64
    assert result.shape == (200,)
    assert result == pytest.approx(sg_deriv1_reference(y, window),
                                   rel=1e-9, abs=1e-9)


@pytest.mark.parametrize("window", [5, 11, 51])
def test_sg_deriv1_interior_matches_scipy_savgol(window):
    y = random_signal(seed=1)
    h = (window - 1) // 2
    expected = savgol_filter(y.astype("f8"), window, 2, deriv=1)
    result = savitzky_golay_order2_deriv1(y, window)
    assert result[h:-h] == pytest.approx(expected[h:-h], abs=1e-9)


@pytest.mark.parametrize("even", [4, 10, 50])
def test_sg_deriv1_even_window_uses_next_odd(even):
    y = random_signal(seed=2)
    np.testing.assert_array_equal(savitzky_golay_order2_deriv1(y, even),
                                  savitzky_golay_order2_deriv1(y, even + 1))


def test_sg_deriv1_float64_input_computed_in_float64():
    """The float32 buffer annotation is not enforced in this pure-Python
    mode module: float64 input is accepted and padded in float64."""
    y = random_signal(seed=8).astype("f8") / 7.0
    result = savitzky_golay_order2_deriv1(y, 11)
    assert result.dtype == np.float64
    assert result == pytest.approx(sg_deriv1_reference(y, 11, "f8"),
                                   rel=1e-9, abs=1e-12)
    interior = savgol_filter(y, 11, 2, deriv=1)
    assert result[5:-5] == pytest.approx(interior[5:-5], abs=1e-12)


def test_sg_deriv1_2d_input_raises_valueerror():
    # the 1-D annotation is not enforced; np.convolve rejects the 2-D
    # padded array
    with pytest.raises(ValueError, match="object too deep for desired array"):
        savitzky_golay_order2_deriv1(np.zeros((4, 5), dtype="f4"), 3)


def test_sg_deriv1_empty_signal_raises_indexerror():
    with pytest.raises(IndexError):
        savitzky_golay_order2_deriv1(f4([]), 5)


# ------------------------------------
# maxima
# ------------------------------------

def test_maxima_two_gaussian_peaks():
    # Peaks centred at 60.3 and 140.6: every symmetric pair around 60
    # (and 140) leans right, around 61 (and 141) leans left, so the
    # smoothed slope turns from + to - between 60/61 and 140/141.
    result = maxima(two_gaussians(), 11)
    assert result.dtype == np.int32
    np.testing.assert_array_equal(result, [60, 140])


@pytest.mark.parametrize("window", [5, 11, 21])
def test_maxima_gaussians_independent_of_window(window):
    np.testing.assert_array_equal(maxima(two_gaussians(), window), [60, 140])


@pytest.mark.parametrize("window", [11, 25, 51])
def test_maxima_matches_reference_on_random_signal(window):
    y = random_signal(500, seed=3)
    np.testing.assert_array_equal(maxima(y, window),
                                  maxima_reference(y, window))


@pytest.mark.parametrize("even", [10, 24, 50])
def test_maxima_even_window_uses_next_odd(even):
    y = random_signal(500, seed=4)
    np.testing.assert_array_equal(maxima(y, even), maxima(y, even + 1))


def test_maxima_default_window_is_51():
    y = random_signal(500, seed=5)
    np.testing.assert_array_equal(maxima(y), maxima(y, 51))


@pytest.mark.parametrize("value", [0.0, 1.0, 1000.0])
def test_maxima_constant_signal_has_none(value):
    # every window sees the same values: one slope everywhere, no sign
    # change
    result = maxima(np.full(100, value, dtype="f4"), 11)
    assert result.dtype == np.int32
    assert result.shape == (0,)


def test_maxima_window_1_has_none():
    # the slope of a one-point fit is 0 everywhere
    assert maxima(two_gaussians(), 1).shape == (0,)
    assert maxima(two_gaussians(), 0).shape == (0,)


def test_maxima_float64_input_accepted():
    # the float32 annotation is not enforced; the smoothing runs in
    # float64 and the result is still int32
    y = random_signal(500, seed=9).astype("f8") / 3.0
    result = maxima(y, 11)
    assert result.dtype == np.int32
    np.testing.assert_array_equal(result, maxima_reference(y, 11, "f8"))


def test_maxima_list_input_raises_typeerror():
    with pytest.raises(TypeError, match="unsupported operand type"):
        maxima([0.0, 1.0, 2.0, 1.0, 0.0], 3)


def test_maxima_empty_signal_raises_indexerror():
    with pytest.raises(IndexError):
        maxima(f4([]), 5)


# ------------------------------------
# enforce_peakyness
# ------------------------------------

def triangle(center, start, stop, top=51):
    """top - |i - center| for i in [start, stop)."""
    return [top - abs(i - center) for i in range(start, stop)]


def test_enforce_peakyness_fewer_than_two_maxima_returned_as_is():
    sig = f4(triangle(50, 0, 101))
    for m in (i4([50]), i4([])):
        assert enforce_peakyness(sig, m) is m


def test_enforce_peakyness_keeps_two_wide_peaks():
    # min 1 at 100. First region [0, 100) minus (2 + sqrt 2) first goes
    # negative at 98 (>= 50 long); last region [100, 201) minus 2 at 100.
    sig = f4(triangle(50, 0, 101) + triangle(150, 101, 201))
    np.testing.assert_array_equal(enforce_peakyness(sig, i4([50, 150])),
                                  [50, 150])


def test_enforce_peakyness_drops_narrow_peak_near_region_start():
    # middle region [100, 125): the 21-bp spike at 105-125 leaves a
    # region of 25 < 50 positions
    sig = (triangle(50, 0, 101) + [1] * 4
           + [1 + 3 * (10 - abs(i - 115)) for i in range(105, 126)]
           + [1] * 24 + triangle(200, 150, 251))
    result = enforce_peakyness(f4(sig), i4([50, 115, 200]))
    np.testing.assert_array_equal(result, [50, 200])


def test_enforce_peakyness_drops_flat_peak():
    # middle region [100, 230) minus 2 holds only {-1, 38}: < 6 values
    sig = (triangle(50, 0, 101) + [1] * 49 + [40] * 80 + [1] * 70
           + triangle(350, 300, 401))
    result = enforce_peakyness(f4(sig), i4([50, 190, 350]))
    np.testing.assert_array_equal(result, [50, 350])


@pytest.mark.parametrize("levels, kept", [(5, False), (6, True)])
def test_enforce_peakyness_six_distinct_values_needed(levels, kept):
    """Re-derived for upstream 9597df4 (#750, issue #748): hard_clip now
    clips on the left too, so the minimum itself (-1 after subtracting the
    threshold) is no longer part of the clipped region. Before, the region
    started at the minimum and held levels + 1 distinct values, so 5 levels
    were enough; now it holds `levels` distinct values and 6 are needed.
    """
    # last region after the minimum (value 1, -1 after subtracting 2):
    # `levels` plateaus of 15 positions at 3, 4, ... (1, 2, ... after
    # subtracting 2); the clipped region is the plateaus only, so its
    # distinct values are `levels`
    stairs = []
    for v in range(3, 3 + levels):
        stairs += [v] * 15
    sig = f4(triangle(50, 0, 101) + stairs)
    top = 101 + 15 * (levels - 1)
    expected = [50, top] if kept else [50]
    np.testing.assert_array_equal(enforce_peakyness(sig, i4([50, top])),
                                  expected)


@pytest.mark.parametrize("plateau, kept", [(1, False), (2, True)])
def test_enforce_peakyness_minimum_width_50(plateau, kept):
    """Re-derived for upstream 9597df4 (#750, issue #748): hard_clip now
    clips on the left too, so the clipped region runs from relative
    position 1 (the first value after the minimum, which is negative) to
    the first negative value after the top, a width of first_neg - 1.
    Before, it started at the minimum (relative position 0) and had width
    first_neg, so one extra top value was enough to reach 50.
    """
    # after the minimum at 100: rise 3..26 (24 values), `plateau` extra
    # top values, fall 25..2, then 1. The first value < 2 after the top is
    # at relative position 49 + plateau, so the clipped region [1,
    # first_neg) is 48 + plateau wide: 49 (dropped) or 50 (kept).
    rise = list(range(3, 27))
    top_extra = [26] * plateau
    region = rise + top_extra + list(range(25, 1, -1)) + [1] * 5
    sig = f4(triangle(50, 0, 101) + region)
    rel = np.arange(1, len(region) + 1)
    first_neg = rel[(np.array(region) < 2) & (rel > 24)][0]
    assert first_neg == 49 + plateau
    assert first_neg - 1 == (50 if kept else 49)
    top = 100 + 24
    expected = [50, top] if kept else [50]
    np.testing.assert_array_equal(enforce_peakyness(sig, i4([50, top])),
                                  expected)


def test_enforce_peakyness_negative_minimum_raises_valueerror():
    sig = f4(triangle(10, 0, 20) + [-1] + triangle(30, 21, 40))
    # the math module's message is worded differently from Python 3.14
    with pytest.raises(ValueError, match="math domain error|expected a (positive|nonnegative) input"):
        enforce_peakyness(sig, i4([10, 30]))


def test_enforce_peakyness_int64_maxima_accepted():
    # the int32 annotation is not enforced; the result keeps the input
    # dtype
    sig = f4(triangle(50, 0, 101) + triangle(150, 101, 201))
    result = enforce_peakyness(sig, np.array([50, 150], dtype="i8"))
    assert result.dtype == np.int64
    np.testing.assert_array_equal(result, [50, 150])


def test_enforce_peakyness_first_peak_uses_same_threshold():
    """Regression test: the first peak subtracted sqrt(threshold) a second
    time, unlike the other peaks.

    Fixed upstream in 9597df4 (#750, issue #748).
    """
    # min 4 at 80, threshold 4 + 2 = 6. The first region stays >= 7 > 6
    # for 80 positions, but 7 < 6 + sqrt(6) cuts it at 41.
    sig = ([i + 10 for i in range(31)] + [40 - 3 * k for k in range(1, 12)]
           + [7] * 38 + [4] + [5 + j for j in range(60)]
           + [64 - j for j in range(1, 61)])
    result = enforce_peakyness(f4(sig), i4([30, 140]))
    np.testing.assert_array_equal(result, [30, 140])


def test_enforce_peakyness_clips_left_side():
    """Regression test: hard_clip never clipped on the left
    (range(right - maximum, 0) is empty), so a narrow peak far from its
    left minimum passed.

    Fixed upstream in 9597df4 (#750, issue #748), which rewrote hard_clip
    to walk outwards from the maximum in both directions.
    """
    # middle region [100, 280): a 21-bp spike at 260-280 after a floor;
    # clipping at the nearest negative left of the summit leaves ~20 bp
    sig = (triangle(50, 0, 101) + [1] * 159
           + [1 + 3 * (10 - abs(i - 270)) for i in range(260, 281)]
           + [1] * 19 + triangle(350, 300, 400))
    result = enforce_peakyness(f4(sig), i4([50, 270, 350]))
    np.testing.assert_array_equal(result, [50, 350])


# ------------------------------------
# enforce_valleys
# ------------------------------------

@pytest.mark.parametrize("signal, summits, expected", [
    # valley 2 < 0.8 * 10
    ([0, 10, 2, 10, 0], [1, 3], [1, 3]),
    # no valley, equal heights: the second summit is dropped
    ([0, 10, 9, 10, 0], [1, 3], [1]),
    # no valley, the second summit is higher: it replaces the first
    ([0, 10, 9, 12, 0], [1, 3], [3]),
    # valley exactly 0.8 * 10 = 8 is not below the requirement
    ([0, 10, 8, 10, 0], [1, 3], [1]),
    ([0, 10, 7.99, 10, 0], [1, 3], [1, 3]),
    # the requirement uses the lower summit: 0.8 * 5 = 4
    ([0, 10, 4.5, 5, 0], [1, 3], [1]),
    ([0, 10, 3.5, 5, 0], [1, 3], [1, 3]),
    # chain: 3 replaces 1, then 5 is separated from 3 by the 2
    ([0, 10, 9, 12, 2, 11, 0], [1, 3, 5], [3, 5]),
])
def test_enforce_valleys(signal, summits, expected):
    result = enforce_valleys(f4(signal), i4(summits))
    assert result.dtype == np.int32
    np.testing.assert_array_equal(result, expected)


def test_enforce_valleys_min_valley_zero_merges_to_highest():
    # nothing is below 0, so summits merge and the higher one is kept
    sig = f4([5, 1, 3, 0.5, 7])
    np.testing.assert_array_equal(enforce_valleys(sig, i4([0, 2, 4]), 0.0),
                                  [4])


def test_enforce_valleys_min_valley_one():
    # 9.5 < 1.0 * 10
    sig = f4([0, 10, 9.5, 10, 0])
    np.testing.assert_array_equal(enforce_valleys(sig, i4([1, 3]), 1.0),
                                  [1, 3])
    np.testing.assert_array_equal(enforce_valleys(sig, i4([1, 3])), [1])


def test_enforce_valleys_default_min_valley_is_0_8():
    sig = random_signal(100, seed=6)
    s = i4([10, 30, 50, 70, 90])
    np.testing.assert_array_equal(enforce_valleys(sig, s),
                                  enforce_valleys(sig, s, 0.8))


def test_enforce_valleys_single_summit_returned_as_is():
    s = i4([3])
    assert enforce_valleys(f4([0, 1, 2, 3, 2]), s) is s


def test_enforce_valleys_does_not_modify_input():
    s = i4([1, 3])
    enforce_valleys(f4([0, 10, 9, 12, 0]), s)
    np.testing.assert_array_equal(s, [1, 3])


def test_enforce_valleys_int64_summits_and_float64_signal_accepted():
    # the dtype annotations are not enforced; the result keeps the
    # summits' dtype
    result = enforce_valleys(np.array([0, 10, 2, 10, 0], dtype="f8"),
                             np.array([1, 3], dtype="i8"))
    assert result.dtype == np.int64
    np.testing.assert_array_equal(result, [1, 3])


# ------------------------------------
# savitzky_golay
# ------------------------------------
# Every call fails before any computation: np.int was removed from
# numpy (1.24), and only ValueError is caught around it. The tests
# below state the intended behaviour: least-squares polynomial
# smoothing over a centred window, with the same end padding as
# savitzky_golay_order2_deriv1, returned as float32.


