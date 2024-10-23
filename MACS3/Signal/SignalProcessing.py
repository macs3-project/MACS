# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-15 11:25:35 Tao Liu>

"""Module Description: functions to find maxima minima or smooth the
signal tracks.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# Python modules
# ------------------------------------
from math import factorial as mathfactorial
from math import sqrt as mathsqrt

# ------------------------------------
# Other modules
# ------------------------------------
# smoothing function
import numpy as np
import cython
import cython.cimports.numpy as cnp
from cython.cimports.cpython import bool


@cython.ccall
def maxima(signal: cnp.ndarray(cython.float, ndim=1),
           window_size: cython.int = 51) -> cnp.ndarray:
    """return the local maxima in a signal after applying a 2nd order
    Savitsky-Golay (polynomial) filter using window_size specified
    """
    m: cnp.ndarray(cython.int, ndim=1)
    smoothed: cnp.ndarray(cython.double, ndim=1)
    sign: cnp.ndarray(cython.double, ndim=1)
    diff: cnp.ndarray(cython.double, ndim=1)

    window_size = window_size//2*2+1  # to make an odd number
    smoothed = savitzky_golay_order2_deriv1(signal, window_size).round(16)
    sign = np.sign(smoothed)
    diff = np.diff(sign)
    m = np.where(diff <= -1)[0].astype("i4")
    return m


@cython.cfunc
def internal_minima(signal: cnp.ndarray(cython.float, ndim=1),
                    maxima: cnp.ndarray(cython.int, ndim=1)) -> cnp.ndarray:
    ret: cnp.ndarray(cython.int, ndim=1)
    n: cython.int = maxima.shape[0]
    i: cython.int

    if n == 0 or n == 1:
        ret = np.ndarray(0, 'i4')
        return ret
    else:
        ret = np.zeros(n - 1, 'i4')
        pos1 = maxima[0]
        for i in range(n - 1):
            pos2 = maxima[i + 1]
            ret[i] = np.where(signal[pos1:pos2] == signal[pos1:pos2].min())[0][0] + pos1
            pos1 = pos2
        return ret


@cython.cfunc
@cython.inline
def sqrt(threshold: cython.float) -> cython.float:
    return mathsqrt(threshold)


@cython.ccall
def enforce_peakyness(signal: cnp.ndarray(cython.float, ndim=1),
                      maxima: cnp.ndarray(cython.int, ndim=1)):
    """requires peaks described by a signal and a set of points where
    the signal is at a maximum to meet a certain set of criteria

    maxima which do not meet the required criteria are discarded

    criteria:
        for each peak:

            calculate a threshold of the maximum of its adjacent two
            minima plus the sqrt of that value

            subtract the threshold from the region bounded by those minima

            clip that region if negative values occur inside it

            require it be > 50 bp in width -- controlled by is_valied_peak()

            require that it not be too flat (< 6 unique values) --
            controlled by is_valid_peak()

    """
    minima: cnp.ndarray(cython.int, ndim=1) = internal_minima(signal, maxima)
    new_signal: cnp.ndarray(cython.float, ndim=1)
    n: cython.int = minima.shape[0]
    threshold: cython.float
    peaky_maxima: cnp.ndarray(cython.int, ndim=1) = maxima.copy()
    j: cython.int = 0

    if n == 0:
        return maxima

    threshold = signal[minima[0]]
    threshold += sqrt(threshold)
    new_signal = signal[0:minima[0]] - threshold - sqrt(threshold)

    if is_valid_peak(new_signal, maxima[0]):
        peaky_maxima[0] = maxima[0]
        j += 1
    for i in range(n - 1):
        threshold = max(signal[minima[i]], signal[minima[i + 1]])
        threshold += sqrt(threshold)
        new_signal = signal[minima[i]:minima[i+1]] - threshold
        new_maximum = maxima[i+1] - minima[i]
        if is_valid_peak(new_signal, new_maximum):
            peaky_maxima[j] = maxima[i + 1]
            j += 1
    threshold = signal[minima[-1]]
    threshold += sqrt(threshold)
    new_signal = signal[minima[-1]:] - threshold
    new_maximum = maxima[-1] - minima[-1]
    if is_valid_peak(new_signal, new_maximum):
        peaky_maxima[j] = maxima[-1]
        j += 1
    peaky_maxima.resize(j, refcheck=False)
    return peaky_maxima


# hardcoded minimum peak width = 50
@cython.cfunc
def is_valid_peak(signal: cnp.ndarray(cython.float, ndim=1),
                  maximum: cython.int) -> bool:
    s: cnp.ndarray
    length: cython.int

    s = hard_clip(signal, maximum)
    length = s.shape[0]
    if length < 50:
        return False
    elif too_flat(s):
        return False
    return True


# require at least 6 different float values -- prevents broad flat peaks
@cython.cfunc
def too_flat(signal: cnp.ndarray(cython.float, ndim=1)) -> bool:
    """return whether signal has at least 6 unique values
    """
    return np.unique(signal).shape[0] < 6


# hard clip a region with negative values
@cython.cfunc
def hard_clip(signal: cnp.ndarray(cython.float, ndim=1),
              maximum: cython.int) -> cnp.ndarray:
    """clip the signal in both directions at the nearest values <= 0
    to position maximum
    """
    i: cython.int
    left: cython.int = 0
    right: cython.int = signal.shape[0]

    # clip left
    for i in range(right - maximum, 0):
        if signal[-i] < 0:
            left = i
            break
    for i in range(maximum, right):
        if signal[i] < 0:
            right = i
            break
    return signal[left:right]


@cython.ccall
def enforce_valleys(signal: cnp.ndarray(cython.float, ndim=1),
                    summits: cnp.ndarray(cython.int, ndim=1),
                    min_valley: cython.float = 0.8) -> cnp.ndarray:
    """require a value of <= min_valley * lower summit
    between each pair of summits
    """
    req_min: cython.float
    v: cython.float
    prev_v: cython.float

    summit_pos: cython.int
    prev_summit_pos: cython.int
    n_summits: cython.int
    n_valid_summits: cython.int

    valid_summits: cnp.ndarray(cython.int, ndim=1)

    n_summits = summits.shape[0]
    n_valid_summits = 1
    valid_summits = summits.copy()
    # Remove peaks that do not have sufficient valleys
    if n_summits == 1:
        return summits
    for i in range(1, n_summits):
        prev_summit_pos = valid_summits[n_valid_summits-1]
        summit_pos = summits[i]
        prev_v = signal[prev_summit_pos]
        v = signal[summit_pos]
        req_min = min_valley * min(prev_v, v)
        if (signal[prev_summit_pos:summit_pos] < req_min).any():
            valid_summits[n_valid_summits] = summit_pos
            n_valid_summits += 1
        elif v > prev_v:
            valid_summits[n_valid_summits-1] = summit_pos
    valid_summits.resize(n_valid_summits, refcheck=False)
    return valid_summits


# Modified from http://www.scipy.org/Cookbook/SavitzkyGolay
# positive window_size not enforced anymore
# needs sane input paramters, window size > 4
# switched to double precision for internal accuracy
@cython.ccall
def savitzky_golay_order2_deriv1(signal: cnp.ndarray(cython.float, ndim=1),
                                 window_size: cython.int) -> cnp.ndarray:
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    half_window: cython.int
    b: cnp.ndarray(cython.long, ndim=2)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals: cnp.ndarray(cython.float, ndim=1)
    lastvals: cnp.ndarray(cython.float, ndim=1)
    m: cnp.ndarray(cython.double, ndim=1)
    ret: cnp.ndarray(cython.double, ndim=1)

    if window_size % 2 != 1:
        window_size += 1
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.array([[1, k, k**2] for k in range(-half_window, half_window+1)],
                 dtype='i8')
    m = np.linalg.pinv(b)[1]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
    lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
    signal = np.concatenate((firstvals, signal, lastvals))
    ret = np.convolve(m[::-1],
                      signal.astype("f8"),
                      mode='valid')
    return ret


# Another modified version from http://www.scipy.org/Cookbook/SavitzkyGolay
@cython.ccall
def savitzky_golay(y: cnp.ndarray(cython.float, ndim=1),
                   window_size: cython.int,
                   order: cython.int,
                   deriv: cython.int = 0,
                   rate: cython.int = 1) -> cnp.ndarray:
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp(-t**2) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    half_window: cython.int
    b: cnp.ndarray(cython.long, ndim=2)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals: cnp.ndarray(cython.float, ndim=1)
    lastvals: cnp.ndarray(cython.float, ndim=1)
    ret: cnp.ndarray(cython.float, ndim=1)
    m: cnp.ndarray(cython.double, ndim=1)

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.array([[k**i
                   for i in range(order + 1)]
                  for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b)[deriv] * rate**deriv * mathfactorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    ret = np.convolve(m[::-1], y, mode='valid').astype("float32")
    return ret
