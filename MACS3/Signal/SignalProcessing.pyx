# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-12-03 20:00:30 Tao Liu>

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
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
from cpython cimport bool


cpdef np.ndarray[int32_t, ndim=1] maxima(np.ndarray[float32_t, ndim=1] signal,
                                            int window_size=51):
    """return the local maxima in a signal after applying a 2nd order
    Savitsky-Golay (polynomial) filter using window_size specified
    """
    cdef:
        np.ndarray[int32_t, ndim=1] m
        np.ndarray[float64_t, ndim=1] smoothed
        np.ndarray[float64_t, ndim=1] sign, diff

    window_size = window_size//2*2+1 # to make an odd number
    smoothed = savitzky_golay_order2_deriv1(signal, window_size).round(16)
    sign = np.sign( smoothed )
    diff = np.diff( sign )
    m = np.where( diff <= -1)[0].astype("int32")
    return m

cdef np.ndarray[int32_t, ndim=1] internal_minima( np.ndarray[float32_t, ndim=1] signal,
                                                  np.ndarray[int32_t, ndim=1] maxima ):
    cdef:
        np.ndarray[int32_t, ndim=1] ret
        int32_t n = maxima.shape[0]
        int32_t i, v, v2
    if n == 0 or n == 1:
        ret = np.ndarray(0, 'int32')
        return ret
    else:
        ret = np.zeros(n - 1, 'int32')
        pos1 = maxima[0]
        for i in range(n - 1):
            pos2 = maxima[i + 1]
            ret[i] = np.where(signal[pos1:pos2] == signal[pos1:pos2].min())[0][0] + pos1
            pos1 = pos2
        return ret

cdef inline float32_t sqrt(float32_t threshold):
    return mathsqrt(threshold)

cpdef enforce_peakyness(np.ndarray[float32_t, ndim=1] signal,
                        np.ndarray[int32_t, ndim=1] maxima):
    """requires peaks described by a signal and a set of points where the signal
    is at a maximum to meet a certain set of criteria

    maxima which do not meet the required criteria are discarded

    criteria:
        for each peak:
            calculate a threshold of the maximum of its adjacent two minima
                plus the sqrt of that value
            subtract the threshold from the region bounded by those minima
            clip that region if negative values occur inside it
            require it be > 50 bp in width
            require that it not be too flat (< 6 unique values)
    """
    cdef:
        np.ndarray[int32_t, ndim=1] minima = internal_minima(signal, maxima)
        np.ndarray[float32_t, ndim=1] new_signal
        int32_t n = minima.shape[0]
        float32_t threshold
        np.ndarray[int32_t, ndim=1] peaky_maxima = maxima.copy()
        int32_t j = 0
    if n == 0: return maxima
#    else:
    threshold = signal[minima[0]]
    threshold += sqrt(threshold)
    new_signal = signal[0:minima[0]] - threshold - sqrt(threshold)
#    assert maxima[0] < minima[0], '%d > %d' % ( maxima[0], minima[0] )
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
    threshold =  signal[minima[-1]]
    threshold += sqrt(threshold)
    new_signal = signal[minima[-1]:] - threshold
    new_maximum = maxima[-1] - minima[-1]
    if is_valid_peak(new_signal, new_maximum):
        peaky_maxima[j] = maxima[-1]
        j += 1
    peaky_maxima.resize(j, refcheck=False)
    return peaky_maxima

# hardcoded minimum peak width = 50
cdef bool is_valid_peak(np.ndarray[float32_t, ndim=1] signal, int maximum):
    cdef:
        np.ndarray s
        int32_t length
    s = hard_clip(signal, maximum)
    length = s.shape[0]
    if length < 50:
        return False
    elif too_flat(s):
        return False
    return True

# require at least 6 different float values -- prevents broad flat peaks
cdef bool too_flat(np.ndarray[float32_t, ndim=1] signal):
    """return whether signal has at least 6 unique values
    """
    return np.unique(signal).shape[0] < 6

# hard clip a region with negative values
cdef np.ndarray[float32_t, ndim=1] hard_clip(np.ndarray[float32_t, ndim=1] signal, int32_t maximum):
    """clip the signal in both directions at the nearest values <= 0
    to position maximum
    """
    cdef:
        int32_t i
        int32_t left = 0
        int32_t right = signal.shape[0]
    # clip left
    for i in range( right - maximum, 0 ):
        if signal[ -i ] < 0:
            left = i
            break
    for i in range(maximum, right):
        if signal[i] < 0:
            right = i
            break
    return signal[ left:right ]

cpdef np.ndarray[ int32_t, ndim=1 ] enforce_valleys(np.ndarray[ float32_t, ndim=1 ] signal,
                                                    np.ndarray[ int32_t, ndim=1 ] summits,
                                                    float32_t min_valley = 0.8 ):
    """require a value of <= min_valley * lower summit
    between each pair of summits
    """
    cdef:
        float32_t req_min, v, prev_v
        int32_t summit_pos, prev_summit_pos
        int32_t n_summits
        int32_t n_valid_summits
        np.ndarray[ int32_t, ndim=1 ] valid_summits
    n_summits = summits.shape[0]
    n_valid_summits = 1    
    valid_summits = summits.copy()        
    # Remove peaks that do not have sufficient valleys
    if n_summits == 1: return summits
    for i in range( 1, n_summits ):
        prev_summit_pos = valid_summits[ n_valid_summits-1 ]
        summit_pos = summits[ i ]
        prev_v = signal[ prev_summit_pos ]
        v = signal[ summit_pos ]
        req_min = min_valley * min( prev_v, v )
        if ( signal[ prev_summit_pos:summit_pos ] < req_min ).any():
            valid_summits[ n_valid_summits ] = summit_pos
            n_valid_summits += 1
        elif v > prev_v:
            valid_summits[ n_valid_summits-1 ] = summit_pos
    valid_summits.resize( n_valid_summits, refcheck=False )
    return valid_summits

# Modified from http://www.scipy.org/Cookbook/SavitzkyGolay
# positive window_size not enforced anymore
# needs sane input paramters, window size > 4
# switched to double precision for internal accuracy
cpdef np.ndarray[float64_t, ndim=1] savitzky_golay_order2_deriv1(np.ndarray[float32_t, ndim=1] signal,
                                                                 int32_t window_size):
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
    cdef:
        int32_t half_window, k
        np.ndarray[int64_t, ndim=2] b
        # pad the signal at the extremes with
        # values taken from the signal itself
        np.ndarray[float32_t, ndim=1] firstvals, lastvals
        np.ndarray[float64_t, ndim=1] m, ret

    if window_size % 2 != 1: window_size += 1
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
    ret = np.convolve( m[::-1], signal.astype("float64"), mode='valid') #.astype("float32").round(8) # round to 8 decimals to avoid signing issue
    return ret

# Another modified version from http://www.scipy.org/Cookbook/SavitzkyGolay
cpdef np.ndarray[float32_t, ndim=1] savitzky_golay( np.ndarray[float32_t, ndim=1] y, int32_t window_size,
                                                    int32_t order, int32_t deriv = 0, int32_t rate = 1 ):
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
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
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
    cdef:
        int32_t half_window, k
        np.ndarray[int64_t, ndim=2] b
        # pad the signal at the extremes with
        # values taken from the signal itself
        np.ndarray[float32_t, ndim=1] firstvals, lastvals, ret
        np.ndarray[float64_t, ndim=1] m

    try:
        window_size = np.abs( np.int( window_size ) )
        order = np.abs( np.int( order ) )
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    half_window = ( window_size -1 ) // 2
    # precompute coefficients
    b = np.array( [ [ k**i for i in range( order + 1 ) ] for k in range( -half_window, half_window+1 ) ] )
    m = np.linalg.pinv( b )[ deriv ] * rate**deriv * mathfactorial( deriv )
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[ 0 ] - np.abs( y[ 1:half_window + 1 ][ ::-1 ] - y[ 0 ] )
    lastvals = y[ -1 ] + np.abs( y[ -half_window - 1:-1 ][ ::-1 ] - y[ -1 ])
    y = np.concatenate( ( firstvals, y, lastvals ) )
    ret = np.convolve( m[ ::-1 ], y, mode = 'valid' ).astype("float32")
    return ret
