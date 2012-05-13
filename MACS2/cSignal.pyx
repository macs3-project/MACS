# smoothing function 
import numpy as np
cimport numpy as np

cpdef maxima(np.ndarray[np.float32_t, ndim=1] signal,
             int window_size=51):
    cdef np.ndarray[np.int32_t, ndim=1] m = np.where(np.diff(np.sign(savitzky_golay_order2(signal, window_size, deriv=1))) == 2)[0].astype('int32')
    return m

cpdef enforce_valleys(np.ndarray[np.float32_t, ndim=1] signal,
                      np.ndarray[np.int32_t, ndim=1] summits,
                      float min_valley = 0.8):
    cdef:
        float req_min, v, prev_v
        int summit_pos, prev_summit_pos
        int n_summits = summits.shape[0]
        int n_valid_summits = 1
        np.ndarray[np.int32_t, ndim=1] valid_summits = summits.copy()
    if n_summits == 1: return summits
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
cpdef savitzky_golay_order2(np.ndarray[np.float32_t, ndim=1] signal,
                     int window_size, int deriv=0):
#    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
#    The Savitzky-Golay filter removes high frequency noise from data.
#    It has the advantage of preserving the original shape and
#    features of the signal better than other types of filtering
#    approaches, such as moving averages techhniques.
#    Parameters
#    ----------
#    y : array_like, shape (N,)
#        the values of the time history of the signal.
#    window_size : int
#        the length of the window. Must be an odd integer number.
#    order : int
#        the order of the polynomial used in the filtering.
#        Must be less then `window_size` - 1.
#    deriv: int
#        the order of the derivative to compute (default = 0 means only smoothing)
#    Returns
#    -------
#    ys : ndarray, shape (N)
#        the smoothed signal (or it's n-th derivative).
#    Notes
#    -----
#    The Savitzky-Golay is a type of low-pass filter, particularly
#    suited for smoothing noisy data. The main idea behind this
#    approach is to make for each point a least-square fit with a
#    polynomial of high order over a odd-sized window centered at
#    the point.
#    Examples
#    --------
#    t = np.linspace(-4, 4, 500)
#    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
#    ysg = savitzky_golay(y, window_size=31, order=4)
#    import matplotlib.pyplot as plt
#    plt.plot(t, y, label='Noisy signal')
#    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
#    plt.plot(t, ysg, 'r', label='Filtered signal')
#    plt.legend()
#    plt.show()
#    References
#    ----------
#    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
#       Data by Simplified Least Squares Procedures. Analytical
#       Chemistry, 1964, 36 (8), pp 1627-1639.
#    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
#       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
#       Cambridge University Press ISBN-13: 9780521880688
#    """
    cdef:
        int half_window
        # precompute coefficients
        np.ndarray[np.int64_t, ndim=2] b
        # pad the signal at the extremes with
        # values taken from the signal itself
        np.ndarray[np.float64_t, ndim=1] y = signal.astype('float64')
        np.ndarray[np.float64_t, ndim=1] m, firstvals, lastvals, ret
    if window_size % 2 != 1: window_size += 1
    half_window = (window_size - 1) / 2
    # precompute coefficients
    b = np.mat([[1, k, k**2, k**3] for k in range(-half_window, half_window+1)],
               dtype='int64')
    m = np.linalg.pinv(b).A[deriv].astype('float64')
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y.astype('float64'), lastvals))
    ret = np.convolve( m, y, mode='valid')
    return ret.astype('float32')

