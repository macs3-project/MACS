# smoothing function 
import numpy as np
cimport numpy as np
from cpython cimport bool
from math import sqrt as mathsqrt

cpdef maxima(np.ndarray[np.float32_t, ndim=1] signal,
             int window_size=51):
    """return the local maxima in a signal after applying a 2nd order
    Savitsky-Golay (polynomial) filter using window_size specified  
    """
    cdef np.ndarray[np.int32_t, ndim=1] m = np.where(np.diff(np.sign(savitzky_golay_order2(signal, window_size, deriv=1))) == 2)[0].astype('int32')
    return m

cdef np.ndarray[np.int32_t, ndim=1] internal_minima(np.ndarray[np.float32_t, ndim=1] signal,
                                                      np.ndarray[np.int32_t, ndim=1] maxima):
    cdef:
        np.ndarray[np.int32_t, ndim=1] ret
        int n = maxima.shape[0]
        int i, v, v2
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

cdef inline float sqrt(float threshold):
    return mathsqrt(threshold)
    
cpdef enforce_peakyness(np.ndarray[np.float32_t, ndim=1] signal,
                        np.ndarray[np.int32_t, ndim=1] maxima):
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
        np.ndarray[np.int32_t, ndim=1] minima = internal_minima(signal, maxima)
        np.ndarray[np.float32_t, ndim=1] new_signal
        int n = minima.shape[0]
        float threshold
        np.ndarray[np.int32_t, ndim=1] peaky_maxima = maxima.copy()
        int j = 0
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
cdef bool is_valid_peak(np.ndarray[np.float32_t, ndim=1] signal, int maximum):
    cdef:
        s = hard_clip(signal, maximum)
        int length = s.shape[0]
    if length < 50: return False
    elif too_flat(s): return False
    return True

# require at least 6 different float values -- prevents broad flat peaks
cdef bool too_flat(np.ndarray[np.float32_t, ndim=1] signal):
#    """return whether signal has at least 6 unique values
#    """
    return np.unique(signal).shape[0] < 6

# hard clip a region with negative values
cdef np.ndarray[np.float32_t, ndim=1] hard_clip(np.ndarray[np.float32_t, ndim=1] signal, int maximum):
#    """clip the signal in both directions at the nearest values <= 0
#    to position maximum
#    """
    cdef:
        int i
        int left = 0
        int right = signal.shape[0]
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

# for all maxima, set min_subpeak_width = 0
#cpdef peak_maxima(np.ndarray[np.float32_t, ndim=1] signal,
#             int window_size=51, int min_subpeak_width=5):
#    cdef:
#        np.ndarray[np.float32_t, ndim=1] D = savitzky_golay_order2(signal, window_size, deriv=1)
#        np.ndarray[np.int32_t, ndim=1] m = np.where(np.diff(np.sign(D)) == 2)[0].astype('int32') + 1
#        np.ndarray[np.int32_t, ndim=1] n = np.zeros_like(m)
#        int i_max
#        int halfw = (min_subpeak_width - 1) / 2
#        int signalw = D.shape[0]
#        int i_new = 0
#        int pos, start, end
#    if m.shape[0] == 0: return m
#    elif m.shape[0] == 1: return m
#    else:
#        i_max = np.where(signal[m] == signal[m].max())[0][0]
#    for pos in m:
#        start = max(pos - halfw, 0)
#        end = min(pos + halfw, signalw)
#        if ((D[start:pos] >= 0).all() and (D[(pos + 1):end] <= 0).all()):
#            n[i_new] = pos
#            i_new += 1
#    if i_new == 0:
#        return m[i_max:(i_max + 1)]
#    else:
#        n.resize(i_new, refcheck=False)
#        return n

cpdef enforce_valleys(np.ndarray[np.float32_t, ndim=1] signal,
                      np.ndarray[np.int32_t, ndim=1] summits,
                      float min_valley = 0.8):
    """require a value of <= min_valley * lower summit
    between each pair of summits
    """
    cdef:
        float req_min, v, prev_v
        int summit_pos, prev_summit_pos
        int n_summits = summits.shape[0]
        int n_valid_summits = 1
        np.ndarray[np.int32_t, ndim=1] valid_summits = summits.copy()
    # Step 1: Remove peaks that do not have sufficient valleys
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
    # Step 2: Re-find peaks from subtracted signal
    # 
    return valid_summits
    
        

# Modified from http://www.scipy.org/Cookbook/SavitzkyGolay
# positive window_size not enforced anymore
# needs sane input paramters, window size > 4
# switched to double precision for internal accuracy
cpdef savitzky_golay_order2(np.ndarray[np.float32_t, ndim=1] signal,
                     int window_size, int deriv=0):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
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
        int half_window, k
        np.ndarray[np.int64_t, ndim=2] b
        # pad the signal at the extremes with
        # values taken from the signal itself
        np.ndarray[np.float32_t, ndim=1] firstvals, lastvals, ret
        np.ndarray[np.float64_t, ndim=1] m
    if window_size % 2 != 1: window_size += 1
    half_window = (window_size - 1) / 2
    # precompute coefficients
    b = np.mat([[1, k, k**2] for k in range(-half_window, half_window+1)],
               dtype='int64')
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
    lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
    signal = np.concatenate((firstvals, signal, lastvals))
    ret = np.convolve( m, signal.astype('float64'), mode='valid').astype('float32')
    return ret

