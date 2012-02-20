"""Data smoothing and local maximum functions from SciPy cookbook


"""
import numpy as np

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'('+str(window_len)+')')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def local_maximum(x):
    index_array = range(len(x))
    f = np.r_[True, x[1:] > x[:-1]] & np.r_[x[:-1] > x[1:], True]
    return (index_array[f],x[f])
    

# the following timeseries codes are from http://www.ieap.uni-kiel.de/et/people/wimmer/teaching/et2/time_series_1.py

"""Try smoothing data with numpy.

This uses the idea that convolving data with a response function
is actually nothing else but smoothing data.


Date: December 7, 2008

Author: Bob Wimmer

"""
class timeseries(object):
    """Time series given by x vlaues (e.g. time) and y values (e.g. temperature).
    Perform various operations on time series wuc as:
    - treat/replace missing values
    - smoothing
    - FFT
    - and more to come (maybe)

    Note: requires numpy

    """
    def __init__(self,x,y):
        if len(x) < 2:
            print "this is not a time series, but a single point!"
        if len(y) != len(x):
            print "len(x) != len(y) -- this needs to be fixed!"

        self.x = x
        self.y = y
        self.beginning = x[0]
        self.end = x[-1]
        
    
    def max(self):
        """returns the maximum of x, y"""
        return self.x.max(), self.y.max()

    def min(self):
        """returns the minimum of x, y"""
        return self.x.min(), self.y.min()

    def cut_ind(self, lo,hi):
        """return subset of time series between lo and hi indices."""
        self.xci = self.x[lo:hi]
        self.yci = self.y[lo:hi]
        return self.xci, self.yci

    def cut_x(self, x_lo,x_hi):
        """return subset of time series between x_lo and x_hi values."""
        """find indices of x_lo and x_hi, then call cut_ind"""
        lo = self.x.searchsorted(x_lo)
        hi = self.x.searchsorted(x_hi)
        print lo, self.x[lo], self.x[lo-1]
        self.xcx = self.x[lo:hi]
        self.ycx = self.y[lo:hi]
        return self.xcx, self.ycx

    def clean(self):
        """remove and/or linearly interpolate appropriately missing data values"""
        """This is not yet implemented"""

    def fft(self):
        self.foury = np.fft.fft(self.y)
        N = len(self.y)
        timestep = self.x[1] - self.x[0]        # if unit=day -> freq unit=cycles/day
        self.fourx = np.fft.fftfreq(N, d=timestep)
        self.fourx = np.fft.fftshift(self.fourx)
        self.foury = np.fft.fftshift(self.foury)
        return self.fourx, self.foury
        
        
    def smooth(self, x,window_len=10,window='hanning'):
        #print "in smooth method"
        #print x, window_len, window
        self.s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
        print(len(self.s))
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('numpy.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),self.s,mode='same')
        dy = self.s.copy()
        dy[window_len-1:-window_len+1] = self.s[window_len-1:-window_len+1] - y[window_len-1:-window_len+1]
        return y[window_len-1:-window_len+1], dy[window_len-1:-window_len+1]

