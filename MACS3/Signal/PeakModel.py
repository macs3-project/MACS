# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-15 10:20:32 Tao Liu>
"""Module Description: Build shifting model

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# Python modules
# ------------------------------------

# ------------------------------------
# MACS3 modules
# ------------------------------------
# from MACS3.Utilities.Constants import *
from MACS3.Signal.Pileup import naive_quick_pileup, naive_call_peaks

# ------------------------------------
# Other modules
# ------------------------------------
import cython
from cython.cimports.cpython import bool
import numpy as np
import cython.cimports.numpy as cnp

# ------------------------------------
# C lib
# ------------------------------------


class NotEnoughPairsException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


@cython.cclass
class PeakModel:
    """Peak Model class.
    """
    # this can be PETrackI or FWTrack
    treatment: object
    # genome size
    gz: cython.double
    max_pairnum: cython.int
    umfold: cython.int
    lmfold: cython.int
    bw: cython.int
    d_min: cython.int
    tag_expansion_size: cython.int

    info: object
    debug: object
    warn: object
    error: object

    summary: str
    max_tags: cython.int
    peaksize: cython.int

    plus_line = cython.declare(cnp.ndarray, visibility="public")
    minus_line = cython.declare(cnp.ndarray, visibility="public")
    shifted_line = cython.declare(cnp.ndarray, visibility="public")
    xcorr = cython.declare(cnp.ndarray, visibility="public")
    ycorr = cython.declare(cnp.ndarray, visibility="public")

    d = cython.declare(cython.int, visibility="public")
    scan_window = cython.declare(cython.int, visibility="public")
    min_tags = cython.declare(cython.int, visibility="public")
    alternative_d = cython.declare(list, visibility="public")

    def __init__(self, opt, treatment, max_pairnum: cython.int = 500):
        # , double gz = 0, int umfold=30, int lmfold=10, int bw=200,
        # int ts = 25, int bg=0, bool quiet=False):
        self.treatment = treatment
        self.gz = opt.gsize
        self.umfold = opt.umfold
        self.lmfold = opt.lmfold
        # opt.tsize| test 10bps. The reason is that we want the best
        # 'lag' between left & right cutting sides. A tag will be
        # expanded to 10bps centered at cutting point.
        self.tag_expansion_size = 10
        # discard any predicted fragment sizes < d_min
        self.d_min = opt.d_min
        self.bw = opt.bw
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn
        self.error = opt.warn
        self.max_pairnum = max_pairnum

    @cython.ccall
    def build(self):
        """Build the model. Main function of PeakModel class.

        1. prepare self.d, self.scan_window, self.plus_line,
        self.minus_line and self.shifted_line.

        2. find paired + and - strand peaks

        3. find the best d using x-correlation
        """
        paired_peakpos: dict
        num_paired_peakpos: cython.long
        c: bytes                # chromosome

        self.peaksize = 2*self.bw
        # mininum unique hits on single strand, decided by lmfold
        self.min_tags = int(round(float(self.treatment.total) *
                                  self.lmfold *
                                  self.peaksize / self.gz / 2))
        # maximum unique hits on single strand, decided by umfold
        self.max_tags = int(round(float(self.treatment.total) *
                                  self.umfold *
                                  self.peaksize / self.gz / 2))
        self.debug(f"#2 min_tags: {self.min_tags}; max_tags:{self.max_tags}; ")
        self.info("#2 looking for paired plus/minus strand peaks...")
        # find paired + and - strand peaks
        paired_peakpos = self.__find_paired_peaks()

        num_paired_peakpos = 0
        for c in list(paired_peakpos.keys()):
            num_paired_peakpos += len(paired_peakpos[c])

        self.info("#2 Total number of paired peaks: %d" % (num_paired_peakpos))

        if num_paired_peakpos < 100:
            self.error(f"#2 MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only find {num_paired_peakpos}! Please make your MFOLD range broader and try again. If MACS3 still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead.")
            self.error("#2 Process for pairing-model is terminated!")
            raise NotEnoughPairsException("No enough pairs to build model")

        # build model, find the best d using cross-correlation
        self.__paired_peak_model(paired_peakpos)

    def __str__(self):
        """For debug...

        """
        return """
Summary of Peak Model:
  Baseline: %d
  Upperline: %d
  Fragment size: %d
  Scan window size: %d
""" % (self.min_tags, self.max_tags, self.d, self.scan_window)

    @cython.cfunc
    def __find_paired_peaks(self) -> dict:
        """Call paired peaks from fwtrackI object.

        Return paired peaks center positions.
        """
        i: cython.int
        chrs: list
        chrom: bytes
        plus_tags: cnp.ndarray(cython.int, ndim=1)
        minus_tags: cnp.ndarray(cython.int, ndim=1)
        plus_peaksinfo: list
        minus_peaksinfo: list
        paired_peaks_pos: dict  # return

        chrs = list(self.treatment.get_chr_names())
        chrs.sort()
        paired_peaks_pos = {}
        for i in range(len(chrs)):
            chrom = chrs[i]
            self.debug(f"Chromosome: {chrom}")
            # extract tag positions
            [plus_tags, minus_tags] = self.treatment.get_locations_by_chr(chrom)
            # look for + strand peaks
            plus_peaksinfo = self.__naive_find_peaks(plus_tags)
            self.debug("Number of unique tags on + strand: %d" % (plus_tags.shape[0]))
            self.debug("Number of peaks in + strand: %d" % (len(plus_peaksinfo)))
            if plus_peaksinfo:
                self.debug(f"plus peaks: first - {plus_peaksinfo[0]} ... last - {plus_peaksinfo[-1]}")
            # look for - strand peaks
            minus_peaksinfo = self.__naive_find_peaks(minus_tags)
            self.debug("Number of unique tags on - strand: %d" % (minus_tags.shape[0]))
            self.debug("Number of peaks in - strand: %d" % (len(minus_peaksinfo)))
            if minus_peaksinfo:
                self.debug(f"minus peaks: first - {minus_peaksinfo[0]} ... last - {minus_peaksinfo[-1]}")
            if not plus_peaksinfo or not minus_peaksinfo:
                self.debug("Chrom %s is discarded!" % (chrom))
                continue
            else:
                paired_peaks_pos[chrom] = self.__find_pair_center(plus_peaksinfo, minus_peaksinfo)
                self.debug("Number of paired peaks in this chromosome: %d" % (len(paired_peaks_pos[chrom])))
        return paired_peaks_pos

    @cython.cfunc
    def __naive_find_peaks(self,
                           taglist: cnp.ndarray(cython.int, ndim=1)) -> list:
        """Naively call peaks based on tags counting.

        Return peak positions and the tag number in peak region by a tuple list[(pos,num)].
        """
        peak_info: list
        pileup_array: list

        # store peak pos in every peak region and unique tag number in
        # every peak region
        peak_info = []

        # less than 2 tags, no need to call peaks, return []
        if taglist.shape[0] < 2:
            return peak_info

        # build pileup by extending both side to half peak size
        pileup_array = naive_quick_pileup(taglist, int(self.peaksize/2))
        peak_info = naive_call_peaks(pileup_array,
                                     self.min_tags,
                                     self.max_tags)

        return peak_info

    @cython.cfunc
    def __paired_peak_model(self, paired_peakpos: dict):
        """Use paired peak positions and treatment tag positions to
        build the model.

        Modify self.(d, model_shift size and scan_window size. and
        extra, plus_line, minus_line and shifted_line for plotting).

        """
        window_size: cython.int
        i: cython.int
        chroms: list
        paired_peakpos_chrom: object

        tags_plus: cnp.ndarray(cython.int, ndim=1)
        tags_minus: cnp.ndarray(cython.int, ndim=1)
        plus_start: cnp.ndarray(cython.int, ndim=1)
        plus_end: cnp.ndarray(cython.int, ndim=1)
        minus_start: cnp.ndarray(cython.int, ndim=1)
        minus_end: cnp.ndarray(cython.int, ndim=1)
        plus_line: cnp.ndarray(cython.int, ndim=1)
        minus_line: cnp.ndarray(cython.int, ndim=1)

        plus_data: cnp.ndarray
        minus_data: cnp.ndarray
        xcorr: cnp.ndarray
        ycorr: cnp.ndarray
        i_l_max: cnp.ndarray

        window_size = 1+2*self.peaksize+self.tag_expansion_size
        # for plus strand pileup
        self.plus_line = np.zeros(window_size, dtype="i4")
        # for minus strand pileup
        self.minus_line = np.zeros(window_size, dtype="i4")
        # for fast pileup
        plus_start = np.zeros(window_size, dtype="i4")
        # for fast pileup
        plus_end = np.zeros(window_size, dtype="i4")
        # for fast pileup
        minus_start = np.zeros(window_size, dtype="i4")
        # for fast pileup
        minus_end = np.zeros(window_size, dtype="i4")
        self.debug("start model_add_line...")
        chroms = list(paired_peakpos.keys())

        for i in range(len(chroms)):
            paired_peakpos_chrom = paired_peakpos[chroms[i]]
            (tags_plus, tags_minus) = self.treatment.get_locations_by_chr(chroms[i])
            # every paired peak has plus line and minus line
            self.__model_add_line(paired_peakpos_chrom,
                                  tags_plus,
                                  plus_start,
                                  plus_end)
            self.__model_add_line(paired_peakpos_chrom,
                                  tags_minus,
                                  minus_start,
                                  minus_end)

        self.__count(plus_start, plus_end, self.plus_line)
        self.__count(minus_start, minus_end, self.minus_line)

        self.debug("start X-correlation...")
        # Now I use cross-correlation to find the best d
        plus_line = self.plus_line
        minus_line = self.minus_line

        # normalize first
        minus_data = (minus_line - minus_line.mean())/(minus_line.std()*len(minus_line))
        plus_data = (plus_line - plus_line.mean())/(plus_line.std()*len(plus_line))

        # cross-correlation
        ycorr = np.correlate(minus_data, plus_data, mode="full")[window_size-self.peaksize:window_size+self.peaksize]
        xcorr = np.linspace(len(ycorr)//2*-1, len(ycorr)//2, num=len(ycorr))

        # smooth correlation values to get rid of local maximums from small fluctuations.
        # window size is by default 11.
        ycorr = smooth(ycorr, window="flat")

        # all local maximums could be alternative ds.
        i_l_max = np.r_[False, ycorr[1:] > ycorr[:-1]] & np.r_[ycorr[:-1] > ycorr[1:], False]
        i_l_max = np.where(i_l_max)[0]
        i_l_max = i_l_max[xcorr[i_l_max] > self.d_min]
        i_l_max = i_l_max[np.argsort(ycorr[i_l_max])[::-1]]

        self.alternative_d = sorted([int(x) for x in xcorr[i_l_max]])
        assert len(self.alternative_d) > 0, "No proper d can be found! Tweak --mfold?"

        self.d = xcorr[i_l_max[0]]

        self.ycorr = ycorr
        self.xcorr = xcorr

        self.scan_window = max(self.d, self.tag_expansion_size)*2

        self.info("#2 Model building with cross-correlation: Done")

        return True

    @cython.cfunc
    def __model_add_line(self,
                         pos1: list,
                         pos2: cnp.ndarray(cython.int, ndim=1),
                         start: cnp.ndarray(cython.int, ndim=1),
                         end: cnp.ndarray(cython.int, ndim=1)):
        """Project each pos in pos2 which is included in
        [pos1-self.peaksize,pos1+self.peaksize] to the line.

        pos1: paired centers -- list of coordinates
        pos2: tags of certain strand -- a numpy.array object
        line: numpy array object where we pileup tags

        """
        i1: cython.int
        i2: cython.int
        i2_prev: cython.int
        i1_max: cython.int
        i2_max: cython.int
        last_p2: cython.int
        psize_adjusted1: cython.int
        p1: cython.int
        p2: cython.int
        max_index: cython.int
        s: cython.int
        e: cython.int

        i1 = 0                  # index for pos1
        i2 = 0                  # index for pos2 index for pos2 in
        # previous pos1 [pos1-self.peaksize,pos1+self.peaksize] region
        i2_prev = 0
        i1_max = len(pos1)
        i2_max = pos2.shape[0]
        flag_find_overlap = False

        max_index = start.shape[0] - 1

        # half window
        psize_adjusted1 = self.peaksize + self.tag_expansion_size // 2

        while i1 < i1_max and i2 < i2_max:
            p1 = pos1[i1]
            p2 = pos2[i2]

            if p1-psize_adjusted1 > p2:
                # move pos2
                i2 += 1
            elif p1+psize_adjusted1 < p2:
                # move pos1
                i1 += 1
                i2 = i2_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    # only the first index is recorded
                    i2_prev = i2
                # project
                s = max(int(p2-self.tag_expansion_size/2-p1+psize_adjusted1), 0)
                start[s] += 1
                e = min(int(p2+self.tag_expansion_size/2-p1+psize_adjusted1), max_index)
                end[e] -= 1
                i2 += 1
        return

    @cython.cfunc
    def __count(self,
                start: cnp.ndarray(cython.int, ndim=1),
                end: cnp.ndarray(cython.int, ndim=1),
                line: cnp.ndarray(cython.int, ndim=1)):
        """
        """
        i: cython.int
        pileup: cython.long

        pileup = 0
        for i in range(line.shape[0]):
            pileup += start[i] + end[i]
            line[i] = pileup
        return

    @cython.cfunc
    def __find_pair_center(self,
                           pluspeaks: list,
                           minuspeaks: list):
        # index for plus peaks
        ip: cython.long = 0
        # index for minus peaks
        im: cython.long = 0
        # index for minus peaks in previous plus peak
        im_prev: cython.long = 0
        pair_centers: list
        ip_max: cython.long
        im_max: cython.long
        flag_find_overlap: bool
        pp: cython.int
        mp: cython.int
        pn: cython.float
        mn: cython.float

        pair_centers = []
        ip_max = len(pluspeaks)
        im_max = len(minuspeaks)
        self.debug(f"ip_max: {ip_max}; im_max: {im_max}")
        flag_find_overlap = False
        while ip < ip_max and im < im_max:
            # for (peakposition, tagnumber in peak)
            (pp, pn) = pluspeaks[ip]
            (mp, mn) = minuspeaks[im]
            if pp-self.peaksize > mp:
                # move minus
                im += 1
            elif pp+self.peaksize < mp:
                # move plus
                ip += 1
                im = im_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    # only the first index is recorded
                    im_prev = im
                # number tags in plus and minus peak region are comparable...
                if pn/mn < 2 and pn/mn > 0.5:
                    if pp < mp:
                        pair_centers.append((pp+mp)//2)
                im += 1
        if pair_centers:
            self.debug(f"Paired centers: first - {pair_centers[0]} ... second - {pair_centers[-1]} ")
        return pair_centers


# smooth function from SciPy cookbook:
# http://www.scipy.org/Cookbook/SignalSmooth
@cython.ccall
def smooth(x,
           window_len: cython.int = 11,
           window: str = 'hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be
                    an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming',
                'bartlett', 'blackman' flat window will produce a
                moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array
          instead of a string

    NOTE: length(output) != length(input), to correct this: return
          y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'flat':        # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[(window_len//2):-(window_len//2)]