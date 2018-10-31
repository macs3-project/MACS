# Time-stamp: <2015-03-10 15:38:17 Tao Liu>

"""Module Description

Copyright (c) 2008,2009 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
import sys, time, random
import numpy as np
cimport numpy as np
from array import array
from MACS2.Constants import *

from cpython cimport bool
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
ctypedef np.float32_t float32_t

cpdef median (nums):
    """Calculate Median.

    Parameters:
    nums:  list of numbers
    Return Value:
    median value
    """
    p = sorted(nums)
    l = len(p)
    if l%2 == 0:
        return (p[l/2]+p[l/2-1])/2
    else:
        return p[l/2]

class NotEnoughPairsException(Exception):
    def __init__ (self,value):
        self.value = value
    def __str__ (self):
        return repr(self.value)

cdef class PeakModel:
    """Peak Model class.
    """
    cdef:
        object treatment
        double gz
        int max_pairnum
        int umfold
        int lmfold
        int bw
        int d_min
        int tag_expansion_size
        object info, debug, warn, error
        str summary
        public np.ndarray plus_line, minus_line, shifted_line
        public int d
        public int scan_window
        public int min_tags
        int max_tags
        int peaksize
        public list alternative_d
        public np.ndarray xcorr, ycorr

    def __init__ ( self, opt , treatment, int max_pairnum=500 ): #, double gz = 0, int umfold=30, int lmfold=10, int bw=200, int ts = 25, int bg=0, bool quiet=False):
        self.treatment = treatment
        #if opt:
        self.gz = opt.gsize
        self.umfold = opt.umfold
        self.lmfold = opt.lmfold
        self.tag_expansion_size = 10         #opt.tsize| test 10bps. The reason is that we want the best 'lag' between left & right cutting sides. A tag will be expanded to 10bps centered at cutting point.
        self.d_min = 20 ### discard any fragment size < d_min
        self.bw = opt.bw
        self.info  = opt.info
        self.debug = opt.debug
        self.warn  = opt.warn
        self.error = opt.warn
        #else:
        #    self.gz = gz
        #    self.umfold = umfold
        #    self.lmfold = lmfold            
        #    self.tag_expansion_size = ts
        #    self.bg = bg
        #    self.bw = bw
        #    self.info  = lambda x: sys.stderr.write(x+"\n")
        #    self.debug = lambda x: sys.stderr.write(x+"\n")
        #    self.warn  = lambda x: sys.stderr.write(x+"\n")
        #    self.error = lambda x: sys.stderr.write(x+"\n")
        #if quiet:
        #    self.info = lambda x: None
        #    self.debug = lambda x: None
        #    self.warn = lambda x: None
        #    self.error = lambda x: None
            
        self.max_pairnum = max_pairnum
        #self.summary = ""
        #self.plus_line = None
        #self.minus_line = None
        #self.shifted_line = None
        #self.d = None
        #self.scan_window = None
        #self.min_tags = None
        #self.peaksize = None
        self.build()
    
    cpdef build (self):
        """Build the model.

        prepare self.d, self.scan_window, self.plus_line,
        self.minus_line and self.shifted_line to use.
        """
        cdef:
            dict paired_peaks
            long num_paired_peakpos, num_paired_peakpos_remained, num_paired_peakpos_picked
            str c
        

        self.peaksize = 2*self.bw
        self.min_tags = int(round(float(self.treatment.total) * self.lmfold * self.peaksize / self.gz /2)) # mininum unique hits on single strand
        self.max_tags = int(round(float(self.treatment.total) * self.umfold * self.peaksize / self.gz /2)) # maximum unique hits on single strand
        #print self.min_tags, self.max_tags
        #print self.min_tags
        #print self.max_tags
        # use treatment data to build model
        self.info("#2 looking for paired plus/minus strand peaks...")
        paired_peakpos = self.__paired_peaks ()
        # select up to 1000 pairs of peaks to build model
        num_paired_peakpos = 0
        num_paired_peakpos_remained = self.max_pairnum
        num_paired_peakpos_picked = 0
        # select only num_paired_peakpos_remained pairs.
        for c in paired_peakpos.keys():
            num_paired_peakpos +=len(paired_peakpos[c])
        # TL: Now I want to use everything

        num_paired_peakpos_picked = num_paired_peakpos

        self.info("#2 number of paired peaks: %d" % (num_paired_peakpos))
        if num_paired_peakpos < 100:
            self.error("Too few paired peaks (%d) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead." % (num_paired_peakpos))
            self.error("Process for pairing-model is terminated!")
            raise NotEnoughPairsException("No enough pairs to build model")
        elif num_paired_peakpos < self.max_pairnum:
            self.warn("Fewer paired peaks (%d) than %d! Model may not be build well! Lower your MFOLD parameter may erase this warning. Now I will use %d pairs to build model!" % (num_paired_peakpos,self.max_pairnum,num_paired_peakpos_picked))
        self.debug("Use %d pairs to build the model." % (num_paired_peakpos_picked))
        self.__paired_peak_model(paired_peakpos)

    def __str__ (self):
        """For debug...

        """
        return """
Summary of Peak Model:
  Baseline: %d
  Upperline: %d
  Fragment size: %d
  Scan window size: %d
""" % (self.min_tags,self.max_tags,self.d,self.scan_window)

    cdef __paired_peak_model (self, paired_peakpos,):
        """Use paired peak positions and treatment tag positions to build the model.

        Modify self.(d, model_shift size and scan_window size. and extra, plus_line, minus_line and shifted_line for plotting).
        """
        cdef:
            int window_size, i
            list chroms
            object paired_peakpos_chrom
            np.ndarray[np.int32_t, ndim=1] tags_plus, tags_minus, plus_start, plus_end, minus_start, minus_end, plus_line, minus_line
            np.ndarray plus_data, minus_data, xcorr, ycorr, i_l_max
        
        window_size = 1+2*self.peaksize+self.tag_expansion_size
        self.plus_line = np.zeros(window_size, dtype="int32") # for plus strand pileup
        self.minus_line = np.zeros(window_size, dtype="int32")# for minus strand pileup
        plus_start = np.zeros(window_size, dtype="int32")     # for fast pileup
        plus_end = np.zeros(window_size, dtype="int32")       # for fast pileup
        minus_start = np.zeros(window_size, dtype="int32")    # for fast pileup
        minus_end = np.zeros(window_size, dtype="int32")      # for fast pileup
        #self.plus_line = [0]*window_size
        #self.minus_line = [0]*window_size        
        self.info("start model_add_line...")
        chroms = paired_peakpos.keys()
        
        for i in range(len(chroms)):
            paired_peakpos_chrom = paired_peakpos[chroms[i]]
            (tags_plus, tags_minus) = self.treatment.get_locations_by_chr(chroms[i])
            # every paired peak has plus line and minus line
            #  add plus_line
            #self.plus_line = self.__model_add_line (paired_peakpos_chrom, tags_plus, self.plus_line) #, plus_strand=1)
            self.__model_add_line (paired_peakpos_chrom, tags_plus, plus_start, plus_end) #, plus_strand=1)
            self.__model_add_line (paired_peakpos_chrom, tags_minus, minus_start, minus_end) #, plus_strand=0)
            #  add minus_line
            #self.minus_line = self.__model_add_line (paired_peakpos_chrom, tags_minus, self.minus_line) #, plus_strand=0)

        self.__count ( plus_start, plus_end, self.plus_line )
        self.__count ( minus_start, minus_end, self.minus_line )

        self.info("start X-correlation...")
        # Now I use cross-correlation to find the best d
        #plus_line = np.asarray(self.plus_line,dtype="int32")
        #minus_line = np.asarray(self.minus_line,dtype="int32")
        plus_line = self.plus_line
        minus_line = self.minus_line
        
        # normalize first
        minus_data = (minus_line - minus_line.mean())/(minus_line.std()*len(minus_line))
        plus_data = (plus_line - plus_line.mean())/(plus_line.std()*len(plus_line))
        #print "plus:",len(plus_data)
        #print "minus:",len(minus_data)

        # cross-correlation
        ycorr = np.correlate(minus_data,plus_data,mode="full")[window_size-self.peaksize:window_size+self.peaksize]
        xcorr = np.linspace(len(ycorr)//2*-1, len(ycorr)//2, num=len(ycorr))

        # smooth correlation values to get rid of local maximums from small fluctuations. 
        ycorr = smooth(ycorr, window="flat") # window size is by default 11.

        # all local maximums could be alternative ds.
        i_l_max = np.r_[False, ycorr[1:] > ycorr[:-1]] & np.r_[ycorr[:-1] > ycorr[1:], False]
        i_l_max = np.where(i_l_max)[0]
        i_l_max = i_l_max[ xcorr[i_l_max] > self.d_min ]
        i_l_max = i_l_max[ np.argsort(ycorr[i_l_max])[::-1]]
#         filter(lambda i: xcorr[i]>self.d_min, i_l_max )
#         i_l_max = sorted(i_l_max,
#                          key=ycorr.__getitem__, 
#                          reverse=True)
        self.alternative_d = sorted(map(int, xcorr[i_l_max]))
        assert len(self.alternative_d) > 0, "No proper d can be found! Tweak --mfold?"
        
        self.d = xcorr[i_l_max[0]]
#         i_l_max = filter(lambda: ycorr)
#         tmp_cor_alternative_d = ycorr[ i_l_max ]
#         tmp_alternative_d = xcorr[ i_l_max ]
#         cor_alternative_d =  tmp_cor_alternative_d [ tmp_alternative_d > 0 ]
#         self.alternative_d = map( int, tmp_alternative_d[ tmp_alternative_d > 0 ] )
        
        
        # best cross-correlation point
        
#         d_cand = xcorr[ np.where( ycorr == sorted( cor_alternative_d [::-1] ) ) ]
#         print (cor_alternative_d)
#         print (d_cand)
#         d_cand = xcorr[ np.where( ycorr== max( cor_alternative_d ) )[0] ]
        
        #### make sure fragment size is not zero
        
#         self.d = [ x for x in d_cand if x > self.d_min ] [0] 
        #self.d = xcorr[np.where(ycorr==max(ycorr))[0][0]] #+self.tag_expansion_size

        # get rid of the last local maximum if it's at the right end of curve.
        
        
        self.ycorr = ycorr
        self.xcorr = xcorr

        #shift_size = self.d/2
        
        self.scan_window = max(self.d,self.tag_expansion_size)*2
        # a shifted model
        #self.shifted_line = [0]*window_size

        self.info("end of X-cor")
        
        return True

    cdef __model_add_line (self, object pos1, np.ndarray[np.int32_t, ndim=1] pos2, np.ndarray[np.int32_t, ndim=1] start, np.ndarray[np.int32_t, ndim=1] end): #, int plus_strand=1):
        """Project each pos in pos2 which is included in
        [pos1-self.peaksize,pos1+self.peaksize] to the line.

        pos1: paired centers -- array.array
        pos2: tags of certain strand -- a numpy.array object
        line: numpy array object where we pileup tags

        """
        cdef:
            int i1, i2, i2_prev, i1_max, i2_max, last_p2, psize_adjusted1, psize_adjusted2, p1, p2, max_index, s, e
        
        i1 = 0                  # index for pos1
        i2 = 0                  # index for pos2
        i2_prev = 0             # index for pos2 in previous pos1
                                # [pos1-self.peaksize,pos1+self.peaksize]
                                # region
        i1_max = len(pos1)
        i2_max = pos2.shape[0]
        last_p2 = -1
        flag_find_overlap = False

        max_index = start.shape[0] - 1

        psize_adjusted1 = self.peaksize + self.tag_expansion_size / 2 # half window

        while i1<i1_max and i2<i2_max:
            p1 = pos1[i1]
            #if plus_strand:
            #    p2 = pos2[i2]
            #else:
            #    p2 = pos2[i2] - self.tag_expansion_size

            p2 = pos2[i2] #- self.tag_expansion_size/2
                
            if p1-psize_adjusted1 > p2: # move pos2
                i2 += 1
            elif p1+psize_adjusted1 < p2: # move pos1
                i1 += 1                 
                i2 = i2_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    i2_prev = i2 # only the first index is recorded
                # project
                #for i in range(p2-p1+self.peaksize,p2-p1+self.peaksize+self.tag_expansion_size):
                s = max(p2-self.tag_expansion_size/2-p1+psize_adjusted1, 0)
                start[s] += 1
                e = min(p2+self.tag_expansion_size/2-p1+psize_adjusted1, max_index)
                end[e] -= 1
                #line[s:e] += 1
                #for i in range(s,e):
                #    #if i>=0 and i<length_l:
                #    line[i]+=1
                i2+=1
        return
    
    cdef __count ( self, np.ndarray[np.int32_t, ndim=1] start, np.ndarray[np.int32_t, ndim=1] end, np.ndarray[np.int32_t, ndim=1] line ):
        """
        """
        cdef:
            int i
            long pileup
        pileup = 0
        for i in range(line.shape[0]):
            pileup += start[i] + end[i]
            line[i] = pileup
        return


    cdef __paired_peaks (self):
        """Call paired peaks from fwtrackI object.

        Return paired peaks center positions.
        """
        cdef:
           int i
           list chrs
           str chrom
           dict paired_peaks_pos
           np.ndarray[np.int32_t, ndim=1] plus_tags, minus_tags

        chrs = self.treatment.get_chr_names()
        chrs.sort()
        paired_peaks_pos = {}
        for i in range( len(chrs) ):
            chrom = chrs[ i ]
            self.debug("Chromosome: %s" % (chrom) )
            [ plus_tags, minus_tags ] = self.treatment.get_locations_by_chr( chrom )
            plus_peaksinfo = self.__naive_find_peaks ( plus_tags, 1 )
            self.debug("Number of unique tags on + strand: %d" % ( plus_tags.shape[0] ) )
            self.debug("Number of peaks in + strand: %d" % ( len(plus_peaksinfo) ) )
            minus_peaksinfo = self.__naive_find_peaks ( minus_tags, 0 )
            self.debug("Number of unique tags on - strand: %d" % ( minus_tags.shape[0] ) )
            self.debug("Number of peaks in - strand: %d" % ( len( minus_peaksinfo ) ) )
            if not plus_peaksinfo or not minus_peaksinfo:
                self.debug("Chrom %s is discarded!" % (chrom) )
                continue
            else:
                paired_peaks_pos[chrom] = self.__find_pair_center (plus_peaksinfo, minus_peaksinfo)
                self.debug("Number of paired peaks: %d" %(len(paired_peaks_pos[chrom])))
        return paired_peaks_pos

    cdef __find_pair_center (self, pluspeaks, minuspeaks):
        ip = 0                  # index for plus peaks
        im = 0                  # index for minus peaks
        im_prev = 0             # index for minus peaks in previous plus peak
        pair_centers = array(BYTE4,[])
        ip_max = len(pluspeaks)
        im_max = len(minuspeaks)
        flag_find_overlap = False
        while ip<ip_max and im<im_max:
            (pp,pn) = pluspeaks[ip] # for (peakposition, tagnumber in peak)
            (mp,mn) = minuspeaks[im]
            if pp-self.peaksize > mp: # move minus
                im += 1
            elif pp+self.peaksize < mp: # move plus
                ip += 1                 
                im = im_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    im_prev = im # only the first index is recorded
                if float(pn)/mn < 2 and float(pn)/mn > 0.5: # number tags in plus and minus peak region are comparable...
                    if pp < mp:
                        pair_centers.append((pp+mp)/2)
                        #self.debug ( "distance: %d, minus: %d, plus: %d" % (mp-pp,mp,pp))
                im += 1
        return pair_centers
            
    cdef __naive_find_peaks (self, np.ndarray[np.int32_t, ndim=1] taglist, int plus_strand=1 ):
        """Naively call peaks based on tags counting. 

        if plus_strand == 0, call peak on minus strand.

        Return peak positions and the tag number in peak region by a tuple list [(pos,num)].
        """
        cdef:
            long i
            int pos
            list peak_info
            int32_t * taglist_ptr
            list current_tag_list
        
        taglist_ptr = <int32_t *> taglist.data

        peak_info = []    # store peak pos in every peak region and
                          # unique tag number in every peak region
        if taglist.shape[0] < 2:
            return peak_info
        pos = taglist[0]

        current_tag_list = [ pos ]

        for i in range( 1, taglist.shape[0] ):
            pos = taglist_ptr[0]
            taglist_ptr += 1

            if ( pos - current_tag_list[0] + 1 ) > self.peaksize: # call peak in current_tag_list when the region is long enough
                # a peak will be called if tag number is ge min tags.
                if len(current_tag_list) >= self.min_tags and len(current_tag_list) <= self.max_tags:
                    peak_info.append( ( self.__naive_peak_pos(current_tag_list,plus_strand), len(current_tag_list) ) )
                #current_tag_list = array(BYTE4, []) # reset current_tag_list
                current_tag_list = []

            current_tag_list.append( pos )   # add pos while 1. no
                                             # need to call peak;
                                             # 2. current_tag_list is []
            
        return peak_info

    cdef __naive_peak_pos (self, pos_list, int plus_strand ):
        """Naively calculate the position of peak.

        plus_strand: 1, plus; 0, minus

        return the highest peak summit position.
        """
        #if plus_strand:
        #    tpos = pos_list + self.tag_expansion_size/2
        #else:
        #    tpos = pos_list - self.tag_expansion_size/2
        cdef:
            int peak_length, start, pos, i, pp, top_p_num, s, e, pileup, j
            np.ndarray[np.int32_t, ndim=1] horizon_line
            list ss, es
            int l_ss, l_es, i_s, i_e

        peak_length = pos_list[-1]+1-pos_list[0]+self.tag_expansion_size
        #if plus_strand:
        #    start = pos_list[0]
        #else:
        #    start = pos_list[0] - self.tag_expansion_size

        start = pos_list[0] - self.tag_expansion_size/2 # leftmost position of project line
        ss = []
        es = []

        #horizon_line = np.zeros(peak_length, dtype="int32") # the line for tags to be projected
        
        horizon_line = np.zeros( peak_length, dtype="int32") # the line for tags to be projected
        #horizon_line = array('i',[0]*peak_length)
        #for pos in pos_list:
        for i in range(len(pos_list)):
            pos = pos_list[i]
            #if plus_strand:
            ss.append( max(pos-start-self.tag_expansion_size/2,0) )
            es.append( min(pos-start+self.tag_expansion_size/2,peak_length) )

        ss.sort()
        es.sort()
        
        pileup = 0

        ls = len( ss )
        le = len( es )

        i_s = 0
        i_e = 0
        
        pre_p = min( ss[ 0 ], es[ 0 ] )
        if pre_p != 0:
            for i in range( pre_p ):
                horizon_line[ i ] = 0

        while i_s < ls and i_e < le:
            if ss[ i_s ] < es[ i_e ]:
                p = ss[ i_s ]
                if p != pre_p:
                    for i in range( pre_p, p ):
                        horizon_line[ i ] = pileup
                    pre_p = p
                pileup += 1
                i_s += 1
            elif ss[ i_s ] > es[ i_e ]:
                p = es[ i_e ]
                if p != pre_p:
                    for i in range( pre_p, p):
                        horizon_line[ i ] = pileup
                    pre_p = p
                pileup -= 1
                i_e += 1
            else:
                i_s += 1
                i_e += 1
        if ( i_e < ls ):
            for j in range( i_e, ls ):
                p = es[ i_e ]
                if p != pre_p:
                    for i in range( pre_p, p ):
                        horizon_line[ i ] = pileup
                    pre_p = p
                pileup -= 1

        # # top indices
        # top_indices = np.where(horizon_line == horizon_line.max())[0]
        # #print top_indices+start
        # print top_indices[ int(top_indices.shape[0]/2) ] + start
        # return top_indices[ int(top_indices.shape[0]/2) ] + start


        top_pos = []            # to record the top positions. Maybe > 1
        top_p_num = 0           # the maximum number of projected points
        for pp in range(peak_length): # find the peak posistion as the highest point
           if horizon_line[pp] > top_p_num:
               top_p_num = horizon_line[pp]
               top_pos = [pp]
           elif horizon_line[pp] == top_p_num:
               top_pos.append(pp)
        
        #print top_pos[int(len(top_pos)/2)]+start
        return (top_pos[int(len(top_pos)/2)]+start)

    cdef __naive_peak_pos2 (self, pos_list, int plus_strand ):
        """Naively calculate the position of peak.

        plus_strand: 1, plus; 0, minus

        return the highest peak summit position.
        """
        cdef int peak_length, start, pos, i, pp, top_p_num
        
        peak_length = pos_list[-1]+1-pos_list[0]+self.tag_expansion_size
        if plus_strand:
            start = pos_list[0] 
        else:
            start = pos_list[0] - self.tag_expansion_size
        horizon_line = np.zeros(peak_length, dtype="int32") # the line for tags to be projected
        for i in range(len(pos_list)):
            pos = pos_list[i]
            if plus_strand:
                for pp in range(int(pos-start),int(pos-start+self.tag_expansion_size)): # projected point
                    horizon_line[pp] += 1
            else:
                for pp in range(int(pos-start-self.tag_expansion_size),int(pos-start)): # projected point
                    horizon_line[pp] += 1

        # top indices
        #print pos_list
        #print horizon_line
        top_indices = np.where(horizon_line == horizon_line.max())[0]
        #print top_indices+start
        return top_indices[ int(top_indices.shape[0]/2) ] + start

# smooth function from SciPy cookbook: http://www.scipy.org/Cookbook/SignalSmooth
cpdef smooth(x, int window_len=11, str window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.
    
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
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
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
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]
