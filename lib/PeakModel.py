# Time-stamp: <2011-02-14 15:52:04 Tao Liu>

"""Module Description

Copyright (c) 2008,2009 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
import sys, time, random

def median (nums):
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

class PeakModel:
    """Peak Model class.
    """
    def __init__ (self, opt=None, treatment=None, max_pairnum=500, gz = 0, umfold=30, lmfold=10, bw=200, ts = 25, bg=0):
        self.treatment = treatment
        if opt:
            self.gz = opt.gsize
            self.umfold = opt.umfold
            self.lmfold = opt.lmfold
            self.tsize = opt.tsize
            self.bw = opt.bw
            self.info  = opt.info
            self.debug = opt.debug
            self.warn  = opt.warn
            self.error = opt.warn
        else:
            self.gz = gz
            self.umfold = umfold
            self.lmfold = lmfold            
            self.tsize = ts
            self.bg = bg
            self.bw = bw
            self.info  = lambda x: sys.stderr.write(x+"\n")
            self.debug = lambda x: sys.stderr.write(x+"\n")
            self.warn  = lambda x: sys.stderr.write(x+"\n")
            self.error = lambda x: sys.stderr.write(x+"\n")
			
        self.max_pairnum = max_pairnum
        self.summary = ""
        self.plus_line = None
        self.minus_line = None
        self.shifted_line = None
        self.d = None
        self.scan_window = None
        self.min_tags = None
        self.peaksize = None
        self.build()
    
    def build (self):
        """Build the model.

        prepare self.d, self.scan_window, self.plus_line,
        self.minus_line and self.shifted_line to use.
        """
        self.peaksize = 2*self.bw
        self.min_tags = float(self.treatment.total) * self.lmfold * self.peaksize / self.gz /2 # mininum unique hits on single strand
        self.max_tags = float(self.treatment.total) * self.umfold * self.peaksize / self.gz /2 # maximum unique hits on single strand
        #print self.min_tags
        #print self.max_tags
        # use treatment data to build model
        paired_peakpos = self.__paired_peaks ()
        # select up to 1000 pairs of peaks to build model
        num_paired_peakpos = 0
        num_paired_peakpos_remained = self.max_pairnum
        num_paired_peakpos_picked = 0
        for c in paired_peakpos.keys():
            num_paired_peakpos +=len(paired_peakpos[c])
            if num_paired_peakpos_remained == 0:
                paired_peakpos.pop(c)
            else:
                paired_peakpos[c] = paired_peakpos[c][:num_paired_peakpos_remained]
                num_paired_peakpos_remained -=  len(paired_peakpos[c])
                num_paired_peakpos_picked += len(paired_peakpos[c])

        self.info("#2 number of paired peaks: %d" % (num_paired_peakpos))
        if num_paired_peakpos < 100:
            self.error("Too few paired peaks (%d) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, please use --nomodel and --shiftsize 100 instead." % (num_paired_peakpos))
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

    def __paired_peak_model (self, paired_peakpos):
        """Use paired peak positions and treatment tag positions to build the model.

        Modify self.(d, model_shift size and scan_window size. and extra, plus_line, minus_line and shifted_line for plotting).
        """
        window_size = 1+2*self.peaksize
        self.plus_line = [0]*window_size
        self.minus_line = [0]*window_size
        for chrom in paired_peakpos.keys():
            paired_peakpos_chrom = paired_peakpos[chrom]
            tags = self.treatment.get_locations_by_chr(chrom)
            tags_plus =  tags[0]
            tags_minus = tags[1]
            # every paired peak has plus line and minus line
            #  add plus_line
            self.plus_line = self.__model_add_line (paired_peakpos_chrom, tags_plus,self.plus_line)
            #  add minus_line
            self.minus_line = self.__model_add_line (paired_peakpos_chrom, tags_minus,self.minus_line)

        # find top 
        plus_tops = []
        minus_tops = []
        plus_max = max(self.plus_line)
        minus_max = max(self.minus_line)
        for i in range(window_size):
            if self.plus_line[i] == plus_max:
                plus_tops.append(i)
            if self.minus_line[i] == minus_max:
                minus_tops.append(i)
        self.d = minus_tops[len(minus_tops)/2] - plus_tops[len(plus_tops)/2] + 1
        shift_size = self.d/2
        # find the median point
        #plus_median = median(self.plus_line) 
        #minus_median = median(self.minus_line)       

        
        self.scan_window = max(self.d,self.tsize)*2
        # a shifted model
        self.shifted_line = [0]*window_size
        plus_shifted = [0]*shift_size
        plus_shifted.extend(self.plus_line[:-1*shift_size])
        minus_shifted = self.minus_line[shift_size:]
        minus_shifted.extend([0]*shift_size)
        #print "d:",self.d,"shift_size:",shift_size
        #print len(self.plus_line),len(self.minus_line),len(plus_shifted),len(minus_shifted),len(self.shifted_line)
        for i in range(window_size):
            self.shifted_line[i]=minus_shifted[i]+plus_shifted[i]
        return True

    def __model_add_line (self, pos1, pos2, line):
        """Project each pos in pos2 which is included in
        [pos1-self.peaksize,pos1+self.peaksize] to the line.

        """
        i1 = 0                  # index for pos1
        i2 = 0                  # index for pos2
        i2_prev = 0             # index for pos2 in previous pos1
                                # [pos1-self.peaksize,pos1+self.peaksize]
                                # region
        i1_max = len(pos1)
        i2_max = len(pos2)
        last_p2 = -1
        flag_find_overlap = False
         
        while i1<i1_max and i2<i2_max:
            p1 = pos1[i1]
            p2 = pos2[i2]
            if p1-self.peaksize > p2: # move pos2
                i2 += 1
            elif p1+self.peaksize < p2: # move pos1
                i1 += 1                 
                i2 = i2_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    i2_prev = i2 # only the first index is recorded
                # project
                for i in range(p2-p1+self.peaksize-self.tsize/2,p2-p1+self.peaksize+self.tsize/2):
                    if i>=0 and i<len(line):
                        line[i]+=1
                i2+=1
        return line
            
    def __paired_peaks (self):
        """Call paired peaks from fwtrackI object.

        Return paired peaks center positions.
        """
        chrs = self.treatment.get_chr_names()
        chrs.sort()
        paired_peaks_pos = {}
        for chrom in chrs:
            self.debug("Chromosome: %s" % (chrom))
            tags = self.treatment.get_locations_by_chr(chrom)
            plus_peaksinfo = self.__naive_find_peaks (tags[0])
            self.debug("Number of unique tags on + strand: %d" % (len(tags[0])))            
            self.debug("Number of peaks in + strand: %d" % (len(plus_peaksinfo)))
            minus_peaksinfo = self.__naive_find_peaks (tags[1])
            self.debug("Number of unique tags on - strand: %d" % (len(tags[1])))            
            self.debug("Number of peaks in - strand: %d" % (len(minus_peaksinfo)))
            if not plus_peaksinfo or not minus_peaksinfo:
                self.debug("Chrom %s is discarded!" % (chrom))
                continue
            else:
                paired_peaks_pos[chrom] = self.__find_pair_center (plus_peaksinfo, minus_peaksinfo)
                self.debug("Number of paired peaks: %d" %(len(paired_peaks_pos[chrom])))
        return paired_peaks_pos

    def __find_pair_center (self, pluspeaks, minuspeaks):
        ip = 0                  # index for plus peaks
        im = 0                  # index for minus peaks
        im_prev = 0             # index for minus peaks in previous plus peak
        pair_centers = []
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
            
    def __naive_find_peaks (self, taglist ):
        """Naively call peaks based on tags counting. 

        Return peak positions and the tag number in peak region by a tuple list [(pos,num)].
        """
        peak_info = []    # store peak pos in every peak region and
                          # unique tag number in every peak region
        if len(taglist)<2:
            return peak_info
        pos = taglist[0]

        current_tag_list = [pos]   # list to find peak pos

        for i in range(1,len(taglist)):
            pos = taglist[i]

            if (pos-current_tag_list[0]+1) > self.peaksize: # call peak in current_tag_list
                # a peak will be called if tag number is ge min tags.
                if len(current_tag_list) >= self.min_tags and len(current_tag_list) <= self.max_tags:
                    peak_info.append((self.__naive_peak_pos(current_tag_list),len(current_tag_list)))
                current_tag_list = [] # reset current_tag_list

            current_tag_list.append(pos)   # add pos while 1. no
                                           # need to call peak;
                                           # 2. current_tag_list is []
        return peak_info

    def __naive_peak_pos (self, pos_list ):
        """Naively calculate the position of peak.

        return the highest peak summit position.
        """
        peak_length = pos_list[-1]+1-pos_list[0]+self.tsize
        start = pos_list[0] -self.tsize/2
        horizon_line = [0]*peak_length # the line for tags to be projected
        for pos in pos_list:
            for pp in range(int(pos-start-self.tsize/2),int(pos-start+self.tsize/2)): # projected point
                horizon_line[pp] += 1

        top_pos = []            # to record the top positions. Maybe > 1
        top_p_num = 0           # the maximum number of projected points
        for pp in range(peak_length): # find the peak posistion as the highest point
            if horizon_line[pp] > top_p_num:
                top_p_num = horizon_line[pp]
                top_pos = [pp]
            elif horizon_line[pp] == top_p_num:
                top_pos.append(pp)
        return (top_pos[int(len(top_pos)/2)]+start)
