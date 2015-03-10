# Time-stamp: <2015-03-05 15:37:59 Tao Liu>

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
#import os
#from array import array
#from copy import deepcopy
from itertools import groupby
from operator import itemgetter
import io
#import subprocess
import gc                               # use garbage collectior

from MACS2.IO.PeakIO import PeakIO
from MACS2.IO.BedGraphIO import bedGraphIO
from MACS2.Constants import *
from MACS2.IO.CallPeakUnit import CallerFromAlignments

cdef str subpeak_letters(short i):
    if i < 26:
        return chr(97+i)
    else:
        return subpeak_letters(i / 26) + chr(97 + (i % 26))

class PeakDetect:
    """Class to do the peak calling.

    e.g
    >>> from MACS2.cPeakDetect import cPeakDetect
    >>> pd = PeakDetect(treat=treatdata, control=controldata, pvalue=pvalue_cutoff, d=100, gsize=3000000000)
    >>> pd.call_peaks()
    """
    def __init__ (self,opt = None,treat = None, control = None, d = None,
                  slocal = None, llocal = None):
        """Initialize the PeakDetect object.

        """
        self.opt = opt
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn

        self.treat = treat
        self.control = control
        self.ratio_treat2control = None
        self.peaks = None
        self.final_peaks = None
        self.PE_MODE = opt.PE_MODE
        self.scoretrack = None

        #self.femax = opt.femax
        #self.femin = opt.femin
        #self.festep = opt.festep
                
        self.log_pvalue = opt.log_pvalue    # -log10pvalue
        self.log_qvalue = opt.log_qvalue    # -log10qvalue
        if d != None:
            self.d = d
        else:
            self.d = self.opt.d
        self.end_shift = self.opt.shift
        self.gsize = opt.gsize
        
        self.nolambda = opt.nolambda

        if slocal != None:
            self.sregion = slocal
        else:
            self.sregion = opt.smalllocal

        if llocal != None:
            self.lregion = llocal
        else:
            self.lregion = opt.largelocal

        if (self.nolambda):
            self.info("#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!")
        #self.diag = opt.diag
        #self.save_score = opt.store_score
        #self.zwig_tr = opt.zwig_tr
        #self.zwig_ctl= opt.zwig_ctl

    def call_peaks (self):
        """Call peaks function.

        Scan the whole genome for peaks. RESULTS WILL BE SAVED IN
        self.final_peaks and self.final_negative_peaks.
        """
        if self.control:                # w/ control
            #if self.opt.broad:
            #    (self.peaks,self.broadpeaks) = self.__call_peaks_w_control()
            #else:
            self.peaks = self.__call_peaks_w_control ()
        else:                           # w/o control
            #if self.opt.broad:
            #    (self.peaks,self.broadpeaks) = self.__call_peaks_wo_control()
            #else:
            self.peaks = self.__call_peaks_wo_control ()
        return self.peaks

    def __call_peaks_w_control (self):
        """To call peaks with control data.

        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment) <-- tuple type

        While calculating pvalue:

        First, t and c will be adjusted by the ratio between total
        reads in treatment and total reads in control, depending on
        --to-small option.

        Then, t and c will be multiplied by the smallest peak size --
        self.d.

        Finally, a poisson CDF is applied to calculate one-side pvalue
        for enrichment.
        """
        cdef:
            int i
            float lambda_bg, effective_depth_in_million
            float treat_scale, d
            list ctrl_scale_s, ctrl_d_s
            long treat_total, control_total
            long treat_sum              # approx sum of treatment pileup values
            long control_sum            # approx sum of control pileup values
            
        treat_total   = self.treat.total
        
        if self.PE_MODE:
            d = self.treat.average_template_length
            control_total = self.control.total * 2 # in PE mode, entire fragment is counted as 1 
                                                   # in treatment whereas both ends of fragment are counted in control/input.
            treat_sum = self.treat.length
            control_sum = control_total * self.treat.average_template_length
            self.ratio_treat2control = float(treat_sum)/control_sum
        else:
            d = self.d
            control_total = self.control.total
            treat_sum = self.treat.total * self.d
            control_sum = self.control.total * self.d
            self.ratio_treat2control = float(treat_sum)/control_sum            

        if self.opt.ratio != 1.0:
            self.ratio_treat2control = self.opt.ratio

        if self.opt.tocontrol:
            # if MACS decides to scale treatment to control data because treatment is bigger
            effective_depth_in_million = control_total / 1000000.0
            lambda_bg = float( control_sum )/ self.gsize
            treat_scale = 1/self.ratio_treat2control
        else:
            # if MACS decides to scale control to treatment because control sample is bigger
            effective_depth_in_million = treat_total / 1000000.0
            lambda_bg = float( treat_sum )/ self.gsize
            treat_scale = 1.0

        # prepare d_s for control data
        if self.sregion:
            assert self.d <= self.sregion, "slocal can't be smaller than d!"
        if self.lregion:
            assert self.d <= self.lregion , "llocal can't be smaller than d!"            
            assert self.sregion <= self.lregion , "llocal can't be smaller than slocal!"

        # Now prepare a list of extension sizes
        ctrl_d_s = [ self.d ]   # note, d doesn't make sense in PE mode.
        # And a list of scaling factors for control
        ctrl_scale_s = []

        # d
        if not self.opt.tocontrol:
            # if user wants to scale everything to ChIP data
            tmp_v = self.ratio_treat2control
        else:
            tmp_v = 1.0
        ctrl_scale_s.append( tmp_v )

        # slocal size local
        if self.sregion:
            ctrl_d_s.append( self.sregion )
            if not self.opt.tocontrol:
                # if user want to scale everything to ChIP data
                tmp_v = float(self.d)/self.sregion*self.ratio_treat2control
            else:
                tmp_v = float(self.d)/self.sregion
            ctrl_scale_s.append( tmp_v )

        # llocal size local
        if self.lregion and self.lregion > self.sregion:
            ctrl_d_s.append( self.lregion )
            if not self.opt.tocontrol:
                # if user want to scale everything to ChIP data
                tmp_v = float(self.d)/self.lregion*self.ratio_treat2control
            else:
                tmp_v = float(self.d)/self.lregion
            ctrl_scale_s.append( tmp_v )                            
            
        #if self.PE_MODE:        # first d/scale are useless in PE mode
        #    ctrl_d_s = ctrl_d_s[1:]
        #    ctrl_scale_s = ctrl_scale_s[1:]
        #    print ctrl_d_s
        #    print ctrl_scale_s
        if self.nolambda:
            ctrl_d_s = []
            ctrl_scale_s = []

        scorecalculator = CallerFromAlignments( self.treat, self.control,
                                                d = d, ctrl_d_s = ctrl_d_s, 
                                                treat_scaling_factor = treat_scale, 
                                                ctrl_scaling_factor_s = ctrl_scale_s,
                                                end_shift = self.end_shift,
                                                lambda_bg = lambda_bg,
                                                save_bedGraph = self.opt.store_bdg,
                                                bedGraph_filename_prefix = self.opt.name,
                                                bedGraph_treat_filename = self.opt.bdg_treat,
                                                bedGraph_control_filename = self.opt.bdg_control,
                                                save_SPMR = self.opt.do_SPMR,
                                                cutoff_analysis_filename = self.opt.cutoff_analysis_file )

        if self.opt.trackline: scorecalculator.enable_trackline()

        # call peaks
        call_summits = self.opt.call_summits
        if call_summits: self.info("#3 Going to call summits inside each peak ...")

        if self.log_pvalue:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10pvalue cutoff and level2: %.5f, %.5f..." % (self.log_pvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['p',], lvl1_cutoff_s=[self.log_pvalue,],lvl2_cutoff_s=[self.opt.log_broadcutoff,],min_length=self.d,
                                                        lvl1_max_gap=self.opt.tsize,lvl2_max_gap=self.d*4,
                                                        auto_cutoff=self.opt.cutoff_analysis )
            else:
                self.info("#3 Call peaks with given -log10pvalue cutoff: %.5f ..." % self.log_pvalue)
                peaks = scorecalculator.call_peaks( ['p',], [self.log_pvalue,],
                                                    min_length=self.d,
                                                    max_gap=self.opt.tsize,
                                                    call_summits=call_summits,
                                                    auto_cutoff=self.opt.cutoff_analysis )
        elif self.log_qvalue:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10qvalue cutoff and level2: %f, %f..." % (self.log_qvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['q',], lvl1_cutoff_s=[self.log_qvalue,],lvl2_cutoff_s=[self.opt.log_broadcutoff,],min_length=self.d,
                                                        lvl1_max_gap=self.opt.tsize,lvl2_max_gap=self.d*4,
                                                        auto_cutoff=self.opt.cutoff_analysis )
            else:
                peaks = scorecalculator.call_peaks( ['q',], [self.log_qvalue,],
                                                    min_length=self.d,
                                                    max_gap=self.opt.tsize,
                                                    call_summits=call_summits,
                                                    auto_cutoff=self.opt.cutoff_analysis )
        scorecalculator.destroy()
        return peaks

    def __call_peaks_wo_control (self):
        """To call peaks without control data.

        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment) <-- tuple type

        While calculating pvalue:

        First, t and c will be adjusted by the ratio between total
        reads in treatment and total reads in control, depending on
        --to-small option.

        Then, t and c will be multiplied by the smallest peak size --
        self.d.

        Finally, a poisson CDF is applied to calculate one-side pvalue
        for enrichment.
        """
        cdef float lambda_bg, effective_depth_in_million
        cdef float treat_scale = 1
        cdef float d
        cdef list ctrl_scale_s, ctrl_d_s

        if self.PE_MODE: d = 0
        else: d = self.d
        treat_length = self.treat.length
        treat_total = self.treat.total
        
        effective_depth_in_million = treat_total / 1000000.0

        # global lambda
        if self.PE_MODE:
        #    # this an estimator, we should maybe test it for accuracy?
            lambda_bg = treat_length / self.gsize
        else:
            lambda_bg = float(d) * treat_total / self.gsize
        treat_scale = 1.0

        # slocal and d-size local bias are not calculated!
        # nothing done here. should this match w control??
        
        if not self.nolambda:
            if self.PE_MODE:
                ctrl_scale_s = [ float(treat_length) / (self.lregion*treat_total*2), ]
            else:
                ctrl_scale_s = [ float(self.d) / self.lregion, ]
            ctrl_d_s     = [ self.lregion, ]
        else:
            ctrl_scale_s = []
            ctrl_d_s     = []

        scorecalculator = CallerFromAlignments( self.treat, None, 
                                                d = d, ctrl_d_s = ctrl_d_s, 
                                                treat_scaling_factor = treat_scale, 
                                                ctrl_scaling_factor_s = ctrl_scale_s,
                                                end_shift = self.end_shift,
                                                lambda_bg = lambda_bg,
                                                save_bedGraph = self.opt.store_bdg,
                                                bedGraph_filename_prefix = self.opt.name,
                                                bedGraph_treat_filename = self.opt.bdg_treat,
                                                bedGraph_control_filename = self.opt.bdg_control,
                                                save_SPMR = self.opt.do_SPMR,
                                                cutoff_analysis_filename = self.opt.cutoff_analysis_file )

        if self.opt.trackline: scorecalculator.enable_trackline()

        # call peaks
        call_summits = self.opt.call_summits
        if call_summits: self.info("#3 Going to call summits inside each peak ...")

        if self.log_pvalue:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10pvalue cutoff and level2: %.5f, %.5f..." % (self.log_pvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['p',], lvl1_cutoff_s=[self.log_pvalue,],lvl2_cutoff_s=[self.opt.log_broadcutoff,],min_length=self.d,
                                                        lvl1_max_gap=self.opt.tsize,lvl2_max_gap=self.d*4)
            else:
                self.info("#3 Call peaks with given -log10pvalue cutoff: %.5f ..." % self.log_pvalue)
                peaks = scorecalculator.call_peaks( ['p',], [self.log_pvalue,],
                                                    min_length=self.d,
                                                    max_gap=self.opt.tsize,
                                                    call_summits=call_summits,
                                                    auto_cutoff=self.opt.cutoff_analysis )
        elif self.log_qvalue:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10qvalue cutoff and level2: %f, %f..." % (self.log_qvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['q',], lvl1_cutoff_s=[self.log_qvalue,],lvl2_cutoff_s=[self.opt.log_broadcutoff,],min_length=self.d,
                                                        lvl1_max_gap=self.opt.tsize,lvl2_max_gap=self.d*4)
            else:
                peaks = scorecalculator.call_peaks( ['q',], [self.log_qvalue,],
                                                    min_length=self.d,
                                                    max_gap=self.opt.tsize,
                                                    call_summits=call_summits,
                                                    auto_cutoff=self.opt.cutoff_analysis )
        scorecalculator.destroy()
        return peaks

    # def __diag_w_control (self):
    #     # sample
    #     sample_peaks = {}
    #     for i in xrange(90,10,-10):
    #         self.info("#3 diag: sample %d%%" % i)
    #         sample_peaks[i]=self.__diag_peakfinding_w_control_sample(float(i)/(i+10))
    #     return self.__overlap (self.final_peaks, sample_peaks,top=90,bottom=10,step=-10)

    # def __diag_peakfinding_w_control_sample (self, percent):
    #     self.treat.sample(percent) # because sampling is after
    #                                # shifting, track.total is used
    #                                # now.
    #     self.control.sample(percent)
    #     ratio_treat2control = float(self.treat.total)/self.control.total

    #     #self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize # bug fixed...
    #     #self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1

    #     self.debug("#3 diag: after shift and merging, treat: %d, control: %d" % (self.treat.total,self.control.total))
    #     self.info("#3 diag: call peak candidates")
    #     peak_candidates = self.__call_peaks_from_trackI (self.treat)

    #     self.info("#3 diag: call negative peak candidates")
    #     negative_peak_candidates = self.__call_peaks_from_trackI (self.control)
        
    #     self.info("#3 diag: use control data to filter peak candidates...")
    #     final_peaks_percent = self.__filter_w_control(peak_candidates,self.treat,self.control, ratio_treat2control)
    #     return final_peaks_percent
        
    # def __diag_wo_control (self):
    #     # sample
    #     sample_peaks = {}
    #     for i in xrange(90,10,-10):
    #         self.info("#3 diag: sample %d%%" % i)
    #         sample_peaks[i]=self.__diag_peakfinding_wo_control_sample(float(i)/(i+10))
    #     return self.__overlap (self.final_peaks, sample_peaks,top=90,bottom=10,step=-10)

    # def __diag_peakfinding_wo_control_sample (self, percent):

    #     #self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize # bug fixed...
    #     #self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1

    #     self.treat.sample(percent)
    #     self.debug("#3 diag: after shift and merging, tags: %d" % (self.treat.total))
    #     self.info("#3 diag: call peak candidates")
    #     peak_candidates = self.__call_peaks_from_trackI (self.treat)
    #     self.info("#3 diag: use self to calculate local lambda and  filter peak candidates...")
    #     final_peaks_percent = self.__filter_w_control(peak_candidates,self.treat,self.treat,pass_sregion=True) # bug fixed...
    #     return final_peaks_percent

    # def __overlap (self, gold_peaks, sample_peaks, top=90,bottom=10,step=-10):
    #     """Calculate the overlap between several fe range for the
    #     golden peaks set and results from sampled data.
        
    #     """
    #     gp = PeakIO()
    #     gp.init_from_dict(gold_peaks)
    #     if self.femax:
    #         femax = min(self.femax, (int(gp.max_fold_enrichment())//self.festep+1)*self.festep)
    #     else:
    #         femax = (int(gp.max_fold_enrichment())//self.festep+1)*self.festep
    #     femin = self.femin
    #     diag_result = []
    #     for f in xrange(femin, femax, self.festep):
            
    #         fe_low = f
    #         fe_up = f + self.festep
    #         self.debug("#3 diag: fe range = %d -- %d" % (fe_low, fe_up))
            
    #         r = self.__overlap_fe(gold_peaks, sample_peaks, fe_low, fe_up, top, bottom, step)
    #         if r:
    #             diag_result.append(r)
    #     return diag_result

    # def __overlap_fe (self, gold_peaks, sample_peaks, fe_low, fe_up, top, bottom, step):
    #     ret = ["%d-%d" % (fe_low,fe_up)]
    #     gp = PeakIO()
    #     gp.init_from_dict(gold_peaks)
    #     gp.filter_fc(fe_low,fe_up)
    #     gptotal =  gp.total()
    #     if gptotal <= 0:
    #         return None

    #     ret.append(gptotal)
    #     for i in xrange(top,bottom,step):
    #         p = PeakIO()
    #         p.init_from_dict(sample_peaks[i])
    #         percent = 100.0*gp.overlap_with_other_peaks(p)/gptotal
    #         ret.append(percent)
    #         del p
    #     return ret


    # def __remove_overlapping_peaks (self, peaks ):
    #     """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

    #     """
    #     new_peaks = {}
    #     chrs = peaks.keys()
    #     chrs.sort()
    #     for chrom in chrs:
    #         new_peaks[chrom]=[]
    #         n_append = new_peaks[chrom].append
    #         prev_peak = None
    #         peaks_chr = peaks[chrom]
    #         for i in xrange(len(peaks_chr)):
    #             if not prev_peak:
    #                 prev_peak = peaks_chr[i]
    #                 continue
    #             else:
    #                 if peaks_chr[i][0] <= prev_peak[1]:
    #                     s_new_peak = prev_peak[0]
    #                     e_new_peak = peaks_chr[i][1]
    #                     l_new_peak = e_new_peak-s_new_peak
    #                     if peaks_chr[i][4] > prev_peak[4]:
    #                         summit_new_peak = peaks_chr[i][3]
    #                         h_new_peak = peaks_chr[i][4]
    #                     else:
    #                         summit_new_peak = prev_peak[3]
    #                         h_new_peak = prev_peak[4]
    #                     prev_peak = (s_new_peak,e_new_peak,l_new_peak,summit_new_peak,h_new_peak,peaks_chr[i][5]+prev_peak[5])
    #                 else:
    #                     n_append(prev_peak)
    #                     prev_peak = peaks_chr[i]
    #         if prev_peak:
    #             n_append(prev_peak)
    #     return new_peaks

