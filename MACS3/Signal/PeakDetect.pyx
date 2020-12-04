# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-11-24 17:39:12 Tao Liu>

"""Module Description: Detect peaks, main module

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# Python modules
# ------------------------------------
from itertools import groupby
from operator import itemgetter
import io
import gc                               # use garbage collectior

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.IO.PeakIO import PeakIO
from MACS3.IO.BedGraphIO import bedGraphIO
from MACS3.Utilities.Constants import *
from MACS3.Signal.CallPeakUnit import CallerFromAlignments

cdef bytes subpeak_letters(short i):
    if i < 26:
        return chr(97+i).encode()
    else:
        return subpeak_letters(i // 26) + chr(97 + (i % 26)).encode()

class PeakDetect:
    """Class to do the peak calling.

    e.g
    >>> from MACS3.cPeakDetect import cPeakDetect
    >>> pd = PeakDetect(treat=treatdata, control=controldata, pvalue=pvalue_cutoff, d=100, gsize=3000000000)
    >>> pd.call_peaks()
    """
    def __init__ (self,opt = None,treat = None, control = None, d = None,
                  maxgap = None, minlen = None, slocal = None, llocal = None):
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

        if opt.maxgap:
            self.maxgap = opt.maxgap
        else:
            self.maxgap = opt.tsize

        if opt.minlen:
            self.minlen = opt.minlen
        else:
            self.minlen = self.d

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
            assert self.d <= self.sregion, f"{self.sregion:} can't be smaller than {self.d:}!"
        if self.lregion:
            assert self.d <= self.lregion , f"{self.lregion:} can't be smaller than {self.d:}!"
            assert self.sregion <= self.lregion , f"{self.lregion:} can't be smaller than {self.sregion:}!"

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

        if self.log_pvalue != None:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10pvalue cutoff and level2: %.5f, %.5f..." % (self.log_pvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['p',],
                                                    lvl1_cutoff_s=[self.log_pvalue,],
                                                    lvl2_cutoff_s=[self.opt.log_broadcutoff,],
                                                    min_length=self.minlen,
                                                    lvl1_max_gap=self.maxgap,
                                                    lvl2_max_gap=self.maxgap*4,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
            else:
                self.info("#3 Call peaks with given -log10pvalue cutoff: %.5f ..." % self.log_pvalue)
                peaks = scorecalculator.call_peaks( ['p',], [self.log_pvalue,],
                                                    min_length=self.minlen,
                                                    max_gap=self.maxgap,
                                                    call_summits=call_summits,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
        elif self.log_qvalue != None:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10qvalue cutoff and level2: %f, %f..." % (self.log_qvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['q',],
                                                    lvl1_cutoff_s=[self.log_qvalue,],
                                                    lvl2_cutoff_s=[self.opt.log_broadcutoff,],
                                                    min_length=self.minlen,
                                                    lvl1_max_gap=self.maxgap,
                                                    lvl2_max_gap=self.maxgap*4,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
            else:
                peaks = scorecalculator.call_peaks( ['q',], [self.log_qvalue,],
                                                    min_length=self.minlen,
                                                    max_gap=self.maxgap,
                                                    call_summits=call_summits,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
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

        if self.log_pvalue != None:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10pvalue cutoff and level2: %.5f, %.5f..." % (self.log_pvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['p',],
                                                    lvl1_cutoff_s=[self.log_pvalue,],
                                                    lvl2_cutoff_s=[self.opt.log_broadcutoff,],
                                                    min_length=self.minlen,
                                                    lvl1_max_gap=self.maxgap,
                                                    lvl2_max_gap=self.maxgap*4,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
            else:
                self.info("#3 Call peaks with given -log10pvalue cutoff: %.5f ..." % self.log_pvalue)
                peaks = scorecalculator.call_peaks( ['p',], [self.log_pvalue,],
                                                    min_length=self.minlen,
                                                    max_gap=self.maxgap,
                                                    call_summits=call_summits,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
        elif self.log_qvalue != None:
            if self.opt.broad:
                self.info("#3 Call broad peaks with given level1 -log10qvalue cutoff and level2: %f, %f..." % (self.log_qvalue,self.opt.log_broadcutoff) )
                peaks = scorecalculator.call_broadpeaks(['q',],
                                                    lvl1_cutoff_s=[self.log_qvalue,],
                                                    lvl2_cutoff_s=[self.opt.log_broadcutoff,],
                                                    min_length=self.minlen,
                                                    lvl1_max_gap=self.maxgap,
                                                    lvl2_max_gap=self.maxgap*4,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
            else:
                peaks = scorecalculator.call_peaks( ['q',], [self.log_qvalue,],
                                                    min_length=self.minlen,
                                                    max_gap=self.maxgap,
                                                    call_summits=call_summits,
                                                    cutoff_analysis=self.opt.cutoff_analysis )
        scorecalculator.destroy()
        return peaks

