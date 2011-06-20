# Time-stamp: <2011-06-19 17:15:00 Tao Liu>

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
import os
from array import array
from copy import deepcopy
import gc                               # use garbage collectior

from MACS2.IO.cPeakIO import PeakIO

from MACS2.Constants import *
from MACS2.cPileup import pileup_bdg
from libc.math cimport log10

class PeakDetect:
    """Class to do the peak calling.

    e.g:
    >>> from MACS2.cPeakDetect import cPeakDetect
    >>> pd = PeakDetect(treat=treatdata, control=controldata, pvalue=pvalue_cutoff, d=100, scan_window=200, gsize=3000000000)
    >>> pd.call_peaks()
    """
    def __init__ (self,opt=None,treat=None, control=None):
        """Initialize the PeakDetect object.

        """
        self.opt  = opt
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn

        self.treat = treat
        self.control = control
        self.ratio_treat2control = None
        self.peaks = None
        self.final_peaks = None

        #self.femax = opt.femax
        #self.femin = opt.femin
        #self.festep = opt.festep
                
        self.pvalue = opt.log_pvalue    # -log10pvalue
        self.d = opt.d
        self.shift_size = self.d/2
        self.scan_window = opt.scanwindow
        self.gsize = opt.gsize
        
        self.nolambda = opt.nolambda

        self.sregion = opt.smalllocal
        self.lregion = opt.largelocal

        if (self.nolambda):
            self.info("#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!")
        #self.diag = opt.diag
        #self.save_score = opt.store_score
        self.zwig_tr = opt.zwig_tr
        self.zwig_ctl= opt.zwig_ctl

    def call_peaks (self):
        """Call peaks function.

        Scan the whole genome for peaks. RESULTS WILL BE SAVED IN
        self.final_peaks and self.final_negative_peaks.
        """
        if self.control:                # w/ control
            self.peaks = self.__call_peaks_w_control ()
        else:                           # w/o control
            self.peaks = self.__call_peaks_wo_control ()
        return self.peaks

    # def diag_result (self):
    #     """Run the diagnosis process on sequencing saturation.
        
    #     """
    #     if not self.diag:
    #         return None
    #     if self.control:                # w/ control
    #         return self.__diag_w_control()
    #     else:                           # w/o control
    #         return self.__diag_wo_control()

    def toxls (self):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        text = ""
        if self.peaks:
            text += "\t".join(("chr","start", "end",  "length",  "abs_summit", "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)"))+"\n"
            #text += "\t".join(("chr","start", "end",  "length",  "abs_summit", "-log10(pvalue)","-log10(qvalue)"))+"\n"
        else:
            return None

        peaks = self.peaks.peaks
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in peaks[chrom]:
                #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                text += "%s\t%d\t%d\t%d" % (chrom,peak[0]+1,peak[1],peak[2])
                text += "\t%d" % (peak[3]+1) # summit position
                text += "\t%.2f" % (peak[5]) # pileup height at summit
                text += "\t%.2f" % (peak[6]) # -log10pvalue at summit
                text += "\t%.2f" % (peak[7]) # fold change at summit                
                text += "\t%.2f" % (peak[8]) # -log10qvalue at summit
                text+= "\n"
        return text

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
        treat_total   = self.treat.total
        control_total = self.control.total
        self.ratio_treat2control = float(treat_total)/control_total

        # Now pileup FWTrackII to form a bedGraphTrackI
        self.info("#3 pileup treatment data by extending tags towards 3' to %d length" % self.d)
        treat_btrack = pileup_bdg(self.treat,self.d,halfextension=self.opt.halfext)

        if self.opt.tocontrol:
            # if user want to scale everything to control data
            treat_btrack.apply_func(lambda x:float(x)/self.ratio_treat2control) 

        # control data needs multiple steps of calculation
        # I need to shift them by 500 bps, then 5000 bps
        assert self.d <= self.sregion, "slocal can't be smaller than d!"
        assert self.sregion <= self.lregion , "llocal can't be smaller than slocal!"

        # d-size local
        self.info("#3 calculate d local lambda for control data")        

        # Now pileup FWTrackII to form a bedGraphTrackI
        c_tmp_btrack = pileup_bdg(self.control,self.d,directional=self.opt.shiftcontrol,halfextension=self.opt.halfext)
        if not self.opt.tocontrol:
            # if user want to scale everything to ChIP data
            tmp_v = self.ratio_treat2control
        else:
            tmp_v = 1
        c_tmp_btrack.apply_func(lambda x:float(x)*tmp_v)
        control_btrack = c_tmp_btrack

        # slocal size local
        if self.sregion:
            self.info("#3 calculate small local lambda for control data")        
            # Now pileup FWTrackII to form a bedGraphTrackI
            c_tmp_btrack = pileup_bdg(self.control,self.sregion,directional=self.opt.shiftcontrol,halfextension=self.opt.halfext)
            if not self.opt.tocontrol:
                # if user want to scale everything to ChIP data
                tmp_v = float(self.d)/self.sregion*self.ratio_treat2control
            else:
                tmp_v = float(self.d)/self.sregion
            c_tmp_btrack.apply_func(lambda x:float(x)*tmp_v)
            control_btrack = control_btrack.overlie(c_tmp_btrack,func=max)

        # llocal size local
        if self.lregion and self.lregion > self.sregion:
            self.info("#3 calculate large local lambda for control data")        
            # Now pileup FWTrackII to form a bedGraphTrackI
            c_tmp_btrack = pileup_bdg(self.control,self.lregion,directional=self.opt.shiftcontrol,halfextension=self.opt.halfext)
            if not self.opt.tocontrol:
                # if user want to scale everything to ChIP data
                tmp_v = float(self.d)/self.lregion*self.ratio_treat2control
            else:
                tmp_v = float(self.d)/self.lregion            
            c_tmp_btrack.apply_func(lambda x:float(x)*tmp_v)
            control_btrack = control_btrack.overlie(c_tmp_btrack,func=max)

        control_btrack.reset_baseline(float(self.d)*treat_total/self.gsize) # set the baseline as lambda_bg

        # calculate pvalue scores
        self.info("#3 Build score track ...")
        score_btrack = treat_btrack.make_scoreTrack_for_macs(control_btrack)
        treat_btrack = None             # clean them
        control_btrack = None
        gc.collect()                    # full collect garbage

        self.info("#3 Calculate qvalues ...")
        pqtable = score_btrack.make_pq_table()
        
        self.info("#3 Saving p-value to q-value table ...")
        pqfhd = open(self.opt.pqtable,"w")
        pqfhd.write( "-log10pvalue\t-log10qvalue\n" )
        for p in sorted(pqtable.keys()):
            q = pqtable[p]
            pqfhd.write("%.2f\t%.2f\n" % (p/100.0,q/100.0))
        pqfhd.close()

        self.info("#3 Assign qvalues ...")
        score_btrack.assign_qvalue( pqtable )
                
        # call peaks
        self.info("#3 Call peaks with given -log10pvalue cutoff: %.2f ..." % self.pvalue)        
        peaks = score_btrack.call_peaks(cutoff=self.pvalue*100,min_length=self.d,max_gap=self.opt.tsize,colname='-100logp')

        if self.opt.store_bdg:
           self.info("#3 save tag pileup into bedGraph file...")
           bdgfhd = open(self.zwig_tr + "_pileup.bdg", "w")
           score_btrack.write_bedGraph( bdgfhd,
                                        self.zwig_tr,
                                        "Fragment pileup at each bp from MACS version %s" % MACS_VERSION,
                                        "sample" )
        
           self.info("#3 save local lambda into bedGraph file...")
           bdgfhd = open(self.zwig_ctl + "_lambda.bdg", "w")
           score_btrack.write_bedGraph( bdgfhd,
                                        self.zwig_ctl,
                                        "Maximum local lambda at each bp from MACS version %s" % MACS_VERSION,
                                        "control" )

           self.info("#3 save the -log10pvalue score track into bedGraph file...")
           bdgfhd = open(self.zwig_tr + "_pvalue.bdg", "w")
           score_btrack.write_bedGraph( bdgfhd,
                                        self.zwig_tr+"_-log10pvalue",
                                        "-log10 pvalue scores at each bp from MACS version %s" % MACS_VERSION,
                                        "-100logp")
            
           self.info("#3 save the -log10qvalue score track into bedGraph file...")
           bdgfhd = open(self.zwig_tr + "_qvalue.bdg", "w")
           score_btrack.write_bedGraph( bdgfhd,
                                        self.zwig_tr+"_-log10qvalue",
                                        "-log10 qvalue scores at each bp from MACS version %s" % MACS_VERSION,
                                        "-100logq")
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
        # llocal size local
        if self.lregion:
            self.info("#3 calculate large local lambda from treatment data")
            # Now pileup FWTrackII to form a bedGraphTrackI
            c_tmp_btrack = pileup_bdg(self.treat,self.lregion,directional=self.opt.shiftcontrol,halfextension=self.opt.halfext)
            tmp_v = float(self.d)/self.lregion
            c_tmp_btrack.apply_func(lambda x:float(x)*tmp_v)
            control_btrack = control_btrack.overlie(c_tmp_btrack,func=max)
 
        # Now pileup FWTrackII to form a bedGraphTrackI
        self.info("#3 pileup treatment data")
        treat_btrack = pileup_bdg(self.treat,self.d,halfextension=self.opt.halfext)

        if self.opt.store_bdg:
            self.info("#3 save tag pileup into bedGraph file...")
            bdgfhd = open(self.zwig_tr + "_pileup.bdg", "w")
            treat_btrack.write_bedGraph(bdgfhd,name=self.zwig_tr,description="Fragment pileup at each bp from MACS version %s" % MACS_VERSION)

        # calculate pvalue scores
        self.info("#3 Calculate pvalue scores...")
        control_btrack.reset_baseline(float(self.d)*self.treat.total/self.gsize) # set the baseline as lambda_bg

        l_bd = float(self.d)*self.treat.total/self.gsize
        treat_btrack.apply_func(lambda x:-10*poisson_cdf(x,l_bd,lower=False,log10=True))
        self.score_btrack = treat_btrack
        if self.opt.store_bdg:
            self.info("#3 save the score track into bedGraph file...")
            bdgfhd = open(self.zwig_tr + "_scores.bdg", "w")
            self.score_btrack.write_bedGraph(bdgfhd,name=self.zwig_tr+"_Scores",description="-10log10 pvalue scores at each bp from MACS version %s" % MACS_VERSION)

        # call peaks
        self.info("#3 Call peaks with given score cutoff: %.2f ..." % self.pvalue)        
        peaks = self.score_btrack.call_peaks(cutoff=self.pvalue,min_length=self.d,max_gap=self.opt.tsize)
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

