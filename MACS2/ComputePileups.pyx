"""Module Description

Copyright (c) 2008,2009 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2013, Tao Liu lab at UB and Xiaole Shirley Liu lab at DFCI
All rights reserved.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

"""Added for WACS by Aseel Awdeh <araed104@uottawa.ca>

ComputePileups allows the user to compute pileups per sample individually
and returns the pileups as position value arrays.

"""

from itertools import groupby
from operator import itemgetter
import io
import gc                               # use garbage collectior

from MACS2.IO.PeakIO import PeakIO
from MACS2.IO.BedGraphIO import bedGraphIO
from MACS2.Constants import *
from MACS2.IO.PileupUnit import PileupUnitFromAlignments

cdef str subpeak_letters(short i):
    if i < 26:
        return chr(97+i)
    else:
        return subpeak_letters(i / 26) + chr(97 + (i % 26))

class ComputePileup:
    """Class to do the pileup computations.

    e.g
    >>> from MACS2.ComputePileups import ComputePileup
    >>> pd = ComputePileup(treat=treatdata, control=controldata, pvalue=pvalue_cutoff, d=100, gsize=3000000000)
    >>> treatment_pileup = pd.pileup_treatment()
    """
    def __init__ (self, opt = None,treat = None, control = None, d = None,
                  slocal = None, llocal = None):
        """Initialize the ComputePileup object.

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

        #Added for the weighted pileups
        self.controlweight = opt.controlweight #weight for one control
        self.pileup_treat_chrom = None
        self.pileup_control_chrom = None
        #----------------------------------

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


    def pileup_treatment (self):
        """
            To find pileups of treatment.
        """
        cdef:
            int i
            float lambda_bg, effective_depth_in_million
            float treat_scale, d
            list ctrl_scale_s, ctrl_d_s
            long treat_total, control_total
            long treat_sum              # approx sum of treatment pileup values
            long control_sum            # approx sum of control pileup values


        print("Compute pileup for treatment only.")
        treat_total   = self.treat.total

        #Aseel
        #For each control in the control list, calculate the control_total and the ratio_treat2control
        if self.PE_MODE:
            d = self.treat.average_template_length
            treat_sum = self.treat.length
            control_total = self.control.total * 2 # in PE mode, entire fragment is counted as 1, in treatment whereas both ends of fragment are counted in control/input.
            control_sum = control_total * self.treat.average_template_length
            self.ratio_treat2control = float(treat_sum)/control_sum
        else:
            d = self.d
            treat_sum = self.treat.total * self.d
            control_total = self.control.total
            control_sum = self.control.total * self.d
            self.ratio_treat2control = float(treat_sum)/control_sum

        if self.opt.ratio != 1.0:
            self.ratio_treat2control = self.opt.ratio

	    # Scaling of controls will always be towards treatment, as there are many controls and one treatment so makes more sense. -- removed if statement
        # if MACS decides to scale control to treatment because control sample is bigger
        effective_depth_in_million = treat_total / 1000000.0
        lambda_bg = float( treat_sum )/ self.gsize
        treat_scale = 1.0

        #Add array of arrays for scaling controls, also for weights
        scorecalculator = PileupUnitFromAlignments( self.treat,
                                                    self.control, #Can be a list.. but one control at a time.
                                                    d = d,
                                                    treat_scaling_factor = treat_scale,
                                                    end_shift = self.end_shift,
                                                    lambda_bg = lambda_bg,
                                                  )

        if self.opt.trackline: scorecalculator.enable_trackline()

        #Compute pileup
        pileup_treat_chrom = scorecalculator.compute_pileup__treatment()
        return pileup_treat_chrom

    def pileup_control (self):
        """
            To find pileups of multiple weighted control data -- one control dataset at a time.
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

        #For each control in the control list, calculate the control_total and the ratio_treat2control
        if self.PE_MODE:
            d = self.treat.average_template_length
            treat_sum = self.treat.length
            control_total = self.control.total * 2 # in PE mode, entire fragment is counted as 1, in treatment whereas both ends of fragment are counted in control/input.
            control_sum = control_total * self.treat.average_template_length
            self.ratio_treat2control = float(treat_sum)/control_sum
        else:
            d = self.d
            treat_sum = self.treat.total * self.d
            control_total = self.control.total
            control_sum = self.control.total * self.d
            self.ratio_treat2control = float(treat_sum)/control_sum

        if self.opt.ratio != 1.0:
            self.ratio_treat2control = self.opt.ratio

	    # Scaling of controls will always be towards treatment, as there are many controls and one treatment so makes more sense. -- removed if statement
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


        #For each control in the control list, prepare a list of extension sizes and
        # a list of scaling factors for control

        # Now prepare a list of extension sizes -- same for all controls
        ctrl_d_s = [ self.d ]   # note, d doesn't make sense in PE mode.

        # slocal size local
        if self.sregion:
            ctrl_d_s.append( self.sregion )

        if self.lregion and self.lregion > self.sregion:
            ctrl_d_s.append( self.lregion )

        # And a list of scaling factors for each control
        ctrl_scale_s = []

        # d
        # if user wants to scale everything to ChIP data
        tmp_v = self.ratio_treat2control
        ctrl_scale_s.append( tmp_v )

        # rlocal size sregion
        # if user want to scale everything to ChIP data
        tmp_v = float(self.d)/self.sregion*self.ratio_treat2control
        ctrl_scale_s.append( tmp_v )

        # llocal size local
        # if user want to scale everything to ChIP data
        tmp_v = float(self.d)/self.lregion*self.ratio_treat2control
        ctrl_scale_s.append( tmp_v )

        #Add array of arrays for scaling controls, also for weights --
        scorecalculator = PileupUnitFromAlignments( self.treat,
                                                    self.control, #Can be a list.. but one control at a time.
                                                    d = d, ctrl_d_s = ctrl_d_s,
                                                    treat_scaling_factor = treat_scale,
                                                    ctrl_scaling_factor_s = ctrl_scale_s, #This changed to take a list of lists (list of controls each having different scaling factors)
                                                    ctrl_weight =  self.controlweight, #This is also added. A weight is associated with each control. In this case, only one weight passed b/c one control.
                                                    end_shift = self.end_shift,
                                                    lambda_bg = lambda_bg,
                                                  )

        if self.opt.trackline: scorecalculator.enable_trackline()

        #Compute pileup
        pileup_control_chrom = scorecalculator.compute_pileup__control()
        return pileup_control_chrom
