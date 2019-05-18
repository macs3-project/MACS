"""Module for Calculate Scores.

Copyright (c) 2013 Tao Liu <vladimir.liu@gmail.com>
Copyright (c) 2013, Tao Liu lab at UB and Xiaole Shirley Liu lab at DFCI
All rights reserved.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: vladimir.liu@gmail.com
"""

"""Added for WACS by Aseel Awdeh <araed104@uottawa.ca>

PileupUnits allows the user to compute pileups per sample individually
and returns the pileups as position value arrays.

"""


# ------------------------------------
# python modules
# ------------------------------------

import numpy as np
cimport numpy as np

from collections import Counter
from copy import copy

from operator import itemgetter
import cPickle
from tempfile import mkstemp
import os

from cpython cimport bool

from MACS2.Signal import maxima, enforce_valleys, enforce_peakyness

#WACS
from MACS2.Pileup import max_over_two_pv_array, sum_over_two_pv_array, multiple_pv_array_w_weight, print_pv

from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
ctypedef np.float32_t float32_t

from libc.math cimport exp,log,log10, M_LN10, log1p, erf, sqrt, floor, ceil

from MACS2.IO.PeakIO import PeakIO, BroadPeakIO, parse_peakname
from MACS2.IO.FixWidthTrack import FWTrack
from MACS2.IO.PairedEndTrack import PETrackI

import logging
from time import time as ttime
from libc.stdio cimport *

# cdef extern from "stdio.h":
#     ctypedef struct FILE
#     FILE *fopen   (const char *filename, const char  *opentype)
#     #FILE * fopen ( const char * filename, const char * mode )
#     int  fclose   (FILE *stream)
#     int fprintf  (FILE *stream, const char *template, ...)

# ------------------------------------
# constants
# ------------------------------------
#__version__ = "scoreCalculate $Revision$"
#__author__ = "Tao Liu <vladimir.liu@gmail.com>"
#__doc__ = "scoreTrackI classes"

# ------------------------------------
# Classes
# ------------------------------------
cdef class PileupUnitFromAlignments:
    """A unit to calculate pileups from alignments --
    FWTrack or PETrackI objects.

    It will compute for each chromosome separately in order to save
    memory usage.
    """
    cdef:
        object treat            # FWTrack or PETrackI object for ChIP
        object ctrl             # FWTrack or PETrackI object for Control

        int  d                           # extension size for ChIP
        list ctrl_d_s                    # extension sizes for Control. Can be multiple values
        float treat_scaling_factor       # scaling factor for ChIP

        #Aseel
        list ctrl_scaling_factor_s       # scaling factor for Control, corresponding to each extension size. -- Needs to be changed for multiple controls
        float ctrl_weight                # weights for multiple controls
        bool multiple_ctrl               # if input is multiple weighted controls, this is set to true. else it is set to false

        float lambda_bg                  # minimum local bias to fill missing values
        list chromosomes                 # name of common chromosomes in ChIP and Control data
        float pseudocount                # the pseudocount used to calcuate logLR, FE or logFE
        str bedGraph_filename_prefix     # prefix will be added to _pileup.bdg for treatment and _lambda.bdg for control

        #SHIFTCONTROL is obsolete
        int  end_shift                   # shift of cutting ends before extension
        bool trackline                   # whether trackline should be saved in bedGraph
        bool save_bedGraph               # whether to save pileup and local bias in bedGraph files
        bool save_SPMR                   # whether to save pileup normalized by sequencing depth in million reads
        bool no_lambda_flag              # whether ignore local bias, and to use global bias instead
        bool PE_mode                     # whether it's in PE mode, will be detected during initiation
        # temporary data buffer
        str chrom                        # name of current chromosome

        dict pileup_data_files           # Record the names of temporary files for storing pileup values of each chromosome


    def __init__ (self, treat, ctrl,
                  int d = 200, list ctrl_d_s = [200, 1000, 10000],
                  float treat_scaling_factor = 1.0,
                  list ctrl_scaling_factor_s = [1.0, 0.2, 0.02],
                  float ctrl_weight = 0.0, #Default case where the controls are not weighted or when its a treatment sample.
                  bool stderr_on = False,
                  float pseudocount = 1.0,
                  int end_shift = 0,
                  float lambda_bg = 0,
                  ):
        """Initialize.

        A calculator is unique to each comparison of treat and
        control. Treat_depth and ctrl_depth should not be changed
        during calculation.

        treat and ctrl are two FWTrack or PETrackI objects.

        treat_depth and ctrl_depth are effective depth in million:
                                    sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        d, sregion, lregion: d is the fragment size, sregion is the
                             small region size, lregion is the large
                             region size

        pseudocount: a pseudocount used to calculate logLR, FE or
                     logFE. Please note this value will not be changed
                     with normalization method. So if you really want
                     to set pseudocount 1 per million reads, set it
                     after you normalize treat and control by million
                     reads by `change_normalizetion_method(ord('M'))`.

        """
        cdef:
            set chr1, chr2
            int i
            char * tmp
            bytes tmp_bytes


        if isinstance(treat, FWTrack):
            self.PE_mode = False
        elif isinstance(treat, PETrackI):
            self.PE_mode = True
        else:
            raise Exception("Should be FWTrack or PETrackI object!")

        self.treat = treat
        if ctrl:
            self.ctrl = ctrl
        else:                   # while there is no control
            self.ctrl = treat
        self.trackline = False
        self.d = d              # note, self.d doesn't make sense in PE mode
        self.ctrl_d_s = ctrl_d_s# note, self.d doesn't make sense in PE mode
        self.treat_scaling_factor = treat_scaling_factor

        #Added for WACS:------------------------------------------------
        self.ctrl_scaling_factor_s = ctrl_scaling_factor_s
        self.ctrl_weight = ctrl_weight
        if not isinstance(self.ctrl, list):
            self.multiple_ctrl = False
        elif len(self.ctrl) > 1 and self.ctrl_weight != 0.0:
            self.multiple_ctrl = True
            logging.info("#3 multiple_ctrl set to True...")
        elif len(self.ctrl) == 1 and self.ctrl_weight != 0.0:
            self.multiple_ctrl = True
            logging.info("#3 One weighted control used...")
        #-----------------------------------------------------------------

        self.end_shift = end_shift
        self.lambda_bg = lambda_bg

        if not self.ctrl_d_s or not self.ctrl_scaling_factor_s:
            self.no_lambda_flag = True
        else:
            self.no_lambda_flag = False

        self.pseudocount = pseudocount

        chr1 = set(self.treat.get_chr_names())

        #Changed to take into account multiple control -- WACS
        if not isinstance(self.ctrl, list):
            chr2 = set(self.ctrl.get_chr_names())
        else:
            chr2 = set()
            for item in self.ctrl:
                chr2 = chr2.union(set(item.get_chr_names()))

        self.chromosomes = list(chr1.intersection(chr2))
        #----------------------------------------------------
        self.pileup_data_files = {}

    cpdef destroy ( self ):
        """Remove temparary files for pileup values of each chromosome.

        Note: This function MUST be called if the class object won't
        be used anymore.

        """
        cdef:
            str f

        for f in self.pileup_data_files.values():
            if os.path.isfile( f ):
                os.unlink( f )
        return

    cpdef set_pseudocount( self, float pseudocount ):
        self.pseudocount = pseudocount

    cpdef enable_trackline( self ):
        """Turn on trackline with bedgraph output
        """
        self.trackline = True


    cdef __pileup_control_weighted_chrom(self, chrom):
        """
            Used to find pileups for weighted multiple controls

            This function is used to compute the pileup for the group of controls:
            - There are 3 scaling factors for each control (sf1, sf2, sf3).
            - For each scaling factor, compute
                = max( W1*DScale1*Spread(Control1,DRegion) + W2*DScale2*Spread(Control2,DRegion),
                       LambdaBG )

            self.ctrl, self.ctrl_weight, chrom, self.ctrl_d_s, self.ctrl_scaling_factor_s, baseline_value = self.lambda_bg, directional = False

            Returns pileup per control per scale factor for chrom
        """
        cdef:
            list pileup_ctrl
            long i
            float t
            object f

        assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom

        pileup_sf = [] #pileup for scalefactor

        for sf_index in range(len(self.ctrl_d_s)):
            index_ctrl = 0
            pileup_ctrl = None #whilst looping through controls one at a time, find max pileup. however n
            ctrl_d_s_index = self.ctrl_d_s[sf_index]
            sf_ctrl = self.ctrl_scaling_factor_s[sf_index]

            if self.ctrl_weight == 0:
                continue

            #print("ctrl_d_s_index: " + str(ctrl_d_s_index) + ", Scale factor: " + str(sf_ctrl) + ", Weight: " + str(self.ctrl_weight))
            if self.PE_mode:
                pileup_ctrl = self.ctrl.pileup_a_chromosome_c( chrom, [ctrl_d_s_index,],
                                                                    [sf_ctrl,],
                                                                    self.ctrl_weight,
                                                                    directional = False )                
            else:
                pileup_ctrl = self.ctrl.pileup_a_chromosome( chrom, [ctrl_d_s_index,],
                                                                    [sf_ctrl,],
                                                                    self.ctrl_weight,
                                                                    directional = False )

            #Pileups for corresponding sf saved
            pileup_sf.append(pileup_ctrl)

        #returns pileups for each SF per chrom
        return  pileup_sf

    cdef __pileup_treatment_weighted_chrom(self, chrom):
        """
            Find pileup for treatment only.
            Used to find pileups for weighted multiple controls, per chromosome. And save it!

            Returns [p,v] arrays for treatment.
        """
        cdef:
            list treat_pv
            long i
            float t
            object f

        assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom

        if self.PE_mode:
            treat_pv = self.treat.pileup_a_chromosome ( chrom, [self.treat_scaling_factor,], baseline_value = 0.0 )
        else:
            treat_pv = self.treat.pileup_a_chromosome( chrom, [self.d,], [self.treat_scaling_factor,], baseline_value = 0.0,
                                                       directional = True,
                                                       end_shift = self.end_shift )
        return treat_pv

    cpdef compute_pileup__control(self):

        """After this function is called, pileup per chromosome per control/treatment is computed.
        All chromosomes will be iterated. So it will take some time.

        """
        cdef:
            str chrom
            dict pileup_control_chrom

        logging.debug ( "Start to pileup for control one chromosome at a time per scale factor..." )

        #pileup_control_chrom = []
        pileup_control_chrom = {}
        for i in range(len( self.chromosomes ) ):
            chrom = self.chromosomes[ i ]
            pileup_control_chrom[chrom] =  self.__pileup_control_weighted_chrom(chrom)

        return pileup_control_chrom

    cpdef compute_pileup__treatment(self):

        """After this function is called, pileup per chromosome per control/treatment is computed.
        All chromosomes will be iterated. So it will take some time.

        """
        cdef:
            str chrom
            dict pileup_treat_chrom

        logging.debug ( "Start to compute pileup for treatment one chromosome at a time..." )

        #pileup_treat_chrom = []
        pileup_treat_chrom = {}

        for i in range(len( self.chromosomes ) ):
            chrom = self.chromosomes[ i ]
            pileup_treat_chrom[chrom] = self.__pileup_treatment_weighted_chrom(chrom)

        #a list of pileups per chromosome for the treatment
        return pileup_treat_chrom
