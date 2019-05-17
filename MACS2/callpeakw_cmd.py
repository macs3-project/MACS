"""Description: MACS 2 main executable

Copyright (c) 2008,2009 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2013, Tao Liu lab at UB and Xiaole Shirley Liu lab at DFCI
All rights reserved.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status: release candidate
@version: $Id$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

"""Added: WACS main executable by Aseel Awdeh <araed104@uottawa.ca>
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
from time import strftime
import tempfile

# ------------------------------------
# own python modules
# ------------------------------------
from MACS2.OptValidator import opt_validate_wacs
from MACS2.OutputWriter import *
from MACS2.Prob import binomial_cdf_inv
from MACS2.PeakModel import PeakModel,NotEnoughPairsException
from MACS2.PeakDetect import PeakDetect
from MACS2.ComputePileups import ComputePileup
from MACS2.Constants import *
from MACS2.Pileup import max_over_two_pv_array, sum_over_two_pv_array, multiple_pv_array_w_weight, print_pv

# ------------------------------------
# Main function
# ------------------------------------
def check_names(treat, control, error_stream, options):
    """
        check common chromosome names -- need to change..
    """
    tchrnames = set(treat.get_chr_names())
    cchrnames = set(control.get_chr_names())
    commonnames = tchrnames.intersection(cchrnames)
    if len(commonnames)==0:
        error_stream("No common chromosome names can be found from treatment and control! Check your input files! MACS will quit...")
        error_stream("Chromosome names in treatment: %s" % ",".join(sorted(tchrnames)))
        error_stream("Chromosome names in control: %s" % ",".join(sorted(cchrnames)))
        sys.exit()

def run( args ):
    """The Main function/pipeline for MACS.
        This is for the weighted option.
    """
    # Parse options...
    options = opt_validate_wacs( args )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    #0 output arguments
    info("\n"+options.argtxt)
    options.PE_MODE = options.format in ('BAMPE','BEDPE')
    if options.PE_MODE: tag = 'fragment' # call things fragments not tags
    else: tag = 'tag'

    tempfile.tempdir = options.tempdir

    #-------------------------------------------------------------------------------------------------------------------
    #-------------------------------------#1 Read tag files-------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------
    info("#1 read %s files...", tag)

    #Load target first -- only once
    if options.PE_MODE: treat = load_frag_treat_files_options (options)
    else:       treat = load_tag_treat_files_options  (options)

    info("#1 %s size = %d", tag, options.tsize)
    tagsinfo  = "# %s size is determined as %d bps\n" % (tag, options.tsize)
    t0 = treat.total
    tagsinfo += "# total %ss in treatment: %d\n" % (tag, t0)
    info("#1  total %ss in treatment: %d", tag, t0)
    info("#1 Redundant reads are not removed. With multiple weighted controls, we assume that files are filtered.")
    t1 = t0

    #-------------------------------------------------------------------------------------------------------------------
    #2 Build Model -- only depends on treatment -- so build model before looping through each control. -- done only once.
    #-------------------------------------------------------------------------------------------------------------------
    info("#2 Build Peak Model...")
    if options.nomodel:
        info("#2 Skipped...")
        if options.PE_MODE:
            #options.shiftsize = 0
            options.d = options.tsize
        else:
            options.d=options.extsize
        if options.shift > 0:
            info("#2 Sequencing ends will be shifted towards 3' by %d bp(s)" % (options.shift))
        elif options.shift < 0:
            info("#2 Sequencing ends will be shifted towards 5' by %d bp(s)" % (options.shift * -1))
        info("#2 Use %d as fragment length" % (options.d))
        options.scanwindow=2*options.d  # remove the effect of --bw
    else:
        try:
            peakmodel = PeakModel(treatment = treat,
                                  max_pairnum = MAX_PAIRNUM,
                                  opt = options
                                  )
            info("#2 finished!")
            debug("#2  Summary Model:")
            debug("#2   min_tags: %d" % (peakmodel.min_tags))
            debug("#2   d: %d" % (peakmodel.d))
            debug("#2   scan_window: %d" % (peakmodel.scan_window))
            info("#2 predicted fragment length is %d bps" % peakmodel.d)
            info("#2 alternative fragment length(s) may be %s bps" % ','.join(map(str,peakmodel.alternative_d)))
            info("#2.2 Generate R script for model : %s" % (options.modelR))
            model2r_script(peakmodel,options.modelR,options.name)
            options.d = peakmodel.d
            options.scanwindow= 2*options.d
            if options.d <= 2*options.tsize:
                warn("#2 Since the d (%.0f) calculated from paired-peaks are smaller than 2*tag length, it may be influenced by unknown sequencing problem!" % (options.d))
                if options.onauto:
                    options.d=options.extsize
                    options.scanwindow=2*options.d
                    warn("#2 MACS will use %d as EXTSIZE/fragment length d. NOTE: if the d calculated is still acceptable, please do not use --fix-bimodal option!" % (options.d))
                else:
                    warn("#2 You may need to consider one of the other alternative d(s): %s" %  ','.join(map(str,peakmodel.alternative_d)))
                    warn("#2 You can restart the process with --nomodel --extsize XXX with your choice or an arbitrary number. Nontheless, MACS will continute computing.")

        except NotEnoughPairsException:
            if not options.onauto:
                sys.exit(1)
            warn("#2 Skipped...")
            options.d=options.extsize
            options.scanwindow=2*options.d
            warn("#2 Since --fix-bimodal is set, MACS will use %d as fragment length" % (options.d))

    #Get chromosomes of treatment
    chr1 = set(treat.get_chr_names())

    #-------------------------------------------------------------------------------------------------------------------
    #------------------------------------Loop through the controls-----------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------
    #Then, load controls one a time.
    if options.cfile:
        prev_pileup_ctrl = None
        treat_pileup = True
        index_ctrl = 0
        chr2 = set()

        #Loop through controls and compute pileups for controls one at a time
        #Step 1. Compute pileups
                #1.1 Compute treatment pileup (once)
                #1.2 Compute control pileup
                #1.3 Sum control pileups
        #Step 2. Replace any pileup value less than lambda with lambda
        #Step 3. Callpeaks
                #3.1
                #self.__chrom_pair_treat_ctrl( treat_pv, ctrl_pv) for each chromosome
                #Construct p-value q-value table for each chromsome
                #self.__chrom_pair_treat_ctrl( treat_pv, ctrl_pv) for each chromosome
                #Peak calling for each chromosome

        for cfile in options.cfile:
            if options.controlweights[index_ctrl] == 0: #dont read controls with weights 0
                index_ctrl = index_ctrl + 1
                continue
            if len(options.cfile) == 1: #if there is only one control, give it a weight of 1 and ignore the weight it already has.
                options.controlweights = [1]

            options.controlweight = options.controlweights[index_ctrl]
            options.info("#3.1 Read multiple input tags for control " +  str(index_ctrl) + " with weight " + str(options.controlweight) + " ...")

            if options.PE_MODE:
                control = options.parser(cfile, buffer_size=options.buffer_size).build_petrack()
            else:
                control = options.parser(cfile, buffer_size=options.buffer_size).build_fwtrack()

            control.finalize() #called from FixedWidthTrack.pyx or PairedEndTrack.pyx

            if control is not None: check_names(treat, control, error, options) #check overlap in chromsomes between control and treatment
            chr2 = chr2.union(set(control.get_chr_names()))

            c0 = control.total
            tagsinfo += "# total %ss in control: %d\n" % (tag, c0)
            info("#3.1  Total %ss in control: %d", tag, c0)
            info("#3.1 Finished for control: %d", index_ctrl)
            info("#3.1 Peak Model already computed.")

            #Compute pileups
            info("#3.2 Compute pileups for control %d", index_ctrl)
            #In case of multiple controls, always scale towards the treatment.. regardless of the controls treatment has less tags than control, since tolarge is false, we will scale control to treatment.
            options.tocontrol = False
            #Send only one control with that weight
            pileup = ComputePileup(treat = treat,
                                   control = control,
                                   opt = options
                                   )

            if treat_pileup:
                print("#3.2 Compute pileup for treatment.")
                pileup_treatment = pileup.pileup_treatment()
                treat_pileup = False

            temp_pileup_ctrl = pileup.pileup_control() #in format, pileups per chrom per scaling factor (24 chromosomes * 3 pileups)

            #prev_pileup_ctrl has pileups for each scale factors.
            #Sum pileups of each sf for for corresponding two controls
            if prev_pileup_ctrl:
                for chr, pileups in prev_pileup_ctrl.iteritems(): #loop through each chr
                    prev_pileup_ctrl_chr = prev_pileup_ctrl.get(chr)
                    temp_pileup_ctrl_chr = temp_pileup_ctrl.get(chr)
                    pileup_chrom_sf = []
                    for sf in range(len(pileups)): #loop through each SF
                        pileup_chrom_sf.append( sum_over_two_pv_array(prev_pileup_ctrl_chr[sf], temp_pileup_ctrl_chr[sf]) )
                    prev_pileup_ctrl[chr] = pileup_chrom_sf

            else:
                prev_pileup_ctrl = temp_pileup_ctrl

            index_ctrl = index_ctrl + 1


        #After looping through all controls, and find control and treatment pileups,
        chromosomes = list(chr1.intersection(chr2))
        treat_sum = t0 * options.d
        lambda_bg = float( treat_sum )/ options.gsize

        #Loop through control pileup if there are any positions with value less than lambda, replace them with lambda
        print("#3.3 Replace values in control pileup less than lambda.")
        for chr, pileups in prev_pileup_ctrl.iteritems():
            prev_pileup_ctrl_chr = prev_pileup_ctrl.get(chr) #3 pileups corresponding to each scale factor
            pileup_chrom_sf = []
            for sf in range(len(pileups)): #loop through each control
                [pos_, value_] = prev_pileup_ctrl_chr[sf]
                for index in range(len(pos_)):
                    if value_[index] < lambda_bg: #If any pileup value less than background, replace with background.
                        value_[index] = lambda_bg
                pileup_chrom_sf.append( [pos_, value_] )
            prev_pileup_ctrl[chr] = pileup_chrom_sf


        #WACS_Combined = max( WACS_D, WACS_S, WACS_L )
        #Now after the summed pileups of controls per scale factor is found, find max of pileups for diff scale factors
        print("#3.4 Find max control pileup per chromosome for different scale factors.")
        for chr, pileups in prev_pileup_ctrl.iteritems():
            prev_pileup_ctrl_chr = prev_pileup_ctrl.get(chr)
            prev_pileup_sf = []
            for sf in range(len(pileups)):
                if prev_pileup_sf:
                    prev_pileup_sf = max_over_two_pv_array ( prev_pileup_sf, prev_pileup_ctrl_chr[sf] )
                else:
                    prev_pileup_sf = prev_pileup_ctrl_chr[sf]

            prev_pileup_ctrl[chr] = prev_pileup_sf


        info("#4 Call peaks ...")
        #Now we have control pileup and treatment pileup -- prev_pileup_sf, pileup_treatment per chromosome
        options.multiple = True
        peakdetect = PeakDetect(treat = treat,
                                control = control,
                                opt = options
                                )

        peakdetect.set_ctrlpileup(pileup_treatment, prev_pileup_ctrl, chromosomes)
        peakdetect.call_peaks()
        # filter out low FE peaks
        #peakdetect.peaks.filter_fc( fc_low = options.fecutoff )

    else:
        info("For the weighted option, there should one or more controls.")
        sys.exit(1)

    #-------------------------------------------------------------------------------------------------------------------
    #4 output
    #4.1 peaks in XLS
    #-------------------------------------------------------------------------------------------------------------------
    info("#4 Write output xls file... %s" % (options.peakxls))
    ofhd_xls = open( options.peakxls, "w" )
    ofhd_xls.write("# This file is generated by MACS version %s\n" % (MACS_VERSION))
    ofhd_xls.write(options.argtxt+"\n")

    ofhd_xls.write(tagsinfo)

    if options.shift > 0:
        ofhd_xls.write("# Sequencing ends will be shifted towards 3' by %d bp(s)\n" % (options.shift))
    elif options.shift < 0:
        ofhd_xls.write("# Sequencing ends will be shifted towards 5' by %d bp(s)\n" % (options.shift * -1))

    ofhd_xls.write("# d = %d\n" % (options.d))
    try:
        ofhd_xls.write("# alternative fragment length(s) may be %s bps\n" % ','.join(map(str,peakmodel.alternative_d)))
    except:
        # when --nomodel is used, there is no peakmodel object. Simply skip this line.
        pass
    if options.nolambda:
        ofhd_xls.write("# local lambda is disabled!\n")
    # pass write method so we can print too, and include name
    peakdetect.peaks.write_to_xls(ofhd_xls, name = options.name)
    ofhd_xls.close()
    #4.2 peaks in BED
    if options.log_pvalue != None:
        score_column = "pscore"
    elif options.log_qvalue != None:
        score_column = "qscore"
    #4.2 peaks in narrowPeak
    if not options.broad:
        #info("#4 Write peak bed file... %s" % (options.peakbed))
        #ofhd_bed = open(options.peakbed,"w")
        #peakdetect.peaks.write_to_bed (ofhd_bed, name_prefix="%s_peak_", name = options.name, description="Peaks for %s (Made with MACS v2, " + strftime("%x") + ")", score_column=score_column, trackline=options.trackline)
        #ofhd_bed.close()
        info("#4 Write peak in narrowPeak format file... %s" % (options.peakNarrowPeak))
        ofhd_bed = open( options.peakNarrowPeak, "w" )
        peakdetect.peaks.write_to_narrowPeak (ofhd_bed, name_prefix="%s_peak_", name=options.name, score_column=score_column, trackline=options.trackline )
        ofhd_bed.close()
        #4.2-2 summits in BED
        info("#4 Write summits bed file... %s" % (options.summitbed))
        ofhd_summits = open( options.summitbed, "w" )
        peakdetect.peaks.write_to_summit_bed (ofhd_summits, name_prefix="%s_peak_", name=options.name,
                                              description="Summits for %s (Made with MACS v2, " + strftime("%x") + ")",
                                              score_column=score_column, trackline=options.trackline )
        ofhd_summits.close()
    #4.2 broad peaks in bed12 or gappedPeak
    else:
        info("#4 Write broad peak in broadPeak format file... %s" % (options.peakBroadPeak))
        ofhd_bed = open( options.peakBroadPeak, "w" )
        peakdetect.peaks.write_to_broadPeak (ofhd_bed, name_prefix="%s_peak_", name=options.name, description=options.name, trackline=options.trackline)
        ofhd_bed.close()
        info("#4 Write broad peak in bed12/gappedPeak format file... %s" % (options.peakGappedPeak))
        ofhd_bed = open( options.peakGappedPeak, "w" )
        peakdetect.peaks.write_to_gappedPeak (ofhd_bed, name_prefix="%s_peak_", name=options.name, description=options.name, trackline=options.trackline)
        ofhd_bed.close()

    info("Done!")

def cal_max_dup_tags ( genome_size, tags_number, p=1e-5 ):
    """Calculate the maximum duplicated tag number based on genome
    size, total tag number and a p-value based on binomial
    distribution. Brute force algorithm to calculate reverse CDF no
    more than MAX_LAMBDA(100000).

    """
    return binomial_cdf_inv(1-p,tags_number,1.0/genome_size)

def load_tag_treat_files_options ( options ):
    """From the options, load treatment tags.
    """
    options.info("#1 read treatment tags...")
    tp = options.parser(options.tfile[0], buffer_size=options.buffer_size)
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize
    treat = tp.build_fwtrack()
    #treat.sort()
    if len(options.tfile) > 1:
        # multiple input
        for tfile in options.tfile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )
            #treat.sort()
    treat.finalize()

    options.info("#1 tag size is determined as %d bps" % options.tsize)
    return treat

def load_frag_treat_files_options ( options ):
    """From the options, load treatment fragments.

    """
    options.info("#1 read treatment fragments...")

    tp = options.parser(options.tfile[0], buffer_size=options.buffer_size)
    treat = tp.build_petrack()
    #treat.sort()
    if len(options.tfile) > 1:
        # multiple input
        for tfile in options.tfile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_petrack( treat )
            #treat.sort()
    treat.finalize()
    options.tsize = tp.d
    options.info("#1 mean fragment size is determined as %d bp from treatment" % options.tsize)
    return treat
