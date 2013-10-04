# Time-stamp: <2013-09-11 17:39:20 Tao Liu>

"""Description: MACS 2 main executable

Copyright (c) 2008,2009 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status: release candidate
@version: $Id$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
from time import strftime

# ------------------------------------
# own python modules
# ------------------------------------
from MACS2.OptValidator import opt_validate
from MACS2.OutputWriter import *
from MACS2.cProb import binomial_cdf_inv
from MACS2.cPeakModel import PeakModel,NotEnoughPairsException
from MACS2.cPeakDetect import PeakDetect
from MACS2.Constants import *
# ------------------------------------
# Main function
# ------------------------------------
def check_names(treat, control, error_stream):
    """check common chromosome names"""
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
    
    """
    # Parse options...
    options = opt_validate( args )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    #0 output arguments
    info("\n"+options.argtxt)
    options.PE_MODE = options.format in ('BAMPE',)
    if options.PE_MODE: tag = 'fragment' # call things fragments not tags
    else: tag = 'tag'
    
    #1 Read tag files
    info("#1 read %s files...", tag)
    if options.PE_MODE: (treat, control) = load_frag_files_options (options)
    else:       (treat, control) = load_tag_files_options  (options)
    if control is not None: check_names(treat, control, error)
    
    info("#1 %s size = %d", tag, options.tsize)
    tagsinfo  = "# %s size is determined as %d bps\n" % (tag, options.tsize)
    
    t0 = treat.total
    tagsinfo += "# total %ss in treatment: %d\n" % (tag, t0)
    info("#1  total %ss in treatment: %d", tag, t0)
    # not ready yet
#    options.filteringmodel = True
#    if options.filteringmodel:
#        treat.separate_dups()
#        t0 = treat.total + treat.dups.total
#        t1 = treat.total
#        info("#1  Redundant rate of treatment: %.2f", float(t0 - t1) / t0)
#        tagsinfo += "# Redundant rate in treatment: %.2f\n" % (float(t0-t1)/t0)
#    elif options.keepduplicates != "all":
    if options.keepduplicates != "all":
        if options.keepduplicates == "auto":
            info("#1 calculate max duplicate %ss in single position based on binomial distribution...", tag)
            treatment_max_dup_tags = cal_max_dup_tags(options.gsize,t0)
            info("#1  max_dup_tags based on binomial = %d" % (treatment_max_dup_tags))
        else:
            info("#1 user defined the maximum %ss...", tag)
            treatment_max_dup_tags = int(options.keepduplicates)
        if options.PE_MODE:
            info("#1 filter out redundant fragments by allowing at most %d identical fragment(s)", treatment_max_dup_tags)
        else:
            info("#1 filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)", treatment_max_dup_tags)
        treat.separate_dups(treatment_max_dup_tags) # changed 5-29
#        treat.filter_dup(treatment_max_dup_tags)
        t1 = treat.total
        info("#1  %ss after filtering in treatment: %d", tag, t1)
        tagsinfo += "# %ss after filtering in treatment: %d\n" % (tag, t1)
        if options.PE_MODE:
            tagsinfo += "# maximum duplicate fragments in treatment = %d\n" % (treatment_max_dup_tags)
        else:
            tagsinfo += "# maximum duplicate tags at the same position in treatment = %d\n" % (treatment_max_dup_tags)
        info("#1  Redundant rate of treatment: %.2f", float(t0 - t1) / t0)
        tagsinfo += "# Redundant rate in treatment: %.2f\n" % (float(t0-t1)/t0)
    else:
        t1 = t0

    if control is not None:
        c0 = control.total
        tagsinfo += "# total %ss in control: %d\n" % (tag, c0)
        info("#1  total %ss in control: %d", tag, c0)
        # not ready yet
        #if options.filteringmodel:
        #    control.separate_dups()
        #    c0 = treat.total + treat.dups.total
        #    c1 = treat.total
        #    info("#1  Redundant rate of treatment: %.2f", float(c0 - c1) / c0)
        #    tagsinfo += "# Redundant rate in treatment: %.2f\n" % (float(c0-c1)/c0)
        #elif options.keepduplicates != "all":
        if options.keepduplicates != "all":
            if options.keepduplicates == "auto":
                info("#1  for control, calculate max duplicate %ss in single position based on binomial distribution...", tag)
                control_max_dup_tags = cal_max_dup_tags(options.gsize,c0)
                info("#1  max_dup_tags based on binomial = %d" % (control_max_dup_tags))
            else:
                info("#1 user defined the maximum %ss...", tag)
                control_max_dup_tags = int(options.keepduplicates)
            if options.PE_MODE:
                info("#1 filter out redundant fragments by allowing at most %d identical fragment(s)", treatment_max_dup_tags)
            else:
                info("#1 filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)", treatment_max_dup_tags)
#            control.filter_dup(treatment_max_dup_tags)
            control.separate_dups(treatment_max_dup_tags) # changed 5-29
            c1 = control.total
            
            info("#1  %ss after filtering in control: %d", tag, c1)
            tagsinfo += "# %ss after filtering in control: %d\n" % (tag, c1)
            if options.PE_MODE:
                tagsinfo += "# maximum duplicate fragments in control = %d\n" % (treatment_max_dup_tags)
            else:
                tagsinfo += "# maximum duplicate tags at the same position in control = %d\n" % (treatment_max_dup_tags)
            
            info("#1  Redundant rate of control: %.2f" % (float(c0-c1)/c0))
            tagsinfo += "# Redundant rate in control: %.2f\n" % (float(c0-c1)/c0)
        else:
            c1 = c0
    info("#1 finished!")

    #2 Build Model
    info("#2 Build Peak Model...")

    if options.nomodel:
        info("#2 Skipped...")
        if options.PE_MODE:
            options.shiftsize = 0
            options.d = options.tsize
        else:
            options.d=options.shiftsize*2
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
                    options.d=options.shiftsize*2
                    options.scanwindow=2*options.d 
                    warn("#2 MACS will use %d as shiftsize, %d as fragment length. NOTE: if the d calculated is still acceptable, please do not use --fix-bimodal option!" % (options.shiftsize,options.d))
                else:
                    warn("#2 You may need to consider one of the other alternative d(s): %s" %  ','.join(map(str,peakmodel.alternative_d)))
                    warn("#2 You can restart the process with --nomodel --shiftsize XXX with your choice or an arbitrary number. Nontheless, MACS will continute computing.")
                
        except NotEnoughPairsException:
            if not options.onauto:
                sys.exit(1)
            warn("#2 Skipped...")
            options.d=options.shiftsize*2
            options.scanwindow=2*options.d 
            warn("#2 Since --fix-bimodal is set, MACS will use %d as shiftsize, %d as fragment length" % (options.shiftsize,options.d))


    #3 Call Peaks
    info("#3 Call peaks...")
    if options.nolambda:
        info("# local lambda is disabled!")

    # decide options.tocontrol according to options.tolarge
    if control:
        if options.downsample:
            # use random sampling to balance treatment and control
            info("#3 User prefers to use random sampling instead of linear scaling.")
            if t1 > c1:
                info("#3 MACS is random sampling treatment %ss...", tag)
                if options.seed < 0:
                    warn("#3 Your results may not be reproducible due to the random sampling!")
                else:
                    info("#3 Random seed (%d) is used." % options.seed)
                treat.sample_num(c1, options.seed)
                info("#3 %d Tags from treatment are kept", treat.total)                
            elif c1 > t1: 
                info("#3 MACS is random sampling control %ss...", tag)
                if options.seed < 0:
                    warn("#3 Your results may not be reproducible due to the random sampling!")
                else:
                    info("#3 Random seed (%d) is used." % options.seed)
                control.sample_num(t1, options.seed)
                info("#3 %d %ss from control are kept", control.total, tag)
            # set options.tocontrol although it would;t matter now
            options.tocontrol = False
        else:
            if options.tolarge:
                if t1 > c1:
                    # treatment has more tags than control, since tolarge is
                    # true, we will scale control to treatment.
                    options.tocontrol = False
                else:
                    # treatment has less tags than control, since tolarge is
                    # true, we will scale treatment to control.
                    options.tocontrol = True
            else:
                if t1 > c1:
                    # treatment has more tags than control, since tolarge is
                    # false, we will scale treatment to control.
                    options.tocontrol = True
                else:
                    # treatment has less tags than control, since tolarge is
                    # false, we will scale control to treatment.
                    options.tocontrol = False

    peakdetect = PeakDetect(treat = treat,
                            control = control,
                            opt = options
                            )
    peakdetect.call_peaks()

    #call refinepeak if needed.
    # if options.refine_peaks:
    #     info("#3 now put back duplicate reads...")
    #     treat.addback_dups()
    #     info("#3 calculate reads balance to refine peak summits...")
    #     refined_peaks = treat.refine_peak_from_tags_distribution ( peakdetect.peaks, options.d, 0 )
    #     info("#3 reassign scores for newly refined peak summits...")
    #     peakdetect.peaks = peakdetect.scoretrack.reassign_peaks( refined_peaks ) # replace
    #     #info("#3 write to file: %s ..." % options.name+"_refined_peaks.encodePeak" )
    #     #refinedpeakfile = open(options.name+"_refined_peaks.encodePeak", "w")
    #     #refined_peaks.write_to_narrowPeak (refinedpeakfile, name_prefix="%s_refined_peak_", name=options.name, score_column=score_column, trackline=options.trackline )

    #diag_result = peakdetect.diag_result()
    #4 output
    #4.1 peaks in XLS
    info("#4 Write output xls file... %s" % (options.peakxls))
    ofhd_xls = open( os.path.join( options.outdir, options.peakxls), "w" )
    ofhd_xls.write("# This file is generated by MACS version %s\n" % (MACS_VERSION))
    ofhd_xls.write(options.argtxt+"\n")

    ofhd_xls.write(tagsinfo)

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
    if options.log_pvalue:
        score_column = "pscore"
    elif options.log_qvalue:
        score_column = "qscore"
    #4.2 peaks in narrowPeak
    if not options.broad:
        #info("#4 Write peak bed file... %s" % (options.peakbed))
        #ofhd_bed = open(options.peakbed,"w")
        #peakdetect.peaks.write_to_bed (ofhd_bed, name_prefix="%s_peak_", name = options.name, description="Peaks for %s (Made with MACS v2, " + strftime("%x") + ")", score_column=score_column, trackline=options.trackline)
        #ofhd_bed.close()
        info("#4 Write peak in narrowPeak format file... %s" % (options.peakNarrowPeak))
        ofhd_bed = open( os.path.join( options.outdir, options.peakNarrowPeak), "w" )
        peakdetect.peaks.write_to_narrowPeak (ofhd_bed, name_prefix="%s_peak_", name=options.name, score_column=score_column, trackline=options.trackline )
        ofhd_bed.close()
        #4.2-2 summits in BED
        info("#4 Write summits bed file... %s" % (options.summitbed))
        ofhd_summits = open( os.path.join( options.outdir, options.summitbed), "w" )
        peakdetect.peaks.write_to_summit_bed (ofhd_summits, name_prefix="%s_peak_", name=options.name,
                                              description="Summits for %s (Made with MACS v2, " + strftime("%x") + ")",
                                              score_column=score_column, trackline=options.trackline )
        ofhd_summits.close()
    #4.2 broad peaks in bed12 or gappedPeak
    else:
        info("#4 Write broad peak in broadPeak format file... %s" % (options.peakBroadPeak))
        ofhd_bed = open( os.path.join( options.outdir, options.peakBroadPeak), "w" )
        peakdetect.peaks.write_to_broadPeak (ofhd_bed, name_prefix="%s_peak_", name=options.name, description=options.name, trackline=options.trackline)
        ofhd_bed.close()
        info("#4 Write broad peak in bed12/gappedPeak format file... %s" % (options.peakGappedPeak))
        ofhd_bed = open( os.path.join( options.outdir, options.peakGappedPeak), "w" )
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

def load_frag_files_options ( options ):
    """From the options, load treatment fragments and control fragments (if available).

    """
    options.info("#1 read treatment fragments...")

    tp = options.parser(options.tfile[0], buffer_size=options.buffer_size)
    treat = tp.build_petrack()
    treat.sort()
    if len(options.tfile) > 1:
        # multiple input
        for tfile in options.tfile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_petrack( treat )
            treat.sort()

    options.tsize = tp.d
    if options.cfile:
        options.info("#1.2 read input fragments...")
        cp = options.parser(options.cfile[0], buffer_size=options.buffer_size)
        control = cp.build_petrack()
        control_d = cp.d
        control.sort()
        if len(options.cfile) > 1:
            # multiple input
            for cfile in options.cfile[1:]:
                cp = options.parser(cfile, buffer_size=options.buffer_size)
                control = cp.append_petrack( control )
                control.sort()
    else:
        control = None
    options.info("#1 mean fragment size is determined as %d bp from treatment" % options.tsize)
#    options.info("#1 fragment size variance is determined as %d bp from treatment" % tp.variance)
    if control is not None:
        options.info("#1 note: mean fragment size in control is %d bp -- value ignored" % control_d)
    return (treat, control)

def load_tag_files_options ( options ):
    """From the options, load treatment tags and control tags (if available).

    """
    options.info("#1 read treatment tags...")
    tp = options.parser(options.tfile[0], buffer_size=options.buffer_size)
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize
    treat = tp.build_fwtrack()
    treat.sort()
    if len(options.tfile) > 1:
        # multiple input
        for tfile in options.tfile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )
            treat.sort()
    
    if options.cfile:
        options.info("#1.2 read input tags...")
        control = options.parser(options.cfile[0], buffer_size=options.buffer_size).build_fwtrack()
        control.sort()
        if len(options.cfile) > 1:
            # multiple input
            for cfile in options.cfile[1:]:
                cp = options.parser(cfile, buffer_size=options.buffer_size)
                control = cp.append_fwtrack( control )
                control.sort()
    else:
        control = None
    options.info("#1 tag size is determined as %d bps" % options.tsize)
    return (treat, control)
