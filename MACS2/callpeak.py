# cython: profile=True
# Time-stamp: <2012-04-13 16:45:26 Tao Liu>

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

import sys
import logging

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
    
    #1 Read tag files
    info("#1 read tag files...")
    (treat, control) = load_tag_files_options (options)
    # check common chromosome names
    if control:
        tchrnames = set(treat.get_chr_names())
        cchrnames = set(control.get_chr_names())
        commonnames = tchrnames.intersection(cchrnames)
        if len(commonnames)==0:
            error("No common chromosome names can be found from treatment and control! Check your input files! MACS will quit...")
            error("Chromosome names in treatment: %s" % ",".join(sorted(tchrnames)))
            error("Chromosome names in control: %s" % ",".join(sorted(cchrnames)))
            sys.exit()
    
    info("#1 tag size = %d" % options.tsize)
    tagsinfo  = "# tag size is determined as %d bps\n" % (options.tsize)

    t0 = treat.total
    tagsinfo += "# total tags in treatment: %d\n" % (t0)
    info("#1  total tags in treatment: %d" % (t0))
    if options.keepduplicates != "all":
        if options.keepduplicates == "auto":
            info("#1 calculate max duplicate tags in single position based on binomal distribution...")
            treatment_max_dup_tags = cal_max_dup_tags(options.gsize,t0)
            info("#1  max_dup_tags based on binomal = %d" % (treatment_max_dup_tags))
            info("#1 filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)" % (treatment_max_dup_tags))
        else:
            info("#1 user defined the maximum tags...")
            treatment_max_dup_tags = int(options.keepduplicates)
            info("#1 filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)" % (treatment_max_dup_tags))

        treat.filter_dup(treatment_max_dup_tags)
        t1 = treat.total
        info("#1  tags after filtering in treatment: %d" % (t1))
        tagsinfo += "# tags after filtering in treatment: %d\n" % (t1)
        tagsinfo += "# maximum duplicate tags at the same position in treatment = %d\n" % (treatment_max_dup_tags)
        info("#1  Redundant rate of treatment: %.2f" % (float(t0-t1)/t0))
        tagsinfo += "# Redundant rate in treatment: %.2f\n" % (float(t0-t1)/t0)
    else:
        t1 = treat.total

    if control:
        c0 = control.total
        tagsinfo += "# total tags in control: %d\n" % (c0)
        info("#1  total tags in control: %d" % (c0))
        if options.keepduplicates != "all":
            if options.keepduplicates == "auto":
                info("#1  for control, calculate max duplicate tags in single position based on binomal distribution...")
                control_max_dup_tags = cal_max_dup_tags(options.gsize,c0)
                info("#1  max_dup_tags based on binomal = %d" % (control_max_dup_tags))
                info("#1 filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)" % (control_max_dup_tags))
            else:
                info("#1 user defined the maximum tags...")
                control_max_dup_tags = int(options.keepduplicates)
                info("#1 filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)" % (treatment_max_dup_tags))
            
            control.filter_dup(control_max_dup_tags)
            c1 = control.total
            info("#1  tags after filtering in control: %d" % (c1))
            tagsinfo += "# tags after filtering in control: %d\n" % (c1)
            tagsinfo += "# maximum duplicate tags at the same position in control = %d\n" % (control_max_dup_tags)
            
            info("#1  Redundant rate of control: %.2f" % (float(c0-c1)/c0))
            tagsinfo += "# Redundant rate in control: %.2f\n" % (float(c0-c1)/c0)
        else:
            c1 = control.total
        
            
    info("#1 finished!")

    #2 Build Model
    info("#2 Build Peak Model...")

    if options.nomodel:
        info("#2 Skipped...")
        options.d=options.shiftsize*2
        info("#2 Use %d as shiftsize, %d as fragment length" % (options.shiftsize,options.d))
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
                info("#3 MACS is random sampling treatment tags...")
                warn("#3 Your results may not be reproducible due to the random sampling!")
                treat.sample_num(c1)
                info("#3 %d tags from treatment are kept" % treat.total)                
            elif c1 > t1: 
                info("#3 MACS is random sampling control tags...")
                warn("#3 Your results may not be reproducible due to the random sampling!")                
                control.sample_num(t1)
                info("#3 %d tags from control are kept" % control.total)
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
    #diag_result = peakdetect.diag_result()
    #4 output
    #4.1 peaks in XLS
    info("#4 Write output xls file... %s" % (options.peakxls))
    ofhd_xls = open(options.peakxls,"w")
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
    ofhd_xls.write(peakdetect.toxls())
    ofhd_xls.close()
    #4.2 peaks in BED
    if options.log_pvalue:
        score_column = "pscore"
    elif options.log_qvalue:
        score_column = "qscore"
    info("#4 Write peak bed file... %s" % (options.peakbed))
    ofhd_bed = open(options.peakbed,"w")
    peakdetect.peaks.write_to_bed (ofhd_bed, name_prefix="MACS_peak_", score_column=score_column)
    ofhd_bed.close()
    #4.2 peaks in narrowPeak
    info("#4 Write peak in narrowPeak format file... %s" % (options.peakNarrowPeak))
    ofhd_bed = open(options.peakNarrowPeak,"w")
    peakdetect.peaks.write_to_narrowPeak (ofhd_bed, name_prefix="MACS_peak_", score_column=score_column)
    ofhd_bed.close()
    #4.2 broad peaks in bed12
    if options.broad:
        info("#4 Write broad peak in bed12 format file... %s" % (options.peakBroadPeak))
        ofhd_bed = open(options.peakBroadPeak,"w")
        peakdetect.broadpeaks.write_to_gappedPeak (ofhd_bed, name_prefix="MACS_peak_", name=options.name, description=options.name)
        ofhd_bed.close()
    #4.2-2 summits in BED
    info("#4 Write summits bed file... %s" % (options.summitbed))
    ofhd_summits = open(options.summitbed,"w")
    peakdetect.peaks.write_to_summit_bed (ofhd_summits, name_prefix="MACS_summit_", score_column=score_column)
    ofhd_summits.close()

def cal_max_dup_tags ( genome_size, tags_number, p=1e-5 ):
    """Calculate the maximum duplicated tag number based on genome
    size, total tag number and a p-value based on binomial
    distribution. Brute force algorithm to calculate reverse CDF no
    more than MAX_LAMBDA(100000).
    
    """
    return binomial_cdf_inv(1-p,tags_number,1.0/genome_size)

def load_tag_files_options ( options ):
    """From the options, load treatment tags and control tags (if available).

    """
    options.info("#1 read treatment tags...")
    tp = options.parser(options.tfile)
        
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize

    treat = tp.build_fwtrack()
    treat.sort()
    if options.cfile:
        options.info("#1.2 read input tags...")
        control = options.parser(options.cfile).build_fwtrack()
        control.sort()
    else:
        control = None
    options.info("#1 tag size is determined as %d bps" % options.tsize)
    return (treat, control)

