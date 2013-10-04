# Time-stamp: <2013-09-11 22:44:34 Tao Liu>

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
from MACS2.IO import cBedGraphIO
from MACS2.IO.cScoreTrack import DiffScoreTrackI
from MACS2.IO.cPeakIO import PeakIO
from MACS2.OptValidator import diff_opt_validate
from MACS2.OutputWriter import *
from MACS2.cProb import binomial_cdf_inv
from MACS2.cPeakModel import PeakModel,NotEnoughPairsException
from MACS2.cPeakDetect import PeakDetect
from MACS2.Constants import *
# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

# ------------------------------------
# Misc functions
# ------------------------------------
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info
# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """The Differential function/pipeline for MACS.
    
    """
    # Parse options...
    options = diff_opt_validate( args )
    #0 output arguments
#    info("\n"+options.argtxt)
 
    ofile_prefix = options.name
    
    # check if tag files exist
    with open(options.t1bdg) as f: pass
    with open(options.c1bdg) as f: pass
    with open(options.t2bdg) as f: pass
    with open(options.c2bdg) as f: pass
    
    if not options.peaks1 == '':
        info("Read peaks for condition 1...")
        p1io = PeakIO()
        with open(options.peaks1, 'rU') as f:
            p1io.read_from_xls(f)

    if not options.peaks2 == '':
        info("Read peaks for condition 2...")
        p2io = PeakIO()
        with open(options.peaks2, 'rU') as f:
            p2io.read_from_xls(f)
    
    #1 Read tag files
    info("Read and build treatment 1 bedGraph...")
    t1bio = cBedGraphIO.bedGraphIO(options.t1bdg)
    t1btrack = t1bio.build_bdgtrack()

    info("Read and build control 1 bedGraph...")
    c1bio = cBedGraphIO.bedGraphIO(options.c1bdg)
    c1btrack = c1bio.build_bdgtrack()

    if len(options.depth) >=2:
        depth1 = options.depth[0]
        depth2 = options.depth[1]
    else:
        depth1 = options.depth[0]
        depth2 = depth1
    
    info("Read and build treatment 2 bedGraph...")
    t2bio = cBedGraphIO.bedGraphIO(options.t2bdg)
    t2btrack = t2bio.build_bdgtrack()

    info("Read and build control 2 bedGraph...")
    c2bio = cBedGraphIO.bedGraphIO(options.c2bdg)
    c2btrack = c2bio.build_bdgtrack()
    
    #3 Call Peaks

    diffscore = DiffScoreTrackI( t1btrack,
                                 c1btrack,
                                 t2btrack,
                                 c2btrack,
                                 depth1, depth2 )
    diffscore.finalize()
    if options.call_peaks:
        diffscore.set_track_score_method(options.track_score_method)
        info("Calling peaks")
        if options.track_score_method == 'p':
            diffscore.call_peaks(cutoff = options.peaks_log_pvalue,
                                 min_length = options.pminlen)
        elif options.track_score_method == 'q':
            diffscore.call_peaks(cutoff = options.peaks_log_qvalue,
                                 min_length = options.pminlen)
        else:
            raise NotImplementedError
    else:
        info("Using existing peaks")
        diffscore.store_peaks(p1io, p2io)
        info("Rebuilding chromosomes")
        diffscore.rebuild_chromosomes()
        diffscore.annotate_peaks()
    
    info("Calling differentially occupied peaks")
    if options.score_method == 'p':
        diffscore.call_diff_peaks(cutoff = options.log_pvalue,
                                  min_length = options.dminlen,
                                  score_method = options.score_method)
    if options.score_method == 'q':
        diffscore.call_diff_peaks(cutoff = options.log_qvalue,
                                  min_length = options.dminlen,
                                  score_method = options.score_method)
#    diffscore.print_some_peaks()
#    diffscore.print_diff_peaks()
    
    info("Write output xls and BED files...")
    ofhd_xls = open( os.path.join( options.outdir, options.peakxls), "w" )
    ofhd_xls.write("# This file is generated by MACS version, using the diffpeak module %s\n" % (MACS_VERSION))
    ofhd_xls.write( options.argtxt+"\n" )
    ofhd_bed = open( os.path.join( options.outdir, options.peakbed), "w" )

    # pass write method so we can print too, and include name
    diffscore.write_peaks(xls=ofhd_xls, bed=ofhd_bed,
                    name = options.name, name_prefix="%s_peak_",
                    description="Peaks for %s (Made with MACS v2, " + strftime("%x") + ")",
                    trackline=options.trackline)
    ofhd_xls.close()
    ofhd_bed.close()
    
    if diffscore.has_peakio():
        info("Write annotated peak xls files...")
        ofhd_xls1 = open( os.path.join( options.outdir, options.peak1xls), "w" )
        ofhd_xls1.write("# This file is generated by MACS version, using the diffpeak module %s\n" % (MACS_VERSION))
        ofhd_xls1.write(options.argtxt+"\n")
        ofhd_xls2 = open( os.path.join( options.outdir, options.peak2xls), "w" )
        ofhd_xls2.write("# This file is generated by MACS version, using the diffpeak module %s\n" % (MACS_VERSION))
        ofhd_xls2.write(options.argtxt+"\n")
        diffscore.write_peaks_by_summit(ofhd_xls1, ofhd_xls2,
                                        name = options.name, name_prefix="%s_peak_")
        ofhd_xls1.close()
        ofhd_xls2.close()
    
    if options.store_bdg:
        info("#4 Write output bedgraph files...")
        ofhd_logLR = open( os.path.join( options.outdir, options.bdglogLR), "w" )
        ofhd_pvalue = open( os.path.join( options.outdir, options.bdgpvalue), "w" )
        ofhd_logFC = open( os.path.join( options.outdir, options.bdglogFC), "w" )
        diffscore.write_bedgraphs(logLR=ofhd_logLR, pvalue=ofhd_pvalue,
                                  logFC=ofhd_logFC, name = options.name,
                                  description=" for %s (Made with MACS v2, " + strftime("%x") + ")",
                                  trackline=options.trackline)
        ofhd_logLR.close()
        ofhd_pvalue.close()
        ofhd_logFC.close()
        
    
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

    tp = options.parser(options.tfile[0])
    treat = tp.build_petrack()
    treat.sort()
    if len(options.tfile) > 1:
        # multiple input
        for tfile in options.tfile[1:]:
            tp = options.parser(tfile)
            treat = tp.append_petrack( treat )
            treat.sort()

    options.tsize = tp.d
    if options.cfile:
        options.info("#1.2 read input fragments...")
        cp = options.parser(options.cfile[0])
        control = cp.build_petrack()
        control_d = cp.d
        control.sort()
        if len(options.cfile) > 1:
            # multiple input
            for cfile in options.cfile[1:]:
                cp = options.parser(cfile)
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
    tp = options.parser(options.tfile[0])
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize
    treat = tp.build_fwtrack()
    treat.sort()
    if len(options.tfile) > 1:
        # multiple input
        for tfile in options.tfile[1:]:
            tp = options.parser(tfile)
            treat = tp.append_fwtrack( treat )
            treat.sort()
    
    if options.cfile:
        options.info("#1.2 read input tags...")
        control = options.parser(options.cfile[0]).build_fwtrack()
        control.sort()
        if len(options.cfile) > 1:
            # multiple input
            for cfile in options.cfile[1:]:
                cp = options.parser(cfile)
                control = cp.append_fwtrack( control )
                control.sort()
    else:
        control = None
    options.info("#1 tag size is determined as %d bps" % options.tsize)
    return (treat, control)
