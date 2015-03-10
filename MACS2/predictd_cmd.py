# Time-stamp: <2015-03-05 13:48:43 Tao Liu>

"""Description: Filter duplicate reads depending on sequencing depth.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
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

# ------------------------------------
# own python modules
# ------------------------------------
from MACS2.OptValidator import opt_validate_predictd as opt_validate
from MACS2.OutputWriter import *
from MACS2.PeakModel import PeakModel,NotEnoughPairsException
from MACS2.Prob import binomial_cdf_inv
from MACS2.Constants import *
# ------------------------------------
# Main function
# ------------------------------------
def run( o_options ):
    """The Main function/pipeline for duplication filter.
    
    """
    # Parse options...
    options = opt_validate( o_options )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    #0 output arguments
    assert options.format != 'BAMPE', "Pair-end data with BAMPE option doesn't work with predictd command. You can pretend your data to be single-end with -f BAM. Please try again!"
    
    #1 Read tag files
    info("# read alignment files...")
    treat = load_tag_files_options  (options)
    
    info("# tag size = %d", options.tsize)
    
    t0 = treat.total
    info("# total tags in alignment file: %d", t0)

    #2 Build Model
    info("# Build Peak Model...")

    try:
        peakmodel = PeakModel(treatment = treat,
                              max_pairnum = MAX_PAIRNUM,
                              opt = options
                              )
        info("# finished!")
        debug("#  Summary Model:")
        debug("#   min_tags: %d" % (peakmodel.min_tags))
        debug("#   d: %d" % (peakmodel.d))
        info("# predicted fragment length is %d bps" % peakmodel.d)
        info("# alternative fragment length(s) may be %s bps" % ','.join(map(str,peakmodel.alternative_d)))
        info("# Generate R script for model : %s" % (options.modelR))
        model2r_script(peakmodel,options.modelR, options.rfile )
        options.d = peakmodel.d

    except NotEnoughPairsException:
        warn("# Can't find enough pairs of symmetric peaks to build model!")

def load_tag_files_options ( options ):
    """From the options, load alignment tags.

    """
    options.info("# read treatment tags...")
    tp = options.parser(options.ifile[0])
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize
    treat = tp.build_fwtrack()
    treat.sort()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(tfile)
            treat = tp.append_fwtrack( treat )
            treat.sort()

    options.info("tag size is determined as %d bps" % options.tsize)
    return treat

