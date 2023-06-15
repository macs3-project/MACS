# Time-stamp: <2020-11-24 16:59:33 Tao Liu>

"""Description: predict fragment size.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_predictd
from MACS3.Signal.PeakModel import PeakModel,NotEnoughPairsException
from MACS3.Signal.Prob import binomial_cdf_inv
from MACS3.IO.OutputWriter import model2r_script
# ------------------------------------
# Main function
# ------------------------------------
def run( o_options ):
    """The Main function/pipeline for duplication filter.

    """
    # Parse options...
    options = opt_validate_predictd( o_options )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    #0 output arguments
    options.PE_MODE = options.format in ('BAMPE','BEDPE')

    #1 Read tag files
    if options.PE_MODE:
        info("# read input file in Paired-end mode.")
        treat = load_frag_files_options ( options ) # return PETrackI object
        t0 = treat.total
        info("# total fragments/pairs in alignment file: %d" % (t0) )
    else:
        info("# read alignment files...")
        treat = load_tag_files_options  (options)
        t0 = treat.total
        info("# tag size = %d" % options.tsize)
        treat.fw = options.tsize
        info("# total tags in alignment file: %d", t0)

    #2 Build Model
    info("# Build Peak Model...")
    if options.PE_MODE:
        d = treat.average_template_length
        info("# Average insertion length of all pairs is %d bps" % d)
        return

    try:
        peakmodel = PeakModel(treatment = treat,
                              max_pairnum = MAX_PAIRNUM,
                              opt = options
                              )
        peakmodel.build()
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
    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize
    treat = tp.build_fwtrack()
    #treat.sort()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )
            #treat.sort()
    treat.finalize()

    options.info("tag size is determined as %d bps" % options.tsize)
    return treat

def load_frag_files_options ( options ):
    """From the options, load treatment fragments and control fragments (if available).

    """
    options.info("# read treatment fragments...")

    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    treat = tp.build_petrack()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(ifile, buffer_size=options.buffer_size)
            treat = tp.append_petrack( treat )
    treat.finalize()
    return treat
