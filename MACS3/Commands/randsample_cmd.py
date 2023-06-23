# Time-stamp: <2020-11-24 17:00:16 Tao Liu>

"""Description: Random sample certain number/percentage of tags.

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
from MACS3.Utilities.OptValidator import opt_validate_randsample

# ------------------------------------
# Main function
# ------------------------------------
def run( options0 ):
    options = opt_validate_randsample( options0 )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    options.PE_MODE = options.format in ('BAMPE','BEDPE')

    #0 check output file
    if options.outputfile:
        outfhd = open( os.path.join( options.outdir, options.outputfile ), "w" )
    else:
        outfhd = sys.stdout

    #1 Read tag files
    if options.PE_MODE:
        info("# read input file in Paired-end mode.")
        treat = load_frag_files_options ( options ) # return PETrackI object
        t0 = treat.total # total fragments
        info("# total fragments/pairs in alignment file: %d" % (t0) )
    else:
        info("read tag files...")
        treat = load_tag_files_options (options)

        info("tag size = %d" % options.tsize)
        treat.fw = options.tsize

        t0 = treat.total
        info(" total tags in alignment file: %d" % (t0))

    if options.number:
        if options.number > t0:
            error(" Number you want is bigger than total number of tags in alignment file! Please specify a smaller number and try again!")
            error(" %.2e > %.2e" % (options.number, t0))
            sys.exit(1)
        info(" Number of tags you want to keep: %.2e" % (options.number))
        options.percentage = float(options.number)/t0*100
    info(" Percentage of tags you want to keep: %.2f%%" % (options.percentage))

    if options.seed >= 0:
        info(" Random seed has been set as: %d" % options.seed )

    treat.sample_percent(options.percentage/100.0, options.seed )

    info(" tags after random sampling in alignment file: %d" % (treat.total))

    info("Write to BED file")
    treat.print_to_bed(fhd=outfhd)
    info("finished! Check %s." % options.outputfile)

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
        for ifile in options.ifile[1:]:
            tp = options.parser(ifile, buffer_size=options.buffer_size)
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
    #treat.sort()
    if len(options.ifile) > 1:
        # multiple input
        for ifile in options.ifile[1:]:
            tp = options.parser(ifile, buffer_size=options.buffer_size)
            treat = tp.append_petrack( treat )
            #treat.sort()
    treat.finalize()
    return treat
