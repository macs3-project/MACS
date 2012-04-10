# Time-stamp: <2012-04-10 17:28:57 Tao Liu>

"""Description: Random sample certain number/percentage of tags.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

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

# ------------------------------------
# own python modules
# ------------------------------------
from MACS2.OptValidator import opt_validate_randsample as opt_validate
from MACS2.Constants import *

# ------------------------------------
# Main function
# ------------------------------------
def run( options0 ):
    options = opt_validate( options0 )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    #0 check output file
    if options.outputfile:
        assert not os.path.exists(options.outputfile), "%s already exists, please check!" % options.outputfile
        outfhd = open(options.outputfile,"w")
    else:
        outfhd = sys.stdout
    
    #1 Read tag files
    info("read tag files...")
    fwtrack = load_tag_files_options (options)
    
    info("tag size = %d" % options.tsize)
    fwtrack.fw = options.tsize

    t0 = fwtrack.total
    info(" total tags in alignment file: %d" % (t0))
    if options.number:
        if options.number > t0:
            error(" Number you want is bigger than total number of tags in alignment file! Please specify a smaller number and try again!")
            error(" %.2e > %.2e" % (options.number, t0))
            sys.exit(1)
        info(" Number of tags you want to keep: %.2e" % (options.number))
        options.percentage = float(options.number)/t0*100
    info(" Percentage of tags you want to keep: %.2f%%" % (options.percentage))
    
    fwtrack.sample_percent(options.percentage/100.0)

    info(" tags after random sampling in alignment file: %d" % (fwtrack.total))

    info("Write to BED file")
    fwtrack.print_to_bed(fhd=outfhd)
    info("finished! Check %s." % options.outputfile)

def load_tag_files_options ( options ):
    """From the options, load alignment tags.

    """
    options.info("read alignment tags...")
    tp = options.parser(options.tfile)

    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize

    treat = tp.build_fwtrack()
    treat.sort()

    options.info("tag size is determined as %d bps" % options.tsize)
    return treat

