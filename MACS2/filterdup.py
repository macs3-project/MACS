# Time-stamp: <2012-10-03 12:18:11 Tao Liu>

"""Description: Filter duplicate reads depending on sequencing depth.

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
from MACS2.OptValidator import opt_validate_filterdup as opt_validate
from MACS2.cProb import binomial_cdf_inv
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

    if options.outputfile != "stdout":
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
    if options.keepduplicates != "all":
        if options.keepduplicates == "auto":
            info("calculate max duplicate tags in single position based on binomal distribution...")
            max_dup_tags = cal_max_dup_tags(options.gsize,t0)
            info(" max_dup_tags based on binomal = %d" % (max_dup_tags))
            info("filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)" % (max_dup_tags))
        else:
            info("user defined the maximum tags...")
            max_dup_tags = int(options.keepduplicates)
            info("filter out redundant tags at the same location and the same strand by allowing at most %d tag(s)" % (max_dup_tags))

        fwtrack = fwtrack.filter_dup(max_dup_tags)
        t1 = fwtrack.total
        info(" tags after filtering in alignment file: %d" % (t1))
        info(" Redundant rate of alignment file: %.2f" % (float(t0-t1)/t0))
    info("Write to BED file")
    fwtrack.print_to_bed(fhd=outfhd)
    info("finished! Check %s." % options.outputfile)

def cal_max_dup_tags ( genome_size, tags_number, p=1e-5 ):
    """Calculate the maximum duplicated tag number based on genome
    size, total tag number and a p-value based on binomial
    distribution. Brute force algorithm to calculate reverse CDF no
    more than MAX_LAMBDA(100000).
    
    """
    return binomial_cdf_inv(1-p,tags_number,1.0/genome_size)

def load_tag_files_options ( options ):
    """From the options, load alignment tags.

    """
    options.info("read alignment tags...")
    tp = options.parser(options.ifile)

    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize

    treat = tp.build_fwtrack()
    treat.sort()

    options.info("tag size is determined as %d bps" % options.tsize)
    return treat

