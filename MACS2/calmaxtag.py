# Time-stamp: <2014-06-14 10:59:30 Tao Liu/JohnUrban>

"""Description: Calculate maximum duplicate (redundant) reads allowed depending on sequencing depth and genome size.

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
from MACS2.OptValidator import opt_validate_calmaxtag as opt_validate
from MACS2.cProb import binomial_cdf_inv
from MACS2.Constants import *
# ------------------------------------
# Main function
# ------------------------------------
def run( o_options ):
    """The calculation based on the binomial distribution for how many tags to allow at one site.
    
    """
    # Parse options...
    options = opt_validate( o_options )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    if options.outputfile != "stdout":
        outfhd = open( os.path.join( options.outdir, options.outputfile ) ,"w" )
    else:
        outfhd = sys.stdout
    
    #1 Read tag files
    if options.ifile:
	if not options.quiet:
            info("counting tags from input files...")
        fwtrack = load_tag_files_options (options)
        t0 = fwtrack.total
	if not options.quiet:
            info(" total tags in alignment file: %d" % (t0))
            info("tag size = %d" % options.tsize)
        fwtrack.fw = options.tsize
    elif options.numTags:
	t0 = options.numTags
    
    if not options.quiet:
    	info("calculate max duplicate tags in single position based on binomal distribution...")
    max_dup_tags = cal_max_dup_tags(options.gsize,t0)

    if not options.quiet:    
	info(" max_dup_tags based on binomal = %d" % (max_dup_tags))
    else:
	print max_dup_tags


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
    if not options.quiet: 
	options.info("read alignment tags...")
    tp = options.parser(options.ifile)

    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize

    treat = tp.build_fwtrack()
    treat.sort()

    if not options.quiet:
	options.info("tag size is determined as %d bps" % options.tsize)
    return treat

