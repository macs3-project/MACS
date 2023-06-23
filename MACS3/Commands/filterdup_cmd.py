# Time-stamp: <2020-11-24 16:49:34 Tao Liu>

"""Description: Filter duplicate reads depending on sequencing depth.

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
from MACS3.Utilities.OptValidator import opt_validate_filterdup
from MACS3.Signal.Prob import binomial_cdf_inv
# ------------------------------------
# Main function
# ------------------------------------
def run( o_options ):
    """The Main function/pipeline for duplication filter.

    """
    # Parse options...
    options = opt_validate_filterdup( o_options )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    options.PE_MODE = options.format in ('BAMPE','BEDPE')

    if options.outputfile != "stdout":
        outfhd = open( os.path.join( options.outdir, options.outputfile ) ,"w" )
    else:
        outfhd = sys.stdout

    #1 Read tag files
    if options.PE_MODE:
        info("# read input file in Paired-end mode.")
        inputtrack = load_frag_files_options ( options ) # return PETrackI object
        t0 = inputtrack.total # total fragments
        info("# total fragments/pairs in alignment file: %d" % (t0) )
    else:
        info("# read tag files...")
        inputtrack = load_tag_files_options (options)

        info("# tag size = %d" % options.tsize)
        inputtrack.fw = options.tsize

        t0 = inputtrack.total
        info("# total tags in alignment file: %d" % (t0))

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

        inputtrack.filter_dup( max_dup_tags )
        t1 = inputtrack.total

        info(" tags after filtering in alignment file: %d" % (t1))
        info(" Redundant rate of alignment file: %.2f" % (float(t0-t1)/t0))

    if not options.dryrun:
        info( "Write to BED file" )
        inputtrack.print_to_bed( fhd=outfhd )
        info( "finished! Check %s." % options.outputfile )
    else:
        info( "Dry-run is finished!" )

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
    options.info("# read treatment tags...")
    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    if not options.tsize:           # override tsize if user specified --tsize
        ttsize = tp.tsize()
        options.tsize = ttsize

    treat = tp.build_fwtrack()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )
            #treat.sort()
    treat.finalize()
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
