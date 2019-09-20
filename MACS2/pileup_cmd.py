# Time-stamp: <2019-09-20 11:36:03 taoliu>

"""Description: Pileup alignment files

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
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
from MACS2.OptValidator import opt_validate_pileup as opt_validate
from MACS2.Pileup import pileup_and_write, pileup_and_write_pe
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
    options.PE_MODE = options.format in ('BAMPE','BEDPE')
    #assert options.format != 'BAMPE', "Pair-end data with BAMPE option currently doesn't work with pileup command. You can pretend your data to be single-end with -f BAM. Please try again!"


    #0 prepare output file
    outfile = os.path.join( options.outdir, options.outputfile )
    if os.path.isfile( outfile ):
        info("# Existing file %s will be replaced!" % outfile )
        os.unlink( outfile )
        
    #1 Read tag files
    info("# read alignment files...")
    if options.PE_MODE:
        info("# read input file in Paired-end mode.")
        treat = load_frag_files_options ( options ) # return PETrackI object
        t0 = treat.total # total fragments
        info("# total fragments/pairs in alignment file: %d" % (t0) )
        info("# Pileup paired-end alignment file.")
        pileup_and_write_pe(treat, outfile )
        
    else:
        (tsize, treat) = load_tag_files_options  (options)
    
        info("# tag size = %d", tsize)
    
        t0 = treat.total
        info("# total tags in alignment file: %d", t0)

        if options.bothdirection:
            info("# Pileup alignment file, extend each read towards up/downstream direction with %d bps" % options.extsize)
            pileup_and_write(treat, outfile, options.extsize * 2, 1, directional=False, halfextension=False)
        else:
            info("# Pileup alignment file, extend each read towards downstream direction with %d bps" % options.extsize)
            pileup_and_write(treat, outfile, options.extsize, 1, directional=True, halfextension=False)

    info("# Done! Check %s" % options.outputfile)

def load_tag_files_options ( options ):
    """From the options, load alignment tags.

    """
    options.info("# read treatment tags...")
    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    tsize = tp.tsize()
    treat = tp.build_fwtrack()
    #treat.sort()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )
            #treat.sort()
    treat.finalize()

    options.info("tag size is determined as %d bps" % tsize)
    return (tsize, treat)

def load_frag_files_options ( options ):
    """From the options, load treatment fragments and control fragments (if available).

    """
    options.info("# read treatment fragments...")

    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    treat = tp.build_petrack()
    #treat.sort()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_petrack( treat )
            #treat.sort()
    treat.finalize()
    return treat
