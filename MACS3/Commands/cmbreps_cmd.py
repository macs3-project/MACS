# Time-stamp: <2020-11-24 16:47:11 Tao Liu>

"""Description: combine replicates

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import os
from math import log as mlog

from MACS3.IO import BedGraphIO
from MACS3.Utilities.OptValidator import opt_validate_cmbreps

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------

def run( options ):
    options = opt_validate_cmbreps( options )
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    
    info("Read and build bedGraph for each replicate...")
    reps = []
    i = 1
    for ifile in options.ifile:
        info("Read file #%d" % i)
        reps.append( BedGraphIO.bedGraphIO( ifile ).build_bdgtrack( ) )
        i += 1

    # first two reps

    info("combining tracks 1-%i with method '%s'" % (i - 1, options.method))
    cmbtrack = reps[ 0 ].overlie( [reps[ j ] for j in range(1, i - 1)], func=options.method )
    ofile = os.path.join( options.outdir, options.ofile )
    info("Write bedGraph of combined scores...")
    ofhd = open(ofile,"w")
    cmbtrack.write_bedGraph(ofhd,name="%s_combined_scores" % (options.method.upper()),description="Scores calculated by %s" % (options.method.upper()))
    info("Finished '%s'! Please check '%s'!" % (options.method, ofile))

