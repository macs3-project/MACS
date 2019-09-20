# Time-stamp: <2019-09-20 11:34:13 taoliu>

"""Description: combine replicates

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import os
import logging

from MACS2.IO import BedGraphIO
from MACS2.OptValidator import opt_validate_cmbreps as opt_validate

from math import log as mlog

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

def run( options ):
    options = opt_validate( options )
    #weights = options.weights

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
    ofhd = open(ofile,"wb")
    cmbtrack.write_bedGraph(ofhd,name="%s_combined_scores" % (options.method.upper()),description="Scores calculated by %s" % (options.method.upper()))
    info("Finished '%s'! Please check '%s'!" % (options.method, ofile))
    
