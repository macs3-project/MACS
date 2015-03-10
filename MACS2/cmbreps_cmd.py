# Time-stamp: <2015-03-05 13:47:22 Tao Liu>

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

    info("combining #1 and #2 with method '%s'" % options.method)
    cmbtrack = reps[ 0 ].overlie( reps[ 1 ], func=options.method )
    ofile = os.path.join( options.outdir, options.ofile )
    info("Write bedGraph of combined scores...")
    ofhd = open(ofile,"wb")
    cmbtrack.write_bedGraph(ofhd,name="%s_combined_scores" % (options.method.upper()),description="Scores calculated by %s" % (options.method.upper()))
    info("Finished '%s'! Please check '%s'!" % (options.method, ofile))
    
