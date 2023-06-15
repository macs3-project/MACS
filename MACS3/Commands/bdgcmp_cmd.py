# Time-stamp: <2020-11-30 13:49:34 Tao Liu>

"""Description: compare bdg files

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import os

from MACS3.IO import BedGraphIO
from MACS3.Utilities.OptValidator import opt_validate_bdgcmp

from math import log as mlog

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
    options = opt_validate_bdgcmp( options )
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    
    scaling_factor = options.sfactor
    pseudo_depth = 1.0/scaling_factor   # not an actual depth, but its reciprocal, a trick to override SPMR while necessary.

    info("Read and build treatment bedGraph...")
    tbio = BedGraphIO.bedGraphIO(options.tfile)
    tbtrack = tbio.build_bdgtrack()

    info("Read and build control bedGraph...")
    cbio = BedGraphIO.bedGraphIO(options.cfile)
    cbtrack = cbio.build_bdgtrack()

    info("Build ScoreTrackII...")
    sbtrack = tbtrack.make_ScoreTrackII_for_macs( cbtrack, depth1 = pseudo_depth, depth2 = pseudo_depth )
    if abs(scaling_factor-1) > 1e-6:
        # Only for the case while your input is SPMR from MACS3 callpeak; Let's override SPMR.
        info("Values in your input bedGraph files will be multiplied by %f ..." % scaling_factor)
        sbtrack.change_normalization_method( ord('M') ) # a hack to override SPMR
    sbtrack.set_pseudocount( options.pseudocount )

    already_processed_method_list = []
    for (i, method) in enumerate(options.method):
        if method in already_processed_method_list:
            continue
        else:
            already_processed_method_list.append( method )

        info("Calculate scores comparing treatment and control by '%s'..." % method)
        if options.ofile:
            ofile = os.path.join( options.outdir, options.ofile[ i ] )
        else:
            ofile = os.path.join( options.outdir, options.oprefix + "_" + method + ".bdg" )
        # build score track
        if method == 'ppois':
            sbtrack.change_score_method( ord('p') )
        elif method == 'qpois':
            sbtrack.change_score_method( ord('q') )
        elif method == 'subtract':
            sbtrack.change_score_method( ord('d') )
        elif method == 'logFE':
            sbtrack.change_score_method( ord('f') )
        elif method == 'FE':
            sbtrack.change_score_method( ord('F') )
        elif method == 'logLR':             # log likelihood
            sbtrack.change_score_method( ord('l') )
        elif method == 'slogLR':             # log likelihood
            sbtrack.change_score_method( ord('s') )
        elif method == 'max':
            sbtrack.change_score_method( ord('M') )
        else:
            raise Exception("Can't reach here!")

        info("Write bedGraph of scores...")
        ofhd = open(ofile,"w")
        sbtrack.write_bedGraph(ofhd,name="%s_Scores" % (method.upper()),description="Scores calculated by %s" % (method.upper()), column = 3)
        info("Finished '%s'! Please check '%s'!" % (method, ofile))
