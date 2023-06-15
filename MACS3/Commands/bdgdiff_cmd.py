# Time-stamp: <2020-11-26 17:04:26 Tao Liu>

"""Description: Naive call differential peaks from 4 bedGraph tracks for scores.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import os
from MACS3.IO import BedGraphIO
from MACS3.Signal import ScoreTrack

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
import logging
import MACS3.Utilities.Logger

logger = logging.getLogger(__name__)
debug   = logger.debug
info    = logger.info
error   = logger.critical
warn    = logger.warning
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( options ):
    if options.maxgap >= options.minlen:
        error("MAXGAP should be smaller than MINLEN! Your input is MAXGAP = %d and MINLEN = %d" % (options.maxgap, options.minlen))

    LLR_cutoff = options.cutoff
    ofile_prefix = options.oprefix

    info("Read and build treatment 1 bedGraph...")
    t1bio = BedGraphIO.bedGraphIO(options.t1bdg)
    t1btrack = t1bio.build_bdgtrack()

    info("Read and build control 1 bedGraph...")
    c1bio = BedGraphIO.bedGraphIO(options.c1bdg)
    c1btrack = c1bio.build_bdgtrack()

    info("Read and build treatment 2 bedGraph...")
    t2bio = BedGraphIO.bedGraphIO(options.t2bdg)
    t2btrack = t2bio.build_bdgtrack()

    info("Read and build control 2 bedGraph...")
    c2bio = BedGraphIO.bedGraphIO(options.c2bdg)
    c2btrack = c2bio.build_bdgtrack()

    depth1 = options.depth1
    depth2 = options.depth2

    if depth1 > depth2:         # scale down condition 1 to size of condition 2
        depth1 = depth2 / depth1
        depth2 = 1.0
    elif depth1 < depth2:       # scale down condition 2 to size of condition 1
        depth2 = depth1/ depth2
        depth1 = 1.0
    else:                       # no need to scale down any
        depth1 = 1.0
        depth2 = 1.0

    twoconditionscore = ScoreTrack.TwoConditionScores( t1btrack,
                                                        c1btrack,
                                                        t2btrack,
                                                        c2btrack,
                                                        depth1,
                                                        depth2 )
    twoconditionscore.build()
    twoconditionscore.finalize()
    (cat1,cat2,cat3) = twoconditionscore.call_peaks(min_length=options.minlen, max_gap=options.maxgap, cutoff=options.cutoff)

    info("Write peaks...")

    if options.ofile:
        ofiles = [os.path.join( options.outdir, x ) for x in options.ofile]
        name_prefix = [ x.encode() for x in options.ofile ]
    else:
        ofiles = [ os.path.join( options.outdir, "%s_c%.1f_cond1.bed" % (options.oprefix,options.cutoff)),
                   os.path.join( options.outdir, "%s_c%.1f_cond2.bed" % (options.oprefix,options.cutoff)),
                   os.path.join( options.outdir, "%s_c%.1f_common.bed" % (options.oprefix,options.cutoff))
                   ]
        name_prefix = [ x.encode() for x in [ options.oprefix+"_cond1_", options.oprefix+"_cond2_", options.oprefix+"_common_" ]]

    nf = open( ofiles[ 0 ], 'w' )
    cat1.write_to_bed(nf, name_prefix=name_prefix[ 0 ], name=b"condition 1", description=b"unique regions in condition 1", score_column="score")
    nf.close()

    nf = open( ofiles[ 1 ], 'w' )
    cat2.write_to_bed(nf, name_prefix=name_prefix[ 1 ], name=b"condition 2", description=b"unique regions in condition 2", score_column="score")
    nf.close()

    nf = open( ofiles[ 2 ], 'w' )
    cat3.write_to_bed(nf, name_prefix=name_prefix[ 2 ], name=b"common", description=b"common regions in both conditions", score_column="score")
    nf.close()
    info("Done")

