# Time-stamp: <2013-07-31 15:23:17 Tao Liu>

"""Description: Naive call differential peaks from 4 bedGraph tracks for scores.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import logging
from MACS2.IO import cBedGraphIO
from MACS2.IO import cScoreTrack

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
    t1bio = cBedGraphIO.bedGraphIO(options.t1bdg)
    t1btrack = t1bio.build_bdgtrack()

    info("Read and build control 1 bedGraph...")
    c1bio = cBedGraphIO.bedGraphIO(options.c1bdg)
    c1btrack = c1bio.build_bdgtrack()

    info("Read and build treatment 2 bedGraph...")
    t2bio = cBedGraphIO.bedGraphIO(options.t2bdg)
    t2btrack = t2bio.build_bdgtrack()

    info("Read and build control 2 bedGraph...")
    c2bio = cBedGraphIO.bedGraphIO(options.c2bdg)
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

    twoconditionscore = cScoreTrack.TwoConditionScores( t1btrack,
                                                        c1btrack,
                                                        t2btrack,
                                                        c2btrack,
                                                        depth1,
                                                        depth2 )
    twoconditionscore.build()
    twoconditionscore.finalize()
    (cat1,cat2,cat3) = twoconditionscore.call_peaks(min_length=options.minlen, max_gap=options.maxgap, cutoff=options.cutoff)

    info("Write peaks...")
    nf = open ("%s_c%.1f_cond1.bed" % (options.oprefix,options.cutoff),"w")        
    cat1.write_to_bed(nf, name_prefix=options.oprefix+"_cond1_", name="condition 1", description="unique regions in condition 1", score_column="score")
    nf = open ("%s_c%.1f_cond2.bed" % (options.oprefix,options.cutoff),"w")        
    cat2.write_to_bed(nf, name_prefix=options.oprefix+"_cond2_", name="condition 2", description="unique regions in condition 2", score_column="score")
    nf = open ("%s_c%.1f_common.bed" % (options.oprefix,options.cutoff),"w")        
    cat3.write_to_bed(nf, name_prefix=options.oprefix+"_common_",name="common", description="common regions in both conditions", score_column="score")
    info("Done")

    

