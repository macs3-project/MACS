# Time-stamp: <2012-09-14 02:27:53 Tao Liu>

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
    LLR_cutoff = options.cutoff
    ofile_prefix = options.oprefix

    info("Read and build treatment 1 bedGraph...")
    t1bio = cBedGraphIO.bedGraphIO(options.t1bdg)
    t1btrack = t1bio.build_bdgtrack()

    info("Read and build control 1 bedGraph...")
    c1bio = cBedGraphIO.bedGraphIO(options.c1bdg)
    c1btrack = c1bio.build_bdgtrack()

    depth1 = options.depth1

    info("Read and build treatment 2 bedGraph...")
    t2bio = cBedGraphIO.bedGraphIO(options.t2bdg)
    t2btrack = t2bio.build_bdgtrack()

    info("Read and build control 2 bedGraph...")
    c2bio = cBedGraphIO.bedGraphIO(options.c2bdg)
    c2btrack = c2bio.build_bdgtrack()

    depth2 = options.depth2

    twoconditionscore = cScoreTrack.TwoConditionScores( t1btrack,
                                                        c1btrack,
                                                        t2btrack,
                                                        c2btrack,
                                                        depth1,
                                                        depth2 )
    twoconditionscore.build()
    twoconditionscore.finalize()
    #(cat1,cat2,cat3,cat4) = twoconditionscore.call_peaks(min_length=options.minlen, cutoff=options.cutoff)
    (cat1,cat2,cat3) = twoconditionscore.call_peaks(min_length=options.minlen, cutoff=options.cutoff)

    info("Write peaks...")
    nf = open ("%s_c%.1f_cat1_peaks.encodePeak" % (options.oprefix,options.cutoff),"w")        
    cat1.write_to_narrowPeak(nf, name_prefix=options.oprefix+"_encodePeak", score_column="score")
    nf = open ("%s_c%.1f_cat2_peaks.encodePeak" % (options.oprefix,options.cutoff),"w")        
    cat2.write_to_narrowPeak(nf, name_prefix=options.oprefix+"_encodePeak", score_column="score")
    nf = open ("%s_c%.1f_cat3_peaks.encodePeak" % (options.oprefix,options.cutoff),"w")        
    cat3.write_to_narrowPeak(nf, name_prefix=options.oprefix+"_encodePeak", score_column="score")
    #nf = open ("%s_c%.1f_cat4_peaks.encodePeak" % (options.oprefix,options.cutoff),"w")        
    #cat4.write_to_narrowPeak(nf, name_prefix=options.oprefix+"_encodePeak", score_column="score")    
    info("Done")

    

