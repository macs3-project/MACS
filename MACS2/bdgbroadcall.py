# Time-stamp: <2012-04-10 17:33:16 Tao Liu>

"""Description: Fine-tuning script to call broad peaks from a single bedGraph track for scores.

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
    info("Read and build bedGraph...")
    bio = cBedGraphIO.bedGraphIO(options.ifile)
    btrack = bio.build_bdgtrack(baseline_value=0)

    info("Call peaks from bedGraph...")
    #(peaks,bpeaks) = btrack.call_broadpeaks (lvl1_cutoff=options.cutoffpeak, lvl2_cutoff=options.cutofflink, min_length=options.minlen, lvl1_max_gap=options.lvl1maxgap, lvl2_max_gap=options.lvl2maxgap)
    bpeaks = btrack.call_broadpeaks (lvl1_cutoff=options.cutoffpeak, lvl2_cutoff=options.cutofflink, min_length=options.minlen, lvl1_max_gap=options.lvl1maxgap, lvl2_max_gap=options.lvl2maxgap)

    info("Write peaks...")
    #nf = open ("%s_c%.1f_l%d_g%d_peaks.encodePeak" % (options.oprefix,options.cutoffpeak,options.minlen,options.lvl1maxgap),"w")
    bf = open ("%s_c%.1f_C%.2f_l%d_g%d_G%d_broad.bed" % (options.oprefix,options.cutoffpeak,options.cutofflink,options.minlen,options.lvl1maxgap,options.lvl2maxgap),"w")        
    #peaks.write_to_narrowPeak(nf, name_prefix=options.oprefix+"_encodePeak", score_column="score")
    bpeaks[1].write_to_gappedPeak(bf, name_prefix=options.oprefix+"_broadRegion")    
    info("Done")
