# Time-stamp: <2019-09-20 11:29:19 taoliu>

"""Description: Fine-tuning script to call broad peaks from a single
bedGraph track for scores.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).

"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import os
import logging
from MACS2.IO import BedGraphIO
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
    bio = BedGraphIO.bedGraphIO(options.ifile)
    btrack = bio.build_bdgtrack(baseline_value=0)

    info("Call peaks from bedGraph...")
    #(peaks,bpeaks) = btrack.call_broadpeaks (lvl1_cutoff=options.cutoffpeak, lvl2_cutoff=options.cutofflink, min_length=options.minlen, lvl1_max_gap=options.lvl1maxgap, lvl2_max_gap=options.lvl2maxgap)
    bpeaks = btrack.call_broadpeaks (lvl1_cutoff=options.cutoffpeak, lvl2_cutoff=options.cutofflink, min_length=options.minlen, lvl1_max_gap=options.lvl1maxgap, lvl2_max_gap=options.lvl2maxgap)

    info("Write peaks...")
    #nf = open ("%s_c%.1f_l%d_g%d_peaks.encodePeak" % (options.oprefix,options.cutoffpeak,options.minlen,options.lvl1maxgap),"w")
    if options.ofile:
        bf = open( os.path.join( options.outdir, options.ofile ), "w" )
        options.oprefix = options.ofile
    else:
        bf = open ( os.path.join( options.outdir, "%s_c%.1f_C%.2f_l%d_g%d_G%d_broad.bed12" % (options.oprefix,options.cutoffpeak,options.cutofflink,options.minlen,options.lvl1maxgap,options.lvl2maxgap)), "w" )
    bpeaks[1].write_to_gappedPeak(bf, name_prefix=options.oprefix+"_broadRegion")
    info("Done")
