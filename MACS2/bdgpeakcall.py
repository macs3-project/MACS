# Time-stamp: <2012-04-10 17:33:25 Tao Liu>

"""Description: Naive call peaks from a single bedGraph track for scores.

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
    peaks = btrack.call_peaks(cutoff=float(options.cutoff),min_length=int(options.minlen),max_gap=int(options.maxgap),call_summits=options.call_summits)

    info("Write peaks...")
    nf = open ("%s_c%.1f_l%d_g%d_peaks.encodePeak" % (options.oprefix,options.cutoff,options.minlen,options.maxgap),"w")        
    peaks.write_to_narrowPeak(nf, name_prefix=options.oprefix+"_encodePeak", score_column="score", trackline=options.trackline)
    info("Done")
    
