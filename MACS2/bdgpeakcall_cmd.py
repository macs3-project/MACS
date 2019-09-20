# Time-stamp: <2019-09-20 11:30:45 taoliu>

"""Description: Naive call peaks from a single bedGraph track for
scores.

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

    if options.cutoff_analysis:
        info("Analyze cutoff vs number of peaks/total length of peaks/average length of peak")
        cutoff_analysis_result = btrack.cutoff_analysis( int(options.maxgap), int(options.minlen), 50 )
        info("Write report...")
        if options.ofile:
            fhd = open( os.path.join( options.outdir, options.ofile ), 'w' )
        else:
            fhd = open ( os.path.join( options.outdir, "%s_l%d_g%d_cutoff_analysis.txt" % (options.oprefix,options.minlen,options.maxgap)), "w" )
        fhd.write( cutoff_analysis_result )
        info("Done")
    else:
        info("Call peaks from bedGraph...")
        peaks = btrack.call_peaks(cutoff=float(options.cutoff),min_length=int(options.minlen),max_gap=int(options.maxgap),call_summits=options.call_summits)

        info("Write peaks...")
        if options.ofile:
            options.oprefix = options.ofile
            nf = open( os.path.join( options.outdir, options.ofile ), 'w' )
        else:
            nf = open ( os.path.join( options.outdir, "%s_c%.1f_l%d_g%d_peaks.narrowPeak" % (options.oprefix,options.cutoff,options.minlen,options.maxgap)), "w" )
        peaks.write_to_narrowPeak(nf, name=options.oprefix, name_prefix=options.oprefix+"_narrowPeak", score_column="score", trackline=options.trackline)
        info("Done")



