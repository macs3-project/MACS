# Time-stamp: <2020-11-24 16:46:33 Tao Liu>

"""Description: Modify bedGraph file

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
from MACS3.Utilities.OptValidator import opt_validate_bdgopt

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( options ):
    options = opt_validate_bdgopt( options )
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    
    info("Read and build bedGraph...")
    bio = BedGraphIO.bedGraphIO(options.ifile)
    btrack = bio.build_bdgtrack(baseline_value=0)

    info("Modify bedGraph...")
    if options.method.lower() == "p2q":
        btrack.p2q()

    else:
        extraparam = float(options.extraparam[0])
        if options.method.lower() == "multiply":
            btrack.apply_func( lambda x: x * extraparam)
        elif options.method.lower() == "add":
            btrack.apply_func( lambda x: x + extraparam)
        elif options.method.lower() == "max":
            btrack.apply_func( lambda x: x if x> extraparam else extraparam )
        elif options.method.lower() == "min":
            btrack.apply_func( lambda x: x if x< extraparam else extraparam )

    ofile = os.path.join( options.outdir, options.ofile )
    info("Write bedGraph of modified scores...")
    ofhd = open(ofile,"w")
    btrack.write_bedGraph(ofhd,name="%s_modified_scores" % (options.method.upper()),description="Scores calculated by %s" % (options.method.upper()))
    info("Finished '%s'! Please check '%s'!" % (options.method, ofile))



