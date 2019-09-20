# Time-stamp: <2019-09-20 11:30:27 taoliu>

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
import logging
from MACS2.IO import BedGraphIO
from MACS2.OptValidator import opt_validate_bdgopt as opt_validate

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
    options = opt_validate( options )
    info("Read and build bedGraph...")
    bio = BedGraphIO.bedGraphIO(options.ifile)
    btrack = bio.build_bdgtrack(baseline_value=0)

    info("Modify bedGraph...")
    if options.method.lower() == "p2q":
        btrack.p2q()
    elif options.method.lower() == "analen":
        btrack.analen()
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
    ofhd = open(ofile,"wb")
    btrack.write_bedGraph(ofhd,name="%s_modified_scores" % (options.method.upper()),description="Scores calculated by %s" % (options.method.upper()))
    info("Finished '%s'! Please check '%s'!" % (options.method, ofile))



