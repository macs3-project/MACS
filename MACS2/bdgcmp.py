# Time-stamp: <2013-07-10 14:28:07 Tao Liu>

import sys
import logging
import io

from MACS2.IO import cBedGraphIO
from MACS2.IO.cBedGraph import scoreTracktoBedGraph
from MACS2.cProb import poisson_cdf
from math import log as mlog

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

LOG10_E = 0.43429448190325176

# ------------------------------------
# Misc functions
# ------------------------------------
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info
# ------------------------------------
# Main function
# ------------------------------------

def run( options ):
    scaling_factor = options.sfactor
    pseudo_depth = 1.0/scaling_factor   # not an actual depth, but its reciprocal, a trick to override SPMR while necessary.

    info("Read and build treatment bedGraph...")
    tbio = cBedGraphIO.bedGraphIO(options.tfile)
    tbtrack = tbio.build_bdgtrack()

    info("Read and build control bedGraph...")
    cbio = cBedGraphIO.bedGraphIO(options.cfile)
    cbtrack = cbio.build_bdgtrack()

    info("Build scoreTrackII...")
    sbtrack = tbtrack.make_scoreTrackII_for_macs( cbtrack, depth1 = pseudo_depth, depth2 = pseudo_depth )
    if abs(scaling_factor-1) > 1e-6:
        # Only for the case while your input is SPMR from MACS2 callpeak; Let's override SPMR.
        info("Values in your input bedGraph files will be multiplied by %f ..." % scaling_factor)
        sbtrack.change_normalization_method( ord('M') ) # a hack to override SPMR
    sbtrack.set_pseudocount( options.pseudocount )
    
    for method in set(options.method):
        info("Calculate scores comparing treatment and control by '%s'..." % method)
        ofile = options.ofile + "_" + method + ".bdg"
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
        else:
            raise Exception("Can't reach here!")
        
        info("Write bedGraph of scores...")
        ofhd = io.open(ofile,"wb")
        sbtrack.write_bedGraph(ofhd,name="%s_Scores" % (method.upper()),description="Scores calculated by %s" % (method.upper()), column = 3)
        info("Finished '%s'! Please check '%s'!" % (method, ofile))
