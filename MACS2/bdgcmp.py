# Time-stamp: <2012-06-06 14:55:12 Tao Liu>

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
    info("Read and build treatment bedGraph...")
    tbio = cBedGraphIO.bedGraphIO(options.tfile)
    tbtrack = tbio.build_bdgtrack()

    treat_depth = options.tdepth

    info("Read and build control bedGraph...")
    cbio = cBedGraphIO.bedGraphIO(options.cfile)
    cbtrack = cbio.build_bdgtrack()

    ctrl_depth = options.cdepth

    info("Build scoreTrackII...")
    sbtrack = tbtrack.make_scoreTrackII_for_macs( cbtrack, depth1 = treat_depth, depth2 = ctrl_depth )
    # normalize by depth
    if abs(treat_depth-1) > 1e-6 or abs(ctrl_depth-1) > 1e-6:
        # if depth of treat and control is 1.0 ( files are generated
        # by MACS2 --SPMR ), no need for the following step.
        info("Normalize by sequencing depth of million reads...")
        sbtrack.change_normalization_method( ord('M') )
    sbtrack.set_pseudocount( options.pseudocount )
    
    #def make_scoreTrackII_for_macs (self, bdgTrack2, float depth1 = 1.0, float depth2 = 1.0 ):
    
    method = options.method

    info("Calculate scores comparing treatment and control by %s..." % method)
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
    else:
        raise Exception("Can't reach here!")

    info("Write bedGraph of scores...")
    ofhd = io.open(options.ofile,"wb")

    #r = sbtrack.get_data_by_chr("chr22")

    #print r

    sbtrack.write_bedGraph(ofhd,name="%s_Scores" % (method.upper()),description="Scores calculated by %s" % (method.upper()), column = 3)
    info("Finished! Please check %s!" % (options.ofile))
