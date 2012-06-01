# Time-stamp: <2012-05-24 13:55:23 Tao Liu>

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
# pscore_dict = {}

# def get_pscore ( observed, expectation ):
#     key_value = (observed, expectation)
#     if pscore_dict.has_key(key_value):
#         return pscore_dict[key_value]
#     else:
#         score = -1*poisson_cdf(observed,expectation,False,True)
#         pscore_dict[(observed,expectation)] = score
#         return score

# logLR_dict = {}

# def logLR ( x, y ):
#     """Calculate log10 Likelihood between H1 ( enriched ) and H0 (
#     chromatin bias ). Then set minus sign for depletion.
    
#     """
#     key_value = ( x, y )
#     if logLR_dict.has_key(key_value):
#         return logLR_dict[key_value]
#     else:
#         if x > y:
#             s = (x*(mlog(x+1)-mlog(y+1))+y-x)*LOG10_E
#         elif x < y:
#             s = (-1*x*(mlog(x+1)-mlog(y+1))-y+x)*LOG10_E
#         else:
#             s = 0
#         logLR_dict[key_value] = s
#         return s

# def logFE ( x, y ):
#     return mlog( (x+1.0)/(y+1.0) )*LOG10_E

def run( options ):
    info("Read and build treatment bedGraph...")
    tbio = cBedGraphIO.bedGraphIO(options.tfile)
    tbtrack = tbio.build_bdgtrack()

    info("Read and build control bedGraph...")
    cbio = cBedGraphIO.bedGraphIO(options.cfile)
    cbtrack = cbio.build_bdgtrack()

    info("Build scoreTrackII...")
    sbtrack = tbtrack.make_scoreTrackII_for_macs( cbtrack )

    method = options.method

    info("Calculate scores comparing treatment and control by %s..." % method)
    # build score track
    if method == 'ppois':
        sbtrack.change_score_method( ord('p') )
    elif method == 'qpois':
        sbtrack.change_score_method( ord('q') )        
    elif method == 'substract':
        sbtrack.change_score_method( ord('d') )        
        sbtrack = tbtrack.overlie(cbtrack,func=lambda x,y:x-y)
    elif method == 'logFE':
        sbtrack.change_score_method( ord('f') )
    elif method == 'FE':
        sbtrack.change_score_method( ord('F') )        
    elif method == 'logLR':             # log likelihood
        sbtrack.change_score_method( ord('l') )
    else:
        raise Exception("Can't reach here!")

    info("Write to output bedGraph...")
    ofhd = io.open(options.ofile,"wb")

    #r = sbtrack.get_data_by_chr("chr22")

    #print r

    sbtrack.write_bedGraph(ofhd,name="%s_Scores" % (method.upper()),description="Scores calculated by %s" % (method.upper()), column = 3)
