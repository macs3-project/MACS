# Time-stamp: <2022-02-23 01:55:55 Tao Liu>

"""Description: Main HMMR command

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see thefile LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
#from typing import Sized

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_hmmratac
from MACS3.IO.PeakIO import PeakIO
from MACS3.IO.Parser import BAMPEParser #BAMaccessor
from MACS3.Signal import HMMR_EM

#from MACS3.IO.BED import BEDreader # this hasn't been implemented yet.

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """The HMMRATAC function/pipeline for MACS.

    """
    #############################################
    # 1. Read the input BAM files
    #############################################
    options = opt_validate_hmmratac( args )
    options.info("\n" + options.argtxt)    
    options.info("# Read fragments from BAM file...")

    bam = BAMPEParser(options.bam_file[0], buffer_size=options.buffer_size)
    petrack = bam.build_petrack()
    if len( options.bam_file ) > 1:
        # multiple input
        for bamfile in options.bam_file[1:]:
            bam = BAMPEParser(bamfile, buffer_size=options.buffer_size)
            petrack = bam.append_petrack( petrack )
    petrack.finalize()

    # read in blacklisted if option entered    
    if options.blacklist:
        options.info("Read blacklist file...")
        peakio = open( options.blacklist )
        blacklist = PeakIO()
        i = 0
        for l in peakio:
            fs = l.rstrip().split()
            i += 1
            blacklist.add( fs[0].encode(), int(fs[1]), int(fs[2]), name=b"%d" % i )
            blacklist.sort()

    #############################################
    # 2. EM - created separate file for HMMR_EM
    #############################################
    # we will use EM to get the best means/stddevs for the mono-, di- and tri- modes of fragment sizes
    em_trainer = HMMR_EM.HMMR_EM( petrack, options.em_means[1:4], options.em_stddevs[1:4] ) # we take the options and initialize the object, then let it run
    # the mean and stddev after EM training
    em_means = [options.em_means[0],]
    em_means.extend(em_trainer.fragMeans)
    em_stddevs = [options.em_stddevs[0],]
    em_stddevs.extend(em_trainer.fragStddevs)    
    options.info( f"The mean and stddevs after EM:")
    options.info( f"Means: {em_means}")
    options.info( f"Stddevs: {em_stddevs}")


#############################################
# 3. Pileup
#############################################
    
    #pileup(SplitBed(GENOME, options.hmm_window), options.bamfile, options.index, options.misc_keep_duplicates)
    #bedGraphMath(pileup.out)
    #MergeBed(bedGraphMath.out, options.hmm_upper, options.hmm_lower) 
    #ExtendBed()  
    #SubtractedBed()
#############################################
# 4. Train HMM
#############################################
    #FragmentPileupGenerator(options.bamfile, options.index, options.training_set, options.em_means, options.em_stddev, options.min_map_quality, options.keep_duplicates)
    #KMeanstoHMM(FragmentPileupGenerator.out, options.hmm_states)
    #BaumWelch(KMeanstoHMM.out)
    #OpdfMultiGaussian(BaumWelch.out)
#############################################
# 5. Predict
#############################################
    #FragPileupGen(options.bamfile, options.index, tempBed, options.em_means, options.em_stddev, options.min_map_quality, options.keep_duplicatess, cpmScale)
    #HMMRTracksToBedgraph(FragPileupGen.out)
    #RobustHMM(FragPileupGen.out, BaumWelch.out)
    #PileupNode(RobustHMM.out)
#############################################
# 6. Output - add to OutputWriter
#############################################
    # bedGraphMath #for scoring
    #from MACS3.IO.OutputWriter import hmmratac_writer
    #hmmratac_writer()
    #
    # 
    #print ( options )
    #return
