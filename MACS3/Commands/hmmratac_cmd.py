# Time-stamp: <2021-12-20 09:19:28 ta32852>

"""Description: Main HMMR command

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
from time import strftime
from typing import Sized

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
# from MACS3.Utilities.OptValidator import opt_validate_hmmratac

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
# 1. Parse Options
#############################################
    from MACS3.Utilities.OptValidator import opt_validate_hmmratac
    # how do we access the functions/classes that are in .c and .pyx formats
    from MACS3.IO import BAM, Parser
    options = opt_validate_hmmratac( args )


    BAM.BAMaccessor.readBamfile(options.bamfile)
    options.info("Read BAM file...")

    BAM.BAIfile.readBAIfile(options.index)
    options.info("Access index file...")

    #what should replace this?
    # replace GenomeFileReader(genomeFile) with index information
    GENOME = GenomeFileReader(genomeFile)

    # read in blacklisted if option entered    
    if options.blacklist:
        Parser.bedParser(options.blacklist)
        options.info("Read in blacklisted...")

#############################################
# 2. EM - created separate file for HMMR_EM
#############################################
    from MACS3.Signal import HMMR_EM

    # these functions should be in their own HMMR_EM file
    # we want to access just the returnable objects from these funcions.
    HMMR_EM(options)

    HMMR_EM.pullLargeLengths(options.bamfile, options.index, options.min_map_quality, GENOME, options.em_means)
    HMMR_EM.HMMR_EM(pullLargeLengths.out, options.em_means, options.em_stddev)

#############################################
# 3. Pileup
#############################################
    
    pileup(SplitBed(GENOME, options.hmm_window), options.bamfile, options.index, options.misc_keep_duplicates)
    bedGraphMath(pileup.out)
    MergeBed(bedGraphMath.out, options.hmm_upper, options.hmm_lower) 
    ExtendBed()  
    SubtractedBed()
#############################################
# 4. Train HMM
#############################################
    FragmentPileupGenerator(options.bamfile, options.index, options.training_set, options.em_means, options.em_stddev, options.min_map_quality, options.keep_duplicates)
    KMeanstoHMM(FragmentPileupGenerator.out, options.hmm_states)
    BaumWelch(KMeanstoHMM.out)
    OpdfMultiGaussian(BaumWelch.out)
#############################################
# 5. Predict
#############################################
    FragPileupGen(options.bamfile, options.index, tempBed, options.em_means, options.em_stddev, options.min_map_quality, options.keep_duplicatess, cpmScale)
    HMMRTracksToBedgraph(FragPileupGen.out)
    RobustHMM(FragPileupGen.out, BaumWelch.out)
    PileupNode(RobustHMM.out)
#############################################
# 6. Output - add to OutputWriter
#############################################
    # bedGraphMath #for scoring
    from MACS3.IO.OutputWriter import hmmratac_writer
    hmmratac_writer()
    
    
    print ( options )
    return
