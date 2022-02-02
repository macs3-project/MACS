# Time-stamp: <2022-02-02 13:50:18 Tao Liu>

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
from time import strftime
#from typing import Sized

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
# from MACS3.Utilities.OptValidator import opt_validate_hmmratac
from MACS3.IO.BAM import BAMaccessor

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
# 1. Parse Options
#############################################
    from MACS3.Utilities.OptValidator import opt_validate_hmmratac
    # how do we access the functions/classes that are in .c and .pyx formats
    from MACS3.IO import BAM, Parser
    options = opt_validate_hmmratac( args )

    bamfile = BAMaccessor( options.bamfile )
    options.info("Read BAM file together with the index BAI file ...")

    #what should replace this?  replace GenomeFileReader(genomeFile)
    # with index information genome file is the file containing
    # chromosome lengths, we get it from .get_rlengths (reference
    # chromosome lengths) function
    genome = bamfile.get_rlengths()
    
    # read in blacklisted if option entered    
    if options.blacklist:
        blacklist = BEDreader( options.blacklist )
        options.info("Read in blacklisted...")

#############################################
# 2. EM - created separate file for HMMR_EM
#############################################
    from MACS3.Signal import HMMR_EM # we have to implement this,
                                     # HMMR_EM would be a python class

    # these functions should be in their own HMMR_EM file
    # we want to access just the returnable objects from these funcions.

    hmmr_em_trainer = HMMR_EM( options, genome ) # we take the options and initialize the object, then let it run
    
    #HMMR_EM.pullLargeLengths(options.bamfile, options.min_map_quality, genome, options.em_means)

    #HMMR_EM.HMMR_EM(pullLargeLengths.out, options.em_means, options.em_stddev)

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
