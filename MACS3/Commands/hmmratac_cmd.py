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
from MACS3.Utilities.OptValidator import opt_validate_hmmratac

# These are just stand-ins... necessary functions need to be added/created still
from MACS3.HMMRATAC.ATACFragments import FragPileupGen
from MACS3.HMMRATAC.FormatConverters import PileupToBedGraph
from MACS3.HMMRATAC.GEMM import HMMR_EM
from MACS3.HMMRATAC.GenomeFileReaders import GenomeFileReader, bedFileReader
from MACS3.HMMRATAC.Node import PileupNode2, ScoreNode, TagNode
from MACS3.HMMRATAC.RobustHMM import KMeansToHMM, RobustHMM
from MACS3.HMMRATAC.WigMath import bedGraphMath, pileup


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
    options = opt_validate_hmmratac( args )
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    # input file
    bamfile = options.bam_file
    # output files
    storebdg = options.store_bdg
    store_bgscore = options.store_bgscore
    store_peaks = options.store_peaks
    # additional options
    print_ex = options.ptint_exclude
    print_train = options.print_train
    em_skip = options.em_skip
    em_means = options.em_means
    em_stddev = options.em_stddev
    states = options.hmm_states
    upper = options.hmm_upper
    lower = options.hmm_lower
    maxTrain = options.hmm_maxTrain
    train_regions = options.hmm_training_regions
    zscore = options.hmm_zscore
    randomSeed = options.hmm_randomSeed
    window = options.hmm_window
    model = options.hmm_file
    modelOnly = options.hmm_modelonly
    minlen = options.call_minlen
    score = options.call_score
    threshold = options.call_threshold
    blacklist = options.misc_blacklist
    keep_dups = options.misc_keep_duplicates
    trim = options.misc_trim
    multiple_processing = options.np
    verbose = options.verbose
    minMapQ = options.min_map_quality

# sections below are meant to be an outline of Main_HMMR_Driver.java
    
    # exit if input file not given
    if not bamfile:
        error( "BAM file must be included!")
        sys.exit( 1 )
    else:
        info("Read BAM file...")

    # print version number ... do we want to list a HMMRATAC version number?

    # print arguments entered

    # read in genome size stats ... is this necessary

    # read in blacklisted if option entered
    if blacklist:
        info("Read in blacklisted")
        #bedFileReader, getData are used in java file

    # set fragment length distribution parameters use input values to set initial values, if provided. else use defaults
    if em_means != [50, 200, 400, 600]:
        info("Input values used for initial mean values for the fragment distribution...")
    else:
        info("Defaults used for initial mean values for the fragment distribution...")
    if em_stddev != [20, 20, 20, 20]:
        info("Input values used for initial standard deviations for the fragment distribution...")
    else:
        info("Defaults used for initial standard deviations for the fragment distribution...")

    # Pull the lengths from the read data. Fragment distribution uses fragments with length > 100 to train the nucleosome distributions. 
    # the short distribution is set to 50 and remains unchanged. Only occurs if EM training occurs, else use the default settings
    if em_skip == False:
        # Perform EM training:
        # pullLargeLengths(bam, index, minMapQ, genomeStats, em_means), length = getSampledLengths, weights = getWeights used in java script
        # HMMR_EM(weights, em_means, em_stddev, length)
        info("Fragment Expectation Maximum performed on the fragment distribution...")
    else:
        info("EM training NOT performed on the fragment distribution...")
    

	# Generate genome wide pileup of read coverage and simultaneously calculate mean and standard deviation to calculate fold change and zscore. Convert pileup to bedgraph for easier storage. Calculate the pileup in chunks for memory efficiency
    # NOTE: this method does not currently exist. Instead we use an inputed big wig file to find fold change training regions and zscored masked regions
	    # various functions used here:
        # pileup, SplitBed(genomeStats, vitWindow, bam, index, rmdup), bedGraphMath(pileupData.getBedGraph)
        # pileupData.getCPMScale, genomeMean, genomeStd, MergeBed(upper, lower)
        # set train size, choose list with randomTrainSeed
    info("Training regions found and Zscore regions for exclusions found...")

    # Create the fragment pileup tracks using the training set and the fragment distribution parameters
    if model:
        info("Binary model file used (generated from previous HMMR run), training skipped... ")
    else:
        # print output to _training.bed
        # FragPileupGen, holder = TrackHolder
        info("Training Fragment Pileup completed...")
    
        # Create the inital model using KMeans and then refine it using Baum-Welch
        KMeansToHMM(holder)
        info("Kmeans Model: ...")

    # Use input model if available

    # Identify peak state as the state with the highest short read signal.
        # Identify flanking nucleosome state as the state witht the second highest mono read signal.

    # Output binary model file
    
    # Stop program if only model is desired  

    # Split the genome file into smaller 25MB chunks 
		#  * Can also split into whatever sized chunks the users prefers
		#  * May be necessary to split into smaller chunks for machines with less memory

    # Subtract excluded regions from the split genome for Viterbi

    # Run Viterbi on the whole genome 

    # Report the final results as peaks, bedgraphs and summits, if desired

        # Execute the scoring commands if the state is a peak or if bedgraph scoring is on

        # report the bedgraph if desired
        # report the peaks and summits if desired

    # print time to run


    print ( options )
    return
