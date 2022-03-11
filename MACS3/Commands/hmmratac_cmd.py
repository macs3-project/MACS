# Time-stamp: <2022-03-10 17:05:08 Tao Liu>

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
from MACS3.Signal.HMMR_EM import HMMR_EM
from MACS3.Signal.HMMR_Signal_Processing import generate_weight_mapping,  generate_digested_signals, extract_signals_from_training_regions
from MACS3.Signal.HMMR_HMM import initial_state_kmeans, hmm_training, hmm_predict


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
    # remember to finalize the petrack
    petrack.finalize()

    # filter duplicates if needed
    if options.misc_keep_duplicates:
        petrack.filter_dup( maxnum=1 )

    # read in blacklisted if option entered    
    if options.blacklist:
        options.info("# Read blacklist file...")
        peakio = open( options.blacklist )
        blacklist = PeakIO()
        i = 0
        for l in peakio:
            fs = l.rstrip().split()
            i += 1
            blacklist.add( fs[0].encode(), int(fs[1]), int(fs[2]), name=b"%d" % i )
            blacklist.sort()

    #############################################
    # 2. EM
    #############################################
    # we will use EM to get the best means/stddevs for the mono-, di- and tri- modes of fragment sizes
    options.info("# Use EM algorithm to estimate means and stddevs of fragment lengths")
    options.info("  for mono-, di-, and tri-nucleosomal signals...") 
    em_trainer = HMMR_EM( petrack, options.em_means[1:4], options.em_stddevs[1:4], seed = options.hmm_randomSeed )
    # the mean and stddev after EM training
    em_means = [options.em_means[0],]
    em_means.extend(em_trainer.fragMeans)
    em_stddevs = [options.em_stddevs[0],]
    em_stddevs.extend(em_trainer.fragStddevs)    
    options.info( f"# The mean and stddevs after EM:")
    options.info( f"         [short,  mono-,  di-,  tri-]")
    options.info( f"  Means: {em_means}")
    options.info( f"  Stddevs: {em_stddevs}")

    #############################################
    # 3. Define training set by peak calling
    #############################################

    # Find regions with fold change within determined range to use as training sites.
    # Find regions with zscore values above certain cutoff to exclude from viterbi.
    # 
    options.info( f"# Look for training set from {petrack.total} fragments" )
    options.info( f"# Pile up all fragments" )
    bdg = petrack.pileup_bdg( [1.0,], baseline_value = 0 )
    (sum_v, n_v, max_v, min_v, mean_v, std_v) = bdg.summary()
    options.info( f"# Call peak above within FC range of {options.hmm_lower} and {options.hmm_upper}" )
    options.info( f"# Convert pileup to fold-change over average signal" )
    bdg.apply_func(lambda x: x/mean_v)
    peaks = bdg.call_peaks (cutoff=options.hmm_lower, up_limit=options.hmm_upper, min_length=petrack.average_template_length, max_gap=petrack.average_template_length, call_summits=False)
    options.info( f"# Total peaks within range: {peaks.total}" )
    # remove peaks overlapping with blacklisted regions
    if options.blacklist:
        peaks.exclude( blacklist )
    options.info( f"# after removing those overlapping with provided blacklisted regions, we have {peaks.total} left" )

#############################################
# 4. Train HMM
#############################################
    options.info( f"# Compute the weights for each fragment length for each of the four signal types")
    fl_dict = petrack.count_fraglengths()
    fl_list = list(fl_dict.keys())
    fl_list.sort()
    # now we will prepare the weights for each fragment length for each of the four distributions based on the EM results
    weight_mapping = generate_weight_mapping( fl_list, em_means, em_stddevs )

    options.info( f"# Generate short, mono-, di-, and tri-nucleosomal signals")    
    digested_atac_signals = generate_digested_signals( petrack, weight_mapping )

    options.info( f"# Extract signals in training regions")
    training_data = extract_signals_from_training_regions( digested_atac_signals, peaks, binsize = 10 )
    #FragmentPileupGenerator(options.bamfile, options.index, options.training_set, options.em_means, options.em_stddev, options.min_map_quality, options.keep_duplicates)
    
    options.info( f"# Use K-means method to build initial states")
    initial_states = initial_state_kmeans( training_data, k=3 )
    #KMeanstoHMM(FragmentPileupGenerator.out, options.hmm_states)
    
    options.info( f"# Use Baum-Welch algorithm to train the HMM")
    hmm_model = hmm_training( training_data, initial_state )

#############################################
# 5. Predict
#############################################
    #FragPileupGen(options.bamfile, options.index, tempBed, options.em_means, options.em_stddev, options.min_map_quality, options.keep_duplicatess, cpmScale)
    #HMMRTracksToBedgraph(FragPileupGen.out)
    
    options.info( f"# Use HMM to predict states")
    predicted_states = hmm_predict( digested_atac_signals, hmm_model, binsize = 10 )

#############################################
# 6. Output - add to OutputWriter
#############################################
    options.info( f"# Write the output...")
    #predicted_states.write_to_bdg( file="" )

