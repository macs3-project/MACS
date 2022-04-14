# Time-stamp: <2022-04-14 16:00:18 Tao Liu>

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
from MACS3.Signal.HMMR_Signal_Processing import generate_weight_mapping, generate_digested_signals, extract_signals_from_training_regions
from MACS3.Signal.HMMR_HMM import hmm_training, hmm_predict


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
    fc_bdg = petrack.pileup_bdg( [1.0,], baseline_value = 0 )
    (sum_v, n_v, max_v, min_v, mean_v, std_v) = fc_bdg.summary()
    options.info( f"# Call peak above within FC range of {options.hmm_lower} and {options.hmm_upper}" )
    options.info( f"# Convert pileup to fold-change over average signal" )
    fc_bdg.apply_func(lambda x: x/mean_v)
    peaks = fc_bdg.call_peaks (cutoff=options.hmm_lower, up_limit=options.hmm_upper, min_length=petrack.average_template_length, max_gap=petrack.average_template_length, call_summits=False)
    options.info( f"# Total peaks within range: {peaks.total}" )

    # we will extend the training regions to both side of 500bps to include background regions ...
    
    # remove peaks overlapping with blacklisted regions
    if options.blacklist:
        peaks.exclude( blacklist )
    options.info( f"# after removing those overlapping with provided blacklisted regions, we have {peaks.total} left" )
    if peaks.total > options.hmm_maxTrain:
        peaks = peaks.randomly_pick( options.hmm_maxTrain, seed = options.hmm_randomSeed )
        options.info( f"# We randomly pick {options.hmm_maxTrain} peaks for training" )

    fhd = open(options.oprefix+"_training_peaks.bed","w")
    peaks.write_to_bed( fhd )
    fhd.close()

#############################################
# 4. Train HMM
#############################################
    options.info( f"# Compute the weights for each fragment length for each of the four signal types")
    fl_dict = petrack.count_fraglengths()
    fl_list = list(fl_dict.keys())
    fl_list.sort()
    # now we will prepare the weights for each fragment length for
    # each of the four distributions based on the EM results
    weight_mapping = generate_weight_mapping( fl_list, em_means, em_stddevs )
    
    options.info( f"# Generate short, mono-, di-, and tri-nucleosomal signals")
    digested_atac_signals = generate_digested_signals( petrack, weight_mapping )
    
    options.info( f"# Saving short, mono-, di-, and tri-nucleosomal signals to bedGraph files")
    
    fhd = open(options.oprefix+"_short.bdg","w")
    digested_atac_signals[ 0 ].write_bedGraph(fhd, "short","short")
    fhd.close()

    fhd = open(options.oprefix+"_mono.bdg","w")
    digested_atac_signals[ 1 ].write_bedGraph(fhd, "mono","mono")
    fhd.close()
    
    fhd = open(options.oprefix+"_di.bdg","w")
    digested_atac_signals[ 2 ].write_bedGraph(fhd, "di","di")
    fhd.close()
    
    fhd = open(options.oprefix+"_tri.bdg","w")
    digested_atac_signals[ 3 ].write_bedGraph(fhd, "tri","tri")
    fhd.close()

    # We first bin the training regions then get four types of signals
    # in the bins, at the same time, we record how many bins for each
    # peak.
    options.info( f"# Extract signals in training regions")
    [ training_data, training_data_lengths ] = extract_signals_from_training_regions( digested_atac_signals, peaks, binsize = 10, flanking = options.hmm_training_flanking )

    f = open(options.oprefix+"_training_data.txt","w")
    for v in training_data:
        f.write( f"{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\n" )
    f.close()
    
    f = open(options.oprefix+"_training_lens.txt","w")
    for v in training_data_lengths:
        f.write( f"{v}\n" )
    f.close()    
    
    options.info( f"# Use Baum-Welch algorithm to train the HMM")
    hmm_model = hmm_training( training_data, training_data_lengths )
    f = open(options.oprefix+"_model.txt","w")
    f.write( str(hmm_model.startprob_)+"\n" )
    f.write( str(hmm_model.transmat_ )+"\n" )
    f.write( str(hmm_model.means_ )+"\n" )
    f.write( str(hmm_model.covars_ )+"\n" )
    f.close()

#############################################
# 5. Predict
#############################################
    # Our prediction strategy will be different with HMMRATAC, we will first ask MACS call peaks with loose cutoff, then for each peak we will run HMM prediction to figure out labels. And for the rest of genomic regions, just mark them as 'background'.
    candidate_peaks = fc_bdg.call_peaks (cutoff=options.hmm_lower/2, up_limit=options.hmm_upper*2, min_length=petrack.average_template_length, max_gap=petrack.average_template_length, call_summits=False)
    options.info( f"# Total candidate peaks : {candidate_peaks.total}" )
    # remove peaks overlapping with blacklisted regions
    if options.blacklist:
        candidate_peaks.exclude( blacklist )
    options.info( f"# after removing those overlapping with provided blacklisted regions, we have {candidate_peaks.total} left" )

    # extract signals
    options.info( f"# Extract signals in candidate regions")
    [ candidate_data, candidate_data_lengths ] = extract_signals_from_training_regions( digested_atac_signals, candidate_peaks, binsize = 10 )
    
    options.info( f"# Use HMM to predict states")
    #predicted_states = hmm_predict( digested_atac_signals, hmm_model, binsize = 10 )
    predicted_states = hmm_predict( candidate_data, candidate_data_lengths, hmm_model, binsize = 10 )
    f = open(options.oprefix+"_predicted.txt","w")
    for l in range(len(predicted_states)):
        f.write ( f"{candidate_data[l]} {predicted_states[l]}\n" )
    f.close()

#############################################
# 6. Output - add to OutputWriter
#############################################
    options.info( f"# Write the output...")
    #predicted_states.write_to_bdg( file="" )

