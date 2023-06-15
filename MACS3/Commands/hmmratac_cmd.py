# Time-stamp: <2023-06-08 11:03:46 Tao Liu>

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
import gc
import numpy as np
import json
from hmmlearn import hmm
#from typing import Sized

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_hmmratac
from MACS3.IO.PeakIO import PeakIO
from MACS3.IO.PeakIO import BroadPeakIO
from MACS3.IO.Parser import BAMPEParser #BAMaccessor
from MACS3.Signal.HMMR_EM import HMMR_EM
from MACS3.Signal.HMMR_Signal_Processing import generate_weight_mapping, generate_digested_signals, extract_signals_from_regions
from MACS3.Signal.HMMR_HMM import hmm_training, hmm_predict, hmm_model_init, hmm_model_save
from MACS3.Signal.Region import Regions
from MACS3.Signal.BedGraph import bedGraphTrackI

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
    options = opt_validate_hmmratac( args )

    #############################################
    # 0. Names of output files
    #############################################
    short_bdgfile = os.path.join( options.outdir, options.name+"_digested_short.bdg" )
    mono_bdgfile = os.path.join( options.outdir, options.name+"_digested_mono.bdg" )
    di_bdgfile = os.path.join( options.outdir, options.name+"_digested_di.bdg" )
    tri_bdgfile = os.path.join( options.outdir, options.name+"_digested_tri.bdg" )
    training_region_bedfile = os.path.join( options.outdir, options.name+"_training_regions.bed" )
    training_datafile = os.path.join( options.outdir, options.name+"_training_data.txt" )
    training_datalengthfile = os.path.join( options.outdir, options.name+"_training_lengths.txt" )
    hmm_modelfile = os.path.join( options.outdir, options.name+"_model.json" )
    open_state_bdgfile = os.path.join( options.outdir, options.name+"_open.bdg" )
    nuc_state_bdgfile = os.path.join( options.outdir, options.name+"_nuc.bdg" )
    bg_state_bdgfile = os.path.join( options.outdir, options.name+"_bg.bdg" )
    states_file = os.path.join( options.outdir, options.name+"_states.bed" )
    accessible_file = os.path.join( options.outdir, options.name+"_accessible_regions.gappedPeak" )
    cutoffanalysis_file = os.path.join( options.outdir, options.name+"_cutoff_analysis.tsv" )
    
    #############################################
    # 1. Read the input BAM files
    #############################################
    options.info("\n" + options.argtxt)
    #options.info("Some important parameters for this run")
    #options.info(f" binsize for HMM: {options.binsize} ()")
    options.info("#1 Read fragments from BAM file...")

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
        options.info("#  Read blacklist file...")
        peakio = open( options.blacklist )
        blacklist = PeakIO()
        i = 0
        for l in peakio:
            fs = l.rstrip().split()
            i += 1
            blacklist.add( fs[0].encode(), int(fs[1]), int(fs[2]), name=b"%d" % i )
            blacklist.sort()
        blacklist_regions = Regions()
        blacklist_regions.init_from_PeakIO( blacklist )

    #############################################
    # 2. EM
    #############################################
    if options.em_skip:
        # Skip EM and use the options.em_means and options.em_stddevs
        em_means = options.em_means
        em_stddevs = options.em_stddevs
        options.info( "#2 EM is skipped. The following means and stddevs will be used:" )
    else:
        # we will use EM to get the best means/stddevs for the mono-, di- and tri- modes of fragment sizes
        options.info("#2 Use EM algorithm to estimate means and stddevs of fragment lengths")
        options.info("#  for mono-, di-, and tri-nucleosomal signals...") 
        em_trainer = HMMR_EM( petrack, options.em_means[1:4], options.em_stddevs[1:4], seed = options.hmm_randomSeed )
        # the mean and stddev after EM training
        em_means = [options.em_means[0],]
        em_means.extend(em_trainer.fragMeans)
        em_stddevs = [options.em_stddevs[0],]
        em_stddevs.extend(em_trainer.fragStddevs)
        # we will round to 1 decimal digit
        for i in range(len(em_means)):
            em_means[ i ] = round(em_means[ i ], 1)
        for i in range(len(em_stddevs)):
            em_stddevs[ i ] = round(em_stddevs[ i ], 1)
        options.info( f"#  The means and stddevs after EM:")

    options.info(  "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s} {0[3]:>10s}".format( ["short", "mono", "di", "tri"] ) )
    options.info(  "#             means: {0[0]:>10.4g} {0[1]:>10.4g} {0[2]:>10.4g} {0[3]:>10.4g}".format( em_means ) )
    options.info(  "#           stddevs: {0[0]:>10.4g} {0[1]:>10.4g} {0[2]:>10.4g} {0[3]:>10.4g}".format( em_stddevs ) )    

    # to finalize the EM training, we will decompose ATAC-seq into four signal tracks
    options.info( f"#  Compute the weights for each fragment length for each of the four signal types")
    fl_dict = petrack.count_fraglengths()
    fl_list = list(fl_dict.keys())
    fl_list.sort()

    # now we will prepare the weights for each fragment length for
    # each of the four distributions based on the EM results
    weight_mapping = generate_weight_mapping( fl_list, em_means, em_stddevs )
    
    options.info( f"#  Generate short, mono-, di-, and tri-nucleosomal signals")
    digested_atac_signals = generate_digested_signals( petrack, weight_mapping )

    # save three types of signals if needed
    if options.save_digested:
        fhd = open(short_bdgfile,"w")
        digested_atac_signals[ 0 ].write_bedGraph(fhd, "short","short")
        fhd.close()
        fhd = open(mono_bdgfile,"w")
        digested_atac_signals[ 0 ].write_bedGraph(fhd, "mono","mono")
        fhd.close()
        fhd = open(di_bdgfile,"w")
        digested_atac_signals[ 0 ].write_bedGraph(fhd, "di","di")
        fhd.close()
        fhd = open(tri_bdgfile,"w")
        digested_atac_signals[ 0 ].write_bedGraph(fhd, "tri","tri")
        fhd.close()        

    minlen = int(petrack.average_template_length)
    # if options.pileup_short is on, we pile up only the short fragments to identify training 
    #  regions and to prescan for candidate regions for decoding.
    if options.pileup_short:
        options.info( f"#  Pile up ONLY short fragments" )
        fc_bdg = digested_atac_signals[ 0 ]
        (sum_v, n_v, max_v, min_v, mean_v, std_v) = fc_bdg.summary()
        print(sum_v, n_v, max_v, min_v, mean_v, std_v)
        options.info( f"#  Convert pileup to fold-change over average signal" )
        fc_bdg.apply_func(lambda x: x/mean_v)
    else:
        options.info( f"#  Pile up all fragments" )
        fc_bdg = petrack.pileup_bdg( [1.0,], baseline_value = 0 )
        (sum_v, n_v, max_v, min_v, mean_v, std_v) = fc_bdg.summary()
        options.info( f"#  Convert pileup to fold-change over average signal" )
        fc_bdg.apply_func(lambda x: x/mean_v)
       
        
    # if cutoff_analysis only, generate and save the report and quit
    if options.cutoff_analysis_only:
        # we will run cutoff analysis only and quit
        options.info( f"#3 Generate cutoff analysis report from {petrack.total} fragments")
        options.info( f"#   Please review the cutoff analysis result in {cutoffanalysis_file}" )

        # Let MACS3 do the cutoff analysis to help decide the lower and upper cutoffs
        with open(cutoffanalysis_file, "w") as ofhd_cutoff:
            ofhd_cutoff.write( fc_bdg.cutoff_analysis( min_length=minlen, max_gap=options.hmm_training_flanking ) )
        #raise Exception("Cutoff analysis only.")
        sys.exit(1)
        
        
    #############################################
    # 3. Define training set by peak calling
    #############################################

    if options.hmm_file:
        # skip this step if hmm_file is given
        options.info( f"#3 Skip this step of looking for training set since a Hidden Markov Model file has been provided!")
    else:
        # Find regions with fold change within determined range to use as training sites.
        # Find regions with zscore values above certain cutoff to exclude from viterbi.
        # 
        options.info( f"#3 Look for training set from {petrack.total} fragments" )
        options.info( f"#  Call peak above within fold-change range of {options.hmm_lower} and {options.hmm_upper}." )
        options.info( f"#   The minimum length of the region is set as the average template/fragment length in the dataset: {minlen}" )
        options.info( f"#   The maximum gap to merge nearby significant regions is set as the flanking size to extend training regions: {options.hmm_training_flanking}" )    
        peaks = fc_bdg.call_peaks (cutoff=options.hmm_lower, min_length=minlen, max_gap=options.hmm_training_flanking, call_summits=False)
        options.info( f"#  Total training regions called after applying the lower cutoff {options.hmm_lower}: {peaks.total}" )
        peaks.filter_score( options.hmm_lower, options.hmm_upper )
        options.info( f"#  Total training regions after filtering with upper cutoff {options.hmm_upper}: {peaks.total}" )

        options.info( f"#  **IMPORTANT**")
        options.info( f"#  Please review the cutoff analysis result in {cutoffanalysis_file} to verify" )
        options.info( f"#   if the choices of lower, upper and prescanning cutoff are appropriate." )
        options.info( f"#   Please read the message in the section 'Choices of cutoff values' by running" )
        options.info( f"#   `macs3 hmmratac -h` for detail." )
        options.info( f"#  ****" )
        
        # Let MACS3 do the cutoff analysis to help decide the lower and upper cutoffs
        with open(cutoffanalysis_file, "w") as ofhd_cutoff:
            ofhd_cutoff.write( fc_bdg.cutoff_analysis( min_length=minlen, max_gap=options.hmm_training_flanking ) )
            
        # we will check if anything left after filtering
        if peaks.total > options.hmm_maxTrain:
            peaks = peaks.randomly_pick( options.hmm_maxTrain, seed = options.hmm_randomSeed )
            options.info( f"#  We randomly pick {options.hmm_maxTrain} regions for training" )
        elif peaks.total == 0:
            options.error( f"# No training regions found. Please adjust the lower or upper cutoff." )
            raise Exception("Not enough training regions!")
        
        # Now we convert PeakIO to Regions and filter blacklisted regions
        training_regions = Regions()
        training_regions.init_from_PeakIO( peaks )
        # We will expand the regions to both directions and merge overlap
        options.info( f"#  We expand the training regions with {options.hmm_training_flanking} basepairs and merge overlap" )
        training_regions.expand( options.hmm_training_flanking )
        training_regions.merge_overlap()
    
        # remove peaks overlapping with blacklisted regions
        if options.blacklist:
            training_regions.exclude( blacklist_regions )
            options.info( f"#  after removing those overlapping with provided blacklisted regions, we have {training_regions.total} left" )
        if options.save_train:
            fhd = open( training_region_bedfile, "w" )
            training_regions.write_to_bed( fhd )
            fhd.close()
            options.info( f"#  Training regions have been saved to `{options.name}_training_regions.bed` " )
    
    #############################################
    # 4. Train HMM
    #############################################
    # if model file is provided, we skip this step
    if options.hmm_file:
        options.info( f"#4 Load Hidden Markov Model from given model file")
        hmm_model, i_open_region, i_background_region, i_nucleosomal_region, options.hmm_binsize = hmm_model_init( options.hmm_file )
    else:
        options.info( f"#4 Train Hidden Markov Model with Multivariate Gaussian Emission" )

        # extract signals within peak using the given binsize
        options.info( f"#  Extract signals in training regions with bin size of {options.hmm_binsize}")
        [ training_bins, training_data, training_data_lengths ] = extract_signals_from_regions( digested_atac_signals, training_regions, binsize = options.hmm_binsize )

        if options.save_train:
            f = open( training_datafile, "w" )
            for i in range( len( training_data ) ):
                v = training_data[ i ]
                p = training_bins[ i ]
                f.write( f"{p[0]}\t{p[1]}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\n" )
            f.close()

            f = open( training_datalengthfile, "w" )
            for v in training_data_lengths:
                f.write( f"{v}\n" )
            f.close()

        options.info( f"#  Use Baum-Welch algorithm to train the HMM")

        hmm_model = hmm_training( training_data, training_data_lengths, random_seed = options.hmm_randomSeed, covar="full" )

        options.info( f"#   HMM converged: {hmm_model.monitor_.converged}")

        # label hidden states
        means_sum = np.sum( hmm_model.means_, axis=1 )

        # first, the state with the highest overall emission is the open state
        i_open_region = np.where( means_sum == max(means_sum) )[0][0]

        # second, the state with lowest overall emission is the bg state 
        i_background_region = np.where( means_sum == min(means_sum) )[0][0]

        # last one is the nuc state (note it may not be accurate though
        i_nucleosomal_region = list(set([0, 1, 2]) - set([i_open_region, i_background_region]))[0]

        # write hmm into model file
        options.info( f"#  Write HMM parameters into JSON: {hmm_modelfile}")
        hmm_model_save( hmm_modelfile, hmm_model, options.hmm_binsize, i_open_region, i_nucleosomal_region, i_background_region )

    # Now tell users the parameters of the HMM
    assignments = [ "", "", "" ]
    assignments[ i_open_region ]        = "open"
    assignments[ i_nucleosomal_region ] = "nuc"
    assignments[ i_background_region ]  = "bg"
    
    options.info( f"#  The Hidden Markov Model for signals of binsize of {options.hmm_binsize} basepairs:")
    options.info( f"#   open state index: state{i_open_region}" )
    options.info( f"#   nucleosomal state index: state{i_nucleosomal_region}" )
    options.info( f"#   background state index: state{i_background_region}" )
    options.info( f"#   Starting probabilities of states:")
    options.info(  "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s}".format( assignments ) )
    options.info(  "#                    {0[0]:>10.4g} {0[1]:>10.4g} {0[2]:>10.4g}".format( hmm_model.startprob_ ) )
    options.info( f"#   HMM Transition probabilities:")
    options.info(  "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s}".format( assignments ) )
    options.info(  "#       {0:>10s}-> {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g}".format(assignments[0], hmm_model.transmat_[0]) )
    options.info(  "#       {0:>10s}-> {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g}".format(assignments[1], hmm_model.transmat_[1]) )
    options.info(  "#       {0:>10s}-> {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g}".format(assignments[2], hmm_model.transmat_[2]) )
    options.info( f"#   HMM Emissions (mean): ")
    options.info(  "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s} {0[3]:>10s}".format( ["short", "mono", "di", "tri"] ) )
    options.info(  "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[0], hmm_model.means_[0]) )
    options.info(  "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[1], hmm_model.means_[1]) )
    options.info(  "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[2], hmm_model.means_[2]) )
    

#############################################
# 5. Predict
#############################################
    # Our prediction strategy will be different with HMMRATAC, we will first ask MACS call peaks with loose cutoff, then for each peak we will run HMM prediction to figure out labels. And for the rest of genomic regions, just mark them as 'background'.
    options.info( f"#5 Decode with Viterbi to predict states" )
    # the following /4 is totally arbitrary, we may need to fix it
    candidate_peaks = fc_bdg.call_peaks (cutoff=options.prescan_cutoff, min_length=minlen, max_gap=options.hmm_training_flanking, call_summits=False)
    options.info( f"#5  Total candidate peaks : {candidate_peaks.total}" )

    # Now we convert PeakIO to Regions and filter blacklisted regions
    candidate_regions = Regions()
    candidate_regions.init_from_PeakIO( candidate_peaks )
    # We will expand the regions to both directions and merge overlap
    options.info( f"#  We expand the candidate regions with {options.hmm_training_flanking} and merge overlap" )
    candidate_regions.expand( options.hmm_training_flanking )
    candidate_regions.merge_overlap()
    
    # remove peaks overlapping with blacklisted regions
    if options.blacklist:
        candidate_regions.exclude( blacklist_regions )
        options.info( f"#  after removing those overlapping with provided blacklisted regions, we have {candidate_regions.total} left" )

    # extract signals
    options.info( f"#  Extract signals in candidate regions and decode with HMM")
    # we will do the extraction and prediction in a step of 10000 regions by default
    
    # Note: we can implement in a different way to extract then predict for each candidate region.
    # predicted_results = hmm_decode_each_region ( digested_atac_signals, candidate_regions, hmm_model, binsize = options.hmm_binsize )
    # Note: we implement in a way that we will decode the candidate regions 10000 regions at a time so 1. we can make it running in parallel in the future; 2. we can reduce the memory usage.
    options.info( f"#  Use HMM to predict states")
    n = 0
    predicted_proba = []
    candidate_bins = []
    while candidate_regions.total != 0:
        n += 1
        cr = candidate_regions.pop( options.decoding_steps )
        options.info( "#    decoding %d..." % ( n*options.decoding_steps ) )        
        [ cr_bins, cr_data, cr_data_lengths ] = extract_signals_from_regions( digested_atac_signals, cr, binsize = options.hmm_binsize )
        #options.info( "#     extracted signals in the candidate regions")
        candidate_bins.extend( cr_bins )
        #options.info( "#     saving information regarding the candidate regions")        
        predicted_proba.extend( hmm_predict( cr_data, cr_data_lengths, hmm_model ) )
        #options.info( "#     prediction done")
        gc.collect()

        
#############################################
# 6. Output - add to OutputWriter
#############################################
    options.info( f"# Write the output...")
    # Now taken the candidate_bins and predicted_proba, we can generate various
    # outputs
    
    # One thing to remember about candidate_bins is that the position
    # in this array is the 'end' of the bin, the actual region is the
    # 'end'-'binsize' to the 'end'.
    
    # First, the likelihoods for each of the three states in a bedGraph
    if options.save_likelihoods:
        options.info( f"# Write the likelihoods for each states into three bedGraph files {options.name}_open.bdg, {options.name}_nuc.bdg, and {options.name}_bg.bdg")
        open_state_bdg_fhd = open( open_state_bdgfile, "w" )
        nuc_state_bdg_fhd = open( nuc_state_bdgfile, "w" )
        bg_state_bdg_fhd = open( bg_state_bdgfile, "w" )
        save_proba_to_bedGraph( candidate_bins, predicted_proba, options.hmm_binsize, open_state_bdg_fhd, nuc_state_bdg_fhd, bg_state_bdg_fhd, i_open_region, i_nucleosomal_region, i_background_region )
        open_state_bdg_fhd.close()
        nuc_state_bdg_fhd.close()
        bg_state_bdg_fhd.close()
    
    # Generate states path:
    states_path = generate_states_path( candidate_bins, predicted_proba, options.hmm_binsize, i_open_region, i_nucleosomal_region, i_background_region )
    
    # Save states path if needed
    # PS: we need to implement extra feature to include those regions NOT in candidate_bins and assign them as 'background state'.
    if options.save_states:
        options.info( f"# Write states assignments in a BED file: {options.name}_states.bed" )
        f = open( states_file, "w" )
        save_states_bed( states_path, f )
        f.close()

    options.info( f"# Write accessible regions in a gappedPeak file: {options.name}_accessible_regions.gappedPeak")
    ofhd = open( accessible_file, "w" )
    save_accessible_regions( states_path, ofhd, options.openregion_minlen )
    ofhd.close()

def save_proba_to_bedGraph( candidate_bins, predicted_proba, binsize, open_state_bdg_file, nuc_state_bdg_file, bg_state_bdg_file, i_open, i_nuc, i_bg ):
    open_state_bdg = bedGraphTrackI( baseline_value = 0 )
    nuc_state_bdg = bedGraphTrackI( baseline_value = 0 )
    bg_state_bdg = bedGraphTrackI( baseline_value = 0 )

    prev_chrom_name = None
    prev_bin_end = None
    for l in range(len(predicted_proba)):
        # note that any region not in the candidate bins will be
        # treated as absolutely the background
        chrname = candidate_bins[l][0]
        end_pos = candidate_bins[l][1]
        start_pos = end_pos - binsize
        if chrname != prev_chrom_name:
            prev_chrom_name = chrname
            # add the first region as background
            if start_pos > 0:
                open_state_bdg.add_loc( chrname, 0, start_pos, 0.0 )
                nuc_state_bdg.add_loc( chrname, 0, start_pos, 0.0 )
                bg_state_bdg.add_loc( chrname, 0, start_pos, 1.0 )
                prev_bin_end = start_pos
            elif start_pos == 0:
                # if start_pos == 0, then the first bin has to be assigned, we set prev_bin_end as 0 
                prev_bin_end = 0
        # now check if the prev_bin_end is start_pos, if not, add a gap of background
        if prev_bin_end < start_pos:
            open_state_bdg.add_loc( chrname, prev_bin_end, start_pos, 0.0 )
            nuc_state_bdg.add_loc( chrname, prev_bin_end, start_pos, 0.0 )
            bg_state_bdg.add_loc( chrname, prev_bin_end, start_pos, 1.0 )

        open_state_bdg.add_loc( chrname, start_pos, end_pos, predicted_proba[l][ i_open ] )
        nuc_state_bdg.add_loc( chrname, start_pos, end_pos, predicted_proba[l][ i_nuc ] )
        bg_state_bdg.add_loc( chrname, start_pos, end_pos, predicted_proba[l][ i_bg ] )
        prev_bin_end = start_pos

    open_state_bdg.write_bedGraph( open_state_bdg_file, "Open States", "Likelihoods of being Open States" )
    nuc_state_bdg.write_bedGraph( nuc_state_bdg_file, "Nucleosomal States", "Likelihoods of being Nucleosomal States" )
    bg_state_bdg.write_bedGraph( bg_state_bdg_file, "Background States", "Likelihoods of being Background States" )
    return

def save_states_bed( states_path, states_bedfile ):
    # we do not need to output background state. 
    for l in range( len( states_path ) ):
        if states_path[l][3] != "bg":
            states_bedfile.write( "%s\t%s\t%s\t%s\n" % states_path[l] )
    return

def generate_states_path(candidate_bins, predicted_proba, binsize, i_open_region, i_nucleosomal_region, i_background_region):
    ret_states_path = []
    labels_list = ["open", "nuc", "bg"]
    start_pos = candidate_bins[0][1] - binsize
    for l in range(1, len(predicted_proba)):
        chromosome = candidate_bins[l][0].decode()
        prev_open, prev_nuc, prev_bg = predicted_proba[l-1][i_open_region], predicted_proba[l-1][i_nucleosomal_region], predicted_proba[l-1][i_background_region]
        curr_open, curr_nuc, curr_bg = predicted_proba[l][i_open_region], predicted_proba[l][i_nucleosomal_region], predicted_proba[l][i_background_region]
        
        label_prev = labels_list[max((prev_open, 0), (prev_nuc, 1), (prev_bg, 2), key=lambda x: x[0])[1]]
        label_curr = labels_list[max((curr_open, 0), (curr_nuc, 1), (curr_bg, 2), key=lambda x: x[0])[1]]

        if candidate_bins[l-1][0] == candidate_bins[l][0]:  # if we are looking at the same chromosome ...
            if label_prev != label_curr:
                end_pos = candidate_bins[l-1][1]
                ret_states_path.append((chromosome, start_pos, end_pos, label_prev))
                start_pos = candidate_bins[l][1] - binsize
            elif l == len(predicted_proba) - 1:
                end_pos = candidate_bins[l][1]
                ret_states_path.append((chromosome, start_pos, end_pos, label_prev))
        else:
            start_pos = candidate_bins[l][1] - binsize
    return ret_states_path

def save_accessible_regions(states_path, accessible_region_file, openregion_minlen):
    # Function to add regions to the list
    def add_regions(i, regions):
        for j in range(i, i+3):
            if not regions or states_path[j][2] != regions[-1][2]:
                regions.append((states_path[j][0], int(states_path[j][1]), int(states_path[j][2]), states_path[j][3]))
        return regions

    # Select only accessible regions from _states.bed, look for nuc-open-nuc pattern
    # This by default is the only final output from HMMRATAC
    accessible_regions = []
    for i in range(len(states_path)-2):
        if (states_path[i][3] == 'nuc' and states_path[i+1][3] == 'open' and states_path[i+2][3] == 'nuc' and 
            states_path[i][2] == states_path[i+1][1] and states_path[i+1][2] == states_path[i+2][1] and 
            states_path[i+1][2] - states_path[i+1][1] > openregion_minlen):
            accessible_regions = add_regions(i, accessible_regions)

    # Group states by region
    list_of_groups = []
    one_group = [accessible_regions[0]]
    for j in range(1, len(accessible_regions)):
        if accessible_regions[j][1] == accessible_regions[j-1][2]:
            one_group.append(accessible_regions[j])
        else:
            list_of_groups.append(one_group)
            one_group = [accessible_regions[j]]
    accessible_regions = list_of_groups

    # Generate broadpeak object
    broadpeak = BroadPeakIO()
    for region in accessible_regions[:-1]:
        block_num = sum('open' in tup for tup in region)
        block_sizes = ','.join(str(region[j][2] - region[j][1]) for j in range(1, len(region) - 1, 2))
        block_starts = ','.join(str(region[j][1] - region[0][1]) for j in range(1, len(region) - 1, 2))
        broadpeak.add(bytes(region[1][0], encoding="raw_unicode_escape"), region[0][1], region[-1][2],
                      thickStart=bytes(str(region[1][1]), encoding="raw_unicode_escape"),
                      thickEnd=bytes(str(region[-2][2]), encoding="raw_unicode_escape"),
                      blockNum=block_num,
                      blockSizes=bytes(block_sizes, encoding="raw_unicode_escape"),
                      blockStarts=bytes(block_starts, encoding="raw_unicode_escape"))
    broadpeak.write_to_gappedPeak(accessible_region_file)
    return
