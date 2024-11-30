# Time-stamp: <2024-11-29 23:19:38 Tao Liu>


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
import tempfile

# ------------------------------------
# own python modules
# ------------------------------------
# from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_hmmratac
from MACS3.IO.PeakIO import PeakIO
from MACS3.IO.Parser import BAMPEParser, BEDPEParser, FragParser  # BAMaccessor
from MACS3.Signal.HMMR_EM import HMMR_EM
from MACS3.Signal.HMMR_Signal_Processing import (generate_weight_mapping,
                                                 generate_digested_signals,
                                                 extract_signals_from_regions)
from MACS3.Signal.HMMR_HMM import (hmm_training,
                                   hmm_predict,
                                   hmm_model_init,
                                   hmm_model_save)
from MACS3.Signal.Region import Regions
from MACS3.IO.BedGraphIO import bedGraphIO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------


def run(args):
    """The HMMRATAC function/pipeline for MACS.

    """
    options = opt_validate_hmmratac(args)

    #############################################
    # 0. Names of output files
    #############################################
    short_bdgfile = os.path.join(options.outdir,
                                 options.name+"_digested_short.bdg")
    mono_bdgfile = os.path.join(options.outdir,
                                options.name+"_digested_mono.bdg")
    di_bdgfile = os.path.join(options.outdir,
                              options.name+"_digested_di.bdg")
    tri_bdgfile = os.path.join(options.outdir,
                               options.name+"_digested_tri.bdg")
    training_region_bedfile = os.path.join(options.outdir,
                                           options.name+"_training_regions.bed")
    training_datafile = os.path.join(options.outdir,
                                     options.name+"_training_data.txt")
    training_datalengthfile = os.path.join(options.outdir,
                                           options.name+"_training_lengths.txt")

    hmm_modelfile = os.path.join(options.outdir,
                                 options.name+"_model.json")

    open_state_bdgfile = os.path.join(options.outdir,
                                      options.name+"_open.bdg")
    nuc_state_bdgfile = os.path.join(options.outdir,
                                     options.name+"_nuc.bdg")
    bg_state_bdgfile = os.path.join(options.outdir,
                                    options.name+"_bg.bdg")

    states_file = os.path.join(options.outdir,
                               options.name+"_states.bed")

    accessible_file = os.path.join(options.outdir,
                                   options.name+"_accessible_regions.narrowPeak")

    cutoffanalysis_file = os.path.join(options.outdir,
                                       options.name+"_cutoff_analysis.tsv")

    #############################################
    # 1. Read the input files
    #############################################
    options.info("\n" + options.argtxt)

    if options.format == "BAMPE":
        options.info("#1 Read fragments from BAMPE file...")
        parser = BAMPEParser
    elif options.format == "BEDPE":
        options.info("#1 Read fragments from BEDPE file...")
        parser = BEDPEParser
    elif options.format == "FRAG":
        options.info("#1 Read fragments from FRAG file...")
        parser = FragParser
    else:
        raise Exception("wrong format")

    alignment = parser(options.input_file[0], buffer_size=options.buffer_size)
    petrack = alignment.build_petrack()
    if len(options.input_file) > 1:
        # multiple input
        for inputfile in options.input_file[1:]:
            alignment = parser(inputfile, buffer_size=options.buffer_size)
            petrack = alignment.append_petrack(petrack)
    # remember to finalize the petrack
    petrack.finalize()

    # filter duplicates if needed
    if options.misc_keep_duplicates:
        petrack.filter_dup(maxnum=1)

    # read in blacklisted if option entered
    if options.blacklist:
        options.info("#  Read blacklist file...")
        peakio = open(options.blacklist)
        blacklist = PeakIO()
        i = 0
        for l_p in peakio:
            fs = l_p.rstrip().split()
            i += 1
            blacklist.add(fs[0].encode(),
                          int(fs[1]),
                          int(fs[2]),
                          name=b"%d" % i)
            blacklist.sort()
        blacklist_regions = Regions()
        blacklist_regions.init_from_PeakIO(blacklist)

    #############################################
    # 2. EM
    #############################################
    if options.em_skip:
        # Skip EM and use the options.em_means and options.em_stddevs
        em_means = options.em_means
        em_stddevs = options.em_stddevs
        options.info("#2 EM is skipped. The following means and stddevs will be used:")
    else:
        # we will use EM to get the best means/stddevs for the mono-, di- and tri- modes of fragment sizes
        options.info("#2 Use EM algorithm to estimate means and stddevs of fragment lengths")
        options.info("#  for mono-, di-, and tri-nucleosomal signals...")
        em_trainer = HMMR_EM(petrack, options.em_means[1:4],
                             options.em_stddevs[1:4], seed=options.hmm_randomSeed)
        # the mean and stddev after EM training
        em_means = [options.em_means[0],]
        em_means.extend(em_trainer.fragMeans)
        em_stddevs = [options.em_stddevs[0],]
        em_stddevs.extend(em_trainer.fragStddevs)
        # we will round to 1 decimal digit
        for i in range(len(em_means)):
            em_means[i] = round(em_means[i], 1)
        for i in range(len(em_stddevs)):
            em_stddevs[i] = round(em_stddevs[i], 1)
        options.info("#  The means and stddevs after EM:")

    options.info("#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s} {0[3]:>10s}".format(["short", "mono", "di", "tri"]))
    options.info("#             means: {0[0]:>10.4g} {0[1]:>10.4g} {0[2]:>10.4g} {0[3]:>10.4g}".format(em_means))
    options.info("#           stddevs: {0[0]:>10.4g} {0[1]:>10.4g} {0[2]:>10.4g} {0[3]:>10.4g}".format(em_stddevs))

    # to finalize the EM training, we will decompose ATAC-seq into four signal tracks
    options.info("#  Compute the weights for each fragment length for each of the four signal types")
    fl_dict = petrack.count_fraglengths()
    fl_list = list(fl_dict.keys())
    fl_list.sort()

    # now we will prepare the weights for each fragment length for
    # each of the four distributions based on the EM results
    weight_mapping = generate_weight_mapping(fl_list, em_means, em_stddevs,
                                             min_frag_p=options.min_frag_p)

    options.info("#  Generate short, mono-, di-, and tri-nucleosomal signals")
    digested_atac_signals = generate_digested_signals(petrack, weight_mapping)

    # save three types of signals if needed
    if options.save_digested:
        bdgshort = bedGraphIO(short_bdgfile, data=digested_atac_signals[0])
        bdgshort.write_bedGraph("short", "short")

        bdgmono = bedGraphIO(mono_bdgfile, data=digested_atac_signals[1])
        bdgmono.write_bedGraph("mono", "mono")

        bdgdi = bedGraphIO(di_bdgfile, data=digested_atac_signals[2])
        bdgdi.write_bedGraph("di", "di")

        bdgtri = bedGraphIO(tri_bdgfile, data=digested_atac_signals[3])
        bdgtri.write_bedGraph("tri", "tri")

    minlen = int(petrack.average_template_length)
    # if options.pileup_short is on, we pile up only the short
    #  fragments to identify training regions and to prescan for
    #  candidate regions for decoding.
    if options.pileup_short:
        options.info("#  Pile up ONLY short fragments")
        fc_bdg = digested_atac_signals[0]
    else:
        options.info("#  Pile up all fragments")
        fc_bdg = petrack.pileup_bdg([1.0, ], baseline_value=0)
    (sum_v, n_v, max_v, min_v, mean_v, std_v) = fc_bdg.summary()
    options.info("#  Convert pileup to fold-change over average signal")
    fc_bdg.apply_func(lambda x: x/mean_v)

    # if cutoff_analysis only, generate and save the report and quit
    if options.cutoff_analysis_only:
        # we will run cutoff analysis only and quit
        options.info(f"#3 Generate cutoff analysis report from {petrack.total} fragments")
        options.info(f"#   Please review the cutoff analysis result in {cutoffanalysis_file}")

        # Let MACS3 do the cutoff analysis to help decide the lower
        # and upper cutoffs
        with open(cutoffanalysis_file, "w") as ofhd_cutoff:
            ofhd_cutoff.write(fc_bdg.cutoff_analysis(min_length=minlen,
                                                     max_gap=options.hmm_training_flanking,
                                                     max_score=options.cutoff_analysis_max,
                                                     steps=options.cutoff_analysis_steps))
        # raise Exception("Cutoff analysis only.")
        sys.exit(1)

    #############################################
    # 3. Define training set by peak calling
    #############################################

    if options.hmm_file:
        # skip this step if hmm_file is given
        options.info("#3 Skip this step of looking for training set since a Hidden Markov Model file has been provided!")
    elif options.hmm_training_regions:
        # if a training region file is provided - need to read in the bedfile and skip the peak calling step
        options.info(f"#3 Read training regions from BED file: {options.hmm_training_regions}")
        # from refinepeak_cmd.py:
        peakio = open(options.hmm_training_regions, "rb")
        peaks = PeakIO()
        for l_p in peakio:
            fs = l_p.rstrip().split()
            peaks.add(chromosome=fs[0], start=int(fs[1]), end=int(fs[2]))  # change based on what expected input file should contain
        peakio.close()
        training_regions = Regions()
        training_regions.init_from_PeakIO(peaks)
        options.info("#  Training regions have been read from bedfile")
    else:
        # Find regions with fold change within determined range to use as training sites.
        # Find regions with zscore values above certain cutoff to exclude from viterbi.
        #
        options.info(f"#3 Look for training set from {petrack.total} fragments")
        options.info(f"#  Call peak above within fold-change range of {options.hmm_lower} and {options.hmm_upper}.")
        options.info(f"#   The minimum length of the region is set as the average template/fragment length in the dataset: {minlen}")
        options.info(f"#   The maximum gap to merge nearby significant regions is set as the flanking size to extend training regions: {options.hmm_training_flanking}")
        peaks = fc_bdg.call_peaks(cutoff=options.hmm_lower, min_length=minlen, max_gap=options.hmm_training_flanking, call_summits=False)
        options.info(f"#  Total training regions called after applying the lower cutoff {options.hmm_lower}: {peaks.total}")
        peaks.filter_score(options.hmm_lower, options.hmm_upper)
        options.info(f"#  Total training regions after filtering with upper cutoff {options.hmm_upper}: {peaks.total}")

        options.info( "#  **IMPORTANT**")
        options.info(f"#  Please review the cutoff analysis result in {cutoffanalysis_file} to verify")
        options.info( "#   if the choices of lower, upper and prescanning cutoff are appropriate.")
        options.info( "#   Please read the message in the section 'Choices of cutoff values' by running")
        options.info( "#   `macs3 hmmratac -h` for detail.")
        options.info( "#  ****")

        # Let MACS3 do the cutoff analysis to help decide the lower and upper cutoffs
        with open(cutoffanalysis_file, "w") as ofhd_cutoff:
            ofhd_cutoff.write(fc_bdg.cutoff_analysis(min_length=minlen, max_gap=options.hmm_training_flanking, max_score=options.cutoff_analysis_max))

        # we will check if anything left after filtering
        if peaks.total > options.hmm_maxTrain:
            peaks = peaks.randomly_pick(options.hmm_maxTrain, seed=options.hmm_randomSeed)
            options.info(f"#  We randomly pick {options.hmm_maxTrain} regions for training")
        elif peaks.total == 0:
            options.error("# No training regions found. Please adjust the lower or upper cutoff.")
            raise Exception("Not enough training regions!")

        # Now we convert PeakIO to Regions and filter blacklisted regions
        training_regions = Regions()
        training_regions.init_from_PeakIO(peaks)
        # We will expand the regions to both directions and merge overlap
        options.info(f"#  We expand the training regions with {options.hmm_training_flanking} basepairs and merge overlap")
        training_regions.expand(options.hmm_training_flanking)
        training_regions.merge_overlap()

        # remove peaks overlapping with blacklisted regions
        if options.blacklist:
            training_regions.exclude(blacklist_regions)
            options.info(f"#  after removing those overlapping with provided blacklisted regions, we have {training_regions.total} left")
        if options.save_train:
            fhd = open(training_region_bedfile, "w")
            training_regions.write_to_bed(fhd)
            fhd.close()
            options.info(f"#  Training regions have been saved to `{options.name}_training_regions.bed` ")

    #############################################
    # 4. Train HMM
    #############################################
    # if model file is provided, we skip this step
    # include options.hmm_type and make it backwards compatible, if no hmm_type default is gaussian
    if options.hmm_file:
        options.info("#4 Load Hidden Markov Model from given model file")
        hmm_model, i_open_region, i_background_region, i_nucleosomal_region, options.hmm_binsize, options.hmm_type = hmm_model_init(options.hmm_file)
    else:
        options.info("#4 Train Hidden Markov Model with Multivariate Gaussian Emission")

        # extract signals within peak using the given binsize
        options.info(f"#  Extract signals in training regions with bin size of {options.hmm_binsize}")
        [training_bins, training_data, training_data_lengths] = extract_signals_from_regions(digested_atac_signals, training_regions, binsize=options.hmm_binsize, hmm_type=options.hmm_type)

        if options.save_train:
            f = open(training_datafile, "w")
            for i in range(len(training_data)):
                v = training_data[i]
                p = training_bins[i]
                f.write(f"{p[0]}\t{p[1]}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\n")
            f.close()

            f = open(training_datalengthfile, "w")
            for v in training_data_lengths:
                f.write(f"{v}\n")
            f.close()

        options.info("#  Use Baum-Welch algorithm to train the HMM")

        hmm_model = hmm_training(training_data,
                                 training_data_lengths,
                                 random_seed=options.hmm_randomSeed,
                                 hmm_type=options.hmm_type)

        options.info(f"#   HMM converged: {hmm_model.monitor_.converged}")

        # label hidden states
        if options.hmm_type == "gaussian":
            means_sum = np.sum(hmm_model.means_, axis=1)
        if options.hmm_type == "poisson":
            means_sum = np.sum(hmm_model.lambdas_, axis=1)

        # first, the state with the highest overall emission is the open state
        i_open_region = np.where(means_sum == max(means_sum))[0][0]

        # second, the state with lowest overall emission is the bg state
        i_background_region = np.where(means_sum == min(means_sum))[0][0]

        # last one is the nuc state (note it may not be accurate though
        i_nucleosomal_region = list(set([0, 1, 2]) - set([i_open_region, i_background_region]))[0]

        # write hmm into model file
        options.info(f"#  Write HMM parameters into JSON: {hmm_modelfile}")
        hmm_model_save(hmm_modelfile, hmm_model, options.hmm_binsize, i_open_region, i_nucleosomal_region, i_background_region, options.hmm_type)

        # if --modelonly option provided, exit script after hmm model is saved
        if options.hmm_modelonly:
            options.info("#  Complete - HMM model was saved, program exited (--modelonly option was provided) ")
            sys.exit()

    # Now tell users the parameters of the HMM
    assignments = ["", "", ""]
    assignments[i_open_region] = "open"
    assignments[i_nucleosomal_region] = "nuc"
    assignments[i_background_region] = "bg"

    options.info(f"#  The Hidden Markov Model for signals of binsize of {options.hmm_binsize} basepairs:")
    options.info(f"#   open state index: state{i_open_region}")
    options.info(f"#   nucleosomal state index: state{i_nucleosomal_region}")
    options.info(f"#   background state index: state{i_background_region}")
    options.info( "#   Starting probabilities of states:")
    options.info( "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s}".format(assignments))
    options.info( "#                    {0[0]:>10.4g} {0[1]:>10.4g} {0[2]:>10.4g}".format(hmm_model.startprob_))
    options.info( "#   HMM Transition probabilities:")
    options.info( "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s}".format(assignments))
    options.info( "#       {0:>10s}-> {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g}".format(assignments[0], hmm_model.transmat_[0]))
    options.info( "#       {0:>10s}-> {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g}".format(assignments[1], hmm_model.transmat_[1]))
    options.info( "#       {0:>10s}-> {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g}".format(assignments[2], hmm_model.transmat_[2]))

    if options.hmm_type == 'gaussian':
        options.info("#   HMM Emissions (means): ")
        options.info( "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s} {0[3]:>10s}".format(["short", "mono", "di", "tri"]))
        options.info( "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[0], hmm_model.means_[0]))
        options.info( "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[1], hmm_model.means_[1]))
        options.info( "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[2], hmm_model.means_[2]))
    if options.hmm_type == 'poisson':
        options.info( "#   HMM Emissions (lambdas): ")
        options.info( "#                    {0[0]:>10s} {0[1]:>10s} {0[2]:>10s} {0[3]:>10s}".format(["short", "mono", "di", "tri"]))
        options.info( "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[0], hmm_model.lambdas_[0]))
        options.info( "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[1], hmm_model.lambdas_[1]))
        options.info( "#       {0:>10s}:  {1[0]:>10.4g} {1[1]:>10.4g} {1[2]:>10.4g} {1[3]:>10.4g}".format(assignments[2], hmm_model.lambdas_[2]))


#############################################
# 5. Predict
#############################################
    # Our prediction strategy will be different with HMMRATAC, we will first ask MACS call peaks with loose cutoff, then for each peak we will run HMM prediction to figure out labels. And for the rest of genomic regions, just mark them as 'background'.
    options.info("#5 Decode with Viterbi to predict states")
    # the following /4 is totally arbitrary, we may need to fix it
    candidate_peaks = fc_bdg.call_peaks(cutoff=options.prescan_cutoff, min_length=minlen, max_gap=options.hmm_training_flanking, call_summits=False)
    options.info(f"#5  Total candidate peaks : {candidate_peaks.total}")

    # Now we convert PeakIO to Regions and filter blacklisted regions
    candidate_regions = Regions()
    candidate_regions.init_from_PeakIO(candidate_peaks)
    # We will expand the regions to both directions and merge overlap
    options.info(f"#  We expand the candidate regions with {options.hmm_training_flanking} and merge overlap")
    candidate_regions.expand(options.hmm_training_flanking)
    candidate_regions.merge_overlap()
    options.info(f"#   after expanding and merging, we have {candidate_regions.total} candidate regions")

    # remove peaks overlapping with blacklisted regions
    if options.blacklist:
        candidate_regions.exclude(blacklist_regions)
        options.info(f"#   after removing those overlapping with provided blacklisted regions, we have {candidate_regions.total} left")

    # extract signals
    options.info("#  Extract signals in candidate regions and decode with HMM")
    # we will do the extraction and prediction in a step of 10000 regions by default

    # Note: we can implement in a different way to extract then predict for each candidate region.
    # predicted_results = hmm_decode_each_region (digested_atac_signals, candidate_regions, hmm_model, binsize = options.hmm_binsize)
    # Note: we implement in a way that we will decode the candidate regions 10000 regions at a time so 1. we can make it running in parallel in the future; 2. we can reduce the memory usage.
    options.info("#  Use HMM to predict states")
    n = 0

    # we create a temporary file to save the proba predicted from hmm
    predicted_proba_file = tempfile.TemporaryFile(mode="w+b")

    while candidate_regions.total != 0:
        n += 1
        # we get DECODING_STEPS number of candidate regions first
        cr = candidate_regions.pop(options.decoding_steps)
        options.info("#    decoding %d..." % (n * options.decoding_steps))

        # then extrac data from digested signals, create cr_bins, cr_data, and cr_data_lengths
        [cr_bins, cr_data, cr_data_lengths] = extract_signals_from_regions(digested_atac_signals, cr, binsize=options.hmm_binsize, hmm_type=options.hmm_type)
        # options.debug("#     extract_signals_from_regions complete")

        prob_data = hmm_predict(cr_data, cr_data_lengths, hmm_model)
        assert len(prob_data) == len(cr_bins)
        for i in range(len(prob_data)):
            predicted_proba_file.write(b"%s,%d" % cr_bins[i])
            predicted_proba_file.write(b",%f,%f,%f\n" % tuple(prob_data[i]))

        cr_data = []
        cr_data_lengths = []
        cr_bins = []
        prob_data = []
        gc.collect()

    predicted_proba_file.seek(0)  # reset
    options.info("# predicted_proba files written...")

#############################################
# 6. Output - add to OutputWriter
#############################################
    options.info("# Write the output...")
    # Now taken the candidate_bins and predicted_proba, we can generate various
    # outputs

    # One thing to remember about candidate_bins is that the position
    # in this array is the 'end' of the bin, the actual region is the
    # 'end'-'binsize' to the 'end'.

    # First, the likelihoods for each of the three states in a bedGraph
    if options.save_likelihoods:
        options.info(f"# Write the likelihoods for each states into three bedGraph files {options.name}_open.bdg, {options.name}_nuc.bdg, and {options.name}_bg.bdg")
        open_state_bdg_fhd = open(open_state_bdgfile, "w")
        nuc_state_bdg_fhd = open(nuc_state_bdgfile, "w")
        bg_state_bdg_fhd = open(bg_state_bdgfile, "w")
        save_proba_to_bedGraph(predicted_proba_file, options.hmm_binsize, open_state_bdg_fhd, nuc_state_bdg_fhd, bg_state_bdg_fhd, i_open_region, i_nucleosomal_region, i_background_region)
        predicted_proba_file.seek(0)  # reset
        open_state_bdg_fhd.close()
        nuc_state_bdg_fhd.close()
        bg_state_bdg_fhd.close()
        options.info("# finished writing proba_to_bedgraph")

    # # Generate states path:
    states_path = generate_states_path(predicted_proba_file, options.hmm_binsize, i_open_region, i_nucleosomal_region, i_background_region)
    options.info("# finished generating states path")
    predicted_proba_file.close()          # kill the temp file
    # Save states path if needed
    # PS: we need to implement extra feature to include those regions NOT in candidate_bins and assign them as 'background state'.
    if options.save_states:
        options.info(f"# Write states assignments in a BED file: {options.name}_states.bed")
        with open(states_file, "w") as f:
            save_states_bed(states_path, f)

    options.info(f"# Write accessible regions in a narrowPeak file: {options.name}_accessible_regions.narrowPeak")
    with open(accessible_file, "w") as ofhd:
        save_accessible_regions(states_path, ofhd, options.openregion_minlen, fc_bdg)

    options.info("# Finished")


def save_proba_to_bedGraph(predicted_proba_file, binsize, open_state_bdg_file, nuc_state_bdg_file, bg_state_bdg_file, i_open, i_nuc, i_bg):

    open_state_bdg_file = bedGraphIO(open_state_bdg_file)
    nuc_state_bdg_file = bedGraphIO(nuc_state_bdg_file)
    bg_state_bdg_file = bedGraphIO(bg_state_bdg_file)

    open_state_bdg = open_state_bdg_file.data
    nuc_state_bdg = nuc_state_bdg_file.data
    bg_state_bdg = bg_state_bdg_file.data

    prev_chrom_name = None
    prev_bin_end = None

    for pp_line in predicted_proba_file:
        pp_data = pp_line.strip().split(b',')

        chrname = pp_data[0]
        end_pos = int(pp_data[1])
        start_pos = end_pos - binsize
        pp_open = float(pp_data[i_open+2])
        pp_nuc = float(pp_data[i_nuc+2])
        pp_bg = float(pp_data[i_bg+2])

        if chrname != prev_chrom_name:
            # we start a new chromosome
            if start_pos > 0:
                # add the first unannotated region as background
                open_state_bdg.add_loc(chrname, 0, start_pos, 0.0)
                nuc_state_bdg.add_loc(chrname, 0, start_pos, 0.0)
                bg_state_bdg.add_loc(chrname, 0, start_pos, 1.0)
            prev_chrom_name = chrname
        else:
            # now check if the prev_bin_end is start_pos, if not, add a gap of background
            if prev_bin_end < start_pos:
                open_state_bdg.add_loc(chrname, prev_bin_end, start_pos, 0.0)
                nuc_state_bdg.add_loc(chrname, prev_bin_end, start_pos, 0.0)
                bg_state_bdg.add_loc(chrname, prev_bin_end, start_pos, 1.0)

        open_state_bdg.add_loc(chrname, start_pos, end_pos, pp_open)
        nuc_state_bdg.add_loc(chrname, start_pos, end_pos, pp_nuc)
        bg_state_bdg.add_loc(chrname, start_pos, end_pos, pp_bg)
        prev_bin_end = end_pos

    open_state_bdg_file.write_bedGraph("Open States", "Likelihoods of being Open States", trackline=False)
    nuc_state_bdg_file.write_bedGraph("Nucleosomal States", "Likelihoods of being Nucleosomal States", trackline=False)
    bg_state_bdg_file.write_bedGraph("Background States", "Likelihoods of being Background States", trackline=False)
    return


def save_states_bed(states_path, states_bedfile):
    # we do not need to output background state.
    for l_len in range(len(states_path)):
        if states_path[l_len][3] != "bg":
            states_bedfile.write("%s\t" % states_path[l_len][0].decode())
            states_bedfile.write("%d\t%d\t%s\n" % states_path[l_len][1:])
    return


def generate_states_path(predicted_proba_file, binsize, i_open, i_nuc, i_bg):
    # predicted_proba_file is a temporary file
    ret_states_path = []
    labels_list = ["open", "nuc", "bg"]

    prev_chrom_name = None
    prev_bin_end = None
    prev_label = None

    for pp_line in predicted_proba_file:
        pp_data = pp_line.strip().split(b',')

        chrname = pp_data[0]
        end_pos = int(pp_data[1])
        start_pos = end_pos - binsize
        pp_open = float(pp_data[i_open+2])
        pp_nuc = float(pp_data[i_nuc+2])
        pp_bg = float(pp_data[i_bg+2])

        # find the best state as label
        label = labels_list[max((pp_open, 0), (pp_nuc, 1), (pp_bg, 2), key=lambda x: x[0])[1]]

        if chrname != prev_chrom_name:
            # we start a new chromosome
            if start_pos > 0:
                # add the first unannotated region as background
                ret_states_path.append((chrname, 0, start_pos, "bg"))
            ret_states_path.append((chrname, start_pos, end_pos, label))
            prev_label = label
            prev_chrom_name = chrname
        else:
            # now check if the prev_bin_end is start_pos, if not, add a gap of background
            if prev_bin_end < start_pos:
                ret_states_path.append((chrname, prev_bin_end, start_pos, "bg"))
                prev_label = "bg"
            # same label, just extend
            if label == prev_label:
                ret_states_path[-1] = (ret_states_path[-1][0], ret_states_path[-1][1], end_pos, label)
            else:
                ret_states_path.append((chrname, start_pos, end_pos, label))
                prev_label = label

        prev_bin_end = end_pos
    return ret_states_path


def save_accessible_regions(states_path, accessible_region_file, openregion_minlen, bdgscore):
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
        if (states_path[i][3] == 'nuc' and
            states_path[i+1][3] == 'open' and
            states_path[i+2][3] == 'nuc' and
            states_path[i][2] == states_path[i+1][1] and
            states_path[i+1][2] == states_path[i+2][1] and
            states_path[i+2][2] - states_path[i][1] > openregion_minlen):  # require nuc-open-nuc entire region start/endpos > openregion_minlen
            accessible_regions = add_regions(i, accessible_regions)

    # remove 'nuc' regions:
    accessible_regions = [tup for tup in accessible_regions if tup[3] != 'nuc']

    # Generate broadpeak object
    openpeak = PeakIO()
    for region in accessible_regions[:-1]:
        openpeak.add(chromosome=region[0], start=region[1], end=region[2])

    # refine peak summit and score using bedGraphTrackI with scores
    openpeak = bdgscore.refine_peaks(openpeak)
    openpeak.write_to_narrowPeak(accessible_region_file)
    return
