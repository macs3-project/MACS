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

    # read in genome size stats ... is this necessary?
    gReader = GenomeFileReader(genomeFile)
    genomeStats = gReader.getMap

    # read in blacklisted if option entered
    if blacklist:
        black = bedFileReader(blacklist)
        info("Read in blacklisted...")
        #bedFileReader, getData are used in java file

    # set fragment length distribution parameters use input values to set initial values, if provided. else use defaults
    #  do we need to include these, since defaults are entered in macs3 file?
    mode = [0.5, 2, 2, 2]
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

        puller = pullLargeLengths(bamfile, index, minMapQ, genomeStats, em_means)
        lengths = puller.getSampledLengths
        weights = puller.getWeights
        # Perform EM training:
        em = HMMR_EM(weights, em_means, em_stddev, lengths)
        info("Fragment Expectation Maximum performed on the fragment distribution...")
    else:
        info("EM training NOT performed on the fragment distribution...")
    

	# Generate genome wide pileup of read coverage and simultaneously calculate mean and standard deviation to calculate fold change and zscore. Convert pileup to bedgraph for easier storage. Calculate the pileup in chunks for memory efficiency
    # NOTE: this method does not currently exist. Instead we use an inputed big wig file to find fold change training regions and zscored masked regions
    pileupData = pileup(SplitBed(genomeStats, window), bamfile, index, keep_dups)
    fc = bedGraphmath(pileupData.getBedgraph)
    cpmScale = pileupData.getCPMScale/1000000
    info("Scaling Factor = %s", cpmScale)
    genomeMean = fc.getMean
    genomeStd = fc.getStd
    
    train = MergeBed(fc.getBetweenRanges(upper, lower))
    
    if len(train) < maxTrain:
        maxTrain = len(train)

    # Shuffle training list before choosing (add randomTrainSeed here)
    newTrain = []
    for i in range(0, maxTrain):
        newTrain.append(train[i])
    
    train = newTrain
    train = ExtendBed(train)
    exclude = MergeBed(fc.getAboveZscore(zscore))
    addback = exclude

    if blacklist:
        exclude.addAll(black)
        exclude = MergeBed(exclude)

    if print_ex:
        ex = "excluded regions to print"

    for i in range(0, len(train)):
        counter = 0
        node1 = train[i]
        for j in range(0, len(exclude)):
            node2 = exclude[j]
            if SubtractedBed.overlap(node1, node2):
                counter += 1
        if counter == 0:
            newTrain.append(node1)

    if train_regions:
        trainReader = bedFileReader(train_regions)
        train = trainReader()
    info("Training regions found and Zscore regions for exclusions found...")

    # Create the fragment pileup tracks using the training set and the fragment distribution parameters
    # Use input model if available
    if model:
        hmm = HMMBinaryReader(model)
        info("Binary model file used (generated from previous HMMR run), training skipped... ")
    else:
        # print output to _training.bed
        gen = FragPileupGen(bamfile, index, train, mode, em_means, em_stddev, minMapQ, keep_dups, cpmScale)
        holder = TrackHolder(gen, trim)
        info("Training Fragment Pileup completed...")
    
        # Create the inital model using KMeans and then refine it using Baum-Welch
        kmeans = KMeansToHMM(holder,states)
        info("Kmeans Model: ...")
        hmm = BaumWelch(kmeans, holder)
    # Identify peak state as the state with the highest short read signal.
    # Identify flanking nucleosome state as the state witht the second highest mono read signal.
    peak = -1
    pmax = 0
    for i in hmm.nbStates:
        pdf = OpdfMultiGaussian(hmm.getOpdf(i))
        sh = pdf.mean()
        if sh > pmax:
            peak = i
            pmax = sh

    pmax = 0
    for i in hmm.nbStates:
        if i != peak:
            pdf = OpdfMultiGaussian(hmm.getOpdf(i))
            mono = pdf.mean()
            if mono > pmax:
                pmax = mono

    # Output binary model file
    outModel = open( modelOnly, "w" )
    HmmBinaryWriter.write(outModel, hmm)  
    outModel = outModel.write(hmm)
    outModel.close()
    info("Model created and refined. See %s.model", output)
    info("Model: \n %s", hmm)

    # Stop program if only model is desired  
    if modelOnly:
        info("HMM Model generated. This can be later applied with `--model`")
        sys.exit(0)

    # Split the genome file into smaller 25MB chunks 
	# Can also split into whatever sized chunks the users prefers
	# May be necessary to split into smaller chunks for machines with less memory
    split = SplitBed(genomeStats, window)

    # Subtract excluded regions from the split genome for Viterbi
    vitBed = SubtractBed(split, exclude)
    info("Genome split and subtracted masked regions...")

    # Run Viterbi on the whole genome 
    if printHMMRTracks: # is print_train = printHMMRTracks?
        pass # should NFR, mono, di, tri files be included in macs3 file?

    genomeAnnotation = []
    for i in range(0, len(vitBed)):
        tempBed = []
        if vitBed[i] >= 10:
            tempBed.append(vitBed[i])
            vGen = FragPileupGen(bamfile, index, tempBed, mode, em_means, em_stddev, minMapQ, keep_dups, cpmScale)
            vHolder = TrackHolder(vGen, trim)
            if printHMMRTracks:
                tracks = HMMRTracksToBedgraph(vHolder, vitBed)
                pass #nfr, mono, di, tri

            HMM = RobustHMM(vHolder, hmm)
            hmm_states = HMM.getStates # is this the same as "states"
            start = vitBed[i].getStart
            remainder = vitBed[i].getLength % 10
            pile = []
            for j in range(0, len(hmm_states)-1):
                pNode = PileupNode2(start+(j*10), hmm_states, vitBed[i])
                pile.append(pNode)
            pNode = PileupNode2((start+((j*10)-remainder)),hmm_states, vitBed[i])
            pile.append(pNode)
            genomeAnnotation.append(PileupToBedGraph(pile).getBedGraph)

            if (i % 50 == 0 or i == len(vitBed)-1):
                info("%d round viterbi done", i)
    if printHMMRTracks:
        pass # close write files

    # Report the final results as peaks, bedgraphs and summits, if desired
    if storebdg:
        pass # print .bedgraph
    if store_peaks:
        pass # print _peaks.gappedPeak, _summits.bed
    
    bdg = fc.getMappedBedgraph

    hmmrBdg = bedGraphMath.toMap(genomeAnnotation)

    # not exactly sure what is happening in java loop here (line 575+):
    for chr in hmmrBdg.keys():
        hmmr = hmmrBdg[chr]
        signal = bdg[chr]
        if signal:
            hmmr.sort() # TagNode.basepairComparator
            signal.sort() # TagNode.basepairComparator
            for i in range(0, len(hmmr)):
                TagNode.temp = hmmr[i]

            # Execute the scoring commands if the state is a peak or if bedgraph scoring is on
                if temp.getScore2 == peak or store_bgscore:
                    overlaps = []
                    hasHadOverlap = False
                    for a in range(0, len(signal)):
                        if SubtractBed.overlap(temp, signal[a]):
                            overlaps.append(signal[a])
                            hasHadOverlap = True
                        else:
                            if hasHadOverlap:
                                index = a
                                break
                    scores = bedGraphMath(temp, overlaps)
                    max_score = scores.getMax
                    mean_score = scores.getMean
                    median_score = scores.getMedian
                    z_score = (scores.getMean - genomeMean)/genomeStd
                    foldchange = scores.getMean/genomeMean
                    
                    if score == "ave":
                        pass
                    elif score == "fc":
                        pass
                    elif score == "zscore":
                        pass
                    elif score == "med":
                        pass
                    elif score == "all":
                        pass
                    else:
                        pass
                if temp.getScore2 == peak:
                    temp = bedGraphMath.setSmooth(temp, overlaps)
                    
                    if i > 0:
                        pass
                    else:
                        pass
                    if i < len(hmmr):
                        pass
            # report the bedgraph if desired
            if storebdg:
                if store_bgscore:
                    pass
                else:
                    pass

            # report the peaks and summits if desired
            if store_peaks and temp.getScore2 == peak and temp.getLength >= minlen and temp.getScore3 >= threshold:
                if temp.getSummit:
                    pass
            
    if storebdg:
        pass #close file
    if store_peaks:
        count = 0
        for i in range(0, len(addback)):
            pass

    # print time to run


    print ( options )
    return
