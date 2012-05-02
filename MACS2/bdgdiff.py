# Time-stamp: <2012-04-30 17:42:07 Tao Liu>

"""Description: Naive call differential peaks from 4 bedGraph tracks for scores.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import logging
from MACS2.IO.cBedGraphIO import bedGraphIO,genericBedIO
from MACS2.IO.cPeakIO import Region
from MACS2.IO.cCompositeScoreTrack import *
from MACS2.cStat import *
#from MACS2.data import PreCompiledGFold as PCGF
from MACS2.data import FakePreCompiledGFold as PCGF
from math import log

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

# ------------------------------------
# Misc functions
# ------------------------------------
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( options ):
    options.do_MCMC = True
    
    # load precompiled matrix
    gfolds_c = PCGF(options.cutoff)

    info("Read peak files...")
    info("Peak of condition 1 treatment...")

    t1_peakio = genericBedIO(options.peak1)
    t1_peak = t1_peakio.build_bedtrack()

    info("Peak of condition 2 treatment...")

    t2_peakio = genericBedIO(options.peak2)
    t2_peak = t2_peakio.build_bedtrack()

    # get union peak regions
    union_peak = t1_peak.overlie(t2_peak)

    info("Read and build bedGraph...")
    info("Pileup of condition 1 treatment...")
    t1_bio = bedGraphIO(options.t1bdg)
    t1_btrack = t1_bio.build_bdgtrack(baseline_value=0)

    info("Pileup of condition 2 treatment...")
    t2_bio = bedGraphIO(options.t2bdg)
    t2_btrack = t2_bio.build_bdgtrack(baseline_value=0)

    # calculate sum of all signals in million
    t1_sum = t1_btrack.summary()[0]
    t2_sum = t2_btrack.summary()[0]    
    n1 = t1_sum/1000000.0               # signal per million
    n2 = t2_sum/1000000.0
    offset = -log(n1,2)+log(n2,2)
    info("t1 sum: %.1f, t2 sum: %.1f, Offset is %.2f" % (t1_sum,t2_sum,offset))

    # combine two tracks
    info("Combine tracks...")
    comb_track = t1_btrack.make_scoreTrack_for_macs2diff(t2_btrack)

    info("Extract average values in union regions...")
    data_in_union = comb_track.extract_average(union_peak) # ([id,...],[count1,...],[count2,...])

    # if n1 > n2:
    #     r1 = n1/n2
    #     r2 = 1
    # else:
    #     r1 = 1
    #     r2 = n2/n1
    for i in xrange(len(data_in_union[0])):
        data_in_union[1][i] = int(data_in_union[1][i]) # actual values are Pileup Per Peak Per Million reads (PPPPM)
        data_in_union[2][i] = int(data_in_union[2][i])

    #info("Convert gfold...")
    info( "Calculate gfold ..." )
    gfolds = convert_gfold(data_in_union, gfolds_c, offset=offset, cutoff=options.cutoff,mcmc=options.do_MCMC)
    
    # sort by gfold
    gfolds.sort(cmp=lambda x,y:cmp(x[1],y[1]))

    # write differential regions with gfold>0

    info( "Write differential regions to %s ..." % (options.oprefix+"_diff.bed") )
    ofhd = open(options.oprefix+"_diff.bed","w")

    for (rid, gf) in gfolds:
        if gf != 0:
            (chrom,start,end) = rid.split('.')
            ofhd.write( "%s\t%s\t%s\t%s\t%.5f\n" % (chrom,start,end,'.',gf) )

    ofhd.close()

    info( "Write gfold values for each region to %s ..." % (options.oprefix+"_diff.txt") )
    ofhd = open(options.oprefix+"_diff.txt","w")

    gf_dict = dict(gfolds)

    for i in xrange(len(data_in_union[0])):
        gftmp = gf_dict[data_in_union[0][i]]
        tmp1 = data_in_union[1][i]
        tmp2 = data_in_union[2][i]
        ofhd.write("%s\t%.5f\t%.5f\t%.5f\n" % (data_in_union[0][i],tmp1/n1,tmp2/n2,gftmp))

    ofhd.close()

