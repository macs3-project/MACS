# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-04-15 00:37:37 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
from math import sqrt
import logging
debug   = logging.debug
info    = logging.info
# ------------------------------------
# Other modules
# ------------------------------------

import numpy as np
cimport numpy as np
from cpython cimport bool

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.Prob import pnorm2
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.Region import Regions

# ------------------------------------
# Misc functions
# ------------------------------------

cdef inline float get_weighted_density( int x, float m, float v, w ):
    """Description:
    
    parameters:
      1. x: the observed value
      2. m: the mean of gaussian
      3. v: the variance of the gaussian
      4. w: the weight
    return value:
    """
    return w * pnorm2( float(x), m, v )

# ------------------------------------
# Classes
# ------------------------------------
# ------------------------------------
# public functions
# ------------------------------------

cpdef list generate_weight_mapping( list fraglen_list, list means, list stddevs ):
    """Generate weights for each fragment length in short, mono, di, and tri-signals track

    return: list of four dictionaries, with key as fraglen and value as the weight.
            ret[0] -- dictionary for short
            ret[1] -- dictionary for mono
            ret[2] -- dictionary for di
            ret[3] -- dictionary for tri
    """
    cdef:
        list ret_mapping
        list variances
        int l
        float m_s, m_m, m_d, m_t
        float v_s, v_m, v_d, v_t
        float p_s, p_m, p_d, p_t
        float w_s, w_m, w_d, w_t
        float s
        int i, j
    assert len(means) == 4
    assert len(stddevs) == 4
    [m_s, m_m, m_d, m_t] = means
    [v_s, v_m, v_d, v_t] = [ x**2 for x in stddevs ]
    ret_mapping = [ {}, {}, {}, {} ]
    for i in range( len(fraglen_list) ):
        l = fraglen_list[ i ]
        p_s = pnorm2( float(l), m_s, v_s )
        p_m = pnorm2( float(l), m_m, v_m )
        p_d = pnorm2( float(l), m_d, v_d )
        p_t = pnorm2( float(l), m_t, v_t )
        s = p_s + p_m + p_d + p_t
        w_s = p_s / s
        w_m = p_m / s
        w_d = p_d / s
        w_t = p_t / s
        ret_mapping[ 0 ][ l ] = w_s
        ret_mapping[ 1 ][ l ] = w_m
        ret_mapping[ 2 ][ l ] = w_d
        ret_mapping[ 3 ][ l ] = w_t
    return ret_mapping

cpdef list generate_digested_signals( object petrack, list weight_mapping ):
    """Generate digested pileup signals (four tracks) using weight mapping 

    return: list of four signals in dictionary, with key as chromosome name and value as a p-v array.
            ret[0] -- dictionary for short
            ret[1] -- dictionary for mono
            ret[2] -- dictionary for di
            ret[3] -- dictionary for tri
    """
    cdef:
        list ret_digested_signals
        list ret_bedgraphs
        object bdg
        int i
        dict certain_signals
        np.ndarray pv
        bytes chrom
    ret_digested_signals = petrack.pileup_bdg_hmmr( weight_mapping )
    ret_bedgraphs = []
    for i in range( 4 ):                  #yes I hardcoded 4!
        certain_signals = ret_digested_signals[ i ]
        bdg = bedGraphTrackI()
        for chrom in certain_signals.keys():
            bdg.add_chrom_data_hmmr_PV( chrom, certain_signals[ chrom ] )
        ret_bedgraphs.append( bdg )
    return ret_bedgraphs

cpdef list extract_signals_from_regions( list signals, object peaks, int binsize = 10, flanking = 500 ):
    # we will take regions in peaks, create a bedGraphTrackI with
    # binned regions in peaks, then let them overlap with signals to
    # create a list (4) of value arrays.  flanking: flanking regions
    # beyond the peak region, we want to include some background
    # regions.
    cdef:
        list extracted_data
        object signaltrack
        object peaksbdg
        bytes chrom
        int i, s, e, tmp_s, tmp_e, tmp_n, n
        list ps
        object p
        list ret_training_data, ret_training_lengths
        object regions

    # peaks.sort()
    # fhd = open("original_peaks.bed","w")
    # peaks.write_to_bed( fhd )
    # fhd.close()
    
    regions = Regions()
    regions.init_from_PeakIO( peaks )

    # fhd = open("before.bed","w")
    # regions.write_to_bed( fhd )
    # fhd.close()

    regions.expand( flanking )
    regions.merge_overlap()

    # fhd = open("after.bed","w")
    # regions.write_to_bed( fhd )
    # fhd.close()

    peaksbdg = bedGraphTrackI(baseline_value=0)

    n = 0
    # here we convert peaks from a PeakIO to BedGraph object with a
    # given binsize.
    for chrom in regions.get_chr_names():
        tmp_p = 0
        ps = regions[ chrom ]
        for i in range( len( ps ) ):
            p = ps[ i ]
            s = p[ 0 ]
            e = p[ 1 ]
            # make bins, no need to be too accurate...
            s = s//binsize*binsize
            e = e//binsize*binsize
            tmp_n = int(( e - s )/binsize)
            for r in range( s, e, binsize ):
                tmp_s = r
                tmp_e = r + binsize
                if tmp_s > tmp_p:
                    peaksbdg.add_loc_wo_merge( chrom, tmp_p, tmp_s, 0 )
                peaksbdg.add_loc_wo_merge( chrom, tmp_s, tmp_e, tmp_n )
                n += 1
                tmp_p = tmp_e
    # we do not merge regions in peaksbdg object so each bin will be seperated.
    debug( f"added {n} bins" )
    #print( peaksbdg.summary() )
    #print( peaksbdg.total() )

    # bfhd = open("a.bdg","w")
    # peaksbdg.write_bedGraph( bfhd, "peaksbdg", "peaksbdg" )
    # bfhd.close()
    # now, let's overlap
    extracted_data = []
    extracted_len = []
    for signaltrack in signals:
        # signaltrack is bedGraphTrackI object
        [ values, lengths ] = signaltrack.extract_value_hmmr( peaksbdg )
        # we only need values
        #print( regions )
        extracted_data.append( values )
        extracted_len.append( lengths )
    ret_training_data = []
    ret_training_lengths = []
    c = 0
    nn = len( extracted_data[0] )
    nnn =len( extracted_len[0] )
    debug( f"{n} bins, {nn}, {nnn}" )
    for i in range( nn ):
        c += 1
        ret_training_data.append(
            [ max( 0.0001, abs(round(extracted_data[0][i], 4))),
              max( 0.0001, abs(round(extracted_data[1][i], 4))),
              max( 0.0001, abs(round(extracted_data[2][i], 4))),
              max( 0.0001, abs(round(extracted_data[3][i], 4))) ] )
        if c == extracted_len[0][i]:
            ret_training_lengths.append( extracted_len[0][i] )
            c = 0
    return [ ret_training_data, ret_training_lengths ]


