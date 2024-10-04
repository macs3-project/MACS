# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-04 10:25:29 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
import cython
from MACS3.Signal.Prob import pnorm2
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.Region import Regions
from MACS3.Utilities.Logger import logging

logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.inline
@cython.cfunc
def get_weighted_density(x: cython.int, m: cython.float, v:
                         cython.float, w: cython.float) -> cython.float:
    """Description:

    parameters:
      1. x: the observed value
      2. m: the mean of gaussian
      3. v: the variance of the gaussian
      4. w: the weight
    return value:
    """
    return w * pnorm2(float(x), m, v)


# ------------------------------------
# Classes
# ------------------------------------
# ------------------------------------
# public functions
# ------------------------------------


@cython.ccall
def generate_weight_mapping(fraglen_list: list, means: list, stddevs:
                            list, min_frag_p: cython.float = 0.001) -> list:
    """Generate weights for each fragment length in short, mono, di,
    and tri-signals track

    return: list of four dictionaries, with key as fraglen and value
    as the weight.

    ret[0] -- dictionary for short
    ret[1] -- dictionary for mono
    ret[2] -- dictionary for di
    ret[3] -- dictionary for tri

    """

    ret_mapping: list
    fl: cython.int
    m_s: cython.float
    m_m: cython.float
    m_d: cython.float
    m_t: cython.float
    v_s: cython.float
    v_m: cython.float
    v_d: cython.float
    v_t: cython.float
    p_s: cython.float
    p_m: cython.float
    p_d: cython.float
    p_t: cython.float
    s: cython.float
    i: cython.int

    assert len(means) == 4
    assert len(stddevs) == 4
    [m_s, m_m, m_d, m_t] = means
    [v_s, v_m, v_d, v_t] = [x**2 for x in stddevs]
    ret_mapping = [{}, {}, {}, {}]
    for i in range(len(fraglen_list)):
        fl = fraglen_list[i]
        p_s = pnorm2(float(fl), m_s, v_s)
        p_m = pnorm2(float(fl), m_m, v_m)
        p_d = pnorm2(float(fl), m_d, v_d)
        p_t = pnorm2(float(fl), m_t, v_t)
        s = p_s + p_m + p_d + p_t
        if p_s < min_frag_p and p_m < min_frag_p and p_d < min_frag_p and p_t < min_frag_p:
            # we exclude the fragment which can't be assigned to
            # short, mono, di-nuc, and tri-nuc (likelihood <
            # min_frag_p, default:0.001) Normally this fragment is too
            # large. We exclude these fragment by setting all weights
            # to zero.
            debug(f"The fragment length {fl} can't be assigned to either distribution so will be excluded!")
            ret_mapping[0][fl] = 0
            ret_mapping[1][fl] = 0
            ret_mapping[2][fl] = 0
            ret_mapping[3][fl] = 0
            continue
        ret_mapping[0][fl] = p_s / s
        ret_mapping[1][fl] = p_m / s
        ret_mapping[2][fl] = p_d / s
        ret_mapping[3][fl] = p_t / s
    return ret_mapping


@cython.ccall
def generate_digested_signals(petrack, weight_mapping: list) -> list:
    """Generate digested pileup signals (four tracks) using weight mapping

    return: list of four signals in dictionary, with key as chromosome name and value as a p-v array.
            ret[0] -- dictionary for short
            ret[1] -- dictionary for mono
            ret[2] -- dictionary for di
            ret[3] -- dictionary for tri
    """
    ret_digested_signals: list
    ret_bedgraphs: list
    bdg: object
    i: int
    certain_signals: dict
    chrom: bytes

    ret_digested_signals = petrack.pileup_bdg_hmmr(weight_mapping)
    ret_bedgraphs = []
    for i in range(4):          # yes I hardcoded 4!
        certain_signals = ret_digested_signals[i]
        bdg = bedGraphTrackI()
        for chrom in sorted(certain_signals.keys()):
            bdg.add_chrom_data_hmmr_PV(chrom, certain_signals[chrom])
        ret_bedgraphs.append(bdg)
    return ret_bedgraphs


@cython.ccall
def extract_signals_from_regions(signals: list, regions, binsize:
                                 cython.int = 10, hmm_type: str = 'gaussian') -> list:
    # we will take regions in peaks, create a bedGraphTrackI with
    # binned regions in peaks, then let them overlap with signals to
    # create a list (4) of value arrays.
    #
    extracted_data: list
    extracted_len: list
    extracted_positions: list
    signaltrack: object
    regionsbdg: object
    i: cython.int
    c: cython.int
    counter: cython.int
    prev_c: cython.int
    ret_training_data: list
    ret_training_lengths: list
    ret_training_bins: list

    regionsbdg = _make_bdg_of_bins_from_regions(regions, binsize)
    debug('#      extract_signals_from_regions: regionsbdg completed')
    # now, let's overlap
    extracted_positions = []
    extracted_data = []
    extracted_len = []
    for signaltrack in signals: # four signal tracks
        # signaltrack is bedGraphTrackI object
        [positions, values, lengths] = signaltrack.extract_value_hmmr(regionsbdg)
        extracted_positions.append(positions)
        extracted_data.append(values)
        extracted_len.append(lengths)
    positions = []
    values = []
    lengths = []
    debug('#      extract_signals_from_regions: extracted positions, data, len')
    ret_training_bins = []
    ret_training_data = []
    ret_training_lengths = []

    nn = len(extracted_data[0])
    assert nn > 0
    assert nn == len(extracted_data[1])
    assert nn == len(extracted_data[2])
    assert nn == len(extracted_data[3])
    counter = 0
    prev_c = extracted_len[0][0]
    c = 0
    if hmm_type == "gaussian":
        for i in range(nn):
            ret_training_bins.append(extracted_positions[0][i])
            ret_training_data.append(
                [max(0.0001, extracted_data[0][i]),
                 max(0.0001, extracted_data[1][i]),
                 max(0.0001, extracted_data[2][i]),
                 max(0.0001, extracted_data[3][i])])
            c = extracted_len[0][i]
            if counter != 0 and c != prev_c:
                ret_training_lengths.append(counter)
                counter = 0
            prev_c = c
            counter += 1
        debug('#      extract_signals_from_regions: ret_training bins, data, lengths - gaussian')
    # poisson can only take int values as input
    if hmm_type == "poisson":
        for i in range(nn):
            ret_training_bins.append(extracted_positions[0][i])
            ret_training_data.append(
                [int(max(0.0001, extracted_data[0][i])),
                 int(max(0.0001, extracted_data[1][i])),
                 int(max(0.0001, extracted_data[2][i])),
                 int(max(0.0001, extracted_data[3][i]))])
            c = extracted_len[0][i]
            if counter != 0 and c != prev_c:
                ret_training_lengths.append(counter)
                counter = 0
            prev_c = c
            counter += 1
        debug('#      extract_signals_from_regions: ret_training bins, data, lengths - poisson')
    # last region
    ret_training_lengths.append(counter)
    assert sum(ret_training_lengths) == len(ret_training_data)
    assert len(ret_training_bins) == len(ret_training_data)
    return [ret_training_bins, ret_training_data, ret_training_lengths]


@cython.cfunc
def _make_bdg_of_bins_from_regions(regions, binsize: cython.int):
    # this function will return a BedGraphTrackI object
    regionsbdg: object
    n: cython.long
    chrom: bytes
    ps: list
    s: cython.int
    e: cython.int
    tmp_p: cython.int
    mark_bin: cython.int
    i: cython.int
    r: cython.int

    assert isinstance(regions, Regions)

    regionsbdg = bedGraphTrackI(baseline_value=-100)

    n = 0
    # here we convert peaks from a PeakIO to BedGraph object with a
    # given binsize.
    mark_bin = 1  # this is to mark the continuous bins in the same region, it will increase by one while moving to the next region
    for chrom in sorted(regions.get_chr_names()):
        tmp_p = 0  # this is to make gap in bedgraph for not covered regions.
        ps = regions[chrom]
        for i in range(len(ps)):
            # for each region
            s = ps[i][0]
            e = ps[i][1]
            # make bins, no need to be too accurate...
            s = s//binsize*binsize
            e = e//binsize*binsize
            # tmp_n = int((e - s)/binsize)
            for r in range(s, e, binsize):
                tmp_s = r
                tmp_e = r + binsize
                if tmp_s > tmp_p:
                    regionsbdg.add_loc_wo_merge(chrom, tmp_p, tmp_s, 0)  # the gap
                regionsbdg.add_loc_wo_merge(chrom, tmp_s, tmp_e, mark_bin)  # the value we put in the bin bedgraph is the number of bins in this region
                n += 1
                tmp_p = tmp_e
            # end of region, we change the mark_bin
            mark_bin += 1
    # we do not merge regions in regionsbdg object so each bin will be separated.
    debug(f"added {n} bins")
    return regionsbdg
