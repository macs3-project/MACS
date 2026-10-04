#!/usr/bin/env python

"""Module Description: Test functions in HMMR_Signal_Processing that
split ATAC-seq fragments into short, mono-, di- and tri-nucleosomal
signals and extract binned signals for the HMMRATAC Hidden Markov
Model.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import math

import numpy as np
import pytest
from scipy.stats import norm

from MACS3.Signal.HMMR_Signal_Processing import (generate_weight_mapping,
                                                 generate_digested_signals,
                                                 extract_signals_from_regions)
from MACS3.Signal.PairedEndTrack import (PETrackI,
                                         PETrackII)
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.Region import Regions

# ------------------------------------
# Helpers
# ------------------------------------

DEFAULT_MEANS = [50, 200, 400, 600]
DEFAULT_SDS = [20, 20, 20, 20]


def ref_weights(fl, means, sds, min_frag_p=0.001):
    """Reference weights for one fragment length from scipy's normal pdf.

    A length whose density is below ``min_frag_p`` under all four
    distributions gets weight 0 everywhere; otherwise the weight of
    each distribution is its density divided by the sum of the four.
    """
    p = [norm.pdf(fl, m, s) for m, s in zip(means, sds)]
    if all(x < min_frag_p for x in p):
        return [0.0, 0.0, 0.0, 0.0]
    s = sum(p)
    return [x / s for x in p]


def make_petrack(frags):
    """PETrackI from (chrom, left, right) tuples."""
    pe = PETrackI()
    for chrom, l, r in frags:
        pe.add_loc(chrom, l, r)
    pe.finalize()
    return pe


def ref_weighted_pileup(frags, weights):
    """Reference pileup for one chromosome as bedGraph (ends, values).

    ``frags`` are (left, right, count); each covers [left, right) with
    ``weights[right - left] * count``. The track starts at 0, ends at
    the largest right end, and neighbours with equal values merge.
    """
    end = max(r for _, r, _ in frags)
    cov = np.zeros(end, dtype=np.float64)
    for l, r, c in frags:
        cov[l:r] += weights[r - l] * c
    pos, val = [], []
    for i in range(end):
        if val and cov[i] == val[-1]:
            pos[-1] = i + 1
        else:
            pos.append(i + 1)
            val.append(float(cov[i]))
    return pos, val


def track_from_segments(segments):
    """bedGraphTrackI from {chrom: [(start, end, value), ...]} (contiguous)."""
    bdg = bedGraphTrackI()
    for chrom, segs in segments.items():
        for s, e, v in segs:
            bdg.add_loc(chrom, s, e, v)
    return bdg


def regions_from(locs):
    """Regions from (chrom, start, end) tuples, sorted."""
    r = Regions()
    for chrom, s, e in locs:
        r.add_loc(chrom, s, e)
    r.sort()
    return r


# ------------------------------------
# generate_weight_mapping
# ------------------------------------

# lengths chosen away from the min_frag_p boundary: 100, 125, 300 and
# 500 are more than 2.45 sd from every mean (density < 0.001 for
# sd=20), the others are close to at least one mean.
FRAGLENS = [30, 50, 64, 100, 125, 150, 200, 250, 300, 350, 400, 450,
            500, 550, 600, 640, 650, 1000]


def test_generate_weight_mapping_structure():
    ret = generate_weight_mapping(FRAGLENS, DEFAULT_MEANS, DEFAULT_SDS)
    assert isinstance(ret, list)
    assert len(ret) == 4
    for d in ret:
        assert isinstance(d, dict)
        assert sorted(d.keys()) == FRAGLENS


@pytest.mark.parametrize("fl", FRAGLENS)
def test_generate_weight_mapping_default_against_normal_pdf(fl):
    ret = generate_weight_mapping(FRAGLENS, DEFAULT_MEANS, DEFAULT_SDS)
    expected = ref_weights(fl, DEFAULT_MEANS, DEFAULT_SDS)
    got = [ret[k][fl] for k in range(4)]
    # float32 arithmetic inside pnorm2: relative 1e-5 is ample
    assert got == pytest.approx(expected, rel=1e-5, abs=1e-12)


@pytest.mark.parametrize("fl, excluded", [
    # density of N(600, 20) at 640 (z=2) is 0.0027 > 0.001: kept
    (640, False),
    # at 650 (z=2.5) it is 0.00088 < 0.001: excluded
    (650, True),
    # 100 is 2.5 sd from the short mean 50: excluded
    (100, True),
    # 125 is 3.75 sd from both 50 and 200: excluded
    (125, True),
    (1000, True),
    (200, False),
])
def test_generate_weight_mapping_min_frag_p_exclusion(fl, excluded):
    ret = generate_weight_mapping([fl], DEFAULT_MEANS, DEFAULT_SDS)
    got = [ret[k][fl] for k in range(4)]
    if excluded:
        assert got == [0, 0, 0, 0]
    else:
        assert sum(got) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("min_frag_p", [1e-6, 1e-4, 0.001, 0.005, 0.019])
def test_generate_weight_mapping_min_frag_p_values(min_frag_p):
    fls = list(range(0, 1001, 5))
    ret = generate_weight_mapping(fls, DEFAULT_MEANS, DEFAULT_SDS,
                                  min_frag_p=min_frag_p)
    for fl in fls:
        expected = ref_weights(fl, DEFAULT_MEANS, DEFAULT_SDS, min_frag_p)
        got = [ret[k][fl] for k in range(4)]
        if expected == [0.0] * 4:
            assert got == [0, 0, 0, 0], fl
        else:
            assert got == pytest.approx(expected, rel=1e-5, abs=1e-12), fl


def test_generate_weight_mapping_larger_min_frag_p_excludes_more():
    fls = list(range(0, 1001))

    def n_kept(p):
        ret = generate_weight_mapping(fls, DEFAULT_MEANS, DEFAULT_SDS,
                                      min_frag_p=p)
        return sum(1 for fl in fls if any(ret[k][fl] for k in range(4)))
    counts = [n_kept(p) for p in (1e-6, 1e-4, 1e-3, 1e-2)]
    assert counts == sorted(counts, reverse=True)
    assert counts[0] > counts[-1]


def test_generate_weight_mapping_min_frag_p_above_max_density():
    # the peak density of N(m, 20) is 1/(20*sqrt(2*pi)) = 0.01995, so a
    # cutoff of 0.02 excludes every fragment length
    fls = [50, 200, 400, 600]
    ret = generate_weight_mapping(fls, DEFAULT_MEANS, DEFAULT_SDS,
                                  min_frag_p=0.02)
    for d in ret:
        assert d == {50: 0, 200: 0, 400: 0, 600: 0}


def test_generate_weight_mapping_custom_means_sds():
    means = [40, 180, 360, 540]
    sds = [10, 25, 30, 35]
    fls = list(range(20, 700, 7))
    ret = generate_weight_mapping(fls, means, sds, min_frag_p=1e-5)
    for fl in fls:
        expected = ref_weights(fl, means, sds, 1e-5)
        got = [ret[k][fl] for k in range(4)]
        assert got == pytest.approx(expected, rel=1e-5, abs=1e-12), fl


def test_generate_weight_mapping_equal_distributions_split_evenly():
    # short and mono identical: a length at their mean is split 50/50
    ret = generate_weight_mapping([100], [100, 100, 400, 600],
                                  [20, 20, 20, 20])
    assert ret[0][100] == pytest.approx(0.5, abs=1e-6)
    assert ret[1][100] == pytest.approx(0.5, abs=1e-6)
    assert ret[2][100] == pytest.approx(0.0, abs=1e-12)
    assert ret[3][100] == pytest.approx(0.0, abs=1e-12)


def test_generate_weight_mapping_negative_stddev_same_as_positive():
    # only the variance (sd squared) is used
    fls = [50, 120, 200, 380]
    a = generate_weight_mapping(fls, DEFAULT_MEANS, [20, 20, 20, 20])
    b = generate_weight_mapping(fls, DEFAULT_MEANS, [-20, 20, -20, 20])
    assert a == b


def test_generate_weight_mapping_empty():
    assert generate_weight_mapping([], DEFAULT_MEANS, DEFAULT_SDS) == \
        [{}, {}, {}, {}]


def test_generate_weight_mapping_single_length():
    ret = generate_weight_mapping([200], DEFAULT_MEANS, DEFAULT_SDS)
    assert [d[200] for d in ret] == pytest.approx(
        ref_weights(200, DEFAULT_MEANS, DEFAULT_SDS), rel=1e-5, abs=1e-12)


def test_generate_weight_mapping_duplicate_lengths():
    ret = generate_weight_mapping([200, 200, 50], DEFAULT_MEANS, DEFAULT_SDS)
    assert sorted(ret[1].keys()) == [50, 200]


def test_generate_weight_mapping_int32_max_length():
    big = 2**31 - 1
    ret = generate_weight_mapping([big], DEFAULT_MEANS, DEFAULT_SDS)
    assert [d[big] for d in ret] == [0, 0, 0, 0]


def test_generate_weight_mapping_length_over_int32_raises():
    with pytest.raises(OverflowError):
        generate_weight_mapping([2**31], DEFAULT_MEANS, DEFAULT_SDS)


@pytest.mark.parametrize("means, sds", [
    ([50, 200, 400], [20, 20, 20, 20]),
    ([50, 200, 400, 600], [20, 20, 20]),
    ([50, 200, 400, 600, 800], [20, 20, 20, 20, 20]),
])
def test_generate_weight_mapping_wrong_number_of_parameters(means, sds):
    with pytest.raises(AssertionError):
        generate_weight_mapping([100], means, sds)


def test_generate_weight_mapping_requires_list():
    with pytest.raises(TypeError):
        generate_weight_mapping((100, 200), DEFAULT_MEANS, DEFAULT_SDS)


# ------------------------------------
# generate_digested_signals
# ------------------------------------

# fragment lengths 50, 200 and 400 with weights that are exact in
# float32, so the reference pileup is exact
HAND_MAPPING = [{50: 1.0, 200: 0.0, 400: 0.0},
                {50: 0.0, 200: 1.0, 400: 0.25},
                {50: 0.0, 200: 0.0, 400: 0.5},
                {50: 0.0, 200: 0.0, 400: 0.25}]

HAND_FRAGS = [(b"chr1", 0, 50), (b"chr1", 10, 60), (b"chr1", 100, 300),
              (b"chr1", 150, 550), (b"chr2", 5, 205), (b"chr2", 30, 80)]


def test_generate_digested_signals_returns_four_bedgraphs():
    ret = generate_digested_signals(make_petrack(HAND_FRAGS), HAND_MAPPING)
    assert len(ret) == 4
    for bdg in ret:
        assert isinstance(bdg, bedGraphTrackI)
        assert bdg.get_chr_names() == {b"chr1", b"chr2"}


@pytest.mark.parametrize("k", [0, 1, 2, 3])
@pytest.mark.parametrize("chrom", [b"chr1", b"chr2"])
def test_generate_digested_signals_against_reference(k, chrom):
    ret = generate_digested_signals(make_petrack(HAND_FRAGS), HAND_MAPPING)
    p, v = ret[k].get_data_by_chr(chrom)
    exp_p, exp_v = ref_weighted_pileup(
        [(l, r, 1) for c, l, r in HAND_FRAGS if c == chrom], HAND_MAPPING[k])
    assert list(p) == exp_p
    assert list(v) == exp_v


def test_generate_digested_signals_hand_values_short_chr1():
    # short weights: only the two 50 bp fragments [0,50) and [10,60)
    # count: [0,10)=1, [10,50)=2, [50,60)=1, then 0 up to 550, the end
    # of the last chr1 fragment
    ret = generate_digested_signals(make_petrack(HAND_FRAGS), HAND_MAPPING)
    p, v = ret[0].get_data_by_chr(b"chr1")
    assert list(p) == [10, 50, 60, 550]
    assert list(v) == [1.0, 2.0, 1.0, 0.0]


def test_generate_digested_signals_hand_values_mono_chr1():
    # mono: [100,300) weight 1 and [150,550) weight 0.25
    ret = generate_digested_signals(make_petrack(HAND_FRAGS), HAND_MAPPING)
    p, v = ret[1].get_data_by_chr(b"chr1")
    assert list(p) == [100, 150, 300, 550]
    assert list(v) == [0.0, 1.0, 1.25, 0.25]


def test_generate_digested_signals_petrackII_counts():
    # FRAG-style track: a fragment with count 3 contributes 3 times
    pe = PETrackII()
    pe.add_loc(b"chr1", 0, 50, b"AAAA", 3)
    pe.add_loc(b"chr1", 20, 70, b"CCCC", 1)
    pe.add_loc(b"chr1", 100, 300, b"AAAA", 2)
    pe.finalize()
    ret = generate_digested_signals(pe, HAND_MAPPING)
    p, v = ret[0].get_data_by_chr(b"chr1")
    assert (list(p), list(v)) == ref_weighted_pileup(
        [(0, 50, 3), (20, 70, 1), (100, 300, 2)], HAND_MAPPING[0])
    assert list(p) == [20, 50, 70, 300]
    assert list(v) == [3.0, 4.0, 1.0, 0.0]
    p, v = ret[1].get_data_by_chr(b"chr1")
    assert list(p) == [100, 300]
    assert list(v) == [0.0, 2.0]


def test_generate_digested_signals_sum_equals_total_pileup():
    # with weights from generate_weight_mapping (which sum to 1 for every
    # kept length), the four tracks add up to the plain fragment pileup
    rs = np.random.RandomState(7)
    frags = []
    for i in range(200):
        ln = int(rs.choice([50, 60, 190, 200, 210, 400, 600]))
        s = int(rs.randint(0, 5000))
        frags.append((b"chr1", s, s + ln))
    pe = make_petrack(frags)
    fls = sorted(set(r - l for _, l, r in frags))
    mapping = generate_weight_mapping(fls, DEFAULT_MEANS, DEFAULT_SDS)
    ret = generate_digested_signals(pe, mapping)
    end = max(r for _, l, r in frags)
    total = np.zeros(end)
    for _, l, r in frags:
        total[l:r] += 1
    summed = np.zeros(end)
    for bdg in ret:
        p, v = bdg.get_data_by_chr(b"chr1")
        pre = 0
        for pp, vv in zip(p, v):
            summed[pre:pp] += vv
            pre = pp
    assert summed == pytest.approx(total, abs=1e-4)


def test_generate_digested_signals_missing_length_raises():
    with pytest.raises(KeyError):
        generate_digested_signals(make_petrack([(b"chr1", 0, 77)]),
                                  HAND_MAPPING)


def test_generate_digested_signals_single_fragment():
    ret = generate_digested_signals(make_petrack([(b"chr1", 30, 230)]),
                                    HAND_MAPPING)
    p, v = ret[1].get_data_by_chr(b"chr1")
    assert list(p) == [30, 230]
    assert list(v) == [0.0, 1.0]
    p, v = ret[0].get_data_by_chr(b"chr1")
    # weight 0 everywhere: one zero block from 0 to the fragment end
    assert list(p) == [230]
    assert list(v) == [0.0]


# ------------------------------------
# extract_signals_from_regions
# ------------------------------------

# signals on chr1 [0, 200):
#   short: [0,25)=1, [25,60)=2, [60,200)=3
#   mono : [0,200)=5
#   di   : [0,100)=0, [100,200)=0.5
#   tri  : [0,200)=0
SIGNALS = [
    {b"chr1": [(0, 25, 1.0), (25, 60, 2.0), (60, 200, 3.0)]},
    {b"chr1": [(0, 200, 5.0)]},
    {b"chr1": [(0, 100, 0.0), (100, 200, 0.5)]},
    {b"chr1": [(0, 200, 0.0)]},
]


def make_signals(spec=SIGNALS):
    return [track_from_segments(s) for s in spec]


def test_extract_signals_gaussian_hand_example():
    # regions [10,50) and [100,135) with binsize 10: both ends are
    # floored to the bin grid, so the bins are [10,20) .. [40,50) and
    # [100,110) .. [120,130). Each bin is reported by its end position
    # and takes the signal value at its last base (end - 1). Zeros are
    # raised to 0.0001 for the gaussian HMM.
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 10, 50), (b"chr1", 100, 135)]),
        binsize=10)
    assert bins == [(b"chr1", 20), (b"chr1", 30), (b"chr1", 40),
                    (b"chr1", 50), (b"chr1", 110), (b"chr1", 120),
                    (b"chr1", 130)]
    assert data == [[1.0, 5.0, 0.0001, 0.0001],
                    [2.0, 5.0, 0.0001, 0.0001],
                    [2.0, 5.0, 0.0001, 0.0001],
                    [2.0, 5.0, 0.0001, 0.0001],
                    [3.0, 5.0, 0.5, 0.0001],
                    [3.0, 5.0, 0.5, 0.0001],
                    [3.0, 5.0, 0.5, 0.0001]]
    assert lengths == [4, 3]


def test_extract_signals_default_binsize_is_10():
    regions = regions_from([(b"chr1", 10, 50), (b"chr1", 100, 135)])
    a = extract_signals_from_regions(make_signals(), regions)
    b = extract_signals_from_regions(make_signals(), regions, binsize=10)
    assert a == b


def test_extract_signals_poisson_truncates_to_int():
    spec = [
        {b"chr1": [(0, 25, 1.0), (25, 60, 2.5), (60, 200, 3.75)]},
        {b"chr1": [(0, 200, 5.0)]},
        {b"chr1": [(0, 100, 0.0), (100, 200, 0.5)]},
        {b"chr1": [(0, 200, 0.0)]},
    ]
    bins, data, lengths = extract_signals_from_regions(
        make_signals(spec),
        regions_from([(b"chr1", 10, 50), (b"chr1", 100, 135)]),
        binsize=10, hmm_type="poisson")
    # int(max(0.0001, v)): 2.5 -> 2, 3.75 -> 3, 0.5 -> 0, 0 -> 0
    assert data == [[1, 5, 0, 0], [2, 5, 0, 0], [2, 5, 0, 0], [2, 5, 0, 0],
                    [3, 5, 0, 0], [3, 5, 0, 0], [3, 5, 0, 0]]
    assert all(isinstance(x, int) for row in data for x in row)
    assert lengths == [4, 3]
    assert len(bins) == 7


def test_extract_signals_binsize_25_floors_region_ends():
    # [10,50) becomes [0,50): bins [0,25), [25,50); [100,135) becomes
    # [100,125): one bin
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 10, 50), (b"chr1", 100, 135)]),
        binsize=25)
    assert bins == [(b"chr1", 25), (b"chr1", 50), (b"chr1", 125)]
    # value at the last base of each bin: short(24)=1, short(49)=2,
    # short(124)=3; di(124)=0.5
    assert data == [[1.0, 5.0, 0.0001, 0.0001],
                    [2.0, 5.0, 0.0001, 0.0001],
                    [3.0, 5.0, 0.5, 0.0001]]
    assert lengths == [2, 1]


@pytest.mark.parametrize("binsize", [1, 5, 10, 20, 40])
def test_extract_signals_bin_count_and_positions(binsize):
    regions = [(b"chr1", 13, 97), (b"chr1", 120, 190)]
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from(regions), binsize=binsize)
    exp_bins = []
    exp_lengths = []
    for _, s, e in regions:
        s0 = s // binsize * binsize
        e0 = e // binsize * binsize
        ends = list(range(s0 + binsize, e0 + 1, binsize))
        exp_bins.extend((b"chr1", x) for x in ends)
        if ends:
            exp_lengths.append(len(ends))
    assert bins == exp_bins
    assert lengths == exp_lengths
    # values: signal at the last base of each bin
    short = SIGNALS[0][b"chr1"]
    for (c, end), row in zip(bins, data):
        val = [v for s, e, v in short if s <= end - 1 < e][0]
        assert row[0] == val
        assert row[1] == 5.0


def test_extract_signals_region_shorter_than_bin_is_skipped():
    # [12,18) has no full bin after flooring (s=10, e=10)
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 12, 18), (b"chr1", 100, 120)]),
        binsize=10)
    assert bins == [(b"chr1", 110), (b"chr1", 120)]
    assert lengths == [2]


def test_extract_signals_adjacent_regions_are_separate_sequences():
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 0, 30), (b"chr1", 30, 50)]),
        binsize=10)
    assert bins == [(b"chr1", 10), (b"chr1", 20), (b"chr1", 30),
                    (b"chr1", 40), (b"chr1", 50)]
    assert lengths == [3, 2]


def test_extract_signals_bins_beyond_signal_end_are_dropped():
    # the signals end at 200: of the bins of [150, 260) only those
    # ending at 160..200 are returned
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 150, 260)]), binsize=10)
    assert bins == [(b"chr1", x) for x in (160, 170, 180, 190, 200)]
    assert lengths == [5]
    assert data == [[3.0, 5.0, 0.5, 0.0001]] * 5


def test_extract_signals_two_chromosomes_order():
    """Pins the current output.

    The chromosomes come out in reverse lexicographic order (chr2 bins
    before chr1 bins) because extract_value_hmmr pops from the end of a
    sorted list; the order is not derivable from a specification, only
    from the implementation.
    """
    spec = [dict(s) for s in SIGNALS]
    for k, v in enumerate([7.0, 6.0, 2.0, 1.0]):
        spec[k][b"chr2"] = [(0, 100, v)]
    bins, data, lengths = extract_signals_from_regions(
        make_signals(spec),
        regions_from([(b"chr1", 10, 30), (b"chr2", 0, 30)]), binsize=10)
    assert bins == [(b"chr2", 10), (b"chr2", 20), (b"chr2", 30),
                    (b"chr1", 20), (b"chr1", 30)]
    assert data == [[7.0, 6.0, 2.0, 1.0]] * 3 + \
        [[1.0, 5.0, 0.0001, 0.0001], [2.0, 5.0, 0.0001, 0.0001]]
    assert lengths == [3, 2]


def test_extract_signals_chromosome_without_signal_is_ignored():
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 10, 30), (b"chrX", 0, 50)]),
        binsize=10)
    assert bins == [(b"chr1", 20), (b"chr1", 30)]
    assert lengths == [2]


def test_extract_signals_no_bins_raises():
    with pytest.raises(AssertionError):
        extract_signals_from_regions(
            make_signals(), regions_from([(b"chrX", 0, 50)]), binsize=10)


def test_extract_signals_tracks_of_unequal_length_raise():
    # the tri track ends at 100, so bins after 100 have no tri value
    spec = list(SIGNALS)
    spec[3] = {b"chr1": [(0, 100, 0.0)]}
    with pytest.raises(AssertionError):
        extract_signals_from_regions(
            make_signals(spec), regions_from([(b"chr1", 50, 150)]),
            binsize=10)


def test_extract_signals_requires_regions_object():
    with pytest.raises(AssertionError):
        extract_signals_from_regions(make_signals(), [(b"chr1", 0, 50)],
                                     binsize=10)


def test_extract_signals_single_bin():
    bins, data, lengths = extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 60, 70)]), binsize=10)
    assert bins == [(b"chr1", 70)]
    assert data == [[3.0, 5.0, 0.0001, 0.0001]]
    assert lengths == [1]


def test_extract_signals_float32_values():
    # the signals are float32 in bedGraphTrackI; 0.1 comes back as the
    # float32 nearest to 0.1
    spec = list(SIGNALS)
    spec[0] = {b"chr1": [(0, 200, 0.1)]}
    _, data, _ = extract_signals_from_regions(
        make_signals(spec), regions_from([(b"chr1", 0, 10)]), binsize=10)
    assert data[0][0] == float(np.float32(0.1))


def test_extract_signals_unknown_hmm_type_returns_no_bins():
    """Pins the current output.

    Only 'gaussian' and 'poisson' fill the result; for any other
    hmm_type no bin is added and the length list holds a single 0 for
    the closing region. The command line limits --hmm-type to those two
    choices, and nothing documents what an unknown type should do.
    """
    assert extract_signals_from_regions(
        make_signals(), regions_from([(b"chr1", 10, 50)]), binsize=10,
        hmm_type="gamma") == [[], [], [0]]


def test_extract_signals_weight_mapping_pipeline():
    # end-to-end: fragments -> weights -> digested signals -> bins; the
    # mono value of each bin equals the reference mono pileup at the
    # last base of the bin
    frags = [(b"chr1", 100 + 7 * i, 300 + 7 * i) for i in range(20)] + \
            [(b"chr1", 150 + 3 * i, 200 + 3 * i) for i in range(20)]
    pe = make_petrack(frags)
    mapping = generate_weight_mapping([50, 200], DEFAULT_MEANS, DEFAULT_SDS)
    signals = generate_digested_signals(pe, mapping)
    bins, data, lengths = extract_signals_from_regions(
        signals, regions_from([(b"chr1", 100, 400)]), binsize=10)
    end = max(r for _, l, r in frags)
    mono = np.zeros(end)
    for _, l, r in frags:
        mono[l:r] += mapping[1][r - l]
    assert lengths == [len(bins)]
    for (c, e), row in zip(bins, data):
        assert row[1] == pytest.approx(max(0.0001, mono[e - 1]), rel=1e-5)
    assert math.isclose(sum(lengths), len(data))
