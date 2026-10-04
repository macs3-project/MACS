#!/usr/bin/env python

"""Module Description: Test the HMMR_EM class, which fits the means and
standard deviations of the mono-, di- and tri-nucleosomal fragment
length distributions for HMMRATAC.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import logging
import re

import numpy as np
import pytest
from scipy.stats import norm

from MACS3.Signal.HMMR_EM import (HMMR_EM)
from MACS3.Signal.PairedEndTrack import (PETrackI,
                                         PETrackII)

# ------------------------------------
# Helpers
# ------------------------------------

EM_LOGGER = "MACS3.Signal.HMMR_EM"


def petrack_from_lengths(lengths, chrom=b"chr1", spacing=2000):
    """PETrackI with one fragment of each length, placed far apart."""
    pe = PETrackI()
    for i, ln in enumerate(lengths):
        pe.add_loc(chrom, i * spacing, i * spacing + int(ln))
    pe.finalize()
    return pe


def em_messages(caplog):
    """Messages logged by HMMR_EM, without the '[N MB] ' prefix."""
    return [re.sub(r"^\[\d+ MB\] ", "", r.getMessage())
            for r in caplog.records if r.name == EM_LOGGER]


def ref_hard_em(data, init_means, init_sds, min_fraglen=100,
                max_fraglen=1000, epsilon=0.05, max_iter=20, jump=1.5):
    """Reference implementation of the HMMR_EM algorithm in float64.

    Hard EM: every fragment length is assigned to the component with
    the largest weighted normal density, then each component's mean
    and variance move ``jump`` times the distance to the mean and
    (population) variance of its assigned lengths, and the weight
    becomes the fraction of assigned lengths. Initial weights split
    the lengths at the midpoints of the initial means. Stops when all
    means, stddevs and weights move less than ``epsilon``.
    Returns (means, stddevs, n_iterations, converged).
    """
    data = np.asarray(data, dtype=np.float64)
    data = data[(data >= min_fraglen) & (data <= max_fraglen)]
    means = np.array(init_means, dtype=np.float64)
    var = np.array(init_sds, dtype=np.float64) ** 2
    sds = np.sqrt(var)
    c1 = (init_means[1] - init_means[0]) / 2 + init_means[0]
    c2 = (init_means[2] - init_means[1]) / 2 + init_means[1]
    n = len(data)
    n1 = int((data < c1).sum())
    n2 = int((data < c2).sum())
    w = np.array([n1, n2 - n1, n - n2], dtype=np.float64) / n
    it = 0
    converged = False
    while True:
        old = (means.copy(), sds.copy(), w.copy())
        dens = w * norm.pdf(data[:, None], means, np.sqrt(var))
        best = dens.argmax(axis=1)
        ties = (dens == dens.max(axis=1, keepdims=True)).sum(axis=1) > 1
        best[ties] = -1
        total = int((best >= 0).sum())
        for j in range(3):
            sel = data[best == j]
            if len(sel) == 0:
                continue
            means[j] += jump * (sel.mean() - means[j])
            var[j] += jump * (sel.var() - var[j])
            sds[j] = np.sqrt(var[j])
            w[j] = len(sel) / total
        it += 1
        if (np.all(np.abs(old[0] - means) < epsilon) and
                np.all(np.abs(old[1] - sds) < epsilon) and
                np.all(np.abs(old[2] - w) < epsilon)):
            converged = True
            break
        if it >= max_iter:
            break
    return means, sds, it, converged


# lengths whose hard assignment reproduces the initial parameters
# exactly: each component gets two lengths at mean -/+ 10, so the
# assigned mean is the initial mean and the population variance is
# 100 (sd 10); the three components get equal weights both from the
# midpoint split (300, 500) and from the assignment.
FIXED_POINT = [190, 210, 390, 410, 590, 610]


# ------------------------------------
# HMMR_EM
# ------------------------------------

def test_hmmr_em_attributes_are_float32_arrays():
    em = HMMR_EM(petrack_from_lengths(FIXED_POINT), [200, 400, 600],
                 [10, 10, 10], sample_percentage=100)
    assert isinstance(em.fragMeans, np.ndarray)
    assert isinstance(em.fragStddevs, np.ndarray)
    assert em.fragMeans.dtype == np.float32
    assert em.fragStddevs.dtype == np.float32
    assert em.fragMeans.shape == (3,)
    assert em.fragStddevs.shape == (3,)


def test_hmmr_em_fixed_point(caplog):
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    em = HMMR_EM(petrack_from_lengths(FIXED_POINT), [200, 400, 600],
                 [10, 10, 10], sample_percentage=100)
    assert em.fragMeans.tolist() == [200.0, 400.0, 600.0]
    assert em.fragStddevs.tolist() == [10.0, 10.0, 10.0]
    msgs = em_messages(caplog)
    assert "# Downsampled 6 fragments will be used for EM training..." in msgs
    assert "# Reached convergence after 1 iterations" in msgs


def test_hmmr_em_one_iteration_by_hand(caplog):
    """Re-derived by hand for the default jump 0.5, which upstream 248fd6d
    ("add cli option for 'jump' for hmmr_em") changed from 1.5 (the update
    rule itself is unchanged). With jump 1.5 the mean was 185."""
    # lengths 180 and 200 go to the first component: their mean is 190
    # and variance 100. With the default jump 0.5 the mean moves from 200
    # to 200 + 0.5 * (190 - 200) = 195; the variance stays 100.
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    em = HMMR_EM(petrack_from_lengths([180, 200, 390, 410, 590, 610]),
                 [200, 400, 600], [10, 10, 10], sample_percentage=100,
                 maxIter=1)
    assert em.fragMeans.tolist() == [195.0, 400.0, 600.0]
    assert em.fragStddevs.tolist() == [10.0, 10.0, 10.0]
    assert "# Reached maximum number (1) of iterations" in em_messages(caplog)


def test_hmmr_em_overrelaxed_convergence_by_hand(caplog):
    """Upstream 248fd6d changed the default jump from 1.5 to 0.5 (and
    added --jump to hmmratac); this test of the over-relaxed update now
    passes jump=1.5 explicitly. The expectation is unchanged."""
    # The first mean starts at 200 and the target is 190; with jump 1.5
    # the error is multiplied by -0.5 each iteration:
    # 185, 192.5, 188.75, 190.625, 189.6875, 190.15625, 189.921875,
    # 190.0390625, 189.98046875, 190.009765625. The step at iteration
    # 10 (0.0293) is the first below epsilon 0.05. All values are exact
    # in float32.
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    em = HMMR_EM(petrack_from_lengths([180, 200, 390, 410, 590, 610]),
                 [200, 400, 600], [10, 10, 10], sample_percentage=100,
                 jump=1.5)
    assert em.fragMeans.tolist() == [190.009765625, 400.0, 600.0]
    assert em.fragStddevs.tolist() == [10.0, 10.0, 10.0]
    assert "# Reached convergence after 10 iterations" in em_messages(caplog)


def test_hmmr_em_jump_one(caplog):
    # jump 1: the mean goes straight to 190, and the second iteration
    # does not move anything
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    em = HMMR_EM(petrack_from_lengths([180, 200, 390, 410, 590, 610]),
                 [200, 400, 600], [10, 10, 10], sample_percentage=100,
                 jump=1.0)
    assert em.fragMeans.tolist() == [190.0, 400.0, 600.0]
    assert "# Reached convergence after 2 iterations" in em_messages(caplog)


def test_hmmr_em_epsilon(caplog):
    """Re-derived by hand for the default jump 0.5 (upstream 248fd6d
    changed it from 1.5; with 1.5 the mean was 189.6875 after 5
    iterations)."""
    # The first mean starts at 200 and the target is 190; with the
    # default jump 0.5 the error halves each iteration: 195, 192.5,
    # 191.25, 190.625 (steps 5, 2.5, 1.25, 0.625). With epsilon 1 the
    # step at iteration 4 (0.625) is the first below 1. All values are
    # exact in float32.
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    em = HMMR_EM(petrack_from_lengths([180, 200, 390, 410, 590, 610]),
                 [200, 400, 600], [10, 10, 10], sample_percentage=100,
                 epsilon=1.0)
    assert em.fragMeans.tolist() == [190.625, 400.0, 600.0]
    assert "# Reached convergence after 4 iterations" in em_messages(caplog)


def test_hmmr_em_variance_update_by_hand():
    # first component gets 170, 190, 210, 230: mean 200, population
    # variance (900+100+100+900)/4 = 500. With jump 1 the variance goes
    # from 100 to 500 in one iteration: sd = sqrt(500).
    em = HMMR_EM(petrack_from_lengths([170, 190, 210, 230, 390, 410,
                                       590, 610]),
                 [200, 400, 600], [10, 10, 10], sample_percentage=100,
                 jump=1.0, maxIter=1)
    assert em.fragMeans.tolist() == [200.0, 400.0, 600.0]
    assert em.fragStddevs[0] == pytest.approx(np.sqrt(500.0), rel=1e-6)
    assert em.fragStddevs[1:].tolist() == [10.0, 10.0]


def test_hmmr_em_component_without_fragments_keeps_initial_values(caplog):
    # max_fraglen 500 removes 590 and 610, so nothing is assigned to the
    # tri-nucleosome component and its parameters stay as given
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    em = HMMR_EM(petrack_from_lengths(FIXED_POINT), [200, 400, 650],
                 [10, 10, 33], sample_percentage=100, max_fraglen=500)
    assert em.fragMeans.tolist() == [200.0, 400.0, 650.0]
    assert em.fragStddevs.tolist() == [10.0, 10.0, 33.0]
    msgs = em_messages(caplog)
    assert "# Downsampled 4 fragments will be used for EM training..." in msgs


@pytest.mark.parametrize("extra, n_used", [
    ([], 6),
    ([50, 99], 6),          # below min_fraglen 100: dropped
    ([100], 7),             # min_fraglen is inclusive
    ([1000], 7),            # max_fraglen is inclusive
    ([1001, 1500], 6),      # above max_fraglen: dropped
])
def test_hmmr_em_fragment_length_range(caplog, extra, n_used):
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    HMMR_EM(petrack_from_lengths(FIXED_POINT + extra), [200, 400, 600],
            [10, 10, 10], sample_percentage=100)
    assert ("# Downsampled %d fragments will be used for EM training..."
            % n_used) in em_messages(caplog)


def test_hmmr_em_out_of_range_lengths_do_not_change_fit():
    em = HMMR_EM(petrack_from_lengths(FIXED_POINT + [20, 60, 1200, 3000]),
                 [200, 400, 600], [10, 10, 10], sample_percentage=100)
    assert em.fragMeans.tolist() == [200.0, 400.0, 600.0]
    assert em.fragStddevs.tolist() == [10.0, 10.0, 10.0]


def test_hmmr_em_custom_fraglen_range(caplog):
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    # jump 1 and one iteration: the single length left in the first
    # component has variance 0, which would break a second iteration
    HMMR_EM(petrack_from_lengths(FIXED_POINT), [200, 400, 600],
            [10, 10, 10], sample_percentage=100, min_fraglen=200,
            max_fraglen=600, jump=1.0, maxIter=1)
    # 210, 390, 410, 590 are inside [200, 600]
    assert "# Downsampled 4 fragments will be used for EM training..." in \
        em_messages(caplog)


@pytest.mark.parametrize("n, pct, n_used", [
    # n * pct / 100, rounded to 5 decimals and truncated to an integer
    (25, 10, 2),
    (30, 10, 3),
    (6, 50, 3),
    (7, 50, 3),
    (40, 25, 10),
])
def test_hmmr_em_sample_percentage(caplog, n, pct, n_used):
    caplog.set_level(logging.INFO, logger=EM_LOGGER)
    lengths = [200 + (i % 5) for i in range(n)]
    # jump 1 keeps the variance non-negative for any sample
    HMMR_EM(petrack_from_lengths(lengths), [200, 400, 600], [10, 10, 10],
            sample_percentage=pct, maxIter=1, jump=1.0)
    assert ("# Downsampled %d fragments will be used for EM training..."
            % n_used) in em_messages(caplog)


def test_hmmr_em_same_seed_same_result():
    rs = np.random.RandomState(3)
    lengths = np.concatenate([rs.normal(190, 20, 300), rs.normal(390, 25, 150),
                              rs.normal(590, 30, 80)]).astype(int)
    pe = petrack_from_lengths(lengths)
    a = HMMR_EM(pe, [200, 400, 600], [20, 20, 20], sample_percentage=50,
                seed=99)
    b = HMMR_EM(pe, [200, 400, 600], [20, 20, 20], sample_percentage=50,
                seed=99)
    assert a.fragMeans.tolist() == b.fragMeans.tolist()
    assert a.fragStddevs.tolist() == b.fragStddevs.tolist()


def test_hmmr_em_seed_irrelevant_without_downsampling():
    rs = np.random.RandomState(4)
    lengths = np.concatenate([rs.normal(190, 20, 200),
                              rs.normal(390, 25, 100),
                              rs.normal(590, 30, 50)]).astype(int)
    pe = petrack_from_lengths(lengths)
    a = HMMR_EM(pe, [200, 400, 600], [20, 20, 20], sample_percentage=100,
                seed=1)
    b = HMMR_EM(pe, [200, 400, 600], [20, 20, 20], sample_percentage=100,
                seed=2)
    assert a.fragMeans.tolist() == b.fragMeans.tolist()


def test_hmmr_em_does_not_modify_petrack():
    pe = petrack_from_lengths(FIXED_POINT * 5)
    before = pe.total
    HMMR_EM(pe, [200, 400, 600], [10, 10, 10], sample_percentage=50)
    assert pe.total == before
    assert sorted(pe.fraglengths().tolist()) == sorted(FIXED_POINT * 5)


# synthetic fragment lengths from three known normals
TRUE_MEANS = [180.0, 370.0, 560.0]
TRUE_SDS = [25.0, 30.0, 35.0]


def synthetic_lengths(seed=1234, n=(1200, 600, 300)):
    rs = np.random.RandomState(seed)
    return np.concatenate([np.rint(rs.normal(m, s, k))
                           for m, s, k in zip(TRUE_MEANS, TRUE_SDS, n)]
                          ).astype(int)


def test_hmmr_em_recovers_known_normals():
    lengths = synthetic_lengths()
    em = HMMR_EM(petrack_from_lengths(lengths), [200, 400, 600],
                 [20, 20, 20], sample_percentage=100)
    # sample means differ from the truth by about sd/sqrt(n) (< 2 bp);
    # the hard assignment trims the tails, which biases the sds down
    assert em.fragMeans.tolist() == pytest.approx(TRUE_MEANS, abs=4)
    assert em.fragStddevs.tolist() == pytest.approx(TRUE_SDS, rel=0.15)


def test_hmmr_em_matches_reference_hard_em():
    lengths = synthetic_lengths()
    em = HMMR_EM(petrack_from_lengths(lengths), [200, 400, 600],
                 [20, 20, 20], sample_percentage=100)
    means, sds, it, converged = ref_hard_em(lengths, [200, 400, 600],
                                            [20, 20, 20])
    assert converged
    # float32 accumulation in HMMR_EM versus float64 here
    assert em.fragMeans.tolist() == pytest.approx(means.tolist(), abs=0.1)
    assert em.fragStddevs.tolist() == pytest.approx(sds.tolist(), abs=0.1)


@pytest.mark.parametrize("init_means, init_sds", [
    ([150, 350, 550], [30, 30, 30]),
    ([200, 400, 600], [20, 20, 20]),
    ([220, 420, 620], [15, 15, 15]),
])
def test_hmmr_em_initial_values(init_means, init_sds):
    # different starting points converge to the same neighbourhood
    lengths = synthetic_lengths(seed=99)
    em = HMMR_EM(petrack_from_lengths(lengths), init_means, init_sds,
                 sample_percentage=100)
    assert em.fragMeans.tolist() == pytest.approx(TRUE_MEANS, abs=5)


def test_hmmr_em_two_chromosomes():
    lengths = FIXED_POINT
    pe = PETrackI()
    for i, ln in enumerate(lengths):
        chrom = b"chr1" if i % 2 else b"chr2"
        pe.add_loc(chrom, i * 1000, i * 1000 + ln)
    pe.finalize()
    em = HMMR_EM(pe, [200, 400, 600], [10, 10, 10], sample_percentage=100)
    assert em.fragMeans.tolist() == [200.0, 400.0, 600.0]
    assert em.fragStddevs.tolist() == [10.0, 10.0, 10.0]


def test_hmmr_em_petrackII_counts():
    # FRAG data: a fragment with count 2 counts as two lengths
    pe = PETrackII()
    pe.add_loc(b"chr1", 0, 190, b"A", 2)
    pe.add_loc(b"chr1", 1000, 1210, b"A", 2)
    pe.add_loc(b"chr1", 2000, 2390, b"B", 1)
    pe.add_loc(b"chr1", 3000, 3410, b"B", 1)
    pe.add_loc(b"chr1", 4000, 4590, b"C", 1)
    pe.add_loc(b"chr1", 5000, 5610, b"C", 1)
    pe.finalize()
    em = HMMR_EM(pe, [200, 400, 600], [10, 10, 10], sample_percentage=100)
    # the online variance of 190, 190, 210, 210 has float32 rounding
    assert em.fragMeans.tolist() == pytest.approx([200.0, 400.0, 600.0],
                                                  rel=1e-6)
    assert em.fragStddevs.tolist() == pytest.approx([10.0, 10.0, 10.0],
                                                    rel=1e-5)


def test_hmmr_em_no_fragment_in_range_raises():
    with pytest.raises(ZeroDivisionError):
        HMMR_EM(petrack_from_lengths([20, 30, 40, 50]), [200, 400, 600],
                [10, 10, 10], sample_percentage=100)


def test_hmmr_em_too_few_initial_means_raises():
    with pytest.raises(IndexError):
        HMMR_EM(petrack_from_lengths(FIXED_POINT), [200, 400], [10, 10],
                sample_percentage=100)


def test_hmmr_em_requires_list_parameters():
    with pytest.raises(TypeError):
        HMMR_EM(petrack_from_lengths(FIXED_POINT), (200, 400, 600),
                [10, 10, 10], sample_percentage=100)


TIGHT = [198, 199, 200, 201, 202] * 20 + [390, 410] * 10 + [590, 610] * 10


def test_hmmr_em_tight_cluster_current_handling(capsys):
    """Upstream 248fd6d changed the default jump from 1.5 to 0.5. With
    jump <= 1 the variance update is a convex combination of two
    variances and cannot go negative, so this test of the negative-variance
    path now passes jump=1.5 explicitly (reachable with hmmratac --jump
    1.5). The expectation is unchanged."""
    # the first component's lengths have variance 2 while the initial
    # variance is 400, so the update 400 + 1.5 * (2 - 400) = -197 is
    # negative; sqrt fails and this message is printed
    HMMR_EM(petrack_from_lengths(TIGHT), [200, 400, 600], [20, 20, 20],
            sample_percentage=100, jump=1.5)
    assert capsys.readouterr().out == \
        " ValueError:  Adjust --means and --stddevs options and re-run " \
        "command\n"


def test_hmmr_em_sparse_tight_component_exits(capsys):
    """Upstream 248fd6d changed the default jump from 1.5 to 0.5, with
    which the variance cannot go negative; this test now passes jump=1.5
    explicitly. The expectation is unchanged."""
    # The over-relaxed update (jump 1.5) is HMMRATAC's published EM step,
    # which HMMRATAC applies to the stddev and MACS3 to the variance.
    # Neither can fit a component whose spread is far below its initial
    # value. Here the tri-nucleosome component gets two lengths 2 bp
    # apart: the variance update 400 + 1.5 * (1 - 400) is negative, and
    # so is HMMRATAC's stddev update for these two lengths,
    # 20 + 1.5 * (sqrt(2) - 20). HMMR_EM prints advice to change the
    # initial values; the other components moved, so EM starts a second
    # iteration, where pnorm2 meets the negative variance and calls
    # sys.exit(1).
    lengths = [170, 185, 200, 215, 230] * 10 + [370, 400, 430] * 5 + \
        [599, 601]
    with pytest.raises(SystemExit) as exc:
        HMMR_EM(petrack_from_lengths(lengths), [200, 400, 600],
                [20, 20, 20], sample_percentage=100, jump=1.5)
    assert exc.value.code == 1
    assert capsys.readouterr().out == \
        " ValueError:  Adjust --means and --stddevs options and re-run " \
        "command\n"


