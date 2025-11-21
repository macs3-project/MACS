# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-11-20 16:57:24 Tao Liu>

"""Scoring utilities for MACS3 signal tracks and peak callers.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from functools import reduce

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.SignalProcessing import maxima, enforce_peakyness
from MACS3.Signal.Prob import poisson_cdf
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO

# ------------------------------------
# Other modules
# ------------------------------------
import cython
import numpy as np
import cython.cimports.numpy as cnp
from cython.cimports.cpython import bool
from cykhash import PyObjectMap, Float32to32Map

# ------------------------------------
# C lib
# ------------------------------------
from cython.cimports.libc.math import (log10,
                                       log)

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.inline
@cython.cfunc
def int_max(a: cython.int, b: cython.int) -> cython.int:
    """Return the larger of ``a`` and ``b``."""
    return a if a >= b else b


@cython.inline
@cython.cfunc
def int_min(a: cython.int, b: cython.int) -> cython.int:
    """Return the smaller of ``a`` and ``b``."""
    return a if a <= b else b


LOG10_E: cython.float = 0.43429448190325176

pscore_dict = PyObjectMap()


@cython.cfunc
def get_pscore(observed: cython.int,
               expectation: cython.float) -> cython.float:
    """Return cached ``-log10`` Poisson tail probability for ``observed``."""
    score: cython.double

    try:
        return pscore_dict[(observed, expectation)]
    except KeyError:
        score = -1 * poisson_cdf(observed,
                                 expectation,
                                 False,
                                 True)
        pscore_dict[(observed, expectation)] = score
        return score


asym_logLR_dict = PyObjectMap()


@cython.cfunc
def logLR_asym(x: cython.float,
               y: cython.float) -> cython.float:
    """Return asymmetric ``log10`` likelihood ratio between ``x`` and ``y``."""
    s: cython.float

    if (x, y) in asym_logLR_dict:
        return asym_logLR_dict[(x, y)]
    else:
        if x > y:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
        elif x < y:
            s = (x*(-log(x)+log(y))-y+x)*LOG10_E
        else:
            s = 0
        asym_logLR_dict[(x, y)] = s
        return s


sym_logLR_dict = PyObjectMap()


@cython.cfunc
def logLR_sym(x: cython.float, y: cython.float) -> cython.float:
    """Return symmetric ``log10`` likelihood ratio between ``x`` and ``y``."""
    s: cython.float

    if (x, y) in sym_logLR_dict:
        return sym_logLR_dict[(x, y)]
    else:
        if x > y:
            s = (x*(log(x)-log(y))+y-x)*LOG10_E
        elif y > x:
            s = (y*(log(x)-log(y))+y-x)*LOG10_E
        else:
            s = 0
        sym_logLR_dict[(x, y)] = s
        return s


@cython.inline
@cython.cfunc
def get_logFE(x: cython.float, y: cython.float) -> cython.float:
    """Return ``log10`` fold enrichment (base-10) for ``x`` over ``y``."""
    return log10(x/y)


@cython.cfunc
def get_subtraction(x: cython.float, y: cython.float) -> cython.float:
    """Return the difference ``x - y``."""
    return x - y

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class ScoreTrackII:
    """Container for treatment/control pileups and derived score tracks."""
    # dictionary for data of each chromosome
    data: dict
    # length of data array of each chromosome
    datalength: dict
    # whether trackline should be saved in bedGraph
    trackline: bool
    # seq depth in million of treatment
    treat_edm: cython.float
    # seq depth in million of control
    ctrl_edm: cython.float
    # method for calculating scores.
    scoring_method: cython.char
    # scale to control? scale to treatment? both scale to 1million reads?
    normalization_method: cython.char
    # the pseudocount used to calcuate logLR, FE or logFE
    pseudocount: cython.float
    # cutoff
    cutoff: cython.float
    # save pvalue<->length dictionary
    pvalue_stat: dict

    def __init__(self,
                 treat_depth: cython.float,
                 ctrl_depth: cython.float,
                 pseudocount: cython.float = 1.0):
        """Initialise score containers with effective library depths.

        Args:
            treat_depth: Effective treatment depth (millions of filtered reads).
            ctrl_depth: Effective control depth (millions of filtered reads).
            pseudocount: Pseudocount added when computing score metrics.
        """
        # for each chromosome, there is a l*4 matrix. First column:
        # end position of a region; Second: treatment pileup; third:
        # control pileup ; forth: score (can be p/q-value/likelihood
        # ratio/fold-enrichment/subtraction depending on -c setting)
        self.data = {}

        self.datalength = {}
        self.trackline = False
        self.treat_edm = treat_depth
        self.ctrl_edm = ctrl_depth

        # scoring_method:  p: -log10 pvalue;
        #                  q: -log10 qvalue;
        #                  l: log10 likelihood ratio (minus for depletion)
        #                  f: log10 fold enrichment
        #                  F: linear fold enrichment
        #                  d: subtraction
        #                  m: fragment pileup per million reads
        #                  N: not set
        self.scoring_method = ord("N")

        # normalization_method: T: scale to depth of treatment;
        #                       C: scale to depth of control;
        #                       M: scale to depth of 1 million;
        #                       N: not set/ raw pileup
        self.normalization_method = ord("N")

        self.pseudocount = pseudocount
        self.pvalue_stat = {}

    @cython.ccall
    def set_pseudocount(self, pseudocount: cython.float):
        """Update the pseudocount used when computing score metrics."""
        self.pseudocount = pseudocount

    @cython.ccall
    def enable_trackline(self):
        """Enable UCSC track line output when exporting bedGraphs."""
        self.trackline = True

    @cython.ccall
    def add_chromosome(self,
                       chrom: bytes,
                       chrom_max_len: cython.int = 200000000):
        """Allocate arrays for ``chrom`` with capacity ``chrom_max_len``."""
        if chrom not in self.data:
            self.data[chrom] = [np.zeros(chrom_max_len, dtype="int32"),  # pos
                                # pileup at each interval, in float32 format
                                np.zeros(chrom_max_len, dtype="float32"),
                                # control at each interval, in float32 format
                                np.zeros(chrom_max_len, dtype="float32"),
                                # score at each interval, in float32 format
                                np.zeros(chrom_max_len, dtype="float32")]
            self.datalength[chrom] = 0

    @cython.ccall
    def add(self,
            chromosome: bytes,
            endpos: cython.int,
            chip: cython.float,
            control: cython.float):
        """Append treatment/control pileup ending at ``endpos`` for ``chromosome``."""
        i: cython.int

        i = self.datalength[chromosome]
        c = self.data[chromosome]
        c[0][i] = endpos
        c[1][i] = chip
        c[2][i] = control
        self.datalength[chromosome] += 1

    @cython.ccall
    def finalize(self):
        """Trim per-chromosome arrays to their populated length."""
        chrom: bytes
        ln: cython.int

        for chrom in sorted(self.data.keys()):
            d = self.data[chrom]
            ln = self.datalength[chrom]
            d[0].resize(ln, refcheck=False)
            d[1].resize(ln, refcheck=False)
            d[2].resize(ln, refcheck=False)
            d[3].resize(ln, refcheck=False)
        return

    @cython.ccall
    def get_data_by_chr(self,
                        chromosome: bytes):
        """Return ``(positions, treatment, control, score)`` arrays for ``chromosome``."""
        if chromosome in self.data:
            return self.data[chromosome]
        else:
            return None

    @cython.ccall
    def get_chr_names(self):
        """Return all the chromosome names stored.

        """
        return set(self.data.keys())

    @cython.ccall
    def change_normalization_method(self,
                                    normalization_method: cython.char):
        """Change/set normalization method. However, I do not
        recommend change this back and forward, since some precision
        issue will happen -- I only keep two digits.

        normalization_method: T: scale to depth of treatment;
                              C: scale to depth of control;
                              M: scale to depth of 1 million;
                              N: not set/ raw pileup
        """
        if normalization_method == ord('T'):
            if self.normalization_method == ord('T'):  # do nothing
                pass
            elif self.normalization_method == ord('C'):
                self.normalize(self.treat_edm/self.ctrl_edm,
                               self.treat_edm/self.ctrl_edm)
            elif self.normalization_method == ord('M'):
                self.normalize(self.treat_edm, self.treat_edm)
            elif self.normalization_method == ord('N'):
                self.normalize(1, self.treat_edm/self.ctrl_edm)
            else:
                raise NotImplementedError
            self.normalization_method = ord('T')
        elif normalization_method == ord('C'):
            if self.normalization_method == ord('T'):
                self.normalize(self.ctrl_edm/self.treat_edm,
                               self.ctrl_edm/self.treat_edm)
            elif self.normalization_method == ord('C'):  # do nothing
                pass
            elif self.normalization_method == ord('M'):
                self.normalize(self.ctrl_edm, self.ctrl_edm)
            elif self.normalization_method == ord('N'):
                self.normalize(self.ctrl_edm/self.treat_edm, 1)
            else:
                raise NotImplementedError
            self.normalization_method = ord('C')
        elif normalization_method == ord('M'):
            if self.normalization_method == ord('T'):
                self.normalize(1/self.treat_edm,
                               1/self.treat_edm)
            elif self.normalization_method == ord('C'):
                self.normalize(1/self.ctrl_edm,
                               1/self.ctrl_edm)
            elif self.normalization_method == ord('M'):  # do nothing
                pass
            elif self.normalization_method == ord('N'):
                self.normalize(1/self.treat_edm,
                               1/self.ctrl_edm)
            else:
                raise NotImplementedError
            self.normalization_method = ord('M')
        elif normalization_method == ord('N'):
            if self.normalization_method == ord('T'):
                self.normalize(self.treat_edm,
                               self.treat_edm)
            elif self.normalization_method == ord('C'):
                self.normalize(self.ctrl_edm,
                               self.ctrl_edm)
            elif self.normalization_method == ord('M'):
                self.normalize(self.treat_edm,
                               self.ctrl_edm)
            elif self.normalization_method == ord('N'):  # do nothing
                pass
            else:
                raise NotImplementedError
            self.normalization_method = ord('N')

    @cython.cfunc
    def normalize(self,
                  treat_scale: cython.float,
                  control_scale: cython.float):
        """Scale treatment and control pileups in-place by the given factors."""
        p: cnp.ndarray
        c: cnp.ndarray
        ln: cython.long
        i: cython.long

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            ln = self.datalength[chrom]
            for i in range(ln):
                p[i] *= treat_scale
                c[i] *= control_scale
        return

    @cython.ccall
    def change_score_method(self,
                            scoring_method: cython.char):
        """
        scoring_method:  p: -log10 pvalue;
                         q: -log10 qvalue;
                         l: log10 likelihood ratio (minus for depletion)
                         s: symmetric log10 likelihood ratio (for comparing two
                            ChIPs)
                         f: log10 fold enrichment
                         F: linear fold enrichment
                         d: subtraction
                         M: maximum
                         m: fragment pileup per million reads
        """
        if scoring_method == ord('p'):
            self.compute_pvalue()
        elif scoring_method == ord('q'):
            # if not already calculated p, compute pvalue first
            if self.scoring_method != ord('p'):
                self.compute_pvalue()
            self.compute_qvalue()
        elif scoring_method == ord('l'):
            self.compute_likelihood()
        elif scoring_method == ord('s'):
            self.compute_sym_likelihood()
        elif scoring_method == ord('f'):
            self.compute_logFE()
        elif scoring_method == ord('F'):
            self.compute_foldenrichment()
        elif scoring_method == ord('d'):
            self.compute_subtraction()
        elif scoring_method == ord('m'):
            self.compute_SPMR()
        elif scoring_method == ord('M'):
            self.compute_max()
        else:
            raise NotImplementedError

    @cython.cfunc
    def compute_pvalue(self):
        """Compute -log_{10}(pvalue)
        """
        p: cnp.ndarray
        c: cnp.ndarray
        v: cnp.ndarray
        pos: cnp.ndarray
        ln: cython.long
        tmp_l: cython.long
        i: cython.long
        prev_pos: cython.long
        chrom: bytes

        for chrom in sorted(self.data.keys()):
            prev_pos = 0
            pos = self.data[chrom][0]
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = get_pscore(cython.cast(cython.int,
                                              (p[i] + self.pseudocount)),
                                  c[i] + self.pseudocount)
                tmp_l = pos[i] - prev_pos
                try:
                    self.pvalue_stat[v[i]] += tmp_l
                except Exception:
                    self.pvalue_stat[v[i]] = tmp_l
                prev_pos = pos[i]

        self.scoring_method = ord('p')
        return

    @cython.cfunc
    def compute_qvalue(self):
        """Compute -log_{10}(qvalue)
        """
        pqtable: object
        i: cython.long
        ln: cython.long
        chrom: bytes
        v: cnp.ndarray

        # pvalue should be computed first!
        assert self.scoring_method == ord('p')
        # make pqtable
        pqtable = self.make_pq_table()

        # convert p to q
        for chrom in sorted(self.data.keys()):
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = pqtable[v[i]]

        self.scoring_method = ord('q')
        return

    @cython.ccall
    def make_pq_table(self):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for
               each pvalue

        Return a dictionary of
        {-log10pvalue:(-log10qvalue,rank,basepairs)} relationships.

        """
        ln: cython.long
        i: cython.long
        j: cython.long
        v: cython.float
        q: cython.float
        pre_q: cython.float     # store the p and q scores
        N: cython.long
        k: cython.float
        f: cython.float
        pvalue2qvalue: object
        pvalue_stat: dict
        unique_values: list

        assert self.scoring_method == ord('p')

        pvalue_stat = self.pvalue_stat

        N = sum(pvalue_stat.values())
        k = 1                           # rank
        f = -log10(N)
        pre_q = 2147483647              # save the previous q-value

        pvalue2qvalue = Float32to32Map(for_int=False)
        unique_values = sorted(list(pvalue_stat.keys()), reverse=True)
        for i in range(len(unique_values)):
            v = unique_values[i]
            ln = pvalue_stat[v]
            q = v + (log10(k) + f)
            if q > pre_q:
                q = pre_q
            if q <= 0:
                q = 0
                break
            pvalue2qvalue[v] = q
            pre_q = q
            k += ln
        # bottom rank pscores all have qscores 0
        for j in range(i, len(unique_values)):
            v = unique_values[j]
            pvalue2qvalue[v] = 0
        return pvalue2qvalue

    @cython.cfunc
    def compute_likelihood(self):
        """Calculate log10 likelihood.

        """
        ln: cython.long
        i: cython.long
        chrom: bytes
        v1: cython.float
        v2: cython.float
        pseudocount: cython.float

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1].flat.__next__  # pileup in treatment
            c = self.data[chrom][2].flat.__next__  # pileup in control
            v = self.data[chrom][3]                # score
            ln = self.datalength[chrom]
            v1 = 2
            v2 = 1
            for i in range(ln):
                v1 = p()
                v2 = c()
                v[i] = logLR_asym(v1 + pseudocount, v2 + pseudocount)
        self.scoring_method = ord('l')
        return

    @cython.cfunc
    def compute_sym_likelihood(self):
        """Calculate symmetric log10 likelihood.

        """
        ln: cython.long
        i: cython.long
        chrom: bytes
        v1: cython.float
        v2: cython.float
        pseudocount: cython.float

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1].flat.__next__
            c = self.data[chrom][2].flat.__next__
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            v1 = 2
            v2 = 1
            for i in range(ln):
                v1 = p()
                v2 = c()
                v[i] = logLR_sym(v1 + pseudocount, v2 + pseudocount)
        self.scoring_method = ord('s')
        return

    @cython.cfunc
    def compute_logFE(self):
        """Calculate log10 fold enrichment (with 1 pseudocount).

        """
        p: cnp.ndarray
        c: cnp.ndarray
        v: cnp.ndarray
        ln: cython.long
        i: cython.long
        pseudocount: cython.float

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = get_logFE(p[i] + pseudocount, c[i] + pseudocount)
        self.scoring_method = ord('f')
        return

    @cython.cfunc
    def compute_foldenrichment(self):
        """Calculate linear scale fold enrichment (with 1 pseudocount).

        """
        p: cnp.ndarray
        c: cnp.ndarray
        v: cnp.ndarray
        ln: cython.long
        i: cython.long
        pseudocount: cython.float

        pseudocount = self.pseudocount

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = (p[i] + pseudocount)/(c[i] + pseudocount)
        self.scoring_method = ord('F')
        return

    @cython.cfunc
    def compute_subtraction(self):
        """Populate scores with treatment minus control pileup."""
        p: cnp.ndarray
        c: cnp.ndarray
        v: cnp.ndarray
        ln: cython.long
        i: cython.long

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = p[i] - c[i]
        self.scoring_method = ord('d')
        return

    @cython.cfunc
    def compute_SPMR(self):
        """Populate scores with treatment pileup per million reads."""
        p: cnp.ndarray
        v: cnp.ndarray
        ln: cython.long
        i: cython.long
        scale: cython.float

        if self.normalization_method == ord('T') or self.normalization_method == ord('N'):
            scale = self.treat_edm
        elif self.normalization_method == ord('C'):
            scale = self.ctrl_edm
        elif self.normalization_method == ord('M'):
            scale = 1

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = p[i] / scale  # two digit precision may not be enough...
        self.scoring_method = ord('m')
        return

    @cython.cfunc
    def compute_max(self):
        """Populate scores with the element-wise maximum of treatment and control."""
        p: cnp.ndarray
        c: cnp.ndarray
        v: cnp.ndarray
        ln: cython.long
        i: cython.long

        for chrom in sorted(self.data.keys()):
            p = self.data[chrom][1]
            c = self.data[chrom][2]
            v = self.data[chrom][3]
            ln = self.datalength[chrom]
            for i in range(ln):
                v[i] = max(p[i], c[i])
        self.scoring_method = ord('M')
        return

    @cython.ccall
    def write_bedGraph(self,
                       fhd,
                       name: str,
                       description: str,
                       column: cython.short = 3):
        """Write all data to fhd in bedGraph Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 1: chip, 2: control, 3: score

        """
        chrom: bytes
        ln: cython.int
        pre: cython.int
        i: cython.int
        p: cython.int
        pre_v: cython.float
        v: cython.float
        chrs: set
        pos: cnp.ndarray
        value: cnp.ndarray

        assert column in range(1, 4), "column should be between 1, 2 or 3."

        write = fhd.write

        if self.trackline:
            # this line is REQUIRED by the wiggle format for UCSC browser
            write("track type=bedGraph name=\"%s\" description=\"%s\"\n" %
                  (name.decode(), description))

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            pos = self.data[chrom][0]
            value = self.data[chrom][column]
            ln = self.datalength[chrom]
            pre = 0
            if pos.shape[0] == 0:
                continue  # skip if there's no data
            pre_v = value[0]
            for i in range(1, ln):
                v = value[i]
                p = pos[i-1]
                if abs(pre_v - v) > 1e-5:  # precision is 5 digits
                    write("%s\t%d\t%d\t%.5f\n" %
                          (chrom.decode(), pre, p, pre_v))
                    pre_v = v
                    pre = p
            p = pos[-1]
            # last one
            write("%s\t%d\t%d\t%.5f\n" %
                  (chrom.decode(), pre, p, pre_v))

        return True

    @cython.ccall
    def cutoff_analysis(self,
                        max_gap: cython.int = 50,
                        min_length: cython.int = 200,
                        steps: cython.int = 100,
                        min_score: cython.float = 0,
                        max_score: cython.float = 1000) -> str:
        """Summarise peak metrics across a range of score thresholds.

        Args:
            max_gap: Maximum distance between merged regions.
            min_length: Minimum peak length to keep.
            steps: Number of cutoff increments between the observed
                minimum and maximum scores.
            min_score: Lower bound for the cutoff sweep.
            max_score: Upper bound for the cutoff sweep.

        Returns:
            str: Tab-delimited report of peak counts and lengths per cutoff.
        """
        chrs: set
        peak_content: list
        ret_list: list
        cutoff_list: list
        cutoff_npeaks: list
        cutoff_lpeaks: list
        chrom: bytes
        ret: str
        cutoff: cython.float
        total_l: cython.long
        total_p: cython.long
        i: cython.long
        n: cython.long
        ts: cython.long
        te: cython.long
        lastp: cython.long
        tl: cython.long
        peak_length: cython.long
        s: cython.float
        ln: cython.long
        chrom_min: cython.float
        chrom_max: cython.float
        obs_min: cython.float
        obs_max: cython.float
        minv: cython.float
        maxv: cython.float

        chrs = self.get_chr_names()

        obs_min = 10000000.0
        obs_max = -10000000.0
        for chrom in chrs:
            ln = self.datalength[chrom]
            if ln == 0:
                continue
            chrom_scores = self.data[chrom][3][:ln]
            if chrom_scores.size == 0:
                continue
            chrom_min = float(np.min(chrom_scores))
            chrom_max = float(np.max(chrom_scores))
            if chrom_min < obs_min:
                obs_min = chrom_min
            if chrom_max > obs_max:
                obs_max = chrom_max

        minv = max(min_score, obs_min)
        maxv = min(obs_max, max_score)
        if steps <= 0 or maxv <= minv:
            return "score\tnpeaks\tlpeaks\tavelpeak\n"

        s = float(maxv - minv)/steps
        if s <= 0:
            return "score\tnpeaks\tlpeaks\tavelpeak\n"

        cutoff_list = [round(value, 3) for value in np.arange(minv, maxv, s)]
        if not cutoff_list:
            return "score\tnpeaks\tlpeaks\tavelpeak\n"

        cutoff_npeaks = [0] * len(cutoff_list)
        cutoff_lpeaks = [0] * len(cutoff_list)

        for chrom in sorted(chrs):
            ln = self.datalength[chrom]
            if ln == 0:
                continue
            pos_array = self.data[chrom][0][:ln]
            score_array = self.data[chrom][3][:ln]

            for n in range(len(cutoff_list)):
                cutoff = cutoff_list[n]
                total_l = 0           # total length of peaks
                total_p = 0           # total number of peaks

                # get the regions with scores above cutoffs. This is
                # not an optimized method. It would be better to store
                # score array in a 2-D ndarray?
                above_cutoff = np.nonzero(score_array > cutoff)[0]
                # end positions of regions where score is above cutoff
                above_cutoff_endpos = pos_array[above_cutoff]
                # start positions of regions where score is above cutoff
                above_cutoff_startpos = pos_array[above_cutoff-1]

                if above_cutoff_endpos.size == 0:
                    continue

                # first bit of region above cutoff
                acs_next = iter(above_cutoff_startpos).__next__
                ace_next = iter(above_cutoff_endpos).__next__

                ts = acs_next()
                te = ace_next()
                peak_content = [(ts, te),]
                lastp = te

                for i in range(1, above_cutoff_startpos.size):
                    ts = acs_next()
                    te = ace_next()
                    tl = ts - lastp
                    if tl <= max_gap:
                        peak_content.append((ts, te))
                    else:
                        peak_length = peak_content[-1][1] - peak_content[0][0]
                        # if the peak is too small, reject it
                        if peak_length >= min_length:
                            total_l += peak_length
                            total_p += 1
                        peak_content = [(ts, te),]
                    lastp = te

                if peak_content:
                    peak_length = peak_content[-1][1] - peak_content[0][0]
                    # if the peak is too small, reject it
                    if peak_length >= min_length:
                        total_l += peak_length
                        total_p += 1
                cutoff_lpeaks[n] += total_l
                cutoff_npeaks[n] += total_p

        # prepare the returnning text
        ret_list = ["score\tnpeaks\tlpeaks\tavelpeak\n"]
        for n in range(len(cutoff_list)-1, -1, -1):
            cutoff = cutoff_list[n]
            if cutoff_npeaks[n] > 0:
                ret_list.append("%.2f\t%d\t%d\t%.2f\n" % (cutoff,
                                                          cutoff_npeaks[n],
                                                          cutoff_lpeaks[n],
                                                          cutoff_lpeaks[n]/cutoff_npeaks[n]))
        ret = ''.join(ret_list)
        return ret

    @cython.ccall
    def call_peaks(self,
                   cutoff: cython.float = 5.0,
                   min_length: cython.int = 200,
                   max_gap: cython.int = 50,
                   call_summits: bool = False):
        """Return peaks where scores remain above ``cutoff``.

        Args:
            cutoff: Minimum score threshold (e.g., ``-log10 p``).
            min_length: Minimum peak length in bases.
            max_gap: Maximum distance between merged segments.
            call_summits: Whether to report all local maxima within peaks.
        """
        i: cython.int
        chrom: bytes
        pos: cnp.ndarray
        sample: cnp.ndarray
        control: cnp.ndarray
        value: cnp.ndarray
        above_cutoff: cnp.ndarray
        above_cutoff_v: cnp.ndarray
        above_cutoff_endpos: cnp.ndarray
        above_cutoff_startpos: cnp.ndarray
        above_cutoff_sv: cnp.ndarray
        peak_content: list

        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks

        self.cutoff = cutoff
        for chrom in sorted(chrs):
            peak_content = []           # to store points above cutoff

            pos = self.data[chrom][0]
            sample = self.data[chrom][1]
            # control = self.data[chrom][2]
            value = self.data[chrom][3]

            # indices where score is above cutoff
            above_cutoff = np.nonzero(value >= cutoff)[0]
            # scores where score is above cutoff
            above_cutoff_v = value[above_cutoff]
            # end positions of regions where score is above cutoff
            above_cutoff_endpos = pos[above_cutoff]
            # start positions of regions where score is above cutoff
            above_cutoff_startpos = pos[above_cutoff-1]
            # sample pileup height where score is above cutoff
            above_cutoff_sv = sample[above_cutoff]
            if above_cutoff_v.size == 0:
                # nothing above cutoff
                continue

            if above_cutoff[0] == 0:
                # first element > cutoff, fix the first point as
                # 0. otherwise it would be the last item in
                # data[chrom]['pos']
                above_cutoff_startpos[0] = 0

            # first bit of region above cutoff
            peak_content.append((above_cutoff_startpos[0],
                                 above_cutoff_endpos[0],
                                 above_cutoff_v[0],
                                 above_cutoff_sv[0],
                                 above_cutoff[0]))
            for i in range(1, above_cutoff_startpos.size):
                if above_cutoff_startpos[i] - peak_content[-1][1] <= max_gap:
                    # append
                    peak_content.append((above_cutoff_startpos[i],
                                         above_cutoff_endpos[i],
                                         above_cutoff_v[i],
                                         above_cutoff_sv[i],
                                         above_cutoff[i]))
                else:
                    # close
                    if call_summits:
                        self.__close_peak2(peak_content,
                                           peaks,
                                           min_length,
                                           chrom,
                                           max_gap//2)
                    else:
                        self.__close_peak(peak_content,
                                          peaks,
                                          min_length,
                                          chrom)
                    peak_content = [(above_cutoff_startpos[i],
                                     above_cutoff_endpos[i],
                                     above_cutoff_v[i],
                                     above_cutoff_sv[i],
                                     above_cutoff[i]),]

            # save the last peak
            if not peak_content:
                continue
            else:
                if call_summits:
                    self.__close_peak2(peak_content,
                                       peaks,
                                       min_length,
                                       chrom,
                                       max_gap//2)
                else:
                    self.__close_peak(peak_content,
                                      peaks,
                                      min_length,
                                      chrom)

        return peaks

    @cython.cfunc
    def __close_peak(self,
                     peak_content: list,
                     peaks: object,
                     min_length: cython.int,
                     chrom: bytes) -> bool:
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        In this function, we define the peak summit as the middle
        point of the region with the highest score, in this peak. For
        example, if the region of the highest score is from 100 to
        200, the summit is 150. If there are several regions of the
        same 'highest score', we will first calculate the possible
        summit for each such region, then pick a position close to the
        middle index (= (len(highest_regions) + 1) / 2) of these
        summits. For example, if there are three regions with the same
        highest scores, [100,200], [300,400], [600,700], we will first
        find the possible summits as 150, 350, and 650, and then pick
        the middle index, the 2nd, of the three positions -- 350 as
        the final summit. If there are four regions, we pick the 2nd
        as well.

        peaks: a PeakIO object

        """
        summit_pos: cython.int
        tstart: cython.int
        tend: cython.int
        summit_index: cython.int
        i: cython.int
        midindex: cython.int
        summit_value: cython.float
        tvalue: cython.float
        tsummitvalue: cython.float

        peak_length = peak_content[-1][1] - peak_content[0][0]
        if peak_length >= min_length:  # if the peak is too small, reject it
            tsummit = []
            summit_pos = 0
            summit_value = 0
            for i in range(len(peak_content)):
                (tstart, tend, tvalue, tsummitvalue, tindex) = peak_content[i]
                #for (tstart,tend,tvalue,tsummitvalue, tindex) in peak_content:
                if not summit_value or summit_value < tsummitvalue:
                    tsummit = [(tend + tstart) / 2,]
                    tsummit_index = [tindex,]
                    summit_value = tsummitvalue
                elif summit_value == tsummitvalue:
                    # remember continuous summit values
                    tsummit.append(int((tend + tstart) / 2))
                    tsummit_index.append(tindex)
            # the middle of all highest points in peak region is defined as summit
            midindex = int((len(tsummit) + 1) / 2) - 1
            summit_pos = tsummit[midindex]
            summit_index = tsummit_index[midindex]
            if self.scoring_method == ord('q'):
                qscore = self.data[chrom][3][summit_index]
            else:
                # if q value is not computed, use -1
                qscore = -1

            peaks.add(chrom,
                      peak_content[0][0],
                      peak_content[-1][1],
                      summit=summit_pos,
                      peak_score=self.data[chrom][3][summit_index],
                      # should be the same as summit_value
                      pileup=self.data[chrom][1][summit_index],  
                      pscore=get_pscore(self.data[chrom][1][summit_index],
                                        self.data[chrom][2][summit_index]),
                      fold_change=(self.data[chrom][1][summit_index] +
                                   self.pseudocount) / (self.data[chrom][2][summit_index] +
                                                        self.pseudocount),
                      qscore=qscore,
                      )
            # start a new peak
            return True

    @cython.cfunc
    def __close_peak2(self,
                      peak_content: list,
                      peaks: object,
                      min_length: cython.int,
                      chrom: bytes,
                      smoothlen: cython.int = 51,
                      min_valley: cython.float = 0.9) -> bool:
        """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object.

        In this function, we use signal processing methods to smooth
        the scores in the peak region, find the maxima and enforce the
        peaky shape, and to define the best maxima as the peak
        summit. The functions used for signal processing is 'maxima'
        (with 2nd order polynomial filter) and 'enfoce_peakyness'
        functions in SignalProcessing.pyx.

        peaks: a PeakIO object

        """
        tstart: cython.int
        tend: cython.int
        tmpindex: cython.int
        summit_index: cython.int
        summit_offset: cython.int
        start: cython.int
        end: cython.int
        i: cython.int
        j: cython.int
        start_boundary: cython.int
        tvalue: cython.float
        peakdata: cnp.ndarray(cython.float, ndim=1)
        peakindices: cnp.ndarray(cython.int, ndim=1)
        summit_offsets: cnp.ndarray(cython.int, ndim=1)

        # Add 10 bp padding to peak region so that we can get true minima
        end = peak_content[-1][1] + 10
        start = peak_content[0][0] - 10
        if start < 0:
            start_boundary = 10 + start
            start = 0
        else:
            start_boundary = 10
        peak_length = end - start
        if end - start < min_length:
            return             # if the region is too small, reject it

        peakdata = np.zeros(end - start, dtype='f4')
        peakindices = np.zeros(end - start, dtype='i4')
        for (tstart, tend, tvalue, tsvalue, tmpindex) in peak_content:
            i = tstart - start + start_boundary
            j = tend - start + start_boundary
            peakdata[i:j] = tsvalue
            peakindices[i:j] = tmpindex
        summit_offsets = maxima(peakdata, smoothlen)
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak(peak_content, peaks, min_length, chrom)
        else:
            # remove maxima that occurred in padding
            i = np.searchsorted(summit_offsets,
                                start_boundary)
            j = np.searchsorted(summit_offsets,
                                peak_length + start_boundary,
                                'right')
            summit_offsets = summit_offsets[i:j]

        summit_offsets = enforce_peakyness(peakdata, summit_offsets)
        if summit_offsets.shape[0] == 0:
            # **failsafe** if no summits, fall back on old approach #
            return self.__close_peak(peak_content, peaks, min_length, chrom)

        summit_indices = peakindices[summit_offsets]
        summit_offsets -= start_boundary

        peak_scores = self.data[chrom][3][summit_indices]
        if not (peak_scores > self.cutoff).all():
            return self.__close_peak(peak_content, peaks, min_length, chrom)
        for summit_offset, summit_index in zip(summit_offsets, summit_indices):
            if self.scoring_method == ord('q'):
                qscore = self.data[chrom][3][summit_index]
            else:
                # if q value is not computed, use -1
                qscore = -1
            peaks.add(chrom,
                      start,
                      end,
                      summit=start + summit_offset,
                      peak_score=self.data[chrom][3][summit_index],
                      # should be the same as summit_value
                      pileup=self.data[chrom][1][summit_index],
                      pscore=get_pscore(self.data[chrom][1][summit_index],
                                        self.data[chrom][2][summit_index]),
                      fold_change=(self.data[chrom][1][summit_index] +
                                   self.pseudocount) / (self.data[chrom][2][summit_index] +
                                                        self.pseudocount),
                      qscore=qscore,
                      )
        # start a new peak
        return True

    @cython.cfunc
    def total(self) -> cython.long:
        """Return the number of regions in this object.

        """
        t: cython.long
        chrom: bytes

        t = 0
        for chrom in sorted(self.data.keys()):
            t += self.datalength[chrom]
        return t

    @cython.ccall
    def call_broadpeaks(self,
                        lvl1_cutoff: cython.float = 5.0,
                        lvl2_cutoff: cython.float = 1.0,
                        min_length: cython.int = 200,
                        lvl1_max_gap: cython.int = 50,
                        lvl2_max_gap: cython.int = 400):
        """Return broad peaks constructed from high- and low-cutoff segments.

        Args:
            lvl1_cutoff: Threshold for core enriched segments.
            lvl2_cutoff: Threshold for linking segments.
            min_length: Minimum peak length to report.
            lvl1_max_gap: Maximum gap when merging level-1 segments.
            lvl2_max_gap: Maximum allowed length for linking segments.
        """
        i: cython.int
        chrom: bytes

        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."
        lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff,
                                     min_length=min_length,
                                     max_gap=lvl1_max_gap)
        lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff,
                                     min_length=min_length,
                                     max_gap=lvl2_max_gap)
        chrs = lvl1_peaks.peaks.keys()
        broadpeaks = BroadPeakIO()
        # use lvl2_peaks as linking regions between lvl1_peaks
        for chrom in sorted(chrs):
            lvl1peakschrom = lvl1_peaks.peaks[chrom]
            lvl2peakschrom = lvl2_peaks.peaks[chrom]
            lvl1peakschrom_next = iter(lvl1peakschrom).__next__
            tmppeakset = []             # to temporarily store lvl1 region inside a lvl2 region
            # our assumption is lvl1 regions should be included in lvl2 regions
            try:
                lvl1 = lvl1peakschrom_next()
                for i in range(len(lvl2peakschrom)):
                    # for each lvl2 peak, find all lvl1 peaks inside
                    # I assume lvl1 peaks can be ALL covered by lvl2 peaks.
                    lvl2 = lvl2peakschrom[i]

                    while True:
                        if lvl2["start"] <= lvl1["start"] and lvl1["end"] <= lvl2["end"]:
                            tmppeakset.append(lvl1)
                            lvl1 = lvl1peakschrom_next()
                        else:
                            # make a hierarchical broad peak
                            #print lvl2["start"], lvl2["end"], lvl2["score"]
                            self.__add_broadpeak(broadpeaks,
                                                 chrom,
                                                 lvl2,
                                                 tmppeakset)
                            tmppeakset = []
                            break
            except StopIteration:
                # no more strong (aka lvl1) peaks left
                self.__add_broadpeak(broadpeaks,
                                     chrom,
                                     lvl2,
                                     tmppeakset)
                tmppeakset = []
                # add the rest lvl2 peaks
                for j in range(i+1, len(lvl2peakschrom)):
                    self.__add_broadpeak(broadpeaks,
                                         chrom,
                                         lvl2peakschrom[j],
                                         tmppeakset)

        return broadpeaks

    def __add_broadpeak(self,
                        bpeaks,
                        chrom: bytes,
                        lvl2peak: dict,
                        lvl1peakset: list):
        """Internal function to create broad peak.
        """

        blockNum: cython.int
        thickStart: cython.int
        thickEnd: cython.int
        start: cython.int
        end: cython.int
        blockSizes: bytes
        blockStarts: bytes

        start = lvl2peak["start"]
        end = lvl2peak["end"]

        # the following code will add those broad/lvl2 peaks with no strong/lvl1 peaks inside
        if not lvl1peakset:
            # will complement by adding 1bps start and end to this region
            # may change in the future if gappedPeak format was improved.
            bpeaks.add(chrom,
                       start,
                       end,
                       score=lvl2peak["score"],
                       thickStart=(b"%d" % start),
                       thickEnd=(b"%d" % end),
                       blockNum=2,
                       blockSizes=b"1,1",
                       blockStarts=(b"0,%d" % (end-start-1)),
                       pileup=lvl2peak["pileup"],
                       pscore=lvl2peak["pscore"],
                       fold_change=lvl2peak["fc"],
                       qscore=lvl2peak["qscore"])
            return bpeaks

        thickStart = b"%d" % lvl1peakset[0]["start"]
        thickEnd = b"%d" % lvl1peakset[-1]["end"]
        blockNum = int(len(lvl1peakset))
        blockSizes = b",".join([b"%d" % x["length"] for x in lvl1peakset])
        blockStarts = b",".join([b"%d" % (x["start"]-start) for x in lvl1peakset])

        if lvl2peak["start"] != thickStart:
            # add 1bp mark for the start of lvl2 peak
            thickStart = b"%d" % start
            blockNum += 1
            blockSizes = b"1,"+blockSizes
            blockStarts = b"0,"+blockStarts
        if lvl2peak["end"] != thickEnd:
            # add 1bp mark for the end of lvl2 peak
            thickEnd = b"%d" % end
            blockNum += 1
            blockSizes = blockSizes+b",1"
            blockStarts = blockStarts + b"," + (b"%d" % (end-start-1))

        # add to BroadPeakIO object
        bpeaks.add(chrom,
                   start,
                   end,
                   score=lvl2peak["score"],
                   thickStart=thickStart,
                   thickEnd=thickEnd,
                   blockNum=blockNum,
                   blockSizes=blockSizes,
                   blockStarts=blockStarts,
                   pileup=lvl2peak["pileup"],
                   pscore=lvl2peak["pscore"],
                   fold_change=lvl2peak["fc"],
                   qscore=lvl2peak["qscore"])
        return bpeaks

# MLX-backed implementation (optional import)
try:  # pragma: no cover - optional dependency bridge
    from MACS3.Signal.ScoreTrackMLX import ScoreTrackMLX  # type: ignore # noqa: E402,F401
except Exception:  # pragma: no cover - ignore when MLX is unavailable
    ScoreTrackMLX = None

@cython.cclass
class TwoConditionScores:
    """Class for saving two condition comparison scores.
    """
    # dictionary for data of each chromosome
    data: dict
    # length of data array of each chromosome
    datalength: dict
    # factor to apply to cond1 pileup values
    cond1_factor: cython.float
    # factor to apply to cond2 pileup values
    cond2_factor: cython.float
    # the pseudocount used to calcuate LLR
    pseudocount: cython.float
    cutoff: cython.float
    t1bdg: object
    c1bdg: object
    t2bdg: object
    c2bdg: object
    pvalue_stat1: dict
    pvalue_stat2: dict
    pvalue_stat3: dict

    def __init__(self,
                 t1bdg,
                 c1bdg,
                 t2bdg,
                 c2bdg,
                 cond1_factor: cython.float = 1.0,
                 cond2_factor: cython.float = 1.0,
                 pseudocount: cython.float = 0.01,
                 proportion_background_empirical_distribution: cython.float = 0.99999):
        """Initialise a differential score track from paired bedGraphs.

        Args:
            t1bdg: Treatment bedGraph for condition 1.
            c1bdg: Control bedGraph for condition 1.
            t2bdg: Treatment bedGraph for condition 2.
            c2bdg: Control bedGraph for condition 2.
            cond1_factor: Scaling factor applied to condition 1 tracks.
            cond2_factor: Scaling factor applied to condition 2 tracks.
            pseudocount: Pseudocount added when comparing conditions.
            proportion_background_empirical_distribution: Fraction of bases
                treated as background when building empirical distributions.
        """
        # for each chromosome, there is a l*4 matrix. First column: end
        # position of a region; Second: treatment pileup; third:
        # control pileup ; forth: score (can be p/q-value/likelihood
        # ratio/fold-enrichment/subtraction depending on -c setting)
        self.data = {}
        self.datalength = {}
        self.cond1_factor = cond1_factor
        self.cond2_factor = cond2_factor
        self.pseudocount = pseudocount
        self.pvalue_stat1 = {}
        self.pvalue_stat2 = {}
        self.t1bdg = t1bdg
        self.c1bdg = c1bdg
        self.t2bdg = t2bdg
        self.c2bdg = c2bdg
        # self.empirical_distr_llr = [] # save all values in histogram

    @cython.ccall
    def set_pseudocount(self, pseudocount: cython.float):
        """Update the pseudocount used for differential scoring."""
        self.pseudocount = pseudocount

    @cython.ccall
    def build(self):
        """Compute scores from 3 types of comparisons and store them
        in self.data.

        """
        common_chrs: set
        chrname: bytes
        chrom_max_len: cython.int
        # common chromosome names
        common_chrs = self.get_common_chrs()
        for chrname in common_chrs:
            (cond1_treat_ps, cond1_treat_vs) = self.t1bdg.get_data_by_chr(chrname)
            (cond1_control_ps, cond1_control_vs) = self.c1bdg.get_data_by_chr(chrname)
            (cond2_treat_ps, cond2_treat_vs) = self.t2bdg.get_data_by_chr(chrname)
            (cond2_control_ps, cond2_control_vs) = self.c2bdg.get_data_by_chr(chrname)
            chrom_max_len = len(cond1_treat_ps) + len(cond1_control_ps) + len(cond2_treat_ps) + len(cond2_control_ps)
            self.add_chromosome(chrname, chrom_max_len)
            self.build_chromosome(chrname,
                                  cond1_treat_ps, cond1_control_ps,
                                  cond2_treat_ps, cond2_control_ps,
                                  cond1_treat_vs, cond1_control_vs,
                                  cond2_treat_vs, cond2_control_vs)

    @cython.cfunc
    def build_chromosome(self, chrname,
                         cond1_treat_ps, cond1_control_ps,
                         cond2_treat_ps, cond2_control_ps,
                         cond1_treat_vs, cond1_control_vs,
                         cond2_treat_vs, cond2_control_vs):
        """Internal function to calculate scores for three types of comparisons.

        cond1_treat_ps, cond1_control_ps: position of treat and control of condition 1
        cond2_treat_ps, cond2_control_ps: position of treat and control of condition 2
        cond1_treat_vs, cond1_control_vs: value of treat and control of condition 1
        cond2_treat_vs, cond2_control_vs: value of treat and control of condition 2

        """
        c1tp: cython.int
        c1cp: cython.int
        c2tp: cython.int
        c2cp: cython.int
        minp: cython.int
        pre_p: cython.int
        c1tv: cython.float
        c1cv: cython.float
        c2tv: cython.float
        c2cv: cython.float

        c1tpn = iter(cond1_treat_ps).__next__
        c1cpn = iter(cond1_control_ps).__next__
        c2tpn = iter(cond2_treat_ps).__next__
        c2cpn = iter(cond2_control_ps).__next__
        c1tvn = iter(cond1_treat_vs).__next__
        c1cvn = iter(cond1_control_vs).__next__
        c2tvn = iter(cond2_treat_vs).__next__
        c2cvn = iter(cond2_control_vs).__next__

        pre_p = 0

        try:
            c1tp = c1tpn()
            c1tv = c1tvn()

            c1cp = c1cpn()
            c1cv = c1cvn()

            c2tp = c2tpn()
            c2tv = c2tvn()

            c2cp = c2cpn()
            c2cv = c2cvn()

            while True:
                minp = min(c1tp, c1cp, c2tp, c2cp)
                self.add(chrname, pre_p, c1tv, c1cv, c2tv, c2cv)
                pre_p = minp
                if c1tp == minp:
                    c1tp = c1tpn()
                    c1tv = c1tvn()
                if c1cp == minp:
                    c1cp = c1cpn()
                    c1cv = c1cvn()
                if c2tp == minp:
                    c2tp = c2tpn()
                    c2tv = c2tvn()
                if c2cp == minp:
                    c2cp = c2cpn()
                    c2cv = c2cvn()
        except StopIteration:
            # meet the end of either bedGraphTrackI, simply exit
            pass
        return

    @cython.cfunc
    def get_common_chrs(self) -> set:
        """Return chromosome names shared across all input bedGraphs."""
        t1chrs: set
        c1chrs: set
        t2chrs: set
        c2chrs: set
        common: set
        t1chrs = self.t1bdg.get_chr_names()
        c1chrs = self.c1bdg.get_chr_names()
        t2chrs = self.t2bdg.get_chr_names()
        c2chrs = self.c2bdg.get_chr_names()
        common = reduce(lambda x, y: x.intersection(y),
                        (t1chrs, c1chrs, t2chrs, c2chrs))
        return common

    @cython.cfunc
    def add_chromosome(self,
                       chrom: bytes,
                       chrom_max_len: cython.int = 200000000):
        """Allocate storage for ``chrom`` with capacity ``chrom_max_len``."""
        if chrom not in self.data:
            self.data[chrom] = [np.zeros(chrom_max_len, dtype="i4"),  # pos
                                np.zeros(chrom_max_len, dtype="f4"),  # LLR t1 vs c1
                                np.zeros(chrom_max_len, dtype="f4"),  # LLR t2 vs c2
                                np.zeros(chrom_max_len, dtype="f4")]  # LLR t1 vs t2
            self.datalength[chrom] = 0

    @cython.cfunc
    def add(self,
            chromosome: bytes,
            endpos: cython.int,
            t1: cython.float,
            c1: cython.float,
            t2: cython.float,
            c2: cython.float):
        """Take chr-endpos-sample1-control1-sample2-control2 and
        compute logLR for t1 vs c1, t2 vs c2, and t1 vs t2, then save
        values.

        chromosome: chromosome name in string
        endpos    : end position of each interval in integer
        t1        : Sample 1 ChIP pileup value of each interval in float
        c1        : Sample 1 Control pileup value of each interval in float
        t2        : Sample 2 ChIP pileup value of each interval in float
        c2        : Sample 2 Control pileup value of each interval in float

        *Warning* Need to add regions continuously.
        """
        i: cython.int
        c: list

        i = self.datalength[chromosome]
        c = self.data[chromosome]
        c[0][i] = endpos
        c[1][i] = logLR_asym((t1+self.pseudocount) * self.cond1_factor,
                             (c1+self.pseudocount) * self.cond1_factor)
        c[2][i] = logLR_asym((t2+self.pseudocount) * self.cond2_factor,
                             (c2+self.pseudocount) * self.cond2_factor)
        c[3][i] = logLR_sym((t1+self.pseudocount) * self.cond1_factor,
                            (t2+self.pseudocount) * self.cond2_factor)
        self.datalength[chromosome] += 1
        return

    @cython.ccall
    def finalize(self):
        """
        Adjust array size of each chromosome.

        """
        chrom: bytes
        ln: cython.int
        d: list

        for chrom in sorted(self.data.keys()):
            d = self.data[chrom]
            ln = self.datalength[chrom]
            d[0].resize(ln, refcheck=False)
            d[1].resize(ln, refcheck=False)
            d[2].resize(ln, refcheck=False)
            d[3].resize(ln, refcheck=False)
        return

    @cython.ccall
    def get_data_by_chr(self,
                        chromosome: bytes):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if chromosome in self.data:
            return self.data[chromosome]
        else:
            return None

    @cython.ccall
    def get_chr_names(self):
        """Return all the chromosome names stored.

        """
        return set(self.data.keys())

    @cython.ccall
    def write_bedGraph(self,
                       fhd,
                       name: str,
                       description: str,
                       column: cython.int = 3):
        """Write all data to fhd in bedGraph Format.

        fhd: a filehandler to save bedGraph.

        name/description: the name and description in track line.

        colname: can be 1: cond1 chip vs cond1 ctrl, 2: cond2 chip vs
        cond2 ctrl, 3: cond1 chip vs cond2 chip

        """
        chrom: bytes
        ln: cython.int
        pre: cython.int
        i: cython.int
        p: cython.int
        pre_v: cython.float
        v: cython.float
        pos: cnp.ndarray
        value: cnp.ndarray

        assert column in range(1, 4), "column should be between 1, 2 or 3."

        write = fhd.write

        # if self.trackline:
        #    # this line is REQUIRED by the wiggle format for UCSC browser
        #    write("track type=bedGraph name=\"%s\" description=\"%s\"\n" % (name.decode(), description))

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            pos = self.data[chrom][0]
            value = self.data[chrom][column]
            ln = self.datalength[chrom]
            pre = 0
            if pos.shape[0] == 0:
                continue        # skip if there's no data
            pre_v = value[0]
            for i in range(1, ln):
                v = value[i]
                p = pos[i-1]
                if abs(pre_v - v) >= 1e-6:
                    write("%s\t%d\t%d\t%.5f\n" %
                          (chrom.decode(), pre, p, pre_v))
                    pre_v = v
                    pre = p
            p = pos[-1]
            # last one
            write("%s\t%d\t%d\t%.5f\n" % (chrom.decode(), pre, p, pre_v))

        return True

    @cython.ccall
    def write_matrix(self,
                     fhd,
                     name: str,
                     description: str):
        """Write all data to fhd into five columns Format:

        col1: chr_start_end
        col2: t1 vs c1
        col3: t2 vs c2
        col4: t1 vs t2

        fhd: a filehandler to save the matrix.

        """
        chrom: bytes
        ln: cython.int
        pre: cython.int
        i: cython.int
        p: cython.int
        v1: cython.float
        v2: cython.float
        v3: cython.float
        pos: cnp.ndarray
        value1: cnp.ndarray
        value2: cnp.ndarray
        value3: cnp.ndarray

        write = fhd.write

        chrs = self.get_chr_names()
        for chrom in sorted(chrs):
            [pos, value1, value2, value3] = self.data[chrom]
            ln = self.datalength[chrom]
            pre = 0
            if pos.shape[0] == 0:
                continue        # skip if there's no data
            for i in range(0, ln):
                v1 = value1[i]
                v2 = value2[i]
                v3 = value3[i]
                p = pos[i]
                write("%s:%d_%d\t%.5f\t%.5f\t%.5f\n" %
                      (chrom.decode(), pre, p, v1, v2, v3))
                pre = p

        return True

    @cython.ccall
    def call_peaks(self,
                   cutoff: cython.float = 3,
                   min_length: cython.int = 200,
                   max_gap: cython.int = 100,
                   call_summits: bool = False) -> tuple:
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        For bdgdiff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 3. For log10 LR, it means 1000 or -1000.
        min_length :  minimum peak length, default 200.
        max_gap   :  maximum gap to merge nearby peaks, default 100.
        ptrack:  an optional track for pileup heights. If it's not None, use it to find summits. Otherwise, use self/scoreTrack.
        """
        chrom: bytes
        pos: cnp.ndarray
        t1_vs_c1: cnp.ndarray
        t2_vs_c2: cnp.ndarray
        t1_vs_t2: cnp.ndarray
        cond1_over_cond2: cnp.ndarray
        cond2_over_cond1: cnp.ndarray
        cond1_equal_cond2: cnp.ndarray
        cond1_sig: cnp.ndarray
        cond2_sig: cnp.ndarray
        cat1: cnp.ndarray
        cat2: cnp.ndarray
        cat3: cnp.ndarray
        cat1_startpos: cnp.ndarray
        cat1_endpos: cnp.ndarray
        cat2_startpos: cnp.ndarray
        cat2_endpos: cnp.ndarray
        cat3_startpos: cnp.ndarray
        cat3_endpos: cnp.ndarray

        chrs = self.get_chr_names()
        cat1_peaks = PeakIO()       # dictionary to save peaks significant at condition 1
        cat2_peaks = PeakIO()       # dictionary to save peaks significant at condition 2
        cat3_peaks = PeakIO()       # dictionary to save peaks significant in both conditions

        self.cutoff = cutoff

        for chrom in sorted(chrs):
            pos = self.data[chrom][0]
            t1_vs_c1 = self.data[chrom][1]
            t2_vs_c2 = self.data[chrom][2]
            t1_vs_t2 = self.data[chrom][3]
            and_ = np.logical_and
            # regions with stronger cond1 signals
            cond1_over_cond2 = t1_vs_t2 >= cutoff
            # regions with stronger cond2 signals
            cond2_over_cond1 = t1_vs_t2 <= -1*cutoff
            cond1_equal_cond2 = and_(t1_vs_t2 >= -1*cutoff, t1_vs_t2 <= cutoff)
            # enriched regions in condition 1
            cond1_sig = t1_vs_c1 >= cutoff
            # enriched regions in condition 2
            cond2_sig = t2_vs_c2 >= cutoff
            # indices where score is above cutoff
            # cond1 stronger than cond2, the indices
            cat1 = np.where(and_(cond1_sig, cond1_over_cond2))[0]
            # cond2 stronger than cond1, the indices
            cat2 = np.where(and_(cond2_over_cond1, cond2_sig))[0]
            # cond1 and cond2 are equal, the indices
            cat3 = np.where(and_(and_(cond1_sig, cond2_sig),
                                 cond1_equal_cond2))[0]

            # end positions of regions where score is above cutoff
            cat1_endpos = pos[cat1]
            # start positions of regions where score is above cutoff            
            cat1_startpos = pos[cat1-1]
            # end positions of regions where score is above cutoff
            cat2_endpos = pos[cat2]
            # start positions of regions where score is above cutoff
            cat2_startpos = pos[cat2-1]
            # end positions of regions where score is above cutoff
            cat3_endpos = pos[cat3]
            # start positions of regions where score is above cutoff
            cat3_startpos = pos[cat3-1]

            # for cat1: condition 1 stronger regions
            self.__add_a_peak(cat1_peaks,
                              chrom,
                              cat1,
                              cat1_startpos,
                              cat1_endpos,
                              t1_vs_t2,
                              max_gap,
                              min_length)
            # for cat2: condition 2 stronger regions
            self.__add_a_peak(cat2_peaks,
                              chrom,
                              cat2,
                              cat2_startpos,
                              cat2_endpos,
                              -1 * t1_vs_t2,
                              max_gap,
                              min_length)
            # for cat3: commonly strong regions
            self.__add_a_peak(cat3_peaks,
                              chrom,
                              cat3,
                              cat3_startpos,
                              cat3_endpos,
                              abs(t1_vs_t2),
                              max_gap,
                              min_length)

        return (cat1_peaks, cat2_peaks, cat3_peaks)

    @cython.cfunc
    def __add_a_peak(self,
                     peaks: object,
                     chrom: bytes,
                     indices: cnp.ndarray,
                     startpos: cnp.ndarray,
                     endpos: cnp.ndarray,
                     score: cnp.ndarray,
                     max_gap: cython.int,
                     min_length: cython.int):
        
        """For a given chromosome, merge nearby significant regions,
        filter out smaller regions, then add regions to PeakIO
        object.

        """
        i: cython.int
        peak_content: list
        mean_logLR: cython.float

        if startpos.size > 0:
            # if it is not empty
            peak_content = []
            if indices[0] == 0:
                # first element > cutoff, fix the first point as
                # 0. otherwise it would be the last item in
                # data[chrom]['pos']
                startpos[0] = 0
            # first bit of region above cutoff
            peak_content.append((startpos[0],
                                 endpos[0],
                                 score[indices[0]]))
            for i in range(1, startpos.size):
                if startpos[i] - peak_content[-1][1] <= max_gap:
                    # append
                    peak_content.append((startpos[i],
                                         endpos[i],
                                         score[indices[i]]))
                else:
                    # close
                    if peak_content[-1][1] - peak_content[0][0] >= min_length:
                        mean_logLR = self.mean_from_peakcontent(peak_content)
                        # if peak_content[0][0] == 22414956:
                        #    print(f"{peak_content} {mean_logLR}")
                        peaks.add(chrom,
                                  peak_content[0][0],
                                  peak_content[-1][1],
                                  summit=-1,
                                  peak_score=mean_logLR,
                                  pileup=0,
                                  pscore=0,
                                  fold_change=0,
                                  qscore=0,
                                  )
                    peak_content = [(startpos[i],
                                     endpos[i],
                                     score[indices[i]]),]

            # save the last peak
            if peak_content:
                if peak_content[-1][1] - peak_content[0][0] >= min_length:
                    mean_logLR = self.mean_from_peakcontent(peak_content)
                    peaks.add(chrom,
                              peak_content[0][0],
                              peak_content[-1][1],
                              summit=-1,
                              peak_score=mean_logLR,
                              pileup=0,
                              pscore=0,
                              fold_change=0,
                              qscore=0,
                              )

        return

    @cython.cfunc
    def mean_from_peakcontent(self,
                              peakcontent: list) -> cython.float:
        """

        """
        tmp_s: cython.int
        tmp_e: cython.int
        ln: cython.int
        tmp_v: cython.long
        sum_v: cython.long      # for better precision
        r: cython.float
        i: cython.int

        ln = 0
        sum_v = 0                         # initialize sum_v as 0
        for i in range(len(peakcontent)):
            tmp_s = peakcontent[i][0]
            tmp_e = peakcontent[i][1]
            tmp_v = peakcontent[i][2]
            sum_v += tmp_v * (tmp_e - tmp_s)
            ln += tmp_e - tmp_s

        r = cython.cast(cython.float, (sum_v / ln))
        return r

    @cython.cfunc
    def total(self) -> cython.long:
        """Return the number of regions in this object.

        """
        t: cython.long
        chrom: bytes

        t = 0
        for chrom in sorted(self.data.keys()):
            t += self.datalength[chrom]
        return t
