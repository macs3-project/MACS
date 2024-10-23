# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-22 16:59:53 Tao Liu>

"""Module for SAPPER PosReadsInfo class.

Copyright (c) 2017 Tao Liu <tliu4@buffalo.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
from MACS3.Signal.VariantStat import (CalModel_Homo,
                                      CalModel_Heter_noAS,
                                      CalModel_Heter_AS)
# calculate_GQ,
# calculate_GQ_heterASsig)
from MACS3.Signal.Prob import binomial_cdf
from MACS3.Signal.PeakVariants import Variant

import cython
import numpy as np
import cython.cimports.numpy as cnp
from cython.cimports.cpython import bool

LN10 = 2.3025850929940458

# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <tliu4@buffalo.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class PosReadsInfo:
    ref_pos: cython.long
    ref_allele: cython.bytes
    alt_allele: cython.bytes
    filterout: bool          # if true, do not output
    bq_set_T: dict       # {A:[], C:[], G:[], T:[], N:[]} for treatment
    bq_set_C: dict
    n_reads_T: dict      # {A:[], C:[], G:[], T:[], N:[]} for treatment
    n_reads_C: dict
    n_reads: dict
    n_strand: list  # [{A:[], C:[], G:[], T:[], N:[]},{A:[], C:[], G:[], T:[], N:[]}] for total appearance on plus strand and minus strand for ChIP sample only
    n_tips: dict  # count of nt appearing at tips
    top1allele: cython.bytes
    top2allele: cython.bytes
    top12alleles_ratio: cython.float
    lnL_homo_major: cython.double
    lnL_heter_AS: cython.double
    lnL_heter_noAS: cython.double
    lnL_homo_minor: cython.double
    BIC_homo_major: cython.double
    BIC_heter_AS: cython.double
    BIC_heter_noAS: cython.double
    BIC_homo_minor: cython.double
    PL_00: cython.double
    PL_01: cython.double
    PL_11: cython.double
    deltaBIC: cython.double
    heter_noAS_kc: cython.int
    heter_noAS_ki: cython.int
    heter_AS_kc: cython.int
    heter_AS_ki: cython.int
    heter_AS_alleleratio: cython.double
    GQ_homo_major: cython.int
    GQ_heter_noAS: cython.int
    GQ_heter_AS: cython.int  # phred scale of prob by standard formular
    GQ_heter_ASsig: cython.int  # phred scale of prob, to measure the difference between AS and noAS
    GQ: cython.double
    GT: str
    type: str
    mutation_type: str       # SNV or Insertion or Deletion
    hasfermiinfor: bool  # if no fermi bam overlap in the position, false; if fermi bam in the position GT: N, false; if anyone of top2allele is not in fermi GT NTs, false;
    fermiNTs: bytearray

    def __cinit__(self):
        self.filterout = False
        self.GQ = 0
        self.GT = "unsure"
        self.alt_allele = b'.'

    def __init__(self,
                 ref_pos: cython.long,
                 ref_allele: cython.bytes):
        self.ref_pos = ref_pos
        self.ref_allele = ref_allele
        self.bq_set_T = {ref_allele: [], b'A': [], b'C': [], b'G': [], b'T': [], b'N': [], b'*': []}
        self.bq_set_C = {ref_allele: [], b'A': [], b'C': [], b'G': [], b'T': [], b'N': [], b'*': []}
        self.n_reads_T = {ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0}
        self.n_reads_C = {ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0}
        self.n_reads = {ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0}
        self.n_strand = [{ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0},
                         {ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0}]
        self.n_tips = {ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0}

    def __getstate__(self):
        return (self.ref_pos,
                self.ref_allele,
                self.alt_allele,
                self.filterout,
                self.bq_set_T,
                self.bq_set_C,
                self.n_reads_T,
                self.n_reads_C,
                self.n_reads,
                self.n_strand,
                self.n_tips,
                self.top1allele,
                self.top2allele,
                self.top12alleles_ratio,
                self.lnL_homo_major,
                self.lnL_heter_AS,
                self.lnL_heter_noAS,
                self.lnL_homo_minor,
                self.BIC_homo_major,
                self.BIC_heter_AS,
                self.BIC_heter_noAS,
                self.BIC_homo_minor,
                self.heter_noAS_kc,
                self.heter_noAS_ki,
                self.heter_AS_kc,
                self.heter_AS_ki,
                self.heter_AS_alleleratio,
                self.GQ_homo_major,
                self.GQ_heter_noAS,
                self.GQ_heter_AS,
                self.GQ_heter_ASsig,
                self.GQ,
                self.GT,
                self.type,
                self.hasfermiinfor,
                self.fermiNTs)

    def __setstate__(self, state):
        (self.ref_pos,
         self.ref_allele,
         self.alt_allele,
         self.filterout,
         self.bq_set_T,
         self.bq_set_C,
         self.n_reads_T,
         self.n_reads_C,
         self.n_reads,
         self.n_strand,
         self.n_tips,
         self.top1allele,
         self.top2allele,
         self.top12alleles_ratio,
         self.lnL_homo_major,
         self.lnL_heter_AS,
         self.lnL_heter_noAS,
         self.lnL_homo_minor,
         self.BIC_homo_major,
         self.BIC_heter_AS,
         self.BIC_heter_noAS,
         self.BIC_homo_minor,
         self.heter_noAS_kc,
         self.heter_noAS_ki,
         self.heter_AS_kc,
         self.heter_AS_ki,
         self.heter_AS_alleleratio,
         self.GQ_homo_major,
         self.GQ_heter_noAS,
         self.GQ_heter_AS,
         self.GQ_heter_ASsig,
         self.GQ,
         self.GT,
         self.type,
         self.hasfermiinfor,
         self.fermiNTs) = state

    @cython.ccall
    def filterflag(self) -> bool:
        return self.filterout

    @cython.ccall
    def apply_GQ_cutoff(self,
                        min_homo_GQ: cython.int = 50,
                        min_heter_GQ: cython.int = 100):
        if self.filterout:
            return
        if self.type.startswith('homo') and self.GQ < min_homo_GQ:
            self.filterout = True
        elif self.type.startswith('heter') and self.GQ < min_heter_GQ:
            self.filterout = True
        return

    @cython.ccall
    def apply_deltaBIC_cutoff(self,
                              min_delta_BIC: cython.float = 10):
        if self.filterout:
            return
        if self.deltaBIC < min_delta_BIC:
            self.filterout = True
        return

    @cython.ccall
    def add_T(self,
              read_index: cython.int,
              read_allele: cython.bytes,
              read_bq: cython.int,
              strand: cython.int,
              tip: bool,
              Q: cython.int = 20):
        """ Strand 0: plus, 1: minus

        Q is the quality cutoff. By default, only consider Q20 or read_bq > 20.
        """
        if read_bq <= Q:
            return
        if not self.n_reads.has_key(read_allele):
            self.bq_set_T[read_allele] = []
            self.bq_set_C[read_allele] = []
            self.n_reads_T[read_allele] = 0
            self.n_reads_C[read_allele] = 0
            self.n_reads[read_allele] = 0
            self.n_strand[0][read_allele] = 0
            self.n_strand[1][read_allele] = 0
            self.n_tips[read_allele] = 0
        self.bq_set_T[read_allele].append(read_bq)
        self.n_reads_T[read_allele] += 1
        self.n_reads[read_allele] += 1
        self.n_strand[strand][read_allele] += 1
        if tip:
            self.n_tips[read_allele] += 1

    @cython.ccall
    def add_C(self,
              read_index: cython.int,
              read_allele: cython.bytes,
              read_bq: cython.int,
              strand: cython.int,
              Q: cython.int = 20):
        if read_bq <= Q:
            return
        if not self.n_reads.has_key(read_allele):
            self.bq_set_T[read_allele] = []
            self.bq_set_C[read_allele] = []
            self.n_reads_T[read_allele] = 0
            self.n_reads_C[read_allele] = 0
            self.n_reads[read_allele] = 0
            self.n_strand[0][read_allele] = 0
            self.n_strand[1][read_allele] = 0
            self.n_tips[read_allele] = 0
        self.bq_set_C[read_allele].append(read_bq)
        self.n_reads_C[read_allele] += 1
        self.n_reads[read_allele] += 1

    @cython.ccall
    def raw_read_depth(self,
                       opt: str = "all") -> cython.int:
        if opt == "all":
            return sum(self.n_reads.values())
        elif opt == "T":
            return sum(self.n_reads_T.values())
        elif opt == "C":
            return sum(self.n_reads_C.values())
        else:
            raise Exception("opt should be either 'all', 'T' or 'C'.")

    @cython.ccall
    def update_top_alleles(self,
                           min_top12alleles_ratio: cython.float = 0.8,
                           min_altallele_count: cython.int = 2,
                           max_allowed_ar: cython.float = 0.95):
        """Identify top1 and top2 NT.  the ratio of (top1+top2)/total
        """
        [self.top1allele, self.top2allele] = sorted(self.n_reads,
                                                    key=self.n_reads_T.get,
                                                    reverse=True)[:2]

        # if top2 allele count in ChIP is lower than
        # min_altallele_count, or when allele ratio top1/(top1+top2)
        # is larger than max_allowed_ar in ChIP, we won't consider
        # this allele at all.  we set values of top2 allele in
        # dictionaries to [] and ignore top2 allele entirely.
        # max(self.n_strand[0][self.top2allele], self.n_strand[1][self.top2allele]) < min_altallele_count
        # if self.ref_pos == 52608504:
        #    prself: cython.int.ref_pos, self.n_reads_T[self.top1allele], self.n_reads_T[self.top2allele], self.n_reads_C[self.top1allele], self.n_reads_C[self.top2allele]
        if self.n_reads_T[self.top1allele] + self.n_reads_T[self.top2allele] == 0:
            self.filterout = True
            return

        if (len(self.top1allele) == 1 and len(self.top2allele) == 1) and (self.top2allele != self.ref_allele and ((self.n_reads_T[self.top2allele] - self.n_tips[self.top2allele]) < min_altallele_count) or self.n_reads_T[self.top1allele]/(self.n_reads_T[self.top1allele] + self.n_reads_T[self.top2allele]) > max_allowed_ar):
            self.bq_set_T[self.top2allele] = []
            self.bq_set_C[self.top2allele] = []
            self.n_reads_T[self.top2allele] = 0
            self.n_reads_C[self.top2allele] = 0
            self.n_reads[self.top2allele] = 0
            self.n_tips[self.top2allele] = 0
            if (self.top1allele != self.ref_allele and (self.n_reads_T[self.top1allele] - self.n_tips[self.top1allele]) < min_altallele_count):
                self.bq_set_T[self.top1allele] = []
                self.bq_set_C[self.top1allele] = []
                self.n_reads_T[self.top1allele] = 0
                self.n_reads_C[self.top1allele] = 0
                self.n_reads[self.top1allele] = 0
                self.n_tips[self.top1allele] = 0

        if self.n_reads_T[self.top1allele] + self.n_reads_T[self.top2allele] == 0:
            self.filterout = True
            return

        self.top12alleles_ratio = (self.n_reads[self.top1allele] + self.n_reads[self.top2allele]) / sum(self.n_reads.values())
        if self.top12alleles_ratio < min_top12alleles_ratio:
            self.filterout = True
            return

        if self.top1allele == self.ref_allele and self.n_reads[self.top2allele] == 0:
            # This means this position only contains top1allele which is the ref_allele. So the GT must be 0/0
            self.type = "homo_ref"
            self.filterout = True
            return
        return

    @cython.ccall
    def top12alleles(self):
        print(self.ref_pos, self.ref_allele)
        print("Top1allele", self.top1allele, "Treatment",
              self.bq_set_T[self.top1allele], "Control",
              self.bq_set_C[self.top1allele])
        print("Top2allele", self.top2allele, "Treatment",
              self.bq_set_T[self.top2allele], "Control",
              self.bq_set_C[self.top2allele])
    
    @cython.ccall
    def call_GT(self, max_allowed_ar: cython.float = 0.99):
        """Require update_top_alleles being called.
        """
        top1_bq_T: cnp.ndarray(cython.int, ndim=1)
        top2_bq_T: cnp.ndarray(cython.int, ndim=1)
        top1_bq_C: cnp.ndarray(cython.int, ndim=1)
        top2_bq_C: cnp.ndarray(cython.int, ndim=1)
        tmp_mutation_type: list
        tmp_alt: cython.bytes

        if self.filterout:
            return

        top1_bq_T = np.array(self.bq_set_T[self.top1allele], dtype="i4")
        top2_bq_T = np.array(self.bq_set_T[self.top2allele], dtype="i4")
        top1_bq_C = np.array(self.bq_set_C[self.top1allele], dtype="i4")
        top2_bq_C = np.array(self.bq_set_C[self.top2allele], dtype="i4")
        (self.lnL_homo_major, self.BIC_homo_major) = CalModel_Homo(top1_bq_T,
                                                                   top1_bq_C,
                                                                   top2_bq_T,
                                                                   top2_bq_C)
        (self.lnL_homo_minor, self.BIC_homo_minor) = CalModel_Homo(top2_bq_T,
                                                                   top2_bq_C,
                                                                   top1_bq_T,
                                                                   top1_bq_C)
        (self.lnL_heter_noAS, self.BIC_heter_noAS) = CalModel_Heter_noAS(top1_bq_T,
                                                                         top1_bq_C,
                                                                         top2_bq_T,
                                                                         top2_bq_C)
        (self.lnL_heter_AS, self.BIC_heter_AS) = CalModel_Heter_AS(top1_bq_T,
                                                                   top1_bq_C,
                                                                   top2_bq_T,
                                                                   top2_bq_C,
                                                                   max_allowed_ar)

        # if self.ref_pos == 71078525:
        #    print "---"
        #    prlen: cython.int(top1_bq_T), len(top1_bq_C), len(top2_bq_T), len(top2_bq_C)
        #    prself: cython.int.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS, self.lnL_heter_AS
        #    prself: cython.int.BIC_homo_major, self.BIC_homo_minor, self.BIC_heter_noAS, self.BIC_heter_AS

        if self.top1allele != self.ref_allele and self.n_reads[self.top2allele] == 0:
            # in this case, there is no top2 nt (or socalled minor
            # allele) in either treatment or control, we should assume
            # it's a 1/1 genotype. We will take 1/1 if it passes BIC
            # test (deltaBIC >=2), and will skip this loci if it can't
            # pass the test.

            self.deltaBIC = min(self.BIC_heter_noAS, self.BIC_heter_AS, self.BIC_homo_minor) - self.BIC_homo_major
            if self.deltaBIC < 2:
                self.filterout = True
                return

            self.type = "homo"
            self.GT = "1/1"

            self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
            self.PL_01 = -10.0 * max(self.lnL_heter_noAS, self.lnL_heter_AS) / LN10
            self.PL_11 = -10.0 * self.lnL_homo_major / LN10

            self.PL_00 = max(0, self.PL_00 - self.PL_11)
            self.PL_01 = max(0, self.PL_01 - self.PL_11)
            self.PL_11 = 0

            self.GQ = min(self.PL_00, self.PL_01)
            self.alt_allele = self.top1allele
        else:
            # assign GQ, GT, and type
            if self.ref_allele != self.top1allele and self.BIC_homo_major + 2 <= self.BIC_homo_minor and self.BIC_homo_major + 2 <= self.BIC_heter_noAS and self.BIC_homo_major + 2 <= self.BIC_heter_AS:
                self.type = "homo"
                self.deltaBIC = min(self.BIC_heter_noAS, self.BIC_heter_AS, self.BIC_homo_minor) - self.BIC_homo_major
                self.GT = "1/1"
                self.alt_allele = self.top1allele

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * max(self.lnL_heter_noAS, self.lnL_heter_AS) / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_11
                self.PL_01 = self.PL_01 - self.PL_11
                self.PL_11 = 0

                self.GQ = min(self.PL_00, self.PL_01)

            elif self.BIC_heter_noAS + 2 <= self.BIC_homo_major and self.BIC_heter_noAS + 2 <= self.BIC_homo_minor and self.BIC_heter_noAS + 2 <= self.BIC_heter_AS:
                self.type = "heter_noAS"
                self.deltaBIC = min(self.BIC_homo_major, self.BIC_homo_minor) - self.BIC_heter_noAS

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * self.lnL_heter_noAS / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_01
                self.PL_11 = self.PL_11 - self.PL_01
                self.PL_01 = 0

                self.GQ = min(self.PL_00, self.PL_11)

            elif self.BIC_heter_AS + 2 <= self.BIC_homo_major and self.BIC_heter_AS + 2 <= self.BIC_homo_minor and self.BIC_heter_AS + 2 <= self.BIC_heter_noAS:
                self.type = "heter_AS"
                self.deltaBIC = min(self.BIC_homo_major, self.BIC_homo_minor) - self.BIC_heter_AS

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * self.lnL_heter_AS / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_01
                self.PL_11 = self.PL_11 - self.PL_01
                self.PL_01 = 0

                self.GQ = min(self.PL_00, self.PL_11)

            elif self.BIC_heter_AS + 2 <= self.BIC_homo_major and self.BIC_heter_AS + 2 <= self.BIC_homo_minor:
                # can't decide if it's noAS or AS
                self.type = "heter_unsure"
                self.deltaBIC = min(self.BIC_homo_major, self.BIC_homo_minor) - max(self.BIC_heter_AS, self.BIC_heter_noAS)

                self.PL_00 = -10.0 * self.lnL_homo_minor / LN10
                self.PL_01 = -10.0 * max(self.lnL_heter_noAS, self.lnL_heter_AS) / LN10
                self.PL_11 = -10.0 * self.lnL_homo_major / LN10

                self.PL_00 = self.PL_00 - self.PL_01
                self.PL_11 = self.PL_11 - self.PL_01
                self.PL_01 = 0

                self.GQ = min(self.PL_00, self.PL_11)

            elif self.ref_allele == self.top1allele and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
                self.type = "homo_ref"
                # we do not calculate GQ if type is homo_ref
                self.GT = "0/0"
                self.filterout = True
            else:
                self.type = "unsure"
                self.filterout = True

            if self.type.startswith("heter"):
                if self.ref_allele == self.top1allele:
                    self.alt_allele = self.top2allele
                    self.GT = "0/1"
                elif self.ref_allele == self.top2allele:
                    self.alt_allele = self.top1allele
                    self.GT = "0/1"
                else:
                    self.alt_allele = self.top1allele+b','+self.top2allele
                    self.GT = "1/2"

        tmp_mutation_type = []
        for tmp_alt in self.alt_allele.split(b','):
            if tmp_alt == b'*':
                tmp_mutation_type.append("Deletion")
            elif len(tmp_alt) > 1:
                tmp_mutation_type.append("Insertion")
            else:
                tmp_mutation_type.append("SNV")
        self.mutation_type = ",".join(tmp_mutation_type)
        return

    @cython.cfunc
    def SB_score_ChIP(self,
                      a: cython.int,
                      b: cython.int,
                      c: cython.int,
                      d: cython.int) -> cython.float:
        """ calculate score for filtering variants with strange strand biases.

        a: top1/major allele plus strand
        b: top2/minor allele plus strand
        c: top1/major allele minus strand
        d: top2/minor allele minus strand

        Return a value: cython.float so that if this value >= 1, the variant will be filtered out.
        """
        p1_l: cython.double
        p1_r: cython.double
        p2_l: cython.double
        p2_r: cython.double

        if a + b == 0 or c + d == 0:
            # if major allele and minor allele both bias to the same strand, allow it
            return 0.0

        # Rule:
        # if there is bias in top2 allele then bias in top1 allele should not be significantly smaller than it.
        # or there is no significant bias (0.5) in top2 allele.

        # pra: cython.int, b, c, d
        p1_l = binomial_cdf(a, (a+c), 0.5, lower=True)      # alternative: less than 0.5
        p1_r = binomial_cdf(c, (a+c), 0.5, lower=True)                   # greater than 0.5
        p2_l = binomial_cdf(b, (b+d), 0.5, lower=True)      # alternative: less than 0.5
        p2_r = binomial_cdf(d, (b+d), 0.5, lower=True)                   # greater than 0.5
        # prp1_l: cython.int, p1_r, p2_l, p2_r

        if (p1_l < 0.05 and p2_r < 0.05) or (p1_r < 0.05 and p2_l < 0.05):
            # we reject loci where the significant biases are inconsistent between top1 and top2 alleles.
            return 1.0
        else:
            # if b<=2 and d=0 or b=0 and d<=2 -- highly possible FPs
            # if (b<=2 and d==0 or b==0 and d<=2):
            #    return 1
            # can't decide
            return 0.0

    @cython.cfunc
    def SB_score_ATAC(self,
                      a: cython.int,
                      b: cython.int,
                      c: cython.int,
                      d: cython.int) -> cython.float:
        """ calculate score for filtering variants with strange strand biases.

        ATAC-seq version

        a: top1/major allele plus strand
        b: top2/minor allele plus strand
        c: top1/major allele minus strand
        d: top2/minor allele minus strand

        Return a value: cython.float so that if this value >= 1, the variant will be filtered out.
        """
        p1_l: cython.double
        p1_r: cython.double
        p2_l: cython.double
        p2_r: cython.double

        if a+b == 0 or c+d == 0:
            # if major allele and minor allele both bias to the same strand, allow it
            return 0.0

        # Rule:
        # if there is bias in top2 allele then bias in top1 allele should not be significantly smaller than it.
        # or there is no significant bias (0.5) in top2 allele.
        # pra: cython.int, b, c, d
        p1_l = binomial_cdf(a, (a+c), 0.5, lower=True)      # alternative: less than 0.5
        p1_r = binomial_cdf(c, (a+c), 0.5, lower=True)      #              greater than 0.5
        p2_l = binomial_cdf(b, (b+d), 0.5, lower=True)      # alternative: less than 0.5
        p2_r = binomial_cdf(d, (b+d), 0.5, lower=True)      #              greater than 0.5
        # prp1_l: cython.int, p1_r, p2_l, p2_r

        if (p1_l < 0.05 and p2_r < 0.05) or (p1_r < 0.05 and p2_l < 0.05):
            # we reject loci where the significant biases are inconsistent between top1 and top2 alleles.
            return 1.0
        else:
            # can't decide
            return 0.0

    @cython.ccall
    def to_vcf(self) -> str:
        """Output REF,ALT,QUAL,FILTER,INFO,FORMAT, SAMPLE columns.
        """
        vcf_ref: str
        vcf_alt: str
        vcf_qual: str
        vcf_filter: str
        vcf_info: str
        vcf_format: str
        vcf_sample: str

        vcf_ref = self.ref_allele.decode()
        vcf_alt = self.alt_allele.decode()
        vcf_qual = "%d" % self.GQ
        vcf_filter = "."
        vcf_info = (b"M=%s;MT=%s;DPT=%d;DPC=%d;DP1T=%d%s;DP2T=%d%s;DP1C=%d%s;DP2C=%d%s;SB=%d,%d,%d,%d;DBIC=%.2f;BICHOMOMAJOR=%.2f;BICHOMOMINOR=%.2f;BICHETERNOAS=%.2f;BICHETERAS=%.2f;AR=%.2f" %
                    (self.type.encode(),
                     self.mutation_type.encode(),
                     sum(self.n_reads_T.values()),
                     sum(self.n_reads_C.values()),
                     self.n_reads_T[self.top1allele],
                     self.top1allele,
                     self.n_reads_T[self.top2allele],
                     self.top2allele,
                     self.n_reads_C[self.top1allele],
                     self.top1allele,
                     self.n_reads_C[self.top2allele],
                     self.top2allele,
                     self.n_strand[0][self.top1allele],
                     self.n_strand[0][self.top2allele],
                     self.n_strand[1][self.top1allele],
                     self.n_strand[1][self.top2allele],
                     self.deltaBIC,
                     self.BIC_homo_major,
                     self.BIC_homo_minor,
                     self.BIC_heter_noAS,
                     self.BIC_heter_AS,
                     self.n_reads_T[self.top1allele]/(self.n_reads_T[self.top1allele]+self.n_reads_T[self.top2allele])
                    )).decode()
        vcf_format = "GT:DP:GQ:PL"
        vcf_sample = "%s:%d:%d:%d,%d,%d" % (self.GT, self.raw_read_depth(opt="all"), self.GQ, self.PL_00, self.PL_01, self.PL_11)
        return "\t".join((vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample))

    @cython.ccall
    def toVariant(self):
        v: Variant

        v = Variant(self.ref_allele.decode(),
                    self.alt_allele.decode(),
                    self.GQ,
                    '.',
                    self.type,
                    self.mutation_type,
                    self.top1allele.decode(),
                    self.top2allele.decode(),
                    sum(self.n_reads_T.values()),
                    sum(self.n_reads_C.values()),
                    self.n_reads_T[self.top1allele],
                    self.n_reads_T[self.top2allele],
                    self.n_reads_C[self.top1allele],
                    self.n_reads_C[self.top2allele],
                    self.n_strand[0][self.top1allele],
                    self.n_strand[0][self.top2allele],
                    self.n_strand[1][self.top1allele],
                    self.n_strand[1][self.top2allele],
                    self.deltaBIC,
                    self.BIC_homo_major,
                    self.BIC_homo_minor,
                    self.BIC_heter_noAS,
                    self.BIC_heter_AS,
                    self.n_reads_T[self.top1allele]/(self.n_reads_T[self.top1allele]+self.n_reads_T[self.top2allele]),
                    self.GT,
                    self.raw_read_depth(opt="all"),
                    self.PL_00,
                    self.PL_01,
                    self.PL_11)
        return v
