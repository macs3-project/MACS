# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-22 16:26:57 Tao Liu>

"""Module for ReadAlignment collection

Copyright (c) 2024 Tao Liu

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: vladimir.liu@gmail.com
"""
# ------------------------------------
# python modules
# ------------------------------------
from collections import Counter
from operator import itemgetter
from copy import copy

from MACS3.Signal.ReadAlignment import ReadAlignment
from MACS3.Signal.PosReadsInfo import PosReadsInfo
from MACS3.Signal.UnitigRACollection import UnitigRAs, UnitigCollection
from MACS3.IO.PeakIO import PeakIO

import cython
from cython.cimports.cpython import bool
# from cython.cimports.cpython.mem import PyMem_Malloc, PyMem_Free

from cython.cimports.libc.stdlib import malloc, free

# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <tliu4@buffalo.edu>"
__doc__ = "All Parser classes"

__DNACOMPLEMENT__ = b'\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQRSAUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff' # A trans table to convert A to T, C to G, G to C, and T to A.

__CIGARCODE__ = "MIDNSHP=X"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class RACollection:
    """A collection of ReadAlignment objects and the corresponding
    PeakIO.

    """
    chrom: bytes
    peak: PeakIO             # A PeakIO object
    RAlists: list           # contain ReadAlignment lists for treatment (0) and control (1)
    left: cython.long               # left position of peak
    right: cython.long              # right position of peak
    length: cython.long             # length of peak
    RAs_left: cython.long           # left position of all RAs in the collection
    RAs_right: cython.long          # right position of all RAs in the collection
    sorted: bool             # if sorted by lpos
    peak_refseq: bytes       # reference sequence in peak region b/w left and right
    peak_refseq_ext: bytes   # reference sequence in peak region with extension on both sides b/w RAs_left and RAs_right

    def __init__(self, chrom: bytes, peak: PeakIO, RAlist_T: list, RAlist_C: list = []):
        """Create RACollection by: object taking:

        1. peak: a PeakIO indicating: object the peak region.

        2. RAlist: a python of: list ReadAlignment objects containing
        all the reads overlapping the peak region. If no RAlist_C
        given, it will be [].

        """
        if len(RAlist_T) == 0:
            # no reads, return None
            raise Exception("No reads from ChIP sample to construct RAcollection!")
        self.chrom = chrom
        self.peak = peak
        # print(len(RAlist_T),"\n")
        # print(len(RAlist_C),"\n")
        self.RAlists = [RAlist_T, RAlist_C]
        self.left = peak["start"]
        self.right = peak["end"]
        self.length = self.right - self.left
        if RAlist_T:
            self.RAs_left = RAlist_T[0]["lpos"]  # initial assignment of RAs_left
            self.RAs_right = RAlist_T[-1]["rpos"]  # initial assignment of RAs_right
            self.sort()                           # it will set self.sorted = True
        else:
            self.RAs_left = -1
            self.RAs_right = -1
        # check RAs_left and RAs_right
        for ra in RAlist_T:
            if ra["lpos"] < self.RAs_left:
                self.RAs_left = ra["lpos"]
            if ra["rpos"] > self.RAs_right:
                self.RAs_right = ra["rpos"]

        for ra in RAlist_C:
            if ra["lpos"] < self.RAs_left:
                self.RAs_left = ra["lpos"]
            if ra["rpos"] > self.RAs_right:
                self.RAs_right = ra["rpos"]
        (self.peak_refseq, self.peak_refseq_ext) = self.__get_peak_REFSEQ()

    def __getitem__(self, keyname):
        if keyname == "chrom":
            return self.chrom
        elif keyname == "left":
            return self.left
        elif keyname == "right":
            return self.right
        elif keyname == "RAs_left":
            return self.RAs_left
        elif keyname == "RAs_right":
            return self.RAs_right
        elif keyname == "length":
            return self.length
        elif keyname == "count":
            return len(self.RAlists[0]) + len(self.RAlists[1])
        elif keyname == "count_T":
            return len(self.RAlists[0])
        elif keyname == "count_C":
            return len(self.RAlists[1])
        elif keyname == "peak_refseq":
            return self.peak_refseq
        elif keyname == "peak_refseq_ext":
            return self.peak_refseq_ext
        else:
            raise KeyError("Unavailable key:", keyname)

    def __getstate__(self):
        #return {"chrom":self.chrom, "peak":self.peak, "RAlists":self.RAlists,
        #        "left":self.left, "right":self.right, "length": self.length,
        #        "RAs_left":self.RAs_left, "RAs_right":self.RAs_right}
        return (self.chrom, self.peak, self.RAlists, self.left, self.right,
                self.length, self.RAs_left, self.RAs_right, self.peak_refseq,
                self.peak_refseq_ext)

    def __setstate__(self, state):
        (self.chrom, self.peak, self.RAlists, self.left, self.right,
         self.length, self.RAs_left, self.RAs_right, self.peak_refseq,
         self.peak_refseq_ext) = state

    @cython.ccall
    def sort(self):
        """Sort RAs according to lpos. Should be used after realignment.

        """
        if self.RAlists[0]:
            self.RAlists[0].sort(key=itemgetter("lpos"))
        if self.RAlists[1]:
            self.RAlists[1].sort(key=itemgetter("lpos"))
        self.sorted = True
        return

    @cython.ccall
    def remove_outliers(self, percent: cython.int = 5):
        """ Remove outliers with too many n_edits. The outliers with
        n_edits in top p% will be removed.

        Default: remove top 5% of reads that have too many differences
        with reference genome.
        """
        n_edits_list: list
        ralist: list
        read: ReadAlignment         # ReadAlignment object
        highest_n_edits: cython.int
        new_RAlist: list
        i: cython.int

        n_edits_list = []
        for ralist in self.RAlists:
            for read in ralist:
                n_edits_list.append(read["n_edits"])
        n_edits_list.sort()
        highest_n_edits = n_edits_list[int(len(n_edits_list) * (1 - percent * .01))]

        for i in (range(len(self.RAlists))):
            new_RAlist = []
            for read in self.RAlists[i]:
                if read["n_edits"] <= highest_n_edits:
                    new_RAlist.append(read)
            self.RAlists[i] = new_RAlist

        return

    @cython.ccall
    def n_edits_sum(self) -> cython.int:
        """
        """
        n_edits_list: list
        ralist: list
        read: ReadAlignment
        c: cython.int
        # highest_n_edits: cython.int

        n_edits_list = []

        for ralist in self.RAlists:
            for read in ralist:
                n_edits_list.append(read["n_edits"])

        n_edits_list.sort()
        # print (n_edits_list)
        c = Counter(n_edits_list)
        return c
        # print(c)

    @cython.cfunc
    def __get_peak_REFSEQ(self) -> tuple:
        """Get the reference sequence within the peak region.

        """
        peak_refseq: bytearray
        # i: cython.int
        # prev_r: cython.long                   #remember the previous filled right end
        start: cython.long
        end: cython.long
        # ind: cython.long
        # ind_r: cython.long
        # read: ReadAlignment
        # read_refseq_ext: bytearray
        # read_refseq: bytearray

        start = min(self.RAs_left, self.left)
        end = max(self.RAs_right, self.right)
        # print ("left",start,"right",end)
        peak_refseq_ext = bytearray(b'N' * (end - start))

        # for treatment.
        peak_refseq_ext = self.__fill_refseq(peak_refseq_ext,
                                             self.RAlists[0])
        # and control if available.
        if self.RAlists[1]:
            peak_refseq_ext = self.__fill_refseq(peak_refseq_ext,
                                                 self.RAlists[1])

        # trim
        peak_refseq = peak_refseq_ext[self.left - start: self.right - start]
        return (bytes(peak_refseq), bytes(peak_refseq_ext))

    @cython.cfunc
    def __fill_refseq(self,
                      seq: bytearray,
                      ralist: list) -> bytearray:
        """Fill refseq sequence of whole peak with refseq sequences of
        each read in ralist.

        """
        prev_r: cython.long         # previous right position of last
        # filled
        ind: cython.long
        ind_r: cython.long
        start: cython.long
        # end: cython.long
        read: ReadAlignment
        read_refseq: bytearray

        start = min(self.RAs_left, self.left)

        # print(len(ralist),"\n")
        prev_r = ralist[0]["lpos"]

        for i in range(len(ralist)):
            read = ralist[i]
            if read["lpos"] > prev_r:
                read = ralist[i - 1]
                read_refseq = read.get_REFSEQ()
                ind = read["lpos"] - start
                ind_r = ind + read["rpos"] - read["lpos"]
                seq[ind: ind_r] = read_refseq
                prev_r = read["rpos"]
        # last
        read = ralist[-1]
        read_refseq = read.get_REFSEQ()
        ind = read["lpos"] - start
        ind_r = ind + read["rpos"] - read["lpos"]
        seq[ind: ind_r] = read_refseq
        return seq

    @cython.ccall
    def get_PosReadsInfo_ref_pos(self,
                                 ref_pos: cython.long,
                                 ref_nt: bytes,
                                 Q: cython.int = 20):
        """Generate a PosReadsInfo for: object a given reference genome
        position.

        Return a PosReadsInfo object.

        """
        s: bytearray
        bq: bytearray
        strand: cython.int
        ra: ReadAlignment
        # bq_list_t: list = []
        # bq_list_c: list = []
        i: cython.int
        pos: cython.int
        tip: bool
        posreadsinfo_p: PosReadsInfo

        posreadsinfo_p = PosReadsInfo(ref_pos, ref_nt)

        # Treatment group
        for i in range(len(self.RAlists[0])):
            ra = self.RAlists[0][i]
            if ra["lpos"] <= ref_pos and ra["rpos"] > ref_pos:
                (s, bq, strand, tip, pos) = ra.get_variant_bq_by_ref_pos(ref_pos)
                posreadsinfo_p.add_T(i, bytes(s), bq[0], strand, tip, Q=Q)

        # Control group
        for i in range(len(self.RAlists[1])):
            ra = self.RAlists[1][i]
            if ra["lpos"] <= ref_pos and ra["rpos"] > ref_pos:
                (s, bq, strand, tip, pos) = ra.get_variant_bq_by_ref_pos(ref_pos)
                posreadsinfo_p.add_C(i, bytes(s), bq[0], strand, Q=Q)

        return posreadsinfo_p

    @cython.ccall
    def get_FASTQ(self) -> bytearray:
        """Get FASTQ file for all reads in RACollection.

        """
        ra: ReadAlignment
        fastq_text: bytearray

        fastq_text = bytearray(b"")

        for ra in self.RAlists[0]:
            fastq_text += ra.get_FASTQ()

        for ra in self.RAlists[1]:
            fastq_text += ra.get_FASTQ()

        return fastq_text

    @cython.cfunc
    def fermi_assemble(self,
                       fermiMinOverlap: cython.int,
                       opt_flag: cython.int = 0x80) -> list:
        """A wrapper function to call Fermi unitig building functions.
        """
        opt: cython.pointer(fml_opt_t)
        # c: cython.int
        n_seqs: cython.int
        n_utg: cython.pointer(cython.int)
        seqs: cython.pointer(bseq1_t)
        utg: cython.pointer(fml_utg_t)
        p: fml_utg_t

        # unitig_k: cython.int
        # merge_min_len: cython.int
        tmps: bytes
        tmpq: bytes
        # ec_k: cython.int = -1
        l: cython.long
        cseq: cython.pointer(cython.char)
        cqual: cython.pointer(cython.char)
        i: cython.int
        j: cython.int
        # tmpunitig: bytes
        # unitig: bytes                 # final unitig
        unitig_list: list             # contain of: list sequences in format: bytes
        # n: cython.pointer(cython.int)

        n_seqs = len(self.RAlists[0]) + len(self.RAlists[1])

        # prn_seqs: cython.int

        # prepare seq and qual, note, we only extract SEQ according to the +
        # strand of reference sequence.
        seqs = cython.cast(cython.pointer(bseq1_t),
                           malloc(n_seqs * cython.sizeof(bseq1_t)))  # we rely on fermi-lite to free this mem

        i = 0
        for ra in self.RAlists[0]:
            tmps = ra["SEQ"]
            tmpq = ra["QUAL"]
            l = len(tmps)
            # we rely on fermi-lite to free this mem
            cseq = cython.cast(cython.pointer(cython.char),
                               malloc((l+1)*cython.sizeof(cython.char)))
            # we rely on fermi-lite to free this mem
            cqual = cython.cast(cython.pointer(cython.char),
                                malloc((l+1)*cython.sizeof(cython.char)))
            for j in range(l):
                cseq[j] = tmps[j]
                cqual[j] = tmpq[j] + 33
            cseq[l] = b'\x00'
            cqual[l] = b'\x00'

            seqs[i].seq = cseq
            seqs[i].qual = cqual
            seqs[i].l_seq = len(tmps)
            i += 1

            # print "@",ra["readname"].decode()
            # prcseq: cython.int.decode()
            # print "+"
            # prcqual: cython.int.decode()

        for ra in self.RAlists[1]:
            tmps = ra["SEQ"]
            tmpq = ra["QUAL"]
            l = len(tmps)
            # we rely on fermi-lite to free this mem
            cseq = cython.cast(cython.pointer(cython.char),
                               malloc((l+1)*cython.sizeof(cython.char)))
            # we rely on fermi-lite to free this mem
            cqual = cython.cast(cython.pointer(cython.char),
                                malloc((l+1)*cython.sizeof(cython.char)))
            for j in range(l):
                cseq[j] = tmps[j]
                cqual[j] = tmpq[j] + 33
            cseq[l] = b'\x00'
            cqual[l] = b'\x00'

            seqs[i].seq = cseq
            seqs[i].qual = cqual
            seqs[i].l_seq = len(tmps)
            i += 1
            # print "@",ra["readname"].decode()
            # prcseq: cython.int.decode()
            # print "+"
            # prcqual: cython.int.decode()

        # if self.RAlists[1]:
        #     unitig_k=int(min(self.RAlists[0][0]["l"],self.RAlists[1][0]["l"])*fermiOverlapMinRatio)

        #     merge_min_len=int(min(self.RAlists[0][0]["l"],self.RAlists[1][0]["l"])*0.5)
        # else:
        #     unitig_k = int(self.RAlists[0][0]["l"]*fermiOverlapMinRatio)

        #     merge_min_len=int(self.RAlists[0][0]["l"]*0.5)
        #fermiMinOverlap = int(self.RAlists[0][0]["l"]*fermiOverlapMinRatio)

        # minimum overlap to merge, default 0
        # merge_min_len= max(25, int(self.RAlists[0][0]["l"]*0.5))
        # merge_min_len= int(self.RAlists[0][0]["l"]*0.5)

        # opt = cython.cast(cython.pointer(fml_opt_t),
        #                   PyMem_Malloc(cython.sizeof(fml_opt_t)))
        # n_utg = cython.cast(cython.pointer(cython.int),
        #                     PyMem_Malloc(cython.sizeof(int)))

        opt = cython.cast(cython.pointer(fml_opt_t),
                          malloc(cython.sizeof(fml_opt_t)))
        n_utg = cython.cast(cython.pointer(cython.int),
                            malloc(cython.sizeof(int)))

        fml_opt_init(opt)
        # k-mer length for error correction (0 for auto; -1 to disable)
        # opt.ec_k = 0

        # min overlap length during initial assembly
        opt.min_asm_ovlp = fermiMinOverlap

        # minimum length to merge, during assembly, don't explicitly merge an overlap if shorter than this value
        # opt.min_merge_len = merge_min_len

        # there are more 'options' for mag clean:
        # flag, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bdiff, max_bvtx, min_merge_len, trim_len, trim_depth, min_dratio1, max_bcov, max_bfrac
        # min_elen (300) will be adjusted
        # min_ensr (4), min_insr (3) will be computed
        # min_merge_len (0) will be updated using opt.min_merge_len

        # We can adjust: flag (0x40|0x80), min_ovlp (0), min_dratio1 (0.7), max_bdiff (50), max_bdist (512), max_bvtx (64), trim_len (0), trim_depth (6), max_bcov (10.), max_bfrac (0.15)

        # 0x20: MAG_F_AGGRESSIVE pop variant bubbles
        # 0x40: MAG_F_POPOPEN aggressive tip trimming
        # 0x80: MAG_F_NO_SIMPL skip bubble simplification
        opt.mag_opt.flag = opt_flag

        # mag_opt.min_ovlp
        #opt.mag_opt.min_ovlp = fermiMinOverlap

        # drop an overlap if its length is below maxOvlpLen*FLOAT
        #opt.mag_opt.min_dratio1 = 0.5

        # retain a bubble if one side is longer than the other side by >INT-bp
        #opt.mag_opt.max_bdiff = 10#merge_min_len

        # trim_len:
        # trim_depth: Parameter used to trim the open end/tip. If trim_len == 0, do nothing

        # max_bdist:
        # max_bvtx: Parameter used to simply bubble while 0x80 flag is set.
        #opt.mag_opt.max_bdist = 1024
        #opt.mag_opt.max_bvtx = 128

        # max_bcov:
        # max_bfrac: Parameter used when aggressive bubble removal is not used. Bubble will be removed if its average coverage lower than max_bcov and fraction (cov1/(cov1+cov2)) is lower than max_bfrac
        #opt.mag_opt.max_bcov = 10.
        #opt.mag_opt.max_bfrac = 0.01

        utg = fml_assemble(opt, n_seqs, seqs, n_utg)
        # get results
        unitig_list = []
        for i in range(n_utg[0]):
            p = utg[i]
            if (p.len < 0):
                continue
            # unitig = b''
            # for j in range(p.len):
            #    unitig += [b'A',b'C',b'G',b'T',b'N'][int(p.seq[j]) - 1]
            # unitig_list.append(unitig)
            unitig_list.append(p.seq)

        fml_utg_destroy(n_utg[0], utg)

        # PyMem_Free(opt)
        # PyMem_Free(n_utg)
        free(opt)
        free(n_utg)

        return unitig_list

    @cython.cfunc
    def align_unitig_to_REFSEQ(self, unitig_list: list) -> tuple:
        """Note: we use smith waterman, but we don't use linear gap
        penalty at this time.

        Also, if unitig is mapped to - strand, we will revcomp the
        unitig. So the unitig_will: list be changed in this case.
        """
        unitig: bytes
        problem: seq_pair_t
        results: cython.pointer(align_t)
        # tmp: cython.pointer(cython.char)
        target: bytes
        reference: bytes
        target_aln_f: bytes
        target_aln_r: bytes
        reference_aln_f: bytes
        reference_aln_r: bytes
        markup_aln_f: bytes
        markup_aln_r: bytes
        score_f: cython.double
        score_r: cython.double
        target_alns: list = []
        reference_alns: list = []
        markup_alns: list = []
        aln_scores: list = []
        i: cython.int

        reference = copy(self.peak_refseq_ext+b'\x00')

        for i in range(len(unitig_list)):
            unitig = unitig_list[i]
            target = copy(unitig + b'\x00')
            # we use swalign.c for local alignment (without affine gap
            # penalty). Will revise later.
            problem.a = target
            problem.alen = len(unitig)
            problem.b = reference
            problem.blen = len(self.peak_refseq_ext)
            results = smith_waterman(cython.address(problem))
            target_aln_f = results.seqs.a
            reference_aln_f = results.seqs.b
            markup_aln_f = results.markup
            score_f = results.score
            free(results.seqs.a)
            free(results.seqs.b)
            free(results.markup)
            free(results)
            # end of local alignment

            # try reverse complement
            target = copy(unitig[::-1] + b'\x00')
            target = target.translate(__DNACOMPLEMENT__)
            problem.a = target
            problem.alen = len(unitig)
            problem.b = reference
            problem.blen = len(self.peak_refseq_ext)
            results = smith_waterman(cython.address(problem))
            target_aln_r = results.seqs.a
            reference_aln_r = results.seqs.b
            markup_aln_r = results.markup
            score_r = results.score
            free(results.seqs.a)
            free(results.seqs.b)
            free(results.markup)
            free(results)
            # end of local alignment

            if score_f > score_r:
                target_alns.append(target_aln_f)
                reference_alns.append(reference_aln_f)
                markup_alns.append(markup_aln_f)
                aln_scores.append(score_f)
            else:
                target_alns.append(target_aln_r)
                reference_alns.append(reference_aln_r)
                markup_alns.append(markup_aln_r)
                aln_scores.append(score_r)
                # we will revcomp unitig
                unitig = unitig[::-1]
                unitig_list[i] = unitig.translate(__DNACOMPLEMENT__)

        return (target_alns, reference_alns, aln_scores, markup_alns)

    @cython.cfunc
    def verify_alns(self, unitig_list, unitig_alns, reference_alns, aln_scores, markup_alns, min_score_100: cython.float = 150):
        """Remove aln/unitig if it contains too many edits in a small region

        default min score is 150, which means under 2/-3/-5/-2 scoring schema, there are 10 mismatches within 100bps region.
        """
        i: cython.int

        for i in range(len(unitig_list)-1, -1, -1):
            # pri: cython.int, aln_scores[i]
            # prunitig_alns: cython.int[i]
            # prmarkup_alns: cython.int[i]
            # prreference_alns: cython.int[i]
            if aln_scores[i] * 100 / len(markup_alns[i]) < min_score_100:
                unitig_list.pop(i)
                unitig_alns.pop(i)
                reference_alns.pop(i)
                aln_scores.pop(i)
                markup_alns.pop(i)
        return

    @cython.cfunc
    def filter_unitig_with_bad_aln(self, unitig_list: list,
                                   target_alns: list,
                                   reference_alns: list,
                                   gratio: float = 0.25) -> tuple:
        """Remove unitigs that has too much gaps (both on target and
        reference) during alignments.

        """
        pass

    @cython.cfunc
    def remap_RAs_w_unitigs(self, unitig_list: list) -> list:
        """Remap RAs to unitigs, requiring perfect match.

        Return RAlists_T, RAlists_C, unmapped_racollection.
        """
        RAlists_T: list = []  # lists of of: list RAs of ChIP mapped to each unitig
        RAlists_C: list = []
        unmapped_RAlist_T: list = []  # of: list RAs of ChIP unmappable to unitigs
        unmapped_RAlist_C: list = []
        # RACollection unmapped_ra_collection
        flag: cython.int = 0
        i: cython.int
        tmp_ra: ReadAlignment
        tmp_ra_seq: bytes
        unitig: bytes

        for i in range(len(unitig_list)):
            RAlists_T.append([])         # for each unitig, there is another of: list RAs
            RAlists_C.append([])

        # assign RAs to unitigs

        for tmp_ra in self.RAlists[0]:
            flag = 0
            tmp_ra_seq = tmp_ra["SEQ"]
            for i in range(len(unitig_list)):
                unitig = unitig_list[i]
                if tmp_ra_seq in unitig:
                    flag = 1
                    RAlists_T[i].append(tmp_ra)
                    break
            if flag == 0:
                unmapped_RAlist_T.append(tmp_ra)
                # print "unmapped:", tmp_ra["SEQ"]

        for tmp_ra in self.RAlists[1]:
            flag = 0
            tmp_ra_seq = tmp_ra["SEQ"]
            for i in range(len(unitig_list)):
                unitig = unitig_list[i]
                if tmp_ra_seq in unitig:
                    flag = 1
                    RAlists_C[i].append(tmp_ra)
                    break
            if flag == 0:
                unmapped_RAlist_C.append(tmp_ra)
                # print "unmapped:", tmp_ra["SEQ"]

        # if unmapped_RAlist_T:
        # unmapped_ra_collection = RACollection(self.chrom, self.peak, unmapped_RAlist_T, unmapped_RAlist_C)
        return [RAlists_T, RAlists_C, unmapped_RAlist_T, unmapped_RAlist_C]

    @cython.cfunc
    def add_to_unitig_list(self, unitig_list, unitigs_2nd) -> list:
        """
        """
        i: cython.int
        j: cython.int
        flag: cython.int
        u0: bytes
        u1: bytes
        new_unitig_list: list

        new_unitig_list = []

        for i in range(len(unitigs_2nd)):
            # initial value: can't be found in unitig_list
            flag = 0
            u0 = unitigs_2nd[i]
            for j in range(len(unitig_list)):
                u1 = unitig_list[j]
                if u1.find(u0) != -1:
                    flag = 1
                    break
                u1 = u1[::-1].translate(__DNACOMPLEMENT__)
                if u1.find(u0) != -1:
                    flag = 1
                    break
            if not flag:
                new_unitig_list.append(u0)
        new_unitig_list.extend(unitig_list)
        return new_unitig_list

    @cython.ccall
    def build_unitig_collection(self, fermiMinOverlap):
        """unitig_and: list tuple_alns are in the same order!

        return UnitigCollection object.

        """
        start: cython.long
        end: cython.long
        unitigs_2nd: list
        # u: bytes
        # target_alns: list
        reference_alns: list
        aln_scores: list
        markup_alns: list
        # target_alns_2nd: list
        # reference_alns_2nd: list
        # aln_scores_2nd: list
        RAlists_T: list = []  # lists of of: list RAs of ChIP mapped to each unitig
        RAlists_C: list = []
        unmapped_RAlist_T: list = []
        unmapped_RAlist_C: list = []
        # tmp_unitig_seq: bytes
        tmp_reference_seq: bytes
        tmp_unitig_aln: bytes
        tmp_reference_aln: bytes

        i: cython.int
        j: cython.int
        left_padding_ref: cython.long
        right_padding_ref: cython.long
        left_padding_unitig: cython.long
        right_padding_unitig: cython.long
        ura_list: list = []
        unmapped_ra_collection: RACollection
        # flag: cython.int = 0
        # n_unmapped: cython.int
        n_unitigs_0: cython.int
        n_unitigs_1: cython.int

        # first round of assembly
        # print (" First round to assemble unitigs")
        unitig_list = self.fermi_assemble(fermiMinOverlap, opt_flag=0x80)
        if len(unitig_list) == 0:
            return 0

        n_unitigs_0 = -1
        n_unitigs_1 = len(unitig_list)
        # print " # of Unitigs:", n_unitigs_1
        # print " Map reads to unitigs"
        (unitig_alns, reference_alns, aln_scores, markup_alns) = self.align_unitig_to_REFSEQ(unitig_list)

        self.verify_alns(unitig_list,
                         unitig_alns,
                         reference_alns,
                         aln_scores,
                         markup_alns)
        if len(unitig_list) == 0:
            # if stop here, it raises a flag that the region may
            # contain too many mismapped reads, we return -1
            return -1
        # print (" # of Unitigs:", n_unitigs_1)

        # assign RAs to unitigs
        [RAlists_T, RAlists_C, unmapped_RAlist_T, unmapped_RAlist_C] = self.remap_RAs_w_unitigs(unitig_list)
        # prunmapped_ra_collection: cython.int.get_FASTQ().decode()

        # n_unmapped = len(unmapped_RAlist_T) + len(unmapped_RAlist_C)

        while len(unmapped_RAlist_T) > 0 and n_unitigs_1 != n_unitigs_0:
            # if there are unmapped reads AND we can get more unitigs
            # from last round of assembly, do assembly again

            # print (" # of RAs not mapped, will be assembled again:", n_unmapped)
            n_unitigs_0 = n_unitigs_1
            # another round of assembly
            unmapped_ra_collection = RACollection(self.chrom,
                                                  self.peak,
                                                  unmapped_RAlist_T,
                                                  unmapped_RAlist_C)
            unitigs_2nd = unmapped_ra_collection.fermi_assemble(fermiMinOverlap,
                                                                opt_flag=0x80)

            if unitigs_2nd:
                unitig_list = self.add_to_unitig_list(unitig_list, unitigs_2nd)
                n_unitigs_1 = len(unitig_list)
                # print " # of Unitigs:", n_unitigs_1
                # print " Map reads to unitigs"
                (unitig_alns, reference_alns, aln_scores, markup_alns) = self.align_unitig_to_REFSEQ(unitig_list)
                self.verify_alns(unitig_list,
                                 unitig_alns,
                                 reference_alns,
                                 aln_scores,
                                 markup_alns)
                [RAlists_T, RAlists_C, unmapped_RAlist_T, unmapped_RAlist_C] = self.remap_RAs_w_unitigs(unitig_list)
                # n_unmapped = len(unmapped_RAlist_T) + len(unmapped_RAlist_C)
            # else:
            #    for r in unmapped_RAlist_T:
            #        prr: cython.int.get_FASTQ().decode().lstrip()

            # print (" # of RAs not mapped, will be assembled again with 1/2 of fermiMinOverlap:", n_unmapped)
            # another round of assembly
            unmapped_ra_collection = RACollection(self.chrom,
                                                  self.peak,
                                                  unmapped_RAlist_T,
                                                  unmapped_RAlist_C)
            unitigs_2nd = unmapped_ra_collection.fermi_assemble(fermiMinOverlap/2,
                                                                opt_flag=0x80)

            if unitigs_2nd:
                unitig_list = self.add_to_unitig_list(unitig_list, unitigs_2nd)
                n_unitigs_1 = len(unitig_list)
                # print " # of Unitigs:", n_unitigs_1
                # print " Map reads to unitigs"
                (unitig_alns, reference_alns, aln_scores, markup_alns) = self.align_unitig_to_REFSEQ(unitig_list)
                self.verify_alns(unitig_list,
                                 unitig_alns,
                                 reference_alns,
                                 aln_scores,
                                 markup_alns)
                [RAlists_T, RAlists_C, unmapped_RAlist_T, unmapped_RAlist_C] = self.remap_RAs_w_unitigs(unitig_list)
                # n_unmapped = len(unmapped_RAlist_T) + len(unmapped_RAlist_C)
            # else:
            #    for r in unmapped_RAlist_T:
            #        prr: cython.int.get_FASTQ().decode().lstrip()
            if len(unitig_list) == 0:
                raise Exception("Shouldn't reach here")
            # print (" # of Unitigs:", n_unitigs_1)

        if len(unitig_list) == 0:
            return None
        # print (" Final round: # of Unitigs:", len(unitig_list))
        # print (" Final round: # of RAs not mapped:", n_unmapped)

        start = min(self.left, self.RAs_left)
        end = max(self.right, self.RAs_right)

        # create UnitigCollection
        for i in range(len(unitig_list)):
            #b'---------------------------AAATAATTTTATGTCCTTCAGTACAAAAAGCAGTTTCAACTAAAACCCAGTAACAAGCTAGCAATTCCTTTTAAATGGTGCTACTTCAAGCTGCAGCCAGGTAGCTTTTTATTACAAAAAATCCCACAGGCAGCCACTAGGTGGCAGTAACAGGCTTTTGCCAGCGGCTCCAGTCAGCATGGCTTGACTGTGTGCTGCAGAAACTTCTTAAATCGTCTGTGTTTGGGACTCGTGGGGCCCCACAGGGCTTTACAAGGGCTTTTTAATTTCCAAAAACATAAAACAAAAAAA--------------'
            #b'GATATAAATAGGATGTTATGAGTTTTCAAATAATTTTATGTCCTTCAGTACAAAAAGCAGTTTCAACTAAAACCCAGTAACAAGCTAGCAATTCCTTTTAAATGGTGCTACTTCAAGCTGCAGCCAGGTAGCTTTTTATTACAAAAA-TCCCACAGGCAGCCACTAGGTGGCAGTAACAGGCTTTTGCCAGCGGCTCCAGTCAGCATGGCTTGACTGTGTGCTGCAGAAACTTCTTAAATCGTCTGTGTTTGGGACTCGTGGGGCCCCACAGGGCTTTACAAGGGCTTTTTAATTTCCAAAAACATAAAACAAAAAAAAATACAAATGTATT'
            tmp_unitig_aln = unitig_alns[i]
            tmp_reference_aln = reference_alns[i]
            # tmp_unitig_seq = tmp_unitig_aln.replace(b'-',b'')
            tmp_reference_seq = tmp_reference_aln.replace(b'-', b'')

            # prtmp_unitig_aln: cython.int
            # prtmp_reference_aln: cython.int
            # prtmp_unitig_seq: cython.int
            # prtmp_reference_aln: cython.int

            # find the position on self.peak_refseq_ext
            left_padding_ref = self.peak_refseq_ext.find(tmp_reference_seq)  # this number of nts should be skipped on refseq_ext from left
            right_padding_ref = len(self.peak_refseq_ext) - left_padding_ref - len(tmp_reference_seq)  # this number of nts should be skipped on refseq_ext from right

            # now, decide the lpos and rpos on reference of this unitig
            # first, trim left padding '-'
            left_padding_unitig = len(tmp_unitig_aln) - len(tmp_unitig_aln.lstrip(b'-'))
            right_padding_unitig = len(tmp_unitig_aln) - len(tmp_unitig_aln.rstrip(b'-'))

            tmp_lpos = start + left_padding_ref
            tmp_rpos = end - right_padding_ref

            for j in range(left_padding_unitig):
                if tmp_reference_aln[j] != b'-':
                    tmp_lpos += 1
            for j in range(1, right_padding_unitig + 1):
                if tmp_reference_aln[-j] != b'-':
                    tmp_rpos -= 1

            tmp_unitig_aln = tmp_unitig_aln[left_padding_unitig:(len(tmp_unitig_aln)-right_padding_unitig)]
            tmp_reference_aln = tmp_reference_aln[left_padding_unitig:(len(tmp_reference_aln)-right_padding_unitig)]

            ura_list.append(UnitigRAs(self.chrom, tmp_lpos, tmp_rpos, tmp_unitig_aln, tmp_reference_aln, [RAlists_T[i], RAlists_C[i]]))

        return UnitigCollection(self.chrom, self.peak, ura_list)
