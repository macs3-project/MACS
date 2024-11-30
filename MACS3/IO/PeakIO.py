# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-15 11:48:33 Tao Liu>

"""Module for PeakIO IO classes.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
from itertools import groupby
from operator import itemgetter
import random
import sys

# ------------------------------------
# MACS3 modules
# ------------------------------------

# from MACS3.Utilities.Constants import *

# ------------------------------------
# Other modules
# ------------------------------------
import cython
from cython.cimports.cpython import bool

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------


@cython.cfunc
def subpeak_letters(i: cython.int) -> str:
    if i < 26:
        return chr(97+i)
    else:
        return subpeak_letters(i // 26) + chr(97 + (i % 26))

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class PeakContent:
    chrom: bytes
    start: cython.int
    end: cython.int
    length: cython.int
    summit: cython.int
    score: cython.float
    pileup: cython.float
    pscore: cython.float
    fc: cython.float
    qscore: cython.float
    name: bytes

    def __init__(self,
                 chrom: bytes,
                 start: cython.int,
                 end: cython.int,
                 summit: cython.int,
                 peak_score: cython.float,
                 pileup: cython.float,
                 pscore: cython.float,
                 fold_change: cython.float,
                 qscore: cython.float,
                 name: bytes = b""):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = end - start
        self.summit = summit
        self.score = peak_score
        self.pileup = pileup
        self.pscore = pscore
        self.fc = fold_change
        self.qscore = qscore
        self.name = name

    def __getitem__(self, a: str):
        if a == "chrom":
            return self.chrom
        elif a == "start":
            return self.start
        elif a == "end":
            return self.end
        elif a == "length":
            return self.length
        elif a == "summit":
            return self.summit
        elif a == "score":
            return self.score
        elif a == "pileup":
            return self.pileup
        elif a == "pscore":
            return self.pscore
        elif a == "fc":
            return self.fc
        elif a == "qscore":
            return self.qscore
        elif a == "name":
            return self.name

    def __setitem__(self, a: str, v):
        if a == "chrom":
            self.chrom = v
        elif a == "start":
            self.start = v
        elif a == "end":
            self.end = v
        elif a == "length":
            self.length = v
        elif a == "summit":
            self.summit = v
        elif a == "score":
            self.score = v
        elif a == "pileup":
            self.pileup = v
        elif a == "pscore":
            self.pscore = v
        elif a == "fc":
            self.fc = v
        elif a == "qscore":
            self.qscore = v
        elif a == "name":
            self.name = v

    def __str__(self):
        return "chrom:%s;start:%d;end:%d;score:%f" % (self.chrom,
                                                      self.start,
                                                      self.end,
                                                      self.score)

    def __getstate__(self):
        return (self.chrom,
                self.start,
                self.end,
                self.length,
                self.summit,
                self.score,
                self.pileup,
                self.pscore,
                self.fc,
                self.qscore,
                self.name)

    def __setstate__(self, state):
        (self.chrom, self.start, self.end, self.length, self.summit,
         self.score, self.pileup, self.pscore, self.fc,
         self.qscore, self.name) = state


@cython.cclass
class PeakIO:
    """IO for peak information.

    """
    # dictionary storing peak contents
    peaks = cython.declare(dict, visibility="public")
    # whether peaks have been sorted by coordinations
    CO_sorted = cython.declare(bool, visibility="public")
    # total number of peaks
    total = cython.declare(cython.long, visibility="public")

    def __init__(self):
        self.peaks = {}
        self.CO_sorted = False
        self.total = 0

    @cython.ccall
    def add(self,
            chromosome: bytes,
            start: cython.int,  # leftmost position
            end: cython.int,    # rightmost position
            summit: cython.int = 0,  # summit position
            peak_score: cython.float = 0,  # score
            pileup: cython.float = 0,      # pileup value
            pscore: cython.float = 0,      # -log10 pvalue
            fold_change: cython.float = 0,  # fold change
            qscore: cython.float = 0,      # -log10 qvalue
            name: bytes = b""):            # peak name
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome] = []
        self.peaks[chromosome].append(PeakContent(chromosome,
                                                  start,
                                                  end,
                                                  summit,
                                                  peak_score,
                                                  pileup,
                                                  pscore,
                                                  fold_change,
                                                  qscore,
                                                  name))
        self.total += 1
        self.CO_sorted = False

    @cython.ccall
    def add_PeakContent(self,
                        chromosome: bytes,
                        peakcontent: PeakContent):
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome] = []
        self.peaks[chromosome].append(peakcontent)
        self.total += 1
        self.CO_sorted = False

    @cython.ccall
    def get_data_from_chrom(self, chrom: bytes) -> list:
        if not self.peaks.has_key(chrom):
            self.peaks[chrom] = []
        return self.peaks[chrom]

    def get_chr_names(self) -> set:
        return set(self.peaks.keys())

    def sort(self):
        chrs: list
        chrom: bytes

        # sort by position
        if self.CO_sorted:
            # if already sorted, quit
            return
        chrs = sorted(list(self.peaks.keys()))
        for chrom in sorted(chrs):
            self.peaks[chrom].sort(key=lambda x: x['start'])
        self.CO_sorted = True
        return

    @cython.ccall
    def randomly_pick(self, n: cython.int, seed: cython.int = 12345):
        """Shuffle the peaks and get n peaks out of it. Return a new
        PeakIO object.
        """
        all_pc: list
        chrs: list
        chrom: bytes
        ret_peakio: PeakIO
        p: PeakContent

        assert n > 0
        chrs = sorted(list(self.peaks.keys()))
        all_pc = []
        for chrom in sorted(chrs):
            all_pc.extend(self.peaks[chrom])
        random.seed(seed)
        random.shuffle(all_pc)
        all_pc = all_pc[:n]
        ret_peakio = PeakIO()
        for p in all_pc:
            ret_peakio.add_PeakContent(p["chrom"], p)
        return ret_peakio

    @cython.ccall
    def filter_pscore(self, pscore_cut: cython.double):
        chrom: bytes
        new_peaks: dict
        chrs: list

        new_peaks = {}
        chrs = sorted(list(self.peaks.keys()))
        self.total = 0
        for chrom in sorted(chrs):
            new_peaks[chrom] = [p for p in self.peaks[chrom] if p['pscore'] >= pscore_cut]
            self.total += len(new_peaks[chrom])
        self.peaks = new_peaks
        self.CO_sorted = True
        self.sort()

    @cython.ccall
    def filter_qscore(self, qscore_cut: cython.double):
        chrom: bytes
        new_peaks: dict
        chrs: list

        new_peaks = {}
        chrs = sorted(list(self.peaks.keys()))
        self.total = 0
        for chrom in sorted(chrs):
            new_peaks[chrom] = [p for p in self.peaks[chrom] if p['qscore'] >= qscore_cut]
            self.total += len(new_peaks[chrom])
        self.peaks = new_peaks
        self.CO_sorted = True
        self.sort()

    @cython.ccall
    def filter_fc(self, fc_low: cython.float, fc_up: cython.float = 0):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)

        """
        chrom: bytes
        new_peaks: dict
        chrs: list

        new_peaks = {}
        chrs = list(self.peaks.keys())
        self.total = 0
        if fc_up > 0 and fc_up > fc_low:
            for chrom in sorted(chrs):
                new_peaks[chrom] = [p for p in self.peaks[chrom] if p['fc'] >= fc_low and p['fc'] < fc_up]
                self.total += len(new_peaks[chrom])
        else:
            for chrom in sorted(chrs):
                new_peaks[chrom] = [p for p in self.peaks[chrom] if p['fc'] >= fc_low]
                self.total += len(new_peaks[chrom])
        self.peaks = new_peaks
        self.CO_sorted = True
        self.sort()

    def filter_score(self, lower_score: cython.float, upper_score: cython.float = 0):
        """Filter peaks in a given score range.

        """
        chrom: bytes
        new_peaks: dict
        chrs: list

        new_peaks = {}
        chrs = list(self.peaks.keys())
        self.total = 0
        if upper_score > 0 and upper_score > lower_score:
            for chrom in sorted(chrs):
                new_peaks[chrom] = [p for p in self.peaks[chrom] if p['score'] >= lower_score and p['score'] < upper_score]
                self.total += len(new_peaks[chrom])
        else:
            for chrom in sorted(chrs):
                new_peaks[chrom] = [p for p in self.peaks[chrom] if p['score'] >= lower_score]
                self.total += len(new_peaks[chrom])
        self.peaks = new_peaks
        self.CO_sorted = True
        self.sort()

    def __str__(self):
        """convert to text -- for debug
        """
        chrs: list
        n_peak: cython.int
        ret: str
        chrom: bytes
        peaks: list

        ret = ""
        chrs = list(self.peaks.keys())
        n_peak = 0
        for chrom in sorted(chrs):
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        ret += "chrom:%s\tstart:%d\tend:%d\tname:peak_%d%s\tscore:%.6g\tsummit:%d\n" % (chrom.decode(), peak['start'], peak['end'], n_peak, subpeak_letters(i), peak["score"], peak["summit"])
                else:
                    peak = peaks[0]
                    ret += "chrom:%s\tstart:%d\tend:%d\tname:peak_%d\tscore:%.6g\tsummit:%d\n" % (chrom.decode(), peak['start'], peak['end'], n_peak, peak["score"], peak["summit"])
        return ret

    @cython.cfunc
    def _to_bed(self,
                name_prefix: bytes = b"%s_peak_",
                name: bytes = b"MACS",
                description: bytes = b"%s",
                score_column: str = "score",
                trackline: bool = False,
                print_func=sys.stdout.write):
        """
        generalization of tobed and write_to_bed
        """
        chrs: list
        n_peak: cython.int
        peakprefix: bytes
        desc: bytes

        chrs = list(self.peaks.keys())
        n_peak = 0
        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix
        try:
            desc = description % name
        except Exception:
            desc = description

        if trackline:
            try:
                print_func('track name="%s (peaks)" description="%s" visibility=1\n' % (name.replace(b"\"", b"\\\"").decode(),
                                                                                        desc.replace(b"\"", b"\\\"").decode()))
            except Exception:
                print_func('track name=MACS description=Unknown\n')
        for chrom in sorted(chrs):
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        print_func("%s\t%d\t%d\t%s%d%s\t%.6g\n" % (chrom.decode(), peak['start'], peak['end'], peakprefix.decode(), n_peak, subpeak_letters(i), peak[score_column]))
                else:
                    peak = peaks[0]
                    print_func("%s\t%d\t%d\t%s%d\t%.6g\n" % (chrom.decode(), peak['start'], peak['end'], peakprefix.decode(), n_peak, peak[score_column]))

    @cython.cfunc
    def _to_summits_bed(self,
                        name_prefix: bytes = b"%s_peak_",
                        name: bytes = b"MACS",
                        description: bytes = b"%s",
                        score_column: str = "score",
                        trackline: bool = False,
                        print_func=sys.stdout.write):
        """
        generalization of to_summits_bed and write_to_summit_bed
        """
        chrs: list
        n_peak: cython.int
        peakprefix: bytes
        desc: bytes

        chrs = list(self.peaks.keys())
        n_peak = 0
        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix
        try:
            desc = description % name
        except Exception:
            desc = description
        if trackline:
            try:
                print_func('track name="%s (summits)" description="%s" visibility=1\n' % (name.replace(b"\"", b"\\\"").decode(),
                                                                                          desc.replace(b"\"", b"\\\"").decode()))
            except Exception:
                print_func('track name=MACS description=Unknown')
        for chrom in sorted(chrs):
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        summit_p = peak['summit']
                        print_func("%s\t%d\t%d\t%s%d%s\t%.6g\n" % (chrom.decode(), summit_p, summit_p+1, peakprefix.decode(), n_peak, subpeak_letters(i), peak[score_column]))
                else:
                    peak = peaks[0]
                    summit_p = peak['summit']
                    print_func("%s\t%d\t%d\t%s%d\t%.6g\n" % (chrom.decode(), summit_p, summit_p+1, peakprefix.decode(), n_peak, peak[score_column]))

    def tobed(self):
        """Print out (stdout) peaks in BED5 format.

        Five columns are chromosome, peak start, peak end, peak name, and peak height.

        start:start
        end:end,
        length:end-start,
        summit:summit,
        score:peak_score,
        pileup:pileup,
        pscore:pvalue,
        fc:fold_change,
        qscore:qvalue
        """
        return self._to_bed(name_prefix=b"%s_peak_", score_column="score", name=self.name, description=b"")

    def to_summits_bed(self):
        """Print out (stdout) peak summits in BED5 format.

        Five columns are chromosome, summit start, summit end, peak name, and peak height.

        """
        return self._to_summits_bed(name_prefix=b"%s_peak_", score_column="score", name=self.name, description=b"")

    # these methods are very fast, specifying types is unnecessary
    def write_to_bed(self, fhd,
                     name_prefix: bytes = b"%s_peak_",
                     name: bytes = b"MACS",
                     description: bytes = b"%s",
                     score_column: str = "score",
                     trackline: bool = True):
        """Write peaks in BED5 format in a file handler. Score (5th
        column) is decided by score_column setting. Check the
        following list. Name column ( 4th column) is made by putting
        name_prefix together with an ascending number.

        Five columns are chromosome, peak start, peak end, peak name,
        and peak score.

        items in peak hash object:

        start:start
        end:end,
        length:end-start,
        summit:summit,
        score:peak_score,
        pileup:pileup,
        pscore:pvalue,
        fc:fold_change,
        qscore:qvalue
        """
        # print(description)
        return self._to_bed(name_prefix=name_prefix,
                            name=name,
                            description=description,
                            score_column=score_column,
                            print_func=fhd.write,
                            trackline=trackline)

    def write_to_summit_bed(self, fhd,
                            name_prefix: bytes = b"%s_peak_",
                            name: bytes = b"MACS",
                            description: bytes = b"%s",
                            score_column: str = "score",
                            trackline: bool = False):
        """Write peak summits in BED5 format in a file handler. Score
        (5th column) is decided by score_column setting. Check the
        following list. Name column ( 4th column) is made by putting
        name_prefix together with an ascending number.

        Five columns are chromosome, summit start, summit end, peak name, and peak score.

        items in peak object:

        start:start
        end:end,
        length:end-start,
        summit:summit,
        score:peak_score,
        pileup:pileup,
        pscore:pvalue,
        fc:fold_change,
        qscore:qvalue
        """
        return self._to_summits_bed(name_prefix=name_prefix, name=name,
                                    description=description, score_column=score_column,
                                    print_func=fhd.write, trackline=trackline)

    def write_to_narrowPeak(self, fhd,
                            name_prefix: bytes = b"%s_peak_",
                            name: bytes = b"MACS",
                            score_column: str = "score",
                            trackline: bool = False):
        """Print out peaks in narrowPeak format.

        This format is designed for ENCODE project, and basically a
        BED6+4 format.

        +-----------+------+----------------------------------------+
        |field      |type  |description                             |
        +-----------+------+----------------------------------------+
        |chrom      |string|Name of the chromosome                  |
        +-----------+------+----------------------------------------+
        |chromStart |int   |The starting position of the feature in |
        |           |      |the chromosome. The first base in a     |
        |           |      |chromosome is numbered 0.               |
        +-----------+------+----------------------------------------+
        |chromEnd   |int   |The ending position of the feature in   |
        |           |      |the chromosome or scaffold. The chromEnd|
        |           |      |base is not included in the display of  |
        |           |      |the feature.  For example, the first 100|
        |           |      |bases of a chromosome are defined as    |
        |           |      |chromStart=0, chromEnd=100, and span the|
        |           |      |bases numbered 0-99.                    |
        +-----------+------+----------------------------------------+
        |name       |string|Name given to a region (preferably      |
        |           |      |unique). Use '.' if no name is assigned.|
        +-----------+------+----------------------------------------+
        |score      |int   |Indicates how dark the peak will be     |
        |(-logpvalue|      |displayed in the browser (1-1000). If   |
        |in MACS3 * |      |'0', the DCC will assign this based on  |
        |10)        |      |signal value. Ideally average           |
        |           |      |signalValue per base spread between     |
        |           |      |100-1000.                               |
        +-----------+------+----------------------------------------+
        |strand     |char  |+/- to denote strand or orientation     |
        |(always .) |      |(whenever applicable). Use '.' if no    |
        |           |      |orientation is assigned.                |
        +-----------+------+----------------------------------------+
        |signalValue|float |Measurement of overall (usually,        |
        |(fc)       |      |average) enrichment for the region.     |
        +-----------+------+----------------------------------------+
        |pValue     |float |Measurement of statistical signficance  |
        |           |      |(-log10). Use -1 if no pValue is        |
        |           |      |assigned.                               |
        +-----------+------+----------------------------------------+
        |qValue     |float |Measurement of statistical significance |
        |           |      |using false discovery rate. Use -1 if no|
        |           |      |qValue is assigned.                     |
        +-----------+------+----------------------------------------+
        |peak       |int   |Point-source called for this peak;      |
        |           |      |0-based offset from chromStart. Use -1  |
        |           |      |if no point-source called.              |
        +-----------+------+----------------------------------------+

        """
        n_peak: cython.int
        chrom: bytes
        s: cython.long
        peakname: str

        chrs = list(self.peaks.keys())
        n_peak = 0
        write = fhd.write
        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix
        if trackline:
            write("track type=narrowPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name.decode(), name.decode()))
        for chrom in sorted(chrs):
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1:  # from call-summits
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix.decode(), n_peak, subpeak_letters(i))
                        if peak['summit'] == -1:
                            s = -1
                        else:
                            s = peak['summit'] - peak['start']
                        fhd.write("%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n" %
                                  (chrom.decode(),
                                   peak['start'],
                                   peak['end'],
                                   peakname,
                                   int(10*peak[score_column]),
                                   peak['fc'],
                                   peak['pscore'],
                                   peak['qscore'],
                                   s))
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix.decode(), n_peak)
                    if peak['summit'] == -1:
                        s = -1
                    else:
                        s = peak['summit'] - peak['start']
                    fhd.write("%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n" %
                              (chrom.decode(),
                               peak['start'],
                               peak['end'],
                               peakname,
                               int(10*peak[score_column]),
                               peak['fc'],
                               peak['pscore'],
                               peak['qscore'],
                               s))
        return

    @cython.ccall
    def write_to_xls(self, ofhd,
                     name_prefix: bytes = b"%s_peak_",
                     name: bytes = b"MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.


        wait... why I have two write_to_xls in this class?

        """
        peakprefix: bytes
        chrs: list
        these_peaks: list
        n_peak: cython.int
        i: cython.int

        write = ofhd.write
        write("\t".join(("chr", "start", "end",  "length",  "abs_summit", "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")

        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix

        peaks = self.peaks
        chrs = list(peaks.keys())
        n_peak = 0
        for chrom in sorted(chrs):
            for end, group in groupby(peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1:
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix.decode(), n_peak, subpeak_letters(i))
                        # [start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                        write("%s\t%d\t%d\t%d" % (chrom.decode(),
                                                  peak['start']+1,
                                                  peak['end'],
                                                  peak['length']))
                        write("\t%d" % (peak['summit']+1))  # summit position
                        write("\t%.6g" % (round(peak['pileup'], 2)))  # pileup height at summit
                        write("\t%.6g" % (peak['pscore']))  # -log10pvalue at summit
                        write("\t%.6g" % (peak['fc']))  # fold change at summit
                        write("\t%.6g" % (peak['qscore']))  # -log10qvalue at summit
                        write("\t%s" % peakname)
                        write("\n")
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix.decode(), n_peak)
                    # [start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                    write("%s\t%d\t%d\t%d" % (chrom.decode(),
                                              peak['start']+1,
                                              peak['end'],
                                              peak['length']))
                    write("\t%d" % (peak['summit']+1))  # summit position
                    write("\t%.6g" % (round(peak['pileup'], 2)))  # pileup height at summit
                    write("\t%.6g" % (peak['pscore']))  # -log10pvalue at summit
                    write("\t%.6g" % (peak['fc']))  # fold change at summit
                    write("\t%.6g" % (peak['qscore']))  # -log10qvalue at summit
                    write("\t%s" % peakname)
                    write("\n")
        return

    @cython.ccall
    def exclude(self, peaksio2: object):
        """ Remove overlapping peaks in peaksio2, another PeakIO object.

        """
        peaks1: dict
        peaks2: dict
        chrs1: list
        chrs2: list
        k: bytes
        ret_peaks: dict
        overlap_found: bool
        r1: PeakContent
        r2: PeakContent
        n_rl1: cython.long
        n_rl2: cython.long

        self.sort()
        peaks1 = self.peaks
        self.total = 0
        assert isinstance(peaksio2, PeakIO)
        peaksio2.sort()
        peaks2 = peaksio2.peaks

        ret_peaks = dict()
        chrs1 = list(peaks1.keys())
        chrs2 = list(peaks2.keys())
        for k in chrs1:
            #print(f"chromosome {k}")
            if not chrs2.count(k):
                # no such chromosome in peaks1, then don't touch the peaks in this chromosome
                ret_peaks[k] = peaks1[k]
                continue
            ret_peaks[k] = []
            n_rl1 = len(peaks1[k])
            n_rl2 = len(peaks2[k])
            rl1_k = iter(peaks1[k]).__next__
            rl2_k = iter(peaks2[k]).__next__
            overlap_found = False
            r1 = rl1_k()
            n_rl1 -= 1
            r2 = rl2_k()
            n_rl2 -= 1
            while (True):
                # we do this until there is no r1 or r2 left.
                if r2["start"] < r1["end"] and r1["start"] < r2["end"]:
                    # since we found an overlap, r1 will be skipped/excluded
                    # and move to the next r1
                    overlap_found = True
                    # print(f"found overlap of {r1['start']} {r1['end']} and {r2['start']} {r2['end']}, move to the next r1")
                    n_rl1 -= 1
                    if n_rl1 >= 0:
                        r1 = rl1_k()
                        # print(f"move to next r1 {r1['start']} {r1['end']}")
                        overlap_found = False
                        continue
                    else:
                        break
                if r1["end"] < r2["end"]:
                    # print(f"now we need to move r1 {r1['start']} {r1['end']}")
                    # in this case, we need to move to the next r1,
                    # we will check if overlap_found is true, if not, we put r1 in a new dict
                    if not overlap_found:
                        # print(f"we add this r1 {r1['start']} {r1['end']} to list")
                        ret_peaks[k].append(r1)
                    n_rl1 -= 1
                    if n_rl1 >= 0:
                        r1 = rl1_k()
                        # print(f"move to next r1 {r1['start']} {r1['end']}")
                        overlap_found = False
                    else:
                        # no more r1 left
                        break
                else:
                    # in this case, we need to move the next r2
                    if n_rl2:
                        r2 = rl2_k()
                        n_rl2 -= 1
                        # print(f"move to next r2 {r2['start']} {r2['end']}")
                    else:
                        # no more r2 left
                        break
            # add the rest of r1
            # print( f"n_rl1: {n_rl1} n_rl2:{n_rl2} last overlap_found is {overlap_found}" )
            # if overlap_found:
            #    n_rl1 -= 1
            if n_rl1 >= 0:
                ret_peaks[k].extend(peaks1[k][-n_rl1-1:])

        for k in ret_peaks.keys():
            self.total += len(ret_peaks[k])

        self.peaks = ret_peaks
        self.CO_sorted = True
        self.sort()
        return

    @cython.ccall
    def read_from_xls(self, ofhd):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.

        """
        line: bytes = b''
        chrom: bytes = b''
        start: cython.int
        end: cython.int
        length: cython.int
        summit: cython.int
        pileup: cython.float
        pscore: cython.float
        fc: cython.float
        qscore: cython.float
        fields: list

        while True:
            if not (line.startswith('#') or line.strip() == ''):
                break
            line = ofhd.readline()

        # sanity check
        columns = line.rstrip().split('\t')
        for a, b in zip(columns, ("chr", "start", "end",  "length", "abs_summit",
                                  "pileup", "-log10(pvalue)", "fold_enrichment",
                                  "-log10(qvalue)", "name")):
            if not a == b:
                raise NotImplementedError('column %s not recognized', a)

        add = self.add
        split = str.split
        rstrip = str.rstrip
        for i, line in enumerate(ofhd.readlines()):
            fields = split(line, '\t')
            chrom = fields[0].encode()
            start = int(fields[1]) - 1
            end = int(fields[2])
            length = int(fields[3])
            if end - start != length:
                raise UserWarning('Malformed peak at line %d:\n%s' % (i, line))
            summit = int(fields[4]) - 1
            pileup = float(fields[5])
            pscore = float(fields[6])
            fc = float(fields[7])
            qscore = float(fields[8])
            peakname = rstrip(fields[9])
            add(chrom, start, end, summit, qscore, pileup, pscore, fc, qscore,
                peakname)


@cython.cclass
class RegionIO:
    """For plain region of chrom, start and end
    """
    regions: dict
    __flag_sorted: bool

    def __init__(self):
        self.regions = {}
        self.__flag_sorted = False

    @cython.ccall
    def add_loc(self, chrom: bytes, start: cython.int, end: cython.int):
        if self.regions.has_key(chrom):
            self.regions[chrom].append((start, end))
        else:
            self.regions[chrom] = [(start, end), ]
        self.__flag_sorted = False
        return

    @cython.ccall
    def sort(self):
        chrom: bytes

        for chrom in sorted(list(self.regions.keys())):
            self.regions[chrom].sort()
        self.__flag_sorted = True

    @cython.ccall
    def get_chr_names(self) -> set:
        return set(sorted(self.regions.keys()))

    @cython.ccall
    def merge_overlap(self):
        """
        merge overlapping regions
        """
        chrom: bytes
        s_new_region: cython.int
        e_new_region: cython.int
        i: cython.int
        regions: dict
        new_regions: dict
        chrs: list
        regions_chr: list
        prev_region: tuple

        if not self.__flag_sorted:
            self.sort()
        regions = self.regions
        new_regions = {}
        chrs = sorted(list(regions.keys()))
        for i in range(len(chrs)):
            chrom = chrs[i]
            new_regions[chrom] = []
            n_append = new_regions[chrom].append
            prev_region = None
            regions_chr = regions[chrom]
            for i in range(len(regions_chr)):
                if not prev_region:
                    prev_region = regions_chr[i]
                    continue
                else:
                    if regions_chr[i][0] <= prev_region[1]:
                        s_new_region = prev_region[0]
                        e_new_region = regions_chr[i][1]
                        prev_region = (s_new_region, e_new_region)
                    else:
                        n_append(prev_region)
                        prev_region = regions_chr[i]
            if prev_region:
                n_append(prev_region)
        self.regions = new_regions
        self.sort()
        return

    @cython.ccall
    def write_to_bed(self, fhd):
        i: cython.int
        chrom: bytes
        chrs: list
        region: tuple

        chrs = sorted(list(self.regions.keys()))
        for i in range(len(chrs)):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                fhd.write("%s\t%d\t%d\n" % (chrom.decode(),
                                            region[0],
                                            region[1]))


@cython.cclass
class BroadPeakContent:
    start: cython.int
    end: cython.int
    length: cython.int
    score: cython.float
    thickStart: bytes
    thickEnd: bytes
    blockNum: cython.int
    blockSizes: bytes
    blockStarts: bytes
    pileup: cython.float
    pscore: cython.float
    fc: cython.float
    qscore: cython.float
    name: bytes

    def __init__(self,
                 start: cython.int,
                 end: cython.int,
                 score: cython.float,
                 thickStart: bytes,
                 thickEnd: bytes,
                 blockNum: cython.int,
                 blockSizes: bytes,
                 blockStarts: bytes,
                 pileup: cython.float,
                 pscore: cython.float,
                 fold_change: cython.float,
                 qscore: cython.float,
                 name: bytes = b"NA"):
        self.start = start
        self.end = end
        self.score = score
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.blockNum = blockNum
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts
        self.length = end - start
        self.pileup = pileup
        self.pscore = pscore
        self.fc = fold_change
        self.qscore = qscore
        self.name = name

    def __getitem__(self, a):
        if a == "start":
            return self.start
        elif a == "end":
            return self.end
        elif a == "length":
            return self.length
        elif a == "score":
            return self.score
        elif a == "thickStart":
            return self.thickStart
        elif a == "thickEnd":
            return self.thickEnd
        elif a == "blockNum":
            return self.blockNum
        elif a == "blockSizes":
            return self.blockSizes
        elif a == "blockStarts":
            return self.blockStarts
        elif a == "pileup":
            return self.pileup
        elif a == "pscore":
            return self.pscore
        elif a == "fc":
            return self.fc
        elif a == "qscore":
            return self.qscore
        elif a == "name":
            return self.name

    def __str__(self):
        return "start:%d;end:%d;score:%f" % (self.start, self.end, self.score)


@cython.cclass
class BroadPeakIO:
    """IO for broad peak information.

    """
    peaks = cython.declare(dict, visibility="public")

    def __init__(self):
        self.peaks = {}

    @cython.ccall
    def add(self,
            chromosome: bytes,
            start: cython.int,
            end: cython.int,
            score: cython.float = 0.0,
            thickStart: bytes = b".",
            thickEnd: bytes = b".",
            blockNum: cython.int = 0,
            blockSizes: bytes = b".",
            blockStarts: bytes = b".",
            pileup: cython.float = 0,
            pscore: cython.float = 0,
            fold_change: cython.float = 0,
            qscore: cython.float = 0,
            name: bytes = b"NA"):
        """items
        chromosome : chromosome name,
        start      : broad region start,
        end        : broad region end,
        score      : average score in all blocks,
        thickStart : start of highly enriched region, # could be b'.'
        thickEnd   : end of highly enriched region,   # could be b'.'
        blockNum   : number of blocks,                # could be 0
        blockSizes : sizes of blocks,                 # could be b'.'
        blockStarts: starts of blocks                 # could be b'.'
        pileup     : median pileup in region          # could be 0
        pscore     : median pvalue score in region    # could be 0
        fold_change: median fold change in region     # could be 0
        qscore     : median pvalue score in region    # could be 0
        name       : peak name                        # could be b'NA'
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome] = []
        self.peaks[chromosome].append(BroadPeakContent(start,
                                                       end,
                                                       score,
                                                       thickStart,
                                                       thickEnd,
                                                       blockNum,
                                                       blockSizes,
                                                       blockStarts,
                                                       pileup,
                                                       pscore,
                                                       fold_change,
                                                       qscore,
                                                       name))

    @cython.ccall
    def filter_pscore(self, pscore_cut: cython.float):
        chrom: bytes
        peaks: dict
        new_peaks: dict
        chrs: list

        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())

        for chrom in sorted(chrs):
            new_peaks[chrom] = [p for p in peaks[chrom] if p['pscore'] >= pscore_cut]
        self.peaks = new_peaks

    @cython.ccall
    def filter_qscore(self, qscore_cut: cython.float):
        chrom: bytes
        peaks: dict
        new_peaks: dict
        chrs: list

        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())

        for chrom in sorted(chrs):
            new_peaks[chrom] = [p for p in peaks[chrom] if p['qscore'] >= qscore_cut]
        self.peaks = new_peaks

    @cython.ccall
    def filter_fc(self, fc_low: float, fc_up: float = -1):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in
        [fc_low,fc_up)

        fc_up has to be a positive number, otherwise it won't be
        applied.

        """
        chrom: bytes
        peaks: dict
        new_peaks: dict
        chrs: list

        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())
        if fc_up >= 0:
            for chrom in sorted(chrs):
                new_peaks[chrom] = [p for p in peaks[chrom] if p['fc'] >= fc_low and p['fc'] < fc_up]
        else:
            for chrom in sorted(chrs):
                new_peaks[chrom] = [p for p in peaks[chrom] if p['fc'] >= fc_low]
        self.peaks = new_peaks

    @cython.ccall
    def total(self):
        chrom: bytes
        peaks: dict
        chrs: list
        x: cython.long = 0

        peaks = self.peaks
        chrs = list(peaks.keys())
        for chrom in sorted(chrs):
            x += len(peaks[chrom])
        return x

    @cython.ccall
    def write_to_gappedPeak(self, fhd,
                            name_prefix: bytes = b"peak_",
                            name: bytes = b'peak',
                            description: bytes = b"%s",
                            score_column: str = "score",
                            trackline: bool = True):
        """Print out peaks in gappedBed format. Only those with stronger enrichment regions are saved.

        This format is basically a BED12+3 format.

        +--------------+------+----------------------------------------+
        |field         |type  |description                             |
        +--------------+------+----------------------------------------+
        |chrom         |string|Name of the chromosome                  |
        +--------------+------+----------------------------------------+
        |chromStart    |int   |The starting position of the feature in |
        |              |      |the chromosome. The first base in a     |
        |              |      |chromosome is numbered 0.               |
        +--------------+------+----------------------------------------+
        |chromEnd      |int   |The ending position of the feature in   |
        |              |      |the chromosome or scaffold. The chromEnd|
        |              |      |base is not included in the display of  |
        |              |      |the feature.  For example, the first 100|
        |              |      |bases of a chromosome are defined as    |
        |              |      |chromStart=0, chromEnd=100, and span the|
        |              |      |bases numbered 0-99.                    |
        +--------------+------+----------------------------------------+
        |name          |string|Name given to a region (preferably      |
        |              |      |unique). Use '.' if no name is assigned.|
        +--------------+------+----------------------------------------+
        |score         |int   |Indicates how dark the peak will be     |
        |(always use   |      |displayed in the browser (1-1000). If   |
        |1000 for      |      |'0', the DCC will assign this based on  |
        |the           |      |signal value. Ideally average           |
        |thickest      |      |signalValue per base spread between     |
        |color)        |      |100-1000.                               |
        +--------------+------+----------------------------------------+
        |strand        |char  |+/- to denote strand or orientation     |
        |(always .)    |      |(whenever applicable). Use '.' if no    |
        |              |      |orientation is assigned.                |
        +--------------+------+----------------------------------------+
        |thickStart    |int   | The starting position at which the     |
        |              |      |feature is drawn thickly. Mark the start|
        |              |      |of highly enriched regions. Not used in |
        |              |      |gappedPeak, so set to 0                 |
        +--------------+------+----------------------------------------+
        |thickEnd      |int   | The ending position at which the       |
        |              |      |feature is drawn thickly. Mark the end  |
        |              |      |of highly enriched regions. Not used, so|
        |              |      |set as 0.                               |
        +--------------+------+----------------------------------------+
        |itemRGB       |string| Not used. Set it as 0.                 |
        +--------------+------+----------------------------------------+
        |blockCounts   |int   | The number of blocks (exons) in the BED|
        |              |      |line.                                   |
        +--------------+------+----------------------------------------+
        |blockSizes    |string| A comma-separated list of the block    |
        |              |      |sizes.                                  |
        +--------------+------+----------------------------------------+
        |blockStarts   |string| A comma-separated list of block starts.|
        +--------------+------+----------------------------------------+
        |signalValue   |float |Measurement of overall (usually,        |
        |(fc)          |      |average) enrichment for the region.     |
        +--------------+------+----------------------------------------+
        |pValue        |float |Measurement of statistical signficance  |
        |              |      |(-log10). Use -1 if no pValue is        |
        |              |      |assigned.                               |
        +--------------+------+----------------------------------------+
        |qValue        |float |Measurement of statistical significance |
        |              |      |using false discovery rate. Use -1 if no|
        |              |      |qValue is assigned.                     |
        +--------------+------+----------------------------------------+

        """
        chrs: list
        n_peak: cython.int = 0
        peak: BroadPeakContent
        desc: bytes
        peakprefix: bytes
        chrom: bytes

        chrs = list(self.peaks.keys())
        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix
        try:
            desc = description % name
        except Exception:
            desc = description
        if trackline:
            fhd.write("track name=\"%s\" description=\"%s\" type=gappedPeak nextItemButton=on\n" % (name.decode(), desc.decode()))
        for chrom in sorted(chrs):
            for peak in self.peaks[chrom]:
                n_peak += 1
                if peak["thickStart"] != b".":
                    fhd.write("%s\t%d\t%d\t%s%d\t%d\t.\t0\t0\t0\t%d\t%s\t%s\t%.6g\t%.6g\t%.6g\n" %
                              (chrom.decode(),
                               peak["start"],
                               peak["end"],
                               peakprefix.decode(),
                               n_peak,
                               int(10*peak[score_column]),
                               peak["blockNum"],
                               peak["blockSizes"].decode(),
                               peak["blockStarts"].decode(),
                               peak['fc'],
                               peak['pscore'],
                               peak['qscore']))

    @cython.ccall
    def write_to_Bed12(self, fhd,
                       name_prefix: bytes = b"peak_",
                       name: bytes = b'peak',
                       description: bytes = b"%s",
                       score_column: str = "score",
                       trackline: bool = True):
        """Print out peaks in Bed12 format.

        +--------------+------+----------------------------------------+
        |field         |type  |description                             |
        +--------------+------+----------------------------------------+
        |chrom         |string|Name of the chromosome                  |
        +--------------+------+----------------------------------------+
        |chromStart    |int   |The starting position of the feature in |
        |              |      |the chromosome. The first base in a     |
        |              |      |chromosome is numbered 0.               |
        +--------------+------+----------------------------------------+
        |chromEnd      |int   |The ending position of the feature in   |
        |              |      |the chromosome or scaffold. The chromEnd|
        |              |      |base is not included in the display of  |
        |              |      |the feature.  For example, the first 100|
        |              |      |bases of a chromosome are defined as    |
        |              |      |chromStart=0, chromEnd=100, and span the|
        |              |      |bases numbered 0-99.                    |
        +--------------+------+----------------------------------------+
        |name          |string|Name given to a region (preferably      |
        |              |      |unique). Use '.' if no name is assigned.|
        +--------------+------+----------------------------------------+
        |score         |int   |Indicates how dark the peak will be     |
        |(always use   |      |displayed in the browser (1-1000). If   |
        |1000 for      |      |'0', the DCC will assign this based on  |
        |the           |      |signal value. Ideally average           |
        |thickest      |      |signalValue per base spread between     |
        |color)        |      |100-1000.                               |
        +--------------+------+----------------------------------------+
        |strand        |char  |+/- to denote strand or orientation     |
        |(always .)    |      |(whenever applicable). Use '.' if no    |
        |              |      |orientation is assigned.                |
        +--------------+------+----------------------------------------+
        |thickStart    |int   | The starting position at which the     |
        |              |      |feature is drawn thickly. Mark the start|
        |              |      |of highly enriched regions.             |
        |              |      |                                        |
        +--------------+------+----------------------------------------+
        |thickEnd      |int   | The ending position at which the       |
        |              |      |feature is drawn thickly. Mark the end  |
        |              |      |of highly enriched regions.             |
        +--------------+------+----------------------------------------+
        |itemRGB       |string| Not used. Set it as 0.                 |
        +--------------+------+----------------------------------------+
        |blockCounts   |int   | The number of blocks (exons) in the BED|
        |              |      |line.                                   |
        +--------------+------+----------------------------------------+
        |blockSizes    |string| A comma-separated list of the block    |
        |              |      |sizes.                                  |
        +--------------+------+----------------------------------------+
        |blockStarts   |string| A comma-separated list of block starts.|
        +--------------+------+----------------------------------------+

        """
        chrs: list
        n_peak: cython.int = 0
        peakprefix: bytes
        peak: BroadPeakContent
        desc: bytes
        peakprefix: bytes
        chrom: bytes

        chrs = list(self.peaks.keys())
        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix
        try:
            desc = description % name
        except Exception:
            desc = description
        if trackline:
            fhd.write("track name=\"%s\" description=\"%s\" type=bed nextItemButton=on\n" % (name.decode(), desc.decode()))
        for chrom in sorted(chrs):
            for peak in self.peaks[chrom]:
                n_peak += 1
                if peak["thickStart"] == b".":
                    # this will violate gappedPeak format, since it's a complement like broadPeak line.
                    fhd.write("%s\t%d\t%d\t%s%d\t%d\t.\n" %
                              (chrom.decode(),
                               peak["start"],
                               peak["end"],
                               peakprefix.decode(),
                               n_peak,
                               int(10*peak[score_column])))
                else:
                    fhd.write("%s\t%d\t%d\t%s%d\t%d\t.\t%s\t%s\t0\t%d\t%s\t%s\n" %
                              (chrom.decode(),
                               peak["start"],
                               peak["end"],
                               peakprefix.decode(),
                               n_peak,
                               int(10*peak[score_column]),
                               peak["thickStart"].decode(),
                               peak["thickEnd"].decode(),
                               peak["blockNum"],
                               peak["blockSizes"].decode(),
                               peak["blockStarts"].decode()))

    @cython.ccall
    def write_to_broadPeak(self, fhd,
                           name_prefix: bytes = b"peak_",
                           name: bytes = b'peak',
                           description: bytes = b"%s",
                           score_column: str = "score",
                           trackline: bool = True):
        """Print out peaks in broadPeak format.

        This format is designed for ENCODE project, and basically a
        BED6+3 format.

        +-----------+------+----------------------------------------+
        |field      |type  |description                             |
        +-----------+------+----------------------------------------+
        |chrom      |string|Name of the chromosome                  |
        +-----------+------+----------------------------------------+
        |chromStart |int   |The starting position of the feature in |
        |           |      |the chromosome. The first base in a     |
        |           |      |chromosome is numbered 0.               |
        +-----------+------+----------------------------------------+
        |chromEnd   |int   |The ending position of the feature in   |
        |           |      |the chromosome or scaffold. The chromEnd|
        |           |      |base is not included in the display of  |
        |           |      |the feature.  For example, the first 100|
        |           |      |bases of a chromosome are defined as    |
        |           |      |chromStart=0, chromEnd=100, and span the|
        |           |      |bases numbered 0-99.                    |
        +-----------+------+----------------------------------------+
        |name       |string|Name given to a region (preferably      |
        |           |      |unique). Use '.' if no name is assigned.|
        +-----------+------+----------------------------------------+
        |score      |int   |Indicates how dark the peak will be     |
        |(-logqvalue|      |displayed in the browser (1-1000). If   |
        |in MACS3 * |      |'0', the DCC will assign this based on  |
        |10)        |      |signal value. Ideally average           |
        |           |      |signalValue per base spread between     |
        |           |      |100-1000.                               |
        +-----------+------+----------------------------------------+
        |strand     |char  |+/- to denote strand or orientation     |
        |(always .) |      |(whenever applicable). Use '.' if no    |
        |           |      |orientation is assigned.                |
        +-----------+------+----------------------------------------+
        |signalValue|float |Measurement of overall (usually,        |
        |(fc)       |      |average) enrichment for the region.     |
        +-----------+------+----------------------------------------+
        |pValue     |float |Measurement of statistical signficance  |
        |           |      |(-log10). Use -1 if no pValue is        |
        |           |      |assigned.                               |
        +-----------+------+----------------------------------------+
        |qValue     |float |Measurement of statistical significance |
        |           |      |using false discovery rate. Use -1 if no|
        |           |      |qValue is assigned.                     |
        +-----------+------+----------------------------------------+

        """
        chrs: list
        n_peak: cython.int = 0
        peakprefix: bytes
        peak: BroadPeakContent
        peakprefix: bytes
        chrom: bytes
        peakname: str

        chrs = list(self.peaks.keys())
        write = fhd.write
        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix
        if trackline:
            write("track type=broadPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name.decode(), name.decode()))
        for chrom in sorted(chrs):
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                peak = these_peaks[0]
                peakname = "%s%d" % (peakprefix.decode(), n_peak)
                fhd.write("%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\n" %
                          (chrom.decode(),
                           peak['start'],
                           peak['end'],
                           peakname,
                           int(10*peak[score_column]),
                           peak['fc'],
                           peak['pscore'],
                           peak['qscore']))
        return

    @cython.ccall
    def write_to_xls(self, ofhd,
                     name_prefix: bytes = b"%s_peak_",
                     name: bytes = b"MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.


        wait... why I have two write_to_xls in this class?

        """
        chrom: bytes
        chrs: list
        peakprefix: bytes
        peaks: dict
        these_peaks: list
        peak: BroadPeakContent
        peakname: str

        write = ofhd.write
        write("\t".join(("chr", "start", "end",  "length",  "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")

        try:
            peakprefix = name_prefix % name
        except Exception:
            peakprefix = name_prefix

        peaks = self.peaks
        chrs = list(peaks.keys())
        n_peak = 0
        for chrom in sorted(chrs):
            for end, group in groupby(peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                peak = these_peaks[0]
                peakname = "%s%d" % (peakprefix.decode(), n_peak)
                write("%s\t%d\t%d\t%d" % (chrom.decode(),
                                          peak['start']+1,
                                          peak['end'],
                                          peak['length']))
                write("\t%.6g" % (round(peak['pileup'], 2)))  # pileup height at summit
                write("\t%.6g" % (peak['pscore']))  # -log10pvalue at summit
                write("\t%.6g" % (peak['fc']))  # fold change at summit
                write("\t%.6g" % (peak['qscore']))  # -log10qvalue at summit
                write("\t%s" % peakname)
                write("\n")
        return
