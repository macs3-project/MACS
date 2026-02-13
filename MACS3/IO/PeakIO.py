# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-11-14 16:52:15 Tao Liu>

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
    """Return the alphabetical label for a zero-based subpeak index.

    Args:
        i: Zero-based subpeak index.

    Returns:
        str: Alphabetical label sequence (``a``, ``b``, ..., ``aa``).
    """
    if i < 26:
        return chr(97+i)
    else:
        return subpeak_letters(i // 26) + chr(97 + (i % 26))

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class PeakContent:
    """Represent a narrow peak and its derived statistics.

    Attributes:
        chrom: Chromosome name (bytes).
        start: 0-based inclusive start coordinate.
        end: 0-based exclusive end coordinate.
        length: Peak length in base pairs.
        summit: 0-based summit position.
        score: Peak score reported by MACS3.
        pileup: Tag pileup at the summit.
        pscore: ``-log10(pvalue)`` at the summit.
        fc: Fold enrichment at the summit.
        qscore: ``-log10(qvalue)`` at the summit.
        name: Optional peak identifier.
    """
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
        """Initialise a peak record with positional and score metrics.

        Args:
            chrom: Chromosome name encoded as bytes.
            start: 0-based inclusive start coordinate.
            end: 0-based exclusive end coordinate.
            summit: 0-based summit position.
            peak_score: Peak score reported by MACS3.
            pileup: Tag pileup at the summit.
            pscore: ``-log10(pvalue)`` at the summit.
            fold_change: Fold enrichment at the summit.
            qscore: ``-log10(qvalue)`` at the summit.
            name: Optional peak identifier.
        """
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
        """Return a peak attribute by symbolic key.

        Args:
            a: Attribute name (for example ``"chrom"`` or ``"score"``).

        Returns:
            Any: Value associated with the requested attribute.
        """
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
        """Assign a new value to a peak attribute.

        Args:
            a: Attribute name to set.
            v: Replacement value.
        """
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
        """Return a concise textual description of the peak."""
        return "chrom:%s;start:%d;end:%d;score:%f" % (self.chrom,
                                                      self.start,
                                                      self.end,
                                                      self.score)

    def __getstate__(self):
        """Serialise the peak content for pickling.

        Returns:
            tuple: Ordered tuple representing the peak state.
        """
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
        """Restore the peak content from a previously serialised tuple.

        Args:
            state: Tuple produced by :meth:`__getstate__`.
        """
        (self.chrom, self.start, self.end, self.length, self.summit,
         self.score, self.pileup, self.pscore, self.fc,
         self.qscore, self.name) = state


@cython.cclass
class PeakIO:
    """Manage in-memory collections of narrow peak intervals.

    Attributes:
        peaks: Dictionary storing peak contents.
        CO_sorted: Whether peaks have been coordinate-sorted.
        total: Total number of peaks.
        name: Collection name used in output.
    """
    # dictionary storing peak contents
    peaks = cython.declare(dict, visibility="public")
    # whether peaks have been sorted by coordinations
    CO_sorted = cython.declare(bool, visibility="public")
    # total number of peaks
    total = cython.declare(cython.long, visibility="public")
    # name
    name = cython.declare(bytes, visibility="public")

    def __init__(self, name: bytes = b"MACS3"):
        """Initialise an empty peak collection."""
        self.peaks = {}
        self.CO_sorted = False
        self.total = 0
        self.name = name

    @cython.ccall
    def add(self,
            chromosome: bytes,
            start: cython.int,
            end: cython.int,
            summit: cython.int = 0,
            peak_score: cython.float = 0,
            pileup: cython.float = 0,
            pscore: cython.float = 0,
            fold_change: cython.float = 0,
            qscore: cython.float = 0,
            name: bytes = b""):
        """Add a peak described by raw coordinates and scores.

        Args:
            chromosome: Chromosome name for the peak.
            start: 0-based inclusive start coordinate.
            end: 0-based exclusive end coordinate.
            summit: 0-based summit position.
            peak_score: Reported peak score.
            pileup: Tag pileup at the summit.
            pscore: ``-log10(pvalue)`` score.
            fold_change: Fold enrichment relative to control.
            qscore: ``-log10(qvalue)`` score.
            name: Optional peak identifier.

        Examples:
            .. code-block:: python

                from MACS3.IO.PeakIO import PeakIO
                peaks = PeakIO()
                peaks.add(b"chr1", 100, 200, summit=150, peak_score=10.0)
        """
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
        """Extend the collection with an existing :class:`PeakContent`.

        Args:
            chromosome: Chromosome name under which to store the peak.
            peakcontent: Peak record to append.

        Examples:
            .. code-block:: python

                from MACS3.IO.PeakIO import PeakIO, PeakContent
                peaks = PeakIO()
                peak = PeakContent(b"chr1", 100, 200, 150, 10.0, 5.0, 3.0, 2.0, 1.0)
                peaks.add_PeakContent(b"chr1", peak)
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome] = []
        self.peaks[chromosome].append(peakcontent)
        self.total += 1
        self.CO_sorted = False

    @cython.ccall
    def get_data_from_chrom(self, chrom: bytes) -> list:
        """Return peaks for ``chrom``, initialising storage if needed.

        Args:
            chrom: Chromosome name to query.

        Returns:
            list: Peaks associated with ``chrom``.

        Examples:
            .. code-block:: python

                peaks = PeakIO()
                chrom_peaks = peaks.get_data_from_chrom(b"chr1")
        """
        if not self.peaks.has_key(chrom):
            self.peaks[chrom] = []
        return self.peaks[chrom]

    def get_chr_names(self) -> set:
        """Return the chromosome names represented in the collection.

        Returns:
            set: Unique chromosome names.

        Examples:
            .. code-block:: python

                peaks = PeakIO()
                names = peaks.get_chr_names()
        """
        return set(self.peaks.keys())

    def sort(self):
        """Sort peaks on each chromosome by ascending start position."""
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
        """Return a new ``PeakIO`` containing ``n`` randomly sampled peaks.

        Args:
            n: Number of peaks to sample.
            seed: RNG seed to ensure reproducibility.

        Returns:
            PeakIO: Fresh instance populated with sampled peaks.

        Examples:
            .. code-block:: python

                sampled = peaks.randomly_pick(100, seed=42)
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
        """Filter peaks by minimum ``-log10(pvalue)``.

        Args:
            pscore_cut: Lower bound (inclusive) for ``-log10(pvalue)``.
        """
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
        """Filter peaks by minimum ``-log10(qvalue)``.

        Args:
            qscore_cut: Lower bound (inclusive) for ``-log10(qvalue)``.
        """
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
        """Filter peaks by fold-change range.

        Args:
            fc_low: Inclusive lower bound on fold change.
            fc_up: Exclusive upper bound; ignored if ``<= 0``.
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
        """Filter peaks by their primary score range.

        Args:
            lower_score: Inclusive lower bound on score.
            upper_score: Exclusive upper bound; if ``<= 0`` the bound is ignored.
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
        """Return a debug-friendly representation of all stored peaks."""
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
        """Render peaks in BED5-compatible format via a callback.

        Args:
            name_prefix: Template used to build peak names.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for optional header line.
            score_column: Peak attribute to emit as the score field.
            trackline: Whether to emit a UCSC ``track`` header line.
            print_func: Callable accepting pre-formatted lines.
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
        """Render peak summits in BED5 format via a callback.

        Args:
            name_prefix: Template used to build summit names.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for optional header line.
            score_column: Peak attribute to emit as the score field.
            trackline: Whether to emit a UCSC ``track`` header line.
            print_func: Callable accepting pre-formatted lines.
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
        """Write peaks in BED5 format to ``stdout``.

        The five columns correspond to chromosome, start, end, name, and
        the attribute selected by ``score_column``.
        """
        return self._to_bed(name_prefix=b"%s_peak_", score_column="score", name=self.name, description=b"")

    def to_summits_bed(self):
        """Write peak summits in BED5 format to ``stdout``.

        Each summit is emitted as a one-base interval with the selected score column.
        """
        return self._to_summits_bed(name_prefix=b"%s_peak_", score_column="score", name=self.name, description=b"")

    # these methods are very fast, specifying types is unnecessary
    def write_to_bed(self, fhd,
                     name_prefix: bytes = b"%s_peak_",
                     name: bytes = b"MACS",
                     description: bytes = b"%s",
                     score_column: str = "score",
                     trackline: bool = True):
        """Write peaks to a file handle in BED5 format.

        Args:
            fhd: Writable file-like object.
            name_prefix: Template used to build peak names.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for optional header line.
            score_column: Peak attribute to emit as the score field.
            trackline: Whether to emit a UCSC ``track`` header line.

        Examples:
            .. code-block:: python

                with open("peaks.bed", "w") as f:
                    peaks.write_to_bed(f)
        """
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
        """Write peak summits to a file handle in BED5 format.

        Args:
            fhd: Writable file-like object.
            name_prefix: Template used to build summit names.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for optional header line.
            score_column: Peak attribute to emit as the score field.
            trackline: Whether to emit a UCSC ``track`` header line.

        Examples:
            .. code-block:: python

                with open("summits.bed", "w") as f:
                    peaks.write_to_summit_bed(f)
        """
        return self._to_summits_bed(name_prefix=name_prefix, name=name,
                                    description=description, score_column=score_column,
                                    print_func=fhd.write, trackline=trackline)

    def write_to_narrowPeak(self, fhd,
                            name_prefix: bytes = b"%s_peak_",
                            name: bytes = b"MACS",
                            score_column: str = "score",
                            trackline: bool = False):
        """Write peaks in the ENCODE narrowPeak (BED6+4) format.

        Args:
            fhd: Writable file-like object.
            name_prefix: Template used to construct peak identifiers.
            name: Dataset label interpolated into ``name_prefix``.
            score_column: Peak attribute mapped to the narrowPeak score field.
            trackline: Whether to emit a UCSC ``track`` header.

        Examples:
            .. code-block:: python

                with open("peaks.narrowPeak", "w") as f:
                    peaks.write_to_narrowPeak(f)
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
        """Export narrow peaks to a tab-delimited ``.xls`` text file.

        Args:
            ofhd: Writable file-like object.
            name_prefix: Template used to build peak identifiers.
            name: Dataset label interpolated into ``name_prefix``.

        Examples:
            .. code-block:: python

                with open("peaks.xls", "w") as f:
                    peaks.write_to_xls(f)
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
        """Remove peaks that overlap any entry in ``peaksio2``.

        Args:
            peaksio2: Another :class:`PeakIO` instance providing exclusion regions.
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
        """Load peak records from a MACS3 ``.xls`` tab-delimited report.

        Args:
            ofhd: Readable file-like object positioned at the beginning of the report.
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
    """Helper for storing and manipulating simple genomic regions.

    Attributes:
        regions: Mapping of chromosome name to a list of ``(start, end)`` tuples.
        __flag_sorted: Whether per-chromosome regions are sorted.
    """
    regions: dict
    __flag_sorted: bool

    def __init__(self):
        """Initialise an empty region collection."""
        self.regions = {}
        self.__flag_sorted = False

    @cython.ccall
    def add_loc(self, chrom: bytes, start: cython.int, end: cython.int):
        """Append a new ``(start, end)`` interval for ``chrom``.

        Examples:
            .. code-block:: python

                regions = RegionIO()
                regions.add_loc(b"chr1", 100, 200)
        """
        if self.regions.has_key(chrom):
            self.regions[chrom].append((start, end))
        else:
            self.regions[chrom] = [(start, end), ]
        self.__flag_sorted = False
        return

    @cython.ccall
    def sort(self):
        """Sort regions for each chromosome by their start coordinate. """
        chrom: bytes

        for chrom in sorted(list(self.regions.keys())):
            self.regions[chrom].sort()
        self.__flag_sorted = True

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return chromosome names present in the region set."""
        return set(sorted(self.regions.keys()))

    @cython.ccall
    def merge_overlap(self):
        """Merge overlapping intervals within each chromosome."""
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
        """Emit regions in BED format to the provided file-like object.

        Examples:
            .. code-block:: python

                with open("regions.bed", "w") as f:
                    regions.write_to_bed(f)
        """
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
    """Container for broad peak metadata used in broadPeak format."""
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
                 name: bytes = b"MACS3"):
        """Initialise a broad peak record with block structure and scores."""
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
        """Provide dict-like read access to stored attributes."""
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
        """Return a compact summary describing the broad peak."""
        return "start:%d;end:%d;score:%f" % (self.start, self.end, self.score)


@cython.cclass
class BroadPeakIO:
    """IO for broad peak information.

    Attributes:
        peaks: Mapping of chromosome name to :class:`BroadPeakContent` list.
    """
    peaks = cython.declare(dict, visibility="public")

    def __init__(self):
        """Create an empty container for :class:`BroadPeakContent` objects."""
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
        """Append a :class:`BroadPeakContent` record.

        Args:
            chromosome: Chromosome name for the region.
            start: 0-based inclusive start coordinate.
            end: 0-based exclusive end coordinate.
            score: Average score across blocks.
            thickStart: Start of the high-enrichment segment or ``b'.'``.
            thickEnd: End of the high-enrichment segment or ``b'.'``.
            blockNum: Number of sub-blocks composing the region.
            blockSizes: Comma-separated block sizes as bytes.
            blockStarts: Comma-separated block starts as bytes.
            pileup: Median pileup within the region.
            pscore: Median ``-log10(pvalue)``.
            fold_change: Median fold-change value.
            qscore: Median ``-log10(qvalue)``.
            name: Optional region identifier.

        Examples:
            .. code-block:: python

                from MACS3.IO.PeakIO import BroadPeakIO
                peaks = BroadPeakIO()
                peaks.add(b"chr1", 100, 500, score=10.0, blockNum=1,
                          blockSizes=b"400", blockStarts=b"0",...)
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
        """Retain broad peaks with ``-log10(pvalue)`` ≥ ``pscore_cut``.

        Args:
            pscore_cut: Inclusive lower bound for ``-log10(pvalue)``."""
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
        """Retain broad peaks with ``-log10(qvalue)`` ≥ ``qscore_cut``.

        Args:
            qscore_cut: Inclusive lower bound for ``-log10(qvalue)``."""
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
        """Filter broad peaks by fold-change range.

        Args:
            fc_low: Inclusive lower bound on fold change.
            fc_up: Exclusive upper bound; ignored when negative.
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
        """Return the total number of broad peaks currently stored.

        Returns:
            int: Number of broad peaks."""
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
        """Write broad peaks in gappedPeak (BED12+3) format.

        Args:
            fhd: Writable file-like object.
            name_prefix: Template used to construct peak identifiers.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for the optional header.
            score_column: Peak attribute mapped to the score column.
            trackline: Whether to emit a UCSC ``track`` header.

        Examples:
            .. code-block:: python

                with open("broad.gappedPeak", "w") as f:
                    peaks.write_to_gappedPeak(f)
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
        """Write broad peaks in BED12 format.

        Args:
            fhd: Writable file-like object.
            name_prefix: Template used to construct peak identifiers.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for the optional header.
            score_column: Peak attribute mapped to the score column.
            trackline: Whether to emit a UCSC ``track`` header.

        Examples:
            .. code-block:: python

                with open("broad.bed12", "w") as f:
                    peaks.write_to_Bed12(f)
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
        """Write broad peaks in the ENCODE broadPeak (BED6+3) format.

        Args:
            fhd: Writable file-like object.
            name_prefix: Template used to construct peak identifiers.
            name: Dataset label interpolated into ``name_prefix``.
            description: Track description for the optional header.
            score_column: Peak attribute mapped to the score column.
            trackline: Whether to emit a UCSC ``track`` header.

        Examples:
            .. code-block:: python

                with open("broad.broadPeak", "w") as f:
                    peaks.write_to_broadPeak(f)
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
        """Export broad peaks to a tab-delimited ``.xls`` text file.

        Args:
            ofhd: Writable file-like object.
            name_prefix: Template used to build peak identifiers.
            name: Dataset label interpolated into ``name_prefix``.

        Examples:
            .. code-block:: python

                with open("broad.xls", "w") as f:
                    peaks.write_to_xls(f)
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
