# cython: language_level=3
# cython: profile=True
# Time-stamp: <2020-11-24 17:11:29 Tao Liu>

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
import re
import sys

# ------------------------------------
# MACS3 modules
# ------------------------------------

from MACS3.Utilities.Constants import *

# ------------------------------------
# Other modules
# ------------------------------------

from cpython cimport bool

# ------------------------------------
# constants
# ------------------------------------
__version__ = "PeakIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "PeakIO class"

# ------------------------------------
# Misc functions
# ------------------------------------
cdef str subpeak_letters( int i):
    if i < 26:
        return chr(97+i)
    else:
        return subpeak_letters(i // 26) + chr(97 + (i % 26))

# ------------------------------------
# Classes
# ------------------------------------

cdef class PeakContent:
    cdef:
        int start
        int end
        int length
        int summit
        float score
        float pileup
        float pscore
        float fc
        float qscore
        bytes name

    def __init__ ( self, int start, int end, int summit,
                   float peak_score, float pileup,
                   float pscore, float fold_change, float qscore,
                   bytes name= b"NA" ):
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

    def __getitem__ ( self, a ):
        if a == "start":
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

    def __setitem__ ( self, a, v ):
        if a == "start":
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

    def __str__ (self):
        return "start:%d;end:%d;score:%f" % ( self.start, self.end, self.score )

cdef class PeakIO:
    """IO for peak information.

    """
    cdef:
        dict peaks

    def __init__ (self):
        self.peaks = {}

    cpdef add (self, bytes chromosome, int start, int end, int summit = 0,
               float peak_score = 0, float pileup = 0,
               float pscore = 0, float fold_change = 0, float qscore = 0,
               bytes name = b"NA"):
        """items:
        start:start
        end:end,
        length:end-start,
        summit:summit,
        score:peak_score,
        pileup:pileup,
        pscore:pscore,
        fc:fold_change,
        qscore:qscore
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append(PeakContent( start, end, summit, peak_score, pileup, pscore, fold_change, qscore, name))

    cpdef add_PeakContent ( self, bytes chromosome, object peakcontent ):
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append(peakcontent)

    def get_data_from_chrom (self, bytes chrom):
        if not self.peaks.has_key( chrom ):
            self.peaks[chrom]= []
        return self.peaks[chrom]

    cpdef set get_chr_names (self):
        return set(sorted(self.peaks.keys()))

    def sort ( self ):
        # sort by position
        chrs = sorted(list(self.peaks.keys()))
        for chrom in chrs:
            self.peaks[chrom].sort(key=lambda x:x['start'])
        return


    def filter_pscore (self, double pscore_cut ):
        cdef bytes chrom

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(list(peaks.keys()))

        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['pscore'] >= pscore_cut]
        self.peaks = new_peaks

    def filter_qscore (self, double qscore_cut ):
        cdef bytes chrom

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(list(peaks.keys()))

        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['qscore'] >= qscore_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low and p['fc']<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x

    # these methods are very fast, specifying types is unnecessary
    def write_to_xls (self, fhd, bytes name_prefix=b"%s_peak_", bytes name=b"MACS"):
        return self._to_xls(name_prefix=name_prefix, name=name,
                            print_func=fhd.write)

    cdef _to_xls (self, bytes name_prefix=b"%s_peak_", bytes name=b"MACS", print_func=sys.stdout.write):

        if self.peaks:
            print_func("\t".join(("chr","start", "end",  "length",  "abs_summit", "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")
        else:
            return

        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix

        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for end, group in groupby(peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1:
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix.decode(), n_peak, subpeak_letters(i))
                        #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                        print_func("%s\t%d\t%d\t%d" % (chrom.decode(),peak['start']+1,peak['end'],peak['length']))
                        print_func("\t%d" % (peak['summit']+1)) # summit position
                        print_func("\t%.6g" % (round(peak['pileup'],2))) # pileup height at summit
                        print_func("\t%.6g" % (peak['pscore'])) # -log10pvalue at summit
                        print_func("\t%.6g" % (peak['fc'])) # fold change at summit
                        print_func("\t%.6g" % (peak['qscore'])) # -log10qvalue at summit
                        print_func("\t%s" % peakname)
                        print_func("\n")
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix.decode(), n_peak)
                    #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                    print_func("%s\t%d\t%d\t%d" % (chrom.decode(),peak['start']+1,peak['end'],peak['length']))
                    print_func("\t%d" % (peak['summit']+1)) # summit position
                    print_func("\t%.6g" % (round(peak['pileup'],2))) # pileup height at summit
                    print_func("\t%.6g" % (peak['pscore'])) # -log10pvalue at summit
                    print_func("\t%.6g" % (peak['fc'])) # fold change at summit
                    print_func("\t%.6g" % (peak['qscore'])) # -log10qvalue at summit
                    print_func("\t%s" % peakname)
                    print_func("\n")
        return

    cdef _to_bed(self, bytes name_prefix=b"%s_peak_", bytes name=b"MACS",
                bytes description=b"%s", str score_column="score",
                trackline=False, print_func=sys.stdout.write):
        """
        generalization of tobed and write_to_bed
        """
        #print(description)
        chrs = list(self.peaks.keys())
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        if trackline:
            try: print_func('track name="%s (peaks)" description="%s" visibility=1\n' % ( name.replace(b"\"", b"\\\"").decode(),\
                                                                                          desc.replace(b"\"", b"\\\"").decode() ) )
            except: print_func('track name=MACS description=Unknown\n')
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        print_func("%s\t%d\t%d\t%s%d%s\t%.6g\n" % (chrom.decode(),peak['start'],peak['end'],peakprefix.decode(),n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    print_func("%s\t%d\t%d\t%s%d\t%.6g\n" % (chrom.decode(),peak['start'],peak['end'],peakprefix.decode(),n_peak,peak[score_column]))

    cdef _to_summits_bed(self, bytes name_prefix=b"%s_peak_", bytes name=b"MACS",
                        bytes description = b"%s", str score_column="score",
                        bool trackline=False, print_func=sys.stdout.write):
        """
        generalization of to_summits_bed and write_to_summit_bed
        """
        chrs = list(self.peaks.keys())
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        if trackline:
            try: print_func('track name="%s (summits)" description="%s" visibility=1\n' % ( name.replace(b"\"", b"\\\"").decode(),\
                                                                                            desc.replace(b"\"", b"\\\"").decode() ) )
            except: print_func('track name=MACS description=Unknown')
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        summit_p = peak['summit']
                        print_func("%s\t%d\t%d\t%s%d%s\t%.6g\n" % (chrom.decode(),summit_p,summit_p+1,peakprefix.decode(),n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    summit_p = peak['summit']
                    print_func("%s\t%d\t%d\t%s%d\t%.6g\n" % (chrom.decode(),summit_p,summit_p+1,peakprefix.decode(),n_peak,peak[score_column]))

    def tobed (self):
        """Print out peaks in BED5 format.

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
        return self._to_bed(name_prefix=b"peak_", score_column="score", name=b"", description=b"")

    def to_summits_bed (self):
        """Print out peak summits in BED5 format.

        Five columns are chromosome, summit start, summit end, peak name, and peak height.

        """
        return self._to_summits_bed(name_prefix=b"peak_", score_column="score", name=b"", description=b"")

    # these methods are very fast, specifying types is unnecessary
    def write_to_bed (self, fhd, bytes name_prefix=b"peak_", bytes name=b"MACS",
                        bytes description = b"%s", str score_column="score", trackline=True):
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
        #print(description)
        return self._to_bed(name_prefix=name_prefix, name=name,
                            description=description, score_column=score_column,
                            print_func=fhd.write, trackline=trackline)

    def write_to_summit_bed (self, fhd, bytes name_prefix = b"peak_", bytes name = b"MACS",
                             bytes description = b"%s", str score_column ="score", trackline=True):
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

    def write_to_narrowPeak (self, fhd, bytes name_prefix = b"peak_", bytes name = b"peak", str score_column="score", trackline=True):
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
        cdef int n_peak
        cdef bytes chrom
        cdef long s
        cdef str peakname

        chrs = list(self.peaks.keys())
        chrs.sort()
        n_peak = 0
        write = fhd.write
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        if trackline:
            write("track type=narrowPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name.decode(), name.decode()))
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1: # from call-summits
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix.decode(), n_peak, subpeak_letters(i))
                        if peak['summit'] == -1:
                            s = -1
                        else:
                            s = peak['summit'] - peak['start']
                        fhd.write( "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n"
                                   %
                                   (chrom.decode(),peak['start'],peak['end'],peakname,int(10*peak[score_column]),
                                    peak['fc'],peak['pscore'],peak['qscore'],s) )
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix.decode(), n_peak)
                    if peak['summit'] == -1:
                        s = -1
                    else:
                        s = peak['summit'] - peak['start']
                    fhd.write( "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n"
                               %
                               (chrom.decode(),peak['start'],peak['end'],peakname,int(10*peak[score_column]),
                                peak['fc'],peak['pscore'],peak['qscore'],s) )
        return

    def write_to_xls (self, ofhd, bytes name_prefix = b"%s_peak_", bytes name = b"MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.


        wait... why I have two write_to_xls in this class?

        """
        write = ofhd.write
        write("\t".join(("chr","start", "end",  "length",  "abs_summit", "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")

        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix

        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for end, group in groupby(peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1:
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix.decode(), n_peak, subpeak_letters(i))
                        #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                        write("%s\t%d\t%d\t%d" % (chrom.decode(),peak['start']+1,peak['end'],peak['length']))
                        write("\t%d" % (peak['summit']+1)) # summit position
                        write("\t%.6g" % (round(peak['pileup'],2))) # pileup height at summit
                        write("\t%.6g" % (peak['pscore'])) # -log10pvalue at summit
                        write("\t%.6g" % (peak['fc'])) # fold change at summit
                        write("\t%.6g" % (peak['qscore'])) # -log10qvalue at summit
                        write("\t%s" % peakname)
                        write("\n")
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix.decode(), n_peak)
                    #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                    write("%s\t%d\t%d\t%d" % (chrom.decode(),peak['start']+1,peak['end'],peak['length']))
                    write("\t%d" % (peak['summit']+1)) # summit position
                    write("\t%.6g" % (round(peak['pileup'],2))) # pileup height at summit
                    write("\t%.6g" % (peak['pscore'])) # -log10pvalue at summit
                    write("\t%.6g" % (peak['fc'])) # fold change at summit
                    write("\t%.6g" % (peak['qscore'])) # -log10qvalue at summit
                    write("\t%s" % peakname)
                    write("\n")
        return


    def overlap_with_other_peaks (self, peaks2, double cover=0):
        """Peaks2 is a PeakIO object or dictionary with can be
        initialized as a PeakIO. check __init__ for PeakIO for detail.

        return how many peaks are intersected by peaks2 by percentage
        coverage on peaks2(if 50%, cover = 0.5).
        """
        cdef int total_num
        cdef list chrs1, chrs2, a
        cdef bytes k

        peaks1 = self.peaks
        if isinstance(peaks2,PeakIO):
            peaks2 = peaks2.peaks
        total_num = 0
        chrs1 = list(peaks1.keys())
        chrs2 = list(peaks2.keys())
        for k in chrs1:
            if not chrs2.count(k):
                continue
            rl1_k = iter(peaks1[k])
            rl2_k = iter(peaks2[k])
            tmp_n = False
            try:
                r1 = rl1_k.next()
                r2 = rl2_k.next()
                while (True):
                    if r2[0] < r1[1] and r1[0] < r2[1]:
                        a = sorted([r1[0],r1[1],r2[0],r2[1]])
                        if float(a[2]-a[1]+1)/r2[2] > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                    if r1[1] < r2[1]:
                        r1 = rl1_k.next()
                        tmp_n = False
                    else:
                        r2 = rl2_k.next()
            except StopIteration:
                continue
        return total_num

    def read_from_xls (self, ofhd):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.

        """
        cdef:
            bytes line = b''
            bytes chrom = b''
            int n_peak = 0
            int start, end, length, summit
            float pileup, pscore, fc, qscore
            list fields
        while True:
            if not (line.startswith('#') or line.strip() == ''): break
            line = ofhd.readline()

        # sanity check
        columns = line.rstrip().split('\t')
        for a,b in zip(columns, ("chr","start", "end",  "length", "abs_summit",
                                 "pileup", "-log10(pvalue)", "fold_enrichment",
                                 "-log10(qvalue)", "name")):
            if not a==b: raise NotImplementedError('column %s not recognized', a)

        add = self.add
        split = str.split
        rstrip = str.rstrip
        for i, line in enumerate(ofhd.readlines()):
            fields = split(line, '\t')
            peak = {}
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

cpdef parse_peakname(peakname):
    """returns peaknumber, subpeak
    """
    cdef:
        bytes peak_id, peaknumber, subpeak
    peak_id = peakname.split(b'_')[-1]
    x = re.split('(\D.*)', peak_id)
    peaknumber = int(x[0])
    try:
        subpeak = x[1]
    except IndexError:
        subpeak = b''
    return (peaknumber, subpeak)

cdef class Region:
    """For plain region of chrom, start and end
    """
    cdef:
        dict regions
        bool __flag_sorted

    def __init__ (self):
        self.regions= {}
        self.__flag_sorted = False

    def add_loc ( self, bytes chrom, int start, int end ):
        if self.regions.has_key(chrom):
            self.regions[chrom].append( (start,end) )
        else:
            self.regions[chrom] = [(start,end), ]
        self.__flag_sorted = False
        return

    def sort (self):
        cdef bytes chrom

        for chrom in list(self.regions.keys()):
            self.regions[chrom].sort()
        self.__flag_sorted = True

    def merge_overlap ( self ):
        cdef bytes chrom
        cdef int s_new_region, e_new_region, i, j

        if not self.__flag_sorted:
            self.sort()
        regions = self.regions
        new_regions = {}
        chrs = list(regions.keys())
        chrs.sort()
        for i in range(len(chrs)):
            chrom = chrs[i]
        #for chrom in chrs:
            new_regions[chrom]=[]
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
                        prev_region = (s_new_region,e_new_region)
                    else:
                        n_append(prev_region)
                        prev_region = regions_chr[i]
            if prev_region:
                n_append(prev_region)
        self.regions = new_regions
        self.sort()
        return True

    def write_to_bed (self, fhd ):
        cdef int i
        cdef bytes chrom

        chrs = list(self.regions.keys())
        chrs.sort()
        for i in range( len(chrs) ):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                fhd.write( "%s\t%d\t%d\n" % (chrom.decode(),region[0],region[1] ) )


cdef class BroadPeakContent:
    cdef:
        long start
        long end
        long length
        float score
        bytes thickStart
        bytes thickEnd
        long blockNum
        bytes  blockSizes
        bytes  blockStarts
        float pileup
        float pscore
        float fc
        float qscore
        bytes name

    def __init__ ( self, long start, long end, float score,
                   bytes thickStart, bytes thickEnd,
                   long blockNum, bytes blockSizes,
                   bytes blockStarts, float pileup,
                   float pscore, float fold_change,
                   float qscore, bytes name = b"NA" ):
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

    def __getitem__ ( self, a ):
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

    def __str__ (self):
        return "start:%d;end:%d;score:%f" % ( self.start, self.end, self.score )


cdef class BroadPeakIO:
    """IO for broad peak information.

    """
    cdef:
        dict peaks

    def __init__ (self):
        self.peaks = {}

    def add (self, char * chromosome, long start, long end, long score = 0,
             bytes thickStart=b".", bytes thickEnd=b".",
             long blockNum=0, bytes blockSizes=b".",
             bytes blockStarts=b".", float pileup = 0,
             float pscore = 0, float fold_change = 0,
             float qscore = 0, bytes name = b"NA" ):
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
        self.peaks[chromosome].append( BroadPeakContent( start, end, score, thickStart, thickEnd,
                                                         blockNum, blockSizes, blockStarts,
                                                         pileup, pscore, fold_change, qscore, name ) )

    def filter_pscore (self, double pscore_cut ):
        cdef:
            bytes chrom
            dict peaks
            dict new_peaks
            list chrs
            BroadPeakContent p

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(list(peaks.keys()))

        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['pscore'] >= pscore_cut]
        self.peaks = new_peaks

    def filter_qscore (self, double qscore_cut ):
        cdef:
            bytes chrom
            dict peaks
            dict new_peaks
            list chrs
            BroadPeakContent p

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(list(peaks.keys()))

        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['qscore'] >= qscore_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)

        """
        cdef:
            bytes chrom
            dict peaks
            dict new_peaks
            list chrs
            BroadPeakContent p

        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low and p['fc']<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        cdef:
            bytes chrom
            dict peaks
            list chrs
            long x

        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x

    def write_to_gappedPeak (self, fhd, bytes name_prefix=b"peak_", bytes name=b'peak', bytes description=b"%s", str score_column="score", trackline=True):
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
        chrs = list(self.peaks.keys())
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        if trackline:
            fhd.write("track name=\"%s\" description=\"%s\" type=gappedPeak nextItemButton=on\n" % (name.decode(), desc.decode()) )
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                if peak["thickStart"] != b".":
                    fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%s\t%s\t0\t%d\t%s\t%s\t%.6g\t%.6g\t%.6g\n"
                               %
                               (chrom.decode(),peak["start"],peak["end"],peakprefix.decode(),n_peak,int(10*peak[score_column]),
                                peak["thickStart"].decode(),peak["thickEnd"].decode(),
                                peak["blockNum"],peak["blockSizes"].decode(),peak["blockStarts"].decode(), peak['fc'], peak['pscore'], peak['qscore'] ) )

    def write_to_Bed12 (self, fhd, bytes name_prefix=b"peak_", bytes name=b'peak', bytes description=b"%s", str score_column="score", trackline=True):
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
        chrs = list(self.peaks.keys())
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        if trackline:
            fhd.write("track name=\"%s\" description=\"%s\" type=bed nextItemButton=on\n" % (name.decode(), desc.decode()) )
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                if peak["thickStart"] == b".":
                    # this will violate gappedPeak format, since it's a complement like broadPeak line.
                    fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\n"
                               %
                               (chrom.decode(),peak["start"],peak["end"],peakprefix.decode(),n_peak,int(10*peak[score_column]) ) )
                else:
                    fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%s\t%s\t0\t%d\t%s\t%s\n"
                               %
                               (chrom.decode(), peak["start"], peak["end"], peakprefix.decode(), n_peak, int(10*peak[score_column]),
                                peak["thickStart"].decode(), peak["thickEnd"].decode(),
                                peak["blockNum"], peak["blockSizes"].decode(), peak["blockStarts"].decode() ))


    def write_to_broadPeak (self, fhd, bytes name_prefix=b"peak_", bytes name=b'peak', bytes description=b"%s", str score_column="score", trackline=True):
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
        cdef int n_peak
        cdef bytes chrom
        cdef long s
        cdef str peakname

        chrs = list(self.peaks.keys())
        chrs.sort()
        n_peak = 0
        write = fhd.write
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        if trackline:
            write("track type=broadPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name.decode(), name.decode()))
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                peak = these_peaks[0]
                peakname = "%s%d" % (peakprefix.decode(), n_peak)
                fhd.write( "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\n" %
                           (chrom.decode(),peak['start'],peak['end'],peakname,int(10*peak[score_column]),
                            peak['fc'],peak['pscore'],peak['qscore'] ) )
        return


    def write_to_xls (self, ofhd, bytes name_prefix=b"%s_peak_", bytes name=b"MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.


        wait... why I have two write_to_xls in this class?

        """
        write = ofhd.write
        write("\t".join(("chr","start", "end",  "length",  "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")

        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix

        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for end, group in groupby(peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                peak = these_peaks[0]
                peakname = "%s%d" % (peakprefix.decode(), n_peak)
                write("%s\t%d\t%d\t%d" % (chrom.decode(),peak['start']+1,peak['end'],peak['length']))
                write("\t%.6g" % (round(peak['pileup'],2))) # pileup height at summit
                write("\t%.6g" % (peak['pscore'])) # -log10pvalue at summit
                write("\t%.6g" % (peak['fc'])) # fold change at summit
                write("\t%.6g" % (peak['qscore'])) # -log10qvalue at summit
                write("\t%s" % peakname)
                write("\n")
        return
