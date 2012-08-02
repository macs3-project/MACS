# cython: profile=True
# Time-stamp: <2012-08-01 18:05:02 Tao Liu>

"""Module for PeakIO IO classes.

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
##############################################################################
# ** NOTE ** THIS MODULE USES python v3-style print() not v2 print keyword   #
##############################################################################
from __future__ import print_function # this line must be first
# ------------------------------------
# python modules
# ------------------------------------
from MACS2.Constants import *
from itertools import groupby
from operator import itemgetter

# ------------------------------------
# constants
# ------------------------------------
__version__ = "PeakIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "PeakIO class"

# ------------------------------------
# Misc functions
# ------------------------------------
cdef subpeak_letters( int i):
    if i < 26:
        return chr(97+i)
    else:
        return subpeak_letters(i / 26) + chr(97 + (i % 26))
    
# ------------------------------------
# Classes
# ------------------------------------
class PeakIO:
    """IO for peak information.

    """

    def __init__ (self):
        self.peaks = {}
    
    def add (self, str chromosome, long start, long end, long summit = 0, 
             double peak_score=0, int pileup=0, 
             double pscore=0, double fold_change=0, double qscore=0):
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
        self.peaks[chromosome].append({"start":start,
                                       "end":end,
                                       "length":end-start,
                                       "summit":summit,
                                       "score":peak_score,
                                       "pileup":pileup,
                                       "pscore":pscore,
                                       "fc":fold_change,
                                       "qscore":qscore})

    def filter_pscore (self, double pscore_cut ):
        cdef str chrom
        
        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p["pscore"] >= pscore_cut]
        self.peaks = new_peaks

    def filter_qscore (self, double qscore_cut ):
        cdef str chrom

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p["qscore"] >= qscore_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p["fc"] >= fc_low and p["fc"]<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p["fc"] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
  
    def _to_bed(self, name_prefix="%s_peak_", name="MACS",
                description="%s", score_column="score",
                print_func=print, trackline=False):
        """
        generalization of tobed and write_to_bed
        """
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        trackcontents = (name.replace("\"", "\\\""), desc.replace("\"", "\\\""))
        if trackline:
            try: print_func('track name="%s (peaks)" description="%s" visibility=1\n' % trackcontents)
            except: print_func('track name=MACS description=Unknown') 
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        print_func("%s\t%d\t%d\t%s%d%s\t%.5f\n" % (chrom,peak["start"],peak["end"],peakprefix,n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    print_func("%s\t%d\t%d\t%s%d\t%.5f\n" % (chrom,peak["start"],peak["end"],peakprefix,n_peak,peak[score_column])) 

    def _to_summits_bed(self, name_prefix="%s_peak_", name="MACS",
                        description = "%s", score_column="score",
                        print_func=print, trackline=False):
        """ 
        generalization of to_summits_bed and write_to_summit_bed
        """
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        trackcontents = (name.replace("\"", "\\\""), desc.replace("\"", "\\\""))
        if trackline:
            try: print_func('track name="%s (summits)" description="%s" visibility=1\n' % trackcontents)
            except: print_func('track name=MACS description=Unknown') 
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        summit_p = peak["summit"]
                        print_func("%s\t%d\t%d\t%s%d%s\t%.5f\n" % (chrom,summit_p,summit_p+1,peakprefix,n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    summit_p = peak["summit"]
                    print_func("%s\t%d\t%d\t%s%d\t%.5f\n" % (chrom,summit_p,summit_p+1,peakprefix,n_peak,peak[score_column]))

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
        return self._to_bed(name_prefix="peak_", score_column="score")

    def to_summits_bed (self):
        """Print out peak summits in BED5 format.

        Five columns are chromosome, summit start, summit end, peak name, and peak height.

        """
        return self._to_summits_bed(name_prefix="peak_", score_column="score")

    # these methods are very fast, specifying types is unnecessary
    def write_to_bed (self, fhd, str name_prefix="peak_", str name="MACS",
                        str description = "%s", str score_column="score", trackline=True):
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
        return self._to_bed(name_prefix=name_prefix, name=name,
                            description=description, score_column=score_column,
                            print_func=fhd.write, trackline=trackline)

    def write_to_summit_bed (self, fhd, name_prefix="peak_", name="MACS",
                             description = "%s", score_column="score", trackline=True):
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

    def write_to_narrowPeak (self, fhd, name_prefix="peak_", name="peak", score_column="score", trackline=True):
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
        |in MACS2 * |      |'0', the DCC will assign this based on  |
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
        cdef str chrom

        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        if trackline:
            fhd.write("track type=narrowPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name, name))
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                # items in peak: (peak start,peak end, peak length,
                # peak summit, peak height, number of tags in peak
                # region, peak pvalue, peak fold_enrichment, qvalue)
                fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%.5f\t%.5f\t%.5f\t%d\n"
                           %
                           (chrom,peak["start"],peak["end"],peakprefix,n_peak,int(10*peak[score_column]),
                            peak["fc"],peak["pscore"],peak["qscore"],peak["summit"]-peak["start"]) )

    def write_to_xls (self, ofhd, name_prefix="%s_peak_", name="MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        write = ofhd.write
        write("\t".join(("chr","start", "end",  "length",  "abs_summit", "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")
        
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix

        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for end, group in groupby(peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1:
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix, n_peak, subpeak_letters(i))
                        #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                        write("%s\t%d\t%d\t%d" % (chrom,peak["start"]+1,peak["end"],peak["length"]))
                        write("\t%d" % (peak["summit"]+1)) # summit position
                        write("\t%.5f" % (peak["pileup"])) # pileup height at summit
                        write("\t%.5f" % (peak["pscore"])) # -log10pvalue at summit
                        write("\t%.5f" % (peak["fc"])) # fold change at summit
                        write("\t%.5f" % (peak["qscore"])) # -log10qvalue at summit
                        write("\t%s" % peakname)
                        write("\n")
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix, n_peak)
                    #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                    write("%s\t%d\t%d\t%d" % (chrom,peak["start"]+1,peak["end"],peak["length"]))
                    write("\t%d" % (peak["summit"]+1)) # summit position
                    write("\t%.5f" % (peak["pileup"])) # pileup height at summit
                    write("\t%.5f" % (peak["pscore"])) # -log10pvalue at summit
                    write("\t%.5f" % (peak["fc"])) # fold change at summit
                    write("\t%.5f" % (peak["qscore"])) # -log10qvalue at summit                    
                    write("\t%s" % peakname)
                    write("\n")
        return


    def overlap_with_other_peaks (self, peaks2, double cover=0):
        """Peaks2 is a PeakIO object or dictionary with can be
        initialzed as a PeakIO. check __init__ for PeakIO for detail.

        return how many peaks are intersected by peaks2 by percentage
        coverage on peaks2(if 50%, cover = 0.5).
        """
        cdef int total_num
        cdef list chrs1, chrs2, a
        cdef str k
        
        peaks1 = self.peaks
        if isinstance(peaks2,PeakIO):
            peaks2 = peaks2.peaks
        total_num = 0
        chrs1 = peaks1.keys()
        chrs2 = peaks2.keys()
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

class Region:
    """For plain region of chrom, start and end
    """
    def __init__ (self):
        self.regions= {}
        self.__flag_sorted = False

    def add_loc ( self, str chrom, int start, int end ):
        if self.regions.has_key(chrom):
            self.regions[chrom].append( (start,end) )
        else:
            self.regions[chrom] = [(start,end), ]
        self.__flag_sorted = False            
        return

    def sort (self):
        cdef str chrom

        for chrom in self.regions.keys():
            self.regions[chrom].sort()
        self.__flag_sorted = True
    
    def merge_overlap ( self ):
        cdef str chrom
        cdef int s_new_region, e_new_region, i, j
        
        if not self.__flag_sorted:
            self.sort()
        regions = self.regions
        new_regions = {}
        chrs = regions.keys()
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
        cdef str chrom
        
        chrs = self.regions.keys()
        chrs.sort()
        for i in range( len(chrs) ):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                fhd.write( "%s\t%d\t%d\n" % (chrom,region[0],region[1] ) )


###
class DiffPeakIO:
    """IO for differential peak information.

    """

    def __init__ (self):
        self.peaks = {}
    
    def add (self, str chromosome, long start, long end, long summit = 0, 
             double diff_score=0, int pileup=0, 
             double pscore=0, double fold_change=0, double qscore=0):
        """items:
        start:start
        end:end,
        length:end-start,
        summit:summit,                  # summit is where the highest pileup is in a differetnial region, or the common region for only condition A.
        score:diff_score,               # diff score is the maximum of qscore in diff/common region, for A vs B and B vs A. If B vs A is bigger, put minus sign later.
        pileup:pileup,                  # the highest pileup the summit.
        pscore:pscore,                  # pscore is the maximum of pscore in diff/common region, for A vs B and B vs A. If B vs A is bigger, put minus sign later.
        fc:fold_change,                 # fc is the maximum of foldchange in diff/common region, for A vs B and B vs A. If B vs A is bigger, put minus sign later.
        qscore:qscore                   # qscore is the maximum of qscore in diff/common region, for A vs B and B vs A. If B vs A is bigger, put minus sign later.
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append({"start":start,
                                       "end":end,
                                       "length":end-start,
                                       "summit":summit,
                                       "score":diff_score,
                                       "pileup":pileup,
                                       "pscore":pscore,
                                       "fc":fold_change,
                                       "qscore":qscore})

    def filter_pscore (self, double pscore_cut ):
        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p["pscore"] >= pscore_cut]
        self.peaks = new_peaks

    def filter_qscore (self, double qscore_cut ):
        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p["qscore"] >= qscore_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p["fc"] >= fc_low and p["fc"]<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p["fc"] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
  
    def _to_bed(self, name_prefix="%s_peak_", name="MACS",
                description="%s", score_column="score",
                print_func=print, trackline=False):
        """
        generalization of tobed and write_to_bed
        """
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        trackcontents = (name.replace("\"", "\\\""), desc.replace("\"", "\\\""))
        if trackline:
            try: print_func('track name="%s (peaks)" description="%s" visibility=1\n' % trackcontents)
            except: print_func('track name=MACS description=Unknown') 
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        print_func("%s\t%d\t%d\t%s%d%s\t%.5f\n" % (chrom,peak["start"],peak["end"],peakprefix,n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    print_func("%s\t%d\t%d\t%s%d\t%.5f\n" % (chrom,peak["start"],peak["end"],peakprefix,n_peak,peak[score_column])) 

    def _to_summits_bed(self, name_prefix="%s_peak_", name="MACS",
                        description = "%s", score_column="score",
                        print_func=print, trackline=False):
        """ 
        generalization of to_summits_bed and write_to_summit_bed
        """
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        trackcontents = (name.replace("\"", "\\\""), desc.replace("\"", "\\\""))
        if trackline:
            try: print_func('track name="%s (summits)" description="%s" visibility=1\n' % trackcontents)
            except: print_func('track name=MACS description=Unknown') 
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                peaks = list(group)
                if len(peaks) > 1:
                    for i, peak in enumerate(peaks):
                        summit_p = peak["summit"]
                        print_func("%s\t%d\t%d\t%s%d%s\t%.5f\n" % (chrom,summit_p,summit_p+1,name_prefix,n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    summit_p = peak["summit"]
                    print_func("%s\t%d\t%d\t%s%d\t%.5f\n" % (chrom,summit_p,summit_p+1,name_prefix,n_peak,peak[score_column]))
  
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
        return self._to_bed(name_prefix="peak_", score_column="score")

    def to_summits_bed (self):
        """Print out peak summits in BED5 format.

        Five columns are chromosome, summit start, summit end, peak name, and peak height.

        """
        return self._to_summits_bed(name_prefix="peak_", score_column="score")

    def write_to_bed (self, fhd, str name_prefix="peak_", str name="MACS",
                      str description = "%s", str score_column="score", trackline=True):
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
        return self._to_bed(name_prefix=name_prefix, score_column=score_column,
                            print_func=fhd.write)

    def write_to_summit_bed (self, fhd, name_prefix="peak_", name="MACS",
                             description = "%s", score_column="score", trackline=True):
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
        return self._to_summits_bed(name_prefix=name_prefix, score_column=score_column,
                            print_func=fhd.write)

    def write_to_narrowPeak (self, fhd, name_prefix="peak_", name="peak", score_column="score", trackline=True):
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
        |in MACS2 * |      |'0', the DCC will assign this based on  |
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
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        if trackline:
            fhd.write("track type=narrowPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name, name))
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                # items in peak: (peak start,peak end, peak length,
                # peak summit, peak height, number of tags in peak
                # region, peak pvalue, peak fold_enrichment, qvalue)
                fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%.5f\t%.5f\t%.5f\t%d\n"
                           %
                           (chrom,peak["start"],peak["end"],peakprefix,n_peak,int(10*peak[score_column]),
                            peak["fc"],peak["pscore"],peak["qscore"],peak["summit"]-peak["start"]) )


#
                
class BroadPeakIO:
    """IO for broad peak information.

    """

    def __init__ (self):
        self.peaks = {}
    
    def add (self, char * chromosome, long start, long end, long score = 0,
             long thickStart=0, long thickEnd=0,
             long blockNum=0, char *blockSizes="", 
             char * blockStarts="" ):
        """items
        chromosome : chromosome name,
        start      : broad region start,
        end        : broad region end,
        score      : average score in all blocks,
        thickStart : start of highly enriched region,
        thickEnd   : end of highly enriched region,
        blockNum   : number of blocks,
        blockSizes : sizes of blocks,
        blockStarts: starts of blocks
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append({"start":start,
                                       "end":end,
                                       "score":score,
                                       "thickStart":thickStart,
                                       "thickEnd":thickEnd,
                                       "blockNum":blockNum,
                                       "blockSizes":blockSizes,
                                       "blockStarts":blockStarts,
                                       }
                                      )

    def total (self):
        cdef str chrom
        cdef long x
        
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
  
    def write_to_gappedPeak (self, fhd, name_prefix="peak_", name='peak', description="%s", trackline=True):
        """Print out peaks in bed12 format.

        This format is basically a BED12 format.

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
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        if trackline:
            fhd.write("track name=\"%s\" description=\"%s\" type=bed nextItemButton=on\n" % (name, desc) )
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%d\t%d\t0\t%d\t%s\t%s\n"
                           %
                           (chrom,peak["start"],peak["end"],peakprefix,n_peak,int(peak["score"]),
                            peak["thickStart"],peak["thickEnd"],
                            peak["blockNum"],peak["blockSizes"],peak["blockStarts"] ) )

