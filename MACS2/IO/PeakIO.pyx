# Time-stamp: <2014-08-22 15:41:27 Tao Liu>

"""Module for PeakIO IO classes.

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
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
import re

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
        str name

    def __init__ ( self, int start, int end, int summit, 
                   float peak_score, float pileup, 
                   float pscore, float fold_change, float qscore,
                   str name="NA" ):
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
    
    cpdef add (self, str chromosome, int start, int end, int summit = 0, 
               float peak_score=0, float pileup=0, 
               float pscore=0, float fold_change=0, float qscore=0,
               str name="NA"):
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

    cpdef add_PeakContent ( self, str chromosome, object peakcontent ):
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append(peakcontent)

    def get_data_from_chrom (self, str chrom):
        return self.peaks[chrom]

    def get_chr_names (self):
        return self.peaks.keys()

    def sort ( self ):
        # sort by position
        chrs = sorted(self.peaks.keys())
        for chrom in chrs:
            self.peaks[chrom].sort(key=lambda x:x['start'])
        return
        

    def filter_pscore (self, double pscore_cut ):
        cdef str chrom
        
        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['pscore'] >= pscore_cut]
        self.peaks = new_peaks

    def filter_qscore (self, double qscore_cut ):
        cdef str chrom

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['qscore'] >= qscore_cut]
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
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low and p['fc']<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
  
    # these methods are very fast, specifying types is unnecessary
    def write_to_xls (self, fhd, str name_prefix="%s_peak_", str name="MACS"):
        return self._to_xls(name_prefix=name_prefix, name=name,
                            print_func=fhd.write)

    def _to_xls (self, name_prefix="%s_peak_", name="MACS", print_func=print):

        if self.peaks:
            print_func("\t".join(("chr","start", "end",  "length",  "abs_summit", "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")
        else:
            return
        
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
                        print_func("%s\t%d\t%d\t%d" % (chrom,peak['start']+1,peak['end'],peak['length']))
                        print_func("\t%d" % (peak['summit']+1)) # summit position
                        print_func("\t%.2f" % (round(peak['pileup'],2))) # pileup height at summit
                        print_func("\t%.5f" % (peak['pscore'])) # -log10pvalue at summit
                        print_func("\t%.5f" % (peak['fc'])) # fold change at summit                
                        print_func("\t%.5f" % (peak['qscore'])) # -log10qvalue at summit
                        print_func("\t%s" % peakname)
                        print_func("\n")
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix, n_peak)
                    #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                    print_func("%s\t%d\t%d\t%d" % (chrom,peak['start']+1,peak['end'],peak['length']))
                    print_func("\t%d" % (peak['summit']+1)) # summit position
                    print_func("\t%.2f" % (round(peak['pileup'],2))) # pileup height at summit
                    print_func("\t%.5f" % (peak['pscore'])) # -log10pvalue at summit
                    print_func("\t%.5f" % (peak['fc'])) # fold change at summit                
                    print_func("\t%.5f" % (peak['qscore'])) # -log10qvalue at summit
                    print_func("\t%s" % peakname)
                    print_func("\n")
        return

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
                        print_func("%s\t%d\t%d\t%s%d%s\t%.5f\n" % (chrom,peak['start'],peak['end'],peakprefix,n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    print_func("%s\t%d\t%d\t%s%d\t%.5f\n" % (chrom,peak['start'],peak['end'],peakprefix,n_peak,peak[score_column])) 

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
                        summit_p = peak['summit']
                        print_func("%s\t%d\t%d\t%s%d%s\t%.5f\n" % (chrom,summit_p,summit_p+1,peakprefix,n_peak,subpeak_letters(i),peak[score_column]))
                else:
                    peak = peaks[0]
                    summit_p = peak['summit']
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
        cdef long s

        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        write = fhd.write
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        if trackline:
            write("track type=narrowPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name, name))
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                if len(these_peaks) > 1: # from call-summits
                    for i, peak in enumerate(these_peaks):
                        peakname = "%s%d%s" % (peakprefix, n_peak, subpeak_letters(i))
                        if peak['summit'] == -1:
                            s = -1
                        else:
                            s = peak['summit'] - peak['start']
                        fhd.write( "%s\t%d\t%d\t%s\t%d\t.\t%.5f\t%.5f\t%.5f\t%d\n"
                                   %
                                   (chrom,peak['start'],peak['end'],peakname,int(10*peak[score_column]),
                                    peak['fc'],peak['pscore'],peak['qscore'],s) )
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix, n_peak)
                    if peak['summit'] == -1:
                        s = -1
                    else:
                        s = peak['summit'] - peak['start']
                    fhd.write( "%s\t%d\t%d\t%s\t%d\t.\t%.5f\t%.5f\t%.5f\t%d\n"
                               %
                               (chrom,peak['start'],peak['end'],peakname,int(10*peak[score_column]),
                                peak['fc'],peak['pscore'],peak['qscore'],s) )
        return

    def write_to_xls (self, ofhd, name_prefix="%s_peak_", name="MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.


        wait... why I have two write_to_xls in this class?
        
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
                        write("%s\t%d\t%d\t%d" % (chrom,peak['start']+1,peak['end'],peak['length']))
                        write("\t%d" % (peak['summit']+1)) # summit position
                        write("\t%.2f" % (round(peak['pileup'],2))) # pileup height at summit
                        write("\t%.5f" % (peak['pscore'])) # -log10pvalue at summit
                        write("\t%.5f" % (peak['fc'])) # fold change at summit
                        write("\t%.5f" % (peak['qscore'])) # -log10qvalue at summit
                        write("\t%s" % peakname)
                        write("\n")
                else:
                    peak = these_peaks[0]
                    peakname = "%s%d" % (peakprefix, n_peak)
                    #[start,end,end-start,summit,peak_height,number_tags,pvalue,fold_change,qvalue]
                    write("%s\t%d\t%d\t%d" % (chrom,peak['start']+1,peak['end'],peak['length']))
                    write("\t%d" % (peak['summit']+1)) # summit position
                    write("\t%.2f" % (round(peak['pileup'],2))) # pileup height at summit
                    write("\t%.5f" % (peak['pscore'])) # -log10pvalue at summit
                    write("\t%.5f" % (peak['fc'])) # fold change at summit
                    write("\t%.5f" % (peak['qscore'])) # -log10qvalue at summit                    
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

    def read_from_xls (self, ofhd):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        cdef:
            str line = ''
            str chrom = ''
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
            chrom = fields[0]
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
        str peak_id, peaknumber, subpeak
    peak_id = peakname.split('_')[-1]
    x = re.split('(\D.*)', peak_id)
    peaknumber = int(x[0])
    try:
        subpeak = x[1]
    except IndexError:
        subpeak = ''
    return (peaknumber, subpeak)

cdef class Region:
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


cdef class BroadPeakContent:
    cdef:
        long start
        long end
        long length
        float score
        str thickStart
        str thickEnd
        long blockNum
        str  blockSizes
        str  blockStarts
        float pileup
        float pscore
        float fc
        float qscore
        str name

    def __init__ ( self, long start, long end, float score,
                   str thickStart, str thickEnd,
                   long blockNum, str blockSizes, 
                   str blockStarts, float pileup, 
                   float pscore, float fold_change, 
                   float qscore, str name = "NA" ):
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
             str thickStart=".", str thickEnd=".",
             long blockNum=0, str blockSizes=".", 
             str blockStarts=".", float pileup = 0,
             float pscore = 0, float fold_change = 0,
             float qscore = 0, str name = "NA" ):
        """items
        chromosome : chromosome name,
        start      : broad region start,
        end        : broad region end,
        score      : average score in all blocks,
        thickStart : start of highly enriched region, # could be '.'
        thickEnd   : end of highly enriched region,   # could be '.'
        blockNum   : number of blocks,                # could be 0
        blockSizes : sizes of blocks,                 # could be '.'
        blockStarts: starts of blocks                 # could be '.'
        pileup     : median pileup in region          # could be 0
        pscore     : median pvalue score in region    # could be 0
        fold_change: median fold change in region     # could be 0
        qscore     : median pvalue score in region    # could be 0
        name       : peak name                        # could be 'NA'
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome] = []
        self.peaks[chromosome].append( BroadPeakContent( start, end, score, thickStart, thickEnd, 
                                                         blockNum, blockSizes, blockStarts, 
                                                         pileup, pscore, fold_change, qscore, name ) )

    def filter_pscore (self, double pscore_cut ):
        cdef str chrom
        
        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['pscore'] >= pscore_cut]
        self.peaks = new_peaks

    def filter_qscore (self, double qscore_cut ):
        cdef str chrom

        peaks = self.peaks
        new_peaks = {}
        chrs = sorted(peaks.keys())
        
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p['qscore'] >= qscore_cut]
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
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low and p['fc']<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p['fc'] >= fc_low]
        self.peaks = new_peaks

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
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        try: desc = description % name
        except: desc = description
        if trackline:
            fhd.write("track name=\"%s\" description=\"%s\" type=gappedPeak nextItemButton=on\n" % (name, desc) )
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                if peak["thickStart"] != ".":
                    fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%s\t%s\t0\t%d\t%s\t%s\t%.5f\t%.5f\t%.5f\n"
                               %
                               (chrom,peak["start"],peak["end"],peakprefix,n_peak,int(10*peak["qscore"]),
                                peak["thickStart"],peak["thickEnd"],
                                peak["blockNum"],peak["blockSizes"],peak["blockStarts"], peak['fc'], peak['pscore'], peak['qscore'] ) )

    def write_to_Bed12 (self, fhd, name_prefix="peak_", name='peak', description="%s", trackline=True):
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
                if peak["thickStart"] == ".":
                    # this will violate gappedPeak format, since it's a complement like broadPeak line.
                    fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\n"
                               %
                               (chrom,peak["start"],peak["end"],peakprefix,n_peak,int(10*peak["qscore"]) ) )
                else:
                    fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%s\t%s\t0\t%d\t%s\t%s\n"
                               %
                               (chrom, peak["start"], peak["end"], peakprefix, n_peak, int(10*peak["qscore"]),
                                peak["thickStart"], peak["thickEnd"],
                                peak["blockNum"], peak["blockSizes"], peak["blockStarts"] ))


    def write_to_broadPeak (self, fhd, name_prefix="peak_", name='peak', description="%s", trackline=True):
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
        
        """
        cdef int n_peak
        cdef str chrom
        cdef long s

        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        write = fhd.write
        try: peakprefix = name_prefix % name
        except: peakprefix = name_prefix
        if trackline:
            write("track type=broadPeak name=\"%s\" description=\"%s\" nextItemButton=on\n" % (name, name))
        for chrom in chrs:
            for end, group in groupby(self.peaks[chrom], key=itemgetter("end")):
                n_peak += 1
                these_peaks = list(group)
                peak = these_peaks[0]
                peakname = "%s%d" % (peakprefix, n_peak)
                fhd.write( "%s\t%d\t%d\t%s\t%d\t.\t%.5f\t%.5f\t%.5f\n" %
                           (chrom,peak['start'],peak['end'],peakname,int(10*peak["qscore"]),
                            peak['fc'],peak['pscore'],peak['qscore'] ) )
        return


    def write_to_xls (self, ofhd, name_prefix="%s_peak_", name="MACS"):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.


        wait... why I have two write_to_xls in this class?
        
        """
        write = ofhd.write
        write("\t".join(("chr","start", "end",  "length",  "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name"))+"\n")
        
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
                peak = these_peaks[0]
                peakname = "%s%d" % (peakprefix, n_peak)
                write("%s\t%d\t%d\t%d" % (chrom,peak['start']+1,peak['end'],peak['length']))
                write("\t%.2f" % (round(peak['pileup'],2))) # pileup height at summit
                write("\t%.5f" % (peak['pscore'])) # -log10pvalue at summit
                write("\t%.5f" % (peak['fc'])) # fold change at summit
                write("\t%.5f" % (peak['qscore'])) # -log10qvalue at summit                    
                write("\t%s" % peakname)
                write("\n")
        return
