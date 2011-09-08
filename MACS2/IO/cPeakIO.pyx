# Time-stamp: <2011-09-07 18:42:01 Tao Liu>

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

# ------------------------------------
# python modules
# ------------------------------------
from MACS2.Constants import *

# ------------------------------------
# constants
# ------------------------------------
__version__ = "PeakIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "PeakIO class"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
class PeakIO:
    """IO for peak information.

    """

    def __init__ (self):
        self.peaks = {}
    
    def add (self, char * chromosome, long start, long end, long summit = 0, 
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
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                text+= "%s\t%d\t%d\tpeak_%d\t%.2f\n" % (chrom,peak["start"],peak["end"],n_peak,peak["score"])
        return text

    def to_summits_bed (self):
        """Print out peak summits in BED5 format.

        Five columns are chromosome, summit start, summit end, peak name, and peak height.

        """
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                summit_p = peak["summit"]
                text+= "%s\t%d\t%d\tpeak_%d\t%.2f\n" % (chrom,summit_p,summit_p+1,n_peak,peak["score"])
        return text

    def write_to_bed (self, fhd, name_prefix="peak_", score_column="score"):
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
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                fhd.write( "%s\t%d\t%d\t%s%d\t%.2f\n" % (chrom,peak["start"],peak["end"],name_prefix,n_peak,peak[score_column]) )


    def write_to_summit_bed (self, fhd, name_prefix="peak_", score_column="score"):
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
        chrs = self.peaks.keys()
        chrs.sort()
        n_peak = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                summit_p = peak["summit"]
                fhd.write( "%s\t%d\t%d\t%s%d\t%.2f\n" % (chrom,summit_p,summit_p+1,name_prefix,n_peak,peak[score_column]) )

    def write_to_narrowPeak (self, fhd, name_prefix="peak_", score_column="score"):
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
        fhd.write("track type=narrowPeak nextItemButton=on\n")
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                # items in peak: (peak start,peak end, peak length,
                # peak summit, peak height, number of tags in peak
                # region, peak pvalue, peak fold_enrichment, qvalue)
                fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%.2f\t%.2f\t%.2f\t%d\n"
                           %
                           (chrom,peak["start"],peak["end"],name_prefix,n_peak,int(10*peak[score_column]),
                            peak["fc"],peak["pscore"],peak["qscore"],peak["summit"]-peak["start"]) )


###
                
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
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
  
    def write_to_gappedPeak (self, fhd, name_prefix="peak_", name="peak", description="peak description"):
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
        fhd.write("track name=\"%s\" description=\"%s\" type=bed nextItemButton=on\n" % (name, description) )
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n_peak += 1
                fhd.write( "%s\t%d\t%d\t%s%d\t%d\t.\t%d\t%d\t0\t%d\t%s\t%s\n"
                           %
                           (chrom,peak["start"],peak["end"],name_prefix,n_peak,int(peak["score"]),
                            peak["thickStart"],peak["thickEnd"],
                            peak["blockNum"],peak["blockSizes"],peak["blockStarts"] ) )

