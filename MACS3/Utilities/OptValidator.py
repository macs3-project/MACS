# Time-stamp: <2022-06-10 10:24:34 Tao Liu>

"""Module Description

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import os
import re
import resource # we turn on memory monitoring inside of logging
from argparse import ArgumentError
from subprocess import Popen, PIPE
from math import log

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.IO.Parser import BEDParser, ELANDResultParser, ELANDMultiParser, \
    ELANDExportParser, SAMParser, BAMParser, BAMPEParser,\
    BEDPEParser, BowtieParser,  guess_parser

from MACS3.Utilities.Constants import EFFECTIVEGS as efgsize
# ------------------------------------
# constants
# ------------------------------------

import logging
import MACS3.Utilities.Logger
logger = logging.getLogger(__name__)

# ------------------------------------
# Misc functions
# ------------------------------------

logger = logging.getLogger(__name__)

def opt_validate_callpeak ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logger.error("Error when interpreting --gsize option: %s" % options.gsize)
            logger.error("Available shortcuts of effective genome sizes are %s" % ",".join(list(efgsize.keys())))
            sys.exit(1)

    # format
    options.gzip_flag = False           # if the input is gzip file

    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
        options.nomodel = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
        options.nomodel = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)

    # duplicate reads
    if options.keepduplicates != "auto" and options.keepduplicates != "all":
        if not options.keepduplicates.isdigit():
            logger.error("--keep-dup should be 'auto', 'all' or an integer!")
            sys.exit(1)

    # shiftsize>0
    #if options.shiftsize:               # only if --shiftsize is set, it's true
    #    options.extsize = 2 * options.shiftsize
    #else:                               # if --shiftsize is not set
    #    options.shiftsize = options.extsize / 2
    if options.extsize < 1 :
        logger.error("--extsize must >= 1!")
        sys.exit(1)

    # refine_peaks, call_summits can't be combined with --broad
    #if options.broad and (options.refine_peaks or options.call_summits):
    #    logger.error("--broad can't be combined with --refine-peaks or --call-summits!")
    #    sys.exit(1)

    if options.broad and options.call_summits:
        logger.error("--broad can't be combined with --call-summits!")
        sys.exit(1)

    if options.pvalue:
        # if set, ignore qvalue cutoff
        options.log_qvalue = None
        options.log_pvalue = log(options.pvalue,10)*-1
    else:
        options.log_qvalue = log(options.qvalue,10)*-1
        options.log_pvalue = None
    if options.broad:
        options.log_broadcutoff = log(options.broadcutoff,10)*-1

    # uppercase the format string
    options.format = options.format.upper()

    # d_min is non-negative
    if options.d_min < 0:
        logger.error("Minimum fragment size shouldn't be negative!" % options.d_min)
        sys.exit(1)

    # upper and lower mfold
    options.lmfold = options.mfold[0]
    options.umfold = options.mfold[1]
    if options.lmfold > options.umfold:
        logger.error("Upper limit of mfold should be greater than lower limit!" % options.mfold)
        sys.exit(1)

    # output filenames
    options.peakxls = os.path.join( options.outdir, options.name+"_peaks.xls" )
    options.peakbed = os.path.join( options.outdir, options.name+"_peaks.bed" )
    options.peakNarrowPeak = os.path.join( options.outdir, options.name+"_peaks.narrowPeak" )
    options.peakBroadPeak = os.path.join( options.outdir, options.name+"_peaks.broadPeak" )
    options.peakGappedPeak = os.path.join( options.outdir, options.name+"_peaks.gappedPeak" )
    options.summitbed = os.path.join( options.outdir, options.name+"_summits.bed" )
    options.bdg_treat = os.path.join( options.outdir, options.name+"_treat_pileup.bdg" )
    options.bdg_control= os.path.join( options.outdir, options.name+"_control_lambda.bdg" )
    if options.cutoff_analysis:
        options.cutoff_analysis_file = os.path.join( options.outdir, options.name+"_cutoff_analysis.txt" )
    else:
        options.cutoff_analysis_file = "None"
    #options.negxls  = os.path.join( options.name+"_negative_peaks.xls" )
    #options.diagxls = os.path.join( options.name+"_diag.xls" )
    options.modelR  = os.path.join( options.outdir, options.name+"_model.r" )
    #options.pqtable  = os.path.join( options.outdir, options.name+"_pq_table.txt" )

    options.argtxt = "\n".join((
        "# Command line: %s" % " ".join(sys.argv[1:]),\
        "# ARGUMENTS LIST:",\
        "# name = %s" % (options.name),\
        "# format = %s" % (options.format),\
        "# ChIP-seq file = %s" % (options.tfile),\
        "# control file = %s" % (options.cfile),\
        "# effective genome size = %.2e" % (options.gsize),\
        #"# tag size = %d" % (options.tsize),\
        "# band width = %d" % (options.bw),\
        "# model fold = %s\n" % (options.mfold),\
        ))

    if options.pvalue:
        if options.broad:
            options.argtxt +=  "# pvalue cutoff for narrow/strong regions = %.2e\n" % (options.pvalue)
            options.argtxt +=  "# pvalue cutoff for broad/weak regions = %.2e\n" % (options.broadcutoff)
            options.argtxt +=  "# qvalue will not be calculated and reported as -1 in the final output.\n"
        else:
            options.argtxt +=  "# pvalue cutoff = %.2e\n" % (options.pvalue)
            options.argtxt +=  "# qvalue will not be calculated and reported as -1 in the final output.\n"
    else:
        if options.broad:
            options.argtxt +=  "# qvalue cutoff for narrow/strong regions = %.2e\n" % (options.qvalue)
            options.argtxt +=  "# qvalue cutoff for broad/weak regions = %.2e\n" % (options.broadcutoff)
        else:
            options.argtxt +=  "# qvalue cutoff = %.2e\n" % (options.qvalue)

    if options.maxgap:
        options.argtxt += "# The maximum gap between significant sites = %d\n" % options.maxgap
    else:
        options.argtxt += "# The maximum gap between significant sites is assigned as the read length/tag size.\n"
    if options.minlen:
        options.argtxt += "# The minimum length of peaks = %d\n" % options.minlen
    else:
        options.argtxt += "# The minimum length of peaks is assigned as the predicted fragment length \"d\".\n"

    if options.downsample:
        options.argtxt += "# Larger dataset will be randomly sampled towards smaller dataset.\n"
        if options.seed >= 0:
            options.argtxt += "# Random seed has been set as: %d\n" % options.seed
    else:
        if options.scaleto == "large":
            options.argtxt += "# Smaller dataset will be scaled towards larger dataset.\n"
        else:
            options.argtxt += "# Larger dataset will be scaled towards smaller dataset.\n"

    if options.ratio != 1.0:
        options.argtxt += "# Using a custom scaling factor: %.2e\n" % (options.ratio)

    if options.cfile:
        options.argtxt += "# Range for calculating regional lambda is: %d bps and %d bps\n" % (options.smalllocal,options.largelocal)
    else:
        options.argtxt += "# Range for calculating regional lambda is: %d bps\n" % (options.largelocal)

    if options.broad:
        options.argtxt += "# Broad region calling is on\n"
    else:
        options.argtxt += "# Broad region calling is off\n"

    if options.fecutoff != 1.0:
        options.argtxt += "# Additional cutoff on fold-enrichment is: %.2f\n" % (options.fecutoff)

    if options.format in ["BAMPE", "BEDPE"]:
        # neutralize SHIFT
        options.shift = 0
        options.argtxt += "# Paired-End mode is on\n"
    else:
        options.argtxt += "# Paired-End mode is off\n"

    #if options.refine_peaks:
    #    options.argtxt += "# Refining peak for read balance is on\n"
    if options.call_summits:
        options.argtxt += "# Searching for subpeak summits is on\n"

    if options.do_SPMR and options.store_bdg:
        options.argtxt += "# MACS will save fragment pileup signal per million reads\n"

    return options

def opt_validate_diffpeak ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # format
    options.gzip_flag = False           # if the input is gzip file

#    options.format = options.format.upper()
    # fox this stuff
#    if True: pass
#    elif options.format == "AUTO":
#        options.parser = guess_parser
#    else:
#        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
#        sys.exit(1)

    if options.peaks_pvalue:
        # if set, ignore qvalue cutoff
        options.peaks_log_qvalue = None
        options.peaks_log_pvalue = log(options.peaks_pvalue,10)*-1
        options.track_score_method = 'p'
    else:
        options.peaks_log_qvalue = log(options.peaks_qvalue,10)*-1
        options.peaks_log_pvalue = None
        options.track_score_method = 'q'

    if options.diff_pvalue:
        # if set, ignore qvalue cutoff
        options.log_qvalue = None
        options.log_pvalue = log(options.diff_pvalue,10)*-1
        options.score_method = 'p'
    else:
        options.log_qvalue = log(options.diff_qvalue,10)*-1
        options.log_pvalue = None
        options.score_method = 'q'

    # output filenames
    options.peakxls = options.name+"_diffpeaks.xls"
    options.peakbed = options.name+"_diffpeaks.bed"
    options.peak1xls = options.name+"_diffpeaks_by_peaks1.xls"
    options.peak2xls = options.name+"_diffpeaks_by_peaks2.xls"
    options.bdglogLR = options.name+"_logLR.bdg"
    options.bdgpvalue = options.name+"_logLR.bdg"
    options.bdglogFC = options.name+"_logLR.bdg"

    options.call_peaks = True
    if not (options.peaks1 == '' or options.peaks2 == ''):
        if options.peaks1 == '':
            raise ArgumentError('peaks1', 'Must specify both peaks1 and peaks2, or neither (to call peaks again)')
        elif options.peaks2 == '':
            raise ArgumentError('peaks2', 'Must specify both peaks1 and peaks2, or neither (to call peaks again)')
        options.call_peaks = False
        options.argtxt = "\n".join((
            "# ARGUMENTS LIST:",\
            "# name = %s" % (options.name),\
#            "# format = %s" % (options.format),\
            "# ChIP-seq file 1 = %s" % (options.t1bdg),\
            "# control file 1 = %s" % (options.c1bdg),\
            "# ChIP-seq file 2 = %s" % (options.t2bdg),\
            "# control file 2 = %s" % (options.c2bdg),\
            "# Peaks, condition 1 = %s" % (options.peaks1),\
            "# Peaks, condition 2 = %s" % (options.peaks2),\
            ""
            ))
    else:
        options.argtxt = "\n".join((
            "# ARGUMENTS LIST:",\
            "# name = %s" % (options.name),\
#            "# format = %s" % (options.format),\
            "# ChIP-seq file 1 = %s" % (options.t1bdg),\
            "# control file 1 = %s" % (options.c1bdg),\
            "# ChIP-seq file 2 = %s" % (options.t2bdg),\
            "# control file 2 = %s" % (options.c2bdg),\
            ""
            ))

        if options.peaks_pvalue:
            options.argtxt +=  "# treat/control -log10(pvalue) cutoff = %.2e\n" % (options.peaks_log_pvalue)
            options.argtxt +=  "# treat/control -log10(qvalue) will not be calculated and reported as -1 in the final output.\n"
        else:
            options.argtxt +=  "# treat/control -log10(qvalue) cutoff = %.2e\n" % (options.peaks_log_qvalue)

    if options.diff_pvalue:
        options.argtxt +=  "# differential pvalue cutoff = %.2e\n" % (options.log_pvalue)
        options.argtxt +=  "# differential qvalue will not be calculated and reported as -1 in the final output.\n"
    else:
        options.argtxt +=  "# differential qvalue cutoff = %.2e\n" % (options.log_qvalue)

    return options

def opt_validate_filterdup ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logger.error("Error when interpreting --gsize option: %s" % options.gsize)
            logger.error("Available shortcuts of effective genome sizes are %s" % ",".join(list(efgsize.keys())))
            sys.exit(1)

    # format

    options.gzip_flag = False           # if the input is gzip file

    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)

    # duplicate reads
    if options.keepduplicates != "auto" and options.keepduplicates != "all":
        if not options.keepduplicates.isdigit():
            logger.error("--keep-dup should be 'auto', 'all' or an integer!")
            sys.exit(1)

    # uppercase the format string
    options.format = options.format.upper()

    return options

def opt_validate_randsample ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # format

    options.gzip_flag = False           # if the input is gzip file

    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)

    # uppercase the format string
    options.format = options.format.upper()

    # percentage or number
    if options.percentage:
        if options.percentage > 100.0:
            logger.error("Percentage can't be bigger than 100.0. Please check your options and retry!")
            sys.exit(1)
    elif options.number:
        if options.number <= 0:
            logger.error("Number of tags can't be smaller than or equal to 0. Please check your options and retry!")
            sys.exit(1)

    return options

def opt_validate_refinepeak ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # format

    options.gzip_flag = False           # if the input is gzip file

    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)

    # uppercase the format string
    options.format = options.format.upper()

    return options

def opt_validate_predictd ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logger.error("Error when interpreting --gsize option: %s" % options.gsize)
            logger.error("Available shortcuts of effective genome sizes are %s" % ",".join(list(efgsize.keys())))
            sys.exit(1)

    # format
    options.gzip_flag = False           # if the input is gzip file

    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
        options.nomodel = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
        options.nomodel = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)

    # uppercase the format string
    options.format = options.format.upper()

    # d_min is non-negative
    if options.d_min < 0:
        logger.error("Minimum fragment size shouldn't be negative!" % options.d_min)
        sys.exit(1)

    # upper and lower mfold
    options.lmfold = options.mfold[0]
    options.umfold = options.mfold[1]
    if options.lmfold > options.umfold:
        logger.error("Upper limit of mfold should be greater than lower limit!" % options.mfold)
        sys.exit(1)

    options.modelR  = os.path.join( options.outdir, options.rfile )

    return options


def opt_validate_pileup ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
    
    # format

    options.gzip_flag = False           # if the input is gzip file

    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    else:
        logger.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)

    # uppercase the format string
    options.format = options.format.upper()

    # extsize
    if options.extsize <= 0 :
        logger.error("--extsize must > 0!")
        sys.exit(1)

    return options

def opt_validate_bdgcmp ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info

    # methods should be valid:

    for method in set(options.method):
        if method not in [ 'ppois', 'qpois', 'subtract', 'logFE', 'FE', 'logLR', 'slogLR', 'max' ]:
            logger.error( "Invalid method: %s" % method )
            sys.exit( 1 )

    # # of --ofile must == # of -m

    if options.ofile:
        if len(options.method) != len(options.ofile):
            logger.error("The number and the order of arguments for --ofile must be the same as for -m.")
            sys.exit(1)

    return options


def opt_validate_cmbreps ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info

    # methods should be valid:


    if options.method not in [ 'fisher', 'max', 'mean']:
        logger.error( "Invalid method: %s" % options.method )
        sys.exit( 1 )

    if len( options.ifile ) < 2:
        logger.error("Combining replicates needs at least two replicates!")
        sys.exit( 1 )

    # # of -i must == # of -w

    # if not options.weights:
    #     options.weights = [ 1.0 ] * len( options.ifile )

    # if len( options.ifile ) != len( options.weights ):
    #     logger.error("Must provide same number of weights as number of input files.")
    #     sys.exit( 1 )

    # if options.method == "fisher" and len( options.ifile ) > 3:
    #     logger.error("NOT IMPLEMENTED! Can't combine more than 3 replicates using Fisher's method.")
    #     sys.exit( 1 )

    return options


def opt_validate_bdgopt ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info

    # methods should be valid:

    if options.method.lower() not in [ 'multiply', 'add', 'p2q', 'max', 'min']:
        logger.error( "Invalid method: %s" % options.method )
        sys.exit( 1 )

    if options.method.lower() in [ 'multiply', 'add' ] and not options.extraparam:
        logger.error( "Need EXTRAPARAM for method multiply or add!")
        sys.exit( 1 )

    return options

def opt_validate_callvar ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info

    # methods should be valid:

    if options.np <= 0:
        options.np = 1
    return options


def opt_validate_hmmratac ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logger.setLevel((4-options.verbose)*10)

    options.error   = logger.critical        # function alias
    options.warn    = logger.warning
    options.debug   = logger.debug
    options.info    = logger.info
   
    # input options.argtxt for hmmratac
    options.argtxt = "# Command line: %s\n" % " ".join(sys.argv[1:])
    #        "# ARGUMENTS LIST:",\
    #        "# outfile = %s" % (options.ofile),\
    #        "# input file = %s\n" % (options.bam_file),\
    # ... add additional

    # Output options
    #if options.store_bdg:
    #    options.argtxt += "# HMMRATAC will report whole genome bedgraph of all state annotations. \n"
    
    #if options.store_bgscore:
    #    options.argtxt += "# HMMRATAC score will be added to each state annotation in bedgraph. \n"

    #if options.store_peaks:
    #    options.argtxt += "# Peaks not reported in bed format\n"

    #if options.print_exclude:
    #    options.print_exclude = os.path.join(options.outdir, options.ofile+"Output_exclude.bed")
    #else:
    #    options.print_exclude = "None"
    
    #if options.print_train:
    #    options.print_train = os.path.join(options.outdir, options.ofile+"Output_training.bed")
    #else:
    #    options.print_train = "None"


    # EM
    # em_skip
    if options.em_skip:
        options.argtxt += "# EM training not performed on fragment distribution. \n"
    # em_means non-negative
    if sum( [ x < 0 for x in options.em_means ] ):
        logger.error(" --means should not be negative! ")
        sys.exit( 1 )
    # em_stddev non-negative
    if sum( [ x < 0 for x in options.em_stddevs ] ):
        logger.error(" --stddev should not be negative! ")
        sys.exit( 1 )


    # HMM
    # hmm_states non-negative int, warn if not k=3
    #if options.hmm_states <=0:
    #    logger.error(" -s, --states must be an integer >= 0.")
    #    sys.exit( 1 )
    #elif options.hmm_states != 3 and options.hmm_states > 0 and options.store_peaks == False:
    #    logger.warn(" If -s, --states not k=3, recommend NOT calling peaks, use bedgraph.")

    # hmm_binsize > 0
    if options.hmm_binsize <=0:
        logger.error(" --binsize must be larger than 0.")
        sys.exit( 1 )

    # hmm_lower less than hmm_upper, non-negative 
    if options.hmm_lower <0:
        logger.error(" -l, --lower should not be negative! ")
        sys.exit( 1 )
    if options.hmm_upper <0:
        logger.error(" -u, --upper should not be negative! ")
        sys.exit( 1 )
    if options.hmm_lower > options.hmm_upper:
        logger.error("Upper limit of fold change range should be greater than lower limit!" % options.mfold)
        sys.exit(1)
    
    # hmm_maxTrain non-negative
    if options.hmm_maxTrain <= 0:
        logger.error(" --maxTrain should be larger than 0!")
        sys.exit( 1 )
    
    # hmm_training_regions
    #if options.hmm_training_regions:
    #    options.argtxt += "# Using -t, --training input to train HMM instead of using fold change settings to select. \n"
    
    # hmm_zscore non-negative
    #if options.hmm_zscore <0:
    #    logger.error(" -z, --zscore should not be negative!")
    #    sys.exit( 1 )
    
    # hmm_randomSeed
    if options.hmm_randomSeed:
        options.argtxt += "# Random seed selected as: %d\n" % options.hmm_randomSeed
    
    # hmm_window non-negative
    #if options.hmm_window <0:
    #    logger.error(" --window should not be negative! ")
    #    sys.exit( 1 )

    # hmm_file
    #if options.hmm_file:
    #    options.argtxt += "# HMM training will be skipped, --model input used instead. \n"

    # hmm_modelonly
    #if options.hmm_modelonly:
    #    options.argtxt += "# Program will stop after generating model, which can be later applied with '--model'. \n"


    # Peak Calling
    if options.prescan_cutoff <= 1:
        logger.error(" In order to use -c or --prescan-cutoff, the cutoff must be larger than 1.")
        sys.exit( 1 )
    
    if options.openregion_minlen < 0: # and options.store_peaks == True:
        logger.error(" In order to use --minlen, the length should not be negative.")
        sys.exit( 1 )

    #if options.call_score.lower() not in [ 'max', 'ave', 'med', 'fc', 'zscore', 'all']:
    #    logger.error( " Invalid method: %s" % options.call_score )
    #    sys.exit( 1 )

    # call_threshold non-negative
    #if options.call_threshold <0:
    #    logger.error(" --threshold should not be negative! ")
    #    sys.exit( 1 )
    

    # Misc
    # misc_blacklist 
    #if options.misc_keep_duplicates:
    #    options.argtxt += "# Duplicate reads from analysis will be stored. \n"

    # misc_trim non-negative
    #if options.misc_trim <0:
    #    logger.error(" --trim should not be negative! ")
    #    sys.exit( 1 )

    # np # should this be mp? non-negative
    #if options.np <0:
    #    logger.error(" -m, --multiple-processing should not be negative! ")
    #    sys.exit( 1 )
    
    # min_map_quality non-negative
    #if options.min_map_quality <0:
    #    logger.error(" -q, --minmapq should not be negative! ")
    #    sys.exit( 1 )

    
    return options


