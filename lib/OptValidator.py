# Time-stamp: <2011-12-15 12:39:02 Tao Liu>

"""Module Description

Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import os
import re
import logging
from subprocess import Popen, PIPE
from math import log
from MACS14.IO.Parser import BEDParser, ELANDResultParser, ELANDMultiParser, ELANDExportParser, PairEndELANDMultiParser, SAMParser, BAMParser, BowtieParser,  guess_parser
# ------------------------------------
# constants
# ------------------------------------

efgsize = {"hs":2.7e9,
           "mm":1.87e9,
           "ce":9e7,
           "dm":1.2e8}

# ------------------------------------
# Misc functions
# ------------------------------------
def opt_validate ( optparser ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()

    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logging.error("Error when interpreting --gsize option: %s" % options.gsize)
            logging.error("Available shortcuts of effective genome sizes are %s" % ",".join(efgsize.keys()))
            sys.exit(1)


    # treatment file
    if not options.tfile:       # only required argument
        optparser.print_help()
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
    elif options.format == "ELANDMULTIPET":
        options.parser = PairEndELANDMultiParser
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
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # for ELANDMULTIPET format
    if options.format == "ELANDMULTIPET":
        # treatment files
        fs = options.tfile.split(',')
        if len(fs) != 2:
            logging.error("Only two filenames are acceptable! But you provided %d:'%s'" % (len(fs),options.tfile))
            sys.exit(1)
        options.tfile = fs
        if not os.path.isfile(options.tfile[0]) or not os.path.isfile(options.tfile[1]):
            logging.error("No such file: %s or %s!" % (options.tfile[0],options.tfile[1]))
            sys.exit(1)
        # input files
        if options.cfile:
            fs = options.cfile.split(',')
            if len(fs) != 2:
                logging.error("Only two filenames are acceptable! But you provided %d:'%s'" % (len(fs),options.cfile))
                sys.exit(1)
            options.cfile = fs
            if not os.path.isfile(options.cfile[0]) or not os.path.isfile(options.cfile[1]):
                logging.error("No such file: %s or %s!" % (options.cfile[0],options.cfile[1]))
                sys.exit(1)
            if not options.largelocal:
                options.largelocal = 10000
        else:
            if not options.largelocal:
                options.largelocal = 10000            
    else:
        if not os.path.isfile (options.tfile):
            logging.error("No such file: %s!" % options.tfile)
            sys.exit(1)

        # input file
        if options.cfile:
            if not os.path.isfile (options.cfile):
                logging.error("No such file: %s!" % options.cfile)
                sys.exit(1)
            if not options.largelocal:
                options.largelocal = 10000
        else:
            if not options.largelocal:
                options.largelocal = 10000                        

    # lambda set
    #lambdaset_hit =  re.match("^(\d+)\,(\d+)\,(\d+)$",options.lambdaset)
    #if lambdaset_hit:
    #    lambdaset = sorted(map(int,lambdaset_hit.groups()))
    #else:
    #    logging.error("Lambda set string must be like \"1000,5000,10000\", which means three integers seperated by commas.")
    #    sys.exit(1)
    #options.lambdaset = lambdaset # overwrite it

    # duplicate reads
    if options.keepduplicates != "auto" and options.keepduplicates != "all":
        if not options.keepduplicates.isdigit():
            logging.error("--keep-dup should be 'auto', 'all' or an integer!")
            sys.exit(1)

    # small sample?

    # shiftsize>0
    if options.shiftsize <=0 :
        logging.error("--shiftsize must > 0!")
        sys.exit(1)

    # -10*log10 pvalue
    options.log_pvalue = log(options.pvalue,10)*-10
    
    # uppercase the format string 
    options.format = options.format.upper()

    # upper and lower mfold
    try:
        (options.lmfold,options.umfold) = map(int, options.mfold.split(","))
    except:
        logging.error("mfold format error! Your input is '%s'. It should be like '10,30'." % options.mfold)
        sys.exit(1)
    
    # check callsubpeaks option should be combined with store_wig
    if options.callsubpeaks:
        if not (options.store_wig or options.store_bdg):
            logging.error("In order to use --call-subpeaks, you need to set --wig or --bdg as well!")
            sys.exit(1)
        try:
            p = Popen(["PeakSplitter","--version"],stdout=PIPE,stderr=PIPE)
        except OSError as (errno,detail):
            logging.error("PeakSplitter can not be accessed through the command line!")
            logging.error("OS gave me this error message: [Errno %d]: %s" % (errno,detail))
            logging.error("Please download PeakSplitter from <http://www.ebi.ac.uk/bertone/software/PeakSplitter_Cpp.tar.gz>.")
            logging.error("After you download and unzip it, read the instruction at <http://www.ebi.ac.uk/bertone/software/PeakSplitter_Cpp_usage.txt>")
            logging.error("Install PeakSplitter so that you can access it through typing 'PeakSplitter' under the command line.")
            sys.exit(1)

    # output filenames
    options.peakxls = options.name+"_peaks.xls"
    options.peakbed = options.name+"_peaks.bed"
    options.summitbed = options.name+"_summits.bed"
    options.zwig_tr = options.name+"_treat_afterfiting"
    options.zwig_ctl= options.name+"_control_afterfiting"
    options.negxls  = options.name+"_negative_peaks.xls"
    options.diagxls = options.name+"_diag.xls"
    options.modelR  = options.name+"_model.r"

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
	
    options.error   = logging.critical		# function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    # options argument text
    #if options.uniquetag:
    #    tmptxt = "ON"
    #else:
    #    tmptxt = "OFF"
    options.argtxt = "\n".join((
        "# ARGUMENTS LIST:",\
        "# name = %s" % (options.name),\
        "# format = %s" % (options.format),\
        "# ChIP-seq file = %s" % (options.tfile),\
        "# control file = %s" % (options.cfile),\
        "# effective genome size = %.2e" % (options.gsize),\
        #"# tag size = %d" % (options.tsize),\
        "# band width = %d" % (options.bw),\
        "# model fold = %s" % (options.mfold),\
        "# pvalue cutoff = %.2e\n" % (options.pvalue),\
        #"# unique tag filter is %s" % (tmptxt),\
        ))

    options.tosmall = not options.tolarge

    if options.tosmall:
        options.argtxt += "# Large dataset will be scaled towards smaller dataset.\n"
    else:
        options.argtxt += "# Small dataset will be scaled towards larger dataset.\n"


    if options.cfile:
        options.argtxt += "# Range for calculating regional lambda is: %d bps and %d bps\n" % (options.smalllocal,options.largelocal)
    else:
        options.argtxt += "# Range for calculating regional lambda is: %d bps\n" % (options.largelocal)

    # wig or bdg subdir
    if options.store_wig or options.store_bdg:
        if options.store_bdg:
            # bdg has higher priority
            subdir = options.name+"_MACS_bedGraph"
        elif options.store_wig:
            subdir = options.name+"_MACS_wiggle"
        # check subdir
        if os.path.exists(subdir):
            options.error("./%s exists! Unable to create directory to store profiles!!" % (subdir))
            sys.exit(1)
        else:
            options.wig_dir_tr = subdir+"/treat/"
            options.wig_dir_ctl= subdir+"/control/"

    return options
