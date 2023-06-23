# Time-stamp: <2020-11-30 16:14:14 Tao Liu>

"""Description: refine peak summits

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
from collections import Counter

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_refinepeak
from MACS3.Signal.Prob import binomial_cdf_inv
from MACS3.IO.PeakIO import PeakIO

# ------------------------------------
# Main function
# ------------------------------------
def run( o_options ):
    """The Main function/pipeline for duplication filter.

    """
    # Parse options...
    options = opt_validate_refinepeak( o_options )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    if options.ofile:
        outputfile = open( os.path.join( options.outdir, options.ofile ), 'w' )
        options.oprefix = options.ofile
    else:
        outputfile = open( os.path.join( options.outdir, "%s_refinepeak.bed" % options.oprefix), "w" )


    peakio = open(options.bedfile,"rb")
    peaks = PeakIO()
    for l in peakio:
        fs = l.rstrip().split()
        peaks.add( fs[0], int(fs[1]), int(fs[2]), name=fs[3] )

    peaks.sort()
    peakio.close()

    #1 Read tag files
    info("read tag files...")
    fwtrack = load_tag_files_options (options)

    retval = fwtrack.compute_region_tags_from_peaks( peaks, find_summit, window_size = options.windowsize, cutoff = options.cutoff )
    outputfile.write( (b"\n".join( [b"%s\t%d\t%d\t%s\t%.2f" % x for x in retval] )).decode() )
    outputfile.close()
    info("Done!")

def find_summit(chrom, plus, minus, peak_start, peak_end, name = b"peak", window_size=100, cutoff = 5):

    left_sum = lambda strand, pos, width = window_size: sum([strand[x] for x in strand if x <= pos and x >= pos - width])
    right_sum = lambda strand, pos, width = window_size: sum([strand[x] for x in strand if x >= pos and x <= pos + width])
    left_forward = lambda strand, pos: strand.get(pos,0) - strand.get(pos-window_size, 0)
    right_forward = lambda strand, pos: strand.get(pos + window_size, 0) - strand.get(pos, 0)

    watson, crick = (Counter(plus), Counter(minus))
    watson_left = left_sum(watson, peak_start)
    crick_left = left_sum(crick, peak_start)
    watson_right = right_sum(watson, peak_start)
    crick_right = right_sum(crick, peak_start)

    wtd_list = []
    for j in range(peak_start, peak_end+1):
        wtd_list.append(2 * (watson_left * crick_right)**0.5 - watson_right - crick_left)
        watson_left += left_forward(watson, j)
        watson_right += right_forward(watson, j)
        crick_left += left_forward(crick, j)
        crick_right += right_forward(crick,j)

    wtd_max_val = max(wtd_list)
    wtd_max_pos = wtd_list.index(wtd_max_val) + peak_start

    #return (chrom, wtd_max_pos, wtd_max_pos+1, wtd_max_val)

    if wtd_max_val > cutoff:
        return (chrom, wtd_max_pos, wtd_max_pos+1, name+b"_R" , wtd_max_val) # 'R'efined
    else:
        return (chrom, wtd_max_pos, wtd_max_pos+1, name+b"_F" , wtd_max_val) # 'F'ailed

def load_tag_files_options ( options ):
    """From the options, load alignment tags.

    """
    options.info("# read treatment tags...")
    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    treat = tp.build_fwtrack()
    #treat.sort()
    if len(options.ifile) > 1:
        # multiple input
        for tfile in options.ifile[1:]:
            tp = options.parser(tfile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )
            #treat.sort()
    treat.finalize()
    return treat

