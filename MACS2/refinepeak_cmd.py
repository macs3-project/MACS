# Time-stamp: <2019-09-20 11:36:55 taoliu>

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
import logging
from collections import Counter

# ------------------------------------
# own python modules
# ------------------------------------
from MACS2.OptValidator import opt_validate_refinepeak as opt_validate
from MACS2.Prob import binomial_cdf_inv
from MACS2.IO.BedGraphIO import bedGraphIO,genericBedIO
from MACS2.IO.PeakIO import PeakIO
from MACS2.Constants import *


# ------------------------------------
# Main function
# ------------------------------------
def run( o_options ):
    """The Main function/pipeline for duplication filter.
    
    """
    # Parse options...
    options = opt_validate( o_options )
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


    peakio = file(options.bedfile)
    peaks = PeakIO()
    for l in peakio:
        fs = l.rstrip().split()
        peaks.add( fs[0], int(fs[1]), int(fs[2]), name=fs[3] )

    peaks.sort()
    
    #1 Read tag files
    info("read tag files...")
    fwtrack = load_tag_files_options (options)
    
    retval = fwtrack.compute_region_tags_from_peaks( peaks, find_summit, window_size = options.windowsize, cutoff = options.cutoff )
    outputfile.write( "\n".join( map(lambda x: "%s\t%d\t%d\t%s\t%.2f" % x , retval) ) )
    info("Done!")

def find_summit(chrom, plus, minus, peak_start, peak_end, name = "peak", window_size=100, cutoff = 5):
    
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
        return (chrom, wtd_max_pos, wtd_max_pos+1, name+"_R" , wtd_max_val) # 'R'efined
    else:
        return (chrom, wtd_max_pos, wtd_max_pos+1, name+"_F" , wtd_max_val) # 'F'ailed

    #return "{}\t{}\t{}\tRefinePeak_summit\t{:.2f}\n".format(chrom,
    #                                                        wtd_max_pos,
    #                                                        wtd_max_pos+1,
    #                                                        wtd_max_val,)



# def find_summit(bed_file, sam_file, window_size, output_file):
#     def count_by_strand(ialign):
#         pred = lambda x:x.is_reverse
#         watson_5_end = lambda x:x.pos
#         crick_5_end = lambda x:x.aend
#         ialign1, ialign2 = tee(ialign)

#         return (Counter(map(watson_5_end,
#                             ifilterfalse(pred, ialign1))),
#                 Counter(map(crick_5_end,
#                             ifilter(pred, ialign2))))
    
#     left_sum = lambda strand, pos, width = window_size: sum([strand[x] for x in strand if x <= pos and x >= pos - width])
#     right_sum = lambda strand, pos, width = window_size: sum([strand[x] for x in strand if x >= pos and x <= pos + width])
#     left_forward = lambda strand, pos: strand.get(pos,0) - strand.get(pos-window_size, 0)
#     right_forward = lambda strand, pos: strand.get(pos + window_size, 0) - strand.get(pos, 0)
#     samfile = pysam.Samfile(sam_file, "rb" )

#     cnt = 0
#     with open(bed_file) as bfile, open(output_file,"w") as ofile:
#         for i in bfile:
#             i = i.split("\t")
#             chrom = i[0]
#             peak_start = int(i[1])
#             peak_end = int(i[2])
            
#             watson, crick = count_by_strand(samfile.fetch(chrom, peak_start-window_size, peak_end+window_size))
#             watson_left = left_sum(watson, peak_start)
#             crick_left = left_sum(crick, peak_start)
#             watson_right = right_sum(watson, peak_start)
#             crick_right = right_sum(crick, peak_start)

#             wtd_list = []
#             for j in range(peak_start, peak_end+1):
#                 wtd_list.append(2 * sqrt(watson_left * crick_right) - watson_right - crick_left)
#                 watson_left += left_forward(watson, j)
#                 watson_right += right_forward(watson, j)
#                 crick_left += left_forward(crick, j)
#                 crick_right += right_forward(crick,j)

#             wtd_max_val = max(wtd_list)
#             wtd_max_pos = wtd_list.index(wtd_max_val) + peak_start
#             cnt += 1

#             ofile.write("{}\t{}\t{}\tSPP_summit_{}\t{:.2f}\n".format(chrom,
#                                                                      wtd_max_pos,
#                                                                      wtd_max_pos+1,
#                                                                      cnt,
#                                                                      wtd_max_val,))
#     samfile.close()



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

