# Time-stamp: <2022-09-15 17:26:49 Tao Liu>

"""Module Description: Functions to write file

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import os
import sys
from array import array as pyarray

from MACS3.Utilities.Constants import MACS_VERSION

# ------------------------------------
# constants
# ------------------------------------
import logging
import MACS3.Utilities.Logger

logger = logging.getLogger(__name__)
debug   = logger.debug
info    = logger.info
# ------------------------------------
# Misc functions
# ------------------------------------
def zwig_write (trackI, subdir, fileprefix, d, log=None,space=10, single=False):
    """Write shifted tags information in wiggle file in a given
    step. Then compress it using 'gzip' program.

    trackI: shifted tags from PeakDetect object
    subdir: directory where to put the wiggle file
    fileprefix: wiggle file prefix
    d     : d length
    log   : logging function, default is sys.stderr.write
    space : space to write tag number on spots, default 10
    """
    chrs = trackI.get_chr_names()
    os.makedirs (subdir)
    step = 10000000 + 2*d

    if single:
        info("write to a wiggle file")
        f = os.path.join(subdir,fileprefix+"_all"+".wig")
        wigfhd = open(f,"w")
        wigfhd.write("track type=wiggle_0 name=\"%s_all\" description=\"Extended tag pileup from MACS version %s for every %d bp\"\n" % (fileprefix.replace('_afterfiting',''), MACS_VERSION, space)) # data type line

    for chrom in sorted(chrs):
        if not single:
            f = os.path.join(subdir,fileprefix+"_"+chrom+".wig")
            info("write to "+f+" for chromosome "+chrom)
            wigfhd = open(f,"w")
            # suggested by dawe
            wigfhd.write("track type=wiggle_0 name=\"%s_%s\" description=\"Extended tag pileup from MACS version %s for every %d bp\"\n" % ( fileprefix.replace('_afterfiting',''), chrom, MACS_VERSION, space)) # data type line
        else:
            info("write data for chromosome "+chrom)

        wigfhd.write("variableStep chrom=%s span=%d\n" % (chrom,space))
        tags = trackI.get_locations_by_chr(chrom)[0]
        l = len(tags)
        window_counts = pyarray('L',[0]*step)
        startp = -1*d
        endp   = startp+step
        index_tag = 0

        while index_tag<l:
            s = tags[index_tag]-d/2     # start of tag
            e = s+d                     # end of tag

            if e < endp:
                # project tag to window_counts line
                ps = s-startp # projection start
                pe = ps+d     # projection end
                for i in range(ps,pe):
                    window_counts[i] += 1
                index_tag += 1
            else:
                # write it to zwig file then reset parameters
                # keep this tag for next window
                for i in range(d,step-d,space):
                    if window_counts[i] == 0:
                        pass
                    else:
                        wigfhd.write("%d\t%d\n" % (i+startp+1,window_counts[i]))
                # reset
                window_counts_next = pyarray('L',[0]*step)
                # copy d values from the tail of previous window to next window
                for n,i in enumerate(range(step-2*d,step)): # debug
                    window_counts_next[n] = window_counts[i]
                window_counts = window_counts_next
                startp = endp - 2*d
                endp = startp+step
        # last window
        for i in range(d,step-d,space):
            if window_counts[i] == 0:
                pass
            else:
                wigfhd.write("%d\t%d\n" % (i+startp+1,window_counts[i]))
        if not single:
            wigfhd.close()
            info("compress the wiggle file using gzip...")
            os.system("gzip "+f)
    if single:
        wigfhd.close()
        info("compress the wiggle file using gzip...")
        os.system("gzip "+f)


def zbdg_write (trackI, subdir, fileprefix, d, log=None, single=False):
    """Write shifted tags information in wiggle file in a given
    step. Then compress it using 'gzip' program.

    trackI: shifted tags from PeakDetect object
    subdir: directory where to put the wiggle file
    fileprefix: wiggle file prefix
    d     : d length
    log   : logging function, default is sys.stderr.write
    space : space to write tag number on spots, default 10
    """
    if not log:
        log = lambda x: sys.stderr.write(x+"\n")
    chrs = trackI.get_chr_names()
    os.makedirs (subdir)
    step = 10000000 + 2*d

    if single:
        info("write to a bedGraph file")
        f = os.path.join(subdir,fileprefix+"_all"+".bdg")
        bdgfhd = open(f,"w")
        bdgfhd.write("track type=bedGraph name=\"%s_all\" description=\"Extended tag pileup from MACS version %s\"\n" % (fileprefix.replace('_afterfiting',''), MACS_VERSION)) # data type line

    for chrom in sorted(chrs):
        if not single:
            f = os.path.join(subdir,fileprefix+"_"+chrom+".bdg")
            info("write to "+f+" for chromosome "+chrom)
            bdgfhd = open(f,"w")
            bdgfhd.write("track type=bedGraph name=\"%s_%s\" description=\"Extended tag pileup from MACS version %s\"\n" % (fileprefix.replace('_afterfiting',''), chrom, MACS_VERSION)) # data type line
        else:
            info("write data for chromosome "+chrom)

        tags = trackI.get_locations_by_chr(chrom)[0]
        l = len(tags)
        window_counts = pyarray('L',[0]*step)
        startp = -1*d
        endp   = startp+step
        index_tag = 0

        while index_tag<l:
            s = tags[index_tag]-d/2     # start of tag
            e = s+d                     # end of tag

            if e < endp:
                # project tag to window_counts line
                ps = s-startp # projection start
                pe = ps+d     # projection end
                for i in range(ps,pe):
                    window_counts[i] += 1
                index_tag += 1
            else:
                # write it to zbdg file then reset parameters
                # keep this tag for next window
                prev = window_counts[d]
                left = startp+d
                right = left+1
                for i in range(d+1,step-d):
                    if window_counts[i] == prev:
                        # same value, extend
                        right += 1
                    else:
                        # diff value, close
                        if prev != 0:
                            bdgfhd.write("%s\t%d\t%d\t%d\n" % (chrom,left,right,prev))
                        prev = window_counts[i]
                        left = right
                        right = left + 1
                # last bin
                if prev != 0:
                    bdgfhd.write("%s\t%d\t%d\t%d\n" % (chrom,left,right,prev))

                # reset
                window_counts_next = pyarray('L',[0]*step)
                # copy d values from the tail of previous window to next window
                for n,i in enumerate(range(step-2*d,step)): # debug
                    window_counts_next[n] = window_counts[i]
                window_counts = window_counts_next
                startp = endp - 2*d
                endp = startp+step
        # last window
        prev = window_counts[d]
        left = startp+d
        right = left+1
        for i in range(d+1,step-d):
            if window_counts[i] == prev:
                # same value, exrend
                right += 1
            else:
                # diff value, close
                if prev != 0:
                    bdgfhd.write("%s\t%d\t%d\t%d\n" % (chrom,left,right,prev))
                prev = window_counts[i]
                left = right
                right = left + 1
        # last bin
        if prev != 0:
            bdgfhd.write("%s\t%d\t%d\t%d\n" % (chrom,left,right,prev))

        if not single:
            bdgfhd.close()
            info("compress the bedGraph file using gzip...")
            os.system("gzip "+f)
    if single:
        bdgfhd.close()
        info("compress the bedGraph file using gzip...")
        os.system("gzip "+f)


def model2r_script(model,filename,name):
    rfhd = open(filename,"w")
    p = model.plus_line
    m = model.minus_line
    ycorr = model.ycorr
    xcorr = model.xcorr
    alt_d = model.alternative_d
    #s = model.shifted_line
    d = model.d
    w = len(p)
    norm_p = [0]*w
    norm_m = [0]*w
    #norm_s = [0]*w
    sum_p = sum(p)
    sum_m = sum(m)
    #sum_s = sum(s)
    for i in range(w):
        norm_p[i] = float(p[i])*100/sum_p
        norm_m[i] = float(m[i])*100/sum_m
        #norm_s[i] = float(s[i])*100/sum_s
    rfhd.write("# R script for Peak Model\n")
    rfhd.write("#  -- generated by MACS\n")

    rfhd.write("""p <- c(%s)
m <- c(%s)
ycorr <- c(%s)
xcorr <- c(%s)
altd  <- c(%s)
x <- seq.int((length(p)-1)/2*-1,(length(p)-1)/2)
pdf('%s_model.pdf',height=6,width=6)
plot(x,p,type='l',col=c('red'),main='Peak Model',xlab='Distance to the middle',ylab='Percentage')
lines(x,m,col=c('blue'))
legend('topleft',c('forward tags','reverse tags'),lty=c(1,1,1),col=c('red','blue'))
plot(xcorr,ycorr,type='l',col=c('black'),main='Cross-Correlation',xlab='Lag between + and - tags',ylab='Correlation')
abline(v=altd,lty=2,col=c('red'))
legend('topleft','alternative lag(s)',lty=2,col='red')
legend('right','alt lag(s) : %s',bty='n')
dev.off()
""" % (','.join(map(str,norm_p)),
       ','.join(map(str,norm_m)),
       ','.join(map(str,ycorr)),
       ','.join(map(str,xcorr)),
       ', '.join(map(str,alt_d)),
       name,
       ','.join(map(str,alt_d))
       ))
    rfhd.close()

def diag_write (filename, diag_result):
    ofhd_diag = open(filename,"w")
    a = diag_result[0]
    l = len(a)-2
    s = [90-x*10 for x in range(l)]
    ofhd_diag.write("FC range\t# of Peaks\tcover by sampling %s\n" % ("%\t".join (map(str,s))+"%"))
    format = "%s\t%d"+"\t%.2f"*l+"\n"
    ofhd_diag.write( "".join( [format % tuple(x) for x in diag_result])  )
    ofhd_diag.close()


def hmmratac_writer():
    pass
# ------------------------------------
# Classes
# ------------------------------------
