# Time-stamp: <2011-07-08 13:54:50 Tao Liu>

"""Module for Feature IO classes.

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
import numpy as np
from numpy import int64,int32,float32

from libc.math cimport sqrt,log10

from MACS2.Constants import *
from MACS2.cProb cimport poisson_cdf
from MACS2.IO.cPeakIO import PeakIO

# ------------------------------------
# constants
# ------------------------------------
__version__ = "scoreTrackI $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "scoreTrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class scoreTrackI:
    """Class for scoreGraph type data. Modified from
    bedGraphTrackI. The only difference is that we store pvalue score,
    qvalue score and foldchange together.

    In bedGraph, data are represented as continuous non-overlapping
    regions in the whole genome. I keep this assumption in all the
    functions. If data has overlaps, some functions will definitely
    give incorrect results.

    1. Continuous: the next region should be after the previous one
    unless they are on different chromosomes;
    
    2. Non-overlapping: the next region should never have overlaps
    with preceding region.

    The way to memorize bedGraph data is to remember the transition
    points together with values of their preceding regions. The last
    data point may exceed chromosome end, unless a chromosome
    dictionary is given. Remember the coordinations in bedGraph and
    this class is 0-indexed and right-open.
    
    """
    def __init__ (self):
        """Different with bedGraphTrackI, missing values are simply
        replaced with 0.
        
        """
        self.data = {}
        self.pointer = {}

    def add_chromosome ( self, chrom, chrom_max_len ):
        if not self.data.has_key(chrom):
            self.data[chrom] = np.zeros(chrom_max_len,dtype=[('pos','int64'),
                                                             ('sample','float32'),
                                                             ('control','float32'),
                                                             ('-100logp','int32'),
                                                             ('-100logq','int32')])
            self.pointer[chrom] = 0

    def add (self,chromosome,endpos,sample,control):
        """Add a chr-endpos-sample-control block into data
        dictionary. At the mean time, calculate pvalues.

        """
        c = self.data[chromosome]
        i = self.pointer[chromosome]
        # get the preceding region
        c[i] = (endpos,sample,control,int(-100*poisson_cdf(sample,control,False,True)),0)
        self.pointer[chromosome] += 1

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([end pos],[value])
        """
        if self.data.has_key(chromosome):
            return self.data[chromosome]
        else:
            return None

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.data.keys())
        return l

    def write_bedGraph (self, fhd, name, description, colname):
        """Write all data to fhd in Wiggle Format.

        fhd: a filehandler to save bedGraph.
        name/description: the name and description in track line.

        colname: can be 'sample','control','-100logp','-100logq'

        shift will be used to shift the coordinates. default: 0
        """
        if colname not in ['sample','control','-100logp','-100logq']:
            raise Exception("%s not supported!" % colname)
        if colname in ['-100logp', '-100logq']:
            flag100 = True              # for pvalue or qvalue, divide them by 100 while writing to bedGraph file
        else:
            flag100 = False
        chrs = self.get_chr_names()
        for chrom in chrs:
            d = self.data[chrom]
            l = self.pointer[chrom]
            pre = 0
            pos   = d['pos']
            if flag100:
                value = d[colname]/100.0
            else:
                value = d[colname]
            for i in xrange( l ):
                fhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,pre,pos[i],value[i]))
                pre = pos[i]

        return True

    def __calculate_fold_change ( self, chrom, index ):
        """From 'sample' and 'control' columns, calculate foldchanges.

        chrom: chromosome name
        index: index in data[chrom]
        
        """
        return self.data[chrom]['sample'][index]/self.data[chrom]['control'][index]

    def make_pq_table ( self ):
        """Make pvalue-qvalue table.

        Step1: get all pvalue and length of block with this pvalue
        Step2: Sort them
        Step3: Apply AFDR method to adjust pvalue and get qvalue for each pvalue

        Return a dictionary of {-100log10pvalue:(-100log10qvalue,rank)} relationships.
        """
        n = self.total()
        value_list = np.empty( n, dtype = [('v', '<f4'), ('l', '<i4')])

        i = 0                           # index for value_list
        for chrom in self.data.keys():
            # for each chromosome
            pre_p  = 0
            pos    = self.data[chrom][ 'pos' ]
            value  = self.data[chrom][ '-100logp' ]
            length = self.pointer[chrom]
            for j in xrange( length ):
                # for each region
                this_l = pos[j]-pre_p
                this_v = value[j]
                value_list[i] = (this_v,this_l)
                pre_p = pos[j]
                i += 1
        # sort
        value_list.sort(order='v')
        
        N = sum(value_list['l'])
        k = 1                           # rank
        #S_q = S_p + log10(k)-log10(N)
        f = -log10(N)
        pre_v = -1e100
        pre_k = 0
        pre_l = 0
        pvalue2qvalue = {pre_v:[0,k,0]}              # pvalue:[qvalue,rank,bp_with_this_pvalue]
        for i in xrange(value_list.size-1,-1,-1):
            (v,l) = value_list[i]
            if v != pre_v:
                # new value
                q = v+int((log10(k)+f)*100) # we save integars here.
                pvalue2qvalue[v] = [max(0,q),k, 0]
                pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
                pre_v = v
            k+=l
        pvalue2qvalue[pre_v][2] = k-pvalue2qvalue[pre_v][1]
        # pop the first -1e100 one
        pvalue2qvalue.pop(-1e100)

        return pvalue2qvalue

    def assign_qvalue ( self , pvalue2qvalue ):
        """Assign -100log10qvalue to every point.

        pvalue2qvalue: a dictionary of -100log10pvalue:-100log10qvalue
        """
        t = 0
        for chrom in self.data.keys():
            pvalue = self.data[chrom]['-100logp']
            qvalue = self.data[chrom]['-100logq']
            for i in xrange( self.pointer[chrom] ):
                qvalue[i] = pvalue2qvalue[pvalue[i]][0]
        return True

    def call_peaks (self, cutoff=50, min_length=200, max_gap=50, colname='-100logp'):
        """This function try to find regions within which, scores
        are continuously higher than a given cutoff.

        This function is NOT using sliding-windows. Instead, any
        regions in bedGraph above certain cutoff will be detected,
        then merged if the gap between nearby two regions are below
        max_gap. After this, peak is reported if its length is above
        min_length.

        cutoff:  cutoff of value, default 1.
        min_length :  minimum peak length, default 200.
        gap   :  maximum gap to merge nearby peaks, default 50.
        colname: can be 'sample','control','-100logp','-100logq'. Cutoff will be applied to the specified column.
        """
        assert (colname in [ 'sample', 'control', '-100logp', '-100logq' ]), "%s not supported!" % colname

        chrs  = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            chrom_pointer = self.pointer[chrom]
            chrom_d       = self.get_data_by_chr( chrom ) # arrays for position and values
            chrom_pos     = chrom_d[ 'pos' ]
            chrom_score   = chrom_d[ colname ]
            chrom_sample  = chrom_d[ 'sample' ]
            chrom_control = chrom_d[ 'control' ]
            chrom_pvalue  = chrom_d[ '-100logp' ]
            chrom_qvalue  = chrom_d[ '-100logq' ]

            x     = 0
            pre_p = 0                   # remember previous position
            peak_content = None         # to store points above cutoff
            
            while True and x < chrom_pointer:
                # find the first region above cutoff
                # try to read the first data range for this chrom
                p = chrom_pos[ x ]
                v = chrom_score[ x ]
                x += 1                  # index for the next point
                if v >= cutoff:
                    peak_content = [ ( pre_p, p, v, x ), ] # remember the index too...
                    pre_p = p
                    break               # found the first range above cutoff
                else:
                    pre_p = p

            for i in xrange( x, chrom_pointer ):
                # continue scan the rest regions
                p = chrom_pos[ i ]
                v = chrom_score[ i ]
                if v < cutoff:
                    pre_p = p
                    continue
                # for points above cutoff
                # if the gap is allowed
                if pre_p - peak_content[ -1 ][ 1 ] <= max_gap:
                    peak_content.append( ( pre_p, p, v, i ) )
                else:
                    # when the gap is not allowed, close this peak
                    peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
                    if peak_length >= min_length: # if the peak is too small, reject it
                        tmpsummit = []
                        summit_pos   = None
                        summit_value = None
                        for (tmpstart,tmpend,tmpvalue,tmpindex) in peak_content:
                            if not summit_value or summit_value < tmpvalue:
                                tmpsummit = [ int(( tmpend+tmpstart )/2), ]
                                tmpsummit_index = [ tmpindex, ]
                                summit_value = tmpvalue
                            elif summit_value == tmpvalue:
                                # remember continuous summit values
                                tmpsummit.append( int( (tmpend+tmpstart)/2 ) )
                                tmpsummit_index.append( tmpindex )
                        middle_summit = int( ( len(tmpsummit)+1 )/2 )-1
                        summit_pos    = tmpsummit[ middle_summit ]
                        summit_index  = tmpsummit_index[ middle_summit ]
                        # char * chromosome, long start, long end, long summit = 0, 
                        # double peak_height=0, int pileup=0, 
                        # double pvalue=0, double fold_change=0, double qvalue=0
                        peaks.add( chrom,
                                   peak_content[0][0],
                                   peak_content[-1][1],
                                   summit      = summit_pos,
                                   peak_score  = summit_value,
                                   pileup      = chrom_sample[ summit_index ],
                                   pscore      = chrom_pvalue[ summit_index ]/100.0,
                                   fold_change = chrom_sample[ summit_index ]/chrom_control[ summit_index ],
                                   qscore      = chrom_qvalue[ summit_index ]/100.0,
                                   )
                    # start a new peak
                    peak_content = [ ( pre_p, p, v, i ), ]
                pre_p = p
                
            # save the last peak
            if not peak_content:
                continue
            peak_length = peak_content[ -1 ][ 1 ] - peak_content[ 0 ][ 0 ]
            if peak_length >= min_length: # if the peak is too small, reject it
                summit_pos = None
                summit_value = None
                for (tmpstart,tmpend,tmpvalue,tmpindex) in peak_content:
                    if not summit_value or summit_value < tmpvalue:
                        tmpsummit = [ int(( tmpend+tmpstart )/2), ]
                        tmpsummit_index = [ tmpindex, ]
                        summit_value = tmpvalue
                    elif summit_value == tmpvalue:
                        # remember continuous summit values
                        tmpsummit.append( int( (tmpend+tmpstart)/2 ) )
                        tmpsummit_index.append( tmpindex )
                middle_summit = int( ( len(tmpsummit)+1 )/2 )-1
                summit_pos    = tmpsummit[ middle_summit ]
                summit_index  = tmpsummit_index[ middle_summit ]

                peaks.add( chrom,
                           peak_content[0][0],
                           peak_content[-1][1],
                           summit      = summit_pos,
                           peak_score  = summit_value,
                           pileup      = chrom_sample[ summit_index ],
                           pscore      = chrom_pvalue[ summit_index ]/100.0,
                           fold_change = chrom_sample[ summit_index ]/chrom_control[ summit_index ],
                           qscore      = chrom_qvalue[ summit_index ]/100.0,
                           )
            
        return peaks

    def total ( self ):
        """Return the number of regions in this object.

        """
        t = 0
        for chrom in self.data.keys():
            t += self.pointer[chrom]
        return t
