# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-09-15 16:13:48 Tao Liu>

"""Module for Region classe.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
#import random
#import re
#import sys

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
__version__ = "Region $Revision$"
__author__ = "Tao Liu <vladimir.liu@gmail.com>"
__doc__ = "Region class"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
cdef class Regions:
    """For plain region of chrom, start and end
    """
    cdef:
        public dict regions
        public int total
        bool __flag_sorted

    def __init__ (self):
        self.regions= {}
        self.__flag_sorted = False
        self.total = 0

    def __getitem__ (self, chrom):
        return self.regions[ chrom ]
        
    cpdef pop ( self, int n ):
        # when called, pop the first n regions in Regions class. Self will be modified.
        cdef:
            list clist # for chromosomes
            int tmp_l
            bytes chrom
            int n_taken # remember the number of regions prepared
            object ret # returned Regions

        if self.total == 0:
            raise Exception("None left")
            
        clist = sorted(list(self.regions.keys()))
        n_taken = n
        ret = Regions()
        for chrom in sorted(clist):
            ret.regions[chrom] = self.regions[chrom][:n_taken]
            self.regions[chrom] = self.regions[chrom][n_taken:]
            if not self.regions[chrom]:
                # remove this chromosome if there is none left
                self.regions.pop( chrom )
            tmp_l = len( ret.regions[chrom] )
            ret.total +=  tmp_l
            # calculate remained
            self.total -= tmp_l
            #print( ret.total, self.total )
            n_taken -= tmp_l
            if not n_taken:
                # when there is no need, quit loop
                break
        return ret
    
    cpdef init_from_PeakIO ( self, object peaks ):
        """Initialize the object with a PeakIO object.

        Note: I intentionally forgot to check if peakio is actually a
        PeakIO...
        """
        cdef:
            bytes chrom
            list ps
            object p
            int i

        peaks.sort()
        self.total = 0
        for chrom in sorted(peaks.get_chr_names()):
            ps = peaks.get_data_from_chrom( chrom )
            self.regions[chrom] = []
            for i in range( len( ps ) ):
                p = ps[ i ]
                self.regions[chrom].append( ( p['start'], p['end'] ) )
                self.total += 1
        self.sort()

    cpdef add_loc ( self, bytes chrom, int start, int end ):
        if self.regions.has_key(chrom):
            self.regions[chrom].append( (start,end) )
        else:
            self.regions[chrom] = [ (start,end), ]
        self.total += 1
        self.__flag_sorted = False
        return

    cpdef sort (self):
        cdef:
            bytes chrom

        if self.__flag_sorted:
            return
        for chrom in sorted(self.regions.keys()):
            self.regions[chrom].sort()
        self.__flag_sorted = True

    cpdef set get_chr_names (self):
        return set( sorted(self.regions.keys()) )

    cpdef expand ( self, int flanking ):
        cdef:
            bytes chrom
            list ps
            int i

        self.sort()
        for chrom in sorted(self.regions.keys()):
            ps = self.regions[ chrom ]
            for i in range( len( ps ) ):
                ps[i] = ( max(0, ps[i][0] - flanking), ps[i][1] + flanking )
            ps.sort()
            self.regions[chrom] = ps
    
    cpdef merge_overlap ( self ):
        """
        merge overlapping regions
        """
        cdef:
            bytes chrom
            int s_new_region, e_new_region, i, j
            dict regions, new_regions
            list chrs
            list regions_chr
            tuple prev_region

        self.sort()
        regions = self.regions
        new_regions = {}
        
        chrs = list(regions.keys())
        chrs.sort()
        self.total = 0
        for i in range(len(chrs)):
            chrom = chrs[i]
            new_regions[chrom]=[]
            n_append = new_regions[chrom].append
            prev_region = ()
            regions_chr = regions[chrom]
            for i in range( len(regions_chr) ):
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
            self.total += len( new_regions[chrom] )
        self.regions = new_regions
        self.sort()
        return True

    cpdef write_to_bed (self, fhd ):
        cdef:
            int i
            bytes chrom
            tuple region

        chrs = list(self.regions.keys())
        chrs.sort()
        for i in range( len(chrs) ):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                fhd.write( "%s\t%d\t%d\n" % (chrom.decode(),region[0],region[1] ) )

    def __str__ ( self ):
        cdef:
            int i
            bytes chrom
            tuple region
            str ret
            
        ret = ""
        chrs = list(self.regions.keys())
        chrs.sort()
        for i in range( len(chrs) ):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                ret += "%s\t%d\t%d\n" % (chrom.decode(),region[0],region[1] )
        return ret
    
    # cpdef object randomly_pick ( self, int n, int seed = 12345 ):
    #     """Shuffle the regions and get n regions out of it. Return a
    #     new Regions object.  
    #     """

    #     cdef:
    #         list all_pc
    #         list chrs
    #         bytes chrom
    #         object ret_peakio, p
    #     assert n > 0
    #     chrs = sorted(list(self.peaks.keys()))
    #     all_pc = []
    #     for chrom in chrs:
    #         all_pc.extend(self.peaks[chrom])
    #     random.seed( seed )
    #     all_pc = random.shuffle( all_pc )[:n]
    #     ret_peakio = PeakIO()
    #     for p in all_pc:
    #         ret_peakio.add_PeakContent ( p["chrom"], p )
    #     return ret_peakio
    

    cpdef void exclude (self, object regionio2):
        """ Remove overlapping regions in regionio2, another Region
        object.

        """
        cdef:
            dict regions1, regions2
            list chrs1, chrs2
            bytes k
            dict ret_regionss
            bool overlap_found
            tuple r1, r2
            long n_rl1, n_rl2

        self.sort()
        regions1 = self.regions
        self.total = 0
        assert isinstance(regionio2,Regions)
        regionio2.sort()
        regions2 = regionio2.regions

        ret_regions = dict()
        chrs1 = list(regions1.keys())
        chrs2 = list(regions2.keys())
        for k in chrs1:
            #print(f"chromosome {k}")
            if not chrs2.count(k):
                # no such chromosome in peaks1, then don't touch the peaks in this chromosome
                ret_regions[ k ] = regions1[ k ]
                self.total += len( ret_regions[ k ] )
                continue
            ret_regions[ k ] = []
            n_rl1 = len( regions1[k] )
            n_rl2 = len( regions2[k] )
            rl1_k = iter( regions1[k] ).__next__
            rl2_k = iter( regions2[k] ).__next__
            overlap_found = False
            r1 = rl1_k()
            n_rl1 -= 1
            r2 = rl2_k()
            n_rl2 -= 1
            while ( True ):
                # we do this until there is no r1 or r2 left.
                if r2[0] < r1[1] and r1[0] < r2[1]:
                    # since we found an overlap, r1 will be skipped/excluded
                    # and move to the next r1
                    overlap_found = True
                    n_rl1 -= 1
                    if n_rl1 >= 0:
                        r1 = rl1_k()
                        overlap_found = False
                        continue
                    else:
                        break
                if r1[1] < r2[1]:
                    # in this case, we need to move to the next r1,
                    # we will check if overlap_found is true, if not, we put r1 in a new dict
                    if not overlap_found:
                        ret_regions[ k ].append( r1 )
                    n_rl1 -= 1
                    if n_rl1 >= 0:
                        r1 = rl1_k()
                        overlap_found = False
                    else:
                        # no more r1 left
                        break
                else:
                    # in this case, we need to move the next r2
                    if n_rl2:
                        r2 = rl2_k()
                        n_rl2 -= 1
                    else:
                        # no more r2 left
                        break
            if n_rl1 >= 0:
                ret_regions[ k ].extend( regions1[ k ][-n_rl1-1:] )
            self.total += len( ret_regions[ k ] )

        self.regions = ret_regions
        self.__flag_sorted = False
        self.sort()
        return
