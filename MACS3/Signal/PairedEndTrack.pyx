# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-09-15 17:07:26 Tao Liu>

"""Module for filter duplicate tags from paired-end data

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# Python modules
# ------------------------------------
import io
import sys
from copy import copy
from array import array as pyarray
from collections import Counter

import logging
import MACS3.Utilities.Logger

logger = logging.getLogger(__name__)
debug   = logger.debug
info    = logger.info
# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Signal.Pileup import quick_pileup, over_two_pv_array, se_all_in_one_pileup
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.PileupV2 import pileup_from_LR_hmmratac
# ------------------------------------
# Other modules
# ------------------------------------
import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, float64_t
from cpython cimport bool
cimport cython


cdef INT_MAX = <int32_t>((<uint32_t>(-1))>>1)

# We don't use the following structs anymore
# cdef packed struct peLoc:
#     int32_t l
#     int32_t r

# cdef class PETrackChromosome:
#     cdef:
#         public np.ndarray locations
#         public uint32_t pointer
#         public uint32_t buffer_size
#         public uint64_t coverage
#         public uint64_t chrlen
#         uint32_t __buffer_increment
#         bool __sorted
#         bool __destroyed

# Let numpy enforce PE-ness using ndarray, gives bonus speedup when sorting
# PE data doesn't have strandedness

cdef class PETrackI:
    """Paired End Locations Track class I along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    cdef:
        public dict __locations
        public dict __size
        public dict __buf_size
        public bool __sorted
        public uint64_t total
        public object annotation
        public dict rlengths
        public int64_t buffer_size
        public int64_t length
        public float32_t average_template_length
        bool   __destroyed

    def __init__ (self, char * anno="", int64_t buffer_size = 100000 ):
        """fw is the fixed-width for all locations.

        """
        self.__locations = {}    # dictionary with chrname as key, nparray with [('l','int32'),('r','int32')] as value
        self.__size = {}      # dictionary with chrname as key, size of the above nparray as value
        self.__buf_size = {}      # dictionary with chrname as key, size of the above nparray as value
        self.__sorted = False
        self.total = 0           # total fragments
        self.annotation = anno   # need to be figured out
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.average_template_length = 0.0

    cpdef void add_loc ( self, bytes chromosome, int32_t start, int32_t end):
        """Add a location to the list according to the sequence name.

        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        """
        cdef:
            int32_t i

        if chromosome not in self.__locations:
            self.__buf_size[chromosome] = self.buffer_size
            self.__locations[chromosome] = np.zeros(shape=self.buffer_size, dtype=[('l','int32'),('r','int32')]) # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            self.__locations[chromosome][0] = ( start, end )
            self.__size[chromosome] = 1
        else:
            i = self.__size[chromosome]
            if self.__buf_size[chromosome] == i:
                self.__buf_size[chromosome] += self.buffer_size
                self.__locations[chromosome].resize((self.__buf_size[chromosome]), refcheck = False )
            self.__locations[chromosome][ i ] = ( start, end )
            self.__size[chromosome] = i + 1
        self.length += end - start
        return

    cpdef void destroy ( self ):
        """Destroy this object and release mem.
        """
        cdef:
            set chrs
            bytes chromosome

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.__locations:
                self.__locations[chromosome].resize( self.buffer_size, refcheck=False )
                self.__locations[chromosome].resize( 0, refcheck=False )
                self.__locations[chromosome] = None
                self.__locations.pop(chromosome)
        self.__destroyed = True
        return

    cpdef bint set_rlengths ( self, dict rlengths ):
        """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this petrack object will be updated.

        If a chromosome in this petrack is not covered by given
        rlengths, and it has no associated length, it will be set as
        maximum integer.
        """
        cdef:
            set valid_chroms, missed_chroms
            bytes chrom

        valid_chroms = set(self.__locations.keys()).intersection(rlengths.keys())
        for chrom in sorted(valid_chroms):
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.__locations.keys()).difference(rlengths.keys())
        for chrom in sorted(missed_chroms):
            self.rlengths[chrom] = INT_MAX
        return True

    cpdef dict get_rlengths ( self ):
        """Get reference chromosome lengths dictionary.

        If self.rlengths is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.__locations.keys()])
        return self.rlengths

    cpdef void finalize ( self ):
        """ Resize np arrays for 5' positions and sort them in place

        Note: If this function is called, it's impossible to append more files to this FWTrack object. So remember to call it after all the files are read!
        """

        cdef:
            int32_t i
            bytes c
            set chrnames

        self.total = 0

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.__locations[c].resize((self.__size[c]), refcheck=False)
            self.__locations[c].sort( order=['l', 'r'] )
            self.total += self.__size[c]

        self.__sorted = True
        self.average_template_length = <float32_t>( self.length ) / self.total
        return

    cpdef get_locations_by_chr ( self, bytes chromosome ):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if chromosome in self.__locations:
            return self.__locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    cpdef set get_chr_names ( self ):
        """Return all the chromosome names in this track object as a python set.
        """
        return set(self.__locations.keys())


    cpdef void sort ( self ):
        """Naive sorting for locations.

        """
        cdef:
            uint32_t i
            bytes c
            set chrnames

        chrnames = self.get_chr_names()

        for c in chrnames:
            #print "before", self.__locations[c][0:100]
            self.__locations[c].sort( order=['l', 'r'] ) # sort by the leftmost location
            #print "before", self.__locations[c][0:100]
        self.__sorted = True
        return

    cpdef dict count_fraglengths ( self ):
        """Return a dictionary of the counts for sizes/fragment lengths of each pair.

        This function is for HMMRATAC.
        """
        cdef:
            np.ndarray[np.int32_t, ndim=1] sizes
            np.int32_t s
            np.ndarray locs
            list chrnames
            int i
            #dict ret_dict
            bytes k

        counter = Counter()
        chrnames = list( self.get_chr_names() )
        for i in range( len(chrnames) ):
            locs = self.__locations[ chrnames[i] ]
            sizes = locs['r'] - locs['l']
            for s in sizes:
                counter[ s ] += 1
        return dict(counter)

    cpdef np.ndarray fraglengths ( self ):
        """Return the sizes/fragment lengths of each pair.

        This function is for HMMRATAC EM training.
        """
        cdef:
            np.ndarray[np.int32_t, ndim=1] sizes
            np.ndarray locs
            list chrnames
            int i

        chrnames = list( self.get_chr_names() )
        locs = self.__locations[ chrnames[ 0 ] ]
        sizes = locs['r'] - locs['l']
        for i in range( 1, len(chrnames) ):
            locs = self.__locations[ chrnames[i] ]
            sizes = np.concatenate( ( sizes, locs['r'] - locs['l'] ) )
        return sizes    
    
    @cython.boundscheck(False) # do not check that np indices are valid
    cpdef void filter_dup ( self, int32_t maxnum=-1):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        cdef:
            int32_t i_chrom, n, start, end
            int32_t loc_start, loc_end, current_loc_start, current_loc_end
            uint64_t i
            bytes k
            np.ndarray locs
            uint64_t locs_size
            set chrnames
            np.ndarray selected_idx

        if maxnum < 0: return # condition to return if not filtering

        if not self.__sorted: self.sort()

        self.total = 0
        #self.length = 0
        self.average_template_length = 0.0
        
        chrnames = self.get_chr_names()

        for k in chrnames: # for each chromosome
            locs = self.__locations[k]
            locs_size = locs.shape[0]
            if locs_size == 1:
                # do nothing and continue
                continue
            # discard duplicate reads and make a new __locations[k]
            # initialize boolean array as all TRUE, or all being kept
            selected_idx = np.ones( locs_size, dtype=bool)
            # get the first loc
            ( current_loc_start, current_loc_end ) = locs[0]
            i = 1 # index of new_locs
            n = 1 # the number of tags in the current genomic location
            for i in range(1, locs_size):
                ( loc_start, loc_end ) = locs[i]
                if loc_start != current_loc_start or loc_end != current_loc_end:
                    # not the same, update currnet_loc_start/end/l, reset n
                    current_loc_start = loc_start
                    current_loc_end = loc_end
                    n = 1
                    continue
                else:
                    # both ends are the same, add 1 to duplicate number n
                    n += 1
                    if n > maxnum:
                        # change the flag to False
                        selected_idx[ i ] = False
                        # subtract current_loc_l from self.length
                        self.length -= current_loc_end - current_loc_start
            self.__locations[k] = locs[ selected_idx ]
            self.__size[k] = self.__locations[k].shape[0]
            self.total += self.__size[k]
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            selected_idx.resize( self.buffer_size, refcheck=False)
            selected_idx.resize( 0, refcheck=False)
        self.average_template_length = self.length / self.total
        return

    cpdef void sample_percent (self, float32_t percent, int32_t seed = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed! If a new PETrackI is wanted, use sample_percent_copy instead.
        """
        cdef:
            uint32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            bytes k
            set chrnames
            object rs, rs_shuffle

        self.total = 0
        self.length = 0
        self.average_template_length = 0.0

        chrnames = self.get_chr_names()

        if seed >= 0:
            info(f"#   A random seed {seed} has been used")
            rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
            rs_shuffle = rs.shuffle
        else:
            rs_shuffle = np.random.shuffle

        for k in sorted(chrnames):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            num = <uint32_t>round(self.__locations[k].shape[0] * percent, 5 )
            rs_shuffle( self.__locations[k] )
            self.__locations[k].resize( num, refcheck = False )
            self.__locations[k].sort( order = ['l', 'r'] ) # sort by leftmost positions
            self.__size[k] = self.__locations[k].shape[0]
            self.length += ( self.__locations[k]['r'] - self.__locations[k]['l'] ).sum()
            self.total += self.__size[k]
        self.average_template_length = <float32_t>( self.length )/ self.total
        return

    cpdef object sample_percent_copy (self, float32_t percent, int32_t seed = -1):
        """Sample the tags for a given percentage. Return a new PETrackI object

        """
        cdef:
            uint32_t num, i_chrom      # num: number of reads allowed on a certain chromosome
            bytes k
            set chrnames
            object ret_petrackI, rs, rs_shuffle
            np.ndarray l

        ret_petrackI = PETrackI( anno=self.annotation, buffer_size = self.buffer_size)
        chrnames = self.get_chr_names()

        if seed >= 0:
            info(f"# A random seed {seed} has been used in the sampling function")
            rs = np.random.default_rng(seed)
        else:
            rs = np.random.default_rng()

        rs_shuffle = rs.shuffle
        for k in sorted(chrnames): # chrnames need to be sorted otherwise we can't assure reproducibility
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            l = np.copy( self.__locations[k] )
            num = <uint32_t>round(l.shape[0] * percent, 5 )
            rs_shuffle( l )
            l.resize( num, refcheck = False )
            l.sort( order = ['l', 'r'] ) # sort by leftmost positions
            ret_petrackI.__locations[ k ] = l
            ret_petrackI.__size[ k ] = l.shape[0]
            ret_petrackI.length += ( l['r'] - l['l'] ).sum()
            ret_petrackI.total += ret_petrackI.__size[ k ]
        ret_petrackI.average_template_length = <float32_t>( ret_petrackI.length )/ ret_petrackI.total
        ret_petrackI.set_rlengths( self.get_rlengths() )
        return ret_petrackI

    cpdef void sample_num (self, uint64_t samplesize, int32_t seed = -1):
        """Sample the tags for a given number.

        Warning: the current object is changed!
        """
        cdef:
            float32_t percent
        percent = <float32_t>(samplesize)/self.total
        self.sample_percent ( percent, seed )
        return

    cpdef object sample_num_copy (self, uint64_t samplesize, int32_t seed = -1):
        """Sample the tags for a given number.

        Warning: the current object is changed!
        """
        cdef:
            float32_t percent
        percent = <float32_t>(samplesize)/self.total
        return self.sample_percent_copy ( percent, seed )

    cpdef void print_to_bed (self, fhd=None):
        """Output to BEDPE format files. If fhd is given, write to a
        file, otherwise, output to standard output.

        """
        cdef:
            int32_t i, i_chrom, s, e
            bytes k
            set chrnames


        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, io.IOBase)

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            locs = self.__locations[k]

            for i in range(locs.shape[0]):
                s, e = locs[ i ]
                fhd.write("%s\t%d\t%d\n" % (k.decode(), s, e))
        return

    cpdef list pileup_a_chromosome ( self, bytes chrom, list scale_factor_s, float32_t baseline_value = 0.0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.

        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float32_t scale_factor

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            scale_factor = scale_factor_s[i]

            tmp_pileup = quick_pileup ( np.sort(self.__locations[chrom]['l']), np.sort(self.__locations[chrom]['r']), scale_factor, baseline_value ) # Can't directly pass partial nparray there since that will mess up with pointer calculation.

            if prev_pileup:
                prev_pileup = over_two_pv_array ( prev_pileup, tmp_pileup, func="max" )
            else:
                prev_pileup = tmp_pileup

        return prev_pileup

    cpdef list pileup_a_chromosome_c ( self, bytes chrom, list ds, list scale_factor_s, float32_t baseline_value = 0.0 ):
        """pileup a certain chromosome, return [p,v] (end position and value) list.

        This function is for control track. Basically, here is a
        simplified function from FixWidthTrack. We pretend the PE is
        SE data and left read is on plus strand and right read is on
        minus strand.

        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float32_t scale_factor
            int64_t d, five_shift, three_shift
            int64_t rlength = self.get_rlengths()[chrom]

        if not self.__sorted: self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            five_shift = d//2
            three_shift= d//2

            tmp_pileup = se_all_in_one_pileup ( self.__locations[chrom]['l'], self.__locations[chrom]['r'], five_shift, three_shift, rlength, scale_factor, baseline_value )

            if prev_pileup:
                prev_pileup = over_two_pv_array ( prev_pileup, tmp_pileup, func="max" )
            else:
                prev_pileup = tmp_pileup

        return prev_pileup


    cpdef object pileup_bdg ( self, list scale_factor_s, float32_t baseline_value = 0.0 ):
        """pileup all chromosomes, and return a bedGraphTrackI object.

        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """
        cdef:
            list tmp_pileup, prev_pileup
            float32_t scale_factor
            bytes chrom
            object bdg
            int32_t prev_s

        #info(f"start to pileup")
        bdg = bedGraphTrackI( baseline_value = baseline_value )

        for chrom in sorted(self.get_chr_names()):
            prev_pileup = None
            for i in range(len(scale_factor_s)):
                scale_factor = scale_factor_s[i]

                tmp_pileup = quick_pileup ( np.sort(self.__locations[chrom]['l']), np.sort(self.__locations[chrom]['r']), scale_factor, baseline_value ) # Can't directly pass partial nparray there since that will mess up with pointer calculation.

                if prev_pileup:
                    prev_pileup = over_two_pv_array ( prev_pileup, tmp_pileup, func="max" )
                else:
                    prev_pileup = tmp_pileup
            # save to bedGraph
            bdg.add_chrom_data( chrom, pyarray('i', prev_pileup[0]), pyarray('f', prev_pileup[1]) )
        return bdg

    cpdef list pileup_bdg_hmmr ( self, list mapping, float32_t baseline_value = 0.0 ):
        """pileup all chromosomes, and return a list of four bedGraphTrackI objects: short, mono, di, and tri nucleosomal signals.

        The idea is that for each fragment length, we generate four bdg using four weights from four distributions. Then we add all sets of four bdgs together.

        Way to generate 'mapping', based on HMMR EM means and stddevs:
        fl_dict = petrack.count_fraglengths()
        fl_list = list(fl_dict.keys())
        fl_list.sort()
        weight_mapping = generate_weight_mapping( fl_list, em_means, em_stddevs )
        """
        cdef:
            list ret_pileup
            set chroms
            bytes chrom
            int i

        ret_pileup = []
        for i in range( len(mapping) ): ret_pileup.append( {} )
        chroms = self.get_chr_names()
        for i in range( len(mapping) ):
            for chrom in sorted(chroms):
                ret_pileup[ i ][ chrom ] = pileup_from_LR_hmmratac( self.__locations[ chrom ], mapping[ i ] )
        return ret_pileup

