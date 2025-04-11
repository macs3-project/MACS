# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-04-11 13:24:36 Tao Liu>

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
from array import array as pyarray
from collections import Counter

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.Pileup import (quick_pileup,
                                 over_two_pv_array,
                                 se_all_in_one_pileup)
from MACS3.Signal.BedGraph import (bedGraphTrackI,
                                   bedGraphTrackII)
from MACS3.Signal.PileupV2 import (pileup_from_LR_hmmratac,
                                   pileup_from_LRC)
from MACS3.Signal.Region import Regions
# ------------------------------------
# Other modules
# ------------------------------------
import cython
import numpy as np
import cython.cimports.numpy as cnp
from cython.cimports.cpython import bool
from cython.cimports.libc.stdint import INT32_MAX as INT_MAX

from MACS3.Utilities.Logger import logging

logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

# Let numpy enforce PE-ness using ndarray, gives bonus speedup when sorting
# PE data doesn't have strandedness


@cython.cclass
class PETrackI:
    """Paired End Locations Track class I along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    locations = cython.declare(dict, visibility="public")
    size = cython.declare(dict, visibility="public")
    buf_size = cython.declare(dict, visibility="public")
    is_sorted = cython.declare(bool, visibility="public")
    total = cython.declare(cython.ulong, visibility="public")
    annotation = cython.declare(str, visibility="public")
    # rlengths: reference chromosome lengths dictionary
    rlengths = cython.declare(dict, visibility="public")
    buffer_size = cython.declare(cython.long, visibility="public")
    length = cython.declare(cython.long, visibility="public")
    average_template_length = cython.declare(cython.float, visibility="public")
    is_destroyed: bool

    def __init__(self, anno: str = "", buffer_size: cython.long = 100000):
        """fw is the fixed-width for all locations.

        """
        # dictionary with chrname as key, nparray with
        # [('l','i4'),('r','i4')] as value
        self.locations = {}
        # dictionary with chrname as key, size of the above nparray as value
        # size is to remember the number of the fragments added to this chromosome
        self.size = {}
        # dictionary with chrname as key, size of the above nparray as value
        self.buf_size = {}
        self.is_sorted = False
        self.total = 0           # total fragments
        self.annotation = anno   # need to be figured out
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.average_template_length = 0.0
        self.is_destroyed = False

    @cython.ccall
    def add_loc(self, chromosome: bytes,
                start: cython.int, end: cython.int):
        """Add a location to the list according to the sequence name.

        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        """
        i: cython.int

        if chromosome not in self.locations:
            self.buf_size[chromosome] = self.buffer_size
            # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            self.locations[chromosome] = np.zeros(shape=self.buffer_size,
                                                  dtype=[('l', 'i4'), ('r', 'i4')])
            self.locations[chromosome][0] = (start, end)
            self.size[chromosome] = 1
        else:
            i = self.size[chromosome]
            if self.buf_size[chromosome] == i:
                self.buf_size[chromosome] += self.buffer_size
                self.locations[chromosome].resize((self.buf_size[chromosome]),
                                                  refcheck=False)
            self.locations[chromosome][i] = (start, end)
            self.size[chromosome] = i + 1
        self.length += end - start
        return

    @cython.ccall
    def destroy(self):
        """Destroy this object and release mem.
        """
        chrs: set
        chromosome: bytes

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.locations:
                self.locations[chromosome].resize(self.buffer_size,
                                                  refcheck=False)
                self.locations[chromosome].resize(0,
                                                  refcheck=False)
                self.locations[chromosome] = None
                self.locations.pop(chromosome)
        self.is_destroyed = True
        return

    @cython.ccall
    def set_rlengths(self, rlengths: dict) -> bool:
        """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this petrack object will be updated.

        If a chromosome in this petrack is not covered by given
        rlengths, and it has no associated length, it will be set as
        maximum integer.
        """
        valid_chroms: set
        missed_chroms: set
        chrom: bytes

        valid_chroms = set(self.locations.keys()).intersection(rlengths.keys())
        for chrom in sorted(valid_chroms):
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.locations.keys()).difference(rlengths.keys())
        for chrom in sorted(missed_chroms):
            self.rlengths[chrom] = INT_MAX
        return True

    @cython.ccall
    def get_rlengths(self) -> dict:
        """Get reference chromosome lengths dictionary.

        If self.rlengths is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

    @cython.ccall
    def finalize(self):
        """Resize np arrays for 5' positions and sort them in place

        Note: If this function is called, it's impossible to append
        more files to this PETrackI object. So remember to call it
        after all the files are read!

        """
        c: bytes
        chrnames: set

        self.total = 0

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c].resize((self.size[c]), refcheck=False)
            self.locations[c].sort(order=['l', 'r'])
            self.total += self.size[c]

        self.is_sorted = True
        self.average_template_length = cython.cast(cython.float, self.length) / self.total
        return

    @cython.ccall
    def get_locations_by_chr(self, chromosome: bytes):
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if chromosome in self.locations:
            return self.locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return all the chromosome names in this track object as a python set.
        """
        return set(self.locations.keys())

    @cython.ccall
    def sort(self):
        """Naive sorting for locations.

        """
        c: bytes
        chrnames: set

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c].sort(order=['l', 'r'])  # sort by the leftmost location
        self.is_sorted = True
        return

    @cython.ccall
    def count_fraglengths(self) -> dict:
        """Return a dictionary of the counts for sizes/fragment
        lengths of each pair.

        This function is for HMMRATAC.

        """
        sizes: cnp.ndarray(cnp.int32_t, ndim=1)
        s: cython.int
        locs: cnp.ndarray
        chrnames: list
        i: cython.int

        counter = Counter()
        chrnames = list(self.get_chr_names())
        for i in range(len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = locs['r'] - locs['l']
            for s in sizes:
                counter[s] += 1
        return dict(counter)

    @cython.ccall
    def fraglengths(self) -> cnp.ndarray:
        """Return the sizes/fragment lengths of each pair.

        This function is for HMMRATAC EM training.
        """
        sizes: cnp.ndarray(np.int32_t, ndim=1)
        locs: cnp.ndarray
        chrnames: list
        i: cython.int

        chrnames = list(self.get_chr_names())
        locs = self.locations[chrnames[0]]
        sizes = locs['r'] - locs['l']
        for i in range(1, len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = np.concatenate((sizes, locs['r'] - locs['l']))
        return sizes

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def exclude(self, regions):
        """Exclude reads overlapping with anything in the regions.

        Run it right after you add all data into this object.
        """
        i: cython.ulong
        k: bytes
        locs: cnp.ndarray
        locs_size: cython.ulong
        chrnames: set
        regions_c: list
        selected_idx: cnp.ndarray
        regions_chrs: list
        r1: cnp.void            # this is the location in numpy.void -- like a tuple
        r2: tuple               # this is the region
        n_rl1: cython.long
        n_rl2: cython.long

        if not self.is_sorted:
            self.sort()

        assert isinstance(regions, Regions)
        regions.sort()
        regions_chrs = list(regions.regions.keys())

        chrnames = self.get_chr_names()

        for k in chrnames:      # for each chromosome
            locs = self.locations[k]
            locs_size = self.size[k]
            # let's check if k is in regions_chr
            if k not in regions_chrs:
                # do nothing and continue
                self.total += locs_size
                continue

            # discard overlapping reads and make a new locations[k]
            # initialize boolean array as all TRUE, or all being kept
            selected_idx = np.ones(locs_size, dtype=bool)

            regions_c = regions.regions[k]

            i = 0
            n_rl1 = locs_size   # the number of locations left
            n_rl2 = len(regions_c)  # the number of regions left
            rl1_k = iter(locs).__next__
            rl2_k = iter(regions_c).__next__
            r1 = rl1_k()        # take the first value
            n_rl1 -= 1          # remaining rl1
            r2 = rl2_k()
            n_rl2 -= 1          # remaining rl2
            while (True):
                # we do this until there is no r1 or r2 left.
                if r2[0] < r1[1] and r1[0] < r2[1]:
                    # since we found an overlap, r1 will be skipped/excluded
                    # and move to the next r1
                    # get rid of this one
                    n_rl1 -= 1
                    self.length -= r1[1] - r1[0]
                    selected_idx[i] = False

                    if n_rl1 > 0:
                        r1 = rl1_k()  # take the next location
                        i += 1
                        continue
                    else:
                        break
                if r1[1] < r2[1]:
                    # in this case, we need to move to the next r1,
                    n_rl1 -= 1
                    if n_rl1 > 0:
                        r1 = rl1_k()  # take the next location
                        i += 1
                    else:
                        # no more r1 left
                        break
                else:
                    # in this case, we need to move the next r2
                    if n_rl2:
                        r2 = rl2_k()  # take the next region
                        n_rl2 -= 1
                    else:
                        # no more r2 left
                        break

            self.locations[k] = locs[selected_idx]
            self.size[k] = self.locations[k].shape[0]
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            selected_idx.resize(self.buffer_size, refcheck=False)
            selected_idx.resize(0, refcheck=False)
        self.finalize()         # use length to set average_template_length and total
        return

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def filter_dup(self, maxnum: cython.int = -1):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        n: cython.int
        loc_start: cython.int
        loc_end: cython.int
        current_loc_start: cython.int
        current_loc_end: cython.int
        i: cython.ulong
        locs_size: cython.ulong
        k: bytes
        locs: cnp.ndarray
        chrnames: set
        selected_idx: cnp.ndarray

        if maxnum < 0:
            return              # condition to return if not filtering

        if not self.is_sorted:
            self.sort()

        self.total = 0
        # self.length = 0
        self.average_template_length = 0.0

        chrnames = self.get_chr_names()

        for k in chrnames:      # for each chromosome
            locs = self.locations[k]
            locs_size = self.size[k]
            if locs_size == 1:
                # do nothing and continue
                self.total += locs_size
                continue
            # discard duplicate reads and make a new locations[k]
            # initialize boolean array as all TRUE, or all being kept
            selected_idx = np.ones(locs_size, dtype=bool)
            # get the first loc
            (current_loc_start, current_loc_end) = locs[0]
            i = 1               # index of new_locs
            n = 1  # the number of tags in the current genomic location
            for i in range(1, locs_size):
                (loc_start, loc_end) = locs[i]
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
                        selected_idx[i] = False
                        # subtract current_loc_l from self.length
                        self.length -= current_loc_end - current_loc_start
            self.locations[k] = locs[selected_idx]
            self.size[k] = self.locations[k].shape[0]
            self.total += self.size[k]
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            selected_idx.resize(self.buffer_size, refcheck=False)
            selected_idx.resize(0, refcheck=False)
        self.average_template_length = self.length / self.total
        return

    @cython.ccall
    def sample_percent(self,
                       percent: cython.float,
                       seed: cython.int = -1):
        """Sample the tags for a given percentage.

        Note that the sampling is on each chromosome using the given
        percentage.

        Warning: the current object is changed! If a new PETrackI is
        wanted, use sample_percent_copy instead.

        """
        # num: number of reads allowed on a certain chromosome
        num: cython.uint
        k: bytes
        chrnames: set

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

            num = cython.cast(cython.uint,
                              round(self.locations[k].shape[0] * percent, 5))
            rs_shuffle(self.locations[k])
            self.locations[k].resize(num, refcheck=False)
            self.locations[k].sort(order=['l', 'r'])  # sort by leftmost positions
            self.size[k] = self.locations[k].shape[0]
            self.length += (self.locations[k]['r'] - self.locations[k]['l']).sum()
            self.total += self.size[k]
        self.average_template_length = cython.cast(cython.float, self.length)/self.total
        return

    @cython.ccall
    def sample_percent_copy(self,
                            percent: cython.float,
                            seed: cython.int = -1):
        """Sample the tags for a given percentage. Return a new
        PETrackI object

        Note that the sampling is on each chromosome using the given
        percentage.

        """
        # num: number of reads allowed on a certain chromosome
        num: cython.uint
        k: bytes
        chrnames: set
        ret_petrackI: PETrackI
        loc: cnp.ndarray

        ret_petrackI = PETrackI(anno=self.annotation,
                                buffer_size=self.buffer_size)
        chrnames = self.get_chr_names()

        if seed >= 0:
            info(f"# A random seed {seed} has been used in the sampling function")
            rs = np.random.default_rng(seed)
        else:
            rs = np.random.default_rng()

        rs_shuffle = rs.shuffle

        # chrnames need to be sorted otherwise we can't assure reproducibility
        for k in sorted(chrnames):
            # for each chromosome.
            # This loop body is too big, I may need to split code later...
            loc = np.copy(self.locations[k])
            num = cython.cast(cython.uint, round(loc.shape[0] * percent, 5))
            rs_shuffle(loc)
            loc.resize(num, refcheck=False)
            loc.sort(order=['l', 'r'])  # sort by leftmost positions
            ret_petrackI.locations[k] = loc
            ret_petrackI.size[k] = loc.shape[0]
            ret_petrackI.length += (loc['r'] - loc['l']).sum()
            ret_petrackI.total += ret_petrackI.size[k]
        ret_petrackI.average_template_length = cython.cast(cython.float, ret_petrackI.length)/ret_petrackI.total
        ret_petrackI.set_rlengths(self.get_rlengths())
        return ret_petrackI

    @cython.ccall
    def sample_num(self,
                   samplesize: cython.ulong,
                   seed: cython.int = -1):
        """Sample the tags for a given number.

        Note that the sampling is on each chromosome using the the
        percentage calculated from the given number.

        Warning: the current object is changed!
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize)/self.total
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def sample_num_copy(self,
                        samplesize: cython.ulong,
                        seed: cython.int = -1):
        """Sample the tags for a given number.

        Note that the sampling is on each chromosome using the the
        percentage calculated from the given number.

        Warning: the current object is changed!
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize)/self.total
        return self.sample_percent_copy(percent, seed)

    @cython.ccall
    def print_to_bed(self, fhd=None):
        """Output to BEDPE format files. If fhd is given, write to a
        file, otherwise, output to standard output.

        """
        i: cython.int
        s: cython.int
        e: cython.int
        k: bytes
        chrnames: set

        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, io.IOBase)

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            locs = self.locations[k]

            for i in range(locs.shape[0]):
                s, e = locs[i]
                fhd.write("%s\t%d\t%d\n" % (k.decode(), s, e))
        return

    @cython.ccall
    def pileup_a_chromosome(self,
                            chrom: bytes,
                            scale_factor: cython.float = 1.0,
                            baseline_value: cython.float = 0.0) -> list:
        """pileup a certain chromosome, return [p,v] (end position and
        value) list.

        scale_factor : linearly scale the pileup value.

        baseline_value : a value to be filled for missing values, and
                         will be the minimum pileup.

        """
        tmp_pileup: list

        tmp_pileup = quick_pileup(np.sort(self.locations[chrom]['l']),
                                  np.sort(self.locations[chrom]['r']),
                                  scale_factor, baseline_value)
        return tmp_pileup

    @cython.ccall
    def pileup_a_chromosome_c(self,
                              chrom: bytes,
                              ds: list,
                              scale_factor_s: list,
                              baseline_value: cython.float = 0.0) -> list:
        """pileup a certain chromosome, return [p,v] (end position and
        value) list.

        This function is for control track. Basically, here is a
        simplified function from FixWidthTrack. We pretend the PE is
        SE data and left read is on plus strand and right read is on
        minus strand.

        ds : tag will be extended to this value to 3' direction,
             unless directional is False. Can contain multiple
             extension values. Final pileup will the maximum.
        scale_factor_s : linearly scale the pileup value applied to
                         each d in ds. The list should have the same
                         length as ds.
        baseline_value : a value to be filled for missing values, and
                         will be the minimum pileup.
        """
        tmp_pileup: list
        prev_pileup: list
        scale_factor: cython.float
        d: cython.long
        five_shift: cython.long
        three_shift: cython.long
        rlength: cython.long = self.get_rlengths()[chrom]

        if not self.is_sorted:
            self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            five_shift = d//2
            three_shift = d//2

            tmp_pileup = se_all_in_one_pileup(self.locations[chrom]['l'],
                                              self.locations[chrom]['r'],
                                              five_shift,
                                              three_shift,
                                              rlength,
                                              scale_factor,
                                              baseline_value)

            if prev_pileup:
                prev_pileup = over_two_pv_array(prev_pileup,
                                                tmp_pileup,
                                                func="max")
            else:
                prev_pileup = tmp_pileup

        with open("tmp_b_c.txt", "w") as f:
            for i in range(len(prev_pileup[0])):
                f.write(f"{prev_pileup[0][i]}\t{prev_pileup[1][i]}\n")

        return prev_pileup

    @cython.ccall
    def pileup_bdg(self,
                   scale_factor: cython.float = 1.0,
                   baseline_value: cython.float = 0.0):
        """pileup all chromosomes, and return a bedGraphTrackI object.

        scale_factor : a value to scale the pileup values.
        baseline_value : a value to be filled for missing values, and
                         will be the minimum pileup.

        """
        tmp_pileup: list
        chrom: bytes
        bdg: bedGraphTrackI

        bdg = bedGraphTrackI(baseline_value=baseline_value)

        for chrom in sorted(self.get_chr_names()):
            tmp_pileup = quick_pileup(np.sort(self.locations[chrom]['l']),
                                      np.sort(self.locations[chrom]['r']),
                                      scale_factor,
                                      baseline_value)

            # save to bedGraph
            bdg.add_chrom_data(chrom,
                               pyarray('i', tmp_pileup[0]),
                               pyarray('f', tmp_pileup[1]))
        return bdg

    @cython.ccall
    def pileup_bdg_hmmr(self,
                        mapping: list,
                        baseline_value: cython.float = 0.0) -> list:
        """pileup all chromosomes, and return a list of four p-v
        ndarray objects: short, mono, di, and tri nucleosomal signals.

        This is specifically designed for hmmratac
        HMM_SignalProcessing.py. Not a general function.

        The idea is that for each fragment length, we generate four
        bdg using four weights from four distributions. Then we add
        all sets of four bdgs together.

        Way to generate 'mapping', based on HMMR EM means and stddevs:
        fl_dict = petrack.count_fraglengths()
        fl_list = list(fl_dict.keys())
        fl_list.sort()
        weight_mapping = generate_weight_mapping(fl_list, em_means, em_stddevs)

        """
        ret_pileup: list
        chroms: set
        chrom: bytes
        i: cython.int

        ret_pileup = []
        for i in range(len(mapping)):
            ret_pileup.append({})
        chroms = self.get_chr_names()
        for i in range(len(mapping)):
            for chrom in sorted(chroms):
                ret_pileup[i][chrom] = pileup_from_LR_hmmratac(self.locations[chrom], mapping[i])
        return ret_pileup


@cython.cclass
class PETrackII:
    """Paired-end track class for fragment files from single-cell
    ATAC-seq experiments. We will store data of start, end, barcode,
    and count from the fragment files.

    * I choose not to inherit PETrackI because there would be a lot of
      differences.

    """
    locations = cython.declare(dict, visibility="public")
    # add another dict for storing barcode for each fragment we will
    # first convert barcode into integer and remember them in the
    # barcode_dict, which will map the rule to numerize
    # key:bytes as value:4bytes_integer
    barcodes = cython.declare(dict, visibility="public")
    barcode_dict = cython.declare(dict, visibility="public")
    # the last number for barcodes, used to map barcode to integer
    barcode_last_n: cython.int

    size = cython.declare(dict, visibility="public")
    buf_size = cython.declare(dict, visibility="public")
    is_sorted = cython.declare(bool, visibility="public")
    total = cython.declare(cython.ulong, visibility="public")
    annotation = cython.declare(str, visibility="public")
    # rlengths: reference chromosome lengths dictionary
    rlengths = cython.declare(dict, visibility="public")
    buffer_size = cython.declare(cython.long, visibility="public")
    length = cython.declare(cython.long, visibility="public")
    average_template_length = cython.declare(cython.float, visibility="public")
    is_destroyed: bool

    def __init__(self, anno: str = "", buffer_size: cython.long = 100000):
        # dictionary with chrname as key, nparray with
        # [('l','i4'),('r','i4'),('c','u2')] as value
        self.locations = {}
        # dictionary with chrname as key, size of the above nparray as value
        # size is to remember the size of the fragments added to this chromosome
        self.size = {}
        # dictionary with chrname as key, size of the above nparray as value
        self.buf_size = {}
        self.is_sorted = False
        self.total = 0           # total fragments
        self.annotation = anno   # need to be figured out
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.average_template_length = 0.0
        self.is_destroyed = False

        self.barcodes = {}
        self.barcode_dict = {}
        self.barcode_last_n = 0

    @cython.ccall
    def add_loc(self,
                chromosome: bytes,
                start: cython.int,
                end: cython.int,
                barcode: bytes,
                count: cython.ushort):
        """Add a location to the list according to the sequence name.

        chromosome: mostly the chromosome name
        start: left position of the fragment
        end: right position of the fragment
        barcode: the barcode of the fragment
        count: the count of the fragment
        """
        i: cython.int
        # bn: the integer in barcode_dict for this barcode
        bn: cython.int

        if barcode not in self.barcode_dict:
            self.barcode_dict[barcode] = self.barcode_last_n
            self.barcode_last_n += 1
        bn = self.barcode_dict[barcode]

        if chromosome not in self.locations:
            self.buf_size[chromosome] = self.buffer_size
            # note: ['l'] is the leftmost end, ['r'] is the rightmost end of fragment.
            # ['c'] is the count number of this fragment
            self.locations[chromosome] = np.zeros(shape=self.buffer_size,
                                                  dtype=[('l', 'i4'), ('r', 'i4'), ('c', 'u2')])
            self.barcodes[chromosome] = np.zeros(shape=self.buffer_size,
                                                 dtype='i4')
            self.locations[chromosome][0] = (start, end, count)
            self.barcodes[chromosome][0] = bn
            self.size[chromosome] = 1
        else:
            i = self.size[chromosome]
            if self.buf_size[chromosome] == i:
                self.buf_size[chromosome] += self.buffer_size
                self.locations[chromosome].resize((self.buf_size[chromosome]),
                                                  refcheck=False)
                self.barcodes[chromosome].resize((self.buf_size[chromosome]),
                                                 refcheck=False)                
            self.locations[chromosome][i] = (start, end, count)
            self.barcodes[chromosome][i] = bn
            self.size[chromosome] = i + 1
        self.length += (end - start) * count
        return

    @cython.ccall
    def destroy(self):
        """Destroy this object and release mem.
        """
        chrs: set
        chromosome: bytes

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.locations:
                self.locations[chromosome].resize(self.buffer_size,
                                                  refcheck=False)
                self.locations[chromosome].resize(0,
                                                  refcheck=False)
                self.locations[chromosome] = None
                self.locations.pop(chromosome)
                self.barcodes.resize(self.buffer_size,
                                     refcheck=False)
                self.barcodes.resize(0,
                                     refcheck=False)
                self.barcodes[chromosome] = None
                self.barcodes.pop(chromosome)
        self.barcode_dict = {}
        self.is_destroyed = True
        return

    @cython.ccall
    def set_rlengths(self, rlengths: dict) -> bool:
        """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this petrack object will be updated.

        If a chromosome in this petrack is not covered by given
        rlengths, and it has no associated length, it will be set as
        maximum integer.
        """
        valid_chroms: set
        missed_chroms: set
        chrom: bytes

        valid_chroms = set(self.locations.keys()).intersection(rlengths.keys())
        for chrom in sorted(valid_chroms):
            self.rlengths[chrom] = rlengths[chrom]
        missed_chroms = set(self.locations.keys()).difference(rlengths.keys())
        for chrom in sorted(missed_chroms):
            self.rlengths[chrom] = INT_MAX
        return True

    @cython.ccall
    def get_rlengths(self) -> dict:
        """Get reference chromosome lengths dictionary.

        If self.rlengths is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

    @cython.ccall
    def finalize(self):
        """Resize np arrays for 5' positions and sort them in place

        Note: If this function is called, it's impossible to append
        more files to this PETrackII object. So remember to call it
        after all the files are read!

        """
        c: bytes
        chrnames: set
        indices: cnp.ndarray

        self.total = 0

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c].resize((self.size[c]), refcheck=False)
            indices = np.argsort(self.locations[c], order=['l', 'r'])
            self.locations[c] = self.locations[c][indices]
            self.barcodes[c] = self.barcodes[c][indices]
            self.total += np.sum(self.locations[c]['c'])  # self.size[c]

        assert self.total > 0, "Error: no fragments in PETrackII"

        self.is_sorted = True
        self.average_template_length = cython.cast(cython.float,
                                                   self.length) / self.total
        return

    @cython.ccall
    def get_locations_by_chr(self, chromosome: bytes):
        """Return a np array of left/right/count for certain chromosome.

        """
        if chromosome in self.locations:
            return self.locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return all the chromosome names in this track object as a
        python set.

        """
        return set(self.locations.keys())

    @cython.ccall
    def sort(self):
        """Naive sorting for locations.

        """
        c: bytes
        chrnames: set
        indices: cnp.ndarray

        chrnames = self.get_chr_names()

        for c in chrnames:
            indices = np.argsort(self.locations[c], order=['l', 'r'])
            self.locations[c] = self.locations[c][indices]
            self.barcodes[c] = self.barcodes[c][indices]
        self.is_sorted = True
        return

    @cython.ccall
    def count_fraglengths(self) -> dict:
        """Return a dictionary of the counts for sizes/fragment
        lengths of each pair.

        This function is for HMMRATAC.

        """
        sizes: cnp.ndarray(cnp.int32_t, ndim=1)
        s: cython.int
        locs: cnp.ndarray
        chrnames: list
        i: cython.int

        counter = Counter()
        chrnames = list(self.get_chr_names())
        for i in range(len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = locs['r'] - locs['l']
            for s in sizes:
                counter[s] += locs['c']
        return dict(counter)

    @cython.ccall
    def fraglengths(self) -> cnp.ndarray:
        """Return the sizes/fragment lengths of each pair.

        This function is for HMMRATAC EM training.
        """
        sizes: cnp.ndarray(np.int32_t, ndim=1)
        t_sizes: cnp.ndarray(np.int32_t, ndim=1)
        locs: cnp.ndarray
        chrnames: list
        i: cython.int

        chrnames = list(self.get_chr_names())
        locs = self.locations[chrnames[0]]
        sizes = locs['r'] - locs['l']
        sizes = [x for x, count in zip(sizes, locs['c']) for _ in range(count)]

        for i in range(1, len(chrnames)):
            locs = self.locations[chrnames[i]]
            t_sizes = locs['r'] - locs['l']
            t_sizes = [x for x, count in zip(t_sizes, locs['c']) for _ in range(count)]
            sizes = np.concatenate((sizes, t_sizes))
        return sizes

    @cython.ccall
    def subset(self, selected_barcodes: set):
        """Make a subset of PETrackII with only the given barcodes.

        Note: the selected_barcodes is a set of barcodes in python
        bytes. For example, {b"ATCTGCTAGTCTACAT", b"ATTCTCGATGCAGTCA"}

        """
        indices: cnp.ndarray
        chrs: set
        selected_barcodes_filtered: list
        selected_barcodes_n: list
        chromosome: bytes
        ret: PETrackII

        ret = PETrackII()
        chrs = self.get_chr_names()

        # first we need to convert barcodes into integers in our
        # barcode_dict
        selected_barcodes_filtered = [b
                                      for b in selected_barcodes
                                      if b in self.barcode_dict]
        ret.barcode_dict = {b: self.barcode_dict[b]
                            for b in selected_barcodes_filtered}
        selected_barcodes_n = [self.barcode_dict[b]
                               for b in selected_barcodes_filtered]
        ret.barcode_last_n = self.barcode_last_n

        # pass some values from self to ret
        ret.annotation = self.annotation
        ret.is_sorted = self.is_sorted
        ret.rlengths = self.rlengths
        ret.buffer_size = self.buffer_size
        # ret.total = 0
        ret.length = 0
        # ret.average_template_length = 0
        ret.is_destroyed = True

        for chromosome in sorted(chrs):
            indices = np.where(np.isin(self.barcodes[chromosome],
                                       list(selected_barcodes_n)))[0]
            ret.barcodes[chromosome] = self.barcodes[chromosome][indices]
            ret.locations[chromosome] = self.locations[chromosome][indices]
            ret.size[chromosome] = len(ret.locations[chromosome])
            ret.buf_size[chromosome] = ret.size[chromosome]
            # ret.total += np.sum(ret.locations[chromosome]['c'])
            ret.length += np.sum((ret.locations[chromosome]['r'] -
                                  ret.locations[chromosome]['l']) *
                                 ret.locations[chromosome]['c'])
        ret.finalize()
        # ret.average_template_length = ret.length / ret.total
        return ret

    @cython.ccall
    def pileup_a_chromosome(self,
                            chrom: bytes,
                            scale_factor: cython.float = 1.0,
                            baseline_value: cython.float = 0.0) -> list:
        """pileup a certain chromosome, return p-v ndarray (end
        position and pileup value).
        """
        pv: cnp.ndarray
        v: cnp.ndarray

        pv = pileup_from_LRC(self.locations[chrom])
        v = pv['v']
        v = v * scale_factor
        v[v < baseline_value] = baseline_value

        return [pv['p'], v]

    @cython.ccall
    def pileup_a_chromosome_c(self,
                              chrom: bytes,
                              ds: list,
                              scale_factor_s: list,
                              baseline_value: cython.float = 0.0) -> list:
        """We will mimic the PETrackI way to pileup PE data as SE to
        build the control signal track.

        """
        prev_pileup: list
        scale_factor: cython.float
        d: cython.long
        five_shift: cython.long

        pv: cnp.ndarray
        v: cnp.ndarray
        tmp_arr_l: cnp.ndarray
        tmp_arr_r: cnp.ndarray
        tmp_arr: cnp.ndarray

        ####
        if not self.is_sorted:
            self.sort()

        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        prev_pileup = None

        for i in range(len(scale_factor_s)):
            d = ds[i]
            scale_factor = scale_factor_s[i]
            five_shift = d//2

            # note, we have to pileup left ends and right ends separately
            tmp_arr_l = self.locations[chrom].copy()
            tmp_arr_l['l'] = tmp_arr_l['l'] - five_shift
            tmp_arr_l['r'] = tmp_arr_l['l'] + d

            tmp_arr_r = self.locations[chrom].copy()
            tmp_arr_r['l'] = tmp_arr_r['r'] - five_shift
            tmp_arr_r['r'] = tmp_arr_r['l'] + d

            tmp_arr = np.concatenate([tmp_arr_l, tmp_arr_r])
            del tmp_arr_l
            del tmp_arr_r

            pv = pileup_from_LRC(tmp_arr)

            v = pv['v']
            v = v * scale_factor
            v[v < baseline_value] = baseline_value

            if prev_pileup:
                prev_pileup = over_two_pv_array(prev_pileup,
                                                [pv['p'], v],
                                                func="max")
            else:
                prev_pileup = [pv['p'], v]

        return prev_pileup

    @cython.ccall
    def pileup_bdg(self):
        """Pileup all chromosome and return a bdg object.
        """
        bdg: bedGraphTrackI
        pv: cnp.ndarray

        bdg = bedGraphTrackI()
        for chrom in self.get_chr_names():
            pv = pileup_from_LRC(self.locations[chrom])
            bdg.add_chrom_data_PV(chrom, pv)
        return bdg

    @cython.ccall
    def pileup_bdg2(self):
        """Pileup all chromosome and return a bdg object.
        """
        bdg: bedGraphTrackII
        pv: cnp.ndarray

        bdg = bedGraphTrackII()
        for chrom in self.get_chr_names():
            pv = pileup_from_LRC(self.locations[chrom])
            bdg.add_chrom_data(chrom, pv)
        # bedGraphTrackII needs to be 'finalized'.
        bdg.finalize()
        return bdg

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def exclude(self, regions):
        """Exclude reads overlapping with anything in the regions.

        Run it right after you add all data into this object.
        """
        i: cython.ulong
        k: bytes
        locs: cnp.ndarray
        locs_size: cython.ulong
        chrnames: set
        regions_c: list
        selected_idx: cnp.ndarray
        regions_chrs: list
        r1: cnp.void
        r2: tuple
        n_rl1: cython.long
        n_rl2: cython.long

        if not self.is_sorted:
            self.sort()

        assert isinstance(regions, Regions)
        regions.sort()
        regions_chrs = list(regions.regions.keys())

        chrnames = self.get_chr_names()

        for k in chrnames:      # for each chromosome
            locs = self.locations[k]
            locs_size = self.size[k]
            # let's check if k is in regions_chr
            if k not in regions_chrs:
                # do nothing and continue
                continue

            # discard overlapping reads and make a new locations[k]
            # initialize boolean array as all TRUE, or all being kept
            selected_idx = np.ones(locs_size, dtype=bool)

            regions_c = regions.regions[k]

            i = 0
            n_rl1 = len(locs)
            n_rl2 = len(regions_c)
            rl1_k = iter(locs).__next__
            rl2_k = iter(regions_c).__next__
            r1 = rl1_k()
            n_rl1 -= 1          # remaining rl1
            r2 = rl2_k()
            n_rl2 -= 1          # remaining rl2
            while (True):
                # we do this until there is no r1 or r2 left.
                if r2[0] < r1[1] and r1[0] < r2[1]:
                    # since we found an overlap, r1 will be skipped/excluded
                    # and move to the next r1
                    # get rid of this one
                    n_rl1 -= 1
                    self.length -= (r1[1] - r1[0])*r1[2]
                    selected_idx[i] = False

                    if n_rl1 >= 0:
                        r1 = rl1_k()
                        i += 1
                        continue
                    else:
                        break
                if r1[1] < r2[1]:
                    # in this case, we need to move to the next r1,
                    n_rl1 -= 1
                    if n_rl1 >= 0:
                        r1 = rl1_k()
                        i += 1
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

            self.locations[k] = locs[selected_idx]
            self.barcodes[k] = self.barcodes[k][selected_idx]
            self.size[k] = self.locations[k].shape[0]
            # free memory?
            # I know I should shrink it to 0 size directly,
            # however, on Mac OSX, it seems directly assigning 0
            # doesn't do a thing.
            selected_idx.resize(self.buffer_size, refcheck=False)
            selected_idx.resize(0, refcheck=False)
        self.finalize()
        return
