# cython: language_level=3
# cython: profile=True
# Time-stamp: <2024-10-14 14:53:06 Tao Liu>

"""Module for FWTrack classes.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import io

# ------------------------------------
# MACS3 modules
# ------------------------------------

from MACS3.IO.PeakIO import PeakIO
from MACS3.Signal.Pileup import se_all_in_one_pileup, over_two_pv_array

# ------------------------------------
# Other modules
# ------------------------------------
import cython
import numpy as np
from cython.cimports.cpython import bool
import cython.cimports.numpy as cnp
from cython.cimports.libc.stdint import INT32_MAX as INT_MAX

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------


@cython.cclass
class FWTrack:
    """Fixed Width Locations Track class  along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    locations: dict
    pointer: dict
    buf_size: dict
    rlengths: dict
    is_sorted: bool
    is_destroyed: bool
    total = cython.declare(cython.ulong, visibility="public")
    annotation = cython.declare(str, visibility="public")
    buffer_size = cython.declare(cython.long, visibility="public")
    length = cython.declare(cython.long, visibility="public")
    fw = cython.declare(cython.int, visibility="public")

    def __init__(self,
                 fw: cython.int = 0,
                 anno: str = "",
                 buffer_size: cython.long = 100000):
        """fw is the fixed-width for all locations.

        """
        self.fw = fw
        self.locations = {}    # location pairs: two strands
        self.pointer = {}      # location pairs
        self.buf_size = {}     # location pairs
        self.is_sorted = False
        self.total = 0           # total tags
        self.annotation = anno   # need to be figured out
        # lengths of reference sequences, e.g. each chromosome in a genome
        self.rlengths = {}
        self.buffer_size = buffer_size
        self.length = 0
        self.is_destroyed = False

    @cython.ccall
    def destroy(self):
        """Destroy this object and release mem.
        """
        chrs: set
        chromosome: bytes

        chrs = self.get_chr_names()
        for chromosome in sorted(chrs):
            if chromosome in self.locations:
                self.locations[chromosome][0].resize(self.buffer_size,
                                                     refcheck=False)
                self.locations[chromosome][0].resize(0,
                                                     refcheck=False)
                self.locations[chromosome][1].resize(self.buffer_size,
                                                     refcheck=False)
                self.locations[chromosome][1].resize(0,
                                                     refcheck=False)
                self.locations[chromosome] = [None, None]
                self.locations.pop(chromosome)
        self.is_destroyed = True
        return

    @cython.ccall
    def add_loc(self,
                chromosome: bytes,
                fiveendpos: cython.int,
                strand: cython.int):
        """Add a location to the list according to the sequence name.

        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for minus strand
        strand     -- 0: plus, 1: minus
        """
        i: cython.int
        b: cython.int
        arr: cnp.ndarray

        if chromosome not in self.locations:
            self.buf_size[chromosome] = [self.buffer_size, self.buffer_size]
            self.locations[chromosome] = [np.zeros(self.buffer_size, dtype='i4'),
                                          np.zeros(self.buffer_size, dtype='i4')]
            self.pointer[chromosome] = [0, 0]
            self.locations[chromosome][strand][0] = fiveendpos
            self.pointer[chromosome][strand] = 1
        else:
            i = self.pointer[chromosome][strand]
            b = self.buf_size[chromosome][strand]
            arr = self.locations[chromosome][strand]
            if b == i:
                b += self.buffer_size
                arr.resize(b, refcheck=False)
                self.buf_size[chromosome][strand] = b
            arr[i] = fiveendpos
            self.pointer[chromosome][strand] += 1
        return

    @cython.ccall
    def finalize(self):
        """ Resize np arrays for 5' positions and sort them in place

        Note: If this function is called, it's impossible to append more files to this FWTrack object. So remember to call it after all the files are read!
        """
        c: bytes
        chrnames: set

        self.total = 0

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c][0].resize(self.pointer[c][0], refcheck=False)
            self.locations[c][0].sort()
            self.locations[c][1].resize(self.pointer[c][1], refcheck=False)
            self.locations[c][1].sort()
            self.total += self.locations[c][0].size + self.locations[c][1].size

        self.is_sorted = True
        self.length = self.fw * self.total
        return

    @cython.ccall
    def set_rlengths(self, rlengths: dict) -> bool:
        """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this fwtrack object will be updated.

        If chromosome in this fwtrack is not covered by given
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

        If self.rlength is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

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
        """Return all the chromosome names stored in this track object.
        """
        return set(sorted(self.locations.keys()))

    @cython.ccall
    def sort(self):
        """Naive sorting for locations.

        """
        c: bytes
        chrnames: set

        chrnames = self.get_chr_names()

        for c in chrnames:
            self.locations[c][0].sort()
            self.locations[c][1].sort()

        self.is_sorted = True
        return

    @cython.boundscheck(False)  # do not check that np indices are valid
    @cython.ccall
    def filter_dup(self, maxnum: cython.int = -1) -> cython.ulong:
        """Filter the duplicated reads.

        Run it right after you add all data into this object.

        Note, this function will *throw out* duplicates
        permenantly. If you want to keep them, use separate_dups
        instead.
        """
        p: cython.int
        n: cython.int
        current_loc: cython.int
        # index for old array, and index for new one
        i_old: cython.ulong
        i_new: cython.ulong
        size: cython.ulong
        k: bytes
        plus: cnp.ndarray(cython.int, ndim=1)
        new_plus: cnp.ndarray(cython.int, ndim=1)
        minus: cnp.ndarray(cython.int, ndim=1)
        new_minus: cnp.ndarray(cython.int, ndim=1)
        chrnames: set

        if maxnum < 0:
            return self.total         # do nothing

        if not self.is_sorted:
            self.sort()

        self.total = 0
        self.length = 0

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            # + strand
            i_new = 0
            plus = self.locations[k][0]
            size = plus.shape[0]
            if len(plus) <= 1:
                new_plus = plus         # do nothing
            else:
                new_plus = np.zeros(self.pointer[k][0] + 1, dtype='i4')
                new_plus[i_new] = plus[i_new]  # first item
                i_new += 1
                # the number of tags in the current location
                n = 1
                current_loc = plus[0]
                for i_old in range(1, size):
                    p = plus[i_old]
                    if p == current_loc:
                        n += 1
                    else:
                        current_loc = p
                        n = 1
                    if n <= maxnum:
                        new_plus[i_new] = p
                        i_new += 1
                new_plus.resize(i_new, refcheck=False)
                self.total += i_new
                self.pointer[k][0] = i_new
                # free memory?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                plus.resize(self.buffer_size, refcheck=False)
                plus.resize(0, refcheck=False)
                # hope there would be no mem leak...

            # - strand
            i_new = 0
            minus = self.locations[k][1]
            size = minus.shape[0]
            if len(minus) <= 1:
                new_minus = minus         # do nothing
            else:
                new_minus = np.zeros(self.pointer[k][1] + 1,
                                     dtype='i4')
                new_minus[i_new] = minus[i_new]  # first item
                i_new += 1
                # the number of tags in the current location
                n = 1
                current_loc = minus[0]
                for i_old in range(1, size):
                    p = minus[i_old]
                    if p == current_loc:
                        n += 1
                    else:
                        current_loc = p
                        n = 1
                    if n <= maxnum:
                        new_minus[i_new] = p
                        i_new += 1
                new_minus.resize(i_new, refcheck=False)
                self.total += i_new
                self.pointer[k][1] = i_new
                # free memory ?
                # I know I should shrink it to 0 size directly,
                # however, on Mac OSX, it seems directly assigning 0
                # doesn't do a thing.
                minus.resize(self.buffer_size, refcheck=False)
                minus.resize(0, refcheck=False)
                # hope there would be no mem leak...

            self.locations[k] = [new_plus, new_minus]

        self.length = self.fw * self.total
        return self.total

    @cython.ccall
    def sample_percent(self, percent: cython.float, seed: cython.int = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        num: cython.int  # num: number of reads allowed on a certain chromosome
        k: bytes
        chrnames: set

        self.total = 0
        self.length = 0

        chrnames = self.get_chr_names()

        if seed >= 0:
            np.random.seed(seed)

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            num = cython.cast(cython.int,
                              round(self.locations[k][0].shape[0] * percent, 5))
            np.random.shuffle(self.locations[k][0])
            self.locations[k][0].resize(num, refcheck=False)
            self.locations[k][0].sort()
            self.pointer[k][0] = self.locations[k][0].shape[0]

            num = cython.cast(cython.int,
                              round(self.locations[k][1].shape[0] * percent, 5))
            np.random.shuffle(self.locations[k][1])
            self.locations[k][1].resize(num, refcheck=False)
            self.locations[k][1].sort()
            self.pointer[k][1] = self.locations[k][1].shape[0]

            self.total += self.pointer[k][0] + self.pointer[k][1]

        self.length = self.fw * self.total
        return

    @cython.ccall
    def sample_num(self, samplesize: cython.ulong, seed: cython.int = -1):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize) / self.total
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def print_to_bed(self, fhd=None):
        """Output FWTrack to BED format files. If fhd is given,
        write to a file, otherwise, output to standard output.

        """
        i: cython.int
        p: cython.int
        k: bytes
        chrnames: set

        if not fhd:
            fhd = sys.stdout
        assert isinstance(fhd, io.IOBase)
        assert self.fw > 0, "FWTrack object .fw should be set larger than 0!"

        chrnames = self.get_chr_names()

        for k in chrnames:
            # for each chromosome.
            # This loop body is too big, I may need to split code later...

            plus = self.locations[k][0]

            for i in range(plus.shape[0]):
                p = plus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k.decode(),
                                                      p,
                                                      p + self.fw,
                                                      "+"))

            minus = self.locations[k][1]

            for i in range(minus.shape[0]):
                p = minus[i]
                fhd.write("%s\t%d\t%d\t.\t.\t%s\n" % (k.decode(),
                                                      p-self.fw,
                                                      p,
                                                      "-"))
        return

    @cython.ccall
    def extract_region_tags(self, chromosome: bytes,
                            startpos: cython.int, endpos: cython.int) -> tuple:
        i: cython.int
        pos: cython.int
        rt_plus: np.ndarray(cython.int, ndim=1)
        rt_minus: np.ndarray(cython.int, ndim=1)
        temp: list
        chrnames: set

        if not self.is_sorted:
            self.sort()

        chrnames = self.get_chr_names()
        assert chromosome in chrnames, "chromosome %s can't be found in the FWTrack object." % chromosome

        (plus, minus) = self.locations[chromosome]

        temp = []
        for i in range(plus.shape[0]):
            pos = plus[i]
            if pos < startpos:
                continue
            elif pos > endpos:
                break
            else:
                temp.append(pos)
        rt_plus = np.array(temp)

        temp = []
        for i in range(minus.shape[0]):
            pos = minus[i]
            if pos < startpos:
                continue
            elif pos > endpos:
                break
            else:
                temp.append(pos)
        rt_minus = np.array(temp)
        return (rt_plus, rt_minus)

    @cython.ccall
    def compute_region_tags_from_peaks(self, peaks: PeakIO,
                                       func,
                                       window_size: cython.int = 100,
                                       cutoff: cython.float = 5.0) -> list:
        """Extract tags in peak, then apply func on extracted tags.

        peaks: redefined regions to extract raw tags in PeakIO type: check cPeakIO.pyx.

        func:  a function to compute *something* from tags found in a predefined region

        window_size: this will be passed to func.

        cutoff: this will be passed to func.

        func needs the fixed number of parameters, so it's not flexible. Here is an example:

        wtd_find_summit(chrom, plus, minus, peak_start, peak_end, name , window_size, cutoff):

        """
        m: cython.int
        i: cython.int
        j: cython.int
        pos: cython.int
        startpos: cython.int
        endpos: cython.int

        plus: cnp.ndarray(cython.int, ndim=1)
        minus: cnp.ndarray(cython.int, ndim=1)
        rt_plus: cnp.ndarray(cython.int, ndim=1)
        rt_minus: cnp.ndarray(cython.int, ndim=1)

        chrom: bytes
        name: bytes

        temp: list
        retval: list
        pchrnames: set
        chrnames: set

        pchrnames = peaks.get_chr_names()
        retval = []

        # this object should be sorted
        if not self.is_sorted:
            self.sort()
        # PeakIO object should be sorted
        peaks.sort()

        chrnames = self.get_chr_names()

        for chrom in sorted(pchrnames):
            assert chrom in chrnames, "chromosome %s can't be found in the FWTrack object." % chrom
            (plus, minus) = self.locations[chrom]
            cpeaks = peaks.get_data_from_chrom(chrom)
            prev_i = 0
            prev_j = 0
            for m in range(len(cpeaks)):
                startpos = cpeaks[m]["start"] - window_size
                endpos = cpeaks[m]["end"] + window_size
                name = cpeaks[m]["name"]

                temp = []
                for i in range(prev_i, plus.shape[0]):
                    pos = plus[i]
                    if pos < startpos:
                        continue
                    elif pos > endpos:
                        prev_i = i
                        break
                    else:
                        temp.append(pos)
                rt_plus = np.array(temp, dtype="i4")

                temp = []
                for j in range(prev_j, minus.shape[0]):
                    pos = minus[j]
                    if pos < startpos:
                        continue
                    elif pos > endpos:
                        prev_j = j
                        break
                    else:
                        temp.append(pos)
                rt_minus = np.array(temp, dtype="i4")

                retval.append(func(chrom, rt_plus, rt_minus, startpos, endpos,
                                   name=name,
                                   window_size=window_size,
                                   cutoff=cutoff))
                # rewind window_size
                for i in range(prev_i, 0, -1):
                    if plus[prev_i] - plus[i] >= window_size:
                        break
                prev_i = i

                for j in range(prev_j, 0, -1):
                    if minus[prev_j] - minus[j] >= window_size:
                        break
                prev_j = j
                # end of a loop

        return retval

    @cython.ccall
    def pileup_a_chromosome(self, chrom: bytes, ds: list,
                            scale_factor_s: list,
                            baseline_value: cython.float = 0.0,
                            directional: bool = True,
                            end_shift: cython.int = 0) -> list:
        """pileup a certain chromosome, return [p,v] (end position and
        value) list.

        ds : tag will be extended to this value to 3' direction,
             unless directional is False. Can contain multiple
             extension values. Final pileup will the maximum.

        scale_factor_s : linearly scale the pileup value applied to
                         each d in ds. The list should have the same
                         length as ds.

        baseline_value : a value to be filled for missing values, and
                         will be the minimum pileup.

        directional : if False, the strand or direction of tag will be
                      ignored, so that extension will be both sides
                      with d/2.

        end_shift : move cutting ends towards 5->3 direction if value
                    is positive, or towards 3->5 direction if
                    negative. Default is 0 -- no shift at all.

        p and v are numpy.ndarray objects.
        """
        d: cython.long
        five_shift: cython.long
        # adjustment to 5' end and 3' end positions to make a fragment
        three_shift: cython.long
        rlength: cython.long
        chrlengths: dict
        five_shift_s: list = []
        three_shift_s: list = []
        tmp_pileup: list
        prev_pileup: list

        chrlengths = self.get_rlengths()
        rlength = chrlengths[chrom]
        assert len(ds) == len(scale_factor_s), "ds and scale_factor_s must have the same length!"

        # adjust extension length according to 'directional' and
        # 'halfextension' setting.
        for d in ds:
            if directional:
                # only extend to 3' side
                five_shift_s.append(- end_shift)
                three_shift_s.append(end_shift + d)
            else:
                # both sides
                five_shift_s.append(d//2 - end_shift)
                three_shift_s.append(end_shift + d - d//2)

        prev_pileup = None

        for i in range(len(ds)):
            five_shift = five_shift_s[i]
            three_shift = three_shift_s[i]
            scale_factor = scale_factor_s[i]
            tmp_pileup = se_all_in_one_pileup(self.locations[chrom][0],
                                              self.locations[chrom][1],
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

        return prev_pileup


@cython.inline
@cython.cfunc
def left_sum(data,
             pos: cython.int,
             width: cython.int) -> cython.int:
    return sum([data[x] for x in data if x <= pos and x >= pos - width])


@cython.inline
@cython.cfunc
def right_sum(data,
              pos: cython.int,
              width: cython.int) -> cython.int:
    return sum([data[x] for x in data if x >= pos and x <= pos + width])


@cython.inline
@cython.cfunc
def left_forward(data,
                 pos: cython.int,
                 window_size: cython.int) -> cython.int:
    return data.get(pos, 0) - data.get(pos-window_size, 0)


@cython.inline
@cython.cfunc
def right_forward(data,
                  pos: cython.int,
                  window_size: cython.int) -> cython.int:
    return data.get(pos + window_size, 0) - data.get(pos, 0)
