# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-11-12 17:20:17 Tao Liu>

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
from collections import Counter,defaultdict
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
    """In-memory paired-end fragment container grouped by chromosome.
    
    The track exposes utilities for sorting, filtering, downsampling, and pileup
    generation on numpy structured arrays of left/right coordinates.
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
    length = cython.declare(cython.ulonglong, visibility="public")
    average_template_length = cython.declare(cython.float, visibility="public")
    is_destroyed: bool

    def __init__(self, anno: str = "", buffer_size: cython.long = 100000):
        """Initialize an empty paired-end track.
        
        Parameters
        ----------
        anno : str, optional
            Annotation label retained with the track metadata.
        buffer_size : int, optional
            Number of fragment slots allocated per growth chunk for each chromosome.
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
        """Append a paired-end fragment to the track.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) that owns the fragment.
        start : int
            Zero-based start coordinate of the fragment (5' end).
        end : int
            Zero-based end coordinate of the fragment (3' end).
        
        Notes
        -----
        Fragments are stored in structured numpy arrays keyed by chromosome, and
        the running fragment count and total template length are updated in place.
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
        """Release numpy buffers held by the track.
        
        All per-chromosome arrays are resized to zero so the memory footprint is
        returned to the allocator, and the track is marked as destroyed.
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
        """Attach reference chromosome lengths to the track.
        
        Parameters
        ----------
        rlengths : dict
            Mapping from chromosome name (bytes) to reference length.
        
        Returns
        -------
        bool
            True when the length mapping has been updated.
        
        Notes
        -----
        Any chromosome stored in the track but missing from ``rlengths`` is
        assigned ``INT_MAX`` so downstream bounds checks can succeed.
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
        """Return the reference chromosome lengths associated with the track.
        
        Returns
        -------
        dict
            Mapping from chromosome name (bytes) to reference length. Chromosomes
            without a recorded length default to ``INT_MAX``.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

    @cython.ccall
    def finalize(self):
        """Shrink backing arrays and sort fragments in place.
        
        Each per-chromosome array is resized to the observed fragment count, sorted
        by the left and right coordinates, and the aggregate counters ``total`` and
        ``average_template_length`` are refreshed. Call this after loading data.
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
        """Return the fragment array for a chromosome.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name, provided as bytes.
        
        Returns
        -------
        numpy.ndarray
            Structured array with ``('l', 'i4')`` and ``('r', 'i4')`` fields.
        
        Raises
        ------
        Exception
            If the chromosome is not present in the track.
        """
        if chromosome in self.locations:
            return self.locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return the set of chromosome names stored in the track.
        
        Returns
        -------
        set
            Chromosome names (bytes) that currently have fragments.
        """
        return set(self.locations.keys())

    @cython.ccall
    def sort(self):
        """Sort fragments for each chromosome by genomic coordinate.
        
        Fragments are ordered first by their left coordinate and then by their right
        coordinate. The ``is_sorted`` flag is set to ``True`` when sorting completes.
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
        """Count observed fragment lengths across the track.
        
        Returns
        -------
        dict
            Mapping from fragment length to observed count, useful for downstream
            models such as HMMRATAC.
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
        """Return all fragment lengths as a single array.
        
        Returns
        -------
        numpy.ndarray
            Concatenated array of ``end - start`` for every stored fragment across
            chromosomes.
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
        """Remove fragments that overlap the provided exclusion regions.
        
        Parameters
        ----------
        regions : MACS3.Signal.Region.Regions
            Sorted region collection whose intervals should be excluded.
        
        Notes
        -----
        The operation mutates the track in place and finishes by calling
        :meth:`finalize` to refresh cached statistics.
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
                    self.length -= cython.cast(cython.ulonglong, r1[1] - r1[0])
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
        """Limit the number of duplicate fragments at identical coordinates.
        
        Parameters
        ----------
        maxnum : int, optional
            Maximum number of fragments allowed per unique ``(start, end)`` pair.
            A negative value disables duplicate filtering.
        
        Notes
        -----
        Fragments exceeding ``maxnum`` for the same coordinates are dropped and the
        aggregate template length is adjusted accordingly.
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
                        self.length -= cython.cast(cython.ulonglong, current_loc_end - current_loc_start)
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
        """Down-sample fragments in place by a fixed percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of fragments to keep per chromosome between 0 and 1 (inclusive).
        seed : int, optional
            Deterministic seed for the RNG; a negative value uses NumPy's global state.
        
        Notes
        -----
        Sampling is performed independently for each chromosome by shuffling the
        fragments, resizing the arrays, and restoring coordinate order.
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
        """Return a down-sampled copy of the track.
        
        Parameters
        ----------
        percent : float
            Fraction of fragments to retain per chromosome between 0 and 1.
        seed : int, optional
            Deterministic seed used when shuffling; a negative value disables seeding.
        
        Returns
        -------
        PETrackI
            New track containing the sampled fragments with metadata copied over.
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
        """Down-sample fragments in place to approximately ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target number of fragments across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent`.
        
        Notes
        -----
        The method converts ``samplesize`` into a sampling fraction using ``self.total``.
        Ensure :meth:`finalize` has been called so counts are up to date.
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize)/self.total
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def sample_num_copy(self,
                        samplesize: cython.ulong,
                        seed: cython.int = -1):
        """Return a down-sampled copy with approximately ``samplesize`` fragments.
        
        Parameters
        ----------
        samplesize : int
            Target number of fragments across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent_copy`.
        
        Returns
        -------
        PETrackI
            New track containing the sampled fragments.
        """
        percent: cython.float

        percent = cython.cast(cython.float, samplesize)/self.total
        return self.sample_percent_copy(percent, seed)

    @cython.ccall
    def print_to_bed(self, fhd=None):
        """Write fragments to a three-column BEDPE-style stream.
        
        Parameters
        ----------
        fhd : io.IOBase, optional
            Writable file-like object. Defaults to ``sys.stdout``.
        
        Notes
        -----
        Each fragment is emitted as ``chrom	start	end`` using decoded chromosome
        names and the stored integer coordinates.
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
        """Compute a coverage pileup for a single chromosome.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        scale_factor : float, optional
            Value used to scale the resulting coverage.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` with numpy arrays describing
            the pileup breakpoints and scaled coverage.
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
        """Project paired-end fragments into pseudo single-end pileups.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        ds : list[int]
            Fragment lengths used to build the projections.
        scale_factor_s : list[float]
            Scale factors paired with each entry in ``ds``.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` representing the merged pileup
            with the maximum value taken across projections.
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

        return prev_pileup

    @cython.ccall
    def pileup_bdg(self,
                   scale_factor: cython.float = 1.0,
                   baseline_value: cython.float = 0.0):
        """Build a ``bedGraphTrackI`` with pileups for every chromosome.
        
        Parameters
        ----------
        scale_factor : float, optional
            Value used to scale the coverage for each chromosome.
        baseline_value : float, optional
            Minimum value enforced on the coverage arrays.
        
        Returns
        -------
        bedGraphTrackI
            BedGraph track populated with per-chromosome pileup data.
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
        """Generate HMMRATAC-style pileups for every chromosome.
        
        Parameters
        ----------
        mapping : list
            Weight mapping produced by HMMRATAC EM training describing the short,
            mono-, di-, and tri-nucleosomal signals.
        baseline_value : float, optional
            Reserved parameter for API compatibility; not currently applied.
        
        Returns
        -------
        list
            List of dictionaries mirroring ``mapping`` where each dictionary maps
            chromosome names to pileup arrays returned by
            :func:`pileup_from_LR_hmmratac`.
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
    """Paired-end track for single-cell ATAC fragments with barcode metadata.
    
    Each chromosome stores a structured array of fragment coordinates and counts
    alongside an integer-encoded barcode array to support barcode-aware analyses.
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
    length = cython.declare(cython.ulonglong, visibility="public")  # total length of all fragments
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
        """Append a fragment together with its barcode and count.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name (as bytes) for the fragment.
        start : int
            Zero-based start coordinate of the fragment.
        end : int
            Zero-based end coordinate of the fragment.
        barcode : bytes
            Raw barcode sequence associated with the fragment.
        count : int
            Number of occurrences represented by the fragment.
        
        Notes
        -----
        Barcodes are interned into integers via ``barcode_dict`` for compact storage
        and the accumulated template length is weighted by ``count``.
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
        """Release fragment and barcode arrays held by the track.
        
        All per-chromosome arrays are resized to zero, barcode mappings are cleared,
        and the track is marked as destroyed.
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
        """Attach reference chromosome lengths to the track.
        
        Parameters
        ----------
        rlengths : dict
            Mapping from chromosome name (bytes) to reference length.
        
        Returns
        -------
        bool
            True when the length mapping has been updated.
        
        Notes
        -----
        Any chromosome stored in the track but missing from ``rlengths`` is assigned
        ``INT_MAX`` so downstream bounds checks can succeed.
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
        """Return the reference chromosome lengths associated with the track.
        
        Returns
        -------
        dict
            Mapping from chromosome name (bytes) to reference length. Chromosomes
            without a recorded length default to ``INT_MAX``.
        """
        if not self.rlengths:
            self.rlengths = dict([(k, INT_MAX) for k in self.locations.keys()])
        return self.rlengths

    @cython.ccall
    def finalize(self):
        """Shrink arrays, sort fragments, and refresh aggregate counters.
        
        Each per-chromosome fragment array is resized to its observed length, sorted
        by ``('l', 'r')``, and the accompanying barcode array is reordered to match.
        The method updates ``total`` and ``average_template_length`` using count
        weights and marks the track as sorted.
        
        Raises
        ------
        AssertionError
            If no fragments are present when finalizing.
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
        """Return the fragment array for a chromosome.
        
        Parameters
        ----------
        chromosome : bytes
            Chromosome name, provided as bytes.
        
        Returns
        -------
        numpy.ndarray
            Structured array with ``('l', 'i4')``, ``('r', 'i4')``, and ``('c', 'u2')`` fields.
        
        Raises
        ------
        Exception
            If the chromosome is not present in the track.
        """
        if chromosome in self.locations:
            return self.locations[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    @cython.ccall
    def get_chr_names(self) -> set:
        """Return the set of chromosome names stored in the track.
        
        Returns
        -------
        set
            Chromosome names (bytes) that currently have fragments.
        """
        return set(self.locations.keys())

    @cython.ccall
    def sort(self):
        """Sort fragments and barcodes for each chromosome.
        
        Fragments are ordered first by their left coordinate and then by their right
        coordinate, and the barcode array is reordered alongside the fragment array.
        The ``is_sorted`` flag is set to ``True`` when sorting completes.
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
        """Count fragment lengths weighted by per-fragment counts.
        
        Returns
        -------
        dict
            Mapping from fragment length to the total count contributed by fragments
            of that length.
        """
        sizes: cnp.ndarray(cnp.int32_t, ndim=1)
        s: cython.int
        locs: cnp.ndarray
        chrnames: list
        i: cython.int
        j: cython.int

        counter = Counter()
        chrnames = list(self.get_chr_names())
        for i in range(len(chrnames)):
            locs = self.locations[chrnames[i]]
            sizes = locs['r'] - locs['l']
            for j in range(len(sizes)):
                s = sizes[j]
                counter[s] += locs['c'][j]
        return dict(counter)

    @cython.ccall
    def fraglengths(self) -> cnp.ndarray:
        """Return all fragment lengths expanded by their counts.
        
        Returns
        -------
        numpy.ndarray
            Array of ``end - start`` values repeated according to the stored counts.
        """
        sizes: cnp.ndarray(np.int32_t, ndim=1)
        chrnames: list
        out: list
        chrom: bytes

        chrnames = list(self.get_chr_names())
        out = []
        for chrom in chrnames:
            locs = self.locations[chrom]
            sizes = locs['r'] - locs['l']
            counts = locs['c']
            out.append(np.repeat(sizes, counts))
        if out:
            return np.concatenate(out).astype(np.int32)
        else:
            return np.array([], dtype=np.int32)

    @cython.ccall
    def subset(self, selected_barcodes: set):
        """Build a new track containing only fragments from selected barcodes.
        
        Parameters
        ----------
        selected_barcodes : set
            Set of barcode byte strings to retain.
        
        Returns
        -------
        PETrackII
            New track restricted to the provided barcodes with metadata preserved.
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
        """Compute a coverage pileup for a single chromosome.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        scale_factor : float, optional
            Value used to scale the resulting coverage.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` with numpy arrays describing
            the pileup breakpoints and scaled coverage.
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
        """Project paired-end fragments into pseudo single-end pileups.
        
        Parameters
        ----------
        chrom : bytes
            Chromosome name to pile up.
        ds : list[int]
            Fragment lengths used to build the projections.
        scale_factor_s : list[float]
            Scale factors paired with each entry in ``ds``.
        baseline_value : float, optional
            Minimum value enforced on the coverage array.
        
        Returns
        -------
        list
            Two-element list ``[positions, values]`` representing the merged pileup
            with the maximum value taken across projections.
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
    def pileup_bdg(self,
                   scale_factor: cython.float = 1.0,
                   baseline_value: cython.float = 0.0):
        """Build a ``bedGraphTrackI`` with pileups for every chromosome.
        
        Parameters
        ----------
        scale_factor : float, optional
            Value used to scale the coverage for each chromosome.
        baseline_value : float, optional
            Minimum value enforced on the coverage arrays.
        
        Returns
        -------
        bedGraphTrackI
            BedGraph track populated with per-chromosome pileup data.
        """
        bdg: bedGraphTrackI
        pv: cnp.ndarray
        chrom: bytes

        bdg = bedGraphTrackI(baseline_value=baseline_value)
        for chrom in sorted(self.get_chr_names()):
            pv = pileup_from_LRC(self.locations[chrom])
            v = pv['v']
            v = v * scale_factor
            v[v < baseline_value] = baseline_value
            bdg.add_chrom_data(chrom,
                               pyarray('i', pv['p']),
                               pyarray('f', v))
        return bdg

    @cython.ccall
    def pileup_bdg2(self):
        """Build a ``bedGraphTrackII`` with pileups for every chromosome.
        
        Returns
        -------
        bedGraphTrackII
            BedGraph track populated with per-chromosome pileup arrays and finalized.
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
        """Remove fragments that overlap the provided exclusion regions.
        
        Parameters
        ----------
        regions : MACS3.Signal.Region.Regions
            Sorted region collection whose intervals should be excluded.
        
        Notes
        -----
        The operation mutates the track in place, adjusts fragment counts and lengths,
        and finishes by calling :meth:`finalize`.
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
                    self.length -= cython.cast(cython.ulonglong, (r1[1] - r1[0])*r1[2])
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

    @cython.ccall
    def sample_percent(self,
                       percent: cython.float,
                       seed: cython.int = -1):
        """Down-sample fragments in place so counts reflect a given percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of total counts to keep per chromosome between 0 and 1 (inclusive).
        seed : int, optional
            Deterministic seed for the RNG; a negative value uses NumPy's global state.
        
        Notes
        -----
        Fragments are sampled proportionally to their counts by expanding to an index
        vector, shuffling, and collapsing counts for the retained entries. Aggregate
        statistics are recomputed and the result is resorted.
        """
        k: bytes
        loc: cnp.ndarray
        bar: cnp.ndarray
        counts: cnp.ndarray
        n: cython.uint
        n_sample: cython.uint
        idx_flat: cnp.ndarray
        unique_idx: cnp.ndarray
        new_counts: cnp.ndarray
        new_locs: cnp.ndarray
        new_bars: cnp.ndarray

        assert 0.0 <= percent <= 1.0, "percent must be in [0, 1]"
        chrnames = sorted(self.get_chr_names())

        # Setup shuffling logic like PETrackI
        if seed >= 0:
            info(f"#   A random seed {seed} has been used")
            rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
            rs_shuffle = rs.shuffle
        else:
            rs_shuffle = np.random.shuffle

        self.length = 0
        self.total = 0
        self.average_template_length = 0.0

        for k in chrnames:
            loc = self.locations[k]
            bar = self.barcodes[k]
            counts = loc['c']
            n = int(counts.sum())
            n_sample = int(round(n * percent))
            if n == 0 or n_sample == 0:
                self.locations[k] = loc[:0]
                self.barcodes[k] = bar[:0]
                self.size[k] = 0
                continue

            # Flatten: build an array of indices into loc, repeated by count
            idx_flat = np.repeat(np.arange(len(loc)), counts)
            rs_shuffle(idx_flat)
            idx_flat = idx_flat[:n_sample]

            # Recount: count how many times each index is chosen
            unique_idx, new_counts = np.unique(idx_flat, return_counts=True)
            # Compose new arrays
            new_locs = loc[unique_idx].copy()
            new_locs['c'] = new_counts
            new_bars = bar[unique_idx].copy()
            self.locations[k] = new_locs
            self.barcodes[k] = new_bars
            self.size[k] = len(new_locs)
            self.length += np.sum((new_locs['r'] - new_locs['l']) * new_locs['c'])
            self.total += np.sum(new_locs['c'])

        if self.total > 0:
            self.average_template_length = float(self.length) / self.total
        else:
            self.average_template_length = 0.0
        self.sort()
        return

    @cython.ccall
    def sample_percent_copy(self,
                            percent: cython.float,
                            seed: cython.int = -1):
        """Return a down-sampled copy whose counts reflect a given percentage.
        
        Parameters
        ----------
        percent : float
            Fraction of total counts to keep per chromosome between 0 and 1 (inclusive).
        seed : int, optional
            Deterministic seed for the RNG; a negative value uses NumPy's global state.
        
        Returns
        -------
        PETrackII
            New track containing the sampled fragments with metadata preserved.
        
        Notes
        -----
        Fragments are sampled proportionally to their counts and the returned track is
        sorted with reference lengths copied from the source track.
        """
        k: bytes
        loc: cnp.ndarray
        bar: cnp.ndarray
        counts: cnp.ndarray
        n: cython.uint
        n_sample: cython.uint
        idx_flat: cnp.ndarray
        unique_idx: cnp.ndarray
        new_counts: cnp.ndarray
        new_locs: cnp.ndarray
        new_bars: cnp.ndarray

        assert 0.0 <= percent <= 1.0, "percent must be in [0, 1]"
        chrnames = sorted(self.get_chr_names())

        # Setup shuffling logic like PETrackI
        if seed >= 0:
            info(f"#   A random seed {seed} has been used")
            rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
            rs_shuffle = rs.shuffle
        else:
            rs_shuffle = np.random.shuffle

        ret = PETrackII(anno=self.annotation, buffer_size=self.buffer_size)
        ret.barcode_dict = dict(self.barcode_dict)
        ret.barcode_last_n = self.barcode_last_n

        ret.length = 0
        ret.total = 0

        for k in chrnames:
            loc = self.locations[k]
            bar = self.barcodes[k]
            counts = loc['c']
            n = int(counts.sum())
            n_sample = int(round(n * percent))
            if n == 0 or n_sample == 0:
                ret.locations[k] = loc[:0]
                ret.barcodes[k] = bar[:0]
                ret.size[k] = 0
                ret.buf_size[k] = 0
                continue

            idx_flat = np.repeat(np.arange(len(loc)), counts)
            rs_shuffle(idx_flat)
            idx_flat = idx_flat[:n_sample]
            unique_idx, new_counts = np.unique(idx_flat, return_counts=True)
            new_locs = loc[unique_idx].copy()
            new_locs['c'] = new_counts
            new_bars = bar[unique_idx].copy()
            ret.locations[k] = new_locs
            ret.barcodes[k] = new_bars
            ret.size[k] = len(new_locs)
            ret.buf_size[k] = len(new_locs)
            ret.length += np.sum((new_locs['r'] - new_locs['l']) * new_locs['c'])
            ret.total += np.sum(new_locs['c'])

        if ret.total > 0:
            ret.average_template_length = float(ret.length) / ret.total
        else:
            ret.average_template_length = 0.0
        ret.set_rlengths(self.get_rlengths())
        ret.sort()
        return ret

    @cython.ccall
    def sample_num(self,
                   samplesize: cython.ulong,
                   seed: cython.int = -1):
        """Down-sample fragments in place so total counts approximate ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target total count across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent`.
        
        Notes
        -----
        The method converts ``samplesize`` into a sampling fraction using the current
        total count and reuses :meth:`sample_percent`.
        """
        chrnames: set
        n_total: cython.uint
        chr_totals: dict
        k: bytes
        percent: cython.float

        chrnames = self.get_chr_names()
        n_total = 0
        chr_totals = {}
        for k in chrnames:
            chr_totals[k] = self.locations[k]['c'].sum()
            n_total += chr_totals[k]
        percent = 0.0 if n_total == 0 else min(samplesize / n_total, 1.0)
        self.sample_percent(percent, seed)
        return

    @cython.ccall
    def sample_num_copy(self,
                        samplesize: cython.ulong,
                        seed: cython.int = -1):
        """Return a down-sampled copy whose total counts approximate ``samplesize``.
        
        Parameters
        ----------
        samplesize : int
            Target total count across all chromosomes.
        seed : int, optional
            Deterministic seed forwarded to :meth:`sample_percent_copy`.
        
        Returns
        -------
        PETrackII
            New track containing the sampled fragments.
        """
        chrnames: set
        n_total: cython.uint
        chr_totals: dict
        k: bytes
        percent: cython.float

        chrnames = self.get_chr_names()
        n_total = 0
        chr_totals = {}
        for k in chrnames:
            chr_totals[k] = self.locations[k]['c'].sum()
            n_total += chr_totals[k]
        percent = 0.0 if n_total == 0 else min(samplesize / n_total, 1.0)
        return self.sample_percent_copy(percent, seed)

    @cython.ccall
    def pileup_bdg_hmmr(self,
                        mapping: list,
                        baseline_value: cython.float = 0.0) -> list:
        """Generate HMMRATAC-style pileups for every chromosome.
        
        Parameters
        ----------
        mapping : list
            Weight mapping produced by HMMRATAC EM training describing the short,
            mono-, di-, and tri-nucleosomal signals.
        baseline_value : float, optional
            Reserved parameter for API compatibility; not currently applied.
        
        Returns
        -------
        list
            List of dictionaries mirroring ``mapping`` where each dictionary maps
            chromosome names to pileup arrays returned by
            :func:`pileup_from_LR_hmmratac`.
        """
        ret_pileup: list
        i: cython.uint
        chroms: set
        chrom: bytes
        locs: cnp.ndarray
        counts: cnp.ndarray
        LR_expanded: cnp.ndarray
        idx: cnp.ndarray

        ret_pileup = []
        for i in range(len(mapping)):
            ret_pileup.append({})

        chroms = self.get_chr_names()
        for i in range(len(mapping)):
            for chrom in sorted(chroms):
                locs = self.locations[chrom]
                counts = locs['c']
                # Efficient numpy "explode"
                if locs.shape[0] == 0 or counts.sum() == 0:
                    LR_expanded = np.zeros((0,),
                                           dtype=[('l', 'i4'), ('r', 'i4')])
                else:
                    idx = np.repeat(np.arange(locs.shape[0]), counts)
                    LR_expanded = np.empty((len(idx),),
                                           dtype=[('l', 'i4'), ('r', 'i4')])
                    LR_expanded['l'] = locs['l'][idx]
                    LR_expanded['r'] = locs['r'][idx]
                ret_pileup[i][chrom] = pileup_from_LR_hmmratac(LR_expanded, mapping[i])
        return ret_pileup
