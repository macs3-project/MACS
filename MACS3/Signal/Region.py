# cython: language_level=3
# cython: profile=True
# Time-stamp: <2025-11-14 16:55:58 Tao Liu>

"""Module for Region classe.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
# import random
# import re
# import sys

# ------------------------------------
# MACS3 modules
# ------------------------------------

# from MACS3.Utilities.Constants import *

# ------------------------------------
# Other modules
# ------------------------------------
import cython
from cython.cimports.cpython import bool

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


@cython.cclass
class Regions:
    """Container for genomic regions organized by chromosome.
    
    This class stores genomic regions (intervals) organized by chromosome
    and provides operations for region manipulation including sorting,
    merging overlaps, expanding, and set operations like intersection
    and exclusion.
    
    Attributes:
        regions (dict): Dictionary mapping chromosome names (bytes) to 
            lists of region tuples (start, end). Regions are 0-based,
            half-open intervals [start, end).
        total (int): Total count of regions across all chromosomes.
    """
    regions = cython.declare(dict, visibility='public')
    total = cython.declare(cython.int, visibility='public')
    __sorted: bool
    __merged: bool

    def __init__(self):
        self.regions = {}
        self.__sorted = False
        self.__merged = False
        self.total = 0

    def __getitem__(self, chrom):
        return self.regions[chrom]

    @cython.ccall
    def pop(self, n: cython.int):
        """Remove and return the first n regions.
        
        Removes the first n regions from the Regions object across 
        chromosomes in sorted order and returns them as a new Regions 
        object. The current object is modified.
        
        Args:
            n (int): Number of regions to pop.
            
        Returns:
            Regions: New Regions object containing the first n regions.
            
        Examples:
            .. code-block:: python

                regions.add_loc(b'chr1', 1000, 2000)
                regions.add_loc(b'chr1', 3000, 4000)
                regions.add_loc(b'chr2', 5000, 6000)
                regions.sort()
                first_two = regions.pop(2)
                print(first_two.total)  # Output: 2
                print(regions.total)    # Output: 1
        """
        # when called, pop the first n regions in Regions class. Self
        # will be modified.
        clist: list             # for chromosomes
        tmp_l: cython.int
        chrom: bytes
        n_taken: cython.int     # remember the number of regions prepared
        ret: object             # returned Regions

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
                self.regions.pop(chrom)
            tmp_l = len(ret.regions[chrom])
            ret.total += tmp_l
            # calculate remained
            self.total -= tmp_l
            # print(ret.total, self.total)
            n_taken -= tmp_l
            if not n_taken:
                # when there is no need, quit loop
                break
        return ret

    @cython.ccall
    def init_from_PeakIO(self, peaks):
        """Initialize the object with a PeakIO object.

        Note: I intentionally forgot to check if peakio is actually a
        PeakIO...
        """
        chrom: bytes
        ps: list
        p: object
        i: int

        peaks.sort()
        self.total = 0
        for chrom in sorted(peaks.get_chr_names()):
            ps = peaks.get_data_from_chrom(chrom)
            self.regions[chrom] = []
            for i in range(len(ps)):
                p = ps[i]
                self.regions[chrom].append((p['start'], p['end']))
                self.total += 1
        self.sort()

    @cython.ccall
    def add_loc(self, chrom: bytes, start: int, end: int):
        """Add a region to a specific chromosome.
        
        Adds a single region interval to the specified chromosome. The 
        region is added without immediate sorting. Call sort() to 
        maintain order.
        
        Args:
            chrom (bytes): Chromosome name.
            start (int): Start position (0-based, inclusive).
            end (int): End position (0-based, exclusive).
            
        Returns:
            None
        """
        if self.regions.has_key(chrom):
            self.regions[chrom].append((start, end))
        else:
            self.regions[chrom] = [(start, end),]
        self.total += 1
        self.__sorted = False
        self.__merged = False
        return

    @cython.ccall
    def sort(self):
        """Sort all regions within each chromosome by start position.
        
        Sorts regions for each chromosome. This method is idempotent—
        calling it multiple times has the same effect as calling once.
        
        Returns:
            None
        """
        chrom: bytes

        if self.__sorted:
            return
        for chrom in sorted(self.regions.keys()):
            self.regions[chrom].sort()
        self.__sorted = True

    @cython.ccall
    def total_length(self) -> cython.long:
        """Calculate the total length covered by all regions.
        
        Returns the sum of lengths of all regions across all chromosomes.
        Overlapping regions are first merged to avoid double-counting.
        
        Returns:
            int: Total nucleotide length covered by all regions.
            
        """
        chrom: bytes
        ps: list
        i: cython.int
        tl: cython.int
        s: cython.int
        e: cython.int
        self.merge_overlap()
        tl = 0
        for chrom in sorted(self.regions.keys()):
            ps = self.regions[chrom]
            for i in range(len(ps)):
                s, e = ps[i]
                tl += e - s
        return tl

    @cython.ccall
    def get_chr_names(self) -> set:
        """Get all chromosome names present in the Regions object.
        
        Returns:
            set: Set of chromosome names (bytes) present in regions.
            
        Examples:
            .. code-block:: python

                from MACS3.Signal.Region import Regions
                regions = Regions()
                regions.add_loc(b'chr1', 1000, 2000)
                regions.add_loc(b'chr1', 3000, 4000)
                regions.add_loc(b'chr2', 5000, 6000)
                print(regions.get_chr_names())
                # Output: {b'chr1', b'chr2'}
        """
        return set(sorted(self.regions.keys()))

    @cython.ccall
    def expand(self, flanking: int):
        """Expand all regions by a fixed distance in both directions.
        
        Extends the start and end positions of each region by 'flanking'
        base pairs. Start positions are expanded leftward but capped at 0.
        The regions are automatically re-sorted after expansion.
        
        Args:
            flanking (int): Number of base pairs to expand in each direction.
            
        Returns:
            None: Modifies the object in place.
            
        Examples:
            .. code-block:: python

                regions.expand(100)
                print(regions.regions[b'chr1'])  # Output: [(900, 2100)]
                
                # With capping at 0 for start positions
                regions2 = Regions()
                regions2.add_loc(b'chr1', 50, 150)
                regions2.expand(100)
                print(regions2.regions[b'chr1'])
                # Output: [(0, 250)]  # Start capped at 0
        """
        chrom: bytes
        ps: list
        i: cython.int

        self.sort()
        for chrom in sorted(self.regions.keys()):
            ps = self.regions[chrom]
            for i in range(len(ps)):
                ps[i] = (max(0, ps[i][0] - flanking), ps[i][1] + flanking)
            ps.sort()
            self.regions[chrom] = ps
        self.__merged = False

    @cython.ccall
    def merge_overlap(self):
        """Merge overlapping or adjacent regions within each chromosome.
        
        Combines regions that overlap or are adjacent. After merging,
        regions are re-sorted. This operation is idempotent and modifies
        the object in place. The total count is updated.
        
        Returns:
            bool: True if merge was performed, False if already merged.
            
        Examples:
            .. code-block:: python
            
                regions = Regions()
                regions.add_loc(b'chr1', 1000, 2000)
                regions.add_loc(b'chr1', 1500, 3000)  # Overlaps with first
                result = regions.merge_overlap()
                print(regions.regions[b'chr1'])  # Output: [(1000, 3000)]
        """
        chrom: bytes
        s_new_region: cython.int
        e_new_region: cython.int
        i: cython.int
        regions: dict
        new_regions: dict
        chrs: list
        regions_chr: list
        prev_region: tuple

        if self.__merged:
            return
        self.sort()
        regions = self.regions
        new_regions = {}

        chrs = list(regions.keys())
        chrs.sort()
        self.total = 0
        for i in range(len(chrs)):
            chrom = chrs[i]
            new_regions[chrom] = []
            n_append = new_regions[chrom].append
            prev_region = ()
            regions_chr = regions[chrom]
            for i in range(len(regions_chr)):
                if not prev_region:
                    prev_region = regions_chr[i]
                    continue
                else:
                    if regions_chr[i][0] <= prev_region[1]:
                        s_new_region = prev_region[0]
                        e_new_region = regions_chr[i][1]
                        prev_region = (s_new_region, e_new_region)
                    else:
                        n_append(prev_region)
                        prev_region = regions_chr[i]
            if prev_region:
                n_append(prev_region)
            self.total += len(new_regions[chrom])
        self.regions = new_regions
        self.sort()
        self.__merged = True
        return True

    @cython.ccall
    def write_to_bed(self, fhd):
        """Write regions to a BED format file.
        
        Args:
            fhd: File handle opened for writing.
            
        Examples:
            .. code-block:: python
            
                with open('output.bed', 'w') as f:
                    regions.write_to_bed(f)
        """
        i: cython.int
        chrom: bytes
        region: tuple

        chrs = list(self.regions.keys())
        chrs.sort()
        for i in range(len(chrs)):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                fhd.write("%s\t%d\t%d\n" % (chrom.decode(),
                                            region[0], region[1]))

    def __str__(self):
        i: cython.int
        chrom: bytes
        region: tuple
        ret: str

        ret = ""
        chrs = list(self.regions.keys())
        chrs.sort()
        for i in range(len(chrs)):
            chrom = chrs[i]
            for region in self.regions[chrom]:
                ret += "%s\t%d\t%d\n" % (chrom.decode(),
                                         region[0], region[1])
        return ret

    # cpdef object randomly_pick (self, int n, int seed = 12345):
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
    #     random.seed(seed)
    #     all_pc = random.shuffle(all_pc)[:n]
    #     ret_peakio = PeakIO()
    #     for p in all_pc:
    #         ret_peakio.add_PeakContent (p["chrom"], p)
    #     return ret_peakio

    @cython.ccall
    def intersect(self, regions_object2):
        """Get the intersection with another Regions object.
        
        Returns a new Regions object containing only the regions that
        overlap between this object and regions_object2. For chromosomes
        present only in this object, all regions are included unchanged.
        
        Args:
            regions_object2 (Regions): Another Regions object to intersect with.
            
        Returns:
            Regions: New object containing intersecting regions.
            
        Examples:
            .. code-block:: python
            
                r1 = Regions()
                r1.add_loc(b'chr1', 1000, 3000)
                
                r2 = Regions()
                r2.add_loc(b'chr1', 2000, 4000)
                
                result = r1.intersect(r2)
                print(result.regions[b'chr1'])  # Output: [(2000, 3000)]
        """
        ret_regions_object: object
        regions1: dict
        regions2: dict
        chrs1: list
        chrs2: list
        four_coords: list
        k: bytes
        ret_regions: dict
        r1: tuple
        r2: tuple
        n_rl1: cython.long
        n_rl2: cython.long

        self.sort()
        regions1 = self.regions
        assert isinstance(regions_object2, Regions)

        regions_object2.sort()
        regions2 = regions_object2.regions

        ret_regions_object = Regions()
        ret_regions = dict()
        ret_regions_object.total = 0
        chrs1 = list(regions1.keys())
        chrs2 = list(regions2.keys())
        for k in chrs1:
            # print(f"chromosome {k}")
            if not chrs2.count(k):
                # no such chromosome in peaks1, then don't touch the
                # peaks in this chromosome
                ret_regions[k] = regions1[k]
                ret_regions_object.total += len(ret_regions[k])
                continue
            ret_regions[k] = []
            n_rl1 = len(regions1[k])    # number of remaining elements in regions1[k]
            n_rl2 = len(regions2[k])    # number of remaining elements in regions2[k]
            rl1_k = iter(regions1[k]).__next__
            rl2_k = iter(regions2[k]).__next__
            r1 = rl1_k()
            n_rl1 -= 1
            r2 = rl2_k()
            n_rl2 -= 1
            while (True):
                # we do this until there is no r1 or r2 left.
                if r2[0] < r1[1] and r1[0] < r2[1]:
                    # We found an overlap, now get the intersecting
                    # region.
                    four_coords = sorted([r1[0], r1[1], r2[0], r2[1]])
                    ret_regions[k].append((four_coords[1], four_coords[2]))
                if r1[1] < r2[1]:
                    # in this case, we need to move to the next r1,
                    if n_rl1:
                        r1 = rl1_k()
                        n_rl1 -= 1
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
            ret_regions_object.total += len(ret_regions[k])

        ret_regions_object.regions = ret_regions
        ret_regions_object.sort()
        return ret_regions_object

    @cython.ccall
    def exclude(self, regions_object2):
        """Remove regions that overlap with another Regions object.
        
        Removes all regions from this object that overlap with regions
        in regions_object2. Regions present only in this object are kept.
        The object is modified in place.
        
        Args:
            regions_object2 (Regions): Another Regions object defining
                regions to exclude.
            
        Returns:
            None: Modifies the object in place.
            
        Examples:
            .. code-block:: python
                
                from MACS3.Signal.Region import Regions
                r1 = Regions()
                r1.add_loc(b'chr1', 1000, 3000)
                r1.add_loc(b'chr1', 5000, 7000)
                
                r2 = Regions()
                r2.add_loc(b'chr1', 2000, 6000)  # Overlaps both regions
                
                r1.exclude(r2)
                print(r1.regions[b'chr1'])
                # Output: [(1000, 2000), (6000, 7000)]
        """
        regions1: dict
        regions2: dict
        chrs1: list
        chrs2: list
        k: bytes
        ret_regions: dict
        overlap_found: bool
        r1: tuple
        r2: tuple
        n_rl1: cython.long
        n_rl2: cython.long

        self.sort()
        regions1 = self.regions
        self.total = 0
        assert isinstance(regions_object2, Regions)
        regions_object2.sort()
        regions2 = regions_object2.regions

        ret_regions = dict()
        chrs1 = list(regions1.keys())
        chrs2 = list(regions2.keys())
        for k in chrs1:
            # print(f"chromosome {k}")
            if not chrs2.count(k):
                # no such chromosome in peaks1, then don't touch the peaks in this chromosome
                ret_regions[k] = regions1[k]
                self.total += len(ret_regions[k])
                continue
            ret_regions[k] = []
            n_rl1 = len(regions1[k])
            n_rl2 = len(regions2[k])
            rl1_k = iter(regions1[k]).__next__
            rl2_k = iter(regions2[k]).__next__
            overlap_found = False
            r1 = rl1_k()
            n_rl1 -= 1
            r2 = rl2_k()
            n_rl2 -= 1
            while (True):
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
                        ret_regions[k].append(r1)
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
                ret_regions[k].extend(regions1[k][-n_rl1-1:])
            self.total += len(ret_regions[k])

        self.regions = ret_regions
        self.__sorted = False
        self.sort()
        return
