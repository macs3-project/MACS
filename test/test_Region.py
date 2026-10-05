#!/usr/bin/env python
# Time-stamp: <2025-09-29 15:10:49 Tao Liu>

import io
import unittest

import pytest

from MACS3.IO.PeakIO import PeakIO
from MACS3.Signal.Region import Regions


class Test_Regions(unittest.TestCase):
    def setUp(self):
        self.test_regions1 = [(b"chrY", 0, 100),
                              (b"chrY", 300, 500),
                              (b"chrY", 700, 900),
                              (b"chrY", 1000, 1200),
                              ]
        self.test_regions2 = [(b"chrY", 100, 200),
                              (b"chrY", 300, 400),
                              (b"chrY", 600, 800),
                              (b"chrY", 1200, 1300),
                              ]
        # when we add test_regions1 and test_region2 into the same
        # Regions, we should get this from 'merge_overlaps'
        self.merge_result_regions = [(b"chrY", 0, 200),
                                     (b"chrY", 300, 500),
                                     (b"chrY", 600, 900),
                                     (b"chrY", 1000, 1300),
                                     ]
        self.merge_result_bedcontent = "chrY\t0\t200\nchrY\t300\t500\nchrY\t600\t900\nchrY\t1000\t1300\n"

        self.test_regions3 = [(b"chr1", 0, 100),
                              (b"chr1", 300, 500),
                              (b"chr1", 700, 900),
                              (b"chr1", 1000, 1200),
                              (b"chrY", 100, 200),
                              (b"chrY", 300, 400),
                              (b"chrY", 600, 800),
                              (b"chrY", 1200, 1300),
                              ]
        self.popped_regions = ["chr1\t0\t100\nchr1\t300\t500\nchr1\t700\t900\n",
                                "chr1\t1000\t1200\nchrY\t100\t200\nchrY\t300\t400\n",
                                "chrY\t600\t800\nchrY\t1200\t1300\n"]

        # test_regions4 and test_regions5 are used to test the intersection
        self.test_regions4 = [(b"chrY", 0, 100),
                              (b"chrY", 300, 500),
                              (b"chrY", 700, 900),
                              (b"chrY", 1000, 1200)
                              ]
        self.test_regions5 = [(b"chrY", 100, 200),
                              (b"chrY", 300, 400),
                              (b"chrY", 600, 800),
                              (b"chrY", 1100, 1150),
                              (b"chrY", 1175, 1300)
                              ]
        # After we add test_regions4 and test_region5 into two Regions
        # objects, we should get this from 'intersect'
        self.intersect_regions_4_vs_5 = {b"chrY": [(300, 400),
                                                   (700, 800),
                                                   (1100, 1150),
                                                   (1175, 1200)]}

    def test_add_loc1(self):
        # make sure the shuffled sequence does not lose any elements
        self.r1 = Regions()
        for a in self.test_regions1:
            self.r1.add_loc(a[0], a[1], a[2])

    def test_add_loc2(self):
        # make sure the shuffled sequence does not lose any elements
        self.r2 = Regions()
        for a in self.test_regions2:
            self.r2.add_loc(a[0], a[1], a[2])

    def test_merge(self):
        self.mr = Regions()
        for a in self.test_regions1:
            self.mr.add_loc(a[0], a[1], a[2])
        for a in self.test_regions2:
            self.mr.add_loc(a[0], a[1], a[2])
        self.mr.merge_overlap()
        self.assertEqual(str(self.mr), self.merge_result_bedcontent)

    def test_pop3(self):
        self.r = Regions()
        for a in self.test_regions3:
            self.r.add_loc(a[0], a[1], a[2])
        # now pop 3 at a time
        ret_list_regions = []
        while self.r.total != 0:
            ret_list_regions.append(str(self.r.pop(3)))
        self.assertEqual(ret_list_regions, self.popped_regions)

    def test_intersect(self):
        self.r4 = Regions()
        for a in self.test_regions4:
            self.r4.add_loc(a[0], a[1], a[2])
        self.r5 = Regions()
        for a in self.test_regions5:
            self.r5.add_loc(a[0], a[1], a[2])
        self.intersect_4_vs_5 = self.r4.intersect(self.r5)
        self.assertEqual(self.intersect_4_vs_5.regions, self.intersect_regions_4_vs_5)
        self.intersect_4_vs_5 = self.r5.intersect(self.r4)
        self.assertEqual(self.intersect_4_vs_5.regions, self.intersect_regions_4_vs_5)


# ------------------------------------
# Reference helpers (independent of MACS3)
# ------------------------------------

def _make_regions(spec):
    """Build a Regions from {chrom: [(start, end), ...]} in the given order."""
    r = Regions()
    for chrom, intervals in spec.items():
        for s, e in intervals:
            r.add_loc(chrom, s, e)
    return r


def _bases(intervals):
    """Set of bases covered by half-open (start, end) intervals."""
    covered = set()
    for s, e in intervals:
        covered.update(range(s, e))
    return covered


def _runs(bases):
    """Maximal runs of consecutive bases as sorted half-open intervals."""
    out = []
    for b in sorted(bases):
        if out and out[-1][1] == b:
            out[-1][1] = b + 1
        else:
            out.append([b, b + 1])
    return [tuple(x) for x in out]


# ------------------------------------
# Regions.add_loc / __getitem__ / get_chr_names
# ------------------------------------

def test_add_loc_keeps_insertion_order_and_counts():
    r = Regions()
    r.add_loc(b"chr1", 50, 60)
    r.add_loc(b"chr2", 0, 5)
    r.add_loc(b"chr1", 10, 20)
    assert r.regions == {b"chr1": [(50, 60), (10, 20)], b"chr2": [(0, 5)]}
    assert r.total == 3
    assert r[b"chr1"] == [(50, 60), (10, 20)]


def test_new_regions_is_empty():
    r = Regions()
    assert r.regions == {}
    assert r.total == 0
    assert str(r) == ""


def test_add_loc_rejects_str_chromosome():
    with pytest.raises(TypeError):
        Regions().add_loc("chr1", 0, 10)


def test_getitem_missing_chromosome_raises_keyerror():
    with pytest.raises(KeyError):
        Regions()[b"chr1"]


def test_add_loc_int32_extremes_merge():
    # (0, 2**31-1) and (2**31-2, 2**31-1) share the last base -> one region
    r = Regions()
    r.add_loc(b"chr1", 2147483646, 2147483647)
    r.add_loc(b"chr1", 0, 2147483647)
    r.merge_overlap()
    assert r.regions == {b"chr1": [(0, 2147483647)]}
    assert r.total == 1


@pytest.mark.parametrize("spec,expected", [
    ({}, set()),
    ({b"chr1": [(0, 1)]}, {b"chr1"}),
    ({b"chr2": [(0, 1)], b"chr1": [(5, 6), (0, 1)], b"chrX": [(1, 2)]},
     {b"chr1", b"chr2", b"chrX"}),
], ids=["empty", "one", "many"])
def test_get_chr_names(spec, expected):
    assert _make_regions(spec).get_chr_names() == expected


# ------------------------------------
# Regions.sort
# ------------------------------------

def test_sort_orders_by_start_then_end_per_chromosome():
    r = _make_regions({b"chr1": [(30, 40), (0, 50), (0, 10)],
                       b"chr2": [(9, 10), (1, 2)]})
    r.sort()
    expected = {b"chr1": [(0, 10), (0, 50), (30, 40)],
                b"chr2": [(1, 2), (9, 10)]}
    assert r.regions == expected
    r.sort()  # idempotent
    assert r.regions == expected
    assert r.total == 5


def test_sort_again_after_add_loc():
    r = _make_regions({b"chr1": [(30, 40)]})
    r.sort()
    r.add_loc(b"chr1", 0, 10)
    r.sort()
    assert r.regions[b"chr1"] == [(0, 10), (30, 40)]


def test_sort_empty():
    r = Regions()
    r.sort()
    assert r.regions == {}


# ------------------------------------
# Regions.merge_overlap
# ------------------------------------

MERGE_CASES = [
    pytest.param([(0, 10), (20, 30)], id="disjoint"),
    pytest.param([(0, 10), (5, 20)], id="overlapping"),
    pytest.param([(0, 10), (10, 20)], id="touching"),
    pytest.param([(0, 10), (11, 20)], id="gap_of_one"),
    pytest.param([(50, 60), (0, 10), (5, 20)], id="unsorted"),
    pytest.param([(0, 5), (4, 9), (8, 12)], id="chain"),
    pytest.param([(3, 7)], id="single"),
    pytest.param([(0, 10), (0, 10)], id="duplicate"),
    pytest.param([(0, 10), (0, 5)], id="same_start"),
]


@pytest.mark.parametrize("intervals", MERGE_CASES)
def test_merge_overlap_matches_base_set_runs(intervals):
    # merged regions = maximal runs of covered bases (touching intervals
    # merge, as documented)
    r = _make_regions({b"chr1": intervals})
    r.merge_overlap()
    expected = _runs(_bases(intervals))
    assert r.regions[b"chr1"] == expected
    assert r.total == len(expected)


def test_merge_overlap_many_chromosomes():
    spec = {b"chr2": [(5, 15), (0, 6)],
            b"chr1": [(100, 200), (200, 300), (400, 500)],
            b"chrX": [(1, 2)]}
    r = _make_regions(spec)
    r.merge_overlap()
    for chrom, intervals in spec.items():
        assert r.regions[chrom] == _runs(_bases(intervals))
    assert r.total == 4
    assert str(r) == "chr1\t100\t300\nchr1\t400\t500\nchr2\t0\t15\nchrX\t1\t2\n"


def test_merge_overlap_empty():
    r = Regions()
    r.merge_overlap()
    assert r.regions == {}
    assert r.total == 0


def test_merge_overlap_returns_true_when_merging():
    r = _make_regions({b"chr1": [(0, 10), (5, 8)]})
    assert r.merge_overlap() is True


def test_merge_overlap_again_after_add_loc():
    r = _make_regions({b"chr1": [(0, 10)]})
    r.merge_overlap()
    r.add_loc(b"chr1", 5, 20)
    r.merge_overlap()
    assert r.regions[b"chr1"] == [(0, 20)]
    assert r.total == 1


# ------------------------------------
# Regions.total_length
# ------------------------------------

@pytest.mark.parametrize("spec", [
    pytest.param({}, id="empty"),
    pytest.param({b"chr1": [(0, 10)]}, id="single"),
    pytest.param({b"chr1": [(0, 10), (5, 20)]}, id="overlap_counted_once"),
    pytest.param({b"chr1": [(0, 10), (10, 20)]}, id="touching"),
    pytest.param({b"chr1": [(0, 10)], b"chr2": [(0, 10)]},
                 id="same_coordinates_two_chromosomes"),
    pytest.param({b"chr1": [(30, 40), (0, 5)], b"chr2": [(7, 9)]},
                 id="unsorted_many"),
])
def test_total_length_counts_covered_bases(spec):
    r = _make_regions(spec)
    expected = sum(len(_bases(v)) for v in spec.values())
    assert r.total_length() == expected


def test_total_length_merges_in_place():
    r = _make_regions({b"chr1": [(5, 20), (0, 10)]})
    assert r.total_length() == 20
    assert r.regions == {b"chr1": [(0, 20)]}
    assert r.total == 1


def test_total_length_int32_max():
    r = _make_regions({b"chr1": [(0, 2147483647)]})
    assert r.total_length() == 2147483647


# ------------------------------------
# Regions.expand
# ------------------------------------

@pytest.mark.parametrize("intervals,flank,expected", [
    pytest.param([(100, 200)], 50, [(50, 250)], id="simple"),
    pytest.param([(30, 200)], 50, [(0, 250)], id="capped_at_zero"),
    pytest.param([(50, 150)], 100, [(0, 250)], id="docstring_example"),
    pytest.param([(0, 10)], 0, [(0, 10)], id="zero_flank"),
    pytest.param([(500, 600), (100, 200)], 10, [(90, 210), (490, 610)],
                 id="unsorted"),
    pytest.param([(100, 200), (250, 300)], 30, [(70, 230), (220, 330)],
                 id="overlap_not_merged"),
    pytest.param([(5, 100), (0, 10)], 20, [(0, 30), (0, 120)],
                 id="both_capped"),
])
def test_expand(intervals, flank, expected):
    # reference: (max(0, s - flank), e + flank), then sorted
    assert sorted((max(0, s - flank), e + flank)
                  for s, e in intervals) == expected
    r = _make_regions({b"chr1": intervals})
    r.expand(flank)
    assert r.regions[b"chr1"] == expected
    assert r.total == len(intervals)


def test_expand_clears_merged_state_many_chromosomes():
    r = _make_regions({b"chr1": [(100, 200), (250, 300)], b"chr2": [(10, 20)]})
    r.merge_overlap()
    r.expand(30)
    r.merge_overlap()
    assert r.regions == {b"chr1": [(70, 330)], b"chr2": [(0, 50)]}
    assert r.total == 2


def test_expand_empty():
    r = Regions()
    r.expand(100)
    assert r.regions == {}
    assert r.total == 0


# ------------------------------------
# Regions.pop
# ------------------------------------

def test_pop_more_than_total_takes_all():
    r = _make_regions({b"chr1": [(20, 30), (0, 10)], b"chr2": [(5, 6)]})
    r.sort()
    got = r.pop(10)
    assert got.total == 3
    assert str(got) == "chr1\t0\t10\nchr1\t20\t30\nchr2\t5\t6\n"
    assert r.total == 0
    assert r.regions == {}
    with pytest.raises(Exception, match="^None left$"):
        r.pop(1)


def test_pop_crosses_chromosome_boundary():
    r = _make_regions({b"chr2": [(0, 1), (2, 3)], b"chr1": [(10, 20), (30, 40)]})
    r.sort()
    got = r.pop(3)
    assert got.regions == {b"chr1": [(10, 20), (30, 40)], b"chr2": [(0, 1)]}
    assert got.total == 3
    assert r.regions == {b"chr2": [(2, 3)]}
    assert r.total == 1


def test_pop_exhausts_one_chromosome():
    r = _make_regions({b"chr1": [(0, 1), (2, 3)], b"chr2": [(4, 5)]})
    got = r.pop(2)
    assert got.regions == {b"chr1": [(0, 1), (2, 3)]}
    assert r.regions == {b"chr2": [(4, 5)]}
    assert (got.total, r.total) == (2, 1)


def test_pop_orders_chromosomes_as_bytes():
    r = _make_regions({b"chr2": [(0, 1)], b"chr10": [(0, 1)], b"chr1": [(0, 1)]})
    order = [list(r.pop(1).regions.keys()) for _ in range(3)]
    # bytes order: b"chr1" < b"chr10" < b"chr2"
    assert order == [[b"chr1"], [b"chr10"], [b"chr2"]]
    assert r.total == 0


def test_pop_zero_takes_nothing():
    r = _make_regions({b"chr1": [(0, 10)]})
    got = r.pop(0)
    assert got.total == 0
    assert str(got) == ""
    assert r.regions == {b"chr1": [(0, 10)]}
    assert r.total == 1


def test_pop_empty_raises():
    with pytest.raises(Exception, match="^None left$"):
        Regions().pop(1)


# ------------------------------------
# Regions.init_from_PeakIO
# ------------------------------------

def test_init_from_PeakIO_sorts_and_counts():
    p = PeakIO()
    p.add(b"chr2", 500, 600)
    p.add(b"chr1", 300, 400)
    p.add(b"chr1", 100, 200)
    p.add(b"chr2", 50, 60)
    r = Regions()
    r.init_from_PeakIO(p)
    assert r.regions == {b"chr1": [(100, 200), (300, 400)],
                         b"chr2": [(50, 60), (500, 600)]}
    assert r.total == 4
    assert r.get_chr_names() == {b"chr1", b"chr2"}


def test_init_from_PeakIO_breaks_start_ties_by_end():
    p = PeakIO()
    p.add(b"chr1", 0, 50)
    p.add(b"chr1", 0, 10)
    r = Regions()
    r.init_from_PeakIO(p)
    assert r.regions == {b"chr1": [(0, 10), (0, 50)]}


def test_init_from_PeakIO_keeps_subpeak_duplicates():
    p = PeakIO()
    p.add(b"chr1", 100, 200, summit=120)
    p.add(b"chr1", 100, 200, summit=180)
    r = Regions()
    r.init_from_PeakIO(p)
    assert r.regions == {b"chr1": [(100, 200), (100, 200)]}
    assert r.total == 2
    r.merge_overlap()
    assert r.regions == {b"chr1": [(100, 200)]}
    assert r.total == 1


def test_init_from_PeakIO_empty():
    r = Regions()
    r.init_from_PeakIO(PeakIO())
    assert r.regions == {}
    assert r.total == 0


# ------------------------------------
# Regions.write_to_bed / __str__
# ------------------------------------

def test_write_to_bed_exact_text(tmp_path):
    r = _make_regions({b"chrX": [(5, 6)], b"chr1": [(100, 200), (0, 50)],
                       b"chr10": [(1, 2)]})
    r.sort()
    path = tmp_path / "r.bed"
    with open(path, "w") as fh:
        r.write_to_bed(fh)
    expected = "chr1\t0\t50\nchr1\t100\t200\nchr10\t1\t2\nchrX\t5\t6\n"
    assert path.read_text() == expected
    assert str(r) == expected


def test_write_to_bed_empty():
    fh = io.StringIO()
    Regions().write_to_bed(fh)
    assert fh.getvalue() == ""


# ------------------------------------
# Regions.intersect
# ------------------------------------

# each operand's intervals are separated by gaps, so the exact pieces are
# the runs of the base-set intersection
INTERSECT_CASES = [
    pytest.param([(0, 100)], [(50, 150)], id="partial"),
    pytest.param([(0, 100)], [(100, 200)], id="touching"),
    pytest.param([(0, 100)], [(10, 20)], id="nested"),
    pytest.param([(0, 10)], [(0, 10)], id="identical"),
    pytest.param([(0, 1000)], [(10, 20), (30, 40), (50, 60)], id="one_vs_many"),
    pytest.param([(10, 20), (30, 40), (50, 60)], [(0, 1000)], id="many_vs_one"),
    pytest.param([(0, 10)], [(20, 30)], id="disjoint"),
    pytest.param([(300, 400), (0, 100)], [(50, 350)], id="unsorted"),
    pytest.param([(0, 10), (20, 30), (40, 50)], [(5, 25), (45, 60)],
                 id="interleaved"),
]


@pytest.mark.parametrize("a,b", INTERSECT_CASES)
def test_intersect_matches_base_sets(a, b):
    got = _make_regions({b"chr1": a}).intersect(_make_regions({b"chr1": b}))
    expected = _runs(_bases(a) & _bases(b))
    assert got.regions.get(b"chr1", []) == expected
    assert got.total == len(expected)


@pytest.mark.parametrize("a,b", [
    pytest.param([(0, 50), (25, 75)], [(40, 60)], id="self_overlapping"),
    pytest.param([(0, 100), (10, 20)], [(5, 15), (90, 120)], id="self_nested"),
    pytest.param([(0, 30), (10, 100)], [(40, 50), (45, 200)],
                 id="both_overlapping"),
])
def test_intersect_covered_bases_with_overlapping_operands(a, b):
    ra, rb = _make_regions({b"chr1": a}), _make_regions({b"chr1": b})
    expected = _bases(a) & _bases(b)
    assert _bases(ra.intersect(rb).regions.get(b"chr1", [])) == expected
    assert _bases(rb.intersect(ra).regions.get(b"chr1", [])) == expected


def test_intersect_ignores_chromosome_only_in_other():
    a = _make_regions({b"chr1": [(0, 100)]})
    b = _make_regions({b"chr1": [(50, 60)], b"chr2": [(0, 100)]})
    got = a.intersect(b)
    assert got.regions == {b"chr1": [(50, 60)]}
    assert got.total == 1


@pytest.mark.parametrize("other", [{b"chr1": [(50, 60)]}, {}],
                         ids=["chromosome_missing", "other_empty"])
def test_intersect_keeps_chromosome_only_in_self(other):
    # Surprising but documented: the intersect docstring says "For
    # chromosomes present only in this object, all regions are included
    # unchanged", so a chromosome absent from the other operand (or an
    # empty other operand) is passed through rather than dropped.
    a = _make_regions({b"chr1": [(0, 100)], b"chr2": [(0, 100)]})
    got = a.intersect(_make_regions(other))
    expected = {b"chr2": [(0, 100)]}
    expected[b"chr1"] = [(50, 60)] if b"chr1" in other else [(0, 100)]
    assert got.regions == expected
    assert got.total == 2


def test_intersect_empty_self():
    got = Regions().intersect(_make_regions({b"chr1": [(0, 10)]}))
    assert got.regions == {}
    assert got.total == 0


def test_intersect_returns_new_object_and_keeps_operands():
    a = _make_regions({b"chr1": [(50, 150), (0, 10)]})
    b = _make_regions({b"chr1": [(100, 200)]})
    got = a.intersect(b)
    assert got is not a and got is not b
    assert a.regions == {b"chr1": [(0, 10), (50, 150)]}
    assert b.regions == {b"chr1": [(100, 200)]}
    assert (a.total, b.total) == (2, 1)


def test_intersect_rejects_non_regions():
    with pytest.raises(AssertionError):
        _make_regions({b"chr1": [(0, 10)]}).intersect(PeakIO())


# ------------------------------------
# Regions.exclude
# ------------------------------------

EXCLUDE_CASES = [
    pytest.param([(0, 100), (200, 300)], [(50, 60)], id="inner_hit"),
    pytest.param([(0, 100)], [(100, 200)], id="touching_right_kept"),
    pytest.param([(100, 200)], [(0, 100)], id="touching_left_kept"),
    pytest.param([(0, 100), (150, 160), (300, 400)], [(90, 155)],
                 id="one_hits_two"),
    pytest.param([(0, 10), (20, 30)], [(40, 50)], id="no_overlap"),
    pytest.param([(0, 10), (20, 30)], [(0, 100)], id="all_removed"),
    pytest.param([(0, 1000), (2000, 2100)], [(10, 20), (30, 40)],
                 id="many_inside_one"),
    pytest.param([(300, 400), (0, 100), (150, 200)], [(350, 360)],
                 id="unsorted"),
    pytest.param([(0, 10), (500, 600), (700, 800)], [(5, 6)],
                 id="tail_after_other_exhausted"),
    pytest.param([(0, 100), (10, 20), (150, 160)], [(15, 16)],
                 id="nested_in_self"),
    pytest.param([(0, 10), (40, 50), (90, 100)], [(5, 45), (20, 30)],
                 id="nested_in_other"),
    pytest.param([(10, 20)], [(0, 5), (25, 30)], id="between_two"),
]


@pytest.mark.parametrize("a,b", EXCLUDE_CASES)
def test_exclude_matches_reference(a, b):
    # reference: keep a region iff it shares no base with the other set
    covered = _bases(b)
    expected = sorted(iv for iv in a if not (_bases([iv]) & covered))
    r = _make_regions({b"chr1": a})
    r.exclude(_make_regions({b"chr1": b}))
    assert r.regions[b"chr1"] == expected
    assert r.total == len(expected)


def test_exclude_removes_whole_regions():
    # The docstring summary says overlapping regions are removed; its
    # example output (trimmed regions) contradicts the summary.
    r = _make_regions({b"chr1": [(1000, 3000), (5000, 7000)]})
    r.exclude(_make_regions({b"chr1": [(2000, 6000)]}))
    assert r.regions == {b"chr1": []}
    assert r.total == 0


def test_exclude_chromosomes_in_one_operand_only():
    r = _make_regions({b"chr1": [(0, 10), (20, 30)], b"chr2": [(5, 9), (0, 3)]})
    r.exclude(_make_regions({b"chr1": [(25, 26)], b"chr3": [(0, 100)]}))
    assert r.regions == {b"chr1": [(0, 10)], b"chr2": [(0, 3), (5, 9)]}
    assert r.total == 3


def test_exclude_with_empty_other_keeps_all():
    r = _make_regions({b"chr1": [(20, 30), (0, 10)], b"chr2": [(0, 1)]})
    r.exclude(Regions())
    assert r.regions == {b"chr1": [(0, 10), (20, 30)], b"chr2": [(0, 1)]}
    assert r.total == 3


def test_exclude_on_empty_self():
    r = Regions()
    r.exclude(_make_regions({b"chr1": [(0, 10)]}))
    assert r.regions == {}
    assert r.total == 0


def test_exclude_rejects_non_regions():
    with pytest.raises(AssertionError):
        _make_regions({b"chr1": [(0, 10)]}).exclude(PeakIO())
