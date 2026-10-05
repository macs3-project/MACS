#!/usr/bin/env python

"""Module Description: Test MACS3.Commands.callvar_cmd and the
``macs3 callvar`` command: check_names, run, call_variants_at_range,
every command-line option and the option validation.

Inputs are small sorted and indexed BAM files written with pysam whose
reads carry MD tags computed against a deterministic pseudo-random
reference. Expected VCF records are derived from the read layout (allele
counts, strands, depths) and from direct numpy/scipy implementations of
the genotype likelihoods (see test_VariantStat.py); the genotype rule
(choose the model with the lowest BIC by a margin of 2) is written out in
``ref_call`` below.

The C-only helpers are covered here through the command as a whole:
ReadAlignment.get_n_edits and relative_ref_pos_to_relative_query_pos
(through RACollection and PosReadsInfo construction), the
GreedyMaxFunction* and calculate_ln functions (through call_GT), and the
fermi-lite/Smith-Waterman methods of RACollection (through ``-F on``).

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import datetime
import math
import re
import signal
import sys

import numpy as np
import pytest
from scipy.stats import binom

from MACS3.Commands.callvar_cmd import (check_names,
                                        run,
                                        call_variants_at_range)
from MACS3.IO.PeakIO import (PeakIO)
from MACS3.Signal.ReadAlignment import (ReadAlignment)
from MACS3.Signal.RACollection import (RACollection)
from MACS3.Signal.PeakVariants import (Variant)
from MACS3.Utilities.Constants import (MACS_VERSION)

LN10 = math.log(10)
CIGAR_OPS = "MIDNSHP=X"
BAM_CODES = "=ACMGRSVTWYHKDBN"
VCF_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"


# ------------------------------------
# synthetic genome and reads
# ------------------------------------

def lcg_seq(n, seed=20261003):
    """Deterministic pseudo-random DNA (linear congruential generator)."""
    x = seed
    out = []
    for _ in range(n):
        x = (1103515245 * x + 12345) % 2**31
        out.append("ACGT"[(x >> 16) & 3])
    return "".join(out)


GENOME = {"chr1": lcg_seq(3000), "chr2": lcg_seq(1000, seed=7)}
REFS = (("chr1", 5000), ("chr2", 2000))


def cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def md_tag(cigar, seq, start, ref):
    md, run_, r, q = "", 0, start, 0
    for n, op in cigar_ops(cigar):
        if op in "M=X":
            for _ in range(n):
                if seq[q] == ref[r]:
                    run_ += 1
                else:
                    md += "%d%s" % (run_, ref[r])
                    run_ = 0
                q += 1
                r += 1
        elif op == "D":
            md += "%d^%s" % (run_, ref[r:r + n])
            run_ = 0
            r += n
        elif op in "IS":
            q += n
    return md + str(run_)


def other_base(b, k=1):
    return "ACGT"[("ACGT".index(b) + k) % 4]


def snv_reads(site, alleles, starts, strands, site_quals, prefix,
              chrom="chr1", readlen=50, md=True):
    """Reads of length ``readlen`` at ``starts`` carrying ``alleles`` at
    ``site``; every base has quality 40 except the site."""
    ref = GENOME[chrom]
    reads = []
    for j, (a, s, strand, q) in enumerate(zip(alleles, starts, strands,
                                              site_quals)):
        seq = list(ref[s:s + readlen])
        seq[site - s] = a
        seq = "".join(seq)
        qual = ["I"] * readlen
        qual[site - s] = chr(33 + q)
        cigar = "%dM" % readlen
        tags = [("MD", md_tag(cigar, seq, s, ref))] if md else []
        reads.append(dict(name="%s%d" % (prefix, j), ref=chrom, pos=s,
                          flag=16 if strand else 0, cigar=cigar, seq=seq,
                          qual="".join(qual), tags=tags))
    return reads


REF1 = GENOME["chr1"]

# heterozygous SNV at chr1:300 (0-based): 5 reference reads (quality 40
# at the site) and 5 alternative reads (quality 35), none at a read end;
# strands (j // 2) % 2 give 3 plus and 2 minus reads for each allele
HET = 300
HET_REF = REF1[HET]
HET_ALT = other_base(HET_REF)
HET_STARTS = [HET - 45 + 4 * j for j in range(10)]
HET_ALLELES = [HET_REF if j % 2 == 0 else HET_ALT for j in range(10)]
HET_STRANDS = [(j // 2) % 2 for j in range(10)]
HET_QUALS = [40 if j % 2 == 0 else 35 for j in range(10)]

# homozygous SNV at chr1:700: 8 alternative reads at quality 40
HOM = 700
HOM_REF = REF1[HOM]
HOM_ALT = other_base(HOM_REF, 2)
HOM_STARTS = [HOM - 40 + 3 * j for j in range(8)]

PEAKS = [("chr1", 200, 400), ("chr1", 600, 800)]


def het_reads(**kw):
    return snv_reads(HET, HET_ALLELES, HET_STARTS, HET_STRANDS, HET_QUALS,
                     "het", **kw)


def hom_reads(**kw):
    return snv_reads(HOM, [HOM_ALT] * 8, HOM_STARTS, [j % 2 for j in range(8)],
                     [40] * 8, "hom", **kw)


def control_reads():
    # 2 reference and 1 alternative read over the heterozygous site
    return snv_reads(HET, [HET_REF, HET_REF, HET_ALT],
                     [HET - 20, HET - 15, HET - 10], [0, 1, 0], [40, 40, 40],
                     "ctl")


def deletion_site(lo):
    """First position >= lo whose base differs from both neighbours (one
    alignment for a 1-bp deletion) and whose two upstream bases differ (so
    an anchor taken one base too far left is visible)."""
    d = lo
    while (REF1[d - 1] == REF1[d] or REF1[d] == REF1[d + 1]
           or REF1[d - 2] == REF1[d - 1]):
        d += 1
    return d


DEL = deletion_site(1100)


def deletion_reads(n=8):
    """Reads with the 1-bp deletion of chr1:DEL, which lies 20, 22, ...,
    34 bases into the reads."""
    reads = []
    for j in range(n):
        s = DEL - 20 - 2 * j
        k = DEL - s
        seq = REF1[s:DEL] + REF1[DEL + 1:s + 51]
        cigar = "%dM1D%dM" % (k, 50 - k)
        reads.append(dict(name="del%d" % j, ref="chr1", pos=s,
                          flag=16 if j % 2 else 0, cigar=cigar, seq=seq,
                          qual="I" * 50,
                          tags=[("MD", md_tag(cigar, seq, s, REF1))]))
    return reads


# ------------------------------------
# reference genotype calls
# ------------------------------------

def err_rate(quals):
    return 10.0 ** (-np.asarray(quals, dtype=float) / 10.0)


def ref_homo(t1T, t1C, t2T, t2C):
    e1 = err_rate(list(t1T) + list(t1C))
    e2 = err_rate(list(t2T) + list(t2C))
    lnL = float(np.log1p(-e1).sum() + np.log(e2).sum())
    return lnL, -2 * lnL


def ref_lnL_k(me, ne, k, p):
    tn = len(me) + len(ne)
    f = k / tn
    em, en = err_rate(me), err_rate(ne)
    return float(binom.logpmf(k, tn, p)
                 + np.log((1 - em) * f + em * (1 - f)).sum()
                 + np.log((1 - en) * (1 - f) + en * f).sum())


def ref_noAS_max(me, ne):
    # MACS3 fixes k at the observed count when one allele has no read
    if len(me) == 0 or len(ne) == 0:
        return ref_lnL_k(me, ne, len(me), 0.5)
    return max(ref_lnL_k(me, ne, k, 0.5) for k in range(len(me) + len(ne) + 1))


def ref_AS_max(me, ne, max_ar):
    tn = len(me) + len(ne)
    mar = float(np.float32(max_ar))
    if len(ne) == 0:
        return ref_lnL_k(me, ne, tn, mar)
    return max(ref_lnL_k(me, ne, k, min(max(k / tn, 1 - mar), mar))
               for k in range(tn + 1))


def phred(lnL):
    return -10.0 * lnL / LN10


def f32(x):
    return float(np.float32(x))


def ref_call(chrom, pos, ref, top1, top2, b1T, b2T, b1C=(), b2C=(),
             strands=(0, 0, 0, 0), other_T=0, other_C=0, max_ar=0.95,
             return_GQ=False):
    """Expected VCF line (list of fields) of one site, or None when the
    site is filtered. b1*/b2* are the base qualities of the top1/top2
    alleles in treatment (T) and control (C); strands = (plus top1, plus
    top2, minus top1, minus top2); other_T/other_C count reads of further
    alleles, which enter the depths but not the likelihoods."""
    tn_T, tn_C = len(b1T) + len(b2T), len(b1C) + len(b2C)
    hM = ref_homo(b1T, b1C, b2T, b2C)
    hm = ref_homo(b2T, b2C, b1T, b1C)
    nA = ref_noAS_max(b1T, b2T)
    A = ref_AS_max(b1T, b2T, max_ar)
    pen_n, pen_a = math.log(tn_T), 2 * math.log(tn_T)
    if tn_C:
        c = ref_noAS_max(b1C, b2C)
        nA += c
        A += c
        pen_n += math.log(tn_C)
        pen_a += math.log(tn_C)
    nA = (nA, -2 * nA + pen_n)
    A = (A, -2 * A + pen_a)
    B = {"hM": hM[1], "hm": hm[1], "nA": nA[1], "A": A[1]}
    if top1 != ref and len(b2T) + len(b2C) == 0:
        dBIC = min(B["nA"], B["A"], B["hm"]) - B["hM"]
        if dBIC < 2:
            return None
        mtype, GT, alt = "homo", "1/1", top1
        PL00 = max(0, phred(hm[0]) - phred(hM[0]))
        PL01 = max(0, phred(max(nA[0], A[0])) - phred(hM[0]))
        PL11 = 0
        GQ = min(PL00, PL01)
    elif top1 != ref and B["hM"] + 2 <= min(B["hm"], B["nA"], B["A"]):
        dBIC = min(B["nA"], B["A"], B["hm"]) - B["hM"]
        mtype, GT, alt = "homo", "1/1", top1
        PL00 = phred(hm[0]) - phred(hM[0])
        PL01 = phred(max(nA[0], A[0])) - phred(hM[0])
        PL11 = 0
        GQ = min(PL00, PL01)
    else:
        if B["nA"] + 2 <= min(B["hM"], B["hm"], B["A"]):
            mtype, lnL01, dBIC = "heter_noAS", nA[0], min(B["hM"], B["hm"]) - B["nA"]
        elif B["A"] + 2 <= min(B["hM"], B["hm"], B["nA"]):
            mtype, lnL01, dBIC = "heter_AS", A[0], min(B["hM"], B["hm"]) - B["A"]
        elif B["A"] + 2 <= B["hM"] and B["A"] + 2 <= B["hm"]:
            mtype, lnL01 = "heter_unsure", max(nA[0], A[0])
            dBIC = min(B["hM"], B["hm"]) - max(B["A"], B["nA"])
        else:
            return None
        PL01 = 0
        PL00 = phred(hm[0]) - phred(lnL01)
        PL11 = phred(hM[0]) - phred(lnL01)
        GQ = min(PL00, PL11)
        if ref == top1:
            GT, alt = "0/1", top2
        elif ref == top2:
            GT, alt = "0/1", top1
        else:
            GT, alt = "1/2", top1 + "," + top2
    mt = ",".join("Deletion" if a == "*" else "Insertion" if len(a) > 1
                  else "SNV" for a in alt.split(","))
    n1T, n2T, n1C, n2C = len(b1T), len(b2T), len(b1C), len(b2C)
    info = ("M=%s;MT=%s;DPT=%d;DPC=%d;DP1T=%d%s;DP2T=%d%s;DP1C=%d%s;DP2C=%d%s;"
            "SB=%d,%d,%d,%d;DBIC=%.2f;BICHOMOMAJOR=%.2f;BICHOMOMINOR=%.2f;"
            "BICHETERNOAS=%.2f;BICHETERAS=%.2f;AR=%.2f"
            % ((mtype, mt, tn_T + other_T, tn_C + other_C, n1T, top1, n2T, top2,
                n1C, top1, n2C, top2) + tuple(strands)
               + (f32(dBIC), f32(B["hM"]), f32(B["hm"]), f32(B["nA"]),
                  f32(B["A"]), f32(n1T / (n1T + n2T)))))
    DP = tn_T + tn_C + other_T + other_C
    sample = "%s:%d:%d:%d,%d,%d" % (GT, DP, int(GQ), int(PL00), int(PL01),
                                    int(PL11))
    fields = [chrom, str(pos + 1), ".", ref, alt, "%d" % int(GQ), ".", info,
              "GT:DP:GQ:PL", sample]
    return (fields, GQ) if return_GQ else fields


def het_expected(control=False, max_ar=0.95):
    kw = {}
    if control:
        kw = dict(b1C=[40, 40], b2C=[40])
    return ref_call("chr1", HET, HET_REF, HET_REF, HET_ALT, [40] * 5, [35] * 5,
                    strands=(3, 3, 2, 2), max_ar=max_ar, **kw)


def hom_expected(max_ar=0.95, return_GQ=False):
    return ref_call("chr1", HOM, HOM_REF, HOM_ALT, HOM_REF, [40] * 8, [],
                    strands=(4, 0, 4, 0), max_ar=max_ar, return_GQ=return_GQ)


# ------------------------------------
# running callvar
# ------------------------------------

def write_peaks(tmp_path, peaks=PEAKS, name="peaks.bed"):
    path = tmp_path / name
    path.write_text("".join("%s\t%d\t%d\n" % p for p in peaks))
    return str(path)


def write_bam(make_alignments, reads, name="treat.bam", refs=REFS):
    return make_alignments(reads, refs=refs, name=name, index=True)


def records(path):
    with open(path) as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#")]


def header(path):
    with open(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.startswith("#")]


def info(rec):
    return dict(x.split("=", 1) for x in rec[7].split(";"))


def callvar(run_macs3, tmp_path, treat, peaks, *extra, out="out.vcf",
            timeout=60):
    out = str(tmp_path / out)
    args = ["callvar", "-b", peaks, "-t", treat, "-o", out] + [str(x) for x in extra]
    res = run_macs3(args, timeout=timeout)
    return res, out, args


@pytest.fixture
def snv_data(tmp_path, make_alignments):
    """Treatment BAM with the heterozygous and homozygous SNVs, a control
    BAM and the two-peak BED file."""
    treat = write_bam(make_alignments, het_reads() + hom_reads())
    ctrl = write_bam(make_alignments, control_reads(), name="ctrl.bam")
    return treat, ctrl, write_peaks(tmp_path)


# ------------------------------------
# check_names
# ------------------------------------

class Track:
    def __init__(self, names):
        self.names = names

    def get_chr_names(self):
        return set(self.names)


def test_check_names_common_names():
    messages = []
    assert check_names(Track([b"chr1", b"chr2"]), Track([b"chr2", b"chr3"]),
                       messages.append) is None
    assert messages == []


def test_check_names_no_common_names_exits():
    messages = []
    with pytest.raises(SystemExit) as excinfo:
        check_names(Track(["chrB", "chrA"]), Track(["1", "2"]), messages.append)
    assert excinfo.value.code is None
    assert messages == [
        "No common chromosome names can be found from treatment and control! "
        "Check your input files! MACS will quit...",
        "Chromosome names in treatment: chrA,chrB",
        "Chromosome names in control: 1,2"]


# ------------------------------------
# call_variants_at_range
# ------------------------------------

def ra_from_dict(r):
    """ReadAlignment from a make_alignments read dict, as BAMaccessor
    builds it."""
    seq = r["seq"]
    codes = [BAM_CODES.index(c) for c in seq] + [0] * (len(seq) % 2)
    packed = bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))
    ops = cigar_ops(r["cigar"])
    span = sum(n for n, op in ops if op in "MDN=X")
    return ReadAlignment(r["name"].encode(), r["ref"].encode(), r["pos"],
                         r["pos"] + span, 1 if r["flag"] & 16 else 0, packed,
                         bytes(ord(c) - 33 for c in r["qual"]),
                         tuple((n << 4) | CIGAR_OPS.index(op) for n, op in ops),
                         dict(r["tags"])["MD"].encode())


def het_collection(control=False):
    peaks = PeakIO()
    peaks.add(b"chr1", 200, 400)
    peak = peaks.get_data_from_chrom(b"chr1")[0]
    reads_C = [ra_from_dict(r) for r in control_reads()] if control else []
    return RACollection(b"chr1", peak, [ra_from_dict(r) for r in het_reads()],
                        reads_C)


CVR_ARGS = dict(top2allelesminr=0.8, max_allowed_ar=0.95, min_altallele_count=2,
                min_homo_GQ=0, min_heter_GQ=0, minQ=20)


def test_call_variants_at_range_finds_het_snv():
    c = het_collection()
    result = call_variants_at_range((c["left"], c["right"]), s=c["peak_refseq"],
                                    collection=c, **CVR_ARGS)
    assert [p for p, v in result] == [HET]
    assert isinstance(result[0][1], Variant)
    assert result[0][1].toVCF().split("\t") == het_expected()[3:]


def test_call_variants_at_range_with_control():
    c = het_collection(control=True)
    result = call_variants_at_range((HET - 5, HET + 5), s=c["peak_refseq"],
                                    collection=c, **CVR_ARGS)
    assert [p for p, v in result] == [HET]
    assert result[0][1].toVCF().split("\t") == het_expected(control=True)[3:]


@pytest.mark.parametrize("lr", [(200, HET), (HET + 1, 400), (HET, HET)])
def test_call_variants_at_range_outside_site(lr):
    c = het_collection()
    assert call_variants_at_range(lr, s=c["peak_refseq"], collection=c,
                                  **CVR_ARGS) == []


def test_call_variants_at_range_skips_N_reference():
    c = het_collection()
    s = bytearray(c["peak_refseq"])
    s[HET - c["left"]] = ord("N")
    assert call_variants_at_range((HET, HET + 1), s=bytes(s), collection=c,
                                  **CVR_ARGS) == []


@pytest.mark.parametrize("change,kept", [
    (dict(minQ=35), False),                 # alternative bases have quality 35
    (dict(min_altallele_count=6), False),
    (dict(min_heter_GQ=10**6), False),
    (dict(top2allelesminr=1.0), True),
])
def test_call_variants_at_range_cutoffs(change, kept):
    c = het_collection()
    kw = dict(CVR_ARGS)
    kw.update(change)
    result = call_variants_at_range((HET, HET + 1), s=c["peak_refseq"],
                                    collection=c, **kw)
    assert len(result) == int(kept)


# ------------------------------------
# run (in process)
# ------------------------------------

def test_run_in_process(snv_data, tmp_path, macs3_argparser, monkeypatch):
    treat, ctrl, peaks = snv_data
    out = str(tmp_path / "inproc.vcf")
    argv = ["callvar", "-b", peaks, "-t", treat, "-o", out]
    monkeypatch.setattr(sys, "argv", ["macs3"] + argv)
    assert run(macs3_argparser.parse_args(argv)) is None
    assert records(out) == [het_expected(), hom_expected()]


# ------------------------------------
# macs3 callvar: default run, -b, -t, -o, header
# ------------------------------------

def test_callvar_default_output(snv_data, tmp_path, run_macs3, parse_log):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks)
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected(), hom_expected()]
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert log == ["Peak: chr1 200 400", " Call variants w/o assembly",
                   "Peak: chr1 600 800", " Call variants w/o assembly"]


def test_callvar_header(snv_data, tmp_path, run_macs3, test_dir):
    """Exact VCF header.

    Pins the current output. The ##INFO/##FORMAT description lines are
    free text, so they are compared with the header of the upstream
    reference file test/standard_results_callvar/PEsample.vcf. The option
    value of --altallele-count is recorded as '--top2allele-count'; it
    is normalised here.
    """
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks)
    assert res.returncode == 0, res.stderr
    # the INFO/FORMAT lines are those of the upstream reference output
    std = header(test_dir / "standard_results_callvar" / "PEsample.vcf")
    std_info = [ln for ln in std if ln.startswith(("##INFO", "##FORMAT"))]
    program_args = " ".join(args[1:] + [
        "-Q", "20", "-D", "1", "--max-ar", "0.95", "--top2alleles-mratio", "0.8",
        "--altallele-count", "2", "-g", "0", "-G", "0",
        " --fermi auto --fermi-overlap 30"])
    got = header(out)
    got[3] = got[3].replace(" --top2allele-count ", " --altallele-count ")
    assert got == (
        ["##fileformat=VCFv4.1",
         "##fileDate=%s" % datetime.date.today().strftime("%Y%m%d"),
         "##source=MACS_V%s" % MACS_VERSION,
         "##Program_Args=callvar " + program_args]
        + std_info
        + ["##contig=<ID=chr1,length=5000,assembly=NA>",
           "##contig=<ID=chr2,length=2000,assembly=NA>",
           VCF_COLUMNS])


def chr2_het_reads(site=800):
    """The heterozygous layout copied to chr2:site."""
    ref2 = GENOME["chr2"]
    alleles = [ref2[site] if j % 2 == 0 else other_base(ref2[site])
               for j in range(10)]
    starts = [site - 45 + 4 * j for j in range(10)]
    reads = snv_reads(site, alleles, starts, HET_STRANDS, HET_QUALS, "c2het",
                      chrom="chr2")
    expected = ref_call("chr2", site, ref2[site], ref2[site], alleles[1],
                        [40] * 5, [35] * 5, strands=(3, 3, 2, 2))
    return reads, expected


def test_callvar_peaks_order_and_empty_peaks(tmp_path, make_alignments, run_macs3,
                                             parse_log):
    reads2, exp2 = chr2_het_reads()
    treat = write_bam(make_alignments, het_reads() + reads2)
    # unsorted peaks, a peak without reads, a chromosome missing from the BAM
    peaks = write_peaks(tmp_path, [("chr2", 700, 900), ("chr1", 1500, 1600),
                                   ("chr1", 200, 400), ("chrX", 0, 100)])
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks)
    assert res.returncode == 0, res.stderr
    recs = records(out)
    assert [(r[0], r[1]) for r in recs] == [("chr1", "301"), ("chr2", "801")]
    assert recs == [het_expected(), exp2]
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert log == ["Peak: chr1 200 400", " Call variants w/o assembly",
                   "Peak: chr1 1500 1600", "No reads found in this peak. Skipped",
                   "Peak: chr2 700 900", " Call variants w/o assembly"]


def test_callvar_narrowpeak_input(snv_data, tmp_path, run_macs3):
    treat, ctrl, _ = snv_data
    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text("chr1\t200\t400\tpeak1\t100\t.\t5.0\t10.0\t8.0\t100\n")
    res, out, args = callvar(run_macs3, tmp_path, treat, str(peaks))
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected()]


def test_callvar_no_peaks(snv_data, tmp_path, run_macs3):
    treat, ctrl, _ = snv_data
    peaks = tmp_path / "empty.bed"
    peaks.write_text("")
    res, out, args = callvar(run_macs3, tmp_path, treat, str(peaks))
    assert res.returncode == 0, res.stderr
    assert records(out) == []
    assert header(out)[-1] == VCF_COLUMNS


def test_callvar_missing_bai(tmp_path, make_alignments, run_macs3):
    treat = make_alignments(het_reads(), refs=REFS, name="noindex.bam",
                            index=False)
    res, out, args = callvar(run_macs3, tmp_path, treat, write_peaks(tmp_path))
    assert res.returncode == 1
    assert ("BAI is not available! Please make sure the `%s.bai` file exists "
            "in the same path" % treat) in res.stderr


def test_callvar_missing_peak_file(snv_data, tmp_path, run_macs3):
    treat, ctrl, _ = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat,
                             str(tmp_path / "nope.bed"))
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr


def test_callvar_output_in_missing_directory(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks,
                             out="no_such_dir/out.vcf")
    assert res.returncode == 1
    assert "FileNotFoundError" in res.stderr


# ------------------------------------
# -c/--control
# ------------------------------------

def test_callvar_control(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-c", ctrl)
    assert res.returncode == 0, res.stderr
    recs = records(out)
    assert recs == [het_expected(control=True), hom_expected()]
    i = info(recs[0])
    assert (i["DPC"], i["DP1C"], i["DP2C"]) == ("3", "2" + HET_REF, "1" + HET_ALT)
    assert recs[0][9].split(":")[1] == "13"


def test_callvar_control_with_other_chromosome_names(snv_data, tmp_path,
                                                     make_alignments, run_macs3):
    treat, _, peaks = snv_data
    reads = [dict(r, ref="1") for r in control_reads()]
    ctrl = write_bam(make_alignments, reads, name="ctrl_other.bam",
                     refs=(("1", 5000),))
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-c", ctrl)
    assert res.returncode == 1
    assert ("It seems Treatment and Control BAM use different naming for "
            "chromosomes! Check headers of both files.") in res.stderr


# ------------------------------------
# --outdir, --verbose
# ------------------------------------

def test_callvar_outdir_is_created(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    res = run_macs3(["callvar", "-b", peaks, "-t", treat, "-o", "rel.vcf",
                     "--outdir", tmp_path / "outdir"], timeout=60)
    assert res.returncode == 0, res.stderr
    assert (tmp_path / "outdir").is_dir()
    assert records(tmp_path / "rel.vcf") == [het_expected(), hom_expected()]


@pytest.mark.parametrize("verbose,n_info", [(0, 0), (1, 0), (2, 4), (3, 4)])
def test_callvar_verbose(snv_data, tmp_path, run_macs3, parse_log, verbose,
                         n_info):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks,
                             "--verbose", verbose)
    assert res.returncode == 0, res.stderr
    log = parse_log(res.stderr)
    assert sum(1 for lvl, m in log if lvl == "INFO") == n_info
    assert all(lvl in ("INFO", "") for lvl, m in log)
    assert records(out) == [het_expected(), hom_expected()]


# ------------------------------------
# -g/--gq-hetero, -G/--gq-homo
# ------------------------------------

def test_callvar_gq_hetero(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    GQ = int(het_expected()[5])
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-g", GQ)
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected(), hom_expected()]
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "--gq-hetero",
                             GQ + 1, out="out2.vcf")
    assert res.returncode == 0, res.stderr
    # the homozygous call is not affected by the heterozygous cutoff
    assert records(out) == [hom_expected()]


def test_callvar_gq_homo(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    GQ = int(hom_expected()[5])
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-G", GQ)
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected(), hom_expected()]
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "--gq-homo",
                             GQ + 1, out="out2.vcf")
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected()]
    assert "-G %s" % float(GQ + 1) in header(out)[3]


# ------------------------------------
# -Q, -D
# ------------------------------------

@pytest.mark.parametrize("Q,het_kept", [(20, True), (34, True), (35, False)])
def test_callvar_Q(snv_data, tmp_path, run_macs3, Q, het_kept):
    # alternative bases at the heterozygous site have quality 35
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-Q", Q)
    assert res.returncode == 0, res.stderr
    expected = ([het_expected()] if het_kept else []) + [hom_expected()]
    assert records(out) == expected


def dup_reads():
    # 4 distinct reference reads and 3 identical alternative reads
    ref_reads = snv_reads(HET, [HET_REF] * 4, [HET - 40, HET - 30, HET - 20,
                                               HET - 12], [0, 1, 0, 1],
                          [40] * 4, "ref")
    alt_reads = snv_reads(HET, [HET_ALT] * 3, [HET - 25] * 3, [0, 0, 0],
                          [40] * 3, "dup")
    return ref_reads + alt_reads


@pytest.mark.parametrize("D,n_alt", [(1, 1), (2, 2), (3, 3), (10, 3)])
def test_callvar_max_duplicate(tmp_path, make_alignments, run_macs3, D, n_alt):
    treat = write_bam(make_alignments, dup_reads())
    res, out, args = callvar(run_macs3, tmp_path, treat, write_peaks(tmp_path),
                             "-D", D)
    assert res.returncode == 0, res.stderr
    expected = ref_call("chr1", HET, HET_REF, HET_REF, HET_ALT, [40] * 4,
                        [40] * n_alt, strands=(2, n_alt, 2, 0))
    # with one alternative read the site fails --altallele-count
    if n_alt >= 2:
        assert expected is not None
    assert records(out) == ([expected] if n_alt >= 2 else [])


# ------------------------------------
# --top2alleles-mratio, --altallele-count, --max-ar
# ------------------------------------

def triallelic_reads():
    alleles = ([HET_REF] * 4 + [other_base(HET_REF, 1)] * 3
               + [other_base(HET_REF, 2)] * 3)
    starts = [HET - 45 + 4 * j for j in range(10)]
    return snv_reads(HET, alleles, starts, [0] * 10, [40] * 10, "tri")


@pytest.mark.parametrize("ratio,kept", [(None, False), (0.8, False), (0.75, False),
                                        (0.65, True)])
def test_callvar_top2alleles_mratio(tmp_path, make_alignments, run_macs3, ratio,
                                    kept):
    # top1 (4) + top2 (3) = 7 of 10 reads at the site
    treat = write_bam(make_alignments, triallelic_reads())
    extra = [] if ratio is None else ["--top2alleles-mratio", ratio]
    res, out, args = callvar(run_macs3, tmp_path, treat, write_peaks(tmp_path),
                             *extra)
    assert res.returncode == 0, res.stderr
    alt1 = other_base(HET_REF, 1)
    alt2 = other_base(HET_REF, 2)
    # top2 is the first of the two tied alternative alleles in ACGT order
    top2 = min(alt1, alt2, key="ACGT".index)
    expected = ref_call("chr1", HET, HET_REF, HET_REF, top2, [40] * 4, [40] * 3,
                        strands=(4, 3, 0, 0), other_T=3)
    assert records(out) == ([expected] if kept else [])


def two_alt_reads():
    alleles = [HET_REF] * 6 + [HET_ALT] * 2
    starts = [HET - 45 + 5 * j for j in range(8)]
    return snv_reads(HET, alleles, starts, [j % 2 for j in range(8)], [40] * 8,
                     "aa")


@pytest.mark.parametrize("count,kept", [(None, True), (2, True), (3, False)])
def test_callvar_altallele_count(tmp_path, make_alignments, run_macs3, count, kept):
    treat = write_bam(make_alignments, two_alt_reads())
    extra = [] if count is None else ["--altallele-count", count]
    res, out, args = callvar(run_macs3, tmp_path, treat, write_peaks(tmp_path),
                             *extra)
    assert res.returncode == 0, res.stderr
    expected = ref_call("chr1", HET, HET_REF, HET_REF, HET_ALT, [40] * 6,
                        [40] * 2, strands=(3, 1, 3, 1))
    assert expected is not None
    assert records(out) == ([expected] if kept else [])


@pytest.mark.parametrize("max_ar", [0.8, 0.99])
def test_callvar_max_ar_changes_homozygous_PL(snv_data, tmp_path, run_macs3,
                                              max_ar):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "--max-ar", max_ar)
    assert res.returncode == 0, res.stderr
    recs = records(out)
    assert recs[-1] == hom_expected(max_ar=max_ar)
    assert recs[-1] != hom_expected()


def test_callvar_max_ar_drops_minor_allele(tmp_path, make_alignments, run_macs3):
    # 6 reference + 2 alternative reads: 6/8 = 0.75 > 0.7, alt is dropped
    treat = write_bam(make_alignments, two_alt_reads())
    res, out, args = callvar(run_macs3, tmp_path, treat, write_peaks(tmp_path),
                             "--max-ar", 0.7)
    assert res.returncode == 0, res.stderr
    assert records(out) == []


# ------------------------------------
# -m/--multiple-processing
# ------------------------------------

@pytest.mark.parametrize("np_", [0, -2])
def test_callvar_np_below_one_runs_single_process(snv_data, tmp_path, run_macs3,
                                                  np_):
    # opt_validate_callvar turns -m <= 0 into 1
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-m", np_)
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected(), hom_expected()]


def test_callvar_two_processes(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks,
                             "--multiple-processing", 2)
    assert res.returncode == 0, res.stderr
    assert records(out) == [het_expected(), hom_expected()]


# ------------------------------------
# -F/--fermi, --fermi-overlap
# ------------------------------------

def deletion_record_without_alleles(rec):
    """A deletion record minus REF, ALT and the allele labels of DP1T..DP2C
    (those are not compared)."""
    i = info(rec)
    for key in ("DP1T", "DP2T", "DP1C", "DP2C"):
        i[key] = re.match(r"\d+", i[key]).group(0)
    return rec[:3] + rec[5:7] + [i] + rec[8:]


def expected_deletion_record():
    # the deleted base is the '*' allele of all 8 reads at quality 93;
    # fix_indels moves the record to the preceding base (POS = DEL, 1-based)
    exp = ref_call("chr1", DEL - 1, REF1[DEL], "*", REF1[DEL], [93] * 8, [],
                   strands=(4, 0, 4, 0))
    return deletion_record_without_alleles(exp)


def test_callvar_fermi_off_deletion(tmp_path, make_alignments, run_macs3,
                                    parse_log):
    treat = write_bam(make_alignments, deletion_reads())
    peaks = write_peaks(tmp_path, [("chr1", 1000, 1200)])
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-F", "off")
    assert res.returncode == 0, res.stderr
    recs = records(out)
    assert len(recs) == 1
    assert deletion_record_without_alleles(recs[0]) == expected_deletion_record()
    assert info(recs[0])["MT"] == "Deletion"
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert log == ["Peak: chr1 1000 1200", " Call variants w/o assembly"]


def test_callvar_fermi_auto_skips_assembly_without_indel(snv_data, tmp_path,
                                                         run_macs3, parse_log):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "--fermi", "auto")
    assert res.returncode == 0, res.stderr
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert " Try to call variants w/ fermi-lite assembly" not in log
    assert records(out) == [het_expected(), hom_expected()]


def test_callvar_fermi_on_assembles(snv_data, tmp_path, run_macs3, parse_log):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-F", "on")
    assert res.returncode == 0, res.stderr
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert log.count(" Try to call variants w/ fermi-lite assembly") == 2
    assert " Call variants w/o assembly" not in log
    assert "--fermi on --fermi-overlap 30" in header(out)[3]
    # every read maps back to an assembled haplotype: same calls as without
    # assembly
    assert records(out) == [het_expected(), hom_expected()]


def test_callvar_fermi_auto_assembles_with_indel(tmp_path, make_alignments,
                                                 run_macs3, parse_log):
    treat = write_bam(make_alignments, deletion_reads())
    peaks = write_peaks(tmp_path, [("chr1", 1000, 1200)])
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks)
    assert res.returncode == 0, res.stderr
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert log == ["Peak: chr1 1000 1200", " Call variants w/o assembly",
                   " Try to call variants w/ fermi-lite assembly"]
    # every read spans the deletion, so the unitig-based call is the same
    recs = records(out)
    assert len(recs) == 1
    assert deletion_record_without_alleles(recs[0]) == expected_deletion_record()


def test_callvar_fermi_overlap_in_header(snv_data, tmp_path, run_macs3):
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks,
                             "--fermi-overlap", 25)
    assert res.returncode == 0, res.stderr
    assert header(out)[3].endswith(" --fermi auto --fermi-overlap 25")
    # no indel: the overlap only matters once fermi-lite runs
    assert records(out) == [het_expected(), hom_expected()]


def test_callvar_fermi_overlap_longer_than_reads_aborts(snv_data, tmp_path,
                                                       run_macs3):
    """--fermi-overlap 60 with 50 bp reads makes fermi-lite abort callvar
    with an assertion in mr_insert_multi (mrope.c).

    Pins the current output. The --fermi-overlap help requires a value
    between 1 and the read length; outside that range the behaviour is not
    documented.
    """
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-F", "on",
                             "--fermi-overlap", 60)
    assert res.returncode == -signal.SIGABRT
    assert "mr_insert_multi: Assertion" in res.stderr
    assert records(out) == []


def test_callvar_fermi_unknown_value_calls_nothing(snv_data, tmp_path,
                                                   run_macs3, parse_log):
    """-F accepts any string; a value other than auto, on and off runs
    neither the calling pass without assembly nor the assembly, so every
    peak is skipped and the VCF holds only its header.

    Pins the current output. The help describes auto, on and off and says
    nothing about other values (the option has no choices list).
    """
    treat, ctrl, peaks = snv_data
    res, out, args = callvar(run_macs3, tmp_path, treat, peaks, "-F", "maybe")
    assert res.returncode == 0, res.stderr
    assert records(out) == []
    assert header(out)[3].endswith(" --fermi maybe --fermi-overlap 30")
    log = [m for lvl, m in parse_log(res.stderr) if lvl == "INFO"]
    assert log == ["Peak: chr1 200 400", "Peak: chr1 600 800"]


# ------------------------------------
# argument errors (argparse exits with 2)
# ------------------------------------

@pytest.mark.parametrize("drop,msg", [
    ("-b", "the following arguments are required: -b/--peak"),
    ("-t", "the following arguments are required: -t/--treatment"),
    ("-o", "the following arguments are required: -o/--ofile"),
])
def test_callvar_required_arguments(tmp_path, run_macs3, drop, msg):
    args = {"-b": "p.bed", "-t": "t.bam", "-o": str(tmp_path / "o.vcf")}
    del args[drop]
    res = run_macs3(["callvar"] + [x for kv in args.items() for x in kv],
                    timeout=30)
    assert res.returncode == 2
    assert msg in res.stderr


@pytest.mark.parametrize("opt,value,msg", [
    ("-Q", "abc", "argument -Q: invalid int value: 'abc'"),
    ("-D", "1.5", "argument -D: invalid int value: '1.5'"),
    ("-g", "x", "argument -g/--gq-hetero: invalid float value: 'x'"),
    ("-G", "x", "argument -G/--gq-homo: invalid float value: 'x'"),
    ("--fermi-overlap", "x", "argument --fermi-overlap: invalid int value: 'x'"),
    ("--top2alleles-mratio", "x",
     "argument --top2alleles-mratio: invalid float value: 'x'"),
    ("--altallele-count", "x", "argument --altallele-count: invalid int value: 'x'"),
    ("--max-ar", "x", "argument --max-ar: invalid float value: 'x'"),
    ("-m", "x", "argument -m/--multiple-processing: invalid int value: 'x'"),
    ("--verbose", "x", "argument --verbose: invalid int value: 'x'"),
])
def test_callvar_invalid_values(tmp_path, run_macs3, opt, value, msg):
    res = run_macs3(["callvar", "-b", "p.bed", "-t", "t.bam", "-o",
                     str(tmp_path / "o.vcf"), opt, value], timeout=30)
    assert res.returncode == 2
    assert msg in res.stderr


def test_callvar_parser_defaults(macs3_argparser):
    a = macs3_argparser.parse_args(["callvar", "-b", "p", "-t", "t", "-o", "o"])
    assert (a.peakbed, a.tfile, a.cfile, a.outdir, a.ofile) == ("p", "t", None,
                                                                "", "o")
    assert (a.verbose, a.GQCutoffHetero, a.GQCutoffHomo, a.Q, a.maxDuplicate) == (
        2, 0, 0, 20, 1)
    assert (a.fermi, a.fermiMinOverlap, a.top2allelesMinRatio,
            a.altalleleMinCount, a.maxAR, a.np) == ("auto", 30, 0.8, 2, 0.95, 1)


# ------------------------------------
# upstream test data
# ------------------------------------

@pytest.mark.slow
def test_callvar_matches_upstream_standard(test_dir, tmp_path, run_macs3):
    """Same command as test/cmdlinetest; compares sorted records with
    test/standard_results_callvar/PEsample.vcf (upstream output).

    Pins the current output. Deriving genotype calls for real ChIP-seq
    reads with fermi-lite assembly by hand is not practical.
    """
    out = tmp_path / "PEsample.vcf"
    res = run_macs3(["callvar", "-b", test_dir / "callvar_testing.narrowPeak",
                     "-t", test_dir / "CTCF_PE_ChIP_chr22_50k.bam",
                     "-c", test_dir / "CTCF_PE_CTRL_chr22_50k.bam",
                     "-o", out], timeout=300)
    assert res.returncode == 0, res.stderr
    std = test_dir / "standard_results_callvar" / "PEsample.vcf"
    assert sorted(map(tuple, records(out))) == sorted(map(tuple, records(std)))


def test_callvar_matches_upstream_standard_two_peaks(test_dir, tmp_path, run_macs3):
    """Peaks 3 (insertions, fermi-lite assembly) and 7a (control reads) of
    test/callvar_testing.narrowPeak; the records must equal those of
    test/standard_results_callvar/PEsample.vcf inside these peaks.

    Pins the current output. Deriving genotype calls for real ChIP-seq
    reads with fermi-lite assembly by hand is not practical.
    """
    lines = (test_dir / "callvar_testing.narrowPeak").read_text().splitlines()
    chosen = [lines[2], lines[6]]
    peaks = tmp_path / "two.narrowPeak"
    peaks.write_text("\n".join(chosen) + "\n")
    out = tmp_path / "two.vcf"
    res = run_macs3(["callvar", "-b", peaks,
                     "-t", test_dir / "CTCF_PE_ChIP_chr22_50k.bam",
                     "-c", test_dir / "CTCF_PE_CTRL_chr22_50k.bam",
                     "-o", out], timeout=120)
    assert res.returncode == 0, res.stderr
    regions = [(f[0], int(f[1]), int(f[2])) for f in (ln.split("\t") for ln in chosen)]
    std = records(test_dir / "standard_results_callvar" / "PEsample.vcf")
    expected = []
    for rec in std:
        inside = any(rec[0] == c and s < int(rec[1]) <= e for c, s, e in regions)
        if inside and rec not in expected:
            expected.append(rec)
    assert len(expected) == 8
    assert records(out) == expected
