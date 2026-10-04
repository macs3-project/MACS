#!/usr/bin/env python

"""Module Description: Test the option validators in
MACS3/Utilities/OptValidator.py.

Options are built with the real ``bin/macs3`` argument parser and then
passed to the validator, as ``bin/macs3`` does. Where a branch can only
be reached with a value that argparse rejects (a format outside the
``choices``), the attribute is set after parsing.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import logging
import math
import os
import re
import sys

import pytest

from MACS3.IO.Parser import (BEDParser,
                             ELANDResultParser,
                             ELANDMultiParser,
                             ELANDExportParser,
                             SAMParser,
                             BAMParser,
                             BAMPEParser,
                             BEDPEParser,
                             BowtieParser,
                             FragParser,
                             guess_parser)
from MACS3.Utilities.OptValidator import (opt_validate_callpeak,
                                          opt_validate_filterdup,
                                          opt_validate_randsample,
                                          opt_validate_refinepeak,
                                          opt_validate_predictd,
                                          opt_validate_pileup,
                                          opt_validate_bdgcmp,
                                          opt_validate_cmbreps,
                                          opt_validate_bdgopt,
                                          opt_validate_callvar,
                                          opt_validate_hmmratac)

# ------------------------------------
# helpers
# ------------------------------------

LOGGER_NAME = "MACS3.Utilities.OptValidator"
LOGGER = logging.getLogger(LOGGER_NAME)
PREFIX = re.compile(r"^\[\d+ MB\] ")

# the smallest command line argparse accepts for each subcommand
MINIMAL_ARGS = {
    "callpeak": ["-t", "t.bed"],
    "filterdup": ["-i", "t.bed"],
    "randsample": ["-i", "t.bed", "-p", "10"],
    "refinepeak": ["-b", "p.bed", "-i", "t.bed", "-o", "o.bed"],
    "predictd": ["-i", "t.bed"],
    # -f BED: opt_validate_pileup rejects the default AUTO
    "pileup": ["-i", "t.bed", "-o", "o.bdg", "-f", "BED"],
    # -m ppois: opt_validate_bdgcmp rejects the default
    "bdgcmp": ["-t", "t.bdg", "-c", "c.bdg", "-m", "ppois",
               "--o-prefix", "x"],
    "cmbreps": ["-i", "a.bdg", "b.bdg", "-o", "o.bdg"],
    "bdgopt": ["-i", "a.bdg", "-o", "o.bdg"],
    "callvar": ["-b", "p.bed", "-t", "t.bam", "-o", "o.vcf"],
    "hmmratac": ["-i", "t.bam"],
}

VALIDATORS = {
    "callpeak": opt_validate_callpeak,
    "filterdup": opt_validate_filterdup,
    "randsample": opt_validate_randsample,
    "refinepeak": opt_validate_refinepeak,
    "predictd": opt_validate_predictd,
    "pileup": opt_validate_pileup,
    "bdgcmp": opt_validate_bdgcmp,
    "cmbreps": opt_validate_cmbreps,
    "bdgopt": opt_validate_bdgopt,
    "callvar": opt_validate_callvar,
    "hmmratac": opt_validate_hmmratac,
}

ALL_SUBCOMMANDS = list(VALIDATORS)

GSIZE_SHORTCUTS = "hs,mm,ce,dm"


@pytest.fixture(autouse=True)
def restore_logger_level():
    """The validators set the module logger level; put it back."""
    level = LOGGER.level
    yield
    LOGGER.setLevel(level)


@pytest.fixture
def parse_argv(macs3_argparser, monkeypatch):
    """``parse_argv([subcommand, ...])`` -> Namespace from the real
    parser, with ``sys.argv`` set to the same command line."""
    def _parse(argv):
        argv = [str(a) for a in argv]
        monkeypatch.setattr(sys, "argv", ["macs3"] + argv)
        return macs3_argparser.parse_args(argv)
    return _parse


@pytest.fixture
def parse(parse_argv):
    """``parse(subcommand, *extra)``: minimal arguments plus ``extra``."""
    def _parse(subcommand, *extra):
        return parse_argv([subcommand] + MINIMAL_ARGS[subcommand]
                          + list(extra))
    return _parse


def logged(caplog):
    """(level, message) of this module's records, memory prefix removed."""
    out = []
    for rec in caplog.records:
        if rec.name == LOGGER_NAME:
            msg = rec.getMessage()
            assert PREFIX.match(msg), msg
            out.append((rec.levelname, PREFIX.sub("", msg, count=1)))
    return out


def exit_messages(validator, options, caplog):
    """Run ``validator`` expecting sys.exit(1); return the log records."""
    with pytest.raises(SystemExit) as excinfo:
        validator(options)
    assert excinfo.value.code == 1
    return logged(caplog)


# ------------------------------------
# all validators
# ------------------------------------

@pytest.mark.parametrize("subcommand", ALL_SUBCOMMANDS)
def test_returns_same_namespace_with_logger_aliases(parse, subcommand):
    options = parse(subcommand)
    assert VALIDATORS[subcommand](options) is options
    assert options.error == LOGGER.critical
    assert options.warn == LOGGER.warning
    assert options.debug == LOGGER.debug
    assert options.info == LOGGER.info


@pytest.mark.parametrize("subcommand", ALL_SUBCOMMANDS)
def test_minimal_options_log_nothing(parse, caplog, subcommand):
    VALIDATORS[subcommand](parse(subcommand))
    assert logged(caplog) == []


@pytest.mark.parametrize("verbose, level", [(0, logging.ERROR),
                                            (3, logging.DEBUG)])
@pytest.mark.parametrize("subcommand", ALL_SUBCOMMANDS)
def test_verbose_sets_logger_level(parse, subcommand, verbose, level):
    VALIDATORS[subcommand](parse(subcommand, "--verbose", verbose))
    assert LOGGER.level == level


@pytest.mark.parametrize("verbose, level", [(0, 40), (1, 30), (2, 20),
                                            (3, 10), (4, 0), (5, -10)])
def test_verbose_level_formula(parse, verbose, level):
    # logger.setLevel((4 - verbose) * 10)
    opt_validate_callpeak(parse("callpeak", "--verbose", verbose))
    assert LOGGER.level == level


def test_default_verbose_is_info(parse):
    opt_validate_callpeak(parse("callpeak"))
    assert LOGGER.level == logging.INFO


def test_aliases_log_through_module_logger(parse, caplog):
    options = opt_validate_callpeak(parse("callpeak", "--verbose", 3))
    options.debug("d %d", 1)
    options.info("i")
    options.warn("w")
    options.error("e")
    assert logged(caplog) == [("DEBUG", "d 1"), ("INFO", "i"),
                              ("WARNING", "w"), ("CRITICAL", "e")]


@pytest.mark.parametrize("verbose, shown", [(0, False), (1, True)])
def test_warning_hidden_at_verbose_0(parse, caplog, verbose, shown):
    opt_validate_callpeak(parse("callpeak", "-f", "FRAG", "--verbose",
                                verbose))
    expected = [("WARNING", "Since the format is 'FRAG', `--keep-dup` "
                            "will be set as 'all'.")] if shown else []
    assert logged(caplog) == expected


# ------------------------------------
# --gsize (callpeak, filterdup, predictd)
# ------------------------------------

GSIZE_SUBCOMMANDS = ["callpeak", "filterdup", "predictd"]


@pytest.mark.parametrize("value, expected", [
    ("hs", 2913022398),
    ("mm", 2652783500),
    ("ce", 100286401),
    ("dm", 142573017),
    ("1e9", 1e9),
    ("2.7e9", 2.7e9),
    ("12345", 12345.0),
    ("1000000000", 1e9),
])
@pytest.mark.parametrize("subcommand", GSIZE_SUBCOMMANDS)
def test_gsize_shortcut_or_number(parse, subcommand, value, expected):
    options = VALIDATORS[subcommand](parse(subcommand, "-g", value))
    assert options.gsize == expected
    # shortcuts give the int from EFFECTIVEGS, numbers a float
    assert type(options.gsize) is type(expected)


@pytest.mark.parametrize("value", ["xx", "HS", "1e9bp"])
@pytest.mark.parametrize("subcommand", GSIZE_SUBCOMMANDS)
def test_gsize_invalid_exits(parse, caplog, subcommand, value):
    options = parse(subcommand, "-g", value)
    assert exit_messages(VALIDATORS[subcommand], options, caplog) == [
        ("ERROR", "Error when interpreting --gsize option: %s" % value),
        ("ERROR", "Available shortcuts of effective genome sizes are %s"
         % GSIZE_SHORTCUTS),
    ]


# ------------------------------------
# --format: parser class and gzip flag
# ------------------------------------

FORMAT_TABLE = {
    # format: (parser, gzip_flag)
    "AUTO": (guess_parser, False),
    "BAM": (BAMParser, True),
    "SAM": (SAMParser, False),
    "BED": (BEDParser, False),
    "ELAND": (ELANDResultParser, False),
    "ELANDMULTI": (ELANDMultiParser, False),
    "ELANDEXPORT": (ELANDExportParser, False),
    "BOWTIE": (BowtieParser, False),
    "BAMPE": (BAMPEParser, True),
    "BEDPE": (BEDPEParser, False),
    "FRAG": (FragParser, False),
}

SE_FORMATS = ["AUTO", "BAM", "SAM", "BED", "ELAND", "ELANDMULTI",
              "ELANDEXPORT", "BOWTIE"]

FORMATS_BY_SUBCOMMAND = {
    "callpeak": SE_FORMATS + ["BAMPE", "BEDPE", "FRAG"],
    "filterdup": SE_FORMATS + ["BAMPE", "BEDPE"],
    "randsample": SE_FORMATS + ["BAMPE", "BEDPE"],
    "refinepeak": SE_FORMATS,
    "predictd": SE_FORMATS + ["BAMPE", "BEDPE"],
    "pileup": SE_FORMATS[1:] + ["BAMPE", "BEDPE", "FRAG"],
}

FORMAT_CASES = [(sub, fmt) for sub, fmts in FORMATS_BY_SUBCOMMAND.items()
                for fmt in fmts]


@pytest.mark.parametrize("subcommand, fmt", FORMAT_CASES)
def test_format_selects_parser_and_gzip_flag(parse, subcommand, fmt):
    options = VALIDATORS[subcommand](parse(subcommand, "-f", fmt))
    parser, gzip_flag = FORMAT_TABLE[fmt]
    assert options.parser is parser
    assert options.gzip_flag is gzip_flag
    assert options.format == fmt


@pytest.mark.parametrize("subcommand", list(FORMATS_BY_SUBCOMMAND))
def test_format_is_uppercased(parse, subcommand):
    options = parse(subcommand)
    options.format = "bed"
    VALIDATORS[subcommand](options)
    assert options.format == "BED"
    assert options.parser is BEDParser


@pytest.mark.parametrize("subcommand, fmt", [
    ("callpeak", "xyz"),
    ("filterdup", "xyz"),
    ("filterdup", "FRAG"),
    ("randsample", "xyz"),
    ("randsample", "FRAG"),
    ("refinepeak", "xyz"),
    ("refinepeak", "BAMPE"),
    ("refinepeak", "BEDPE"),
    ("refinepeak", "FRAG"),
    ("predictd", "xyz"),
    ("predictd", "FRAG"),
    ("pileup", "xyz"),
])
def test_unrecognized_format_exits(parse, caplog, subcommand, fmt):
    options = parse(subcommand)
    options.format = fmt        # outside argparse's choices
    assert exit_messages(VALIDATORS[subcommand], options, caplog) == [
        ("ERROR", 'Format "%s" cannot be recognized!' % fmt.upper())]


@pytest.mark.parametrize("subcommand, fmt", [("refinepeak", "BAMPE"),
                                             ("filterdup", "FRAG"),
                                             ("callpeak", "bed")])
def test_argparse_rejects_formats_outside_choices(parse, subcommand, fmt):
    with pytest.raises(SystemExit) as excinfo:
        parse(subcommand, "-f", fmt)
    assert excinfo.value.code == 2


# ------------------------------------
# --keep-dup (callpeak, filterdup)
# ------------------------------------

@pytest.mark.parametrize("value", ["auto", "all", "0", "1", "5", "007"])
@pytest.mark.parametrize("subcommand", ["callpeak", "filterdup"])
def test_keepdup_valid_values_unchanged(parse, subcommand, value):
    options = VALIDATORS[subcommand](parse(subcommand, "--keep-dup", value))
    assert options.keepduplicates == value


@pytest.mark.parametrize("value", ["-1", "1.5", "AUTO", "x", ""])
@pytest.mark.parametrize("subcommand", ["callpeak", "filterdup"])
def test_keepdup_invalid_exits(parse, caplog, subcommand, value):
    options = parse(subcommand, "--keep-dup", value)
    assert exit_messages(VALIDATORS[subcommand], options, caplog) == [
        ("ERROR", "--keep-dup should be 'auto', 'all' or an integer!")]


@pytest.mark.parametrize("subcommand, default", [("callpeak", "1"),
                                                 ("filterdup", "auto")])
def test_keepdup_defaults(parse, subcommand, default):
    options = VALIDATORS[subcommand](parse(subcommand))
    assert options.keepduplicates == default


# ------------------------------------
# callpeak
# ------------------------------------

def test_callpeak_defaults(parse):
    o = opt_validate_callpeak(parse("callpeak"))
    assert o.gsize == 2913022398
    assert o.format == "AUTO"
    assert o.parser is guess_parser
    assert o.gzip_flag is False
    assert o.nomodel is False
    assert o.keepduplicates == "1"
    assert o.log_qvalue == pytest.approx(-math.log10(0.05), rel=1e-12)
    assert o.log_pvalue is None
    assert not hasattr(o, "log_broadcutoff")
    assert (o.lmfold, o.umfold) == (5, 50)
    assert o.shift == 0
    assert o.tsize is None
    assert o.cutoff_analysis_file == "None"


def test_callpeak_tsize_passes_through(parse):
    assert opt_validate_callpeak(parse("callpeak", "-s", 36)).tsize == 36


CALLPEAK_OUTPUTS = [
    ("peakxls", "_peaks.xls"),
    ("peakbed", "_peaks.bed"),
    ("peakNarrowPeak", "_peaks.narrowPeak"),
    ("peakBroadPeak", "_peaks.broadPeak"),
    ("peakGappedPeak", "_peaks.gappedPeak"),
    ("summitbed", "_summits.bed"),
    ("bdg_treat", "_treat_pileup.bdg"),
    ("bdg_control", "_control_lambda.bdg"),
    ("modelR", "_model.r"),
]


@pytest.mark.parametrize("outdir, name", [("", "NA"), ("out", "exp1"),
                                          ("/abs/dir", "x.y"),
                                          ("a/b/", "n")])
def test_callpeak_output_file_names(parse, outdir, name):
    extra = ["-n", name] + (["--outdir", outdir] if outdir else [])
    o = opt_validate_callpeak(parse("callpeak", *extra))
    for attr, suffix in CALLPEAK_OUTPUTS:
        assert getattr(o, attr) == os.path.join(outdir, name + suffix)


def test_callpeak_output_file_names_literal(parse):
    o = opt_validate_callpeak(parse("callpeak", "-n", "exp1", "--outdir",
                                    "out"))
    assert o.peakxls == "out/exp1_peaks.xls"
    assert o.peakNarrowPeak == "out/exp1_peaks.narrowPeak"
    assert o.summitbed == "out/exp1_summits.bed"
    assert o.bdg_treat == "out/exp1_treat_pileup.bdg"
    assert o.bdg_control == "out/exp1_control_lambda.bdg"
    assert o.modelR == "out/exp1_model.r"
    assert o.cutoff_analysis_file == "None"


@pytest.mark.parametrize("outdir, expected", [
    ("", "exp1_cutoff_analysis.txt"),
    ("out", "out/exp1_cutoff_analysis.txt"),
])
def test_callpeak_cutoff_analysis_file(parse, outdir, expected):
    extra = ["-n", "exp1", "--cutoff-analysis"]
    if outdir:
        extra += ["--outdir", outdir]
    o = opt_validate_callpeak(parse("callpeak", *extra))
    assert o.cutoff_analysis_file == expected


@pytest.mark.parametrize("fmt, pe", [("BAMPE", True), ("BEDPE", True),
                                     ("FRAG", True), ("BED", False),
                                     ("BAM", False), ("AUTO", False)])
def test_callpeak_pe_formats_force_nomodel_and_zero_shift(parse, fmt, pe):
    o = opt_validate_callpeak(parse("callpeak", "-f", fmt, "--shift", 50))
    assert o.nomodel is pe
    assert o.shift == (0 if pe else 50)
    assert ("# Paired-End mode is on\n" in o.argtxt) is pe
    assert ("# Paired-End mode is off\n" in o.argtxt) is (not pe)


def test_callpeak_nomodel_flag_kept_for_se(parse):
    assert opt_validate_callpeak(parse("callpeak", "--nomodel")).nomodel


def test_callpeak_frag_sets_keepdup_all_with_warning(parse, caplog):
    o = opt_validate_callpeak(parse("callpeak", "-f", "FRAG"))
    assert o.keepduplicates == "all"
    assert logged(caplog) == [
        ("WARNING",
         "Since the format is 'FRAG', `--keep-dup` will be set as 'all'.")]


@pytest.mark.parametrize("value", ["auto", "5"])
def test_callpeak_frag_overrides_any_keepdup(parse, caplog, value):
    o = opt_validate_callpeak(parse("callpeak", "-f", "FRAG", "--keep-dup",
                                    value))
    assert o.keepduplicates == "all"
    assert [lv for lv, _ in logged(caplog)] == ["WARNING"]


def test_callpeak_frag_keepdup_all_no_warning(parse, caplog):
    o = opt_validate_callpeak(parse("callpeak", "-f", "FRAG", "--keep-dup",
                                    "all"))
    assert o.keepduplicates == "all"
    assert logged(caplog) == []


def test_callpeak_frag_invalid_keepdup_exits_before_override(parse, caplog):
    options = parse("callpeak", "-f", "FRAG", "--keep-dup", "x")
    assert exit_messages(opt_validate_callpeak, options, caplog) == [
        ("ERROR", "--keep-dup should be 'auto', 'all' or an integer!")]


@pytest.mark.parametrize("subcommand", ["callpeak", "pileup"])
def test_frag_negative_max_count_exits(parse, caplog, subcommand):
    options = parse(subcommand, "-f", "FRAG", "--max-count", -1)
    assert exit_messages(VALIDATORS[subcommand], options, caplog) == [
        ("ERROR", "--max-count can't be a negative value")]


@pytest.mark.parametrize("subcommand", ["callpeak", "pileup"])
@pytest.mark.parametrize("value", [0, 1, 70000])
def test_frag_nonnegative_max_count_accepted(parse, subcommand, value):
    o = VALIDATORS[subcommand](parse(subcommand, "-f", "FRAG",
                                     "--max-count", value))
    assert o.maxcount == value


@pytest.mark.parametrize("subcommand", ["callpeak", "pileup"])
def test_max_count_only_checked_for_frag(parse, subcommand):
    o = VALIDATORS[subcommand](parse(subcommand, "-f", "BED",
                                     "--max-count", -1))
    assert o.maxcount == -1


@pytest.mark.parametrize("value", [0, -5])
def test_callpeak_extsize_below_1_exits(parse, caplog, value):
    options = parse("callpeak", "--extsize", value)
    assert exit_messages(opt_validate_callpeak, options, caplog) == [
        ("ERROR", "--extsize must >= 1!")]


def test_callpeak_extsize_1_accepted(parse):
    o = opt_validate_callpeak(parse("callpeak", "--extsize", 1))
    assert o.extsize == 1


def test_callpeak_broad_with_call_summits_exits(parse, caplog):
    options = parse("callpeak", "--broad", "--call-summits")
    assert exit_messages(opt_validate_callpeak, options, caplog) == [
        ("ERROR", "--broad can't be combined with --call-summits!")]


@pytest.mark.parametrize("p", [0.01, 1e-5, 0.5])
def test_callpeak_pvalue_sets_log_pvalue(parse, p):
    o = opt_validate_callpeak(parse("callpeak", "-p", p))
    assert o.log_pvalue == pytest.approx(-math.log10(p), rel=1e-12)
    assert o.log_qvalue is None


@pytest.mark.parametrize("q", [0.05, 0.01, 1.0])
def test_callpeak_qvalue_sets_log_qvalue(parse, q):
    o = opt_validate_callpeak(parse("callpeak", "-q", q))
    assert o.log_qvalue == pytest.approx(-math.log10(q), rel=1e-12,
                                         abs=1e-15)
    assert o.log_pvalue is None


def test_callpeak_pvalue_zero_falls_back_to_qvalue(parse):
    """`if options.pvalue` is false for 0, so -p 0 behaves as if -p were
    not given and the default q-value cutoff is used."""
    o = opt_validate_callpeak(parse("callpeak", "-p", 0))
    assert o.log_pvalue is None
    assert o.log_qvalue == pytest.approx(-math.log10(0.05), rel=1e-12)
    assert "# qvalue cutoff = 5.00e-02\n" in o.argtxt


def test_callpeak_p_and_q_are_mutually_exclusive(parse):
    with pytest.raises(SystemExit) as excinfo:
        parse("callpeak", "-p", 0.01, "-q", 0.05)
    assert excinfo.value.code == 2


@pytest.mark.parametrize("cutoff", [0.1, 0.05, 0.5])
def test_callpeak_broad_cutoff_log(parse, cutoff):
    o = opt_validate_callpeak(parse("callpeak", "--broad", "--broad-cutoff",
                                    cutoff))
    assert o.log_broadcutoff == pytest.approx(-math.log10(cutoff), rel=1e-12)


def test_callpeak_broad_cutoff_ignored_without_broad(parse):
    o = opt_validate_callpeak(parse("callpeak", "--broad-cutoff", 0.05))
    assert not hasattr(o, "log_broadcutoff")


@pytest.mark.parametrize("args", [["-q", "0"], ["-q", "-0.1"],
                                  ["-p", "-0.1"],
                                  ["--broad", "--broad-cutoff", "0"]])
def test_callpeak_nonpositive_cutoff_raises_valueerror(parse, args):
    """No check covers the range of -p/-q/--broad-cutoff; math.log raises
    ValueError for values <= 0 (worded differently from Python 3.14)."""
    with pytest.raises(ValueError, match="math domain error|expected a (positive|nonnegative) input"):
        opt_validate_callpeak(parse("callpeak", *args))


@pytest.mark.parametrize("subcommand", ["callpeak", "predictd"])
def test_d_min_zero_accepted(parse, subcommand):
    assert VALIDATORS[subcommand](parse(subcommand, "--d-min", 0)).d_min == 0


@pytest.mark.parametrize("subcommand", ["callpeak", "predictd"])
@pytest.mark.parametrize("low, high", [(10, 30), (7, 7), (1, 1000)])
def test_mfold_sets_lower_and_upper(parse, subcommand, low, high):
    o = VALIDATORS[subcommand](parse(subcommand, "-m", low, high))
    assert (o.lmfold, o.umfold) == (low, high)


@pytest.mark.parametrize("subcommand", ["callpeak", "predictd"])
def test_mfold_lower_above_upper_exits(parse, caplog, subcommand):
    # The message is %-formatted with the mfold list; str % list treats
    # the list as a mapping, so no TypeError (unlike --d-min's int).
    options = parse(subcommand, "-m", 50, 5)
    assert exit_messages(VALIDATORS[subcommand], options, caplog) == [
        ("ERROR", "Upper limit of mfold should be greater than lower "
                  "limit!")]


# ------------------------------------
# callpeak: argtxt
# ------------------------------------

def test_callpeak_argtxt_defaults(parse):
    o = opt_validate_callpeak(parse("callpeak"))
    assert o.argtxt == (
        "# Command line: callpeak -t t.bed\n"
        "# ARGUMENTS LIST:\n"
        "# name = NA\n"
        "# format = AUTO\n"
        "# ChIP-seq file = ['t.bed']\n"
        "# control file = None\n"
        "# effective genome size = 2.91e+09\n"
        "# band width = 300\n"
        "# model fold = [5, 50]\n"
        "# qvalue cutoff = 5.00e-02\n"
        "# The maximum gap between significant sites is assigned as the "
        "read length/tag size.\n"
        "# The minimum length of peaks is assigned as the predicted "
        "fragment length \"d\".\n"
        "# Larger dataset will be scaled towards smaller dataset.\n"
        "# Range for calculating regional lambda is: 10000 bps\n"
        "# Broad region calling is off\n"
        "# Paired-End mode is off\n")


def test_callpeak_argtxt_command_line_from_sys_argv(parse):
    o = opt_validate_callpeak(parse("callpeak", "-q", "0.01", "-n", "x"))
    assert o.argtxt.splitlines()[0] == ("# Command line: callpeak -t t.bed "
                                        "-q 0.01 -n x")


QVALUE_NOT_CALCULATED = ("# qvalue will not be calculated and reported as "
                         "-1 in the final output.\n")

ARGTXT_CASES = [
    # extra arguments, lines present, substrings absent
    (["-p", "0.01"],
     ["# pvalue cutoff = 1.00e-02\n", QVALUE_NOT_CALCULATED],
     ["# qvalue cutoff"]),
    (["-p", "0.01", "--broad"],
     ["# pvalue cutoff for narrow/strong regions = 1.00e-02\n",
      "# pvalue cutoff for broad/weak regions = 1.00e-01\n",
      QVALUE_NOT_CALCULATED, "# Broad region calling is on\n"],
     ["# pvalue cutoff = "]),
    (["--broad"],
     ["# qvalue cutoff for narrow/strong regions = 5.00e-02\n",
      "# qvalue cutoff for broad/weak regions = 1.00e-01\n",
      "# Broad region calling is on\n"],
     ["# qvalue cutoff = ", "# Broad region calling is off",
      "qvalue will not be calculated"]),
    (["--broad", "--broad-cutoff", "0.2", "-q", "0.01"],
     ["# qvalue cutoff for narrow/strong regions = 1.00e-02\n",
      "# qvalue cutoff for broad/weak regions = 2.00e-01\n"],
     []),
    (["--max-gap", "30"],
     ["# The maximum gap between significant sites = 30\n"],
     ["assigned as the read length"]),
    (["--min-length", "100"],
     ["# The minimum length of peaks = 100\n"],
     ["predicted fragment length"]),
    (["--down-sample"],
     ["# Larger dataset will be randomly sampled towards smaller "
      "dataset.\n"],
     ["# Random seed", "scaled towards"]),
    (["--down-sample", "--seed", "0"],
     ["# Random seed has been set as: 0\n"],
     []),
    (["--down-sample", "--seed", "12"],
     ["# Random seed has been set as: 12\n"],
     []),
    (["--seed", "12"],
     [],
     ["# Random seed"]),
    (["--scale-to", "large"],
     ["# Smaller dataset will be scaled towards larger dataset.\n"],
     ["# Larger dataset will be scaled"]),
    (["--scale-to", "small"],
     ["# Larger dataset will be scaled towards smaller dataset.\n"],
     ["# Smaller dataset"]),
    (["--ratio", "0.5"],
     ["# Using a custom scaling factor: 5.00e-01\n"],
     []),
    (["--ratio", "1.0"],
     [],
     ["custom scaling factor"]),
    (["-c", "c.bed"],
     ["# control file = ['c.bed']\n",
      "# Range for calculating regional lambda is: 1000 bps and 10000 "
      "bps\n"],
     []),
    (["-c", "c.bed", "--slocal", "500", "--llocal", "20000"],
     ["# Range for calculating regional lambda is: 500 bps and 20000 "
      "bps\n"],
     []),
    (["--llocal", "5000"],
     ["# Range for calculating regional lambda is: 5000 bps\n"],
     []),
    (["-c"],
     ["# control file = []\n",
      "# Range for calculating regional lambda is: 10000 bps\n"],
     []),
    (["--fe-cutoff", "2"],
     ["# Additional cutoff on fold-enrichment is: 2.00\n"],
     []),
    (["--call-summits"],
     ["# Searching for subpeak summits is on\n"],
     []),
    (["--SPMR", "-B"],
     ["# MACS will save fragment pileup signal per million reads\n"],
     []),
    (["--SPMR"],
     [],
     ["per million reads"]),
    (["-f", "FRAG", "--max-count", "3"],
     ["# Maximum count in fragment file is set as 3\n",
      "# Paired-End mode is on\n"],
     []),
    (["-f", "FRAG", "--max-count", "0"],
     [],
     ["# Maximum count"]),
    (["-f", "BED", "--max-count", "3"],
     [],
     ["# Maximum count"]),
    (["-g", "1e6"],
     ["# effective genome size = 1.00e+06\n"],
     []),
    (["-g", "dm"],
     ["# effective genome size = 1.43e+08\n"],
     []),
    (["-n", "myexp", "-f", "BED"],
     ["# name = myexp\n", "# format = BED\n"],
     []),
    (["--bw", "150"],
     ["# band width = 150\n"],
     []),
    (["-m", "3", "30"],
     ["# model fold = [3, 30]\n"],
     []),
    (["-t", "a.bed", "b.bed"],
     ["# ChIP-seq file = ['a.bed', 'b.bed']\n"],
     []),
]


@pytest.mark.parametrize("extra, present, absent", ARGTXT_CASES,
                         ids=[" ".join(c[0]) for c in ARGTXT_CASES])
def test_callpeak_argtxt_lines(parse, extra, present, absent):
    o = opt_validate_callpeak(parse("callpeak", *extra))
    for line in present:
        assert line in o.argtxt
    for text in absent:
        assert text not in o.argtxt
    assert o.argtxt.endswith("\n")


# ------------------------------------
# filterdup
# ------------------------------------

def test_filterdup_defaults(parse):
    o = opt_validate_filterdup(parse("filterdup"))
    assert o.gsize == 2913022398
    assert o.format == "AUTO"
    assert o.parser is guess_parser
    assert o.gzip_flag is False
    assert o.keepduplicates == "auto"
    assert not hasattr(o, "nomodel")


# ------------------------------------
# randsample
# ------------------------------------

@pytest.mark.parametrize("value", [100, 50, 0.001])
def test_randsample_percentage_accepted(parse, value):
    o = opt_validate_randsample(parse("randsample", "-p", value))
    assert o.percentage == value


@pytest.mark.parametrize("value", [100.5, 1000])
def test_randsample_percentage_above_100_exits(parse, caplog, value):
    options = parse("randsample", "-p", value)
    assert exit_messages(opt_validate_randsample, options, caplog) == [
        ("ERROR", "Percentage can't be bigger than 100.0. Please check your "
                  "options and retry!")]


@pytest.mark.parametrize("value", [1, 1e6])
def test_randsample_number_accepted(parse_argv, value):
    o = opt_validate_randsample(parse_argv(["randsample", "-i", "t.bed",
                                            "-n", value]))
    assert o.number == value
    assert o.percentage is None


@pytest.mark.parametrize("value", [-1, -0.5])
def test_randsample_negative_number_exits(parse_argv, caplog, value):
    options = parse_argv(["randsample", "-i", "t.bed", "-n", value])
    assert exit_messages(opt_validate_randsample, options, caplog) == [
        ("ERROR", "Number of tags can't be smaller than or equal to 0. "
                  "Please check your options and retry!")]


def test_randsample_needs_percentage_or_number(parse_argv):
    with pytest.raises(SystemExit) as excinfo:
        parse_argv(["randsample", "-i", "t.bed"])
    assert excinfo.value.code == 2


def test_randsample_defaults(parse):
    o = opt_validate_randsample(parse("randsample"))
    assert o.format == "AUTO"
    assert o.parser is guess_parser
    assert o.gzip_flag is False
    assert not hasattr(o, "gsize")


# ------------------------------------
# refinepeak
# ------------------------------------

def test_refinepeak_defaults_and_outputs_untouched(parse):
    o = opt_validate_refinepeak(parse("refinepeak"))
    assert o.parser is guess_parser
    assert o.gzip_flag is False
    assert o.ofile == "o.bed"
    assert o.oprefix is None
    assert (o.cutoff, o.windowsize) == (5, 200)


# ------------------------------------
# predictd
# ------------------------------------

@pytest.mark.parametrize("outdir, rfile, expected", [
    ("", None, "predictd_model.R"),
    ("out", None, "out/predictd_model.R"),
    ("", "m.R", "m.R"),
    ("/tmp/x", "m.R", "/tmp/x/m.R"),
])
def test_predictd_model_r_path(parse, outdir, rfile, expected):
    extra = (["--outdir", outdir] if outdir else []) + (
        ["--rfile", rfile] if rfile else [])
    o = opt_validate_predictd(parse("predictd", *extra))
    assert o.modelR == expected


@pytest.mark.parametrize("fmt", ["BAMPE", "BEDPE"])
def test_predictd_pe_formats_set_nomodel(parse, fmt):
    assert opt_validate_predictd(parse("predictd", "-f", fmt)).nomodel is True


@pytest.mark.parametrize("fmt", ["BED", "BAM", "AUTO"])
def test_predictd_se_formats_leave_nomodel_unset(parse, fmt):
    # predictd has no --nomodel option; only PE formats add the attribute
    assert not hasattr(opt_validate_predictd(parse("predictd", "-f", fmt)),
                       "nomodel")


def test_predictd_defaults(parse):
    o = opt_validate_predictd(parse("predictd"))
    assert o.gsize == 2913022398
    assert (o.lmfold, o.umfold) == (5, 50)
    assert o.d_min == 20


# ------------------------------------
# pileup
# ------------------------------------

@pytest.mark.parametrize("value", [0, -1])
def test_pileup_extsize_not_positive_exits(parse, caplog, value):
    options = parse("pileup", "--extsize", value)
    assert exit_messages(opt_validate_pileup, options, caplog) == [
        ("ERROR", "--extsize must > 0!")]


@pytest.mark.parametrize("value", [1, 200])
def test_pileup_extsize_positive_accepted(parse, value):
    assert opt_validate_pileup(parse("pileup", "--extsize",
                                     value)).extsize == value


def test_pileup_frag_does_not_force_nomodel_or_keepdup(parse):
    o = opt_validate_pileup(parse("pileup", "-f", "FRAG"))
    assert o.parser is FragParser
    assert not hasattr(o, "nomodel")
    assert not hasattr(o, "keepduplicates")


# ------------------------------------
# bdgcmp
# ------------------------------------

ALL_BDGCMP_METHODS = ["ppois", "qpois", "subtract", "logFE", "FE", "logLR",
                      "slogLR", "max"]


@pytest.mark.parametrize("methods", [["ppois"], ["max"], ALL_BDGCMP_METHODS,
                                     ["FE", "FE"]])
def test_bdgcmp_valid_methods(parse, methods):
    o = opt_validate_bdgcmp(parse("bdgcmp", "-m", *methods))
    assert o.method == methods


def test_bdgcmp_invalid_method_exits(parse, caplog):
    options = parse("bdgcmp")
    options.method = ["ppois", "bogus"]     # outside argparse's choices
    assert exit_messages(opt_validate_bdgcmp, options, caplog) == [
        ("ERROR", "Invalid method: bogus")]


@pytest.mark.parametrize("methods, ofiles", [
    (["ppois", "FE"], ["a.bdg"]),
    (["ppois"], ["a.bdg", "b.bdg"]),
    (["FE", "FE"], ["a.bdg"]),          # counted with repeats
])
def test_bdgcmp_ofile_count_must_match_methods(parse_argv, caplog, methods,
                                               ofiles):
    options = parse_argv(["bdgcmp", "-t", "t.bdg", "-c", "c.bdg", "-m"]
                         + methods + ["-o"] + ofiles)
    assert exit_messages(opt_validate_bdgcmp, options, caplog) == [
        ("ERROR", "The number and the order of arguments for --ofile must "
                  "be the same as for -m.")]


def test_bdgcmp_matching_ofiles_accepted(parse_argv):
    o = opt_validate_bdgcmp(parse_argv(["bdgcmp", "-t", "t.bdg", "-c",
                                        "c.bdg", "-m", "ppois", "FE", "-o",
                                        "a.bdg", "b.bdg"]))
    assert o.ofile == ["a.bdg", "b.bdg"]


def test_bdgcmp_needs_ofile_or_oprefix(parse_argv):
    with pytest.raises(SystemExit) as excinfo:
        parse_argv(["bdgcmp", "-t", "t.bdg", "-c", "c.bdg"])
    assert excinfo.value.code == 2


# ------------------------------------
# cmbreps
# ------------------------------------

@pytest.mark.parametrize("method", ["fisher", "max", "mean"])
def test_cmbreps_methods_accepted(parse, method):
    o = opt_validate_cmbreps(parse("cmbreps", "-m", method))
    assert o.method == method


def test_cmbreps_default_method_fisher(parse):
    assert opt_validate_cmbreps(parse("cmbreps")).method == "fisher"


def test_cmbreps_invalid_method_exits(parse, caplog):
    options = parse("cmbreps")
    options.method = "median"           # outside argparse's choices
    assert exit_messages(opt_validate_cmbreps, options, caplog) == [
        ("ERROR", "Invalid method: median")]


def test_cmbreps_single_replicate_exits(parse_argv, caplog):
    options = parse_argv(["cmbreps", "-i", "a.bdg", "-o", "o.bdg"])
    assert exit_messages(opt_validate_cmbreps, options, caplog) == [
        ("ERROR", "Combining replicates needs at least two replicates!")]


def test_cmbreps_method_checked_before_replicate_count(parse_argv, caplog):
    options = parse_argv(["cmbreps", "-i", "a.bdg", "-o", "o.bdg"])
    options.method = "median"
    assert exit_messages(opt_validate_cmbreps, options, caplog) == [
        ("ERROR", "Invalid method: median")]


def test_cmbreps_three_replicates_accepted(parse_argv):
    o = opt_validate_cmbreps(parse_argv(["cmbreps", "-i", "a", "b", "c",
                                         "-o", "o.bdg"]))
    assert o.ifile == ["a", "b", "c"]


# ------------------------------------
# bdgopt
# ------------------------------------

@pytest.mark.parametrize("method", ["p2q", "max", "min"])
def test_bdgopt_methods_without_extra_param(parse, method):
    o = opt_validate_bdgopt(parse("bdgopt", "-m", method))
    assert o.method == method
    assert o.extraparam is None


@pytest.mark.parametrize("method", ["multiply", "add"])
def test_bdgopt_multiply_add_with_extra_param(parse, method):
    o = opt_validate_bdgopt(parse("bdgopt", "-m", method, "-p", 2))
    assert o.extraparam == [2.0]


@pytest.mark.parametrize("method", ["multiply", "add"])
@pytest.mark.parametrize("extra", [[], ["-p"]])
def test_bdgopt_multiply_add_need_extra_param(parse, caplog, method, extra):
    options = parse("bdgopt", "-m", method, *extra)
    assert exit_messages(opt_validate_bdgopt, options, caplog) == [
        ("ERROR", "Need EXTRAPARAM for method multiply or add!")]


def test_bdgopt_invalid_method_exits(parse, caplog):
    options = parse("bdgopt")
    options.method = "Divide"           # outside argparse's choices
    assert exit_messages(opt_validate_bdgopt, options, caplog) == [
        ("ERROR", "Invalid method: Divide")]


@pytest.mark.parametrize("method", ["MULTIPLY", "Add", "P2Q"])
def test_bdgopt_method_check_is_case_insensitive(parse, method):
    options = parse("bdgopt", "-p", 1)
    options.method = method             # outside argparse's choices
    assert opt_validate_bdgopt(options).method == method


def test_bdgopt_default_method_p2q(parse):
    assert opt_validate_bdgopt(parse("bdgopt")).method == "p2q"


# ------------------------------------
# callvar
# ------------------------------------

@pytest.mark.parametrize("value, expected", [(0, 1), (-2, 1), (1, 1),
                                             (4, 4)])
def test_callvar_np_at_least_1(parse, value, expected):
    assert opt_validate_callvar(parse("callvar", "-m", value)).np == expected


def test_callvar_np_default_1(parse):
    assert opt_validate_callvar(parse("callvar")).np == 1


@pytest.mark.parametrize("extra, expected", [([], "auto"),
                                             (["-F", "on"], "on"),
                                             (["--fermi", "off"], "off")])
def test_callvar_fermi_passes_through(parse, extra, expected):
    assert opt_validate_callvar(parse("callvar", *extra)).fermi == expected


# ------------------------------------
# hmmratac
# ------------------------------------

HMM_TYPE_LINE = ("# Use --hmm-type to select a Gaussian ('gaussian') or "
                 "Poisson ('poisson') model for the hidden markov model in "
                 "HMMRATAC. Default: 'gaussian'. \n")


def test_hmmratac_argtxt_defaults(parse):
    o = opt_validate_hmmratac(parse("hmmratac"))
    assert o.argtxt == ("# Command line: hmmratac -i t.bam\n"
                        "# Random seed selected as: 10151\n"
                        + HMM_TYPE_LINE)


@pytest.mark.parametrize("extra, line", [
    (["--no-fragem"],
     "# EM training not performed on fragment distribution. \n"),
    (["-t", "train.bed"],
     "# Using -t, --training input to train HMM instead of using fold "
     "change settings to select. \n"),
    (["--randomSeed", "12345"],
     "# Random seed selected as: 12345\n"),
    (["--modelonly"],
     "# Program will stop after generating model, which can be later "
     "applied with '--model'. \n"),
    (["--hmm-type", "poisson"], HMM_TYPE_LINE),
])
def test_hmmratac_argtxt_lines(parse, extra, line):
    o = opt_validate_hmmratac(parse("hmmratac", *extra))
    assert line in o.argtxt


def test_hmmratac_argtxt_order(parse):
    o = opt_validate_hmmratac(parse("hmmratac", "--no-fragem", "-t", "r.bed",
                                    "--modelonly"))
    assert o.argtxt == (
        "# Command line: hmmratac -i t.bam --no-fragem -t r.bed "
        "--modelonly\n"
        "# EM training not performed on fragment distribution. \n"
        "# Using -t, --training input to train HMM instead of using fold "
        "change settings to select. \n"
        "# Random seed selected as: 10151\n"
        "# Program will stop after generating model, which can be later "
        "applied with '--model'. \n"
        + HMM_TYPE_LINE)


def test_hmmratac_defaults_unchanged(parse):
    o = opt_validate_hmmratac(parse("hmmratac"))
    assert o.em_means == [50, 200, 400, 600]
    assert o.em_stddevs == [20, 20, 20, 20]
    assert (o.hmm_lower, o.hmm_upper) == (10, 20)
    assert o.min_frag_p == 0.001
    assert o.prescan_cutoff == 1.2
    assert o.format == "BAMPE"
    assert not hasattr(o, "parser")


@pytest.mark.parametrize("extra, message", [
    (["--means", "-1", "200", "400", "600"],
     " `--means` should not be negative! "),
    (["--means", "50", "200", "400", "-0.5"],
     " `--means` should not be negative! "),
    (["--stddevs", "20", "-20", "20", "20"],
     " `--stddev` should not be negative! "),
    (["--min-frag-p", "0"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--min-frag-p", "1"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--min-frag-p", "-0.1"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--min-frag-p", "1.5"],
     " `--min-frag-p` should be larger than 0 and smaller than 1!"),
    (["--binsize", "0"], " `--binsize` must be larger than 0."),
    (["--binsize", "-10"], " `--binsize` must be larger than 0."),
    (["-l", "-1"], " `-l` or `--lower` should not be negative! "),
    (["-u", "-1"], " `-u` or `--upper` should not be negative! "),
    (["--maxTrain", "0"], " `--maxTrain` should be larger than 0!"),
    (["--maxTrain", "-5"], " `--maxTrain` should be larger than 0!"),
    (["-c", "1"],
     " In order to use -c or --prescan-cutoff, the cutoff must be larger "
     "than 1."),
    (["-c", "0.5"],
     " In order to use -c or --prescan-cutoff, the cutoff must be larger "
     "than 1."),
    (["--minlen", "-1"],
     " In order to use --minlen, the length should not be negative."),
])
def test_hmmratac_invalid_values_exit(parse, caplog, extra, message):
    options = parse("hmmratac", *extra)
    assert exit_messages(opt_validate_hmmratac, options, caplog) == [
        ("ERROR", message)]


@pytest.mark.parametrize("extra", [
    ["--means", "0", "0", "0", "0"],
    ["--stddevs", "0", "1", "2", "3"],
    ["--min-frag-p", "0.5"],
    ["--min-frag-p", "0.999"],
    ["--binsize", "1"],
    ["-l", "0", "-u", "0"],
    ["-l", "20", "-u", "20"],
    ["--maxTrain", "1"],
    ["-c", "1.01"],
    ["--minlen", "0"],
])
def test_hmmratac_boundary_values_accepted(parse, caplog, extra):
    opt_validate_hmmratac(parse("hmmratac", *extra))
    assert logged(caplog) == []


def test_hmmratac_means_checked_before_stddevs(parse, caplog):
    options = parse("hmmratac", "--means", "-1", "1", "1", "1",
                    "--stddevs", "-1", "1", "1", "1")
    assert exit_messages(opt_validate_hmmratac, options, caplog) == [
        ("ERROR", " `--means` should not be negative! ")]


