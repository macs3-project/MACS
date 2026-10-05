#!/usr/bin/env python

"""Module Description: Test the constants in
MACS3/Utilities/Constants.py and the modules that rely on them.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import re

import pytest

from MACS3.Utilities import Constants
from MACS3.Utilities.Constants import (MACS_VERSION,
                                       MAX_PAIRNUM,
                                       MAX_LAMBDA,
                                       FESTEP,
                                       BUFFER_SIZE,
                                       READ_BUFFER_SIZE,
                                       N_MP,
                                       EFFECTIVEGS)

# ------------------------------------
# helpers
# ------------------------------------

BUFFERED_SUBCOMMANDS = ["callpeak", "filterdup", "predictd", "pileup",
                        "randsample", "refinepeak", "hmmratac"]


def subparser(argparser, name):
    """The argparse sub-parser of ``macs3 <name>``."""
    for action in argparser._subparsers._group_actions:
        if name in action.choices:
            return action.choices[name]
    raise KeyError(name)


def option(parser, dest):
    for action in parser._actions:
        if action.dest == dest:
            return action
    raise KeyError(dest)


# ------------------------------------
# values
# ------------------------------------

@pytest.mark.parametrize("name, value", [
    # 3.0.4 -> 3.0.5: upstream bumped the version in 6cefb42 and d8918ba
    # (3.0.5 release, #744).
    ("MACS_VERSION", "3.0.5"),
    ("MAX_PAIRNUM", 1000),
    ("MAX_LAMBDA", 100000),
    ("FESTEP", 20),
    ("BUFFER_SIZE", 100000),
    ("READ_BUFFER_SIZE", 10000000),
    ("N_MP", 2),
])
def test_constant_values(name, value):
    got = getattr(Constants, name)
    assert type(got) is type(value)
    assert got == value


def test_imported_names_are_the_module_attributes():
    assert (MACS_VERSION, MAX_PAIRNUM, MAX_LAMBDA, FESTEP, BUFFER_SIZE,
            READ_BUFFER_SIZE, N_MP) == (
        Constants.MACS_VERSION, Constants.MAX_PAIRNUM, Constants.MAX_LAMBDA,
        Constants.FESTEP, Constants.BUFFER_SIZE, Constants.READ_BUFFER_SIZE,
        Constants.N_MP)
    assert EFFECTIVEGS is Constants.EFFECTIVEGS


def test_effectivegs_exact():
    # deeptools effective genome sizes: GRCh38, GRCm38, WBcel235, dm6
    assert EFFECTIVEGS == {"hs": 2913022398, "mm": 2652783500,
                           "ce": 100286401, "dm": 142573017}
    assert list(EFFECTIVEGS) == ["hs", "mm", "ce", "dm"]
    assert all(type(v) is int for v in EFFECTIVEGS.values())


# ------------------------------------
# MACS_VERSION
# ------------------------------------

def test_macs_version_is_semantic_version():
    assert re.fullmatch(r"\d+\.\d+\.\d+", MACS_VERSION)


def test_macs_version_matches_cli_version(run_macs3):
    proc = run_macs3(["--version"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == "macs3 %s\n" % MACS_VERSION
    assert proc.stderr == ""


def test_macs_version_matches_changelog_head(test_dir):
    # the first release entry of ChangeLog names the current version
    with open(test_dir.parent / "ChangeLog") as fh:
        head = fh.read(2000)
    m = re.search(r"MACS (\d+\.\d+\.\d+)", head)
    assert m is not None
    assert m.group(1) == MACS_VERSION


@pytest.mark.parametrize("module", ["MACS3.IO.OutputWriter",
                                    "MACS3.Commands.callpeak_cmd",
                                    "MACS3.Commands.callvar_cmd"])
def test_macs_version_used_by_writers(module):
    import importlib
    mod = importlib.import_module(module)
    assert mod.MACS_VERSION is MACS_VERSION


# ------------------------------------
# EFFECTIVEGS: -g shortcuts
# ------------------------------------

@pytest.mark.parametrize("subcommand", ["callpeak", "filterdup", "predictd"])
def test_gsize_help_lists_effectivegs_values(macs3_argparser, subcommand):
    # the -g help text quotes each shortcut with its comma-grouped size
    helptext = option(subparser(macs3_argparser, subcommand), "gsize").help
    for key, value in EFFECTIVEGS.items():
        assert "'%s'" % key in helptext
        assert "(%s)" % format(value, ",") in helptext


@pytest.mark.parametrize("subcommand", ["callpeak", "filterdup", "predictd"])
def test_gsize_default_is_a_shortcut(macs3_argparser, subcommand):
    assert option(subparser(macs3_argparser, subcommand),
                  "gsize").default == "hs"
    assert "hs" in EFFECTIVEGS


def test_optvalidator_uses_effectivegs():
    from MACS3.Utilities import OptValidator
    assert OptValidator.efgsize is EFFECTIVEGS


# ------------------------------------
# buffer sizes and limits used elsewhere
# ------------------------------------

@pytest.mark.parametrize("subcommand", BUFFERED_SUBCOMMANDS)
def test_buffer_size_option_default_matches_constant(macs3_argparser,
                                                     subcommand):
    args = {"callpeak": ["-t", "t.bed"],
            "filterdup": ["-i", "t.bed"],
            "predictd": ["-i", "t.bed"],
            "pileup": ["-i", "t.bed", "-o", "o.bdg"],
            "randsample": ["-i", "t.bed", "-p", "10"],
            "refinepeak": ["-b", "p.bed", "-i", "t.bed", "-o", "o.bed"],
            "hmmratac": ["-i", "t.bam"]}[subcommand]
    options = macs3_argparser.parse_args([subcommand] + args)
    assert options.buffer_size == BUFFER_SIZE


@pytest.mark.parametrize("module", ["MACS3.IO.Parser", "MACS3.IO.BAM"])
def test_read_buffer_size_used_by_readers(module):
    import importlib
    mod = importlib.import_module(module)
    assert mod.READ_BUFFER_SIZE == READ_BUFFER_SIZE


@pytest.mark.parametrize("module", ["MACS3.Commands.callpeak_cmd",
                                    "MACS3.Commands.predictd_cmd"])
def test_max_pairnum_used_by_model_building(module):
    import importlib
    mod = importlib.import_module(module)
    assert mod.MAX_PAIRNUM == MAX_PAIRNUM
