#!/usr/bin/env python
"""Module Description: Test the top-level ``bin/macs3`` script: main(),
prepare_argparser() and the add_*_parser functions.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import argparse
import importlib
import importlib.machinery
import importlib.util
import os
import re
import sys
from pathlib import Path

import pytest

from MACS3.Utilities.Constants import MACS_VERSION

# ------------------------------------
# local helpers
# ------------------------------------

MACS3_SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "macs3"

SUBCOMMANDS = ["callpeak", "bdgpeakcall", "bdgbroadcall", "bdgcmp", "bdgopt",
               "cmbreps", "bdgdiff", "filterdup", "predictd", "pileup",
               "randsample", "refinepeak", "callvar", "hmmratac"]

# every option of every subcommand (from the add_*_parser functions)
OPTIONS = {
    "callpeak": "-t --treatment -c --control -f --format -g --gsize -s "
                "--tsize --keep-dup --barcodes --max-count --outdir -n --name "
                "-B --bdg --verbose --trackline --SPMR --nomodel --shift "
                "--extsize --bw --d-min -m --mfold --fix-bimodal -q --qvalue "
                "-p --pvalue --scale-to --down-sample --seed --tempdir "
                "--nolambda --slocal --llocal --max-gap --min-length --broad "
                "--broad-cutoff --cutoff-analysis --call-summits --fe-cutoff "
                "--to-large --ratio --buffer-size",
    "bdgpeakcall": "-i --ifile -c --cutoff -l --min-length -g --max-gap "
                   "--cutoff-analysis --cutoff-analysis-max "
                   "--cutoff-analysis-steps --no-trackline --verbose "
                   "--outdir -o --ofile --o-prefix",
    "bdgbroadcall": "-i --ifile -c --cutoff-peak -C --cutoff-link -l "
                    "--min-length -g --lvl1-max-gap -G --lvl2-max-gap "
                    "--no-trackline --verbose --outdir -o --ofile --o-prefix",
    "bdgcmp": "-t --tfile -c --cfile -S --scaling-factor -p --pseudocount "
              "-m --method --verbose --outdir --o-prefix -o --ofile",
    "bdgopt": "-i --ifile -m --method -p --extra-param --outdir -o --ofile "
              "--verbose",
    "cmbreps": "-i -m --method --outdir -o --ofile --verbose",
    "bdgdiff": "--t1 --t2 --c1 --c2 -C --cutoff -l --min-len -g --max-gap "
               "--d1 --depth1 --d2 --depth2 --verbose --outdir --o-prefix "
               "-o --ofile",
    "filterdup": "-i --ifile -f --format -g --gsize -s --tsize -p --pvalue "
                 "--keep-dup --buffer-size --verbose --outdir -o --ofile "
                 "-d --dry-run",
    "predictd": "-i --ifile -f --format -g --gsize -s --tsize --bw --d-min "
                "-m --mfold --outdir --rfile --buffer-size --verbose",
    "pileup": "-i --ifile -o --ofile --outdir -f --format --barcodes "
              "--max-count -B --both-direction --extsize --buffer-size "
              "--verbose",
    "randsample": "-i --ifile -p --percentage -n --number --seed -o --ofile "
                  "--outdir -s --tsize -f --format --buffer-size --verbose",
    "refinepeak": "-b -i --ifile -f --format -c --cutoff -w --window-size "
                  "--buffer-size --verbose --outdir -o --ofile --o-prefix",
    "callvar": "-b --peak -t --treatment -c --control --outdir -o --ofile "
               "--verbose -g --gq-hetero -G --gq-homo -Q -D -F --fermi "
               "--fermi-overlap --top2alleles-mratio --altallele-count "
               "--max-ar -m --multiple-processing",
    "hmmratac": "-i --input -f --format --barcodes --max-count --outdir -n "
                "--name --cutoff-analysis-only --cutoff-analysis-max "
                "--cutoff-analysis-steps --save-digested --save-states "
                "--save-likelihoods --save-training-data --no-fragem --means "
                "--stddevs --min-frag-p --binsize -u --upper -l --lower "
                "--maxTrain --training-flanking -t --training --model "
                "--modelonly --hmm-type -c --prescan-cutoff --minlen "
                "--pileup-short --randomSeed --decoding-steps -e --blacklist "
                "--remove-dup --verbose --buffer-size",
}

# required options, as argparse lists them when missing
REQUIRED = {
    "callpeak": "-t/--treatment",
    "bdgpeakcall": "-i/--ifile",
    "bdgbroadcall": "-i/--ifile",
    "bdgcmp": "-t/--tfile, -c/--cfile",
    "bdgopt": "-i/--ifile, -o/--ofile",
    "cmbreps": "-i, -o/--ofile",
    "bdgdiff": "--t1, --t2, --c1, --c2",
    "filterdup": "-i/--ifile",
    "predictd": "-i/--ifile",
    "pileup": "-i/--ifile, -o/--ofile",
    "randsample": "-i/--ifile",
    "refinepeak": "-b, -i/--ifile",
    "callvar": "-b/--peak, -t/--treatment, -o/--ofile",
    "hmmratac": "-i/--input",
}

# a minimal valid command line for each subcommand
MINIMAL = {
    "callpeak": ["-t", "t.bed"],
    "bdgpeakcall": ["-i", "x.bdg", "-o", "y"],
    "bdgbroadcall": ["-i", "x.bdg", "-o", "y"],
    "bdgcmp": ["-t", "t.bdg", "-c", "c.bdg", "-o", "y"],
    "bdgopt": ["-i", "x.bdg", "-o", "y"],
    "cmbreps": ["-i", "a.bdg", "b.bdg", "-o", "y"],
    "bdgdiff": ["--t1", "a", "--t2", "b", "--c1", "c", "--c2", "d",
                "--o-prefix", "p"],
    "filterdup": ["-i", "x.bed"],
    "predictd": ["-i", "x.bed"],
    "pileup": ["-i", "x.bed", "-o", "y"],
    "randsample": ["-i", "x.bed", "-p", "10"],
    "refinepeak": ["-b", "p.bed", "-i", "x.bed", "-o", "y"],
    "callvar": ["-b", "p.bed", "-t", "t.bam", "-o", "y"],
    "hmmratac": ["-i", "x.bam"],
}

OPTION_TOKEN = re.compile(r"(?<![\w-])(--?[A-Za-z][\w-]*)")


@pytest.fixture(scope="module")
def macs3_module():
    """``bin/macs3`` loaded as a module (it has no .py suffix)."""
    loader = importlib.machinery.SourceFileLoader("macs3_script_under_test",
                                                  str(MACS3_SCRIPT))
    spec = importlib.util.spec_from_loader("macs3_script_under_test", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


@pytest.fixture
def run_macs3(run_macs3):
    """The conftest runner with single-threaded OpenMP/BLAS pools and a
    fixed help width (argparse wraps at $COLUMNS)."""
    def _run(args, env=None, **kwargs):
        one = {v: "1" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                "MKL_NUM_THREADS")}
        one["COLUMNS"] = "100"
        one.update(env or {})
        return run_macs3(args, env=one, **kwargs)
    return _run


def write_tiny_bed(path):
    path.write_text("".join("chr1\t%d\t%d\tr\t0\t+\n" % (100 * i, 100 * i + 36)
                            for i in range(1, 6)))
    return path


# ------------------------------------
# --version, no subcommand, unknown subcommand, -h
# ------------------------------------

def test_cli_version(run_macs3):
    proc = run_macs3(["--version"], timeout=60)
    assert proc.returncode == 0
    assert proc.stdout == "macs3 %s\n" % MACS_VERSION
    assert proc.stderr == ""


def test_version_constant():
    """Upstream bumped the version from 3.0.4 to 3.0.5 in 6cefb42 and
    d8918ba (3.0.5 release, #744)."""
    assert MACS_VERSION == "3.0.5"


def test_cli_no_subcommand(run_macs3):
    proc = run_macs3([], timeout=60)
    assert proc.returncode == 2
    assert proc.stdout == ""
    assert proc.stderr.startswith("usage: macs3 [-h] [--version]")
    assert proc.stderr.rstrip("\n").endswith(
        "macs3: error: the following arguments are required: subcommand")


def test_cli_unknown_subcommand(run_macs3):
    proc = run_macs3(["peakcall"], timeout=60)
    assert proc.returncode == 2
    assert "macs3: error: argument subcommand: invalid choice: 'peakcall'" \
        in proc.stderr


def test_cli_top_level_help(run_macs3):
    proc = run_macs3(["-h"], timeout=60)
    assert proc.returncode == 0
    assert proc.stdout.startswith("usage: macs3 [-h] [--version]")
    assert "macs3 -- Model-based Analysis for ChIP-Sequencing" in proc.stdout
    assert "For command line options of each command, type: macs3 COMMAND -h" \
        in proc.stdout
    assert "{%s}" % ",".join(SUBCOMMANDS) in proc.stdout
    for sub in SUBCOMMANDS:
        assert re.search(r"^\s+%s\s" % sub, proc.stdout, re.M), sub


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_cli_subcommand_help_lists_every_option(run_macs3, sub):
    proc = run_macs3([sub, "-h"], timeout=60)
    assert proc.returncode == 0
    assert proc.stdout.startswith("usage: macs3 %s [-h]" % sub)
    found = set(OPTION_TOKEN.findall(proc.stdout))
    expected = set(OPTIONS[sub].split()) | {"-h", "--help"}
    if sub == "bdgpeakcall":
        expected.discard("--call-summits")
    assert expected <= found, sorted(expected - found)


def test_cli_bdgpeakcall_hides_call_summits(run_macs3):
    proc = run_macs3(["bdgpeakcall", "-h"], timeout=60)
    assert proc.returncode == 0
    assert "--call-summits" not in proc.stdout


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_subcommand_missing_required_arguments(macs3_argparser, capsys, sub):
    with pytest.raises(SystemExit) as exc:
        macs3_argparser.parse_args([sub])
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert ("error: the following arguments are required: %s" %
            REQUIRED[sub]) in err


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_subcommand_minimal_arguments_parse(macs3_argparser, sub):
    options = macs3_argparser.parse_args([sub] + MINIMAL[sub])
    assert options.subcommand == sub
    assert options.outdir == ""


# ------------------------------------
# prepare_argparser and the add_*_parser helpers
# ------------------------------------

def test_prepare_argparser_subcommands(macs3_module):
    parser = macs3_module.prepare_argparser()
    assert isinstance(parser, argparse.ArgumentParser)
    sub_actions = [a for a in parser._actions
                   if isinstance(a, argparse._SubParsersAction)]
    assert len(sub_actions) == 1
    assert list(sub_actions[0].choices) == SUBCOMMANDS
    assert sub_actions[0].required is True
    assert sub_actions[0].dest == "subcommand"


def test_prepare_argparser_returns_new_parser(macs3_module):
    assert macs3_module.prepare_argparser() is not \
        macs3_module.prepare_argparser()


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_add_parser_functions(macs3_module, sub):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    assert getattr(macs3_module, "add_%s_parser" % sub)(subparsers) is None
    options = parser.parse_args([sub] + MINIMAL[sub])
    assert options.subcommand == sub


def test_add_outdir_option(macs3_module):
    parser = argparse.ArgumentParser()
    macs3_module.add_outdir_option(parser)
    assert parser.parse_args([]).outdir == ""
    assert parser.parse_args(["--outdir", "d"]).outdir == "d"


def test_add_output_group_required(macs3_module, capsys):
    parser = argparse.ArgumentParser()
    macs3_module.add_output_group(parser)
    assert parser.parse_args(["-o", "f"]).ofile == "f"
    assert parser.parse_args(["--o-prefix", "p"]).oprefix == "p"
    with pytest.raises(SystemExit):
        parser.parse_args([])
    assert "one of the arguments -o/--ofile --o-prefix is required" in \
        capsys.readouterr().err
    with pytest.raises(SystemExit):
        parser.parse_args(["-o", "f", "--o-prefix", "p"])
    assert "not allowed with argument" in capsys.readouterr().err


def test_add_output_group_optional(macs3_module):
    parser = argparse.ArgumentParser()
    macs3_module.add_output_group(parser, required=False)
    options = parser.parse_args([])
    assert options.ofile is None and options.oprefix is None


# ------------------------------------
# main(): dispatch to the subcommand's run()
# ------------------------------------

# subcommands whose run() lives in a module not named after them
DISPATCH_MODULE = {"pileup": "pileup_v2"}


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_main_dispatches_to_run(macs3_module, monkeypatch, sub):
    """``pileup`` dispatches to ``pileup_v2_cmd.run`` since upstream 88af6b3
    ("Use v2 implementations for pileup from now"); before it went to
    ``pileup_cmd.run``."""
    module = importlib.import_module("MACS3.Commands.%s_cmd"
                                     % DISPATCH_MODULE.get(sub, sub))
    calls = []
    monkeypatch.setattr(module, "run", calls.append)
    monkeypatch.setattr(sys, "argv", ["macs3", sub] + MINIMAL[sub])
    assert macs3_module.main() is None
    assert len(calls) == 1
    assert calls[0].subcommand == sub


def test_main_creates_outdir_before_run(macs3_module, monkeypatch, tmp_path):
    module = importlib.import_module("MACS3.Commands.filterdup_cmd")
    seen = []
    monkeypatch.setattr(module, "run",
                        lambda o: seen.append(os.path.isdir(o.outdir)))
    outdir = tmp_path / "x" / "y"
    monkeypatch.setattr(sys, "argv", ["macs3", "filterdup", "-i", "a.bed",
                                      "--outdir", str(outdir)])
    macs3_module.main()
    assert seen == [True]


def test_main_outdir_error_raises_system_exit(macs3_module, monkeypatch,
                                              tmp_path):
    blocker = tmp_path / "file.txt"
    blocker.write_text("x")
    target = str(blocker / "sub")
    monkeypatch.setattr(sys, "argv", ["macs3", "filterdup", "-i", "a.bed",
                                      "--outdir", target])
    with pytest.raises(SystemExit) as exc:
        macs3_module.main()
    assert exc.value.code == ("Output directory (%s) could not be created. "
                              "Terminate program." % target)


# ------------------------------------
# main(): --outdir through the CLI
# ------------------------------------

def test_cli_outdir_nested_is_created(run_macs3, tmp_path):
    bed = write_tiny_bed(tmp_path / "r.bed")
    outdir = tmp_path / "a" / "b" / "c"
    proc = run_macs3(["filterdup", "-i", bed, "--keep-dup", "all",
                      "--outdir", outdir, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert os.listdir(outdir) == ["out.bed"]


def test_cli_outdir_existing_is_used(run_macs3, tmp_path):
    bed = write_tiny_bed(tmp_path / "r.bed")
    outdir = tmp_path / "exists"
    outdir.mkdir()
    (outdir / "keep.txt").write_text("keep")
    proc = run_macs3(["filterdup", "-i", bed, "--keep-dup", "all",
                      "--outdir", outdir, "-o", "out.bed"], timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert sorted(os.listdir(outdir)) == ["keep.txt", "out.bed"]
    assert (outdir / "keep.txt").read_text() == "keep"


def test_cli_outdir_dangling_symlink(run_macs3, tmp_path):
    """os.path.exists() is False for a dangling link, and makedirs then
    raises FileExistsError."""
    link = tmp_path / "link"
    link.symlink_to(tmp_path / "missing")
    bed = write_tiny_bed(tmp_path / "r.bed")
    proc = run_macs3(["filterdup", "-i", bed, "--outdir", link], timeout=60)
    assert proc.returncode == 1
    assert proc.stderr == ("Output directory (%s) could not be created since "
                           "it already exists. Terminate program.\n" % link)


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root ignores directory permissions")
def test_cli_outdir_permission_denied(run_macs3, tmp_path):
    parent = tmp_path / "ro"
    parent.mkdir()
    bed = write_tiny_bed(tmp_path / "r.bed")
    parent.chmod(0o555)
    try:
        proc = run_macs3(["filterdup", "-i", bed,
                          "--outdir", parent / "child"], timeout=60)
    finally:
        parent.chmod(0o755)
    assert proc.returncode == 1
    assert proc.stderr == ("Output directory (%s) could not be created due "
                           "to permission. Terminate program.\n" %
                           (parent / "child"))


def test_cli_outdir_under_a_file(run_macs3, tmp_path):
    blocker = tmp_path / "file.txt"
    blocker.write_text("x")
    bed = write_tiny_bed(tmp_path / "r.bed")
    proc = run_macs3(["filterdup", "-i", bed, "--outdir", blocker / "sub"],
                     timeout=60)
    assert proc.returncode == 1
    assert proc.stderr == ("Output directory (%s) could not be created. "
                           "Terminate program.\n" % (blocker / "sub"))


def test_cli_outdir_checked_before_options(run_macs3, tmp_path):
    """The directory is handled in main() before any subcommand
    validation: an invalid --gsize is not reported."""
    blocker = tmp_path / "file.txt"
    blocker.write_text("x")
    proc = run_macs3(["filterdup", "-i", "missing.bed", "-g", "bad",
                      "--outdir", blocker / "sub"], timeout=60)
    assert proc.returncode == 1
    assert "gsize" not in proc.stderr
