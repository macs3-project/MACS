"""Pileup alignment files using NumPy-backed PileupV2 routines.

This command mirrors :mod:`MACS3.Commands.pileup_cmd` but routes
through the PileupV2 implementation.
"""

import os

from MACS3.Utilities.OptValidator import opt_validate_pileup
from MACS3.Signal import PileupV2

# Reuse the existing loaders so option handling remains consistent.
from MACS3.Commands.pileup_cmd import load_tag_files_options, load_frag_files_options  # noqa: E402


def run(o_options):
    """Main entry for the pileup_v2 command."""
    options = opt_validate_pileup(o_options)
    info = options.info

    options.PE_MODE = options.format in ("BAMPE", "BEDPE", "FRAG")

    outfile = os.path.join(options.outdir, options.outputfile)
    if os.path.isfile(outfile):
        info("# Existing file %s will be replaced!" % outfile)
        os.unlink(outfile)

    info("# read alignment files...")
    if options.PE_MODE:
        info("# read input file in Paired-end mode.")
        treat = load_frag_files_options(options)
        t0 = treat.total
        info("# total fragments/pairs in alignment file: %d" % t0)

        if options.format == "FRAG" and options.barcodefile:
            info("# extract fragments with given barcodes")
            barcodes_set = set()
            with open(options.barcodefile, "r") as bfhd:
                for l in bfhd:
                    barcodes_set.add(l.rstrip().encode())
            treat = treat.subset(barcodes_set)
            info("#   extracted %d fragments", treat.total)

        info("# Pileup paired-end alignment file with PileupV2.")
        PileupV2.pileup_and_write_pe(treat,
                                     outfile.encode(),
                                     scale_factor=1,
                                     baseline_value=0.0)
    else:
        (tsize, treat) = load_tag_files_options(options)
        info("# tag size = %d", tsize)
        t0 = treat.total
        info("# total tags in alignment file: %d", t0)

        if options.bothdirection:
            info("# Pileup alignment file (bidirectional), extend each read +/- %d bps" % options.extsize)
            PileupV2.pileup_and_write_se(
                treat,
                outfile.encode(),
                options.extsize * 2,
                1,
                directional=False,
                halfextension=False,
            )
        else:
            info("# Pileup alignment file, extend each read downstream by %d bps" % options.extsize)
            PileupV2.pileup_and_write_se(
                treat,
                outfile.encode(),
                options.extsize,
                1,
                directional=True,
                halfextension=False,
            )

    info("# Done! Check %s" % options.outputfile)
