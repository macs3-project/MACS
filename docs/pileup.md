# Pileup

## Overview
The `pileup` command is part of the MACS3 suite of tools and is used to pile up alignment files. It is particularly useful in ChIP-Seq analysis where summarizing the read depth at each genomic location is required.

## Detailed Description

The `pileup` command takes in one or multiple input files and produces an output file with the piled-up alignments. It uses an efficient algorithm to pile up the alignments, improving the quality of your data for further analysis.

Pileup aligned reads with a given extension size (fragment size or d in MACS language). Note there will be no step for duplicate reads filtering or sequencing depth scaling, so you may need to do certain pre/post-processing.

## Command Line Options

The command line options for `pileup` are defined in `/MACS3/Commands/pileup_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: The input alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. Note that pair-end data is not supposed to work with this command. REQUIRED.
- `-o` or `--ofile`: The output bedGraph file name. If not specified, will write to standard output. REQUIRED
- `--outdir`: If specified all output files will be written to that directory. Default: the current working directory
- `-f` or `--format`: Format of tag file, "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE", "BAMPE", or "BEDPE". The default AUTO option will let 'macs3 pileup' decide which format the file is. DEFAULT: "AUTO", MACS3 will pick a format from "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM" and "BOWTIE". If the format is BAMPE or BEDPE, please specify it explicitly. Please note that when the format is BAMPE or BEDPE, the -B and --extsize options would be ignored.
- `-B` or `--both-direction`: By default, any read will be extended towards downstream direction by extension size. So it's [0,size-1] (1-based index system) for plus strand read and [-size+1,0] for minus strand read where position 0 is 5' end of the aligned read. Default behavior can simulate MACS3 way of piling up ChIP sample reads where extension size is set as fragment size/d. If this option is set as on, aligned reads will be extended in both upstream and downstream directions by extension size. It means [-size,size] where 0 is the 5' end of a aligned read. It can partially simulate MACS3 way of piling up control reads. However MACS3 local bias is calculated by maximizing the expected pileup over a ChIP fragment size/d estimated from 10kb, 1kb, d and whole genome background. This option will be ignored when the format is set as BAMPE or BEDPE. DEFAULT: False
- `--extsize`: The extension size in bps. Each alignment read will become a EXTSIZE of fragment, then be piled up. Check description for -B for detail. It's twice the `shiftsize` in old MACSv1 language. This option will be ignored when the format is set as BAMPE or BEDPE. DEFAULT: 200
- `--buffer-size`: Buffer size for incrementally increasing internal array size to store reads alignment information. In most cases, you don't have to change this parameter. However, if there are large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
- `--verbose`: Set verbose level. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show debug messages. If you want to know where are the duplicate reads, use 3. DEFAULT:2

## Example Usage

Here is an example of how to use the `pileup` command:

```bash
macs3 pileup -i treatment.bam -o piledup.bedGraph -f BAM -g hs -n experiment1
```

In this example, the program will pile up the alignments in the `treatment.bam` file and write the result to `piledup.bedGraph`. The input file is in BAM format, the genome size is set to 'hs' (human), and the name of the experiment is 'experiment1'.
