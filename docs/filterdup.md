# Filterdup

## Overview
The `filterdup` command is part of the MACS3 suite of tools and is used to filter duplicate reads from your data. It is particularly useful in sequencing analysis where duplicate reads can bias the downstream analysis.

## Detailed Description

The `filterdup` command takes an input file and produces an output file with duplicate reads removed. It uses an efficient algorithm to detect and filter duplicate reads, greatly improving the quality of your data for further analysis.

The `filteredup` command removes duplicate reads at the same position, then saves the remaining alignments to a BED or BEDPE file. If you use '--keep-dup all option', this script can be utilized to convert any acceptable format into BED or BEDPE format.

## Command Line Options

The command line options for `filterdup` are defined in `/MACS3/Commands/filterdup_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: The input alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. Note that pair-end data is not supposed to work with this command. REQUIRED.
- `-f`or `--format`: The format of the tag file. Options include: "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or "BAMPE" or "BEDPE". The default AUTO option will let 'macs3 filterdup' decide which format the file is. Please check the definition in README file if you choose ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE or BAMPE/BEDPE. DEFAULT: "AUTO"
- `-g` or `--gsize`: The effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8). DEFAULT:hs
- `-s` or `--tsize`: The tag size. This will override the auto detected tag size. DEFAULT: Not set
- `-p` or `--pvalue`: The pvalue cutoff for binomial distribution test. DEFAULT:1e-5
- `--keep-dup`: The number of duplicates to keep. It controls the 'macs3 filterdup' behavior towards duplicate tags/pairs at the exact same location -- the same coordination and the same strand. The 'auto' option makes 'macs3 filterdup' calculate the maximum tags at the exact same location based on binomal distribution using given -p as pvalue cutoff; and the 'all' option keeps every tags (useful if you only want to convert formats). If an integer is given, at most this number of tags will be kept at the same location. Note, MACS3 callpeak function uses KEEPDUPLICATES=1 as default. Note, if you've used samtools or picard to flag reads as 'PCR/Optical duplicate' in bit 1024, MACS3 will still read them although the reads may be decided by MACS3 as duplicate later. Default: auto
- `--buffer-size`: The buffer size for incrementally increasing internal array size to store reads alignment information. In most cases, you don't have to change this parameter. However, if there are large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
- `--verbose`: The verbose level. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show debug messages. If you want to know where are the duplicate reads, use 3. DEFAULT:2
- `--outdir`: If specified all output files will be written to that directory. Default: the current working directory
- `-o` or `--ofile`: The output BED file name. If not specified, will write to standard output. Note, if the input format is BAMPE or BEDPE, the output will be in BEDPE format. DEFAULT: stdout
- `-d` or `--dry-run`: When set, filterdup will only output numbers instead of writing output files, including maximum allowable duplicates, total number of reads before filtering, total number of reads after filtering, and redundant rate. Default: not set

## Example Usage

Here is an example of how to use the `filterdup` command:

```bash
macs3 filterdup -i input.bam -o output.bam --gsize hs --format AUTO --keep-dup 1 --buffer-size 100000
```

In this example, the program will remove duplicate reads from the `input.bam` file and write the result to `output.bam`. The mappable genome size is set to `hs` (Homo Sapiens), the format of the input file is determined automatically, and the program keeps only one duplicate.

