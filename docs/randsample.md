# Randsample

## Overview
The `randsample` command is part of the MACS3 suite of tools and is used to randomly sample a certain number or percentage of tags from alignment files. This can be useful in ChIP-Seq analysis where a subset of the data is required for downstream analysis.

## Detailed Description

The `randsample` command takes in one or multiple input alignment files and produces an output file with the randomly sampled tags. It uses an efficient algorithm to sample the tags, improving the quality of your data for further analysis.

## Command Line Options

The command line options for `randsample` are defined in `/MACS3/Commands/randsample_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: Alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. Note that pair-end data is not supposed to work with this command. REQUIRED.
- `-p` or `--percentage`: Percentage of tags you want to keep. Input 80.0 for 80%. This option can't be used at the same time with -n/--num. REQUIRED
- `-n` or `--number`: Number of tags you want to keep. Input 8000000 or 8e+6 for 8 million. This option can't be used at the same time with -p/--percent. Note that the number of tags in the output is approximate as the number specified here. REQUIRED
- `--seed`: Set the random seed while downsampling data. Must be a non-negative integer in order to be effective. DEFAULT: not set
- `-o` or `--ofile`: Output BED file name. If not specified, will write to standard output. Note, if the input format is BAMPE or BEDPE, the output will be in BEDPE format. DEFAULT: stdout
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `-s` or `--tsize`: Tag size. This will override the auto-detected tag size. DEFAULT: Not set
- `-f` or `--format`: Format of the tag file.
  - `AUTO`: MACS3 will pick a format from "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE", "BAMPE", and "BEDPE". Please check the definition in the README file if you choose ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE or BAMPE/BEDPE. DEFAULT: "AUTO"
- `--buffer-size`: Buffer size for incrementally increasing the internal array size to store read alignment information. In most cases, you don't have to change this parameter. However, if there are a large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
- `--verbose`: Set the verbose level. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. If you want to know where the duplicate reads are, use 3. DEFAULT: 2


## Example Usage

Here is an example of how to use the `randsample` command:

```bash
macs3 randsample -i treatment.bam -o sampled.bed -f BAM -g hs -n experiment1 --mfold 10 30
```

In this example, the program will randomly sample tags from the `treatment.bam` file and write the result to `sampled.bed`. The input file is in BAM format, the genome size is set to 'hs' (human), and the name of the experiment is 'experiment1'.
