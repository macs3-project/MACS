# refinepeak

## Overview
The `refinepeak` command is part of the MACS3 suite of tools and is
used to refine peak summits. It is particularly useful in ChIP-Seq
analysis where refining the peak summits can lead to more accurate
results. 

## Detailed Description

The `refinepeak` command takes in a BED file containing peaks and raw
reads alignment, then produces an output BED file with refined peak
summits.  It will refine peak summits and give scores measuring the
balance of Watson/Crick tags, inspired by SPP. Basically, we assume
that a good summit in a peak region should have balanced Watson/Crick
tags around.

## Command Line Options

Here is a brief overview of the `refinepeak` options:

- `-b`: Candidate peak file in BED format. REQUIRED.
- `-i` or `--ifile`: ChIP-seq alignment file. If multiple files are
  given as '-t A B C', then they will all be read and combined. Note
  that pair-end data is not supposed to work with this
  command. REQUIRED. 
- `-f` or `--format`: Format of the tag file.
  - `AUTO`: MACS3 will pick a format from "AUTO", "BED", "ELAND",
    "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE". Please check
    the definition in the README file if you choose
    ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE. DEFAULT: "AUTO" 
- `-c` or `--cutoff`: Cutoff. Regions with SPP wtd score lower than
  cutoff will not be considered. DEFAULT: 5 
- `-w` or `--window-size`: Scan window size on both sides of the
  summit (default: 100bp) 
- `--buffer-size`: Buffer size for incrementally increasing the
  internal array size to store read alignment information. In most
  cases, you don't have to change this parameter. However, if there
  are a large number of chromosomes/contigs/scaffolds in your
  alignment, it's recommended to specify a smaller buffer size in
  order to decrease memory usage (but it will take longer time to read
  alignment files). Minimum memory requested for reading an alignment
  file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT:
  100000 
- `--verbose`: Set the verbose level. 0: only show critical messages,
  1: show additional warning messages, 2: show process information, 3:
  show debug messages. If you want to know where the duplicate reads
  are, use 3. DEFAULT: 2 
- `--outdir`: If specified, all output files will be written to that
  directory. Default: the current working directory 
- `-o` or `--ofile`: Output file name. Mutually exclusive with
  --o-prefix. 
- `--o-prefix`: Output file prefix. Mutually exclusive with
  -o/--ofile. 

## Example Usage

Here is an example of how to use the `refinepeak` command:

```bash
macs3 refinepeak -b peaks.bed -i alignment.bam -o refined_peaks.bed
```

In this example, the program will refine the peak summits in the
`peaks.bed` file taking in the alignment file `alignment.bam`, and
write the result to `refined_peaks.bed`. 
