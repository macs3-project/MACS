# Refinepeak

## Overview
The `refinepeak` command is part of the MACS3 suite of tools and is used to refine peak summits. It is particularly useful in ChIP-Seq analysis where refining the peak summits can lead to more accurate results.

## Detailed Description

The `refinepeak` command takes in a BED file containing peaks and produces an output BED file with refined peak summits. It uses an efficient algorithm to refine the peak summits, improving the quality of your data for further analysis.


(Experimental) Take raw reads alignment, refine peak summits and give scores measuring balance of waston/crick tags. Inspired by
                        SPP.
## Command Line Options

The command line options for `refinepeak` are defined in `/MACS3/Commands/refinepeak_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-b`: Candidate peak file in BED format. REQUIRED.
- `-i` or `--ifile`: ChIP-seq alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. Note that pair-end data is not supposed to work with this command. REQUIRED.
- `-f` or `--format`: Format of the tag file.
  - `AUTO`: MACS3 will pick a format from "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE". Please check the definition in the README file if you choose ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE. DEFAULT: "AUTO"
- `-c` or `--cutoff`: Cutoff. Regions with SPP wtd score lower than cutoff will not be considered. DEFAULT: 5
- `-w` or `--window-size`: Scan window size on both sides of the summit (default: 100bp)
- `--buffer-size`: Buffer size for incrementally increasing the internal array size to store read alignment information. In most cases, you don't have to change this parameter. However, if there are a large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
- `--verbose`: Set the verbose level. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. If you want to know where the duplicate reads are, use 3. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `-o` or `--ofile`: Output file name. Mutually exclusive with --o-prefix.
- `--o-prefix`: Output file prefix. Mutually exclusive with -o/--ofile.


## Example Usage

Here is an example of how to use the `refinepeak` command:

```bash
macs3 refinepeak -i peaks.bed -o refined_peaks.bed -g hs -n experiment1 --bw 300 --nomodel
```

In this example, the program will refine the peak summits in the `peaks.bed` file and write the result to `refined_peaks.bed`. The genome size is set to 'hs' (human), the name of the experiment is 'experiment1', the bandwidth is set to 300, and the shifting model is not built.
