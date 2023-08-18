# Bdgpeakcall

## Overview
The `predictd` command is part of the MACS3 suite of tools and is used to call peaks from a single bedGraph track for scores. It is particularly useful in ChIP-Seq analysis for identifying peaks in the data.

## Detailed Description

The `predictd` command takes an input bedGraph file and ...

Predict d or fragment size from alignment results. In case of PE data, report the average insertion/fragment size from all pairs. *Will NOT filter duplicates*

## Command Line Options

The command line options for `predictd` are defined in `/MACS3/Commands/predictd_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i IFILE [IFILE ...]` or `--ifile IFILE [IFILE ...]`: ChIP-seq alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. Note that pair-end data is not supposed to work with this command. REQUIRED.
- `-f` or `--format`: Format of the tag file.
  - `AUTO`: MACS3 will pick a format from "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE", "BAMPE", and "BEDPE". However, if you want to decide the average insertion size/fragment size from PE data such as BEDPE or BAMPE, please specify the format as BAMPE or BEDPE since MACS3 won't automatically recognize these two formats with -f AUTO. Please be aware that in PE mode, -g, -s, --bw, --d-min, -m, and --rfile have NO effect. DEFAULT: "AUTO"
- `-g` or `--gsizeE`: Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7), and 'dm' for fruit fly (1.2e8). Default: hs
- `-s` or `--tsize`: Tag size. This will override the auto-detected tag size. DEFAULT: Not set
- `--bw`: Bandwidth for picking regions to compute the fragment size. This value is only used while building the shifting model. DEFAULT: 300
- `--d-min`: Minimum fragment size in base pairs. Any predicted fragment size less than this will be excluded. DEFAULT: 20
- `-m` or `--mfoldD`: Select the regions within MFOLD range of high-confidence enrichment ratio against background to build the model. Fold-enrichment in regions must be lower than the upper limit and higher than the lower limit. Use as "-m 10 30". DEFAULT: 5 50
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `--rfile`: PREFIX of the filename of the R script for drawing the X-correlation figure. DEFAULT: 'predictd' and the R file will be predicted_model.R
- `--buffer-size`: Buffer size for incrementally increasing the internal array size to store read alignment information. In most cases, you don't have to change this parameter. However, if there is a large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
- `--verbose`: Set the verbose level of runtime messages. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. DEFAULT: 2



## Example Usage

Here is an example of how to use the `predictd` command:

```bash
macs3 predictd -i input.bedGraph -o output.narrowPeak 
```

