# Cmbreps

## Overview
The `cmbreps` command is part of the MACS3 suite of tools and is used to combine replicate bedGraph files. It is particularly useful in ChIP-Seq analysis where multiple replicates of the same experiment are performed.

## Detailed Description

The `cmbreps` command takes a list of input bedGraph files (replicates) and produces an output file with combined scores. Note: All regions on the same chromosome in the bedGraph file should be continuous so only bedGraph files from MACS3 are acceptable.
## Command Line Options

The command line options for `cmbreps` are defined in `/MACS3/Commands/cmbreps_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: The input files. These should be in bedGraph format. This option requires at least 2 files such as '-i A B C D'. REQUIRED
- `-m` or `--method`: The method for combining scores from replicates. Options include 'fisher', 'max', or 'mean'. 1) fisher: Fisher's combined probability test. It requires scores in ppois form (-log10 pvalues) from bdgcmp. Other types of scores for this method may cause cmbreps unexpected errors. 2) max: take the maximum value from replicates for each genomic position. 3) mean: take the average value. Note, except for Fisher's method, max or mean will take scores AS IS which means they won't convert scores from log scale to linear scale or vice versa.
- `--outdir`: If specified all output files will be written to that directory. Default: the current working directory
- `-o` or `--ofile`: The output filename for combined scores in BEDGraph format.
- `--verbose`: Set verbose level of runtime message. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show debug messages. DEFAULT:2

## Example Usage

Here is an example of how to use the `cmbreps` command:

```bash
macs3 cmbreps -i replicate1.bedGraph replicate2.bedGraph replicate3.bedGraph -o combined.bedGraph --method mean
```

In this example, the program will combine the scores in the `replicate1.bedGraph`, `replicate2.bedGraph`, and `replicate3.bedGraph` files and write the result to `combined.bedGraph`. The method used for combining scores is `mean`.

