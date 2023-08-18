# Bdgdiff

## Overview
The `bdgdiff` command is part of the MACS3 suite of tools and is used to call differential peaks from four bedGraph tracks for scores. It is particularly useful in ChIP-Seq analysis for identifying differential peaks in the data.

## Detailed Description

The `bdgdiff` command takes four input bedGraph files (two treatment and two control files) and produces three output files with differential peaks called. It uses an efficient algorithm to detect and call differential peaks, greatly improving the quality of your data for further analysis.

Differential peak detection based on paired four bedgraph files. Note: All regions on the same chromosome in the bedGraph file
                        should be continuous so only bedGraph files from MACS3 are accpetable.

## Command Line Options

The command line options for `bdgdiff` are defined in `/MACS3/Commands/bdgdiff_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `--t1`: MACS pileup bedGraph for condition 1. Incompatible with callpeak --SPMR output. REQUIRED
- `--t2`: MACS pileup bedGraph for condition 2. Incompatible with callpeak --SPMR output. REQUIRED
- `--c1`: MACS control lambda bedGraph for condition 1. Incompatible with callpeak --SPMR output. REQUIRED
- `--c2`: MACS control lambda bedGraph for condition 2. Incompatible with callpeak --SPMR output. REQUIRED
- `-C` or `--cutoff`: logLR cutoff. Regions with signals lower than the cutoff will not be considered as enriched regions. DEFAULT: 3 (likelihood ratio=1000)
- `-l` or `--min-len`: Minimum length of the differential region. Try a bigger value to remove small regions. DEFAULT: 200
- `-g` or `--max-gap`: Maximum gap to merge nearby differential regions. Consider a wider gap for broad marks. The maximum gap should be smaller than the minimum length (-g). DEFAULT: 100
- `--d1` or `--depth1`: Sequencing depth (# of non-redundant reads in million) for condition 1. It will be used together with --d2. See the description for --d2 below for how to assign them. Default: 1
- `--d2` or `--depth2`: Sequencing depth (# of non-redundant reads in million) for condition 2. It will be used together with --d1. DEPTH1 and DEPTH2 will be used to calculate the scaling factor for each sample, to down-scale the larger sample to the level of the smaller one. For example, while comparing 10 million condition 1 and 20 million condition 2, use --d1 10 --d2 20, then the pileup value in bedGraph for condition 2 will be divided by 2. Default: 1
- `--verbose`: Set the verbose level of runtime messages. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `--o-prefix`: Output file prefix. Actual files will be named as PREFIX_cond1.bed, PREFIX_cond2.bed, and PREFIX_common.bed. Mutually exclusive with -o/--ofile.
- `-o` or `--ofile`: Output filenames. Must give three arguments in order: 1. file for unique regions in condition 1; 2. file for unique regions in condition 2; 3. file for common regions in both conditions. Note: mutually exclusive with --o-prefix.


## Example Usage

Here is an example of how to use the `bdgdiff` command:

```bash
macs3 bdgdiff -t1 treatment1.bedGraph -c1 control1.bedGraph -t2 treatment2.bedGraph -c2 control2.bedGraph --depth1 1.0 --depth2 1.0 -o output.bedGraph --minlen 500 --maxgap 1000 --cutoff 1.0
```

In this example, the program will call differential peaks from the two pairs of treatment and control bedGraph files and write the result to `output.bedGraph`. The depth of the first and second condition is set to 1.0, the minimum length of differential peaks is set to 500, the maximum gap between differential peaks is set to 1000, and the cutoff for LLR to call differential peaks is set to 1.0.

