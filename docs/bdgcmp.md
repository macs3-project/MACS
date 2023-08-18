# Bdgcmp

## Overview
The `bdgcmp` command is part of the MACS3 suite of tools and is used to compare bedGraph files. It is particularly useful in ChIP-Seq analysis for comparing different bedGraph files.

## Detailed Description

The `bdgcmp` command takes two input bedGraph files (a control and a treatment bedgraph) and produces an output bedGraph of comparison scores. The `bdgcmp` command is used to deduct noise by comparing two signal tracks in bedGraph. Note: All regions on the same chromosome in the bedGraph file should be continuous so only bedGraph files from MACS3 are acceptable.

## Command Line Options

The command line options for `bdgcmp` are defined in `/MACS3/Commands/bdgcmp_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-t` or `--tfile`: Treatment bedGraph file, e.g. *_treat_pileup.bdg from MACSv2. REQUIRED
- `-c` or `--cfile`: Control bedGraph file, e.g. *_control_lambda.bdg from MACSv2. REQUIRED
- `-S` or `--scaling-factor`: Scaling factor for treatment and control track. Keep it as 1.0 or default in most cases. Set it ONLY while you have SPMR output from MACS3 callpeak, and plan to calculate scores as MACS3 callpeak module. If you want to simulate 'callpeak' w/o '--to-large', calculate effective smaller sample size after filtering redundant reads in million (e.g., put 31.415926 if effective reads are 31,415,926) and input it for '-S'; for 'callpeak --to-large', calculate effective reads in a larger sample. DEFAULT: 1.0
- `-p` or `--pseudocount`: The pseudocount used for calculating logLR, logFE or FE. The count will be applied after normalization of sequencing depth. DEFAULT: 0.0, no pseudocount is applied.
- `-m` or `--method`: Method to use while calculating a score in any bin by comparing the treatment value and control value. Available choices are: ppois, qpois, subtract, logFE, logLR, and slogLR. They represent Poisson P-value (-log10(pvalue) form) using control as lambda and treatment as observation, q-value through a BH process for Poisson p-values, subtraction from treatment, linear scale fold enrichment, log10 fold enrichment (need to set pseudocount), log10 likelihood between ChIP-enriched model and open chromatin model (need to set pseudocount), symmetric log10 likelihood between two ChIP-enrichment models, or the maximum value between the two tracks. The default option is ppois.
- `--verbose`: Set the verbose level of runtime messages. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `--o-prefix`: The PREFIX of the output bedGraph file to write scores. If it is given as A, and the method is 'ppois', the output file will be A_ppois.bdg. Mutually exclusive with -o/--ofile.
- `-o` or `--ofile`: Output filename. Mutually exclusive with --o-prefix. The number and the order of arguments for --ofile must be the same as for -m.

## Example Usage

Here is an example of how to use the `bdgcmp` command:

```bash
macs3 bdgcmp -t treatment.bedGraph -c control.bedGraph -m ppois -p 1.0 -S 1.0 -o output.bedGraph
```

In this example, the program will compare the `treatment.bedGraph` file and the `control.bedGraph` file and write the result to `output.bedGraph`. The method used for comparison is `ppois`, the pseudo-count is set to 1.0, and the scaling factor is set to 1.0.
