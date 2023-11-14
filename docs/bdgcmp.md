# bdgcmp

## Overview
The `bdgcmp` command is part of the MACS3 suite of tools and is used
to compare two bedGraph files in each basepair that are commonly
covered by the two files. The typical use case is to calculate pvalue
or qvalue using Poisson model for each basepair given a treatment
pileup signal file in bedGraph format and a control lambda bedGraph
file. But we provides more functions rather than pvalue and qvalue,
including subtract, division (FE) and more.

## Detailed Description

The `bdgcmp` command takes two input bedGraph files (e.g. a control
and a treatment bedgraph) and produces an output bedGraph of
comparison scores for each genomic position involved in the bedGraph
files. The `bdgcmp` command normally is used to deduct noise from a
signal track in bedGraph (e.g. ChIP treatment) over another signal
track in bedGraph (e.g. control). Note: All regions on the same
chromosome in the bedGraph file should be continuous so we recommand
you use the bedGraph files from MACS3. We provide the following
function to 'compare two tracks':

- `ppois` Poisson p-value (-log10(pvalue) form) using the second file
  (-c) as lambda and treatment (-t) as observation
- `qpoi`s The q-value through a BH process for poisson pvalues
- `subtract` Subtraction from treatment
- `FE` linear scale fold enrichment, or the score from file A divided
  by the score from file B
- `logFE` log10 fold enrichment(need to set pseudocount)
- `logLR` log10 likelihood between ChIP-enriched model and open
  chromatin model (need to set pseudocount)
- `slogLR` symmetric log10 likelihood between two ChIP-enrichment
  models using Poison distribution, and this can be used to compare
  ChIP signals from two differen conditions (differential binding)
- `max` Maximum value between the two tracks.

## Command Line Options

Here is a brief description of the command line options for `bdgcmp` :

- `-t` or `--tfile`: Treatment bedGraph file, e.g. *_treat_pileup.bdg
  from MACS. REQUIRED
- `-c` or `--cfile`: Control bedGraph file, e.g. *_control_lambda.bdg
  from MACS. REQUIRED
- `-S` or `--scaling-factor`: Scaling factor for treatment and control
  track. Keep it as 1.0 or default in most cases. Set it ONLY while
  you have SPMR output from MACS3 callpeak, and plan to calculate
  scores as MACS3 callpeak module. If you want to simulate 'callpeak'
  w/o '--to-large', calculate effective smaller sample size after
  filtering redundant reads in million (e.g., put 31.415926 if
  effective reads are 31,415,926) and input it for '-S'; for 'callpeak
  --to-large', calculate effective reads in a larger sample. DEFAULT:
  1.0
- `-p` or `--pseudocount`: The pseudocount used for calculating logLR,
  logFE or FE. The count will be applied after normalization of
  sequencing depth. DEFAULT: 0.0, no pseudocount is applied.
- `-m` or `--method`: Method to use while calculating a score in any
  bin by comparing the treatment value and control value. Available
  choices are: `ppois`, `qpois`, `subtract`, `logFE`,` logLR`,
  `slogLR`, and `max`. They represent Poisson P-value (-log10(pvalue)
  form) using control as lambda and treatment as observation, q-value
  through a BH process for Poisson p-values, subtraction from
  treatment, linear scale fold enrichment, log10 fold enrichment (need
  to set pseudocount), log10 likelihood between ChIP-enriched model
  and open chromatin model (need to set pseudocount), symmetric log10
  likelihood between two ChIP-enrichment models, or the maximum value
  between the two tracks. The default option is ppois.
- `--verbose`: Set the verbose level of runtime messages. 0: only show
  critical messages, 1: show additional warning messages, 2: show
  process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that
  directory. Default: the current working directory
- `--o-prefix`: The PREFIX of the output bedGraph file to write
  scores. If it is given as A, and the method is 'ppois', the output
  file will be A_ppois.bdg. Mutually exclusive with -o/--ofile.
- `-o` or `--ofile`: Output filename. Mutually exclusive with
  --o-prefix. The number and the order of arguments for --ofile must
  be the same as for -m.

## Example Usage

Here is an example of how to use the `bdgcmp` command:

```bash
macs3 bdgcmp -t treatment.bedGraph -c control.bedGraph -m ppois -p 1.0 -S 1.0 -o output.bedGraph
```

In this example, the program will compare the `treatment.bedGraph`
file and the `control.bedGraph` file and write the result to
`output.bedGraph`. The method used for comparison is `ppois`, the
pseudo-count is set to 1.0, and the scaling factor is set to 1.0.
