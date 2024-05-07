# cmbreps

## Overview
The `cmbreps` command is part of the MACS3 suite of tools and is used
to combine bedGraph files from replicates. It is particularly useful
in ChIP-Seq analysis where multiple replicates of the same experiment
are performed. 

## Detailed Description

The `cmbreps` command takes a list of input bedGraph files
(replicates) and produces an output file with combined scores. Note:
All regions on the same chromosome in the bedGraph file should be
continuous so bedGraph files from MACS3 are recommended. 

The `cmbreps` command provides different way to combine replicates,
compared with the `callpeak` command where all replicates will be
simply pooled. A possible usage is that: for each replicate, we can
first follow the instructions in the [Advanced Step-by-step Peak
Calling](./Advanced_Step-by-step_Peak_Calling.md) to generate the
p-value scores through `bdgcmp -m ppois`, use `cmbreps -m fisher` to
use Fisher's combined probability test to combine all the p-value
score tracks and generate a single BedGraph, then call peaks using
`bdgpeakcall`.

## Command Line Options

Here is a brief overview of command line options:

- `-i IFILE1 IFILE2 [IFILE3 ...]`: MACS score in bedGraph for each
  replicate. Require at least 2 files. REQUIRED 
- `-m` or `--method`: Method to use while combining scores from
  replicates. 
  - `fisher`: Fisher's combined probability test. It requires scores
    in ppois form (-log10 pvalues) from `bdgcmp`. Other types of
    scores for this method may cause cmbreps unexpected errors. 
  - `max`: Take the maximum value from replicates for each genomic
    position. 
  - `mean`: Take the average value. Note, except for Fisher's method,
    max or mean will take scores AS IS which means they won't convert
    scores from log scale to linear scale or vice versa. 
- `--outdir`: If specified, all output files will be written to that
  directory. Default: the current working directory 
- `-o` or `--ofile`: Output BEDGraph filename for combined scores. 
- `--verbose`: Set the verbose level of runtime messages. 0: only show
  critical messages, 1: show additional warning messages, 2: show
  process information, 3: show debug messages. DEFAULT: 2 


## Example Usage

Here is an example of how to use the `cmbreps` command:

```bash
macs3 cmbreps -i replicate1.bedGraph replicate2.bedGraph replicate3.bedGraph -o combined.bedGraph --method mean
```

In this example, the program will combine the scores in the
`replicate1.bedGraph`, `replicate2.bedGraph`, and
`replicate3.bedGraph` files and write the result to
`combined.bedGraph`. The method used for combining scores is `mean` so
it will take the average score from the three replicates at each
genomic location.

