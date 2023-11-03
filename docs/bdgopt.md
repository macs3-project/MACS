# bdgopt

## Overview
The `bdgopt` command is part of the MACS3 suite of tools and is used
to modify a single bedGraph file. It provides various operations to
modify the value in the fourth column of the bedGraph file -- the
score column. 

## Detailed Description

The `bdgopt` command takes an input bedGraph file and produces an
output file with modified scores. It uses various methods to modify
the scores in the bedGraph files, greatly improving the flexibility of
your data for further analysis. Operations on score column of bedGraph
file include multiplication, addition, maximization with a given
value, minimization with a given value, and pvalue-to-qvalue
conversion (-log10 form).  Note: All regions on the same chromosome in
the bedGraph file should be continuous. We recommend to use the
bedGraph files from MACS3. 

## Command Line Options

Here is a brief overview of the commandline options:

- `-i` or `--ifile`: A bedGraph file containing scores. Note: this
  must be a bedGraph file covering the ENTIRE genome. REQUIRED
- `-m` or `--method`: Method to modify the score column of the
  bedGraph file. Available choices are: multiply, add, max, min, or
  p2q. 
    - `multiply`: The EXTRAPARAM is required and will be multiplied to
      the score column. If you intend to divide the score column by X,
      use the value of 1/X as EXTRAPARAM. 
    - `add`: The EXTRAPARAM is required and will be added to the score
      column. If you intend to subtract the score column by X, use the
      value of -X as EXTRAPARAM. 
    - `max`: The EXTRAPARAM is required and will take the maximum
      value between the score and the EXTRAPARAM. 
    - `min`: The EXTRAPARAM is required and will take the minimum
      value between the score and the EXTRAPARAM. 
    - `p2q`: This will convert p-value scores to q-value scores using
      the Benjamini-Hochberg process. The EXTRAPARAM is not
      required. This method assumes the scores are -log10 p-value from
      MACS3. Any other types of scores will cause unexpected errors. 
- `-p` or `--extra-param`: The extra parameter for METHOD. Check the
  detail of the -m option. 
- `--outdir`: If specified, all output files will be written to that
  directory. Default: the current working directory 
- `-o` or `--ofile`: Output BEDGraph filename.
- `--verbose`: Set the verbose level of runtime messages. 0: only show
  critical messages, 1: show additional warning messages, 2: show
  process information, 3: show debug messages. DEFAULT: 2 


## Example Usage

Here is an example of how to use the `bdgopt` command:

```bash
macs3 bdgopt -i input.bedGraph -o output.bedGraph --method multiply --extraparam 2.0
```

In this example, the program will modify the scores in the
`input.bedGraph` file and write the result to `output.bedGraph`. The
method used for modification is `multiply`, and the extra parameter is
set to 2.0, meaning that all scores will be multiplied by 2.0. 

Some use cases for `bdgopt`:

1. If you plan to scale up or down the scores in the bedGraph file,
   you can use `-m multiply` with a larger than 1 (>1) EXTRAPARAM in
   `-p` to scale up, or a positive value smaller than 1 (>0 and <1)
   EXTRAPARAM in `-p` to scale up; or `-m add` with a positive value
   in `-p` to increase the scores by a fixed amount or a negative
   value to decrease the scores.
2. If you want to cap the score in the bedGraph, you can use `-m max`
   with the upper limit score you want to use in `-p`. If you want to
   set the minimum score in the bedGraph, for example to set the whole
   genome background signal in the MACS contral lambda track, you can
   use `-m min` with the value in `-p`.

