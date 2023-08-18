# Bdgpeakcall

## Overview
The `bdgpeakcall` command is part of the MACS3 suite of tools and is used to call peaks from a single bedGraph track for scores. It is particularly useful in ChIP-Seq analysis for identifying peaks in the data.

## Detailed Description

The `bdgpeakcall` command takes an input bedGraph file and produces an output file with peaks called. It uses an efficient algorithm to detect and call peaks, greatly improving the quality of your data for further analysis.



Call peaks from bedGraph output. Note: All regions on the same chromosome in the bedGraph file should be continuous so only
                        bedGraph files from MACS3 are accpetable.

## Command Line Options

The command line options for `bdgpeakcall` are defined in `/MACS3/Commands/bdgpeakcall_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: MACS score in bedGraph. REQUIRED
- `-c` or `--cutoff`: Cutoff depends on which method you used for the score track. If the file contains p-value scores from MACS3, score 5 means pvalue 1e-5. Regions with signals lower than the cutoff will not be considered as enriched regions. DEFAULT: 5
- `-l` or `--min-length`: Minimum length of peak, better to set it as d value. DEFAULT: 200
- `-g` or `--max-gap`: Maximum gap between significant points in a peak, better to set it as the tag size. DEFAULT: 30
- `--cutoff-analysis`: While set, bdgpeakcall will analyze the number or total length of peaks that can be called by different cutoff then output a summary table to help the user decide a better cutoff. Note, minlen and maxgap may affect the results. DEFAULT: False
- `--no-trackline`: Tells MACS not to include a trackline with bedGraph files. The trackline is required by UCSC.
- `--verbose`: Set the verbose level of runtime messages. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `-o` or `--ofile`: Output file name. Mutually exclusive with --o-prefix.
- `--o-prefix`: Output file prefix. Mutually exclusive with -o/--ofile.


## Example Usage

Here is an example of how to use the `bdgpeakcall` command:

```bash
macs3 bdgpeakcall -i input.bedGraph -o output.narrowPeak --cutoff 1.0 --minlen 500 --maxgap 1000 --cutoff-analysis
```

In this example, the program will call peaks from the `input.bedGraph` file and write the result to `output.narrowPeak`. The cutoff for calling peaks is set to 1.0, the minimum length of peaks is set to 500, the maximum gap between peaks is set to 1000, and the cutoff-analysis option is enabled.
