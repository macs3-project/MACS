# Bdgbroadcall

## Overview
The `bdgbroadcall` command is part of the MACS3 suite of tools and is used to call broad peaks from a single bedGraph track for scores. It is particularly useful in ChIP-Seq analysis for identifying broad peaks in the data.

## Detailed Description

The `bdgbroadcall` command takes an input bedGraph file and produces an output file with broad peaks called. It is important to note: only bedGraph files from MACS3 are acceptable to use in the `bdgbroadcall` command, as All regions on the same chromosome in the bedGraph file should be continuous.

## Command Line Options

The command line options for `bdgbroadcall` are defined in `/MACS3/Commands/bdgbroadcall_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: MACS score in bedGraph. REQUIRED
- `-c` or `--cutoff-peak`: Cutoff for peaks depending on which method you used for the score track. If the file contains qvalue scores from MACS3, score 2 means qvalue 0.01. Regions with signals lower than the cutoff will not be considered as enriched regions. DEFAULT: 2
- `-C` or `--cutoff-link`: Cutoff for linking regions/low abundance regions depending on which method you used for the score track. If the file contains qvalue scores from MACS3, score 1 means qvalue 0.1, and score 0.3 means qvalue 0.5. DEFAULT: 1
- `-l` or `--min-length`: Minimum length of peak, better to set it as d value. DEFAULT: 200
- `-g` or `--lvl1-max-gap`: Maximum gap between significant peaks, better to set it as the tag size. DEFAULT: 30
- `-G` or `--lvl2-max-gap`: Maximum linking between significant peaks, better to set it as 4 times the d value. DEFAULT: 800
- `--no-trackline`: Tells MACS not to include a trackline with bedGraph files. The trackline is required by UCSC.
- `--verbose`: Set verbose level of runtime messages. 0: only show critical messages, 1: show additional warning messages, 2: show process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that directory. Default: the current working directory
- `-o` or `--ofile`: Output file name. Mutually exclusive with --o-prefix.
- `--o-prefix`: Output file prefix. Mutually exclusive with -o/--ofile.

## Example Usage

Here is an example of how to use the `bdgbroadcall` command:

```bash
macs3 bdgbroadcall -i input.bedGraph -o output.bedGraph --cutoff 1.0 --minlen 500 --maxgap 1000 --cutofflink 1.5 --maxgaplink 500
```

In this example, the program will call broad peaks from the `input.bedGraph` file and write the result to `output.bedGraph`. The cutoff value for calling peaks is set to 1.0, the minimum length of peaks is set to 500, the maximum gap between peaks is set to 1000, the cutoff value for linking peaks is set to 1.5, and the maximum gap for linking peaks is set to 500.

