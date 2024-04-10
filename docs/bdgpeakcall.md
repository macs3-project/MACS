# bdgpeakcall

## Overview
The `bdgpeakcall` command is part of the MACS3 suite of tools and is
used to call peaks from a single bedGraph track for scores.

## Detailed Description
The `bdgpeakcall` command takes an input bedGraph file, a cutoff
value, and the options to control peak width, then produces an output
file with peaks called. This tool can be used as a generic peak caller
that can be applied to any scoring system in a BedGraph file, no
matter the score is the pileup, the p-value, or the fold change -- or
any other measurement to decide the 'significancy' of the genomic
loci.

Note: All regions on the same chromosome in the bedGraph file should
be continuous.

## Command Line Options

Here is a brief overview of the command line options for `bdgpeakcall`:

- `-i` or `--ifile`: MACS score, or any numerical measurement score in
   bedGraph. The only assumption on the score is that higher the score
   is, more significant the genomic loci is. REQUIRED
- `-c` or `--cutoff`: Cutoff depends on which method you used for the
   score track. If the file contains -log10(p-value)  scores from
   MACS3, score 5 means pvalue 1e-5. Regions with signals lower than
   the cutoff will not be considered as enriched regions. DEFAULT: 5 
- `-l` or `--min-length`: Minimum length of peak, better to set it as
   d value if we will deal with MACS3 generated scores. DEFAULT: 200 
- `-g` or `--max-gap`: Maximum gap between significant points in a 
   peak, better to set it as the tag size if we will deal with MACS3
   generated scores. DEFAULT: 30 
- `--cutoff-analysis`: While set, bdgpeakcall will analyze the number
   or total length of peaks that can be called by different cutoff
   then output a summary table to help the user decide a better
   cutoff. Note, minlen and maxgap may affect the results. Also, if
   this option is on, `bdgpeakcall` will analyze the outcome of
   different cutoff values and then quit without calling any peaks.
   DEFAULT: False
- `--cutoff-analysis-steps`: Steps for performing cutoff analysis. It
   will be used to decide which cutoff value should be included in the
   final report. Larger the value, higher resolution the cutoff
   analysis can be. The cutoff analysis function will first find the
   smallest (at least 0) and the largest (at most 1,000) score in the
   bedGraph, then break the range of score into
   `CUTOFF_ANALYSIS_STEPS` intervals. It will then use each score as
   cutoff to call peaks and calculate the total number of candidate
   peaks, the total basepairs of peaks, and the average length of peak
   in basepair. Please note that the final report ideally should
   include `CUTOFF_ANALYSIS_STEPS` rows, but in practice, if the
   cutoff yield zero peak, the row for that value won't be included.
   DEFAULT: 100", default = 100 )
- `--no-trackline`: Tells MACS not to include a trackline with
   bedGraph files. The trackline is used by UCSC for the options for
   display.
- `--verbose`: Set the verbose level of runtime messages. 0: only show
   critical messages, 1: show additional warning messages, 2: show
   process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that
   directory. Default: the current working directory 
- `-o` or `--ofile`: Output file name. Mutually exclusive with
   `--o-prefix`. 
- `--o-prefix`: Output file prefix. Mutually exclusive with
   `-o/--ofile`. 


## Example Usage

Here is an example of how to use the `bdgpeakcall` command:

```bash
macs3 bdgpeakcall -i input.bedGraph -o output.narrowPeak --cutoff 1.0 --minlen 500 --maxgap 1000
```

In this example, the program will call peaks from the `input.bedGraph`
file and write the result to `output.narrowPeak`. The cutoff for
calling peaks is set to 1.0, the minimum length of peaks is set to
500, and the maximum gap between peaks is set to 1000.

## Cutoff Analysis

The cutoff analysis function is provided by `--cutoff-analysis` option
in `callpeak`, `bdgpeakcall`, and `hmmratac`. However, the function is
`bdgpeakcall` is more flexible and can be applied on any scoring
scheme. We will separate this function into a dedicated subcommand in
the future.

Please note that if this `--cutoff-anlaysis` option is on, the
`bdgpeakcall` won't write any results of the peaks into narrowPeak
format file, ignoring `-c` you specified. Instead, it will write a
cutoff analysis report (`-o`) and quit.

When the option is on, we will first calculate the window of step for
increasing the cutoff from the lowest (`min_v`) to the highest
(`max_v`), using the `--cutoff-analysis-steps`:

`win_step = (max_v-min_v)/steps`. 

Then for each cutoff we plan to investigate, we will check the number
of peaks that can be called, their average peak length, and their
total length. The smallest cutoff (close to `min_v`) and any cutoff
that can't lead to any peak will be excluded in the final report.

The report consists of four columns:

1. score: the possible fold change cutoff value.
2. npeaks: the number of peaks under this cutoff.
3. lpeaks: the total length of all peaks.
4. avelpeak: the average length of peaks.

While there's no universal rule to suggest the best cutoff, here are a
few suggestions:

- You can use elbow analysis to find the cutoff that dramatically
   change the trend of npeaks, lpeaks, or avelpeak. But you need to
   think about how to define 'dramatical change'.
- You can use some common expectation to decide the cutoff. For
   example, the number of peaks should be thousands/ or the avelpeak
   should be around 500bps. Of course, it's arbitrary but the table
   will give you some insight.
