# bdgdiff

## Overview
The `bdgdiff` command is part of the MACS3 suite of tools and is used
to call differential peaks from four bedGraph tracks of scores,
including treatment and control track for each condition.


## Detailed Description

The `bdgdiff` command takes four input bedGraph files (two treatment
and two control files) and produces three output files with
differential peaks called. Users should provide paired four bedgraph
files: for each condition, a treatment pileup signal track in bedGraph
format, and a control lambda track in bedGraph format. This
differential calling can only handle one replicate per condition, so
if you have multiple replicates, you should either combine the
replicates when `callpeak`, or choose other tool that can consider
within-group variation (such as DESeq2 or edgeR). The method we use to
define the differential peaks is based on multiple likelihood tests,
based on the poisson distribution. Suppose that we have two conditions
A and B, the unique binding sites in condition A over condition B
should be *more likely* to be a binding event in treatment A over
treatment B, and also *more likely* to be a real binding site in
condition A while comparing treatment A over control A; the unique
binding sites in condition B is defined in the same way; the common
peaks of both condition should be *more likely* to be a real binding
sites in condition A while comparing treatment A and control A, and in
condition B while comparing treatment B over control B, and also the
likelihood test while comparing treatment A and treatment B can't
decide which condition is stronger.

The likelihood function we used while comparing two conditions: ChIP
(enrichment) or control (chromatin bias) is:

```math
ln(LR) = x*(ln(x)-ln(y)) + y - x
```

Here $`LR`$ is the likelihood ratio, x is the signal (fragment pileup)
we observed in condition 1, and y is the signal in condition
2. And $`ln`$ is the natural logarithm.

Note: All regions on the same chromosome in the bedGraph file should
be continuous so only bedGraph files from MACS3 are acceptable.

## Command Line Options

The command line options for `bdgdiff` are defined in `/MACS3/Commands/bdgdiff_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `--t1`: MACS pileup bedGraph for condition 1. Incompatible with
  callpeak --SPMR output. REQUIRED
- `--t2`: MACS pileup bedGraph for condition 2. Incompatible with
  callpeak --SPMR output. REQUIRED
- `--c1`: MACS control lambda bedGraph for condition 1. Incompatible
  with callpeak --SPMR output. REQUIRED
- `--c2`: MACS control lambda bedGraph for condition 2. Incompatible
  with callpeak --SPMR output. REQUIRED
- `-C` or `--cutoff`: log10LR cutoff. Regions with signals lower than
  the cutoff will not be considered as enriched regions. DEFAULT: 3
  (likelihood ratio=1000)
- `-l` or `--min-len`: Minimum length of the differential region. Try
  a bigger value to remove small regions. DEFAULT: 200
- `-g` or `--max-gap`: Maximum gap to merge nearby differential
  regions. Consider a wider gap for broad marks. The maximum gap
  should be smaller than the minimum length (-g). DEFAULT: 100
- `--d1` or `--depth1`: Sequencing depth (# of non-redundant reads in
  million) for condition 1. It will be used together with --d2. See
  the description for --d2 below for how to assign them. Default: 1
- `--d2` or `--depth2`: Sequencing depth (# of non-redundant reads in
  million) for condition 2. It will be used together with --d1. DEPTH1
  and DEPTH2 will be used to calculate the scaling factor for each
  sample, to down-scale the larger sample to the level of the smaller
  one. For example, while comparing 10 million condition 1 and 20
  million condition 2, use --d1 10 --d2 20, then the pileup value in
  bedGraph for condition 2 will be divided by 2. Default: 1
- `--verbose`: Set the verbose level of runtime messages. 0: only show
  critical messages, 1: show additional warning messages, 2: show
  process information, 3: show debug messages. DEFAULT: 2
- `--outdir`: If specified, all output files will be written to that
  directory. Default: the current working directory
- `--o-prefix`: Output file prefix. Actual files will be named as
  PREFIX_cond1.bed, PREFIX_cond2.bed, and PREFIX_common.bed. Mutually
  exclusive with -o/--ofile.
- `-o` or `--ofile`: Output filenames. Must give three arguments in
  order: 1. file for unique regions in condition 1; 2. file for unique
  regions in condition 2; 3. file for common regions in both
  conditions. Note: mutually exclusive with --o-prefix.


## Example Usage

Here is an example of how to use the `bdgdiff` command:

```bash
macs3 bdgdiff -t1 treatment1.bedGraph -c1 control1.bedGraph -t2 treatment2.bedGraph -c2 control2.bedGraph --depth1 1.0 --depth2 1.0 -o output.bedGraph --minlen 500 --maxgap 1000 --cutoff 1.0
```

In this example, the program will call differential peaks from the two
pairs of treatment and control bedGraph files and write the result to
`output.bedGraph`. The depth of the first and second condition is set
to 1.0, the minimum length of differential peaks is set to 500, the
maximum gap between differential peaks is set to 1000, and the cutoff
for log10LR to call differential peaks is set to 1.0 (or likelihood
ratio 10).

## Step-by-step Instruction for calling differential peaks

In this chatper, we will describe how to use MACS3 to identify
differential regions by comparing pileup tracks of two conditions,
starting from the alignment files. Two modules will be involved:
`callpeak` and `bdgdiff` ( `predictd` is optional ). We will use human
ChIP-seq data as example, and filenames of raw data are:
cond1_ChIP.bam and cond1_Control.bam for condition 1; cond2_ChIP.bam
and cond2_Control.bam for condition 2.

### Step 1: Generate pileup tracks using callpeak module

Purpose of this step is to use `callpeak` with -B option to generate
bedGraph files for both conditions. There are several things to be
remember: 1. `--SPMR` is not compatible with `bdgdiff`, so avoid using
it; 2. prepare a pen to write down the number of non-redundant reads
of both conditions -- you will find such information in runtime
message or xls output from `callpeak`; 3. keep using the same
`--extsize` for both conditions ( you can get it from `predictd`
module).

To get a uniform extension size for running `callpeak`, run `predictd`:

```
 $ macs3 predictd -i cond1_ChIP.bam

 $ macs3 predictd -i cond2_ChIP.bam
```

An easy solution is to use the average of two 'fragment size'
predicted in `callpeak`, however any reasonable value will work. For
example, you can use `200` for most ChIP-seq datasets for
transcription factors, or ''147'' for most histone modification
ChIP-seq. The only requirement is that you have to keep using the same
extsize for the following commands:

```
 $ macs3 callpeak -B -t cond1_ChIP.bam -c cond1_Control.bam -n cond1 --nomodel --extsize 120
 
 $ macs3 callpeak -B -t cond2_ChIP.bam -c cond2_Control.bam -n cond2 --nomodel --extsize 120
```

Pay attention to runtime message, or extract the "tags after filtering in treatment" and "tags after filtering in control" lines from xls to see the effective sequencing depths for both conditions. In our previous command lines, '--to-large' is not used, so the effective sequencing depth is the smaller number of treatment and control. For example:

```
 $ egrep "tags after filtering in treatment|tags after filtering in control" cond1_peaks.xls
 # tags after filtering in treatment: 19291269
 # tags after filtering in control: 12914669

 $ egrep "tags after filtering in treatment|tags after filtering in control" cond2_peaks.xls
 # tags after filtering in treatment: 19962431
 # tags after filtering in control: 14444786
```

Then actual effective depths of condition 1 and 2 are: 12914669
and 14444786. Keep record of these two numbers and we will use them
later. After successfully running '''callpeak''', you will have
''cond1_treat_pileup.bdg'', ''cond1_control_lambda.bdg'',
''cond2_treat_pileup.bdg'', and ''cond2_control_lambda.bdg'' in the
working directory.

### Step 2: Call differential regions

The purpose of this step is to do a three ways comparisons to find out
where in the genome has differential enrichment between two
conditions. A basic requirement is that this region should be at least
enriched in either condition. A log10 likelihood ratio cutoff (C) will
be applied in this step. Three types of differential regions will be
reported: 1. those having more enrichment in condition 1 over
condition 2 ( cond1_ChIP > cond1_Control and cond1_ChIP > cond2_ChIP
); 2. those having more enrichment in condition 2 over condition 1 (
cond2_ChIP > cond2_Control and cond2_ChIP > cond1_ChIP ); those having
similar enrichment in both conditions ( cond1_ChIP > cond1_Control and
cond2_ChIP > cond2_Control and cond1_ChIP ≈ cond1_ChIP ).

Run this:

```
 $ macs3 bdgdiff --t1 cond1_treat_pileup.bdg --c1 cond1_control_lambda.bdg --t2 cond2_treat_pileup.bdg\
   --c2 cond2_control_lambda.bdg --d1 12914669 --d2 14444786 -g 60 -l 120 --o-prefix diff_c1_vs_c2
```

You will get the following three files in working directory:

 1. `diff_c1_vs_c2_c3.0_cond1.bed` This file stores regions that are
 highly enriched in condition 1 comparing to condition 2. The last
 column in the file represent the log10 likelihood ratio to show how
 likely the observed signal in condition 1 in this region is from
 condition 1 comparing to condition 2. Higher the value, bigger the
 difference.

 2. `diff_c1_vs_c2_c3.0_cond2.bed` This file stores regions that are
 highly enriched in condition 2 comparing to condition 1. The last
 column in the file represent the log10 likelihood ratio to show how
 likely the observed signal in condition 2 in this region is from
 condition 2 comparing to condition 1. Higher the value, bigger the
 difference.

 3. `diff_c1_vs_c2_c3.0_common.bed` This file stores regions that are
 highly enriched in both condition 1 and condition 2, and the
 difference between condition 1 and condition 2 is not
 significant. The last column in the file represent the difference
 between condition 1 and condition 2 in log10 likelihood ratios.
