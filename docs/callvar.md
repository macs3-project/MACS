# Callvar

## Overview
The `callvar` command is part of the MACS3 suite of tools and is used to call variants in given peak regions from the alignment BAM files. It is particularly useful in ChIP-Seq analysis where the identification of genomic variants is required.

## Detailed Description

The `callvar` command takes in treatment and control BAM files along with a bed file containing peak regions. The command identifies variants in these regions using a multi-process approach, greatly improving the speed and efficiency of variant calling.

The `callvar` command assumes you have two types of BAM files. The first type, what we call `TREAT`, is from DNA enrichment assay such as ChIP-seq or ATAC-seq where the DNA fragments in the sequencing library are enriched in certain genomics regions with potential allele biases; the second type, called `CTRL` for control, is from genomic assay in which the DNA enrichment is less biased in multiploid chromosomes and more uniform across the whole genome (the later one is optional). In order to run `callvar`, please sort (by coordinates) and index the BAM files.
Example:

1. Sort the BAM file:
    `$ samtools sort TREAT.bam -o TREAT_sorted.bam`
    `$ samtools sort CTRL.bam -o CTRL_sorted.bam`
2. Index the BAM file:
    `$ samtools index TREAT_sorted.bam`
    `$ samtools index CTRL_sorted.bam`
3. Make sure .bai files are available:
    `$ ls TREAT_sorted.bam.bai`
    `$ ls CTRL_sorted.bam.bai`

To call variants:
    `$ macs3 callvar -b peaks.bed -t TREAT_sorted.bam -c CTRL_sorted.bam -o peaks.vcf`

## Command Line Options

The command line options for `callvar` are defined in `/MACS3/Commands/callvar_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

### Input files arguments:
- `-b` or `--peak`: The peak regions in BED format, sorted by coordinates. This option is required.
- `-t` or `--treatment`: The ChIP-seq/ATAC-seq treatment file in BAM format, sorted by coordinates. Make sure the .bai file is avaiable in the same directory. This option is required.
- `-c` or `--control`: Optional control file in BAM format, sorted by coordinates. Make sure the .bai file is avaiable in the same directory.

### Output arguments:
- `--outdir`: The directory for all output files to be written to. Default: writes output files to the current working directory.
- `-o` or `--ofile`: The output VCF file name.
- `--verbose`: The directory for all output files to be written to. Default: writes output files to the current working directory.

### Variant calling arguments: 
- `-g` or `--gq-hetero`: The Genotype Quality score (-10log10((L00+L11)/(L01+L00+L11))) cutoff for Heterozygous allele type. Default is 0, or there is no cutoff on GQ.
- `-G` or `--gq-homo`: The Genotype Quality score (-10log10((L00+L01)/(L01+L00+L11))) cutoff for Homozygous allele (not the same as reference) type. Default is 0, or there is no cutoff on GQ.
- `-Q`: The cutoff for the quality score. Only consider bases with quality score greater than this value. Default is 20, which means Q20 or 0.01 error rate.
- `-F` or `--fermi`: The option to control when to apply local assembly through Fermi. By default (set as 'auto'), while SAPPER detects any INDEL variant in a peak region, it will utilize Fermi to recover the actual DNA sequences to refine the read alignments. If set as 'on', Fermi will be always invoked. It can increase specificity however sensivity and speed will be significantly lower. If set as 'off', Fermi won't be invoked at all. If so, speed and sensitivity can be higher but specificity will be significantly lower.
- `--fermi-overlap`: The minimal overlap for fermi to initially assemble two reads. Must be between 1 and read length. A longer fermiMinOverlap is needed while read length is small (e.g. 30 for 36bp read, but 33 for 100bp read may work). Default is 30.
- `--top2alleles-mratio`: The reads for the top 2 most frequent alleles (e.g. a ref allele and an alternative allele) at a loci shouldn't be too few comparing to total reads mapped. The minimum ratio is set by this optoin. Must be a float between 0.5 and 1. Default:0.8 which means at least 80% of reads contain the top 2 alleles.
- `--altallele-count`: The count of the alternative (non-reference) allele at a loci shouldn't be too few. By default, we require at least two reads support the alternative allele. Default:2
- `--max-ar`: The maximum Allele-Ratio allowed while calculating likelihood for allele-specific binding. If we allow higher maxAR, we may mistakenly assign some homozygous loci as heterozygous. Default:0.95

### Misc arguments:
- `-m` or `--multiple-processing`: The CPU used for mutliple processing. Please note that, assigning more CPUs does not guarantee the process being faster. Creating too many parrallel processes need memory operations and may negate benefit from multi processing. Default: 1


## Example Usage

Here is an example of how to use the `callvar` command:

```bash
macs3 callvar -b peaks.bed -t treatment.bam -c control.bam -o experiment1
```

In this example, the program will identify variants in the `treatment.bam` file relative to the `control.bam` file. The name of the experiment is 'experiment1'. All tags that pass quality filters will be stored in a BAM file.