# callvar

## Overview
The `callvar` command is part of the MACS3 suite of tools and is used
to call variants in giveqn peak regions from the alignment BAM
files. It is particularly useful in ChIP-Seq analysis where the
identification of genomic variants is required. 

## Detailed Description of usage

The `callvar` command takes in treatment and control BAM files along
with a bed file containing peak regions. The command identifies
variants in these regions using a multi-process approach, greatly
improving the speed and efficiency of variant calling. Please check
the section *Callvar Algorithm* for detail on this variant calling
algorithm. 

The `callvar` command assumes you have two types of BAM files. The
first type, what we call `TREAT`, is from DNA enrichment assay such as
ChIP-seq or ATAC-seq where the DNA fragments in the sequencing library
are enriched in certain genomics regions with potential allele biases;
the second type, called `CTRL` for control, is from genomic assay in
which the DNA enrichment is less biased in multiploid chromosomes and
more uniform across the whole genome (the later one is optional). In
order to run `callvar`, please sort (by coordinates) and index the BAM
files. 

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

Here is a brief overview of these options:

### Input files Options:

- `-b` or `--peak`: The peak regions in BED format, sorted by
  coordinates. This option is required. 
- `-t` or `--treatment`: The ChIP-seq/ATAC-seq treatment file in BAM
  format, sorted by coordinates. Make sure the .bai file is avaiable
  in the same directory. This option is required. 
- `-c` or `--control`: Optional control file in BAM format, sorted by
  coordinates. Make sure the .bai file is avaiable in the same
  directory. 

### Output Options:
- `--outdir`: The directory for all output files to be written
  to. Default: writes output files to the current working directory. 
- `-o` or `--ofile`: The output VCF file name. 
- `--verbose`: Set the verbose level of runtime messages. 0: only show
  critical messages, 1: show additional warning messages, 2: show
  process information, 3: show debug messages. DEFAULT: 2 

### Variant calling Options: 
- `-g` or `--gq-hetero`: The Genotype Quality score
  (-10log10((L00+L11)/(L01+L00+L11))) cutoff for Heterozygous allele
  type. Default is 0, or there is no cutoff on GQ. 
- `-G` or `--gq-homo`: The Genotype Quality score
  (-10log10((L00+L01)/(L01+L00+L11))) cutoff for Homozygous allele
  (not the same as reference) type. Default is 0, or there is no
  cutoff on GQ. 
- `-Q`: The cutoff for the quality score. Only consider bases with
  quality score greater than this value. Default is 20, which means
  Q20 or 0.01 error rate. 
- `-F` or `--fermi`: The option to control when to apply local
  assembly through Fermi. By default (set as 'auto'), while `callvar`
  detects any INDEL variant in a peak region, it will utilize Fermi to
  recover the actual DNA sequences to refine the read alignments. If
  set as 'on', Fermi will always be invoked. It can increase
  specificity, however sensivity and speed will be significantly
  lower. If set as 'off', Fermi won't be invoked at all. If so, speed
  and sensitivity can be higher but specificity will be significantly
  lower. 
- `--fermi-overlap`: The minimal overlap for fermi to initially
  assemble two reads. Must be between 1 and read length. A longer
  fermiMinOverlap is needed while read length is small (e.g. 30 for
  36bp read, but 33 for 100bp read may work). Default is 30. 
- `--top2alleles-mratio`: The reads for the top 2 most frequent
  alleles (e.g. a ref allele and an alternative allele) at a loci
  shouldn't be too few comparing to total reads mapped. The minimum
  ratio is set by this optoin. Must be a float between 0.5
  and 1. Default:0.8 which means at least 80% of reads contain the top
  2 alleles. 
- `--altallele-count`: The count of the alternative (non-reference)
  allele at a loci shouldn't be too few. By default, we require at
  least two reads support the alternative allele. Default:2 
- `--max-ar`: The maximum Allele-Ratio allowed while calculating
  likelihood for allele-specific binding. If we allow higher maxAR, we
  may mistakenly assign some homozygous loci as
  heterozygous. Default:0.95 

### Misc Options:
- `-m` or `--multiple-processing`: The CPU used for mutliple
  processing. Please note that, assigning more CPUs does not guarantee
  the process being faster. Creating too many parrallel processes need
  memory operations and may negate benefit from multi
  processing. Default: 1 

## Example Usage

Here is an example of how to use the `callvar` command:

```
macs3 callvar -b peaks.bed -t treatment.bam -c control.bam -o experiment1
```

In this example, the program will identify variants in the
`treatment.bam` file relative to the `control.bam` file. The name of
the experiment is 'experiment1'. All tags that pass quality filters
will be stored in a BAM file. 

## `callvar` Algorithm

![Callvar Algorithm](./callvar_algorithm.jpeg)

Functional sequencing assays which targeted at particular sequences,
such as ChIP-Seq, were thought to be unsuitable for *de novo*
variation predictions because their genome-wide sequencing coverage is
not as uniform as Whole Genome Sequencing (WGS). However, if we aim at
discover the variations and allele usage at the targeted genomic
regions, the coverage should be much higher and sufficient. We
therefore proposed a noval method to call the variants directly at the
called peaks by MACS3.

At each peak region, we extracted the reads and assembled the DNA
sequences using [fermi-lite](https://github.com/lh3/fermi-lite) , a
unitig graph based assembly algorithm developed by Heng Li. Then, we
aligned the unitigs (i.e., assembled short DNA sequences) to the
reference genome sequence using Smith-Waterman algorithm. Differences
between the reference sequence and the unitigs revealed possible SNVs
and INDELs. For each possible SNV or INDEL, we built a statistical
model incorporating the sequences and sequencing errors (base
qualities) from both treatment (ChIP) and control (genomic input) to
predict the most likely genotype using Bayesian Information Criterion
(BIC) among four allele types: homozygous loci (genotype 1/1),
heterozygous loci (genotype 0/1 or 1/2) with allele bias, and
heterozygous loci without allele bias. The detailed explanation of our
statistical model is as the following.  We retrieved the base quality
scores $\epsilon$, which represents sequencing errors, then we
calculated the likelihoods of each of the four types. We assumed the
independence of ChIP and control experiments so that the generalized
likelihood function was the product of the likelihood functions of
ChIP and control data:

$$L(\omega,\phi,g\_c,g\_i:D)=L(\omega,g\_c:D\_c)L(\phi,g\_i:D\_i)$$

where $`D_c`$ and $`D_i`$ represent the ChIP-Seq and control (e.g.,
genomic input) data observed at the position including base coverage
and base qualities. The parameter $\omega$ stands for the allele ratio
of allele A (chosen as the more abundant or stronger allele compared
with the others) from the ChIP-Seq data and $\phi$ represents the
allele ratio in the control. The parameter $`g_c`$ represents the
actual number of ChIPed DNA fragments containing allele A, which could
differ from the observed count $`r_(c,A)`$ considering that some
observations could be due to sequencing errors. The symbol $`g_i`$
represents the control analogously to $`g_c`$. We used $`r_c`$ to
denote the total number of observed allele A ($`r_(c,A)`$) and allele
B ($`r_(c,B)`$). We assumed the occurrence of the allele A ($`g_c`$)
is from a Bernoulli trial from $`r_c`$ with the allele ratio
$\omega$. The probability of observing the ChIP-Seq data at a certain
position under a given type is as follows:

```math
Pr(D_c|g_c,\omega) = Pr(D_c|g_c) =
 \sum^{r_{c,A}}_{j=1}\left((1-\epsilon_j)g_c/r_c+\epsilon_j(1-g_c/r_c)\right)\sum_{j=1}^{r_{c,B}}\left((1-\epsilon_j)(1-g_c/r_c)+\epsilon_j
 g_c/r_c\right)
 ```

where ϵ_j represents the sequencing error of the base showing
difference with reference genome in case of mismatch (corresponding to
SNV) and insertion. In case of deletion, the sequencing errors from
the two bases on sequenced read surrounding the deletion would be
considered. We modeled the control data in the similar way. We
assessed the likelihood functions of the 4 major type using the
following parameters: ω=1,φ=1,g_c=r_(c,0),g_i=r_(i,0) for A/A
genotype; ω=0,φ=0,g_c=0,g_i=0 for B/B genotype, ω=0.5,φ=0.5 and
g_c,g_i as free variables for A/B genotype with unbiased binding;
φ=0.5 and 〖ω,g〗_c,g_i as free variables for A/B genotype with biased
binding or allele usage. Next, we applied the Bayesian Information
Criterion (BIC) to select the best type as our prediction with the
minimal BIC value among the 4 models. If the best type was either
“A/B, noAS” or “A/B, AS”, we concluded that the genotype was
heterozygous (A/B). We consider two types of data from the same assay
independently: ChIP sample that can have biased allele usage, and
control sample that won’t have biased allele usage. So that in case
control is not available, such as in ATAC-Seq assay, our model can
still work. Furthermore, in case a good quality WGS is available, it
can be regarded as the control sample and be inserted into our
calculation to further increase the sensitivity. 
