# MACS: Model-based Analysis for ChIP-Seq

![Status](https://img.shields.io/pypi/status/macs2.svg) [![License](https://img.shields.io/github/license/taoliu/MACS)](https://github.com/taoliu/MACS/blob/master/LICENSE) ![Programming languages](https://img.shields.io/github/languages/top/taoliu/MACS) ![Commit activity](https://img.shields.io/github/commit-activity/m/taoliu/MACS) [![TravisCI Build Status](https://img.shields.io/travis/com/taoliu/MACS/master)](https://travis-ci.com/taoliu/MACS)

[![PyPI download](https://img.shields.io/pypi/dm/macs2?label=pypi%20downloads)](https://pypistats.org/packages/macs2) [![Bioconda download](https://img.shields.io/conda/dn/bioconda/macs2?label=bioconda%20downloads)](https://anaconda.org/bioconda/macs2)

Latest Release:
* Github: [![Github Release](https://img.shields.io/github/v/release/taoliu/MACS)](https://github.com/taoliu/MACS/releases)
* PyPI: [![PyPI Release](https://img.shields.io/pypi/v/macs2.svg) ![PyPI Python Version](https://img.shields.io/pypi/pyversions/MACS2) ![PyPI Format](https://img.shields.io/pypi/format/macs2)](https://pypi.org/project/macs2/)
* Bioconda: [![Bioconda Release](https://img.shields.io/conda/v/bioconda/macs2) ![Bioconda Platform](https://img.shields.io/conda/pn/bioconda/macs2)](https://anaconda.org/bioconda/macs2)
* Debian Med: [![Debian Stable](https://img.shields.io/debian/v/macs/stable?label=debian%20stable)](https://packages.debian.org/stable/macs) [![Debian Unstable](https://img.shields.io/debian/v/macs/sid?label=debian%20sid)](https://packages.debian.org/sid/macs)

## Introduction

With the improvement of sequencing techniques, chromatin
immunoprecipitation followed by high throughput sequencing (ChIP-Seq)
is getting popular to study genome-wide protein-DNA interactions. To
address the lack of powerful ChIP-Seq analysis method, we presented
the **M**odel-based **A**nalysis of **C**hIP-**S**eq (MACS), for
identifying transcript factor binding sites. MACS captures the
influence of genome complexity to evaluate the significance of
enriched ChIP regions and MACS improves the spatial resolution of
binding sites through combining the information of both sequencing tag
position and orientation. MACS can be easily used for ChIP-Seq data
alone, or with a control sample with the increase of
specificity. Moreover, as a general peak-caller, MACS can also be
applied to any "DNA enrichment assays" if the question to be asked is
simply: *where we can find significant reads coverage than the random
background*.

## Recent Changes for MACS (2.2.8)

### 2.2.8
   * Bug fix:

	MACS2 typo in 'setup.py': 'numpy>=>=1.17' -> 'numpy>=1.17'
	PR #543, issues #535, #541, #544

	Now test on and support Python 3.6/3.7/3.8/3.9/3.10/3.11

### 2.2.7.1
	* hotfix: 

	Add 'wheel' and 'pip' to pyproject.toml so that `pip install` can 
	work.


### 2.2.7
	* Bugs fixed

	1) MACS2 has been tested on multiple architectures to make sure it
	can successfully generate consistent results. Currently the
	supported architectures are: AMD64, ARM64, i386, PPC64LE, and
	S390X. Thanks to @mr-c, @junaruga, and @tillea! Related to issue
	#340, #349, #351, and #359; to PR #348, #350, #360, #361, #367,
	and #370. The lesson is that if the project is built on Cython and
	is aimed at memory efficiency, we should specifically define all
	int/float types in pyx files such as int8_t or uint32_t using
	either libc or numpy (c version) instead of relying on Cython
	types such as short, long, double.

	2) MACS2 setup script will check numpy and install numpy if
	necessary. PR #378, issue #364

	3) `bdgbroadcall` command will correctly add the score column (5th
	column). The score (5th) column contains 10 times of the average
	score in the broad region. PR #373, issue #362

	4) The missing test on `bdgopt` subcommand has been added. PR #363

	5) The obsolete option `--ratio` from `callpeak` subcommand has
	been removed. PR #369, issue #366

	6) Fixed the incorrect description in README on the 'maximum
	length of broad region is 4 times of d' to 'maximum gap for
	merging broad regions is 4 times of tag size by default'. PR #380,
	issue #365.

	* Other

	1) CODE OF CONDUCT document has been added to MACS2 github
	repository. PR #358

### 2.2.6
	* New Features 

	1) Speed up MACS2. Some programming tricks and code cleanup. The 
	filter_dup function replaces separate_dups. The later one was 
	implemented for potentially putting back duplicate reads in 
	certain downstream analysis. However such analysis hasn't been 
	implemented. Optimize speed of writing bedGraph files. Optimize 
	BAM and BAMPE parsing with pointer casting instead of python 
	unpack. 

	2) The comment lines in the headers of BED or SAM files will be
	correctly skipped. However, MACS2 won't check comment lines in the
	middle of the file.

	* Bugs fixed 

	1) Cutoff-analysis in callpeak command. #341
	
	2) Issues related to SAMParser and three ELAND Parsers are
	fixed. #347
	
	* Other 

	1) cmdlinetest script in test/ folder has been updated to: 1. test 
	cutoff-analysis with callpeak cmd; 2. output the 2 lines before 
	and after the error or warning message during tests; 3. output 
	only the first 10 lines if the difference between test result and 
	standard result can be found; 4. prockreport monitor CPU time and 
	memory usage in 1 sec interval -- a bit more accurate.
	
	2) Python3.5 support is removed. Now MACS2 requires Python>=3.6.

## Install

Please check the file 'INSTALL.md' in the distribution.

## Usage

```
macs2 [-h] [--version]
    {callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak}
```

Example for regular peak calling: `macs2 callpeak -t ChIP.bam -c
Control.bam -f BAM -g hs -n test -B -q 0.01`

Example for broad peak calling: `macs2 callpeak -t ChIP.bam -c
Control.bam --broad -g hs --broad-cutoff 0.1`

There are twelve functions available in MAC2S serving as sub-commands.

Subcommand | Description
-----------|----------
`callpeak` | Main MACS2 Function to call peaks from alignment results.
`bdgpeakcall` | Call peaks from bedGraph output.
`bdgbroadcall` | Call broad peaks from bedGraph output.
`bdgcmp` | Comparing two signal tracks in bedGraph format.
`bdgopt` | Operate the score column of bedGraph file.
`cmbreps` | Combine BEDGraphs of scores from replicates.
`bdgdiff` | Differential peak detection based on paired four bedGraph files.
`filterdup` | Remove duplicate reads, then save in BED/BEDPE format.
`predictd` | Predict d or fragment size from alignment results.
`pileup` | Pileup aligned reads (single-end) or fragments (paired-end)
`randsample` | Randomly choose a number/percentage of total reads.
`refinepeak` | Take raw reads alignment, refine peak summits.

We only cover `callpeak` subcommand in this document. Please use
`macs2 COMMAND -h` to see the detail description for each option of
each subcommand.

### Call peaks

This is the main function in MACS2. It can be invoked by `macs2
callpeak` . If you type this command with `-h`, you will see a full
description of command-line options. Here we only list the essentials.

#### Essential Options

##### `-t`/`--treatment FILENAME`

This is the only REQUIRED parameter for MACS. The file can be in any
supported format -- see detail in the `--format` option. If you have
more than one alignment file, you can specify them as `-t A B C`. MACS
will pool up all these files together.

##### `-c`/`--control`

The control, genomic input or mock IP data file. Please follow the
same direction as for `-t`/`--treatment`.

##### `-n`/`--name`

The name string of the experiment. MACS will use this string NAME to
create output files like `NAME_peaks.xls`, `NAME_negative_peaks.xls`,
`NAME_peaks.bed` , `NAME_summits.bed`, `NAME_model.r` and so on. So
please avoid any confliction between these filenames and your existing
files.

##### `--outdir`

MACS2 will save all output files into the specified folder for this
option. A new folder will be created if necessary.

##### `-f`/`--format FORMAT`

Format of tag file can be `ELAND`, `BED`, `ELANDMULTI`, `ELANDEXPORT`,
`SAM`, `BAM`, `BOWTIE`, `BAMPE`, or `BEDPE`. Default is `AUTO` which
will allow MACS to decide the format automatically. `AUTO` is also
useful when you combine different formats of files. Note that MACS
can't detect `BAMPE` or `BEDPE` format with `AUTO`, and you have to
implicitly specify the format for `BAMPE` and `BEDPE`.

Nowadays, the most common formats are `BED` or `BAM` (including
`BEDPE` and `BAMPE`). Our recommendation is to convert your data to
`BED` or `BAM` first.

Also, MACS2 can detect and read gzipped file. For example, `.bed.gz`
file can be directly used without being uncompressed with `--format
BED`.

Here are detailed explanation of the recommanded formats:

###### `BED`

The BED format can be found at [UCSC genome browser
website](http://genome.ucsc.edu/FAQ/FAQformat#format1).

The essential columns in BED format input are the 1st column
`chromosome name`, the 2nd `start position`, the 3rd `end position`,
and the 6th, `strand`.

Note that, for `BED` format, the 6th column of strand information is
required by MACS. And please pay attention that the coordinates in BED
format are zero-based and half-open. See more detail at
[UCSC site](http://genome.ucsc.edu/FAQ/FAQtracks#tracks1).

###### `BAM`/`SAM`

If the format is `BAM`/`SAM`, please check the definition in
(http://samtools.sourceforge.net/samtools.shtml).  If the `BAM` file is
generated for paired-end data, MACS will only keep the left mate(5'
end) tag. However, when format `BAMPE` is specified, MACS will use the
real fragments inferred from alignment results for reads pileup.

###### `BEDPE` or `BAMPE`

A special mode will be triggered while the format is specified as
`BAMPE` or `BEDPE`. In this way, MACS2 will process the `BAM` or `BED`
files as paired-end data. Instead of building a bimodal distribution
of plus and minus strand reads to predict fragment size, MACS2 will
use actual insert sizes of pairs of reads to build fragment pileup.

The `BAMPE` format is just a `BAM` format containing paired-end alignment
information, such as those from `BWA` or `BOWTIE`.

The `BEDPE` format is a simplified and more flexible `BED` format,
which only contains the first three columns defining the chromosome
name, left and right position of the fragment from Paired-end
sequencing. Please note, this is NOT the same format used by
`BEDTOOLS`, and the `BEDTOOLS` version of `BEDPE` is actually not in a
standard `BED` format. You can use MACS2 subcommand `randsample` to
convert a `BAM` file containing paired-end information to a `BEDPE`
format file:

```
macs2 randsample -i the_BAMPE_file.bam -f BAMPE -p 100 -o the_BEDPE_file.bed
```

##### `-g`/`--gsize`

PLEASE assign this parameter to fit your needs!

It's the mappable genome size or effective genome size which is
defined as the genome size which can be sequenced. Because of the
repetitive features on the chromosomes, the actual mappable genome
size will be smaller than the original size, about 90% or 70% of the
genome size. The default *hs* -- 2.7e9 is recommended for human
genome. Here are all precompiled parameters for effective genome size:

 * hs: 2.7e9
 * mm: 1.87e9
 * ce: 9e7
 * dm: 1.2e8

Users may want to use k-mer tools to simulate mapping of Xbps long
reads to target genome, and to find the ideal effective genome
size. However, usually by taking away the simple repeats and Ns from
the total genome, one can get an approximate number of effective
genome size. A slight difference in the number won't cause a big
difference of peak calls, because this number is used to estimate a
genome-wide noise level which is usually the least significant one
compared with the *local biases* modeled by MACS.

##### `-s`/`--tsize`

The size of sequencing tags. If you don't specify it, MACS will try to
use the first 10 sequences from your input treatment file to determine
the tag size. Specifying it will override the automatically determined
tag size.

##### `-q`/`--qvalue`

The q-value (minimum FDR) cutoff to call significant regions. Default
is 0.05. For broad marks, you can try 0.05 as the cutoff. Q-values are
calculated from p-values using the Benjamini-Hochberg procedure.

##### `-p`/`--pvalue`

The p-value cutoff. If `-p` is specified, MACS2 will use p-value instead
of q-value.

##### `--min-length`, `--max-gap`

These two options can be used to fine-tune the peak calling behavior
by specifying the minimum length of a called peak and the maximum
allowed a gap between two nearby regions to be merged. In other words,
a called peak has to be longer than `min-length`, and if the distance
between two nearby peaks is smaller than `max-gap` then they will be
merged as one. If they are not set, MACS2 will set the DEFAULT value
for `min-length` as the predicted fragment size `d`, and the DEFAULT
value for `max-gap` as the detected read length. Note, if you set a
`min-length` value smaller than the fragment size, it may have NO
effect on the result. For broad peak calling with `--broad` option
set, the DEFAULT `max-gap` for merging nearby stronger peaks will be
the same as narrow peak calling, and 4 times of the `max-gap` will be
used to merge nearby weaker (broad) peaks. You can also use
`--cutoff-analysis` option with the default setting, and check the
column `avelpeak` under different cutoff values to decide a reasonable
`min-length` value.

##### `--nolambda`

With this flag on, MACS will use the background lambda as local
lambda. This means MACS will not consider the local bias at peak
candidate regions.

##### `--slocal`, `--llocal`

These two parameters control which two levels of regions will be
checked around the peak regions to calculate the maximum lambda as
local lambda. By default, MACS considers 1000bp for small local
region(`--slocal`), and 10000bps for large local region(`--llocal`)
which captures the bias from a long-range effect like an open
chromatin domain. You can tweak these according to your
project. Remember that if the region is set too small, a sharp spike
in the input data may kill a significant peak.

##### `--nomodel`

While on, MACS will bypass building the shifting model.

##### `--extsize`

While `--nomodel` is set, MACS uses this parameter to extend reads in
5'->3' direction to fix-sized fragments. For example, if the size of
the binding region for your transcription factor is 200 bp, and you
want to bypass the model building by MACS, this parameter can be set
as 200. This option is only valid when `--nomodel` is set or when MACS
fails to build model and `--fix-bimodal` is on.

##### `--shift`

Note, this is NOT the legacy `--shiftsize` option which is replaced by
`--extsize`! You can set an arbitrary shift in bp here. Please Use
discretion while setting it other than the default value (0). When
`--nomodel` is set, MACS will use this value to move cutting ends (5')
then apply `--extsize` from 5' to 3' direction to extend them to
fragments. When this value is negative, ends will be moved toward
3'->5' direction, otherwise 5'->3' direction. Recommended to keep it
as default 0 for ChIP-Seq datasets, or -1 * half of *EXTSIZE* together
with `--extsize` option for detecting enriched cutting loci such as
certain DNAseI-Seq datasets. Note, you can't set values other than 0
if the format is BAMPE or BEDPE for paired-end data. The default is 0.

Here are some examples for combining `--shift` and `--extsize`:

1. To find enriched cutting sites such as some DNAse-Seq datasets. In
this case, all 5' ends of sequenced reads should be extended in both
directions to smooth the pileup signals. If the wanted smoothing
window is 200bps, then use `--nomodel --shift -100 --extsize 200`.

2. For certain nucleosome-seq data, we need to pile up the centers of
nucleosomes using a half-nucleosome size for wavelet analysis
(e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about
147bps, this option can be used: `--nomodel --shift 37 --extsize 73`.

##### `--keep-dup`

It controls the MACS behavior towards duplicate tags at the exact same
location -- the same coordination and the same strand. The default
`auto` option makes MACS calculate the maximum tags at the exact same
location based on binomial distribution using 1e-5 as p-value cutoff;
and the `all` option keeps every tag.  If an integer is given, at most
this number of tags will be kept at the same location. The default is
to keep one tag at the same location. Default: 1

##### `--broad`

When this flag is on, MACS will try to composite broad regions in
BED12 ( a gene-model-like format ) by putting nearby highly enriched
regions into a broad region with loose cutoff. The broad region is
controlled by another cutoff through `--broad-cutoff`. Please note
that, the `max-gap` value for merging nearby weaker/broad peaks is 4
times of `max-gap` for merging nearby stronger peaks. The later one
can be controlled by `--max-gap` option, and by default it is the
average fragment/insertion length in the PE data. DEFAULT: False

##### `--broad-cutoff`

Cutoff for the broad region. This option is not available unless
`--broad` is set. If `-p` is set, this is a p-value cutoff, otherwise,
it's a q-value cutoff.  DEFAULT: 0.1

##### `--scale-to <large|small>`

When set to `large`, linearly scale the smaller dataset to the same
depth as the larger dataset. By default or being set as `small`, the
larger dataset will be scaled towards the smaller dataset. Beware, to
scale up small data would cause more false positives.

##### `-B`/`--bdg`

If this flag is on, MACS will store the fragment pileup, control
lambda in bedGraph files. The bedGraph files will be stored in the
current directory named `NAME_treat_pileup.bdg` for treatment data,
`NAME_control_lambda.bdg` for local lambda values from control.

##### `--call-summits`

MACS will now reanalyze the shape of signal profile (p or q-score
depending on the cutoff setting) to deconvolve subpeaks within each
peak called from the general procedure. It's highly recommended to
detect adjacent binding events. While used, the output subpeaks of a
big peak region will have the same peak boundaries, and different
scores and peak summit positions.

##### `--buffer-size`

MACS uses a buffer size for incrementally increasing internal array
size to store reads alignment information for each chromosome or
contig. To increase the buffer size, MACS can run faster but will
waste more memory if certain chromosome/contig only has very few
reads. In most cases, the default value 100000 works fine. However, if
there are a large number of chromosomes/contigs in your alignment and
reads per chromosome/contigs are few, it's recommended to specify a
smaller buffer size in order to decrease memory usage (but it will
take longer time to read alignment files). Minimum memory requested
for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE *
8 Bytes. DEFAULT: 100000

#### Output files

1. `NAME_peaks.xls` is a tabular file which contains information about
   called peaks. You can open it in excel and sort/filter using excel
   functions. Information include:
   
    - chromosome name
    - start position of peak
    - end position of peak
    - length of peak region
    - absolute peak summit position
    - pileup height at peak summit
    - -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then
      this value should be 10)
    - fold enrichment for this peak summit against random Poisson
      distribution with local lambda,
    - -log10(qvalue) at peak summit
   
   Coordinates in XLS is 1-based which is different from BED
   format. When `--broad` is enabled for broad peak calling, the
   pileup, p-value, q-value, and fold change in the XLS file will be
   the mean value across the entire peak region, since peak summit
   won't be called in broad peak calling mode.

2. `NAME_peaks.narrowPeak` is BED6+4 format file which contains the
   peak locations together with peak summit, p-value, and q-value. You
   can load it to the UCSC genome browser. Definition of some specific
   columns are:
   
   - 5th: integer score for display. It's calculated as
     `int(-10*log10pvalue)` or `int(-10*log10qvalue)` depending on
     whether `-p` (pvalue) or `-q` (qvalue) is used as score
     cutoff. Please note that currently this value might be out of the
     [0-1000] range defined in [UCSC ENCODE narrowPeak
     format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12). You
     can let the value saturated at 1000 (i.e. p/q-value = 10^-100) by
     using the following 1-liner awk: `awk -v OFS="\t"
     '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak`
   - 7th: fold-change at peak summit
   - 8th: -log10pvalue at peak summit
   - 9th: -log10qvalue at peak summit
   - 10th: relative summit position to peak start
   
   The file can be loaded directly to the UCSC genome browser. Remove
   the beginning track line if you want to analyze it by other tools.

3. `NAME_summits.bed` is in BED format, which contains the peak
   summits locations for every peak. The 5th column in this file is
   the same as what is in the `narrowPeak` file. If you want to find
   the motifs at the binding sites, this file is recommended. The file
   can be loaded directly to the UCSC genome browser. Remove the
   beginning track line if you want to analyze it by other tools.

4. `NAME_peaks.broadPeak` is in BED6+3 format which is similar to
   `narrowPeak` file, except for missing the 10th column for
   annotating peak summits. This file and the `gappedPeak` file will
   only be available when `--broad` is enabled. Since in the broad
   peak calling mode, the peak summit won't be called, the values in
   the 5th, and 7-9th columns are the mean value across all positions
   in the peak region. Refer to `narrowPeak` if you want to fix the
   value issue in the 5th column.

5. `NAME_peaks.gappedPeak` is in BED12+3 format which contains both
   the broad region and narrow peaks. The 5th column is the score for
   showing grey levels on the UCSC browser as in `narrowPeak`. The 7th
   is the start of the first narrow peak in the region, and the 8th
   column is the end. The 9th column should be RGB color key, however,
   we keep 0 here to use the default color, so change it if you
   want. The 10th column tells how many blocks including the starting
   1bp and ending 1bp of broad regions. The 11th column shows the
   length of each block and 12th for the start of each block. 13th:
   fold-change, 14th: *-log10pvalue*, 15th: *-log10qvalue*. The file can
   be loaded directly to the UCSC genome browser. Refer to
   `narrowPeak` if you want to fix the value issue in the 5th column.

6. `NAME_model.r` is an R script which you can use to produce a PDF
   image of the model based on your data. Load it to R by:

   `$ Rscript NAME_model.r`

   Then a pdf file `NAME_model.pdf` will be generated in your current
   directory. Note, R is required to draw this figure.

7. The `NAME_treat_pileup.bdg` and `NAME_control_lambda.bdg` files are
   in bedGraph format which can be imported to the UCSC genome browser
   or be converted into even smaller bigWig files. The
   `NAME_treat_pielup.bdg` contains the pileup signals (normalized
   according to `--scale-to` option) from ChIP/treatment sample. The
   `NAME_control_lambda.bdg` contains local biases estimated for each
   genomic location from the control sample, or from treatment sample
   when the control sample is absent. The subcommand `bdgcmp` can be
   used to compare these two files and make a bedGraph file of scores
   such as p-value, q-value, log-likelihood, and log fold changes.

## Other useful links

 * [Cistrome](http://cistrome.org/ap/)
 * [bedTools](http://code.google.com/p/bedtools/)
 * [UCSC toolkits](http://hgdownload.cse.ucsc.edu/admin/exe/)

## Tips of fine-tuning peak calling

There are several subcommands within MACSv2 package to fine-tune or
customize your analysis:

1. `bdgcmp` can be used on `*_treat_pileup.bdg` and
   `*_control_lambda.bdg` or bedGraph files from other resources to
   calculate the score track.

2. `bdgpeakcall` can be used on `*_treat_pvalue.bdg` or the file
   generated from bdgcmp or bedGraph file from other resources to call
   peaks with given cutoff, maximum-gap between nearby mergeable peaks
   and a minimum length of peak. bdgbroadcall works similarly to
   bdgpeakcall, however, it will output `_broad_peaks.bed` in BED12
   format.

3. Differential calling tool -- `bdgdiff`, can be used on 4 bedGraph
   files which are scores between treatment 1 and control 1, treatment
   2 and control 2, treatment 1 and treatment 2, treatment 2 and
   treatment 1. It will output consistent and unique sites according
   to parameter settings for minimum length, the maximum gap and
   cutoff.

4. You can combine subcommands to do a step-by-step peak calling. Read
   detail at [MACS2
   wikipage](https://github.com/taoliu/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands)
