# MACS: Model-based Analysis for ChIP-Seq

![Status](https://img.shields.io/pypi/status/macs3.svg) ![License](https://img.shields.io/github/license/macs3-project/MACS) ![Programming languages](https://img.shields.io/github/languages/top/macs3-project/MACS) ![CI x64](https://github.com/macs3-project/MACS/workflows/CI%20x64/badge.svg) ![CI non x64](https://github.com/macs3-project/MACS/workflows/CI%20non%20x64,%20python%203.7/badge.svg)

[![PyPI download](https://img.shields.io/pypi/dm/macs3?label=pypi%20downloads)](https://pypistats.org/packages/macs3) [![Bioconda download](https://img.shields.io/conda/dn/bioconda/macs3?label=bioconda%20downloads)](https://anaconda.org/bioconda/macs3)

Latest Release:
* Github: [![Github Release](https://img.shields.io/github/v/release/macs3-project/MACS)](https://github.com/macs3-project/MACS/releases)
* PyPI: [![PyPI Release](https://img.shields.io/pypi/v/macs3.svg) ![PyPI Python Version](https://img.shields.io/pypi/pyversions/MACS3) ![PyPI Format](https://img.shields.io/pypi/format/macs3)](https://pypi.org/project/macs3/)
* Bioconda: [![Bioconda Release](https://img.shields.io/conda/v/bioconda/macs3) ![Bioconda Platform](https://img.shields.io/conda/pn/bioconda/macs3)](https://anaconda.org/bioconda/macs3)
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

**Please note that current MACS3 is still in alpha stage, although we
utilize Github Action to implement the CI (Continous Integration) to
make sure that the main branch passes unit testing on certain
functions and subcommands. More new featuer will be added soon.**

## Recent Changes for MACS (3.0.0a1)

### 3.0.0a1
	* Features
	
	1) Speed/memory optimization, including using the cykhash to
    replace python dictionary

	2) Code cleanup

	3) Unit testing

	4) R wrappers for MACS

    5) Switching to Github Action for CI, support multi-arch testing.

## Install

The common way to install MACS is through
[PYPI](https://pypi.org/project/macs3/)) or
[conda](https://anaconda.org/bioconda/macs3). Please check the
[INSTALL](./docs/INSTALL.md) document for detail.

## Usage

There are currently twelve functions available in MAC2S serving as sub-commands.

```
macs3 [-h] [--version]
    {callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak}
```

Example for regular peak calling: `macs2 callpeak -t ChIP.bam -c
Control.bam -f BAM -g hs -n test -B -q 0.01`

Example for broad peak calling: `macs2 callpeak -t ChIP.bam -c
Control.bam --broad -g hs --broad-cutoff 0.1`

Please click on the link to see the detail description of the subcommands.

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

For advanced usage, for example, to run `macs3` in a modular way,
please read the [advanced usage](./docs/advanced_usage.md). There is a
[Q&A](./docs/qa.md) document where we collected some common questions
from users.

## Contribute

Please read our [CODE OF CONDUCT](./CODE_OF_CONDUCT.md) and
[How to contribute](./CONTRIBUTING.md) documents.

## Ackowledgement

MACS3 project is sponsored by
[CZI EOSS](https://chanzuckerberg.com/eoss/). And we particularly want
to thank the user community for their supports, feedbacks and
contributions over the years.

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
