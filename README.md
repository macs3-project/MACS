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

**Please note that current MACS3 is still in alpha stage. However, we
utilize Github Action to implement the CI (Continous Integration) to
make sure that the main branch passes unit testing on certain
functions and subcommands to reproduce the correct outputs. We will
add more new features in the future.**

## Recent Changes for MACS (3.0.0b1)

### 3.0.0b1
        The first beta version of MACS3, with HMMRATAC feature recently added.
	   
	* New features from alpha7:

	1) HMMRATAC module is added
	HMMRATAC is a dedicated software to analyze ATAC-seq data. The
	basic idea behind HMMRATAC is to digest ATAC-seq data according to
	the fragment length of read pairs into four signal tracks: short
	fragments, mononucleosomal fragments, di-nucleosomal fragments and
	tri-nucleosomal fragments. Then integrate the four tracks again
	using Hidden Markov Model to consider three hidden states: open
	region, nucleosomal region, and background region. The orginal
	paper was published in 2019 written in JAVA, by Evan Tarbell. We
	implemented it in Python/Cython and optimize the whole process
	using existing MACS functions and hmmlearn. Now it can run much
	faster than the original JAVA version. Note: evaluation of the
	peak calling results is underway.
	
	2) Multiple updates regarding dependencies, anaconda built, CI/CD
	process.

## Install

The common way to install MACS is through
[PYPI](https://pypi.org/project/macs3/)) or
[conda](https://anaconda.org/bioconda/macs3). Please check the
[INSTALL](./docs/INSTALL.md) document for detail.

MACS3 has been tested in CI for every push and PR in the following
architectures:

 * x86_64
 * aarch64
 * armv7
 * ppc64le
 * s390x 

## Usage

Example for regular peak calling on TF ChIP-seq:

`macs3 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01`

Example for broad peak calling on Histone Mark ChIP-seq:

`macs3 callpeak -t ChIP.bam -c Control.bam --broad -g hs --broad-cutoff 0.1`

Example for peak calling on ATAC-seq (paired-end mode):

`macs3 callpeak -f BAMPE -t ATAC.bam -g hs -n test -B -q 0.01`

There are currently twelve functions available in MAC3 serving as
sub-commands. Please click on the link to see the detail description
of the subcommands.

Subcommand | Description
-----------|----------
[`callpeak`](./docs/callpeak.md) | Main MACS3 Function to call peaks from alignment results.
[`bdgpeakcall`](./docs/bdgpeakcall.md) | Call peaks from bedGraph output.
[`bdgbroadcall`](./docs/bdgbroadcall.md) | Call broad peaks from bedGraph output.
[`bdgcmp`](./docs/bdgcmp.md) | Comparing two signal tracks in bedGraph format.
[`bdgopt`](./docs/bdgopt.md) | Operate the score column of bedGraph file.
[`cmbreps`](./docs/cmbreps.md) | Combine BEDGraphs of scores from replicates.
[`bdgdiff`](./docs/bdgdiff.md) | Differential peak detection based on paired four bedGraph files.
[`filterdup`](./docs/filterdup.md) | Remove duplicate reads, then save in BED/BEDPE format.
[`predictd`](./docs/predictd.md) | Predict d or fragment size from alignment results.
[`pileup`](./docs/pileup.md) | Pileup aligned reads (single-end) or fragments (paired-end)
[`randsample`](./docs/randsample.md) | Randomly choose a number/percentage of total reads.
[`refinepeak`](./docs/refinepeak.md) | Take raw reads alignment, refine peak summits.
[`callvar`](./docs/callvar.md) | Call variants in given peak regions from the alignment BAM files.
[`hmmratac`](./docs/hmmratac.md) | Dedicated peak calling based on Hidden Markov Model for ATAC-seq data.

For advanced usage, for example, to run `macs3` in a modular way,
please read the [advanced usage](./docs/advanced_usage.md). There is a
[Q&A](./docs/qa.md) document where we collected some common questions
from users.

## Contribute

Please read our [CODE OF CONDUCT](./CODE_OF_CONDUCT.md) and
[How to contribute](./CONTRIBUTING.md) documents. If you have any
questions, suggestion/ideas, or just want to have conversions with
developers and other users in the community, we recommand you use the
[MACS Discussions](https://github.com/macs3-project/MACS/discussions)
instead of posting to our
[Issues](https://github.com/macs3-project/MACS/issues) page.

## Ackowledgement

MACS3 project is sponsored by
[CZI EOSS](https://chanzuckerberg.com/eoss/). And we particularly want
to thank the user community for their supports, feedbacks and
contributions over the years.

## Other useful links

 * [Cistrome](http://cistrome.org/)
 * [bedTools](http://code.google.com/p/bedtools/)
 * [UCSC toolkits](http://hgdownload.cse.ucsc.edu/admin/exe/)
 * [deepTools](https://github.com/deeptools/deepTools/)

