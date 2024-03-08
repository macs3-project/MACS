# MACS: Model-based Analysis for ChIP-Seq

![Status](https://img.shields.io/pypi/status/macs3.svg) ![License](https://img.shields.io/github/license/macs3-project/MACS) ![Programming languages](https://img.shields.io/github/languages/top/macs3-project/MACS) ![CI x64](https://github.com/macs3-project/MACS/workflows/MACS3%20CI%20x64/badge.svg) ![CI non x64](https://github.com/macs3-project/MACS/workflows/MACS3%20CI%20non%20x64/badge.svg) ![CI MacOS](https://github.com/macs3-project/MACS/workflows/MACS3%20CI%20Mac%20OS/badge.svg)

[![PyPI
download](https://img.shields.io/pypi/dm/macs3?label=pypi%20downloads)](https://pypistats.org/packages/macs3)

Latest Release:
* Github: [![Github Release](https://img.shields.io/github/v/release/macs3-project/MACS)](https://github.com/macs3-project/MACS/releases) 
* PyPI: [![PyPI Release](https://img.shields.io/pypi/v/macs3.svg)![PyPI Python Version](https://img.shields.io/pypi/pyversions/MACS3)![PyPI Format](https://img.shields.io/pypi/format/macs3)](https://pypi.org/project/macs3/) 
* Anaconda: [![Anaconda-Server Badge](https://anaconda.org/macs3/macs3/badges/version.svg)](https://anaconda.org/macs3/macs3) 
* Debian Med: [![Debian Stable](https://img.shields.io/debian/v/macs/stable?label=debian%20stable)](https://packages.debian.org/stable/macs)[![Debian Unstable](https://img.shields.io/debian/v/macs/sid?label=debian%20sid)](https://packages.debian.org/sid/macs)

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

## Changes for MACS (3.0.1) 

*Bugs fixed*

1) Fixed a bug that the `hmmatac` can't correctly save the digested
	signal
	files. [#605](https://github.com/macs3-project/MACS/issues/605)
	[#611](https://github.com/macs3-project/MACS/pull/611)

2) Applied a patch to remove cython requirement from the installed
	system. (it's needed for building the
	package). [#606](https://github.com/macs3-project/MACS/issues/606)
	[#612](https://github.com/macs3-project/MACS/pull/612)

3) Relax the testing script while comparing the peaks called from
	current codes and the standard peaks. To implement this, we added
	'intersection' function to 'Regions' class to find the
	intersecting regions of two Regions object (similar to PeakIO but
	only recording chromosome, start and end positions). And we
	updated the unit test 'test_Region.py' then implemented a script
	'jaccard.py' to compute the Jaccard Index of two peak files. If
	the JI > 0.99 we would think the peaks called and the standard
	peaks are similar. This is to avoid the problem caused by
	different Numpy/SciPy/sci-kit learn libraries, when certain peak
	coordinates may have 10bps
	difference. [#615](https://github.com/macs3-project/MACS/issues/615)
	[#619](https://github.com/macs3-project/MACS/pull/619)
	
4) Due to [the changes in scikit-learn
	1.3.0](https://scikit-learn.org/1.3/whats_new/v1.3.html), the way
	hmmlearn 0.3 uses Kmeans will end up with inconsistent results
	between sklearn <1.3 and sklearn >=1.3. Therefore, we patched the
	class hmm.GaussianHMM and adjusted the standard output from
	`hmmratac` subcommand. The change is based on [hmmlearn
	PR#545](https://github.com/hmmlearn/hmmlearn/pull/545). The idea
	is to do the random seeding of KMeans 10 times. Now the `hmmratac`
	results should be more consistent (at least
	JI>0.99). [#615](https://github.com/macs3-project/MACS/issues/615)
	[#620](https://github.com/macs3-project/MACS/pull/620)

*Other*
	
1) We added some dependencies to MACS3. `hmmratc` subcommand needs 
	`hmmlearn` library, `hmmlearn` needs `scikit-learn` and 
	`scikit-learn` needs `scipy`. Since major releases have happened 
	for both`scipy` and `scikit-learn`, we have to set specific 
	version requirements for them in order to make sure the output 
	results from `hmmratac` are consistent. 

2) We updated our documentation website using 
    Sphinx. https://macs3-project.github.io/MACS/

## Changes for MACS (3.0.0)

1) Call variants in peak regions directly from BAM files. The
	function was originally developed under code name SAPPER. Now
	SAPPER has been merged into MACS as the `callvar` command. It can
	be used to call SNVs and small INDELs directly from alignment
	files for ChIP-seq or ATAC-seq. We call `fermi-lite` to assemble
	the DNA sequence at the enriched genomic regions (binding sites or
	accessible DNA) and to refine the alignment when necessary. We
	added `simde` as	a submodule in order to support fermi-lite
	library under non-x64 architectures.

2) HMMRATAC module is added as subcommand `hmmratac`. HMMRATAC is a
	dedicated software to analyze ATAC-seq data. The basic idea behind
	HMMRATAC is to digest ATAC-seq data according to the fragment
	length of read pairs into four signal tracks: short fragments,
	mono-nucleosomal fragments, di-nucleosomal fragments and
	tri-nucleosomal fragments. Then integrate the four tracks again
	using Hidden Markov Model to consider three hidden states: open
	region, nucleosomal region, and background region. The orginal
	paper was published in 2019 written in JAVA, by Evan Tarbell. We
	implemented it in Python/Cython and optimize the whole process
	using existing MACS functions and hmmlearn. Now it can run much
	faster than the original JAVA version. Note: evaluation of the
	peak calling results is still underway.

3) Speed/memory optimization.  Use the cykhash to replace python
	dictionary. Use buffer (10MB) to read and parse input file (not
	available for BAM file parser). And many optimization tweaks. We
	added memory monitoring to the runtime messages.

4) R wrappers for MACS -- MACSr for bioconductor.

5) Code cleanup. Reorganize source codes. 

6) Unit testing. 

7) Switch to Github Action for CI, support multi-arch testing
	including x64, armv7, aarch64, s390x and ppc64le. We also test on
	Mac OS 12.

8) MACS tag-shifting model has been refined. Now it will use a naive
	peak calling approach to find ALL possible paired peaks at + and -
	strand, then use all of them to calculate the
	cross-correlation. (a related bug has been fix
	[#442](https://github.com/macs3-project/MACS/issues/442))

9) BAI index and random access to BAM file now is
	supported. [#449](https://github.com/macs3-project/MACS/issues/449).

10) Support of Python > 3.10 [#498](https://github.com/macs3-project/MACS/issues/498)

11) The effective genome size parameters have been updated
	according to deeptools. [#508](https://github.com/macs3-project/MACS/issues/508)

12) Multiple updates regarding dependencies, anaconda built, CI/CD
	process.

13) Cython 3 is supported.

14) Documentations for each subcommand can be found under /docs

*Other*

1) Missing header line while no peaks can be called
[#501](https://github.com/macs3-project/MACS/issues/501)
[#502](https://github.com/macs3-project/MACS/issues/502)

2) Note: different numpy, scipy, sklearn may give slightly
	different results for hmmratac results. The current standard
	results for automated testing in `/test` directory are from Numpy
	1.25.1, Scipy 1.11.1, and sklearn 1.3.0.

## Install

The common way to install MACS is through
[PYPI](https://pypi.org/project/macs3/)) or
[conda](https://anaconda.org/macs3/macs3). Please check the
[INSTALL](docs/INSTALL.md) document for detail.

MACS3 has been tested using GitHub Actions for every push and PR in
the following architectures:

 * x86_64 (Ubuntu 22, Python 3.9, 3.10, 3.11)
 * aarch64 (Ubuntu 22, Python 3.10)
 * armv7 (Ubuntu 22, Python 3.10)
 * ppc64le (Ubuntu 22, Python 3.10)
 * s390x (Ubuntu 22, Python 3.10)
 * Apple chips (Mac OS 13, Python 3.11)

In general, you can install through PyPI as `pip install macs3`.  To
use virtual environment is highly recommended. Or you can install
after unzipping the released package downloaded from Github, then use
`pip install .` command. Please note that, we haven't tested
installation on any Windows OS, so currently only Linux and Mac OS
systems are supported. Also, for aarch64, armv7, ppc64le and s390x,
due to some unknown reason potentially related to the scientific
calculation libraries MACS3 depends on, such as Numpy, Scipy,
hmm-learn, scikit-learn, the results from `hmmratac` subcommand may
not be consistent with the results from x86 or Apple chips. Please be
aware.

## Usage

Example for regular peak calling on TF ChIP-seq:

`macs3 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01`

Example for broad peak calling on Histone Mark ChIP-seq:

`macs3 callpeak -t ChIP.bam -c Control.bam --broad -g hs --broad-cutoff 0.1`

Example for peak calling on ATAC-seq (paired-end mode):

`macs3 callpeak -f BAMPE -t ATAC.bam -g hs -n test -B -q 0.01`

There are currently 14 functions available in MACS3 serving as
sub-commands. Please click on the link to see the detail description
of the subcommands.

Subcommand | Description
-----------|----------
[`callpeak`](docs/callpeak.md) | Main MACS3 Function to call peaks from alignment results.
[`bdgpeakcall`](docs/bdgpeakcall.md) | Call peaks from bedGraph file.
[`bdgbroadcall`](docs/bdgbroadcall.md) | Call nested broad peaks from bedGraph file.
[`bdgcmp`](docs/bdgcmp.md) | Comparing two signal tracks in bedGraph format.
[`bdgopt`](docs/bdgopt.md) | Operate the score column of bedGraph file.
[`cmbreps`](docs/cmbreps.md) | Combine bedGraph files of scores from replicates.
[`bdgdiff`](docs/bdgdiff.md) | Differential peak detection based on paired four bedGraph files.
[`filterdup`](docs/filterdup.md) | Remove duplicate reads, then save in BED/BEDPE format file.
[`predictd`](docs/predictd.md) | Predict d or fragment size from alignment results. In case of PE data, report the average insertion/fragment size from all pairs.
[`pileup`](docs/pileup.md) | Pileup aligned reads (single-end) or fragments (paired-end)
[`randsample`](docs/randsample.md) | Randomly choose a number/percentage of total reads, then save in BED/BEDPE format file.
[`refinepeak`](docs/refinepeak.md) | Take raw reads alignment, refine peak summits.
[`callvar`](docs/callvar.md) | Call variants in given peak regions from the alignment BAM files.
[`hmmratac`](docs/hmmratac.md) | Dedicated peak calling based on Hidden Markov Model for ATAC-seq data.

For advanced usage, for example, to run `macs3` in a modular way,
please read the [advanced usage](docs/Advanced_Step-by-step_Peak_Calling.md). There is a
[Q&A](docs/qa.md) document where we collected some common questions
from users.

## Contribute

Please read our [CODE OF CONDUCT](CODE_OF_CONDUCT.md) and [How to
contribute](CONTRIBUTING.md) documents. If you have any questions,
suggestion/ideas, or just want to have conversions with developers and
other users in the community, we recommend using the [MACS
Discussions](https://github.com/macs3-project/MACS/discussions)
instead of posting to our
[Issues](https://github.com/macs3-project/MACS/issues) page.

## Ackowledgement

MACS3 project is sponsored by
[CZI EOSS](https://chanzuckerberg.com/eoss/). And we particularly want
to thank the user community for their supports, feedbacks and
contributions over the years.

## Citation

2008: [Model-based Analysis of ChIP-Seq
(MACS)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137)

## Other useful links

 * [Cistrome](http://cistrome.org/)
 * [bedTools](http://code.google.com/p/bedtools/)
 * [UCSC toolkits](http://hgdownload.cse.ucsc.edu/admin/exe/)
 * [deepTools](https://github.com/deeptools/deepTools/)


```{toctree}
:maxdepth: 2
:hidden:

docs/INSTALL.md
docs/index.md
docs/Advanced_Step-by-step_Peak_Calling.md
docs/qa.md
docs/tutorial.md
CODE_OF_CONDUCT.md
CONTRIBUTING.md
