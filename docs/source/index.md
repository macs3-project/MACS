# MACS: Model-based Analysis for ChIP-Seq

![Status](https://img.shields.io/pypi/status/macs3.svg) ![License](https://img.shields.io/github/license/macs3-project/MACS) ![Programming languages](https://img.shields.io/github/languages/top/macs3-project/MACS) [![CI x64](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-x64.yml/badge.svg)](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-x64.yml) [![CI non x64](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-non-x64.yml/badge.svg)](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-non-x64.yml) [![CI Mac OS](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-macos.yml/badge.svg)](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-macos.yml) [![PyPI download](https://img.shields.io/pypi/dm/macs3?label=pypi%20downloads)](https://pypistats.org/packages/macs3)

Latest Release:
* Github: [![Github Release](https://img.shields.io/github/v/release/macs3-project/MACS)](https://github.com/macs3-project/MACS/releases)
* PyPI: [![PyPI Release](https://img.shields.io/pypi/v/macs3.svg)](https://pypi.org/project/MACS3/)
* Bioconda:[![Bioconda Badge](https://anaconda.org/bioconda/macs3/badges/version.svg)](https://anaconda.org/bioconda/macs3)
* Debian Med: [![Debian Stable](https://img.shields.io/debian/v/macs/stable?label=debian%20stable)](https://packages.debian.org/stable/macs)[![Debian Unstable](https://img.shields.io/debian/v/macs/sid?label=debian%20sid)](https://packages.debian.org/sid/macs3)

## Introduction

With the advancement of sequencing technologies, Chromatin
Immunoprecipitation followed by high-throughput sequencing (ChIP-Seq)
has become a popular method for studying genome-wide protein-DNA
interactions. With the purpose of addressing the need for a robust
ChIP-Seq analysis tool, we introduce **M**odel-based **A**nalysis of
**C**hIP-**S**eq (MACS), a powerful tool for identifying transcription
factor binding sites. MACS accounts for the complexity of the genome
to assess the significance of enriched ChIP regions and enhances the
spatial resolution of binding sites by integrating both sequencing tag
position and orientation. MACS can be readily applied to ChIP-Seq data
alone, or in conjunction with a control sample, thus enhancing
specificity. Furthermore, as a versatile peak-caller, MACS can be
employed in any "DNA enrichment assay" to answer the fundamental
question: *Where are the regions with significant read coverage
compared to random background?*

## Changes for MACS 3.0.5

### Features added

1. Added `PETrackII.return_anndata` to build a sparse barcode-by-peak
   AnnData count matrix directly from single-cell fragment data. The
   method preserves fragment counts and merges overlapping or adjacent
   input regions.

2. Added `hmmratac --jump` to control the update factor used by the
   fragment-length EM algorithm. The default is 0.5.

### Performance improvements

1. Reworked and vectorized the NumPy-backed `PileupV2` implementation.
   `callpeak`, `pileup`, peak-model construction, and the single-end
   and paired-end track classes now use the optimized routines by
   default. Now it has 1.5x speedup over 3.0.4.

2. Improved the performance of `PETrackII.exclude` and sparse
   barcode-by-peak matrix construction.

3. Replaced `cykhash` caches with Python dictionaries and removed the
   `cykhash` dependency.

### Bugs fixed

1. Fixed `bdgdiff` region scores being truncated to integers before
   calculating their length-weighted mean. Decimal scores are now
   preserved, including with cutoffs below 1 (#715).

2. Fixed `callpeak --call-summits` assigning the first above-cutoff
   chunk's score to a smoothed maximum in a below-cutoff gap. Invalid
   gap maxima are now discarded (#741).

3. Fixed an `IndexError` in bedGraph peak refinement when a peak has no
   overlapping bedGraph content, as encountered in `hmmratac` (#735).

4. `hmmratac --cutoff-analysis-only` now exits with status 0 after
   successfully writing its report (#704).

### Compatibility changes

1. MACS3 now requires Python 3.12 or later and declares support for
   Python 3.12, 3.13, and 3.14.

2. `pandas` and `anndata` are now runtime dependencies for the AnnData
   export API.

### Documentation

1. Added automated major-version and release-to-release performance
   benchmarking and their documentation.

2. Corrected the cutoff-analysis documentation for `callpeak`,
   `bdgpeakcall`, and `hmmratac`, and updated the recommended MACS3
   citation.
	
## Install

The common way to install MACS is through
[PYPI](https://pypi.org/project/macs3/)) or
[conda](https://anaconda.org/macs3/macs3). Please check the
[INSTALL](docs/INSTALL.md) document for detail.

MACS3 has been tested using GitHub Actions for every push and PR in
the following architectures:

 * x86_64 (Ubuntu 22, Python 3.9, 3.10, 3.11, 3.12, 3.13)
 * aarch64 (Ubuntu 22, Python 3.10)
 * armv7 (Ubuntu 22, Python 3.10)
 * ppc64le (Ubuntu 22, Python 3.10)
 * s390x (Ubuntu 22, Python 3.10)
 * Apple chips (Mac OS 13, Python 3.9, 3.10, 3.11, 3.12, 3.13)

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

Example for peak calling on ATAC-seq with HMMATAC:

`macs3 hmmratac -i ATAC.bam -f BAMPE -n test`

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
[`hmmratac`](docs/hmmratac.md) | Dedicated peak calling based on Hidden Markov Model for ATAC-seq or scATAC-seq data.

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

## Support MACS3

I maintain MACS3 in my spare time. If you find the project useful and
would like to support its continued development, you can
[![Buy Me a Coffee](https://img.shields.io/badge/Buy_Me_a_Coffee-support-FFDD00?style=for-the-badge&logo=buy-me-a-coffee&logoColor=000000)](https://buymeacoffee.com/taoliu). Your contribution
will help cover my ever-growing consumption of coffee and tokens.

## Ackowledgement

MACS3 project is sponsored by [![CZI's Essential Open Source Software for Science](https://chanzuckerberg.github.io/open-science/badges/CZI-EOSS.svg)](https://czi.co/EOSS) through EOSS2 (2020-2022) and EOSS4 (2021-2025). And we particularly want to thank the user community for their supports, feedbacks and contributions over the years.

## Citation

For MACS version 2 and 3, please cite our 2026 paper [MACS3: A Peak-calling Platform for Bulk and Single-cell Regulatory Genomics](https://academic.oup.com/gpb/advance-article/doi/10.1093/gpbjnl/qzag097/8802115)

If you are using MACS version 1, please cite our 2008 paper [Model-based Analysis of ChIP-Seq
(MACS)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137)

## Other useful links

 * [Cistrome](http://cistrome.org/)
 * [bedTools](http://code.google.com/p/bedtools/)
 * [UCSC toolkits](http://hgdownload.cse.ucsc.edu/admin/exe/)
 * [deepTools](https://github.com/deeptools/deepTools/)


```{toctree}
:maxdepth: 3
:hidden:

docs/INSTALL.md
docs/version_benchmark.md
docs/subcommands_index.md
docs/fileformats_index.md
docs/tutorial.md
docs/qa.md
CODE_OF_CONDUCT.md
CONTRIBUTING.md
docs/api/index
