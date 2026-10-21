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

## Changes for MACS 3.1.0

See the [full 3.1.0 release notes](docs/release_notes/v3.1.0.md).

### Features and performance

- **Major `callpeak` speedup:** Jacob Schreiber's auto-optimize work in
  [PR #756](https://github.com/macs3-project/MACS/pull/756) measured 41–64×
  faster results across four default workflows on an 8-core host. Speedups
  vary by workload, and large inputs can use more peak memory.
- Added `--threads` to supported commands for a per-run CPU budget. Set
  `MACS3_NUM_THREADS` to choose a default, or pass an execution policy from
  the Python API. See [CPU budgets and process backends](docs/execution.md).
- Added `--process-backend {auto,serial}` to `callpeak` and `callvar`.

### Fixes

- Automatic duplicate filtering now uses the control sample's own threshold
  for control reads.
- Corrected lower-tail Poisson probabilities for very large expected counts.
- FRAG input may end without a trailing newline. Truncated BAM alignment
  records are rejected as invalid input.
- Improved reliability and cleanup when peak-calling workers fail or multiple
  operations run concurrently.

### Documentation

- Clarified the Poisson tail convention used for reported p-values.

### Acknowledgements

Thanks to Jacob Schreiber ([PR #756](https://github.com/macs3-project/MACS/pull/756)),
Cindy Krafft ([PR #753](https://github.com/macs3-project/MACS/pull/753),
[PR #777](https://github.com/macs3-project/MACS/pull/777), and
[PR #789](https://github.com/macs3-project/MACS/pull/789)), Michael R. Crusoe
([PR #754](https://github.com/macs3-project/MACS/pull/754)), and Apoorva Garg
([PR #793](https://github.com/macs3-project/MACS/pull/793)) for their
contributions to this release.
	
## Install

Install MACS3 from [PyPI](https://pypi.org/project/macs3/) or, when the
desired version is available, [Bioconda](https://bioconda.github.io/recipes/macs3/README.html).
MACS3 3.1.0 requires Python 3.12 or later. See the
[installation guide](docs/INSTALL.md) for supported environments, build
requirements, and source-install instructions.

GitHub Actions tests each code change on x86_64 Linux and macOS with
Python 3.12, 3.13, and 3.14, and on aarch64 Ubuntu 24.04 with Python
3.12. The armv7, ppc64le, and s390x Ubuntu 24.04 tests use Python
3.12 and run through a separate manually triggered workflow.

We recommend a virtual environment. Install from PyPI with
`python -m pip install macs3`, or install from an unpacked source distribution
with `python -m pip install .`. Windows is not currently supported. On
aarch64, armv7, ppc64le, and s390x, `hmmratac` results may differ from
results on x86 or Apple Silicon; validate cross-architecture comparisons.

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

Please read our {doc}`Code of Conduct </CODE_OF_CONDUCT>` and
{doc}`Contributing guide </CONTRIBUTING>`. For questions, suggestions, or
ideas, we recommend using the [MACS
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
docs/release_notes/v3.1.0.md
docs/execution.md
docs/subcommands_index.md
docs/fileformats_index.md
docs/tutorial.md
docs/qa.md
CODE_OF_CONDUCT.md
CONTRIBUTING.md
docs/api/index
