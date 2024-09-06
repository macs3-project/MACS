# MACS: Model-based Analysis for ChIP-Seq

![Status](https://img.shields.io/pypi/status/macs3.svg) ![License](https://img.shields.io/github/license/macs3-project/MACS) ![Programming languages](https://img.shields.io/github/languages/top/macs3-project/MACS) ![CI x64](https://github.com/macs3-project/MACS/workflows/MACS3%20CI%20x64/badge.svg?branch=master) ![CI non x64](https://github.com/macs3-project/MACS/workflows/MACS3%20CI%20non%20x64/badge.svg?branch=master) ![CI Mac OS](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-macos.yml/badge.svg?branch=master) [![CZI's Essential Open Source Software for Science](https://chanzuckerberg.github.io/open-science/badges/CZI-EOSS.svg)](https://czi.co/EOSS)


[![PyPI
download](https://img.shields.io/pypi/dm/macs3?label=pypi%20downloads)](https://pypistats.org/packages/macs3)

Latest Release:
* Github: [![Github Release](https://img.shields.io/github/v/release/macs3-project/MACS)](https://github.com/macs3-project/MACS/releases)
* PyPI: [![PyPI Release](https://img.shields.io/pypi/v/macs3.svg)](https://pypi.org/project/MACS3/)
* Bioconda:[![Bioconda Badge](https://anaconda.org/bioconda/macs3/badges/version.svg)](https://anaconda.org/bioconda/macs3)
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

Please find MACS3 documentations through [MACS3
website](https://macs3-project.github.io/MACS/).

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

