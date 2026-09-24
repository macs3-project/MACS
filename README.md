# MACS: Model-based Analysis for ChIP-Seq

![Status](https://img.shields.io/pypi/status/macs3.svg) ![License](https://img.shields.io/github/license/macs3-project/MACS) ![Programming languages](https://img.shields.io/github/languages/top/macs3-project/MACS) [![CI x64](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-x64.yml/badge.svg)](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-x64.yml) [![CI non x64](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-non-x64.yml/badge.svg)](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-non-x64.yml) [![CI Mac OS](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-macos.yml/badge.svg)](https://github.com/macs3-project/MACS/actions/workflows/build-and-test-MACS3-macos.yml)

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

## What's new in MACS 3.0.5

- Added `PETrackII.return_anndata` for efficient creation of sparse
  barcode-by-peak AnnData matrices from single-cell fragment data.
- Switched peak calling, pileup output, and peak-model construction to
  faster NumPy-backed `PileupV2` routines, and optimized fragment
  exclusion and score caches. Now it has 1.5x speedup over 3.0.4.
- Added `hmmratac --jump` to control fragment-length EM updates.
- Fixed decimal-score truncation in `bdgdiff` (#715), incorrect summit
  scores for maxima in below-cutoff gaps (#741), and an `IndexError`
  during `hmmratac` peak refinement (#735).
- Corrected summit coordinates and peak boundaries around the temporary
  search padding (#747), and fixed shape filtering of summit candidates
  on both sides of a local maximum (#748).
- Made successful `hmmratac --cutoff-analysis-only` runs exit with
  status 0 (#704).
- MACS3 now requires Python 3.12 or later. The `cykhash` dependency was
  removed, while `pandas` and `anndata` were added for AnnData export.

See the [ChangeLog](ChangeLog) for the complete release notes.

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

## Note

Now the default branch of MACS3 has been renamed to 'main'. If you are still using old 'master' branch in your cloned Git repository for MACS3, please use the following command to rename it:

```
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```
