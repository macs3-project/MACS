# Contributing to MACS project

First of all, thanks for taking the time to contribute! :tada:

The following is a set of guidelines for contributing to MACS, which is hosted in the [MACS project](https://github.com/macs3-project) on GitHub. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

[What should I know before I get started?](#what-should-i-know-before-i-get-started)
  * [MACS documentations](#macs-documentations)
  * [MACS Slack Community](#macs-slack-community)

[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Your First Code Contribution](#your-first-code-contribution)
  * [Pull Requests](#pull-requests)

[Styleguides](#styleguides)
  <!-- * [Git Commit Messages](#git-commit-messages) this doesnt link anywhere -->
  * [Python Styleguide](#python-styleguide)
  * [Documentation Styleguide](#documentation-styleguide)

[Additional Notes](#additional-notes)
  * [Issue and Pull Request Labels](#issue-and-pull-request-labels)

## Code of Conduct

This project and everyone participating in it is governed by the [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [macs3-project@google.com](mailto:macs3-project@google.com).

## What Should I Know Before I Get Started

:see_no_evil: **Note:** We do not prohibit issues regarding general questions (e.g. usage of MACS or suggestions on data analysis). However, please label such issues with the tag  [General Question][search-macs-label-general-question]. But, before you do so, you'll get faster results by using the resources below.

### MACS Documentations

We have an official tutorial and a detailed FAQ compiled from the community. 

* [MACS documentations](https://macs3-project.github.io/MACS/)
  * [MACS tutorial](https://macs3.github.io) <!-- change link? -->
  * [MACS FAQ](https://macs3.github.io) <!-- change link? -->
  * [MACS API](https://macs3.github.io) <!-- change link? -->

### MACS Slack Community

If chat is more your speed, you can join the [MACS3 Slack community](https://macs3.slack.com/):
* Even though Slack is a chat service, sometimes it takes several hours for community members to respond &mdash; please be patient!
* Use the `#general` channel for general questions or discussion about MACS3
* Use the `#data-analysis` channel for detailed questions about data-analysis. For example, how to analyze your own data using MACS.
* Use the `#feature-request` channel for suggestions on new features.
* Use the `#bug-reports` channel for reporting and discussing bugs you find.
* There are many other channels available, check the channel list

## How Can I Contribute?

### Reporting Bugs

This section guides you through submitting a bug report for MACS. Following these guidelines helps maintainers and the community understand your report :pencil:, reproduce the behavior :computer: :computer:, and find related reports :mag_right:.

Before creating bug reports, please check [this list](#before-submitting-a-bug-report) as you might find out that you don't need to create one. When you are creating a bug report, please [include as many details as possible](#how-do-i-submit-a-good-bug-report). Fill out [the required template](https://github.com/macs3-project/macs/.github/blob/master/.github/ISSUE_TEMPLATE/bug_report.md), the information it asks for helps us resolve issues faster.

**Note:** If you find a **Closed** issue that seems like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one.

#### Before Submitting A Bug Report

* Search [previous bugs](https://github.com/macs3-project/MACS/issues) that other people have submitted. Your issue may have been addressed, someone else may have had a similar experience, or a solution may have been implemented. Consider searching for issues that are closed or are tagged with the label "Bug Report".
* Similarly, search through the [discussion page](https://github.com/macs3-project/MACS/discussions). There may be a discussion already started about your issue.

#### How Do I Submit A (Good) Bug Report?

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/). 

Explain the problem and include additional details to help maintainers reproduce the problem:

...

Provide more context by answering these questions:

...

Provide details about your configuration and environment:

...

### Suggesting Enhancements

...

#### Before Submitting An Enhancement Suggestion

...

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/). 
...

### Your First Code Contribution

#### Local development

MACS can be developed locally. 

* [Setting up local testing environment through Docker](./docs/testing_in_docker.md)

### Pull Requests

The process described here has several goals:

- Maintain MACS quality
- Fix problems that are important to users
- Engage the community in working toward the best possible MACS
- Enable a sustainable system for MACS maintainers to review contributions

Please follow these steps to have your contribution considered by the maintainers:

1. Follow all instructions in [the template](PULL_REQUEST_TEMPLATE.md)
2. Follow the [styleguides](#styleguides)
3. After you submit your pull request, verify that all [status checks](https://help.github.com/articles/about-status-checks/) are passing <details><summary>What if the status checks are failing?</summary>If a status check is failing, and you believe that the failure is unrelated to your change, please leave a comment on the pull request explaining why you believe the failure is unrelated. A maintainer will re-run the status check for you. If we conclude that the failure was a false positive, then we will open an issue to track that problem with our status check suite.</details>

While the prerequisites above must be satisfied prior to having your pull request reviewed, the reviewer(s) may ask you to complete additional design work, tests, or other changes before your pull request can be ultimately accepted.

## Styleguides

### Python Styleguide

All Python/Cython codes must adhere to PEP.


### Documentation Styleguide

* Use [Markdown](https://daringfireball.net/projects/markdown).

## Additional Notes

### Issue and Pull Request Labels

This section lists the labels we use to help us track and manage issues and pull requests. 

[GitHub search](https://help.github.com/articles/searching-issues/) makes it easy to use labels for finding groups of issues or pull requests you're interested in. We  encourage you to read about [other search filters](https://help.github.com/articles/searching-issues/) which will help you write more focused queries.

The labels are loosely grouped by their purpose, but it's not required that every issue have a label from every group or that an issue can't have more than one label from the same group.

Please open an issue on `macs3-project/MACS` if you have suggestions for new labels, and if you notice some labels are missing on some repositories, then please open an issue on that repository.

#### Type of Issue and Issue State

#### Topic Categories


#### Pull Request Labels

