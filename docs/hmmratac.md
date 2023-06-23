# Hidden Markov Model peak caller for ATAC-seq -- HMMRATAC

HMMRATAC (`macs3 hmmratac`) is a dedicated peak calling algorithm
based on Hidden Markov Model for ATAC-seq data. The basic idea behind
HMMRATAC is to digest ATAC-seq data according to the fragment length
of read pairs into four signal tracks: short fragments,
mono-nucleosomal fragments, di-nucleosomal fragments and
tri-nucleosomal fragments. Then integrate the four tracks using Hidden
Markov Model to consider three hidden states: open region, nucleosomal
region, and background region. The [orginal
paper](https://academic.oup.com/nar/article/47/16/e91/5519166) was
published in 2019, and the original software was written in JAVA, by
the then PhD student Evan Tarbell, a mohawk bioinformatician. In MACS3
project, we implemented HMMRATAC idea in Python/Cython and optimize
the whole process using existing MACS functions and hmmlearn.

Here's an example of runing the `hmmratac` command:

```
$ macs3 hmmratac -b test.bam -n test_result
```

Please use `macs3 hmmratac --help` to see all the options. Here we
list the essential ones.

## Essential Options

### `-b BAM_FILE [BAM_FILE ...]` / `--bam BAM_FILE [BAM_FILE ...]`

This is the only REQUIRED parameter for `hmmratac`. The file can
should be in BAMPE format. If multiple files are given as '-b A B C',
then they will all be read and pooled together. REQUIRED.

### `--outdir OUTDIR`

If specified all output files will be written to that directory. Default: the current working directory

### `-n NAME`/ `--name NAME`
Name for this experiment, which will be used as a prefix to generate output file names. DEFAULT: "NA"

### `--cutoff-analysis-only`

 Only run the cutoff analysis and output a report. After generating
 the report, the process will stop. The report will help user decide
 the three crucial parameters for `-l`, `-u`, and `-c`. So it's highly
 recommanded to run this first! Please read the report and
 instructions in `Choices of cutoff values` on how to decide the three
 crucial parameters

### `-u HMM_UPPER` / `--upper HMM_UPPER`

Upper limit on fold change range for choosing training sites. This is
an important parameter for training so please read. The purpose of
this parameter is to EXCLUDE those unusually highly enriched chromatin
regions so we can get training samples in 'ordinary' regions
instead. It's highly recommended to run the `--cutoff-analysis-only`
first to decide the lower cutoff `-l`, the upper cutoff `-u`, and the
pre-scanning cutoff `-c`. The upper cutoff should be the cutoff in the
cutoff analysis result that can capture some (typically hundreds of)
extremely high enrichment and unusually wide peaks. Default: 20

### `-l HMM_LOWER` / `--lower HMM_LOWER`
Lower limit on fold change range for choosing training sites. This is
an important parameter for training so please read. The purpose of
this parameter is to ONLY INCLUDE those chromatin regions having
ordinary enrichment so we can get training samples to learn the common
features through HMM. It's highly recommended to run the
`--cutoff-analysis-only` first to decide the lower cutoff `-l`, the
upper cutoff `-u`, and the pre-scanning cutoff `-c`. The lower cutoff
should be the cutoff in the cutoff analysis result that can capture
moderate number ( about 10k ) of peaks with normal width ( average
length 500-1000bps long). Default: 10

### `-c PRESCAN_CUTOFF` / `--prescan-cutoff PRESCAN_CUTOFF`

The fold change cutoff for prescanning candidate regions in the whole
dataset. Then we will use HMM to predict/decode states on these
candidate regions. Higher the prescan cutoff, fewer regions will be
considered. Must > 1. This is an important parameter for decoding so
please read. The purpose of this parameter is to EXCLUDE those
chromatin regions having noises/random enrichment so we can have a
large number of possible regions to predict the HMM states. It's
highly recommended to run the `--cutoff-analysis-only` first to decide
the lower cutoff `-l`, the upper cutoff `-u`, and the pre-scanning
cutoff `-c`. The pre-scanning cutoff should be the cutoff close to the
BOTTOM of the cutoff analysis result that can capture large number of
possible peaks with normal length (average length 500-1000bps). In
most cases, please do not pick a cutoff too low that capture almost
all the background noises from the data. Default: 1.2


## Choices of cutoff values

Before you proceed, it's highly recommended to run with
`--cutoff-analysis-only` for the initial attempt. When this option is
activated, `hmmratac` will use EM to estimate the best parameters for
fragment sizes of short fragments, mono-, di-, and tri-nucleosomes,
pileup fragments, convert the pileup values into fold-change, and
analyze each possible cutoff. This analysis includes the number of
peaks that can be called, their average peak length, and their total
length. After the report is generated, you can review its contents and
decide on the optimal `-l`, `-u`, and `-c`.

The report consists of four columns:

1. Score: the possible fold change cutoff value.
2. npeaks: the number of peaks.
3. lpeaks: the total length of all peaks.
4. avelpeak: the average length of peaks.

While there's no universal rule, here are a few suggestions:

- The lower cutoff should be the cutoff in the report that captures a
  moderate number (about 10k) of peaks with a normal width (average
  length 500-1000bps long).
- The upper cutoff should capture some (typically hundreds of)
  extremely high enrichment and unusually wide peaks in the
  report. The aim here is to exclude abnormal enrichment caused by
  artifacts such as repetitive regions.
- The pre-scanning cutoff should be the cutoff close to the BOTTOM of
  the report that can capture a large number of potential peaks with a
  normal length (average length 500-1000bps). However, it's
  recommended not to use the lowest cutoff value in the report as this
  may include too much noise from the genome.
  
## Tune the HMM model

It's highly recommended to check the runtime message of the HMM model
after training. An example is like this:

```
#4 Train Hidden Markov Model with Multivariate Gaussian Emission
#  Extract signals in training regions with bin size of 10
#  Use Baum-Welch algorithm to train the HMM
#   HMM converged: True
#  Write HMM parameters into JSON: test_model.json
#  The Hidden Markov Model for signals of binsize of 10 basepairs:
#   open state index: state2
#   nucleosomal state index: state1
#   background state index: state0
#   Starting probabilities of states:
#                            bg        nuc       open
#                        0.7994     0.1312    0.06942
#   HMM Transition probabilities:
#                            bg        nuc       open
#               bg->     0.9842    0.01202   0.003759
#              nuc->    0.03093     0.9562    0.01287
#             open->   0.007891    0.01038     0.9817
#   HMM Emissions (mean):
#                         short       mono         di        tri
#               bg:      0.2551      1.526     0.4646    0.07071
#              nuc:       6.538      17.94      3.422    0.05819
#             open:       5.016      17.47      6.897      2.121
```

We will 'guess' which hidden state is for the open region, which is
the nucleosomal region, and which is the background. We compute from
the HMM Emission matrix to pick the state with the highest sum of mean
signals as the open state, the lowest as the backgound state, then the
rest is the nucleosomal state. However it may not work in every
case. In the above example, it may be tricky to call the second row as
'nuc' and the third as 'open'. If the users want to exchange the state
assignments of the 'nuc' and 'open', they can modify the state
assignment in the HMM model file (e.g. test_model.json). For the above
example, the model.json looks like this (we skipped some detail):

```
{"startprob": [...], "transmat": [...], "means": [...], "covars": [...], 
"covariance_type": "full", "n_features": 4, 
"i_open_region": 2, "i_background_region": 0, "i_nucleosomal_region": 1,
"hmm_binsize": 10}
```

We can modify the assignment of: `"i_open_region": 2,
"i_background_region": 0, "i_nucleosomal_region": 1,` by assigning `1`
to open, and `2` to nucleosomal as: `"i_open_region": 1,
"i_background_region": 0, "i_nucleosomal_region": 2,` Then save the
HMM in a new model file such as `new_model.json`.

Then next, we can re-run `macs3 hmmratac` with the same parameters
plus an extra option for the HMM model file like `macs3 hmmratac
--model new_model.json`

