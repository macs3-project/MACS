# hmmratac

## Description

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
the whole process using existing MACS functions and hmmlearn, and
integrate HMMRATAC with all other modules in MACS3 such as the support
for single-cell assay (FRAG format support).

Here's an example of how to run the `hmmratac` command:

```
$ macs3 hmmratac -i yeast.bam -n yeast
```

or with the BEDPE format of a much smaller size:

```
$ macs3 hmmratac -i yeast.bedpe.gz -f BEDPE -n yeast
```

You can convert BAMPE to BEDPE by using

```
$ macs3 filterdup --keep-dup all -f BAMPE -i yeast.bam -o yeast.bedpe
```

You can also call accessible regions on the fragment files from
scATAC-seq analysis:

```
$ macs3 hmmratac -i yeast.scATAC.frag.gz -f FRAG --barcodes selected_barcodes.txt --max-count 1
```

Please note that in order to save memory usage and fasten the process,
`hmmratac` will save intermediate temporary file to the disk. The file
size can range from megabytes to gigabytes, depending on how many
candidate regions `hmmratac` needs to decode. The temporary file will
be removed after the job is done. So please make sure there is enough
space in the 'tmp' directory of your system. 

Please use `macs3 hmmratac --help` to see all the options. Here we
list the essential ones.

## Output

The final output file from `hmmratac` is in narrowPeak format
containing the accessible regions (open state in `hmmratac` HMM). The
columns are:

 1. chromosome name
 2. start position of the accessible region
 3. end position of the accesssible region
 4. peak name
 5. peak score. The score is the 10times the maximum foldchange 
    (signal/average signal) within the peak. By default, the 'signal'
    used to calculate foldchange is the total pileup of all types of
    fragments from short to tri-nuc size fragments.
 7. Not used
 8. Not used
 9. Not used
 10. peak summit position. It's the relative position from the start
    position to the peak summit which is defined as the position with
    the maximum foldchange score.

## Essential Options

### `-i INPUT_FILE [INPUT_FILE ...]` / `--input INPUT_FILE [INPUT_FILE ...]`

This is the only REQUIRED parameter for `hmmratac`. Input files
containing the aligment results for ATAC-seq paired end reads. If
multiple files are given as '-t A B C', then they will all be read and
pooled together. The file should be in BAMPE, BEDPE format (aligned in
paired end mode) or FRAG format (from scATAC analysis). Files can be
gzipped. Note: all files should be in the same format. REQUIRED.

### `-f {BAMPE,BEDPE,FRAG}` / `--format {BAMPE,BEDPE,FRAG}`

Format of input files, "BAMPE", "BEDPE", or "FRAG". If there are
multiple files, they should be in the same format -- either BAMPE,
BEDPE, or FRAG. Please note that the BEDPE only contains three columns
-- chromosome, left position of the whole pair, right position of the
whole pair-- and is NOT the same BEDPE format used by BEDTOOLS. To
convert BAMPE to BEDPE, you can use this command `macs3 filterdup
--keep-dup all -f BAMPE -i input.bam -o output.bedpe`. And the FRAG
format is like BEDPE but with two extra columns -- barcode and count
of the fragment. Please note that if you plan to analyze single-cell
ATAC-seq data on a specific list of barcodes, you can only use FRAG
format. For other formats, the barcode information will be
ignored. DEFAULT: "BAMPE".

### `--barcodes` and `--max-count`

These two options are for `-f FRAG` format only. You can let
`hmmratac` work on a selected set of barcodes by specifying
`--barcodes barcode.txt`. The `barcode.txt` should contain list of
selected barcodes and each row represents a specific barcode, like:

```
ATCGATCGATCGATCG
GCTAGCTAGCTAGCTA
...
```

The `--max-count` option is recommended for scATAC-seq since for each
single cell, Tn5 can only cut the same location once for each DNA
molecule. So theoratically, if you are studying a haploid genome, the
maximum count of the same fragment can't be larger than 2 (i.e. should
use `--max-count 1` or `--max-count 2`).

### `--outdir OUTDIR`

If specified all output files will be written to that
directory. Default: the current working directory 

### `-n NAME`/ `--name NAME`

Name for this experiment, which will be used as a prefix to generate
output file names. DEFAULT: "NA" 

### `-e BLACKLIST`/`--blacklist BLACKLIST`

Filename of the file containing the blacklisted regions to exclude
from the process. Any fragments overlapping with blacklisted regions
are excluded. An example of such file can be found from the ENCODE
project at: https://github.com/Boyle-Lab/Blacklist/. Alternatively, if
you wish to exclude centromeres and telomeres, you can find their
genomic coordinates and write them to a BED format file. By default,
there is no blacklist file in use.

### `--modelonly`

This option will only generate the HMM model as a JSON file and
quit. This model can then be applied using the `--model`
option. Default: False

### `--model`

If provided, HMM training will be skipped and a JSON file generated
from a previous HMMRATAC run will be used instead of creating new
one. Default: NA 
   
### `-t HMM_TRAINING_REGIONS` / `--training HMM_TRAINING_REGIONS`

Customized training regions can be provided through this option. `-t`
takes the filename of training regions (previously was BED_file) to
use for training HMM, instead of using foldchange settings to
select. Default: NA 

### `--min-frag-p MIN_FRAG_P`

We will exclude the abnormal fragments that can't be assigned to any
of the four signal tracks. After we use EM to find the means and
stddevs of the four distributions, we will calculate the likelihood
that a given fragment length fit any of the four using normal
distribution. The criteria we will use is that if a fragment length
has less than MIN_FRAG_P probability to be like either of short, mono,
di, or tri-nuc fragment, we will exclude it while generating the four
signal tracks for later HMM training and prediction. The value should
be between 0 and 1. Larger the value, more abnormal fragments will be
allowed. So if you want to include more 'ideal' fragments, make this
value smaller. Default = 0.001

### `--cutoff-analysis-only`

Only run the cutoff analysis and output a report. After generating the
report, the whole process will stop. By default, the cutoff analysis
will be included in the whole process, but won't quit after the report
is generated. The report will help user decide the three crucial
parameters for `-l`, `-u`, and `-c`. So it's highly recommanded to run
this first!  Please read the report and instructions in [Choices of
cutoff values](#choices-of-cutoff-values) on how to decide the three
crucial parameters. The resolution of cutoff analysis can be
controlled by `--cutoff-analysis-max` and `--cutoff-analysis-steps`
options.

### `--cutoff-analysis-max`

The maximum cutoff score for performing cutoff analysis. Together with
`--cutoff-analysis-steps`, the resolution in the final report can be
controlled. Please check the description in `--cutoff-analysis-steps`
for detail. The default value is 100.

### `--cutoff-analysis-steps`

Steps for performing cutoff analysis. It will be used to decide which
cutoff value should be included in the final report. Larger the value,
higher resolution the cutoff analysis can be. The cutoff analysis
function will first find the smallest (at least 0) and the largest (at
most 100, and controlled by --cutoff-analysis-max) foldchange score in
the data, then break the range of foldchange score into
`CUTOFF_ANALYSIS_STEPS` intervals. It will then use each foldchange
score as cutoff to call peaks and calculate the total number of
candidate peaks, the total basepairs of peaks, and the average length
of peak in basepair. Please note that the final report ideally should
include `CUTOFF_ANALYSIS_STEPS` rows, but in practice, if the
foldchange cutoff yield zero peak, the row for that foldchange value
won't be included.  The default is 100.

### `--hmm-type`

We provide two types of emissions for the Hidden Markov Model -- the
Gaussian model and the Poisson model. By default, the Gaussian
emission will be used (as `--hmm-type gaussian`). To choose Poisson
emission, use `--hmm-type poisson`. The Gaussian emission can be
described by mean and variance for each state, while the simpler
Poisson only needs the lambda value. The difference can be found in
the saved json file for HMM.

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
candidate regions. The higher the prescan cutoff, the fewer regions
will be considered. Must be > 1. This is an important parameter for
decoding so please read. The purpose of this parameter is to EXCLUDE
those chromatin regions having noises/random enrichment so we can have
a large number of possible regions to predict the HMM states. It's
highly recommended to run the `--cutoff-analysis-only` first to decide
the lower cutoff `-l`, the upper cutoff `-u`, and the pre-scanning
cutoff `-c`. The pre-scanning cutoff should be the cutoff close to the
BOTTOM of the cutoff analysis result that can capture a large number
of possible peaks with normal length (average length 500-1000bps). In
most cases, please do not pick a cutoff too low that captures almost
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

