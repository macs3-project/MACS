# Advanced Step-by-step peak calling using MACS3 commands

Over the years, many users have asked us questions about whether they
can analyze their X-Seq (not ChIP-Seq) data using MACS, or customize
the `callpeak` function to suit their specific needs. Typically, we
recommend finding a tool that's more suited to their data type, as
`callpeak` is specifically optimized for ChIP-Seq. However, MACS3 does
offer a range of subcommands that allow you to customize every step of
your analysis. In this tutorial, I will demonstrate how to break down
the main `callpeak` function into a pipeline using various MACS3
subcommands, such as `filterdup`, `predictd`, `pileup`, `bdgcmp`,
`bdgopt`, and `bdgpeakcall` (or `bdgbroadcall` for broad marks). This
approach allows you to adjust or skip steps and modify parameters to
analyze your data in a highly customized way. Additionally, you will
have a complete idea of how `callpeak` works. We'll use two test
files, `CTCF_ChIP_200K.bed.gz` and `CTCF_Control_200K.bed.gz`, which
you can find in the MACS3 GitHub repository in the `test` directory.

*Note, currently this tutorial is mainly for single-end
datasets. Please read the tutorial carefully and modify the command
line for paired-end data accordingly.*

## Step 1: Filter duplicates

In the initial step of ChIP-Seq analysis with `callpeak`, we read both
ChIP and control data and remove redundant reads from each genomic
location. While we won't delve into the reasons behind this process,
we'll explain how you can accomplish this using the `filterdup`
subcommand. By default, `callpeak` permits only one duplicate read per
location, which is set with `--keep-dup=1`. To replicate this setting,
you can do the following:

```
$ macs3 filterdup -i CTCF_ChIP_200K.bed.gz --keep-dup=1 -o CTCF_ChIP_200K_filterdup.bed`
$ macs3 filterdup -i CTCF_Control_200K.bed.gz --keep-dup=1 -o CTCF_Control_200K_filterdup.bed`
```

You can choose a different setting for `--keep-dup` or allow MACS3 to
automatically determine the maximum number of allowed duplicated reads
for each genomic loci for ChIP and control separately. For more
details, check `macs3 filterdup -h`. If you opt for the `--keep-dup
auto` setting, ensure the genome size is set appropriately. MACS3 will
report the final number of reads retained after filtering. It’s
crucial to record these numbers as we need them to normalize the ChIP
and control signals to the same depth. In this example, the numbers
are 199,583 for ChIP and 199,867 for control, resulting in a ratio of
approximately 0.99858. The output files from `filterdup` are in BED
format. Please note that if you are using Paired-end data in BAM
format as input and specify the `-f BAMPE` option, you will get the
output in `BEDPE` format. Check [this document for detail](./BEDPE.md)
on this `BEDPE` format.

## Step 2: Decide the fragment length `d`

This is a crucial step for analyzing ChIP-Seq with MACS3, as well as
other types of data. The location of a sequenced read typically
indicates only the end of a DNA fragment of interest originated from
such as transcription factor binding sites (TFBS) or DNA
hypersensitive sites. To accurately determine the enrichment, you
must estimate the actual length of this DNA fragment. This process can
also be seen as a data smoothing technique. With `macs3 callpeak`, the
output will typically include the identification of a certain number
of peak pairs and the predicted fragment length, denoted as `d` in
MACS3 terminology. This can also be accomplished using the `predictd`
subcommand, which we need to apply only to ChIP data:

```
$ macs3 predictd -i CTCF_ChIP_200K_filterdup.bed -g hs -m 5 50
```

In this process, the `-g` (genome size) needs to be set based on your
sample, and the mfold parameters for `-m` must be reasonably
established. To mimic the default behavior of macs3 `callpeak`, use
`-m 5 50`. Of course, you can adjust these parameters as needed. The
output from `predictd` will provide the fragment length `d`, which in
this example is `254`. Make sure to write this number down, as we'll
need it in the next step. However, if you prefer not to extend the
reads or if you have a more accurate estimate of the fragment length,
you can skip this step. 

Please note that for paired-end data, this step is still necessary
since this function will tell us the average fragment length defined
by the mapping locations of read pairs. For example, if you run this
on the `CTCF_PE_ChIP_chr22_50k.bedpe.gz` file in the test directory:

```
$ macs3 predictd -f BEDPE -i CTCF_PE_ChIP_chr22_50k.bedpe.gz
```

This function will print out that `# Average insertion length of all
pairs is 253 bps`. You should write down this number `253` for
building the control bias track for pair-end data in Step 4.

## Step 3: Extend ChIP sample to get ChIP coverage track

Now that you've estimated the fragment length, we can proceed to
generate a pileup track for the ChIP sample using the MACS3 `pileup`
subcommand. The following command replicates the behavior of callpeak:

`
$ macs3 pileup -i CTCF_ChIP_200K_filterdup.bed -o CTCF_ChIP_200K_filterdup.pileup.bdg --extsize 254
`

This command produces a file in BEDGRAPH format,
`CTCF_ChIP_200K_filterdup.pileup.bdg`, which contains the fragment
pileup signals for the ChIP sample. We use `--extsize 254`, where
`254` is the number obtained from the `predictd` step. In ChIP-Seq
data processing, as demonstrated in this tutorial, we extend reads in
the 5' -> 3' direction, which is the default behavior of the `pileup`
function.

If you are analyzing DNAse-Seq data, or if the cutting site detected
by short read sequencing (the 5' mapping location of the read) is
believed to be centrally located within the DNA region of interest,
you should use the `-B` option. This setting allows for the extension
of the 5' location in both directions. For example, `-B 100` would
extend the 5' cutting site on each side by 100 base pairs (200bps in
total) before generating the pileup signal.

For paired-end data, there is no need to specify `--extsize` as
`pileup` will automatically pile up the entire region between the two
5' mapping locations of the read pair. However, it is crucial to use
`-f BEDPE` or `-f BAMPE` to indicate to MACS3 that the input file
contains paired-end data.

## Step 4: Build local bias track from control

By default, the MACS3 `callpeak` function calculates local bias by
considering the maximum bias from the surrounding 1kb (set by
`--slocal`), 10kb (set by `--llocal`), the fragment length `d`
(predicted from the `predictd` function), and the whole genome
background. In this section, we demonstrate how each bias is
calculated and how they can be combined using the subcommands. 

Please note that, for paired-data data, we should treat the control as
single-end data by NOT specifying `-f BAMPE` or `-f BEDPE`. Also, we
should use the average fragment length from Step 2 as `d` in the
following steps.

### The `d` background

To create the background noise track, extend the control read to both
sides using the `-B` option in the `pileup` function. This approach
assumes that the cutting site from the control sample represent
certain noise from the exact same location. To accomplish this, use
half of `d` obtained from the `predictd` module. For instance, with
`d` equaling 254, use 127 (half of 254) as follows:

```
$ macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 127 -o d_bg.bdg
```

The BEDGRAPH file `d_bg.bdg` contains the background noise data (`d`
background) from the control sample.

### The slocal background

Next, you can create a background noise track for the `slocal`bps
local window, or a 1kb window by default. Imagine that each cutting
site represents a surrounding noise of 1kb (this value is
adjustable). Use the following command:

```
$ macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 500 -o 1k_bg.bdg
```

Here, 500 represents half of 1k, as the default `slocal` for
`callpeak`. Since the ChIP signal track was constructed by extending
reads to `d` size fragments, we need to normalize the 1kb noise by
multiplying the values by `d/slocal`, which is 254/1000 = 0.254 in our
example. To do this, apply the `bdgopt` subcommand:

```
$ macs3 bdgopt -i 1k_bg.bdg -m multiply -p 0.254 -o 1k_bg_norm.bdg
```

The file `1k_bg_norm.bdg` contains the normalized slocal background
from the control sample. Note, normalization for the `d` background is
not necessary since the multiplier is simply 1.

### The llocal background

The background noise from a larger region can be generated similarly
to the previous approach for the slocal background, with the only
difference being the extension size. By default, MACS3 `callpeak` uses
a 10kb surrounding window, but this value can be adjusted:

```
$ macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 5000 -o 10k_bg.bdg
```

Here, the `extsize` should be set as half of the llocal, which is 5000
in this case. The appropriate multiplier is `d/llocal`, or 0.0254 in
our example:

```
$ macs3 bdgopt -i 10k_bg.bdg -m multiply -p 0.0254 -o 10k_bg_norm.bdg
```

The file `10k_bg_norm.bdg` now contains the normalized llocal
background from the control.

### The genome background

The whole genome background is calculated using the formula:
`number_of_control_reads * fragment_length / genome_size`. In our
example, this calculation would be:

```
199867 * 254 / 2700000000 ≈ 0.0188023
```

You don't need to execute any subcommands to create a genome
background track, as it's represented by just a single value.

### Combine and generate the maximum background noise

To compute the maximum bias for each genomic location, you can follow
the default behavior of MACS3 `callpeak` or customize your pipeline to
include additional noise backgrounds (such as 5k or 50k). Here's how
to combine the background noises and determine the maximum bias:

First, take the maximum between the slocal (1k) and llocal (10k)
backgrounds:

```
macs3 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
```

Next, compare this maximum with the `d` background:

```
macs3 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
```

Finally, combine this with the genome-wide background using the
`bdgopt` subcommand:

```
macs3 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p .0188023 -o local_bias_raw.bdg
```

The resulting file `local_bias_raw.bdg` is a BEDGRAPH file containing
the maximum local bias from control data.

## Step 5: Scale the ChIP and control to the same sequencing depth

To ensure accurate comparison between ChIP and control signals, both
must be scaled to the same sequencing depth. The `callpeak` module
typically scales down the larger sample to match the smaller one,
preventing inflation of smaller values and enhancing the specificity
of the results. In our example, after duplicate filtering, the final
read counts are 199,583 for ChIP and 199,867 for control. Therefore,
the control bias needs to be scaled down using the ratio between ChIP
and control, which is 199,583/199,867 = 0.99858. To accomplish this:


```
$ macs3 bdgopt -i local_bias_raw.bdg -m multiply -p .99858 -o local_lambda.bdg
```

The output file is named `local_lambda.bdg`, as it contains values
that represent the lambda (or expected value), which can be compared
with ChIP signals using the local Poisson test.

## Step 6: Compare ChIP and local lambda to get the scores in pvalue or qvalue

To identify enriched regions and predict peaks, the ChIP signals and
local lambda stored in the BEDGRAPH file must be compared using a
statistical model. This is done using the `bdgcmp` module, which
outputs a score for each base pair in the genome. Despite potentially
large data sizes, the BEDGRAPH format efficiently manages this by
merging nearby regions with identical scores. Theoretically, the size
of the output file for scores depends on the complexity of your
data. The maximum number of data points, when using `d`, `slocal`, and
`llocal` backgrounds, is the minimum between the genome size and
approximately `(number_of_ChIP_reads + number_of_control_reads*3)*2`,
which in our case is about 1.6 million. 

The commands to generate score tracks are:

```
$ macs3 bdgcmp -t CTCF_ChIP_200K_filterdup.pileup.bdg -c local_lambda.bdg -m qpois -o CTCF_ChIP_200K_qvalue.bdg
```

or

```
$ macs3 bdgcmp -t CTCF_ChIP_200K_filterdup.pileup.bdg -c local_lambda.bdg -m ppois -o CTCF_ChIP_200K_pvalue.bdg
```

The `CTCF_ChIP_200K_pvalue.bdg` or `CTCF_ChIP_200K_qvalue.bdg` file
contains the `-log10(p-values)` or `-log10(q-values)` for each base
pair, derived through a local Poisson test. This means the ChIP signal
at each base pair will be compared against the corresponding local
lambda from the control using a Poisson model. *Note*, following this
tutorial ensures that there are no zeros in the local lambda track, as
the smallest value is the whole genome background. However, if the
genome background is not included, many zeros in the local lambda can
disrupt the Poisson test. In such cases, you need to set the
`pseudocount` for `bdgcmp` using the `-p` option. This pseudocount
will be added to both ChIP and local lambda values before the
test. The choice of pseudocount is largely arbitrary, and you may find
various discussions online. Generally, a higher pseudocount increases
specificity but decreases sensitivity.

## Step 7: Call peaks on score track using a cutoff

The final step in peak calling is to identify regions that surpass a
specific score cutoff using the `bdgpeakcall` function for narrow peak
calling. Although this process may seem straightforward, it involves
two additional parameters:

1. **Merging Nearby Regions**: If two regions exceed the cutoff and
   are separated by a smaller, lower-scoring region, they should be
   merged into a single larger region to account for
   fluctuations. This merging threshold is set as the read length in
   the MACS3 `callpeak` function, reflecting the dataset's
   resolution. For `bdgpeakcall`, you must obtain the read length from
   Step 1 or by examining the raw fastq file and set this with the
   `-g` option.

2. **Minimum Peak Length**: To avoid calling excessively small peaks,
   a minimum peak length is required. The MACS3 `callpeak` function
   automatically uses the fragment size `d` for this purpose. In
   `bdgpeakcall`, set the `-l` option to the `d` value determined in
   Step 2.

Finally, set the cutoff value. Remember, the scores from `bdgcmp` are
in -log10 format. For example, to set a cutoff of 0.05, use a -log10
value of approximately 1.3. The command is as follows:


```
$ macs3 bdgpeakcall -i CTCF_ChIP_200K_qvalue.bdg -c 1.301 -l 245 -g 100 -o CTCF_ChIP_200K_peaks.bed
```

The output is essentially a narrowPeak format file (a type of BED
file), which includes the locations of peaks with the summit location
noted in the last column. If you want to explore how to better decide
the cutoff, please read the [tutorial for cutoff analysis](./cutoffanalysis.md)

