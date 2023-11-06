# Advanced Step-by-step peak calling using MACS3 commands

Over the years, I have got many emails from users asking if they can
analyze their X-Seq (not ChIP-Seq) data using MACS, or if they can
turn on or off some features in `callpeak` for their special needs. In
most of cases, I would simply reply that they may have to find more
dedicated tool for the type of your data, because the `callpeak`
module is specifically designed and tuned for ChIP-Seq data. However,
MACS3 in fact contains a suite of subcommands and if you can design a
pipeline to combine them, you can control every single step and
analyze your data in a highly customized way. In this tutorial, I show
how the MACS3 main function `callpeak` can be decomposed into a
pipeline containing MACS3 subcommands,
including `filterdup`, `predictd`, `pileup`, `bdgcmp`, `bdgopt`,
and `bdgpeakcall` (or `bdgbroadcall` in case of broad mark).  To
analyze your special data in a special way, you may need to skip some
of the steps or tweak some of the parameters of certain steps. Now
let\'s suppose we are dealing with the two testing
files `CTCF_ChIP_200K.bed.gz` and `CTCF_Control_200K.bed.gz`, that you
can find in MACS3 github repository. 

*Note, currently this tutorial is for single-end datasets. Please
modify the instructions for paired-end data by yourself.*

## Step 1: Filter duplicates

In the first step of ChIP-Seq analysis by `callpeak`, ChIP and control
data need to be read and the redundant reads at each genomic loci have
to be removed. I won\'t go over the rationale, but just tell you how
this can be done by `filterdup` subcommand. By default, the maximum
number of allowed duplicated reads is 1, or `--keep-dup=1` for
`callpeak`. To simulate this behavior, do the following:

`$ macs3 filterdup -i CTCF_ChIP_200K.bed.gz --keep-dup=1 -o CTCF_ChIP_200K_filterdup.bed`
`$ macs3 filterdup -i CTCF_Control_200K.bed.gz --keep-dup=1 -o CTCF_Control_200K_filterdup.bed`

You can set different number for `--keep-dup` or let MACS3
automatically decide the maximum allowed duplicated reads for each
genomic loci for ChIP and control separately. Check `macs3 filterdup
-h` for detail, and remember if you go with auto way, the genome size
should be set accordingly. *Note*, in the output, MACS3 will give
you the final number of reads kept after filtering, you\'d better
write the numbers down since we need them when we have to scale the
ChIP and control signals to the same depth. In this case, the number
is 199583 for ChIP and 199867 for control, and the ratio between them
is 199583/199867=.99858

## Step 2: Decide the fragment length `d`

This is an important step for MACS3 to analyze ChIP-Seq and also for
other types of data since the location of sequenced read may only tell
you the end of a DNA fragment that you are interested in (such as TFBS
or DNA hypersensitive regions), and you have to estimate how long this
DNA fragment is in order to recover the actual enrichment. You can also
regard this as a data smoothing technic. You know that `macs3 callpeak`
will output something like, it can identify certain number of pairs of
peaks and it can predict the fragment length, or `d` in MACS3
terminology, using cross-correlation. In fact, you can also do this
using `predictd` module. Normally, we only need to do this for ChIP
data:

`
$ macs3 predictd -i CTCF_ChIP_200K_filterdup.bed -g hs -m 5 50
`

Here the `-g` (the genome size) need to be set according to your sample,
and the mfold parameters have to be set reasonably. To simulate the
default behavior of `macs3 callpeak`, set `-m 5 50`. Of course, you can
tweak it. The output from `predictd` will tell you the fragment length
`d`, and in this example, it is *254*. Write the number down on your
notebook since we will need it in the next step. Of course, if you do
not want to extend the reads or you have a better estimation on fragment
length, you can simply skip this step.

## Step 3: Extend ChIP sample to get ChIP coverage track

Now you have estimated the fragment length, next, we can use MACS3
`pileup` subcommand to generate a pileup track in BEDGRAPH format for
ChIP sample. Since we are dealing with ChIP-Seq data in this tutorial,
we need to extend reads in 5\' to 3\' direction which is the default
behavior of `pileup` function. If you are dealing with some DNAse-Seq
data or you think the cutting site, that is detected by short read
sequencing, is just in the `middle` of the fragment you are interested
in, you need to use `-B` option to extend the read in both direction.
Here is the command to simulate `callpeak` behavior:

`
$ macs3 pileup -i CTCF_ChIP_200K_filterdup.bed -o CTCF_ChIP_200K_filterdup.pileup.bdg --extsize 254
`

The file `CTCF_ChIP_200K_filterdup.pileup.bdg` now contains the
fragment pileup signals for ChIP sample.

## Step 4: Build local bias track from control

By default, MACS3 `callpeak` function computes the local bias by taking
the maximum bias from surrounding 1kb (set by `--slocal`), 10kb (set by
`--llocal`), the size of fragment length `d` (predicted as what you got
from `predictd`), and the whole genome background. Here I show you how
each of the bias is calculated and how they can be combined using the
subcommands.

### The `d` background

Basically, to create the background noise track, you need to extend the
control read to both sides (-B option) using `pileup` function. The idea
is that the cutting site from control sample contains the noise
representing a region surrounding it. To do this, take half of `d` you
got from `predictd`, 127 (1/2 of 254) for our example, then:

`
$ macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 127 -o d_bg.bdg
`

The file `d_bg.bdg` contains the `d` background from control.

### The slocal background

Next, you can create a background noise track of slocal local window, or
1kb window by default. Simply imagine that each cutting site (sequenced
read) represent a 1kb (default, you can tweak it) surrounding noise. So:

`
$ macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 500 -o 1k_bg.bdg
`

Note, here 500 is the 1/2 of 1k. Because the ChIP signal track was built
by extending reads into `d` size fragments, we have to normalize the 1kb
noise by multiplying the values by `d/slocal`, which is 254/1000=0.254
in our example. To do so, use the `bdgopt` subcommand:

`
$ macs3 bdgopt -i 1k_bg.bdg -m multiply -p 0.254 -o 1k_bg_norm.bdg
`

The file`1k_bg_norm.bdg` contains the slocal background from control.
Note, we don\'t have to do this for `d` background because the
multiplier is simply 1.

### The llocal background

The background noise from larger region can be generated in the same way
as slocal backgound. The only difference is the size for extension.
MACS3 `callpeak` by default asks for 10kb (you can change this value)
surrounding window, so:

`
$ macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 5000 -o 10k_bg.bdg
`

The extsize has to be set as 1/2 of llocal. Then, the multiplier now is
`d/llocal`, or 0.0254 in our example.

`
$ macs3 bdgopt -i 10k_bg.bdg -m multiply -p 0.0254 -o 10k_bg_norm.bdg
`

The file `10k_bg_norm.bdg` now contains the slocal background from
control.

### The genome background

The whole genome background can be calculated as
`the_number_of_control_reads\fragment_length/genome_size`, and in our
example, it is $`199867*254/2700000000 ~= .0188023`$. You don\'t need to
run subcommands to build a genome background track since it\'s just a
single value.

### Combine and generate the maximum background noise

Now all the above background noises have to be combined and the maximum
bias for each genomic location need be computed. This is the default
behavior of MACS3 `callpeak`, but you can have your own pipeline to
include some of them or even make more noise (such as 5k or 50k
background) then include more tracks. Here is the way to combine and get
the maximum bias.

Take the maximum between slocal (1k) and llocal (10k) background:

`
macs3 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
`

Then, take the maximum then by comparing with `d` background:

`
macs3 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
`

Finally, combine with the genome wide background using `bdgopt`
subcommand

`
macs3 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p .0188023 -o local_bias_raw.bdg
`

Now the file `local_bias_raw.bdg` is a BEDGRAPH file containing the
raw local bias from control data.

## Step 5: Scale the ChIP and control to the same sequencing depth

In order to compare ChIP and control signals, the ChIP pileup and
control lambda have to be scaled to the same sequencing depth. The
default behavior in `callpeak` module is to scale down the larger sample
to the smaller one. This will make sure the noise won\'t be amplified
through scaling and improve the specificity of the final results. In our
example, the final number of reads for ChIP and control are 199583 and
199867 after filtering duplicates, so the control bias have to be scaled
down by multiplying with the ratio between ChIP and control which is
199583/199867=.99858. To do so:

`
$ macs3 bdgopt -i local_bias_raw.bdg -m multiply -p .99858 -o local_lambda.bdg
`

Now, I name the output file as `local_lambda.bdg` since the values in
the file can be regarded as the lambda (or expected value) and can be
compared with ChIP signals using the local Poisson test.

## Step 6: Compare ChIP and local lambda to get the scores in pvalue or qvalue

Next, to find enriched regions and predict the so-called \'peaks\',
the ChIP signals and local lambda stored in BEDGRAPH file have to be
compared using certain statistic model. To do so, you need to use
`bdgcmp` module, which will output score for each basepair in the
genome. You may wonder it may need a huge file to save score for each
basepair in the genome, however the format BEDGRAPH can deal with the
problem by merging nearby regions with the same score. So
theoratically, the size of the output file for scores depends on the
complexity of your data, and the maximum number of data points, if
`d`, `slocal`, and `llocal` background are all used, is the minimum
value of the genome size and approximately
`(#read_ChIP+#reads_control\*3)\*2`, in our case about 1.6 million.
The command to generate score tracks is

`
$ macs3 bdgcmp -t CTCF_ChIP_200K_filterdup.pileup.bdg -c local_lambda.bdg -m qpois -o CTCF_ChIP_200K_qvalue.bdg
`
or

`
$ macs3 bdgcmp -t CTCF_ChIP_200K_filterdup.pileup.bdg -c local_bias.bdg -m ppois -o CTCF_ChIP_200K_pvalue.bdg
`

The `CTCF_ChIP_200K_pvalue.bdg` or `CTCF_ChIP_200K_qvalue.bdg` file
contains the -log10(p-value)s or -log10(q-value)s for each basepair
through local Poisson test, which means the ChIP signal at each basepair
will be tested against the corresponding local lambda derived from
control with Poisson model. *Note*, if you follow this tutorial, then
you won\'t get any `0` in the local lambda track because the smallest
value is the whole genome background; however if you do not include
genome background, you will see many `0`s in local lambda which will
crash the Poisson test. In this case, you need to set the
`pseudocount` for `bdgcmp` through `-p` option. The pseudocount will
be added to both ChIP and local lambda values before Poisson test. The
choice of pseudocount is mainly arbitrary and you can search on the web
to see some discussion. But in general, higher the pseudocount, higher
the specificity and lower the sensitivity.

## Step 7: Call peaks on score track using a cutoff

The final task of peak calling is to just take the scores and call those
regions higher than certain cutoff. We can use the `bdgpeakcall`
function for narrow peak calling and `bdgrroadcall` for broad peak
calling, and I will only cover the usage of `bdgpeakcall` in this
tutorial. It looks simple however there are extra two parameters for the
task. First, if two nearby regions are both above cutoff but the region
in-between is lower, and if the region in-between is small enough, we
should merge the two nearby regions together into a bigger one and
tolerate the fluctuation. This value is set as the read length in MACS3
`callpeak` function since the read length represent the resolution of
the dataset. As for `bdgpeakcall`, you need to get the read length
information from Step 1 or by checking the raw fastq file and set the
`-g` option. Secondly, we don\'t want to call too many small `peaks`
so we need to have a minimum length for the peak. MACS3 `callpeak`
function automatically uses the fragment size `d` as the parameter of
minimum length of peak, and as for `bdgpeakcall`, you need to set the
`-l` option as the `d` you got from Step 2. Last, you have to set the
cutoff value. Remember the scores in the output from `bdgcmp` are in
-log10 form, so if you need the cutoff as 0.05, the -log10 value is
about 1.3. The final command is like:

`
$ macs3 bdgpeakcall -i CTCF_ChIP_200K_qvalue.bdg -c 1.301 -l 245 -g 100 -o CTCF_ChIP_200K_peaks.bed
`

The output is in fact a narrowPeak format file (a type of BED file)
which contains locations of peaks and the summit location in the last
column.

 
