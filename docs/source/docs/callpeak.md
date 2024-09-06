# callpeak

## Overview
This is the main function in MACS3. It will take alignment files in
various format (please check the detail below) and call the
significantly enriched regions in the genome as 'peaks'.  It can be
invoked by `macs3 callpeak` . If you type this command with `-h`, you
will see a full description of command-line options. Here we only list
the essentials.

## Essential Commandline Options

### Input and Output

- `-t`/`--treatment`

  This is the only REQUIRED parameter for MACS3. The file can be in
  any supported format -- see detail in the `--format` option. If you
  have more than one alignment file, you can specify them as `-t A B
  C`.  MACS3 will pool up all these files together.

- `-c`/`--control`

  The control, genomic input or mock IP data file. Please follow the
  same direction as for `-t`/`--treatment`.

- `-n`/`--name`

  The name string of the experiment. MACS3 will use this string NAME
  to create output files like `NAME_peaks.xls`,
  `NAME_negative_peaks.xls`, `NAME_peaks.bed` , `NAME_summits.bed`,
  `NAME_model.r` and so on. So please avoid any confliction between
  these filenames and your existing files.

- `-f`/`--format FORMAT`

  Format of tag file can be `ELAND`, `BED`, `ELANDMULTI`,
  `ELANDEXPORT`, `SAM`, `BAM`, `BOWTIE`, `BAMPE`, or `BEDPE`. Default
  is `AUTO` which will allow MACS3 to decide the format
  automatically. `AUTO` is also useful when you combine different
  formats of files. Note that MACS3 can't detect `BAMPE` or `BEDPE`
  format with `AUTO`, and you have to implicitly specify the format
  for `BAMPE` and `BEDPE`.

  Nowadays, the most common formats are `BED` or `BAM` (including
  `BEDPE` and `BAMPE`). Our recommendation is to convert your data to
  `BED` or `BAM` first.

  Also, MACS3 can detect and read gzipped file. For example, `.bed.gz`
  file can be directly used without being uncompressed with `--format
  BED`.

  Here are detailed explanation of the recommended formats:

  - `BED`

    The BED format can be found at [UCSC genome browser
    website](http://genome.ucsc.edu/FAQ/FAQformat#format1).

    The essential columns in BED format input are the 1st column
    `chromosome name`, the 2nd `start position`, the 3rd `end
    position`, and the 6th, `strand`.

    Note that, for `BED` format, the 6th column of strand information
    is required by MACS3. And please pay attention that the
    coordinates in BED format are zero-based and half-open. See more
    detail at [UCSC
    site](http://genome.ucsc.edu/FAQ/FAQtracks#tracks1).

  - `BAM`/`SAM`

    If the format is `BAM`/`SAM`, please check the definition in
    [samtools](https://samtools.github.io/hts-specs/SAMv1.pdf).  If
    the `BAM` file is generated for paired-end data, MACS3 will only
    keep the left mate(5' end) tag. However, when format `BAMPE` is
    specified, MACS3 will use the real fragments inferred from
    alignment results for reads pileup.

  - `BEDPE` or `BAMPE`

    A special mode will be triggered while the format is specified as
    `BAMPE` or `BEDPE`. In this way, MACS3 will process the `BAM` or
    `BED` files as paired-end data. Instead of building a bimodal
    distribution of plus and minus strand reads to predict fragment
    size, MACS3 will use actual insert sizes of pairs of reads to
    build fragment pileup.

    The `BAMPE` format is just a `BAM` format containing paired-end
    alignment information, such as those from `BWA` or `BOWTIE`.

    The `BEDPE` format is a simplified and more flexible `BED` format,
    which only contains the first three columns defining the
    chromosome name, left and right position of the fragment from
    Paired-end sequencing. Please note, this is NOT the same format
    used by `bedtools`, and the `bedtools` version of `BEDPE` is
    actually not in a standard `BED` format. You can use MACS3
    subcommand [`randsample`](./randsample.md) or
    [`filterdup`](./filterdup.md) to convert a `BAMPE` file containing
    paired-end information to a `BEDPE` format file:
    
    ```
    macs3 randsample -i the_BAMPE_file.bam -f BAMPE -p 100 -o the_BEDPE_file.bed
    ```
	or
	
	```
    macs3 filterdup -i the_BAMPE_file.bam -f BAMPE --keep-dup all -o the_BEDPE_file.bed
    ```
	  

- `--outdir`

  MACS3 will save all output files into the specified folder for this 
  option. A new folder will be created if necessary. 

- `-B`/`--bdg`

  If this flag is on, MACS3 will store the fragment pileup, control
  lambda in bedGraph files. The bedGraph files will be stored in the
  current directory named `NAME_treat_pileup.bdg` for treatment data,
  `NAME_control_lambda.bdg` for local lambda values from control.

- `--trackline`

  MACS3 will include the trackline in the header of output files,
  including the bedGraph, narrowPeak, gappedPeak, BED format files. To
  include this trackline in the header is necessary while uploading
  them to the UCSC genome browser. You can also mannually add these
  trackline to corresponding output files. For example, in order to
  upload narrowPeak file to UCSC browser, add this to as the first
  line -- `track type=narrowPeak name=`"my_peaks`" description=\"my
  peaks\"`. Default: Not to include any trackline.

### Options controling peak calling behaviors

- `-g`/`--gsize`

  It's the mappable genome size or effective genome size which is
  defined as the genome size which can be sequenced. Because of the
  repetitive features on the chromosomes, the actual mappable genome
  size will be smaller than the original size, about 90% or 70% of the
  genome size. The default *hs* ~2.9e9 is recommended for human
  genome. Here are all precompiled parameters for effective genome
  size from
  [deeptools](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html):

  * hs: 2,913,022,398 for GRCh38
  * mm: 2,652,783,500 for GRCm38
  * ce: 100,286,401 for WBcel235
  * dm: 142,573,017 for dm6

  Please check deeptools webpage to find the appropriate effective
  genome size if you want a more accurate estimation regarding
  specific assembly and read length.

  Users may want to use k-mer tools to simulate mapping of Xbps long
  reads to target genome, and to find the ideal effective genome
  size. However, usually by taking away the simple repeats and Ns from
  the total genome, one can get an approximate number of effective
  genome size. A slight difference in the number won't cause a big
  difference of peak calls, because this number is used to estimate a
  genome-wide noise level which is usually the least significant one
  compared with the *local biases* modeled by MACS3.

- `-s`/`--tsize`

  The size of sequencing tags. If you don't specify it, MACS3 will try
  to use the first 10 sequences from your input treatment file to
  determine the tag size. Specifying it will override the
  automatically determined tag size.

- `-q`/`--qvalue`

  The q-value (minimum FDR) cutoff to call significant
  regions. Default is 0.05. For broad marks, you can try 0.01 as the
  cutoff. The q-values are calculated from p-values using the
  [Benjamini-Hochberg
  procedure](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure).

- `-p`/`--pvalue`

  The p-value cutoff. If `-p` is specified, MACS3 will use p-value
  instead of q-value.

- `--min-length`, `--max-gap`

  These two options can be used to fine-tune the peak calling behavior
  by specifying the minimum length of a called peak and the maximum
  allowed a gap between two nearby regions to be merged. In other
  words, a called peak has to be longer than `min-length`, and if the
  distance between two nearby peaks is smaller than `max-gap` then
  they will be merged as one. If they are not set, MACS3 will set the
  DEFAULT value for `min-length` as the predicted fragment size `d`,
  and the DEFAULT value for `max-gap` as the detected read
  length. Note, if you set a `min-length` value smaller than the
  fragment size, it may have NO effect on the result. For broad peak
  calling with `--broad` option set, the DEFAULT `max-gap` for merging
  nearby stronger peaks will be the same as narrow peak calling, and 4
  times of the `max-gap` will be used to merge nearby weaker (broad)
  peaks. You can also use `--cutoff-analysis` option with the default
  setting, and check the column `avelpeak` under different cutoff
  values to decide a reasonable `min-length` value.

- `--nolambda`

  With this flag on, MACS3 will use the background lambda as local
  lambda. This means MACS3 will not consider the local bias at peak
  candidate regions. It is particularly recommended while calling
  peaks without control sample.

- `--slocal`, `--llocal`

  These two parameters control which two levels of regions will be
  checked around the peak regions to calculate the maximum lambda as
  local lambda. By default, MACS3 considers 1000bp for small local
  region(`--slocal`), and 10000bps for large local region(`--llocal`)
  which captures the bias from a long-range effect like an open
  chromatin domain. You can tweak these according to your
  project. Remember that if the region is set too small, a sharp spike
  in the input data may kill a significant peak.

- `--nomodel`

  While on, MACS3 will bypass building the shifting model. Please
  combine the usage of `--extsize` and `--shift` to achieve the effect
  you expect.

- `--extsize`

  While `--nomodel` is set, MACS3 uses this parameter to extend reads
  in 5'->3' direction to fix-sized fragments. For example, if the size
  of the binding region for your transcription factor is 200 bp, and
  you want to bypass the model building by MACS3, this parameter can
  be set as 200. This option is only valid when `--nomodel` is set or
  when MACS3 fails to build model and `--fix-bimodal` is on.

- `--shift`

  Note, this is NOT the legacy `--shiftsize` option which is replaced
  by `--extsize` from MACS version 2! You can set an arbitrary shift
  in bp here to adjust the alignment positions of reads in the whole
  library. Please use discretion while setting it other than the
  default value (0). When `--nomodel` is set, MACS3 will use this
  value to move cutting ends (5') then apply `--extsize` from 5' to 3'
  direction to extend them to fragments. When this value is negative,
  the cutting ends (5') will be moved toward 3'->5' direction,
  otherwise 5'->3' direction. Recommended to keep it as default 0 for
  ChIP-Seq datasets, or -1 * half of *EXTSIZE* together with
  `--extsize` option for detecting enriched cutting loci such as
  certain DNAseI-Seq datasets. Note, you can't set values other than 0
  if the format is BAMPE or BEDPE for paired-end data. The default is
  0.

  Here are some examples for combining `--shift` and `--extsize`:

  1. To find enriched cutting sites such as some DNAse-Seq
    datasets. In this case, all 5' ends of sequenced reads should be
    extended in both directions to smooth the pileup signals. If the
    wanted smoothing window is 200bps, then use `--nomodel --shift
    -100 --extsize 200`.

  2. For certain nucleosome-seq data, we need to pile up the centers
    of nucleosomes using a half-nucleosome size for wavelet analysis
    (e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about
    147bps, this option can be used: `--nomodel --shift 37 --extsize
    73`.

- `--keep-dup`

  It controls the MACS3 behavior towards duplicate tags at the exact
  same location -- the same coordination and the same strand. You can
  set this as `auto`, `all`, or an integer value. The `auto` option
  makes MACS3 calculate the maximum tags at the exact same location
  based on binomial distribution using 1e-5 as p-value cutoff; and the
  `all` option keeps every tag.  If an integer is given, at most this
  number of tags will be kept at the same location. The default is to
  keep one tag at the same location. Default: 1

- `--broad`
  
  This option, along with the `bdgbroadcall` command, facilitates
  broad peak calling, producing results in the UCSC gappedPeak format
  which encapsulates a nested structure of peaks. To conceptualize
  'nested' peaks, picture a gene structure housing regions analogous
  to exons (strong peaks) and introns coupled with UTRs (weak
  peaks). The broad peak calling process utilizes two distinct cutoffs
  to discern broader, weaker peaks (`--broad-cutoff`) and narrower,
  stronger peaks (`-p` or `-q`), which are subsequently nested to
  provide a detailed peak landscape. Please note that, the `max-gap`
  value for merging nearby weaker/broad peaks is 4 times of `max-gap`
  for merging nearby stronger peaks. The later one can be controlled
  by `--max-gap` option, and by default it is the average
  fragment/insertion length in the PE data. DEFAULT: False

  Please note that, if you only want to call 'broader' peak and not
  interested in the nested peak structure, please simply use `-p` or
  `-q` with weaker cutoff instead of using `--broad` option.

- `--broad-cutoff`

  Cutoff for the broad region. This option is not available unless
  `--broad` is set. Please note that if `-p` is set, this is a p-value
  cutoff, otherwise, it's a q-value cutoff.  DEFAULT: 0.1

- `--scale-to <large|small>`

  When set to `large`, linearly scale the smaller dataset to the same
  depth as the larger dataset. By default or being set as `small`, the
  larger dataset will be scaled towards the smaller dataset. Beware,
  to scale up small data would cause more false positives. So the
  default behavior `small` is recommended.

- `--call-summits`

  MACS3 will now reanalyze the shape of signal profile (p or q-score
  depending on the cutoff setting) to deconvolve subpeaks within each
  peak called from the general procedure. It's highly recommended to
  detect adjacent binding events. While used, the output subpeaks of a
  big peak region will have the same peak boundaries, and different
  scores and peak summit positions.

### Other options

- `--buffer-size`

  MACS3 uses a buffer size for incrementally increasing internal array
  size to store reads alignment information for each chromosome or
  contig. To increase the buffer size, MACS3 can run faster but will
  waste more memory if certain chromosome/contig only has very few
  reads. In most cases, the default value 100000 works fine. However,
  if there are a large number of chromosomes/contigs in your alignment
  and reads per chromosome/contigs are few, it's recommended to
  specify a smaller buffer size in order to decrease memory usage (but
  it will take longer time to read alignment files). Minimum memory
  requested for reading an alignment file is about # of CHROMOSOME *
  BUFFER_SIZE * 8 Bytes. DEFAULT: 100000

- `--cutoff-analysis`

  While set, MACS3 will analyze the number or total length of peaks
  that can be called by different cutoff then output a summary table
  to help the user decide a better cutoff. Note, minlen and maxgap may
  affect the results. DEFAULT: False
  
  Different with the option in `bdgpeakcall`, `callpeak` will perform
  both tasks to call peaks and to generate a report for cutoff
  analysis. Please check the section *Cutoff Analysis* for more
  detail.

## Output files

1. `NAME_peaks.xls` is a tabular file which contains information about
   called peaks. You can open it in Excel and sort/filter using excel
   functions. Information include:
   
   - chromosome name
   - start position of peak
   - end position of peak
   - length of peak region
   - absolute peak summit position
   - pileup height at peak summit
   - -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then
      this value should be 10)
   - fold enrichment for this peak summit against random Poisson
      distribution with local lambda,
   - -log10(qvalue) at peak summit
   
   Coordinates in XLS is 1-based which is different from BED
   format. When `--broad` is enabled for broad peak calling, the
   pileup, p-value, q-value, and fold change in the XLS file will be
   the mean value across the entire peak region, since peak summit
   won't be called in broad peak calling mode.

2. `NAME_peaks.narrowPeak` is BED6+4 format file which contains the
   peak locations together with peak summit, p-value, and q-value. If you
   plan to load it to the UCSC genome browser, please make sure that
   you turn on `--trackline` option. Definition of some specific
   columns are:
   
   - 5th: integer score for display. It's calculated as
     `int(-10*log10pvalue)` or `int(-10*log10qvalue)` depending on
     whether `-p` (pvalue) or `-q` (qvalue) is used as score
     cutoff. Please note that currently this value might be out of the
     [0-1000] range defined in [UCSC ENCODE narrowPeak
     format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12). You
     can let the value saturated at 1000 (i.e. p/q-value = 10^-100) by
     using the following 1-liner awk: `awk -v OFS="\t"
     '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak`
   - 7th: fold-change at peak summit
   - 8th: -log10pvalue at peak summit
   - 9th: -log10qvalue at peak summit
   - 10th: relative summit position to peak start
   
   Remove the beginning track line if you want to analyze it by other
   tools.

3. `NAME_summits.bed` is in BED format, which contains the peak
   summits locations for every peak. The 5th column in this file is
   the same as what is in the `narrowPeak` file. If you want to find
   the motifs at the binding sites, this file is recommended. The file
   can be loaded directly to the UCSC genome browser with
   `--trackline` option on. Remove the beginning track line if you
   want to analyze it by other tools.

4. `NAME_peaks.broadPeak` is in BED6+3 format which is similar to
   `narrowPeak` file, except for missing the 10th column for
   annotating peak summits. This file and the `gappedPeak` file will
   only be available when `--broad` is enabled. Since in the broad
   peak calling mode, the peak summit won't be called, the values in
   the 5th, and 7-9th columns are the mean value across all positions
   in the peak region. Refer to `narrowPeak` if you want to fix the
   value issue in the 5th column. The file can be loaded directly to
   the UCSC genome browser with `--trackline` option on.

5. `NAME_peaks.gappedPeak` is in BED12+3 format which contains both
   the broad region and narrow peaks. The 5th column is the score for
   showing grey levels on the UCSC browser as in `narrowPeak`. The 7th
   is the start of the first narrow peak in the region, and the 8th
   column is the end. The 9th column should be RGB color key, however,
   we keep 0 here to use the default color, so change it if you
   want. The 10th column tells how many blocks including the starting
   1bp and ending 1bp of broad regions. The 11th column shows the
   length of each block and 12th for the start of each block. 13th:
   fold-change, 14th: *-log10pvalue*, 15th: *-log10qvalue*. The file
   can be loaded directly to the UCSC genome browser. Refer to
   `narrowPeak` if you want to fix the value issue in the 5th
   column. The file can be loaded directly to the UCSC genome browser
   with `--trackline` option on.

6. `NAME_model.r` is an R script which you can use to produce a PDF
   image of the model based on your data. Load it to R by:

   `$ Rscript NAME_model.r`

   Then a pdf file `NAME_model.pdf` will be generated in your current
   directory. Note, R is required to draw this figure.

7. The `NAME_treat_pileup.bdg` and `NAME_control_lambda.bdg` files are
   in bedGraph format which can be imported to the UCSC genome browser
   or be converted into even smaller bigWig files. The
   `NAME_treat_pielup.bdg` contains the pileup signals (normalized
   according to `--scale-to` option) from ChIP/treatment sample. The
   `NAME_control_lambda.bdg` contains local biases estimated for each
   genomic location from the control sample, or from treatment sample
   when the control sample is absent. The subcommand `bdgcmp` can be
   used to compare these two files and make a bedGraph file of scores
   such as p-value, q-value, log-likelihood, and log fold changes.

## Cutoff Analysis

Since cutoff can be an arbitrary value during peak calling, there are
many approaches proposed in the community to guide the cutoff
selection such as the [IDR
approach](https://doi.org/doi:10.1214%2F11-AOAS466). In MACS3, we
provide a simple way to do the cutoff analysis. The cutoff analysis
function is provided by `--cutoff-analysis` option in `callpeak`,
`bdgpeakcall`, and `hmmratac`. Among them, the function in
`bdgpeakcall` is more flexible and can be applied on any scoring
scheme. We will sperate this function into a dedicated subcommand in
the future.

Please note that if this `--cutoff-anlaysis` option is on, the report
will be written into a file named `NAME_cutoff_analysis.txt`.

When the option is on, we will generate a list of possible pvalue
cutoffs to check from pscore cutoff from 0.3 to 10, with a step of
0.3. When -log10(pvalue) is 0.3, it represents an extremely loose
cutoff pvalue 0.5; and when it's 10, it represents an extremely
strigent cutoff pvalue 1e-10. Please note that the is different with
`bdgpeakcall` where users can control how the cutoff should be
calculated.

Then for each cutoff we plan to investigate, we will check the number
of peaks that can be called, their average peak length, and their
total length.

The report consists of four columns:

1. score: the possible fold change cutoff value.
2. npeaks: the number of peaks under this cutoff.
3. lpeaks: the total length of all peaks.
4. avelpeak: the average length of peaks.

While there's no universal rule to suggest the best cutoff, here are a
few suggestions:

- You can use elbow analysis to find the cutoff that dramatically
   change the trend of npeaks, lpeaks, or avelpeak. But you need to
   think about how to define 'dramatical change'.
- You can use some common expectation to decide the cutoff. For
   example, the number of peaks should be thousands/ or the avelpeak
   should be around 500bps. Of course, it's arbitrary but the table
   will give you some insight.
