========================
README for MACS (2.1.0)
========================
Time-stamp: <2015-04-20 15:37:09 Tao Liu>

Introduction
============

With the improvement of sequencing techniques, chromatin
immunoprecipitation followed by high throughput sequencing (ChIP-Seq)
is getting popular to study genome-wide protein-DNA interactions. To
address the lack of powerful ChIP-Seq analysis method, we present a
novel algorithm, named Model-based Analysis of ChIP-Seq (MACS), for
identifying transcript factor binding sites. MACS captures the
influence of genome complexity to evaluate the significance of
enriched ChIP regions, and MACS improves the spatial resolution of
binding sites through combining the information of both sequencing tag
position and orientation. MACS can be easily used for ChIP-Seq data
alone, or with control sample with the increase of specificity.

Install
=======

Please check the file 'INSTALL' in the distribution.

Usage of MACS2
==============

::

  macs2 [-h] [--version]
             {callpeak,filterdup,bdgpeakcall,bdgcmp,randsample,bdgdiff,bdgbroadcall}

:Example for regular peak calling: ``macs2 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01``

:Example for broad peak calling: ``macs2 callpeak -t ChIP.bam -c Control.bam --broad -g hs --broad-cutoff 0.1``

There are seven major functions available in MACS serving as sub-commands.

:callpeak:            Main MACS2 Function to `Call peaks`_ from alignment results.
:bdgpeakcall:         Call peaks from bedGraph output.
:bdgbroadcall:        Call broad peaks from bedGraph output.
:bdgcmp:              Deduct noise by comparing two signal tracks in bedGraph.
:bdgdiff:             Differential peak detection based on paired four bedgraph files.
:filterdup:           Remove duplicate reads at the same position, then convert acceptable format to BED format.
:predictd:            Predict d or fragment size from alignment results.
:pileup:              Pileup aligned reads with a given extension
                      size (fragment size or d in MACS language). Note there will be no
                      step for duplicate reads filtering or sequencing depth scaling, so you may need to do certain post-
                      processing.
:randsample:          Randomly sample number/percentage of total reads.
:refinepeak:          (Experimental) Take raw reads alignment, refine peak
                          summits and give scores measuring balance of forward-
                          backward tags. Inspired by SPP.


We only cover 'callpeak' module in this document. Please use 'macs2
COMMAND -h' to see the detail description for each option of each
module.

Call peaks
~~~~~~~~~~

This is the main function in MACS2. It can be invoked by 'macs2
callpeak' command. If you type this command without parameters, you
will see a full description of commandline options. Here we only list
commonly used ones.

Options
--------------

-t/--treatment FILENAME
```````````````````````

This is the only REQUIRED parameter for MACS. File can be in any
supported format specified by --format option. Check --format for
detail. If you have more than one alignment files, you can specify
them as ```-t A B C```. MACS will pool up all these files together.

-c/--control
````````````

The control or mock data file. Please follow the same direction as for
-t/--treatment.

-n/--name
`````````

The name string of the experiment. MACS will use this string NAME to
create output files like 'NAME_peaks.xls', 'NAME_negative_peaks.xls',
'NAME_peaks.bed' , 'NAME_summits.bed', 'NAME_model.r' and so on. So
please avoid any confliction between these filenames and your
existing files.

--outdir
````````

MACS2 will save all output files into speficied folder for this
option.

-f/--format FORMAT
``````````````````

Format of tag file, can be "ELAND", "BED", "ELANDMULTI",
"ELANDEXPORT", "ELANDMULTIPET" (for pair-end tags), "SAM", "BAM",
"BOWTIE" or "BAMPE". Default is "AUTO" which will allow MACS to decide
the format automatically. "AUTO" is also usefule when you combine
different formats of files.

The BED format can be found at `UCSC genome browser website <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_.

If the format is ELAND, the file must be ELAND result output file,
each line MUST represents only ONE tag, with fields of:

1. Sequence name (derived from file name and line number if format is not Fasta)
2. Sequence
3. Type of match:

:NM: no match found.
:QC: no matching done: QC failure (too many Ns basically).
:RM: no matching done: repeat masked (may be seen if repeatFile.txt was specified).
:U0: Best match found was a unique exact match.
:U1: Best match found was a unique 1-error match. 
:U2: Best match found was a unique 2-error match. 
:R0: Multiple exact matches found.
:R1: Multiple 1-error matches found, no exact matches.
:R2: Multiple 2-error matches found, no exact or 1-error matches.

4. Number of exact matches found.
5. Number of 1-error matches found.
6. Number of 2-error matches found.  
   Rest of fields are only seen if a unique best match was found
   (i.e. the match code in field 3 begins with "U").
7. Genome file in which match was found.
8. Position of match (bases in file are numbered starting at 1).
9. Direction of match (F=forward strand, R=reverse).
10. How N characters in read were interpreted: ("."=not applicable,
    "D"=deletion, "I"=insertion). Rest of fields are only seen in
    the case of a unique inexact match (i.e. the match code was U1 or
    U2).
11. Position and type of first substitution error (e.g. 12A: base 12
    was A, not whatever is was in read).
12. Position and type of first substitution error, as above. 

If the format is ELANDMULTI, the file must be ELAND output file from
multiple-match mode, each line MUST represents only ONE tag, with
fields of:

1. Sequence name 
2. Sequence 
3. Either NM, QC, RM (as described above) or the following: 
4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches found
5. Blank, if no matches found or if too many matches found, or the following:
   BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 This says
   there are two matches to BAC_plus_vector.fa: one in the reverse
   direction starting at position 160322 with one error, one in the
   forward direction starting at position 170128 with two
   errors. There is also a single-error match to E_coli.fa.

If the format is BAM/SAM, please check the definition in
(http://samtools.sourceforge.net/samtools.shtml).  Pair-end mapping
results can be saved in a single BAM file, if so, MACS will
automatically keep the left mate(5' end) tag. However, when format
BAMPE is specified, MACS will use the real fragments inferred from
alignment results for reads pileup.

If the format is BOWTIE, you need to provide the ASCII bowtie output
file with the suffix '.map'. Please note that, you need to make sure
that in the bowtie output, you only keep one location for one
read. Check the bowtie manual for detail if you want at
(http://bowtie-bio.sourceforge.net/manual.shtml)

Here is the definition for Bowtie output in ASCII characters I copied
from the above webpage:

1. Name of read that aligned
2. Orientation of read in the alignment, '-' for reverse complement, '+'
   otherwise
3. Name of reference sequence where alignment occurs, or ordinal ID
   if no name was provided
4. 0-based offset into the forward reference strand where leftmost
   character of the alignment occurs
5. Read sequence (reverse-complemented if orientation is -)
6. ASCII-encoded read qualities (reversed if orientation is -). The
   encoded quality values are on the Phred scale and the encoding is
   ASCII-offset by 33 (ASCII char !).
7. Number of other instances where the same read aligns against the
   same reference characters as were aligned against in this
   alignment. This is not the number of other places the read aligns
   with the same number of mismatches. The number in this column is
   generally not a good proxy for that number (e.g., the number in
   this column may be '0' while the number of other alignments with
   the same number of mismatches might be large). This column was
   previously described as "Reserved".
8. Comma-separated list of mismatch descriptors. If there are no
   mismatches in the alignment, this field is empty. A single
   descriptor has the format offset:reference-base>read-base. The
   offset is expressed as a 0-based offset from the high-quality (5')
   end of the read.

Notes:

1) For BED format, the 6th column of strand information is required by
MACS. And please pay attention that the coordinates in BED format is
zero-based and half-open
(http://genome.ucsc.edu/FAQ/FAQtracks#tracks1).

2) For plain ELAND format, only matches with match type U0, U1 or U2
is accepted by MACS, i.e. only the unique match for a sequence with
less than 3 errors is involed in calculation. If multiple hits of a
single tag are included in your raw ELAND file, please remove the
redundancy to keep the best hit for that sequencing tag.

3) For the experiment with several replicates, it is recommended to
concatenate several ChIP-seq treatment files into a single file. To
do this, under Unix/Mac or Cygwin (for windows OS), type:

```$ cat replicate1.bed replicate2.bed replicate3.bed > all_replicates.bed```

For BAM or SAM files, samtools can be used to combine replicates.

4) ELAND export format support sometimes may not work on your
datasets, because people may mislabel the 11th and 12th column. MACS
uses 11th column as the sequence name which should be the chromosome
names.

5) A special mode will be triggered while format is specified as
'BAMPE'. In this way, MACS2 will process the BAM files as paired-end
data. Instead of building bimodal distribution of plus and minus
strand reads to predict fragment size, MACS2 now will use actual
insert sizes of pairs of reads to build fragment pileup. 


-g/--gsize
``````````

PLEASE assign this parameter to fit your needs!

It's the mappable genome size or effective genome size which is
defined as the genome size which can be sequenced. Because of the
repetitive features on the chromsomes, the actual mappable genome size
will be smaller than the original size, about 90% or 70% of the genome
size. The default hs -- 2.7e9 is recommended for UCSC human hg18
assembly. Here are all precompiled parameters for effective genome
size:

:hs: 2.7e9
:mm: 1.87e9
:ce: 9e7
:dm: 1.2e8

-s/--tsize
``````````

The size of sequencing tags. If you don't specify it, MACS will try to
use the first 10 sequences from your input treatment file to determine
the tag size. Specifying it will override the automatically determined
tag size.

--bw
````

The band width which is used to scan the genome ONLY for model
building. You can set this parameter as the sonication fragment size
expected from wet experiment. The previous side effect on the peak
detection process has been removed. So this parameter only affects the
model building.

-q/--qvalue
```````````

The qvalue (minimum FDR) cutoff to call significant regions. Default
is 0.01. For broad marks, you can try 0.05 as cutoff. Q-values are
calculated from p-values using Benjamini-Hochberg procedure.

-p/--pvalue
```````````

The pvalue cutoff. If -p is specified, MACS2 will use pvalue instead
of qvalue.

-m/--mfold
``````````

This parameter is used to select the regions within MFOLD range of
high-confidence enrichment ratio against background to build
model. The regions must be lower than upper limit, and higher than
the lower limit of fold enrichment. DEFAULT:10,30 means using all
regions not too low (>10) and not too high (<30) to build
paired-peaks model. If MACS can not find more than 100 regions to
build model, it will use the --extsize parameter to continue the
peak detection ONLY if --fix-bimodal is set.


--nolambda
``````````

With this flag on, MACS will use the background lambda as local
lambda. This means MACS will not consider the local bias at peak
candidate regions.

--slocal, --llocal
``````````````````

These two parameters control which two levels of regions will be
checked around the peak regions to calculate the maximum lambda as
local lambda. By default, MACS considers 1000bp for small local
region(--slocal), and 10000bps for large local region(--llocal) which
captures the bias from a long range effect like an open chromatin
domain. You can tweak these according to your project. Remember that
if the region is set too small, a sharp spike in the input data may
kill the significant peak.

--fix-bimodal
`````````````

Whether turn on the auto paired-peak model process. If it's set, when
MACS failed to build paired model, it will use the nomodel settings,
the '--extsize' parameter to extend each tags. If set, MACS will be
terminated if paried-peak model is failed.

--nomodel
`````````

While on, MACS will bypass building the shifting model.

--extsize
`````````

While '--nomodel' is set, MACS uses this parameter to extend reads in
5'->3' direction to fix-sized fragments. For example, if the size of
binding region for your transcription factor is 200 bp, and you want
to bypass the model building by MACS, this parameter can be set
as 200. This option is only valid when --nomodel is set or when MACS
fails to build model and --fix-bimodal is on.

--shift
```````

Note, this is NOT the legacy --shiftsize option which is replaced by
--extsize! You can set an arbitrary shift in bp here. Please Use
discretion while setting it other than default value (0). When
--nomodel is set, MACS will use this value to move cutting ends (5')
then apply --extsize from 5' to 3' direction to extend them to
fragments. When this value is negative, ends will be moved toward
3'->5' direction, otherwise 5'->3' direction. Recommended to keep it
as default 0 for ChIP-Seq datasets, or -1 * half of EXTSIZE together
with --extsize option for detecting enriched cutting loci such as
certain DNAseI-Seq datasets. Note, you can't set values other than 0
if format is BAMPE for paired-end data. Default is 0.

Here are some examples for combining --shift and --extsize:

1. To find enriched cutting sites such as some DNAse-Seq datasets. In
this case, all 5' ends of sequenced reads should be extended in both
direction to smooth the pileup signals. If the wanted smoothing window
is 200bps, then use '--nomodel --shift -100 --extsize 200'.

2. For certain nucleosome-seq data, we need to pileup the centers of
nucleosomes using a half-nucleosome size for wavelet analysis
(e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about
147bps, this option can be used: '--nomodel --shift 37 --extsize 73'.

--keep-dup
``````````

It controls the MACS behavior towards duplicate tags at the exact same
location -- the same coordination and the same strand. The default
'auto' option makes MACS calculate the maximum tags at the exact same
location based on binomal distribution using 1e-5 as pvalue cutoff;
and the 'all' option keeps every tags.  If an integer is given, at
most this number of tags will be kept at the same location. The
default is to keep one tag at the same location. Default: 1

--broad
```````

When this flag is on, MACS will try to composite broad regions in
BED12 ( a gene-model-like format ) by putting nearby highly enriched
regions into a broad region with loose cutoff. The broad region is
controlled by another cutoff through --broad-cutoff. The maximum
length of broad region length is 4 times of d from MACS. DEFAULT:
False

--broad-cutoff
``````````````

Cutoff for broad region. This option is not available unless --broad
is set. If -p is set, this is a pvalue cutoff, otherwise, it's a
qvalue cutoff.  DEFAULT: 0.1

--to-large
``````````

When set, linearly scale the smaller dataset to the same depth as
larger dataset, by default, the smaller dataset will be scaled
towards the larger dataset. Beware, to scale up small data would
cause more false positives.

--down-sample
`````````````

When set, random sampling method will scale down the bigger
sample. By default, MACS uses linear scaling. This option will make
the results unstable and irreproducible since each time, random reads
would be selected, especially the numbers (pileup, pvalue, qvalue)
would change. Consider to use 'randsample' script before MACS2 runs
instead.

-B/--bdg
````````

If this flag is on, MACS will store the fragment pileup, control
lambda, -log10pvalue and -log10qvalue scores in bedGraph files. The
bedGraph files will be stored in current directory named
NAME+'_treat_pileup.bdg' for treatment data,
NAME+'_control_lambda.bdg' for local lambda values from control,
NAME+'_treat_pvalue.bdg' for Poisson pvalue scores (in -log10(pvalue)
form), and NAME+'_treat_qvalue.bdg' for q-value scores from
Benjamini–Hochberg–Yekutieli procedure
<http://en.wikipedia.org/wiki/False_discovery_rate#Dependent_tests>

--call-summits
``````````````

MACS will now reanalyze the shape of signal profile (p or q-score
depending on cutoff setting) to deconvolve subpeaks within each peak
called from general procedure. It's highly recommended to detect
adjacent binding events. While used, the output subpeaks of a big
peak region will have the same peak boundaries, and different scores
and peak summit positions.

--verbose
`````````

If you don't want to see any message during the running of MACS, set
it to 0. But the CRITICAL messages will never be hidden. If you want
to see rich information like how many peaks are called for every
chromosome, you can set it to 3 or larger than 3.

Output files
~~~~~~~~~~~~

1. NAME_peaks.xls is a tabular file which contains information about
   called peaks. You can open it in excel and sort/filter using excel
   functions. Information include: chromosome name, start position of
   peak, end position of peak, length of peak region, absolute peak
   summit position, pileup height at peak summit, -log10(pvalue) for
   the peak summit (e.g. pvalue =1e-10, then this value should be
   10), fold enrichment for this peak summit against random Poisson
   distribution with local lambda, -log10(qvalue) at peak
   summit. Coordinates in XLS is 1-based which is different with BED
   format.

2. NAME_peaks.narrowPeak is BED6+4 format file which contains the
   peak locations together with peak summit, pvalue and qvalue. You
   can load it to UCSC genome browser. Definition of some specific
   columns are: 5th: integer score for display, 7th: fold-change,
   8th: -log10pvalue, 9th: -log10qvalue, 10th: relative summit
   position to peak start. The file can be loaded directly to UCSC
   genome browser. Remove the beginning track line if you want to
   analyze it by other tools.

3. NAME_summits.bed is in BED format, which contains the peak summits
   locations for every peaks. The 5th column in this file is
   -log10pvalue the same as NAME_peaks.bed. If you want to find the
   motifs at the binding sites, this file is recommended. The file
   can be loaded directly to UCSC genome browser. Remove the
   beginning track line if you want to analyze it by other tools.

4. NAME_peaks.broadPeak is in BED6+3 format which is similar to
   narrowPeak file, except for missing the 10th column for annotating
   peak summits.

4. NAME_peaks.gappedPeak is in BED12+3 format which contains both the
   broad region and narrow peaks. The 5th column is 10*-log10qvalue,
   to be more compatible to show grey levels on UCSC browser. Tht 7th
   is the start of the first narrow peak in the region, and the 8th
   column is the end. The 9th column should be RGB color key, however,
   we keep 0 here to use the default color, so change it if you
   want. The 10th column tells how many blocks including the starting
   1bp and ending 1bp of broad regions. The 11th column shows the
   length of each blocks, and 12th for the starts of each blocks. 13th:
   fold-change, 14th: -log10pvalue, 15th: -log10qvalue. The file can be
   loaded directly to UCSC genome browser. 

5. NAME_model.r is an R script which you can use to produce a PDF
   image about the model based on your data. Load it to R by:

   ```$ Rscript NAME_model.r```

   Then a pdf file NAME_model.pdf will be generated in your current
   directory. Note, R is required to draw this figure.

6. The .bdg files are in bedGraph format which can be imported to
   UCSC genome browser or be converted into even smaller bigWig
   files. There are two kinds of bdg files: treat_pileup, and
   control_lambda.

Other useful links
==================

:Cistrome: http://cistrome.org/ap/
:bedTools: http://code.google.com/p/bedtools/
:UCSC toolkits: http://hgdownload.cse.ucsc.edu/admin/exe/

Tips of fine-tuning peak calling
================================

Check the three scripts within MACSv2 package:

1. bdgcmp can be used on ```*_treat_pileup.bdg``` and
   ```*_control_lambda.bdg``` or bedGraph files from other resources
   to calculate score track.

2. bdgpeakcall can be used on ```*_treat_pvalue.bdg``` or the file
   generated from bdgcmp or bedGraph file from other resources to
   call peaks with given cutoff, maximum-gap between nearby mergable
   peaks and minimum length of peak. bdgbroadcall works similarly to
   bdgpeakcall, however it will output _broad_peaks.bed in BED12
   format.

3. Differential calling tool -- bdgdiff, can be used on 4 bedgraph
   files which are scores between treatment 1 and control 1,
   treatment 2 and control 2, treatment 1 and treatment 2, treatment
   2 and treatment 1. It will output the consistent and unique sites
   according to parameter settings for minimum length, maximum gap
   and cutoff.
