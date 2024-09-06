# gappedPeak format

You can find the official definition of gappedPeak format
[here](https://genome.ucsc.edu/FAQ/FAQformat.html#format14). We chose
this format for broad peak calling, since it can represent a nested
peak calling result where there are some highlighted narrower regions
inside of a wider genomic region.

The fields/columns in the file are defined as followed (based on UCSC
definition, with specific notes on MACS3 output).

1. chrom - Name of the chromosome (or contig, scaffold, etc.), as used
   in the input alignment file for MACS3.
2. chromStart - The starting position of the feature in the chromosome
   or scaffold. The first base in a chromosome is numbered **0**.
3. chromEnd - The ending position of the feature in the chromosome or
   scaffold. The chromEnd base **is not included** in the display of
   the feature. For example, the first 100 bases of a chromosome are
   defined as chromStart=0, chromEnd=100, and span the bases numbered
   0-99. This BED-style coordination system is machine-friendly but
   not human-friendly. As a note, this is different with the
   coordination system used in GFF format file, where the same region
   in the above example will be defined as chromStart=1, chromEnd=100.
4. name - Name given to a region (preferably unique). In MACS3, it
   will be the experiment name plus a number such as 'FoxA1_99'.
5. score - Indicates how dark the peak will be displayed in the
   browser (0-1000). Thus, it's for the purpose of displaying on
   genome browser. In MACS3 `callpeak` output, we use the
   -log10qvalue*10. However, it may happen when the value in this
   column goes above 1000, and cause trouble while loading it in
   genome browsers. In this case, use the following `awk` command to
   fix: `awk -F'\t' '{ if ($5 > 1000) $5=1000; OFS="\t"; print }' peak.broadPeak`
6. strand - `+`/`-` to denote strand or orientation (whenever
   applicable). We use `.` since we can't decide the strand of DNA
   where the binding happens.
7. thickStart - The starting position at which the feature is drawn
   thickly. Not used in gappedPeak type, set to 0.
8. thickEnd - The ending position at which the feature is drawn
   thickly. Not used in gappedPeak type, set to 0.
9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). Not used
   in gappedPeak type, set to 0.
10. blockCount - The number of blocks (exons) in the BED line. In
    MACS3 broad peak calling mode, the 'blocks' are not exons but the
    stronger peaks within this peak region. Because of the restriction
    in this format -- see the blockStarts, we have to make the first
    base and the last base in the region in the 'blocks' if they are
    not. 
11. blockSizes - A comma-separated list of the block sizes. The number
    of items in this list should correspond to blockCount.
12. blockStarts - A comma-separated list of block starts. The first
	value must be 0 and all of the blockStart positions should be
	calculated relative to chromStart. The number of items in this
	list should correspond to blockCount.
13. signalValue - Measurement of overall (usually, average) enrichment
   for the region. In MACS3, we use the fold-enrichment values
   here. Please note that, same as the columns 8 and 9, this value
   corresponds to the mean value (in this case, the fold-enrichment)
   across the whole broad peak region.
14. pValue - Measurement of statistical significance (**-log10**). We
   are using a 'local Poisson' model where the Poisson lambda is
   estimated (by default) as the average input signal from surrounding
   region.
15. qValue - Measurement of statistical significance using false
   discovery rate (**-log10**). We used Benjamini-horchberg process
   over the whole genome, treating the p-value calculation on each
   basepair as an independent test, to get the FDR and q-value.

Please also note that, by default, we don't include a header line
(called trackline) in the gappedPeak output from MACS3. However, this
line will be necessary while uploading the gappedPeak file to UCSC for
display. Therefore, if you plan to upload the gappedPeak file, either
you turn on `--trackline` option in corresponding MACS3 subcommands,
or add the trackline mannually to the beginning of the file. A minimal
trackline is like:

`track type=gappedPeak name="track name" description="track description"`
