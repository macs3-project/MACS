# broadPeak format

You can find the official definition of broadPeak format
[here](https://genome.ucsc.edu/FAQ/FAQformat.html#format13). We chose
this format for broad peak calling results. It is a BED variant, and
comparing with narrowPeak format, it just misses the last peak summit
column.

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
7. signalValue - Measurement of overall (usually, average) enrichment
   for the region. In MACS3, we use the fold-enrichment values
   here. Please note that, same as the columns 8 and 9, this value
   corresponds to the mean value (in this case, the fold-enrichment)
   across the whole broad peak region.
8. pValue - Measurement of statistical significance (**-log10**). We
   are using a 'local Poisson' model where the Poisson lambda is
   estimated (by default) as the average input signal from surrounding
   region.
9. qValue - Measurement of statistical significance using false
   discovery rate (**-log10**). We used Benjamini-horchberg process
   over the whole genome, treating the p-value calculation on each
   basepair as an independent test, to get the FDR and q-value.

But please also note that, by default, we don't include a header line
(called trackline) in the broadPeak output from MACS3. However, this
line will be necessary while uploading the broadPeak file to UCSC for
display. Therefore, if you plan to upload the broadPeak file, either
you turn on `--trackline` option in corresponding MACS3 subcommands,
or add the trackline mannually to the beginning of the file. A minimal
trackline is like:

`track type=broadPeak name="track name" description="track description"`
