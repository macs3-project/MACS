# SAM/BAM/BAMPE format

SAM and BAM formats are the most common format used to store alignment
information of sequencing reads. The BAM format is in fact the binary
version of SAM which is a text file. Please refer to [SAM format
specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for
detail. As in MACS3, we call the BAM file from paired-end reads
mapping as `BAMPE` -- BAM file for Paired-End data. When you specify
the input format as `BAMPE` such as `-f BAMPE` in `callpeaks`, you
will trigger the PE mode for MACS3 where we will consider the whole
fragment between read pair as the DNA fragment bound by the protein of
interest.  On the contrast, even if the BAM storing PE mapping
information, if you specify `-f BAM`, MACS3 will treat the input as
single-end data and use only the 1st mate of read pair to estimate the
fragment length using the cross-correlation model.

Most of MACS3 modules only take the mapping location of reads from
SAM/BAM/BAMPE file. If necessary, you can use `macs3 filterdup
--keep-dup all` or `macs3 randsample -p 100` to convert the
SAM/BAM/BAMPE into BED or [BEDPE](./BEDPE.md) format, then gzip them
to save storage. 

The only exception is that the `callvar` command in MACS3 will need to
read the actual sequences and mapping information such as the
mismatch, insertion, deletion, etc from BAM. Also, please make sure
that the BAM file for `callvar` has to be sorted and indexed.
