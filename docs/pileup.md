# Pileup alignment file and generate signal track

...

To see a full description of command-line options for `pileup`,
type `macs3 pileup -h`.

## Options


### `-i IFILE [IFILE ...]` or `--ifile IFILE [IFILE ...]`

Required option. Alignment file. If multiple files are given as '-t A
B C', then they will all be read and combined. 

### `-B` or `--both-direction`

By default, any read will be extended towards downstream direction by
extension size. So it's [0,size-1] (1-based index system) for plus
strand read and [-size+1,0] for minus strand read where position 0 is
5' end of the aligned read. Default behavior can simulate MACS3 way of
piling up ChIP sample reads where extension size is set as fragment
size or `d`. If this option is set as on, aligned reads will be
extended in both upstream and downstream directions by extension
size. It means [-size,size] where 0 is the 5' end of a aligned
read. It can partially simulate MACS3 way of piling up control
reads. However MACS3 local bias is calculated by maximizing the
expected pileup over a ChIP fragment size/d estimated from 10kb, 1kb,
d and whole genome background. This option will be ignored when the
format is set as BAMPE or BEDPE. Because in the case of `-f BAMPE` or
`-f BEDPE`, the regions defined by left and rightmost positions of
fragments will be piled up. DEFAULT: False

###  `--extsize EXTSIZE`

The extension size in bps. Each alignment read will become a EXTSIZE
of fragment, then be piled up. Check description for -B for
detail. This option will be ignored when the format is set as BAMPE or
BEDPE. DEFAULT: 200

###  `-o OUTPUTFILE` or  `--ofile OUTPUTFILE`

Output BedGraph file name. Required.

### `--outdir OUTDIR`

If specified all output files will be written to that
directory. Default: the current working directory

### `-f` or `--format`

Format of tag file. It can be "BED" or "ELAND" or "ELANDMULTI" or
"ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or "BAMPE" or "BEDPE". Or
you can set it as "AUTO". The default "AUTO" setting will let `macs3
randsample` decide which format the file is. Please check the
definition in [`callpeak` function](./callpeak.md) if you choose
ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE or BAMPE/BEDPE. DEFAULT:
"AUTO".

### `--buffer-size BUFFER_SIZE`
Buffer size for incrementally increasing internal array size to store
reads alignment information. In most cases, you don't have to change
this parameter.  However, if there are large number of
chromosomes/contigs/scaffolds in your alignment, it's recommended to
specify a smaller buffer size in order to decrease memory usage (but
it will take longer time to read alignment files). Minimum memory
requested for reading an alignment file is about # of CHROMOSOME *
BUFFER_SIZE * 8 Bytes. DEFAULT: 100000

## Example of usage

## Pileup aligned reads

Assuming aligned reads are all 50bps long:

```
$ macs3 pileup -i alignment.bed.gz --extsize 50 -o output.bdg
```

## Pileup the 5' sites of aligned reads

We need certain extension to smooth the signal. Otherwise, a single
base pileup will cause the signal track too sparse. Assuming we want a
100bp smoothing window, so that the value (4th column) in the
resulting bedGraph contains the total number of 5' sites within a
100bp window centered at each genomic location.

```
$ macs3 pileup -i alignment.bed.gz --extsize 50 -B -o output.bdg
```

## Pileup fragments defined by pairs of reads

## Pileup any defined genomic regions

## Generate a motif occurrence track from motif hits in BED format

