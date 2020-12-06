## Getting Started
```sh
git clone https://github.com/lh3/fermi-lite
cd fermi-lite && make
./fml-asm test/MT-simu.fq.gz > MT.fq
# to compile your program:
gcc -Wall -O2 prog.c -o prog -L/path/to/fermi-lite -lfml -lz -lm -lpthread
```

## Introduction

Fermi-lite is a standalone C library as well as a command-line tool for
assembling Illumina short reads in regions from 100bp to 10 million bp in size.
It is largely a light-weight in-memory version of [fermikit][fk] without
generating any intermediate files. It inherits the performance, the relatively
small memory footprint and the features of fermikit. In particular, fermi-lite
is able to retain heterozygous events and thus can be used to assemble diploid
regions for the purpose of variant calling. It is one of the limited choices
for local re-assembly and arguably the easiest to interface.

## Usage

For now, see [example.c][example] for the basic use of the library. Here is a
sketch of the example:
```cpp
#include <stdio.h>                      // for printf()
#include "fml.h"                        // only one header file required

int main(int argc, char *argv[])
{
	int i, n_seqs, n_utgs;
	bseq1_t *seqs;                      // array of input sequences
	fml_utg_t *utgs;                    // array of output unitigs
	fml_opt_t opt;
	if (argc == 1) return 1;            // do nothing if there is no input file
	seqs = bseq_read(argv[1], &n_seqs); // or fill the array with callers' functions
	fml_opt_init(&opt);                 // initialize parameters
	utgs = fml_assemble(&opt, n_seqs, seqs, &n_utgs); // assemble!
	for (i = 0; i < n_utgs; ++i)        // output in fasta
		printf(">%d\n%s\n", i+1, utgs[i].seq);
	fml_utg_destroy(n_utgs, utgs);      // deallocate unitigs
	return 0;
}
```
The `fml_assemble()` output is in fact a graph. You may have a look at the
`fml_utg_print_gfa()` function in [misc.c][misc] about how to derive a
[GFA][gfa] representation from an array of `fml_utg_t` objects.

## Overview of the Assembly Algorithm

Fermi-lite is an overlap-based assembler. Given a set of input reads, it counts
*k*-mers, estimates the *k*-mer coverage, sets a threshold on *k*-mer
occurrences to determine solid *k*-mers and then use them correct sequencing
errors ([Li, 2015][bfc-paper]). After error correction, fermi-lite trims a read
at an *l*-mer unique to the read. It then constructs an FM-index for trimmed
reads ([Li, 2014][rb2-paper]) and builds a transitively reduced overlap graph from the
FM-index ([Simpson and Durbin, 2010][sga-paper]; [Li, 2012][fm1-paper]),
requiring at least *l*-bp overlaps. In this graph, fermi-lite trims tips and
pops bubbles caused by uncorrected errors. If a sequence in the graph has
multiple overlaps, fermi-lite discards overlaps significantly shorter than the
longest overlap -- this is a technique applied to overlap graph only. The graph
after these procedure is the final output. Sequences in this graph are unitigs.

## Limitations

1. Fermi-lite can efficiently assemble bacterial genomes. However, it has not
   been carefully tuned for this type of assembly. While on a few GAGE-B data
   sets fermi-lite appears to work well, it may not compete with recent
   mainstream assemblers in general.

2. Fermi-lite does not work with genomes more than tens of megabases as a
   whole. It would take too much memory to stage all data in memory. For large
   genomes, please use [fermikit][fk] instead.

3. This is the first iteration of fermi-lite. It is still immarture. In
   particular, I hope fermi-lite can be smart enough to automatically figure
   out various parameters based on input, which is very challenging given the
   high variability of input data.

[sga-paper]: http://www.ncbi.nlm.nih.gov/pubmed/20529929
[bfc-paper]: http://www.ncbi.nlm.nih.gov/pubmed/25953801
[rb2-paper]: http://www.ncbi.nlm.nih.gov/pubmed/25107872
[fm1-paper]: http://www.ncbi.nlm.nih.gov/pubmed/22569178
[bfc]: http://github.com/lh3/bfc
[rb2]: http://github.com/lh3/ropebwt2
[fm2]: http://github.com/lh3/fermi2
[fk]: http://github.com/lh3/fermikit
[example]: https://github.com/lh3/fermi-lite/blob/master/example.c
[header]: https://github.com/lh3/fermi-lite/blob/master/fml.h
[misc]: https://github.com/lh3/fermi-lite/blob/master/misc.c
[gfa]: https://github.com/pmelsted/GFA-spec
