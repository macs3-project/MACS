# bedGraph format

The official definition of bedGraph can be found
[here](https://genome.ucsc.edu/goldenPath/help/bedgraph.html). As for
all other files generated directly from MACS3, you may need to add a
trackline to the beginning of the file when you upload the file to
UCSC genome browser for display. A minimal trackline for bedGraph is
`track type=bedGraph`. Please read the above link to the official
definition page for more instruction. The bedGraph format file can be
further converted into a binary format named 'bigWig' for smaller file
size and for fast and random access. Check this
[gist](https://gist.github.com/taoliu/2469050) for a script to convert
the plain text bedGraph to binary bigWig.
