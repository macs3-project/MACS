# BEDPE format

The BEDPE format is specifically designed for keeping the alignment
locations of each read pair from Paired-End library. This is not a
general format but only a format for MACS3. It only contains three
columns -- the chromosome, the leftmost position of read pair, and the
rightmost position of the read pair. These three columns All other information from
alignment will not be kept in this format, such as the length of the
read, the mismatches/gaps in the alignment, and etc. An example is as
followed:

```
chrXIII	0	60
chrXIII	1	64
chrXIII	1	211
chrXIII	2	46
chrXIII	3	154
chrXIII	3	209
chrXIII	9	71
chrXIII	11	67
chrXIII	11	71
chrXIII	14	71
...
```


