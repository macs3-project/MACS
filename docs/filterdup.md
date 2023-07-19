# Filterdup

## Overview
The `filterdup` command is part of the MACS3 suite of tools and is used to filter duplicate reads from your data. It is particularly useful in sequencing analysis where duplicate reads can bias the downstream analysis.

## Detailed Description

The `filterdup` command takes an input file and produces an output file with duplicate reads removed. It uses an efficient algorithm to detect and filter duplicate reads, greatly improving the quality of your data for further analysis.

## Command Line Options

The command line options for `filterdup` are defined in `/MACS3/Commands/filterdup_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-i` or `--ifile`: The input file. This should be in BAM format. This option is required.
- `-o` or `--ofile`: The output file. This will be in the same format as the input file. This option is required.
- `-g` or `--gsize`: The mappable genome size. This can be an integer or a string. If it's an integer, it represents the genome size. If it's a string, it should be one of the following: `hs` (for Homo Sapiens), `mm` (for Mus musculus), `ce` (for Caenorhabditis elegans), or `dm` (for Drosophila melanogaster). The default is `hs`.
- `-s` or `--format`: The format of the input file. This can be `AUTO`, `BED`, `ELAND`, `ELANDMULTI`, `ELANDEXPORT`, `SAM`, `BAM`, or `BOWTIE`. The default is `AUTO`, which means the program will try to determine the format automatically.
- `--keep-dup`: The number of duplicates to keep. This can be `all`, `auto`, or an integer. If it's `all`, all duplicates will be kept. If it's `auto`, the program will use an algorithm to determine the number of duplicates to keep. If it's an integer, it represents the number of duplicates to keep. The default is `1`.
- `--buffer-size`: Buffer size for incrementally increasing internal array size to store reads alignment information. Default: `100000`.

## Example Usage

Here is an example of how to use the `filterdup` command:

```bash
macs3 filterdup -i input.bam -o output.bam --gsize hs --format AUTO --keep-dup 1 --buffer-size 100000
```

In this example, the program will remove duplicate reads from the `input.bam` file and write the result to `output.bam`. The mappable genome size is set to `hs` (Homo Sapiens), the format of the input file is determined automatically, and the program keeps only one duplicate.

## FAQs about `filterdup`

Q: What does `filterdup` do?
A: `filterdup` is a tool in the MACS3 suite that filters duplicate reads from your data. Duplicate reads can bias your analysis, so it's important to remove them.

Q: How do I use `filterdup`?
A: You can use `filterdup` by providing it with an input file, an output file, and optional parameters that control its behavior. See the [Example Usage](#example-usage) section above for an example of how to use `filterdup`.

## Troubleshooting `filterdup`

If you're having trouble using `filterdup`, here are some things to try:

- Make sure your input file is in the correct format. `filterdup` can handle several formats, and it can try to determine the format automatically, but it's best to specify the format explicitly if possible.
- Make sure the `--gsize` option matches your genome. If you're not working with one of the four supported genomes (`hs`, `mm`, `ce`, `dm`), you will need to specify the genome size as an integer.
- If the program is running out of memory or taking too long, try adjusting the `--buffer-size` option.

## Known issues or limitations of `filterdup`

As of now, there are no known issues or limitations with `filterdup`. If you encounter a problem, please submit an issue on the MACS3 GitHub page.