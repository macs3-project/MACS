# Bdgdiff

## Overview
The `bdgdiff` command is part of the MACS3 suite of tools and is used to call differential peaks from four bedGraph tracks for scores. It is particularly useful in ChIP-Seq analysis for identifying differential peaks in the data.

## Detailed Description

The `bdgdiff` command takes four input bedGraph files (two treatment and two control files) and produces three output files with differential peaks called. It uses an efficient algorithm to detect and call differential peaks, greatly improving the quality of your data for further analysis.

## Command Line Options

The command line options for `bdgdiff` are defined in `/MACS3/Commands/bdgdiff_cmd.py` and `/bin/macs3` files. Here is a brief overview of these options:

- `-t1` or `--t1bdg`: The first treatment file. This should be in bedGraph format. This option is required.
- `-c1` or `--c1bdg`: The first control file. This should be in bedGraph format. This option is required.
- `-t2` or `--t2bdg`: The second treatment file. This should be in bedGraph format. This option is required.
- `-c2` or `--c2bdg`: The second control file. This should be in bedGraph format. This option is required.
- `-d1` or `--depth1`: The depth of the first condition. This is a float value.
- `-d2` or `--depth2`: The depth of the second condition. This is a float value.
- `-o` or `--ofile`: The output file. This will be in bedGraph format.
- `-m` or `--minlen`: The minimum length of differential peaks. This is an integer value.
- `-g` or `--maxgap`: The maximum gap between differential peaks. This is an integer value.
- `-c` or `--cutoff`: The cutoff for LLR to call differential peaks. This is a float value.

## Example Usage

Here is an example of how to use the `bdgdiff` command:

```bash
macs3 bdgdiff -t1 treatment1.bedGraph -c1 control1.bedGraph -t2 treatment2.bedGraph -c2 control2.bedGraph --depth1 1.0 --depth2 1.0 -o output.bedGraph --minlen 500 --maxgap 1000 --cutoff 1.0
```

In this example, the program will call differential peaks from the two pairs of treatment and control bedGraph files and write the result to `output.bedGraph`. The depth of the first and second condition is set to 1.0, the minimum length of differential peaks is set to 500, the maximum gap between differential peaks is set to 1000, and the cutoff for LLR to call differential peaks is set to 1.0.

## FAQs about `bdgdiff`

Q: What does `bdgdiff` do?
A: `bdgdiff` is a tool in the MACS3 suite that calls differential peaks from four bedGraph tracks for scores. It is useful in ChIP-Seq analysis for identifying differential peaks in the data.

Q: How do I use `bdgdiff`?
A: You can use `bdgdiff` by providing it with four input files (two treatment and two control), an output file, and optional parameters that control its behavior. See the [Example Usage](#example-usage) section above for an example of how to use `bdgdiff`.

## Troubleshooting `bdgdiff`

If you're having trouble using `bdgdiff`, here are some things to try:

- Make sure your input files are in the correct format. `bdgdiff` requires four bedGraph format files.
- Make sure the values you provide for depth, minlen, maxgap, and cutoff are appropriate for your data. If the values are too small or too large, you may get inaccurate results.

## Known issues or limitations of `bdgdiff`

As of now, there are no known issues or limitations with `bdgdiff`. If you encounter a problem, please submit an issue on the MACS3 GitHub page.