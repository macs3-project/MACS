# Cutoff Analysis

Since cutoff can be an arbitrary value during peak calling, there are
many approaches proposed in the community to guide the cutoff
selection such as the [IDR
approach](https://doi.org/doi:10.1214%2F11-AOAS466). In MACS3, we
provide a simple way to do the cutoff analysis. The cutoff analysis
function is provided by `--cutoff-analysis` option in `callpeak`,
`bdgpeakcall`, and `hmmratac`. Among them, the function in
`bdgpeakcall` is more flexible and can be applied on any scoring
scheme. We will use `bdgpeakcall` as example. 

The cutoff analysis will generate a list of possible cutoffs to check
from a loose cutoff to a stringent cutoff, with a given step. For
example, if the `--cutoff-analysis-min` is set as 0, the
`--cutoff-analysis-max` is set as 10, and the
`--cutoff-analysis-steps` is set as 20, then the list of cutoff values
to be investigated will range from 0 to 10, with an increase of 0.5
for each step.

Then for each cutoff we plan to investigate, we will check the number
of peaks that can be called, their average peak length, and their
total length. The exact report format depends on the subcommand.

The `bdgpeakcall` report used in this example consists of four columns:

1. score: the possible cutoff value for the scores in the input bedGraph.
2. npeaks: the number of peaks called at this cutoff.
3. lpeaks: the total length of all peaks.
4. avelpeak: the average length of peaks.

The `hmmratac` report uses the same four columns, with `score`
representing a fold-change cutoff. The `callpeak` report instead
contains five columns:

1. pscore: the possible `-log10(p-value)` cutoff.
2. qscore: the corresponding `-log10(q-value)` cutoff.
3. npeaks: the number of peaks called at this cutoff.
4. lpeaks: the total length of all peaks.
5. avelpeak: the average length of peaks.

While there's no universal rule to suggest the best cutoff, here are a
few suggestions:

- You can use elbow analysis to find the cutoff that dramatically
   change the trend of npeaks, lpeaks, or avelpeak. But you need to
   think about how to define 'dramatical change'.
- You can use some common expectation to decide the cutoff. For
   example, the number of peaks should be thousands/ or the avelpeak
   should be around 500bps. Of course, it's arbitrary but the table
   will give you some insight.

Please check the document for [`hmmratac`](./hmmratac.md) function for
choosing lower, upper and pre-scan cutoff values from cutoff-analysis
result.
