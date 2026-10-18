MACS3.Signal.CallPeakUnit module
================================

``CallerFromAlignments`` uses the shared invocation execution policy. Process
phases also require at least 262,144 stored treatment and distinct control
locations across two nonempty chromosomes. FRAG counts do not multiply the
stored-location cost. This minimum is an internal startup heuristic.

Chromosomes empty in both inputs are ignored; an empty comparison returns an
empty peak collection. One-sided empty chromosomes retain zero treatment or
background control coverage. Large chromosomes remain individual tasks and
small ones are batched, with at most 32 tasks per eligible worker. Results
retain serial chromosome and peak order. See :doc:`../../execution` for
resource limits and empty output artifacts.

.. automodule:: MACS3.Signal.CallPeakUnit
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: bool, erf, exp, float, int, log, math, np, os, pd, re, sys, fopen, fprintf, get_logFE, 
                     get_logLR_asym,get_pscore, chi2_k1_cdf, chi2_k2_cdf, chi2_k4_cdf, mean_from_value_length,
                     median_from_value_length, right_forward, right_sum, apply_multiple_cutoffs, get_subtraction,
                     left_forward,  left_sum, log10_chi2_k1_cdf, log10_chi2_k2_cdf, log10_chi2_k4_CDF, log10, log1p, sqrt