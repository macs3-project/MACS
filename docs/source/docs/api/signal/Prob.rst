MACS3.Signal.Prob module
========================

.. automodule:: MACS3.Signal.Prob
   :members:
   :undoc-members:
   :show-inheritance:

Concurrent Poisson scoring
--------------------------

``poisson_cdf`` (log10 upper tail), ``poisson_cdf_Q_log10_many``,
``poisson_cdf_Q_log10_prepare``, and ``poisson_cdf_Q_log10_range`` may run
concurrently. Each scoring call retains an immutable snapshot of the
logarithm tables until its native computation finishes. Range and batch
computation release the GIL. Use separate output buffers or disjoint output
ranges, and keep input arrays unchanged while calls are running.

Preparation is optional and can reduce later allocation work. Growth creates
both replacement tables before publishing them under the GIL. If allocation
fails, the previous tables remain valid and scoring uses the ordered uncached
calculation for entries outside that snapshot. Preparation returns ``False``
when its requested table prefix could not be allocated. These paths preserve
the scalar calculation's operation order and five-decimal rounding.

At the default limit of 2**22 entries, the two float64 tables occupy at most
64 MiB per generation. The current generation remains cached for reuse.
Retired generations are released when their last active reader finishes;
construction and overlapping readers can temporarily use more than 64 MiB.
