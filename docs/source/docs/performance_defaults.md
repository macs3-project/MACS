# Performance defaults (developer reference)

`MACS3.Utilities.Performance.DEFAULTS` is an immutable, slotted record of
internal policy defaults. These values retain the tuning introduced in
PR #756. They are not CLI options or per-invocation settings. Every
compiled consumer copies the relevant field into an explicitly typed
module value at import, outside native and `nogil` loops. Changing a
module attribute cannot consistently configure already imported consumers.

`Threads.MAX_THREADS` and `Threads.MIN_THREADED_SIZE` remain derived
compatibility aliases, as do `CallPeakUnit.MAX_WORKERS` and
`CallPeakUnit.PSCORE_CACHE_INIT_SIZE`. Compiled consumers use typed copies.
CPU detection and backend eligibility are unchanged.
The thread caps include the calling thread; the process cap counts child
processes. A shared configurable CPU budget is tracked in #768, and the
positive scheduling target/minimum-work policy in #769.

## Python policy inventory

Sizes and multipliers are positive native-size integers. Hash capacities
are powers of two; CPU and FRAG task counts and the Poisson cache bound fit
signed 32-bit integers. The Poisson minimum does not exceed its maximum.
The record validates these constraints once at construction with explicit
errors, including under `python -O`. All reserve growth factors and load
factors remain algorithm-owned and unchanged.

| DEFAULTS field | Existing value and units | Current definition / consumers | Rationale / constraint |
| --- | --- | --- | --- |
| cpu_limit | 8 execution participants | Threads.MAX_THREADS, Parser._BAM_MAX_THREADS, PeakModel._MODEL_MAX_THREADS, CallPeakUnit.MAX_WORKERS | Preserve each resolver's accounting: threads include caller; process limit counts children. Shared request/accounting belongs to #768. |
| min_parallel_track_elements | 65,536 elements per chromosome | Threads.MIN_THREADED_SIZE; typed copies in FixWidthTrack and PairedEndTrack | Avoid thread handoff on small chromosomes; preserve >= comparison. |
| bam_serial_input_bytes | 1,048,576 compressed bytes | Parser._BAM_IN_SIZE / shared BAM/FRAG stream | Positive initial read buffer; retain stream growth/fallback. |
| bam_serial_output_bytes | 4,194,304 decompressed bytes | Parser._BAM_OUT_SIZE / shared BAM/FRAG stream | Initial window, not a record-size limit. |
| bgzf_batch_input_bytes | 4,194,304 compressed bytes | Parser._BAM_PAR_IN | Parallel compressed batch capacity. |
| bgzf_batch_output_bytes | 16,777,216 decompressed bytes | Parser._BAM_PAR_OUT | Parallel output budget, unrelated to float32 exactness. |
| bgzf_tail_reserve_bytes | 262,144 decompressed bytes | Parser._BAM_PAR_HEAD | Prefix space for residual records; longer records still grow/fallback. |
| bgzf_max_batch_jobs | 8,192 BGZF blocks/jobs | Parser._BAM_PAR_MAXJOBS / bgzf_split job allocation and bound | Allocate and dispatch with the same capacity; preserve stop/fallback behavior. |
| frag_append_batch_records | 262,144 fragments | Parser._FRAG_BATCH / _FragBatch | Positive append staging capacity, not bytes. |
| frag_parallel_min_bytes | 65,536 decompressed bytes | Parser._FRAG_PAR_MIN | Existing threaded line-walk eligibility threshold. |
| frag_tasks_per_thread | 2 chunks per inflate thread | Parser._FRAG_TASKS_PER_THREAD / chunk sizing | Positive multiplier; preserve current pool reuse. |
| frag_interner_initial_slots | 1,024 hash slots | Parser._Interner.__cinit__ / new_slots | Positive power of two; retain quarter-full growth and mask/shift rules. |
| frag_interner_initial_arena_bytes | 16,384 key bytes | Parser._Interner.add | Initial long-key arena; retain doubling. |
| frag_initial_chromosomes | 64 chromosome entries | Parser._FragBatch.__init__ / size_chroms | Positive initial per-chromosome capacity; retain growth. |
| model_initial_midpoints | 64 midpoint entries | PeakModel._scan_push_mid | Initial reserve; preserve doubling and failure handling. |
| model_initial_peaks | 256 peak entries | PeakModel._scan_close | Initial paired summit/height reserve; preserve common capacity. |
| model_initial_pairs | 1,024 position-pair entries | PeakModel pair collector | Initial reserve; preserve geometric growth. |
| pscore_initial_slots | 65,536 hash slots | CallPeakUnit.PSCORE_CACHE_INIT_SIZE / PScoreCache | Positive power of two; mask derives from actual capacity. |
| pscore_initial_pending_keys | 4,096 pending keys | CallPeakUnit.PScoreCache.reserve_pending | Initial reserve; preserve doubling. |
| pscore_pending_chunk_keys | 8,192 keys per scheduling claim | CallPeakUnit.PENDING_CHUNK / score_pending | Positive chunk; preserve serial threshold and chunk boundaries. |
| qscore_min_slots | 16 hash slots | CallPeakUnit.QScoreTable.__cinit__ | Positive power of two; retain hint-based growth and half-full policy. |
| qtable_prefetch_ahead_slots | 16 upcoming slots | CallPeakUnit.QTABLE_PREFETCH_AHEAD / q-table loops | Independent from qscore_min_slots despite equal value. |
| chromosome_tasks_per_worker | 32 tasks per worker target | CallPeakUnit._chromosome_batches denominator | Name existing multiplier only; positive-target fix belongs to #769. |
| deferred_close_max_files | 256 descriptors per destroy operation | CallPeakUnit.CLOSE_LATER_MAX_FILES / destroy | Preserve current asynchronous cleanup selection and fallback. |
| deferred_close_min_bytes | 1,048,576 bytes per eligible file | CallPeakUnit.CLOSE_LATER_MIN_BYTES / destroy | Independent from serial compressed read capacity. |
| peak_writer_min_bytes | 65,536 UTF-8 bytes | PeakIO._TextBuf.reserve | Initial buffer floor; preserve doubling and chromosome flush behavior. |
| poisson_log_cache_max_entries | 4,194,304 entries per float64 table | Prob.LN_SUM_MAX / acquire_log_tables and scoring wrappers | Positive int32-safe bound; two current tables use up to 64 MiB, with temporary overlap possible. |
| poisson_log_cache_min_entries | 1,024 entries per table | Prob.acquire_log_tables | Preserve current minimum allocation and growth; <= maximum. |


## Native policy and structural bounds

`MACS3/IO/cparser_layout.pxd` defines the following compile-time bounds.
They are not arbitrary Python settings:

| Native definition | Value and units | Consumers and constraints |
| --- | --- | --- |
| FRAG_LOCAL_CHROMOSOMES | 8 chromosomes per chunk | Chunk allocation/stride/bounds and the fixed local ID map use this one enum. Must be positive. |
| FRAG_LOCAL_RUNS | 8 runs per chunk | Chunk allocation, stride and fallback bounds use this separate enum. Must be positive. |
| INTERNER_INLINE_KEY_BYTES | 36 bytes | Both the slot key array and inline/arena comparison bound use this enum; at least 24 bytes for the three-word lookup. |
| INTERNER_ALIGNMENT_BYTES | 64 bytes | Slot-pointer alignment uses its derived mask; must be a power of two. Slot size and allocation slack are checked at import on both 32- and 64-bit builds. |

The interner validates power-of-two capacities before allocation. Slot size
cannot exceed the alignment, and one extra slot must cover the maximum
alignment adjustment given malloc's native alignment. The inline-key bound
is tied to the array declaration instead of a second Python value. The
64-bit slot remains one 64-byte cache line; native type sizes on 32-bit
builds can produce smaller slots.

Allocator-only policy stays in its native owner:

| Definition | Value and units | Owner / rationale |
| --- | --- | --- |
| MACS3_HUGE_ADVICE_MIN_BYTES | 4 MiB of fresh allocation | Signal/cmalloc.pxd: smaller allocations skip huge-page advice. |
| MACS3_HP_MAX_HUGE | 32 MiB per huge page | Signal/chugepages.pxd: maximum supported kernel huge-page size. |
| MACS3_HP_RESERVE_DIVISOR | 4 (extra fraction is 1/4) | Signal/chugepages.pxd: reserve growth headroom. |
| MACS3_HP_RESERVE_PAGES | 1 huge page | Signal/chugepages.pxd: additional headroom before rounding. |

Huge-page reservation still rounds `need + need / 4 + huge_page` and
retains its original mapping fallback, page-size checks and overflow bound.
These native-only policies have no duplicate Python definitions. The IO
layout enums are shared by all their native consumers via the `.pxd`;
new shared native policy belongs in a native definition, not a runtime
Python constant. Source distributions include these `.pxd` files.

## Values that stay with their domain

BAM core length, BGZF block-size limits and framing, gzip window bits,
flags, FRAG grammar/count widths and overflow boundaries, hash multipliers
and byte masks, IEEE masks, integer bounds and the float32 `2^24` exactness
limit are format, layout or arithmetic contracts. In particular, a 16 MiB
buffer and `2^24` exact integer accumulation express unrelated constraints.

Probability convergence and rounding, log constants, model gap/minimum
length, output precision, and `PeakIO.LINE_SLACK` (256 bytes of proven
formatting overhead) retain their domain meaning. Allocator headers,
magic values, overflow bounds, sort implementation details and vendored
libdeflate constants also stay local. Existing scientific/CLI defaults and
legacy `Utilities.Constants` settings are outside this extraction.
