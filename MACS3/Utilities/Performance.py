"""Immutable internal performance defaults, with units in field names.

These are development defaults, not per-job settings. Compiled consumers
copy them into typed values once at import, before entering native loops.
Rebinding module attributes cannot configure an already imported extension.
"""

from dataclasses import dataclass, fields
import sys
from typing import Final


@dataclass(frozen=True, slots=True)
class _PerformanceDefaults:
    # Default invocation budget, including the caller; explicit requests may exceed it.
    cpu_limit: int = 8
    # Small chromosomes do not pay for handing work to another thread.
    min_parallel_track_elements: int = 1 << 16

    # Initial stream buffers can grow for long records; these are not limits
    # on record length. The stream also reads FRAG data.
    bam_serial_input_bytes: int = 1 << 20
    bam_serial_output_bytes: int = 1 << 22
    bgzf_batch_input_bytes: int = 1 << 22
    bgzf_batch_output_bytes: int = 1 << 24
    bgzf_tail_reserve_bytes: int = 1 << 18
    bgzf_max_batch_jobs: int = 8192

    # Append staging, parallel line-walk eligibility and scheduling.
    frag_append_batch_records: int = 1 << 18
    frag_parallel_min_bytes: int = 1 << 16
    frag_tasks_per_thread: int = 2
    # Hash slots must be a power of two; storage grows geometrically.
    frag_interner_initial_slots: int = 1024
    frag_interner_initial_arena_bytes: int = 16384
    frag_initial_chromosomes: int = 64

    # Initial native reserves, doubled when full.
    model_initial_midpoints: int = 64
    model_initial_peaks: int = 256
    model_initial_pairs: int = 1024
    pscore_initial_slots: int = 1 << 16
    pscore_initial_pending_keys: int = 4096
    pscore_pending_chunk_keys: int = 8192
    qscore_min_slots: int = 16
    # A prefetch hint in the q-table insertion loops, independent of capacity.
    qtable_prefetch_ahead_slots: int = 16
    chromosome_tasks_per_worker: int = 32
    # Stored treatment + distinct control records, not FRAG multiplicities.
    # Conservative startup heuristic; not a measured Linux break-even point.
    min_parallel_peak_locations: int = 1 << 18

    # Defer closing at most this many sufficiently large unlinked files.
    deferred_close_max_files: int = 256
    deferred_close_min_bytes: int = 1 << 20
    peak_writer_min_bytes: int = 1 << 16

    # Two float64 tables: up to 64 MiB for the current generation; readers
    # retaining older generations can temporarily increase total memory use.
    poisson_log_cache_max_entries: int = 1 << 22
    poisson_log_cache_min_entries: int = 1024

    def __post_init__(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if type(value) is not int or not 0 < value <= sys.maxsize:
                raise ValueError(f"{field.name} must be a positive native-size integer")
        for name in ("frag_interner_initial_slots", "pscore_initial_slots",
                     "qscore_min_slots"):
            value = getattr(self, name)
            if value & (value - 1):
                raise ValueError(f"{name} must be a power of two")
        for name in ("cpu_limit", "frag_tasks_per_thread",
                     "poisson_log_cache_max_entries"):
            if getattr(self, name) > (1 << 31) - 1:
                raise ValueError(f"{name} must fit a signed 32-bit integer")
        if self.poisson_log_cache_min_entries > self.poisson_log_cache_max_entries:
            raise ValueError("Poisson log-cache minimum exceeds its maximum")


DEFAULTS: Final = _PerformanceDefaults()
