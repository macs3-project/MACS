"""Internal, invocation-owned peak-caller scheduling and pool lifetime.

No scoring types are imported here. The caller passes typed worker adapters
and the pool initializer; only children bind a caller in that initializer.
Each PeakExecution owns a workload snapshot and an immutable ExecutionPolicy,
but never retains the caller or live pool handles between operations.

The shared lock is synchronization only: as before, it serializes fork pools
through shutdown so a new pool cannot inherit another pool's descriptors.
"""

from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import signal
import threading

from MACS3.Utilities.Execution import (resolve_execution_policy,
                                      process_backend_reason)
from MACS3.Utilities.Logger import logging
from MACS3.Utilities.Performance import DEFAULTS


execution_info = logging.getLogger("MACS3.Utilities.Execution").info
_worker_pool_lock = threading.Lock()


def peak_processes_allowed(process_backend: str) -> bool:
    """Whether the current host supports the automatic Linux fork policy.

    Python thread detection cannot detect native library threads or prevent
    a host from starting a thread after this check. Embedded callers should
    select ``serial`` when they cannot establish fork safety themselves.
    """
    return process_backend_reason(process_backend) is None


def n_available_workers(execution_policy=None) -> int:
    """Allowed child processes, reserving one budget slot for the parent."""
    policy = resolve_execution_policy(execution_policy=execution_policy)
    return (policy.process_workers
            if peak_processes_allowed(policy.process_backend) else 0)


def n_locations(track, chrom: bytes) -> int:
    """Return the number of reads or fragments of a track on a chromosome."""
    locs = track.get_locations_by_chr(chrom)
    if isinstance(locs, (list, tuple)):
        return sum([x.shape[0] for x in locs])
    return locs.shape[0]


class _PeakProcessContext:
    """Own process handles even if executor startup stops before its manager.

    ProcessPoolExecutor.shutdown cannot join a partially launched fork pool
    before its manager thread exists. Record handles at the context boundary
    so failure cleanup does not depend on executor implementation attributes.
    """

    def __init__(self):
        self.context = multiprocessing.get_context('fork')
        self.processes = []

    def __getattr__(self, name):
        return getattr(self.context, name)

    def Process(self, *args, **kwargs):
        process = self.context.Process(*args, **kwargs)
        self.processes.append(process)
        return process


def shutdown_peak_pool(pool, processes=()):
    """Defer SIGINT until joining finishes, then restore the host's handler.

    On Python 3.12 an interrupted Thread.join can mark a live executor manager
    thread stopped. Retrying shutdown then returns before its workers finish.
    Record SIGINT without raising inside join, and deliver it after shutdown.
    """
    interrupted = None
    pending = []
    handler = signal.getsignal(signal.SIGINT)
    defer_signal = (threading.current_thread() is threading.main_thread() and
                    callable(handler))

    def defer(signum, frame):
        pending.append((signum, frame))

    if defer_signal:
        signal.signal(signal.SIGINT, defer)
    try:
        try:
            while True:
                try:
                    pool.shutdown(wait=True, cancel_futures=True)
                    break
                except KeyboardInterrupt as exc:
                    interrupted = exc
        finally:
            for process in processes:
                if process.pid is not None:
                    if process.is_alive():
                        process.terminate()
                    process.join()
    finally:
        if defer_signal:
            signal.signal(signal.SIGINT, handler)
    if pending:
        handler(*pending[-1])
    if interrupted is not None:
        raise interrupted


class PeakExecution:
    """One caller's workload and synchronous execution scopes.

    Call prepare once per public call. Work indices retain serial chromosome
    order for histogram ties; batches are submitted largest first. Safety is
    rechecked under the fork lock immediately before pool creation. A scope
    yields None for serial fallback, or an executor joined before scope exit,
    including submission/result failure, interrupted join, and partial startup.
    """

    def __init__(self, execution_policy):
        self.policy = execution_policy
        self.chromosome_indices = []
        self.active_chromosomes = []
        self.chromosome_work = []
        self.chromosome_work_total = 0

    def prepare(self, treat, ctrl, chromosomes):
        """Snapshot stored-location counts without changing serial indices."""

        self.chromosome_indices = []
        self.active_chromosomes = []
        self.chromosome_work = []
        self.chromosome_work_total = 0
        for i in range(len(chromosomes)):
            w = n_locations(treat, chromosomes[i])
            if ctrl is not treat:
                w += n_locations(ctrl, chromosomes[i])
            if w == 0:
                continue
            self.chromosome_indices.append(i)
            self.active_chromosomes.append(chromosomes[i])
            self.chromosome_work.append((-w, i))
            self.chromosome_work_total += w
        self.chromosome_work.sort()

    def workers(self):
        """Apply the same budget, host and minimum-work gate in both phases."""
        if not peak_processes_allowed(self.policy.process_backend):
            execution_info("Peak process backend: serial (%s)",
                           process_backend_reason(self.policy.process_backend))
            return 0
        if (len(self.chromosome_indices) < 2 or
                self.chromosome_work_total < DEFAULTS.min_parallel_peak_locations):
            execution_info(
                "Peak process backend: serial (workload: %d locations, %d nonempty chromosomes; minimum %d locations)",
                self.chromosome_work_total, len(self.chromosome_indices),
                DEFAULTS.min_parallel_peak_locations)
            return 0
        return min(self.policy.process_workers, len(self.chromosome_indices))

    def batches(self, n_workers):
        """Largest first, at most tasks-per-worker * workers nonempty tasks.

        Every task except a possible trailing batch weighs at least the
        positive ceiling target. Zero-weight chromosomes never enter tasks.
        """
        batches: list = []
        batch: list = []

        task_limit = DEFAULTS.chromosome_tasks_per_worker * n_workers
        target = max(1, self.chromosome_work_total // task_limit +
                     (self.chromosome_work_total % task_limit != 0))
        batch_weight = 0
        for (w, i) in self.chromosome_work:
            batch.append(i)
            batch_weight -= w
            if batch_weight >= target:
                batches.append(batch)
                batch = []
                batch_weight = 0
        if batch:
            batches.append(batch)
        return batches

    @contextmanager
    def pool(self, caller, n_workers, n_tasks, initializer):
        """Own pool handles locally; release the lock only after every join."""
        with _worker_pool_lock:
            if not peak_processes_allowed(self.policy.process_backend):
                execution_info(
                    "Peak process backend: serial after safety recheck (%s)",
                    process_backend_reason(self.policy.process_backend))
                yield None
                return
            workers = min(n_workers, n_tasks)
            execution_info("Peak process backend: fork, workers=%d", workers)
            context = _PeakProcessContext()
            pool = ProcessPoolExecutor(max_workers=workers, mp_context=context,
                                       initializer=initializer, initargs=(caller,))
            try:
                yield pool
            finally:
                shutdown_peak_pool(pool, context.processes)

    def results(self, pool, worker, tasks):
        """Submit a phase and yield results in completion order inside pool()."""
        futures = [pool.submit(worker, task) for task in tasks]
        for future in as_completed(futures):
            yield future.result()
