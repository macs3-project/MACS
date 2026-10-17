# Threads.py: run per-chromosome work of a track on threads
"""Run the independent per-chromosome steps of a track on threads.

FWTrack, PETrackI and PETrackII run finalize and filter_dup chromosome by
chromosome. A chromosome of fewer than MIN_THREADED_SIZE elements they
run on the spot; the others they hand to map_chromosomes, which runs them
on threads. That work releases the GIL in its C loops and in the numpy
calls it makes (sort, take, boolean indexing, sum), so chromosomes run in
parallel. Every thread is joined before map_chromosomes returns.
"""
import threading
from MACS3.Utilities.Performance import DEFAULTS
from MACS3.Utilities.Execution import resolve_execution_policy

# Informational compatibility alias; jobs use an explicit ExecutionPolicy.
MAX_THREADS = DEFAULTS.cpu_limit

# The track methods run a chromosome with fewer elements than this on the
# calling thread rather than through map_chromosomes. Its work is too short
# to pay for handing the GIL between threads, and a track of many small
# contigs then runs as it did before.
MIN_THREADED_SIZE = DEFAULTS.min_parallel_track_elements


def n_threads(execution_policy=None) -> int:
    """Return the invocation's total thread budget (caller included)."""
    return resolve_execution_policy(execution_policy=execution_policy).threads


def run_on_threads(work, n_workers, on_error, args=()):
    """Run work on the caller and n_workers additional threads.

    Report startup/caller failures through on_error before joining every
    successfully started worker, then propagate the original exception.
    on_error must stop further task dispatch and must not raise.
    """
    started = []
    try:
        for _ in range(n_workers):
            thread = threading.Thread(target=work, args=args)
            thread.start()
            started.append(thread)
        work(*args)
    except BaseException as error:
        on_error(error)
        raise
    finally:
        for thread in started:
            thread.join()


def map_chromosomes(fn, chroms: list, sizes: list, *args,
                    execution_policy=None) -> list:
    """Return [fn(c, *args) for c in chroms], calling fn once per
    chromosome.

    sizes[i] is the number of elements of chroms[i]. With two or more
    chromosomes and a budget above 1, the chromosomes run largest
    first on min(budget, len(chroms)) threads, the calling one
    included; otherwise they run on the calling thread in the order
    given. The calls must not depend on each other. The first exception
    raised stops the calls not yet started and is re-raised once every
    thread has stopped.
    """
    results = [None] * len(chroms)
    budget = n_threads(execution_policy)
    n = min(budget, len(chroms)) if len(chroms) >= 2 else 1
    if n < 2:
        for i in range(len(chroms)):
            results[i] = fn(chroms[i], *args)
        return results

    lock = threading.Lock()
    queue = iter(sorted(range(len(chroms)), key=lambda i: -sizes[i]))
    errors = []

    def work():
        while True:
            with lock:
                if errors:
                    return
                i = next(queue, None)
            if i is None:
                return
            try:
                results[i] = fn(chroms[i], *args)
            except BaseException as e:
                with lock:
                    errors.append(e)
                return

    def stop(error):
        with lock:
            errors.append(error)

    run_on_threads(work, n - 1, stop)
    if errors:
        raise errors[0]
    return results
