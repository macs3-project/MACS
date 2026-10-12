# Threads.py: run per-chromosome work of a track on threads
"""Run the independent per-chromosome steps of a track on threads.

FWTrack, PETrackI and PETrackII run finalize and filter_dup chromosome by
chromosome. A chromosome of fewer than MIN_THREADED_SIZE elements they
run on the spot; the others they hand to map_chromosomes, which runs them
on threads. That work releases the GIL in its C loops and in the numpy
calls it makes (sort, take, boolean indexing, sum), so chromosomes run in
parallel. Every thread is joined before map_chromosomes returns.
"""
import os
import threading

# at most this many threads, the calling one included
MAX_THREADS = 8

# The track methods run a chromosome with fewer elements than this on the
# calling thread rather than through map_chromosomes. Its work is too short
# to pay for handing the GIL between threads, and a track of many small
# contigs then runs as it did before.
MIN_THREADED_SIZE = 65536


def n_threads() -> int:
    """min(MAX_THREADS, the cores this process may run on)."""
    if hasattr(os, "sched_getaffinity"):
        return min(MAX_THREADS, len(os.sched_getaffinity(0)))
    return min(MAX_THREADS, os.cpu_count() or 1)


def map_chromosomes(fn, chroms: list, sizes: list, *args) -> list:
    """Return [fn(c, *args) for c in chroms], calling fn once per
    chromosome.

    sizes[i] is the number of elements of chroms[i]. With two or more
    chromosomes and n_threads() above 1, the chromosomes run largest
    first on min(n_threads(), len(chroms)) threads, the calling one
    included; otherwise they run on the calling thread in the order
    given. The calls must not depend on each other. The first exception
    raised stops the calls not yet started and is re-raised once every
    thread has stopped.
    """
    results = [None] * len(chroms)
    n = min(n_threads(), len(chroms)) if len(chroms) >= 2 else 1
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

    threads = [threading.Thread(target=work) for _ in range(n - 1)]
    for t in threads:
        t.start()
    try:
        work()
    finally:
        for t in threads:
            t.join()
    if errors:
        raise errors[0]
    return results
