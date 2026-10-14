"""Subprocess driver: exercise only the public PeakModel workflow."""

import importlib
import json
import os
from pathlib import Path
import sys
import threading
from types import SimpleNamespace

import numpy as np

from MACS3.Signal.FixWidthTrack import FWTrack


def treatment(kind, two_strands=False):
    track = FWTrack(fw=36)
    if kind in ("ties", "normal_ties"):
        # 69 equal-height maxima in one connected peak grow the midpoint
        # buffer from 64 to 128 (200-bp tag extensions overlap by 1 bp).
        positions = 1000 + np.arange(70) * 199
        if kind == "normal_ties":
            positions = (positions + np.arange(120)[:, None] * 20000).ravel()
    else:
        n = {"last": 1, "between": 2, "growth_last": 257,
             "growth_between": 258, "normal": 258}[kind]
        positions = np.repeat(1000 + np.arange(n) * 1000, 5)
    for pos in positions:
        track.add_loc(b"chr1", int(pos), 0)
        if two_strands:
            track.add_loc(b"chr1", int(pos) + 100, 1)
    track.finalize()
    return track


def build(module, track):
    log = []
    opt = SimpleNamespace(gsize=track.total * 1000, lmfold=0,
                          umfold=200, bw=100, d_min=20,
                          info=log.append, debug=log.append, warn=log.append)
    model = module.PeakModel(opt, track)
    try:
        model.build()
    except module.NotEnoughPairsException:
        # Single-strand fault cases intentionally stop before pairing.
        return {"log": log, "outcome": "not_enough_pairs"}
    attrs = ("plus_line", "minus_line", "shifted_line", "xcorr", "ycorr",
             "d", "scan_window", "alternative_d")
    if not all(hasattr(model, attr) for attr in attrs):
        return {"outcome": "unavailable"}
    result = {attr: (getattr(model, attr).tolist()
                     if isinstance(getattr(model, attr), np.ndarray)
                     else getattr(model, attr)) for attr in attrs}
    return dict(result, log=log, outcome="built")


def main():
    # Disable core files before loading the extension under test.
    import resource
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    directory, module_name, kind, target, worker_only = sys.argv[1:]
    sys.path.insert(0, str(Path(directory).resolve()))
    module = importlib.import_module(module_name)
    target, worker_only = int(target), bool(int(worker_only))
    # Control CPU detection rather than calling the private scheduler.
    os.cpu_count = lambda: 2
    if hasattr(os, "sched_getaffinity"):
        os.sched_getaffinity = lambda pid: {0, 1}
    track = treatment(kind, two_strands=worker_only or kind.startswith("normal"))
    instrumented = hasattr(module, "allocation_reset")
    before = set(threading.enumerate())
    if instrumented:
        module.allocation_reset(target, worker_only)
    if target:
        try:
            build(module, track)
        except MemoryError:
            pass
        else:
            raise AssertionError("PeakModel.build() did not raise MemoryError")
        stats = module.allocation_stats()
        assert stats[1] == 1, stats
        assert stats[2:4] == (0, 0), stats
        assert set(threading.enumerate()) == before
        if worker_only:
            assert stats[6] == 1, stats
        else:
            assert stats[7] == 0, stats  # failed scan made no further allocations
        module.allocation_reset(0, False)
        recovered = build(module, track)
        assert module.allocation_stats()[2:4] == (0, 0)
        # Compare recovery with a fresh, uninjected run of the same input.
        module.allocation_reset(0, False)
        assert recovered == build(module, track)
        assert module.allocation_stats()[2:4] == (0, 0)
        print(json.dumps({"stats": stats, "recovered": recovered["outcome"]}))
    else:
        result = build(module, track)
        if instrumented:
            assert module.allocation_stats()[2:4] == (0, 0)
        print(json.dumps(result))


if __name__ == "__main__":
    main()
