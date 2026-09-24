import sys
import tempfile
import types
import pytest
import numpy as np

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.CallPeakUnit import CallerFromAlignments
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.IO.PeakIO import PeakIO


# ---------------------------------------------------------------------------
# Provide tiny Cython stubs so the Python sources import cleanly.
# ---------------------------------------------------------------------------
cython_stub = sys.modules.get("cython")
if cython_stub is None:
    cython_stub = types.ModuleType("cython")
    sys.modules["cython"] = cython_stub


def _identity_decorator(*args, **kwargs):
    def decorate(func):
        return func
    if args and callable(args[0]) and len(args) == 1 and not kwargs:
        return args[0]
    return decorate


for name in ("cfunc", "ccall", "cclass", "locals", "inline", "returns",
             "boundscheck", "wraparound"):
    setattr(cython_stub, name, getattr(cython_stub, name, _identity_decorator))

for attr, default in [
    ("declare", lambda *a, **k: None),
    ("short", int),
    ("float", float),
    ("double", float),
    ("int", int),
    ("long", int),
    ("ulong", int),
    ("bint", bool),
]:
    setattr(cython_stub, attr, getattr(cython_stub, attr, default))

cimports_mod = sys.modules.get("cython.cimports")
if cimports_mod is None:
    cimports_mod = types.ModuleType("cython.cimports")
    sys.modules["cython.cimports"] = cimports_mod
cython_stub.cimports = cimports_mod

if "cython.cimports.cpython" not in sys.modules:
    cpython_mod = types.ModuleType("cython.cimports.cpython")
    cpython_mod.bool = bool
    sys.modules["cython.cimports.cpython"] = cpython_mod
    cimports_mod.cpython = cpython_mod

if "cython.cimports.numpy" not in sys.modules:
    numpy_mod = types.ModuleType("cython.cimports.numpy")
    numpy_mod.ndarray = lambda *args, **kwargs: None
    sys.modules["cython.cimports.numpy"] = numpy_mod
    cimports_mod.numpy = numpy_mod

def make_fwtrack(layout, fw=50):
    track = FWTrack(fw=fw)
    for chrom, positions in layout.items():
        for pos in positions:
            track.add_loc(chrom, pos, strand=0)
    track.finalize()
    return track


def make_tracks():
    treat = make_fwtrack({b"chr1": [10, 30]})
    ctrl = make_fwtrack({b"chr1": [15], b"chr2": [5]})
    return treat, ctrl


def make_two_summit_profile(length=260, right_apex_distance=15):
    """Return a peak profile with a summit close to the right boundary."""
    signal = np.full(length, 10.0)
    half_width = 70
    for center in (60, length - right_apex_distance):
        positions = np.arange(max(center - half_width + 1, 0),
                              min(center + half_width, length))
        signal[positions] = np.maximum(
            signal[positions],
            10 + 30 * (1 - np.abs(positions - center) / half_width),
        )
    return np.rint(signal).astype(int)


def make_petrack_from_profile(signal, start):
    """Encode an integer coverage profile as a paired-end track."""
    track = PETrackI(buffer_size=200000)
    for level in range(1, int(signal.max()) + 1):
        edges = np.diff(np.r_[False, signal >= level, False].astype(int))
        for left, right in zip(np.flatnonzero(edges == 1),
                               np.flatnonzero(edges == -1)):
            track.add_loc(b"chrSynthetic", start + int(left),
                          start + int(right))
    track.finalize()
    return track


def test_constructor_rejects_unknown_track_type():
    with pytest.raises(Exception):
        CallerFromAlignments(object(), object())


def test_destroy_accepts_empty_state(tmp_path):
    calc = CallerFromAlignments(*make_tracks())
    if hasattr(calc, "pileup_data_files"):
        dummy = tmp_path / "pile.tmp"
        dummy.write_text("tmp")
        calc.pileup_data_files = {b"chr1": str(dummy)}
    calc.destroy()  # should not raise even if files missing or attribute hidden


def test_set_pseudocount_no_exception():
    calc = CallerFromAlignments(*make_tracks())
    calc.set_pseudocount(5.0)
    if hasattr(calc, "pseudocount"):
        assert calc.pseudocount == pytest.approx(5.0)


def test_enable_trackline_idempotent():
    calc = CallerFromAlignments(*make_tracks())
    calc.enable_trackline()
    calc.enable_trackline()
    if hasattr(calc, "trackline"):
        assert calc.trackline is True


def test_call_peaks_with_no_chromosomes_returns_peakio():
    treat, ctrl = make_tracks()
    calc = CallerFromAlignments(treat, ctrl)
    if hasattr(calc, "chromosomes"):
        calc.chromosomes = []
    if not hasattr(calc, "pqtable"):
        pytest.skip("pqtable not exposed in current build")
    calc.pqtable[0.0] = 0.0
    peaks = calc.call_peaks(['p'], [1.0], min_length=10, max_gap=5, call_summits=False, cutoff_analysis=False)
    assert isinstance(peaks, PeakIO)


def test_call_summits_rejects_maximum_in_below_cutoff_gap(tmp_path,
                                                           monkeypatch):
    """A smoothed summit must map to an above-cutoff signal chunk."""
    signal = np.rint(np.interp(np.arange(123),
                               np.linspace(0, 122, 8),
                               [27, 1, 27, 8, 14, 10, 27, 12])).astype(int)
    track = PETrackI(buffer_size=10000)

    for level in range(1, int(signal.max()) + 1):
        edges = np.diff(np.r_[False, signal >= level, False].astype(int))
        starts = np.flatnonzero(edges == 1)
        ends = np.flatnonzero(edges == -1)
        for start, end in zip(starts, ends):
            track.add_loc(b"chrSynthetic", 1000 + int(start),
                          1000 + int(end))
    track.finalize()

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    caller = CallerFromAlignments(track, None, ctrl_d_s=[],
                                  ctrl_scaling_factor_s=[], lambda_bg=1.0)
    try:
        peaks = caller.call_peaks(["f"], [9.0], min_length=75,
                                  max_gap=100, call_summits=True)
    finally:
        caller.destroy()

    peak = peaks.peaks[b"chrSynthetic"][0]
    actual_pileup = signal[peak["summit"] - 1000]
    assert peak["summit"] == 1035
    assert peak["pileup"] == actual_pileup == 27
    assert peak["fc"] == pytest.approx((actual_pileup + 1) / 2)


@pytest.mark.parametrize("start", [0, 1, 5, 9, 10, 1000])
def test_call_summits_keeps_right_edge_candidate(tmp_path, monkeypatch,
                                                  start):
    """Regression test for the padding-coordinate mismatch in issue #747."""
    signal = make_two_summit_profile()
    track = make_petrack_from_profile(signal, start)

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    caller = CallerFromAlignments(track, None, ctrl_d_s=[],
                                  ctrl_scaling_factor_s=[], lambda_bg=1.0)
    try:
        peaks = caller.call_peaks(["f"], [2.0], min_length=50,
                                  max_gap=100, call_summits=True)
    finally:
        caller.destroy()

    rows = peaks.peaks[b"chrSynthetic"]
    relative_summits = [row["summit"] - start for row in rows]
    assert len(relative_summits) == 2
    assert abs(relative_summits[0] - 60) <= 1
    assert abs(relative_summits[1] - 236) <= 1
    assert all(row["start"] == start for row in rows)
    assert all(row["end"] == start + len(signal) for row in rows)
    assert all(row["start"] <= row["summit"] < row["end"] for row in rows)
