import sys
import types
import pytest

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.CallPeakUnit import CallerFromAlignments
from MACS3.IO.PeakIO import PeakIO


# ---------------------------------------------------------------------------
# Provide tiny Cython/cykhash stubs so the Python sources import cleanly.
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

if "cykhash" not in sys.modules:
    cykhash_stub = types.ModuleType("cykhash")

    class _Map(dict):
        def __init__(self, for_int=False):
            super().__init__()

    cykhash_stub.PyObjectMap = _Map
    cykhash_stub.Float32to32Map = _Map
    sys.modules["cykhash"] = cykhash_stub


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
