import sys
import types
from types import SimpleNamespace

import numpy as np
import pytest


# Ensure a minimal cython stub is present so the module can import.
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


cython_stub.cfunc = getattr(cython_stub, "cfunc", _identity_decorator)
cython_stub.ccall = getattr(cython_stub, "ccall", _identity_decorator)
cython_stub.cclass = getattr(cython_stub, "cclass", _identity_decorator)
cython_stub.locals = getattr(cython_stub, "locals", _identity_decorator)
cython_stub.inline = getattr(cython_stub, "inline", _identity_decorator)
cython_stub.returns = getattr(cython_stub, "returns", _identity_decorator)
cython_stub.declare = getattr(cython_stub, "declare", lambda *a, **k: None)
cython_stub.short = getattr(cython_stub, "short", int)
cython_stub.float = getattr(cython_stub, "float", float)
cython_stub.double = getattr(cython_stub, "double", float)
cython_stub.int = getattr(cython_stub, "int", int)
cython_stub.long = getattr(cython_stub, "long", int)
cython_stub.bint = getattr(cython_stub, "bint", bool)

cimports_mod = sys.modules.get("cython.cimports")
if cimports_mod is None:
    cimports_mod = types.ModuleType("cython.cimports")
    sys.modules["cython.cimports"] = cimports_mod
cython_stub.cimports = cimports_mod

cpython_mod = sys.modules.get("cython.cimports.cpython")
if cpython_mod is None:
    cpython_mod = types.ModuleType("cython.cimports.cpython")
    cpython_mod.bool = bool
    sys.modules["cython.cimports.cpython"] = cpython_mod
cimports_mod.cpython = cpython_mod

numpy_mod = sys.modules.get("cython.cimports.numpy")
if numpy_mod is None:
    numpy_mod = types.ModuleType("cython.cimports.numpy")
    numpy_mod.ndarray = np.ndarray
    sys.modules["cython.cimports.numpy"] = numpy_mod
cimports_mod.numpy = numpy_mod

from MACS3.Signal.PeakModel import PeakModel, NotEnoughPairsException


class DummyTreatment:
    def __init__(self, total=0, length=0, chromosomes=None, chr_locations=None):
        self.total = total
        self.length = length
        self._chromosomes = chromosomes or [b"chr1"]
        self._chr_locations = chr_locations or {
            self._chromosomes[0]: (np.array([], dtype="int32"), np.array([], dtype="int32"))
        }

    def get_chr_names(self):
        return list(self._chromosomes)

    def get_locations_by_chr(self, chrom):
        return self._chr_locations.get(chrom, (np.array([], dtype="int32"), np.array([], dtype="int32")))


def make_opt(**overrides):
    capture = overrides.pop("warn_capture", None)
    noop = lambda *args, **kwargs: None
    warn_fn = (capture.append if capture is not None else noop)
    params = {
        "gsize": 1000000,
        "umfold": 30,
        "lmfold": 10,
        "d_min": 20,
        "bw": 100,
        "info": noop,
        "debug": noop,
        "warn": warn_fn,
    }
    params.update(overrides)
    return SimpleNamespace(**params)


def test_build_raises_when_not_enough_pairs():
    warnings = []
    opt = make_opt(warn_capture=warnings)
    treatment = DummyTreatment(total=1000)
    model = PeakModel(opt, treatment)

    with pytest.raises(NotEnoughPairsException):
        model.build()

    assert warnings, "Expected build to emit warning messages"


def test_build_computes_thresholds_before_pairing():
    warnings = []
    opt = make_opt(warn_capture=warnings, lmfold=5, umfold=15, bw=120, gsize=500000)
    treatment = DummyTreatment(total=300)
    model = PeakModel(opt, treatment)

    if not hasattr(model, "peaksize"):
        pytest.skip("PeakModel attributes unavailable in current build")

    with pytest.raises(NotEnoughPairsException):
        model.build()

    expected_peaksize = 2 * opt.bw
    assert model.peaksize == expected_peaksize
    expected_min = int(round(float(treatment.total) * opt.lmfold * expected_peaksize / opt.gsize / 2))
    expected_max = int(round(float(treatment.total) * opt.umfold * expected_peaksize / opt.gsize / 2))
    assert model.min_tags == expected_min
    assert model.max_tags == expected_max


def test_str_representation_reflects_summary_fields():
    opt = make_opt()
    treatment = DummyTreatment(total=100)
    model = PeakModel(opt, treatment)

    if not all(hasattr(model, attr) for attr in ("min_tags", "max_tags", "d", "scan_window")):
        pytest.skip("PeakModel attributes unavailable in current build")

    model.min_tags = 3
    model.max_tags = 9
    model.d = 147
    model.scan_window = 300

    summary = str(model)

    assert "Baseline: 3" in summary
    assert "Upperline: 9" in summary
    assert "Fragment size: 147" in summary
    assert "Scan window size: 300" in summary
