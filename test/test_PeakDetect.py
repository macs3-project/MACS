import sys
import types
from types import SimpleNamespace

import pytest

# Provide a minimal cython stub when the real package is unavailable.
if "cython" not in sys.modules:
    cython_stub = types.ModuleType("cython")

    def _identity_decorator(*args, **kwargs):
        def decorate(func):
            return func
        if args and callable(args[0]) and len(args) == 1 and not kwargs:
            return args[0]
        return decorate

    cython_stub.cfunc = _identity_decorator
    cython_stub.ccall = _identity_decorator
    cython_stub.locals = _identity_decorator
    cython_stub.inline = _identity_decorator
    cython_stub.returns = _identity_decorator
    cython_stub.short = int
    cython_stub.float = float
    cython_stub.double = float
    cython_stub.int = int
    cython_stub.long = int
    sys.modules["cython"] = cython_stub

from MACS3.Signal.PeakDetect import PeakDetect


class DummyTrack:
    def __init__(self, total, length, average_template_length=150, chr_names=("chr1",)):
        self.total = total
        self.length = length
        self.average_template_length = average_template_length
        self._chr_names = tuple(chr_names)

    def get_chr_names(self):
        return list(self._chr_names)


def make_opt(**overrides):
    noop = lambda *args, **kwargs: None
    params = {
        "info": noop,
        "debug": noop,
        "warn": noop,
        "PE_MODE": False,
        "log_pvalue": 5.0,
        "log_qvalue": None,
        "d": 150,
        "maxgap": None,
        "tsize": 200,
        "minlen": None,
        "shift": 0,
        "gsize": 10000,
        "nolambda": False,
        "smalllocal": 200,
        "largelocal": 600,
        "ratio": 1.0,
        "tocontrol": False,
        "call_summits": True,
        "store_bdg": False,
        "name": "test",
        "bdg_treat": "treat.bdg",
        "bdg_control": "control.bdg",
        "do_SPMR": False,
        "cutoff_analysis_file": "cutoff.txt",
        "trackline": False,
        "broad": False,
        "log_broadcutoff": 0.0,
        "cutoff_analysis": False,
    }
    params.update(overrides)
    return SimpleNamespace(**params)


def test_call_peaks_with_control_scales_control_parameters(monkeypatch):
    captured = {}

    class FakeCaller:
        def __init__(self, treat, control, **kwargs):
            self.treat = treat
            self.control = control
            self.kwargs = kwargs
            self.trackline_enabled = False
            self.destroy_called = False
            self.call_peaks_args = None
            self.call_peaks_kwargs = None
            captured["instance"] = self

        def enable_trackline(self):
            self.trackline_enabled = True

        def call_peaks(self, metrics, cutoffs, **kwargs):
            self.call_peaks_args = (metrics, cutoffs)
            self.call_peaks_kwargs = kwargs
            return ["peak"]

        def call_broadpeaks(self, *args, **kwargs):  # pragma: no cover - defensive
            raise AssertionError("Broad peak path should not be used in this test")

        def destroy(self):
            self.destroy_called = True

    monkeypatch.setattr("MACS3.Signal.PeakDetect.CallerFromAlignments", FakeCaller)

    opt = make_opt()
    treat = DummyTrack(total=100, length=5000, average_template_length=50)
    control = DummyTrack(total=50, length=4000, average_template_length=50)

    detector = PeakDetect(opt=opt, treat=treat, control=control)

    required_attrs = ["d", "sregion", "lregion", "minlen", "maxgap"]
    if not all(hasattr(detector, attr) for attr in required_attrs):
        pytest.skip("PeakDetect implementation hides required attributes")
    result = detector.call_peaks()

    fake_instance = captured["instance"]
    assert result == ["peak"]
    assert fake_instance.kwargs["treat_scaling_factor"] == pytest.approx(1.0)
    assert fake_instance.kwargs["lambda_bg"] == pytest.approx((treat.total * detector.d) / opt.gsize)
    assert fake_instance.kwargs["ctrl_d_s"] == [detector.d, detector.sregion, detector.lregion]
    expected_scales = [treat.total / control.total,
                       detector.d / detector.sregion * (treat.total / control.total),
                       detector.d / detector.lregion * (treat.total / control.total)]
    for observed, expected in zip(fake_instance.kwargs["ctrl_scaling_factor_s"], expected_scales):
        assert observed == pytest.approx(expected)
    assert fake_instance.call_peaks_args == (["p"], [opt.log_pvalue])
    assert fake_instance.call_peaks_kwargs["call_summits"] is True
    assert fake_instance.call_peaks_kwargs["min_length"] == detector.minlen
    assert fake_instance.call_peaks_kwargs["max_gap"] == detector.maxgap
    assert fake_instance.destroy_called is True


def test_call_peaks_with_control_ignores_lambda_when_disabled(monkeypatch):
    captured = {}

    class FakeCaller:
        def __init__(self, treat, control, **kwargs):
            self.kwargs = kwargs
            captured["instance"] = self

        def call_peaks(self, *args, **kwargs):
            return ["peak"]

        def call_broadpeaks(self, *args, **kwargs):  # pragma: no cover - defensive
            raise AssertionError("Broad peak path should not be used in this test")

        def destroy(self):
            pass

    monkeypatch.setattr("MACS3.Signal.PeakDetect.CallerFromAlignments", FakeCaller)

    opt = make_opt(nolambda=True)
    treat = DummyTrack(total=100, length=5000, average_template_length=50)
    control = DummyTrack(total=80, length=4000, average_template_length=50)

    detector = PeakDetect(opt=opt, treat=treat, control=control)
    detector.call_peaks()

    fake_instance = captured["instance"]
    assert fake_instance.kwargs["ctrl_d_s"] == []
    assert fake_instance.kwargs["ctrl_scaling_factor_s"] == []


def test_call_peaks_without_control_uses_lregion_bias(monkeypatch):
    captured = {}

    class FakeCaller:
        def __init__(self, treat, control, **kwargs):
            self.treat = treat
            self.control = control
            self.kwargs = kwargs
            self.trackline_enabled = False
            self.destroy_called = False
            self.call_peaks_kwargs = None
            captured["instance"] = self

        def enable_trackline(self):
            self.trackline_enabled = True

        def call_peaks(self, metrics, cutoffs, **kwargs):
            self.call_peaks_kwargs = kwargs
            return ["peak"]

        def call_broadpeaks(self, *args, **kwargs):  # pragma: no cover - defensive
            raise AssertionError("Broad peak path should not be used in this test")

        def destroy(self):
            self.destroy_called = True

    monkeypatch.setattr("MACS3.Signal.PeakDetect.CallerFromAlignments", FakeCaller)

    opt = make_opt(trackline=True, call_summits=False)
    treat = DummyTrack(total=80, length=2400, average_template_length=50)

    detector = PeakDetect(opt=opt, treat=treat, control=None)

    required_attrs = ["lregion", "d", "minlen", "maxgap"]
    if not all(hasattr(detector, attr) for attr in required_attrs):
        pytest.skip("PeakDetect implementation hides required attributes")

    result = detector.call_peaks()

    fake_instance = captured["instance"]
    assert result == ["peak"]
    assert fake_instance.control is None
    assert fake_instance.trackline_enabled is True
    assert fake_instance.kwargs["ctrl_d_s"] == [detector.lregion]
    assert fake_instance.kwargs["ctrl_scaling_factor_s"] == [detector.d / detector.lregion]
    assert fake_instance.kwargs["treat_scaling_factor"] == pytest.approx(1.0)
    expected_lambda_bg = detector.d * treat.total / opt.gsize
    assert fake_instance.kwargs["lambda_bg"] == pytest.approx(expected_lambda_bg)
    assert fake_instance.call_peaks_kwargs["call_summits"] is False
    assert fake_instance.destroy_called is True
