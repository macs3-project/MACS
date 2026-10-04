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


# ===========================================================================
# Comprehensive tests for PeakDetect (added below the original tests,
# which are kept unchanged).
#
# PeakDetect turns an options namespace into the arguments of
# CallerFromAlignments and picks the calling method.  The parameter tests
# replace CallerFromAlignments with a recorder and compare its arguments
# with values computed from these definitions:
#   SE with control: treat_sum = treat.total * d, control_sum =
#       control.total * d
#   PE with control: treat_sum = treat.length, control_sum =
#       2 * control.total * treat.average_template_length truncated to an
#       integer (PeakDetect keeps both sums in C long variables), and the caller's
#       d is the average template length
#   ratio = treat_sum / control_sum unless --ratio is not 1
#   control scaled to treatment (tocontrol False): lambda_bg = treat_sum /
#       gsize, treatment scale 1, control windows (d, slocal, llocal) with
#       scales ratio * (1, d/slocal, d/llocal)
#   treatment scaled to control (tocontrol True): lambda_bg = control_sum
#       / gsize, treatment scale 1/ratio, scales (1, d/slocal, d/llocal)
#   slocal is used when non-zero, llocal when larger than slocal
#   SE without control: lambda_bg = d * total / gsize, one llocal window
#       with scale d / llocal
#   PE without control: caller d 0, lambda_bg = length / gsize, one llocal
#       window with scale length / (llocal * total * 2)
#   nolambda: no control windows at all
# The end-to-end tests run the real caller on small tracks and check the
# written pileup/lambda bedGraphs against a per-base numpy reference.
# ===========================================================================

import logging

import numpy as np
from scipy.stats import poisson

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI
from MACS3.Signal.CallPeakUnit import CallerFromAlignments
from MACS3.IO.PeakIO import (PeakIO,
                             BroadPeakIO)

F32 = np.float32


class Messages:
    """Collects what PeakDetect sends to opt.info / opt.debug / opt.warn."""

    def __init__(self):
        self.info = []
        self.debug = []
        self.warn = []


def opt_and_log(**overrides):
    log = Messages()
    opt = make_opt(info=log.info.append, debug=log.debug.append,
                   warn=log.warn.append, **overrides)
    return opt, log


class RecordingCaller:
    """Stands in for CallerFromAlignments and records how it is used."""
    created = []

    def __init__(self, treat, control, **kwargs):
        self.treat = treat
        self.control = control
        self.kwargs = kwargs
        self.calls = []
        RecordingCaller.created.append(self)

    def enable_trackline(self):
        self.calls.append(("enable_trackline",))

    def call_peaks(self, syms, cutoffs, **kwargs):
        self.calls.append(("call_peaks", syms, cutoffs, kwargs))
        return "narrow peaks"

    def call_broadpeaks(self, syms, **kwargs):
        self.calls.append(("call_broadpeaks", syms, kwargs))
        return "broad peaks"

    def destroy(self):
        self.calls.append(("destroy",))


@pytest.fixture
def recorder(monkeypatch):
    RecordingCaller.created = []
    monkeypatch.setattr("MACS3.Signal.PeakDetect.CallerFromAlignments",
                        RecordingCaller)
    return RecordingCaller.created


def expected_caller_args(opt, treat, control):
    """Arguments of CallerFromAlignments from the definitions above."""
    d = opt.d
    s, l = opt.smalllocal, opt.largelocal
    if control is not None:
        if opt.PE_MODE:
            d_arg = treat.average_template_length
            treat_sum = treat.length
            control_sum = int(control.total * 2 * treat.average_template_length)
        else:
            d_arg = d
            treat_sum = treat.total * d
            control_sum = control.total * d
        ratio = treat_sum / control_sum if opt.ratio == 1.0 else opt.ratio
        if opt.tocontrol:
            lambda_bg, treat_scale, base = control_sum / opt.gsize, 1 / ratio, 1.0
        else:
            lambda_bg, treat_scale, base = treat_sum / opt.gsize, 1.0, ratio
        ds, scales = [d], [base]
        if s:
            ds.append(s)
            scales.append(d / s * base)
        if l and l > s:
            ds.append(l)
            scales.append(d / l * base)
    else:
        treat_scale = 1.0
        if opt.PE_MODE:
            d_arg = 0
            lambda_bg = treat.length / opt.gsize
            scales = [treat.length / (l * treat.total * 2)]
        else:
            d_arg = d
            lambda_bg = d * treat.total / opt.gsize
            scales = [d / l]
        ds = [l]
    if opt.nolambda:
        ds, scales = [], []
    return dict(d=d_arg, ctrl_d_s=ds, ctrl_scaling_factor_s=scales,
                treat_scaling_factor=treat_scale, lambda_bg=lambda_bg,
                end_shift=opt.shift, save_bedGraph=opt.store_bdg,
                bedGraph_filename_prefix=opt.name,
                bedGraph_treat_filename=opt.bdg_treat,
                bedGraph_control_filename=opt.bdg_control,
                save_SPMR=opt.do_SPMR,
                cutoff_analysis_filename=opt.cutoff_analysis_file)


def assert_caller_args(got, exp):
    """lambda_bg, the treatment scale and d pass through float32 locals in
    PeakDetect, hence rel=1e-6 for those."""
    assert sorted(got) == sorted(exp)
    for k, v in exp.items():
        if k in ("lambda_bg", "treat_scaling_factor", "d"):
            assert got[k] == pytest.approx(v, rel=1e-6), k
        elif k == "ctrl_scaling_factor_s":
            assert got[k] == pytest.approx(v, rel=1e-12), k
        else:
            assert got[k] == v, k


def se_dummy():
    return (DummyTrack(total=100, length=5000, average_template_length=50),
            DummyTrack(total=40, length=2000, average_template_length=50))


def pe_dummy():
    return (DummyTrack(total=100, length=25000, average_template_length=250.0),
            DummyTrack(total=80, length=20000, average_template_length=250.0))


# ------------------------------------
# PeakDetect.__init__
# ------------------------------------

def test_init_reads_options():
    opt, log = opt_and_log(d=150, maxgap=70, minlen=90, shift=-5,
                           gsize=12345, smalllocal=300, largelocal=900,
                           log_pvalue=3.0, log_qvalue=None, PE_MODE=False)
    treat, control = se_dummy()
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    assert (pd.d, pd.maxgap, pd.minlen, pd.end_shift, pd.gsize) == \
        (150, 70, 90, -5, 12345)
    assert (pd.sregion, pd.lregion) == (300, 900)
    assert (pd.log_pvalue, pd.log_qvalue) == (3.0, None)
    assert pd.PE_MODE is False
    assert pd.treat is treat and pd.control is control
    assert pd.opt is opt
    assert pd.ratio_treat2control is None
    assert pd.peaks is None and pd.final_peaks is None
    assert pd.scoretrack is None
    assert log.info == []


@pytest.mark.parametrize("maxgap,minlen,d_arg,exp_maxgap,exp_minlen", [
    (None, None, None, 200, 150),     # max gap = tsize, min length = d
    (0, 0, None, 200, 150),           # 0 counts as unset
    (30, None, None, 30, 150),
    (None, 40, None, 200, 40),
    (None, None, 99, 200, 99),        # explicit d also sets min length
    (25, 35, 99, 25, 35),
])
def test_init_maxgap_minlen_defaults(maxgap, minlen, d_arg, exp_maxgap,
                                     exp_minlen):
    opt, _ = opt_and_log(d=150, tsize=200, maxgap=maxgap, minlen=minlen)
    pd = PeakDetect(opt=opt, treat=None, control=None, d=d_arg)
    assert pd.maxgap == exp_maxgap
    assert pd.minlen == exp_minlen
    assert pd.d == (150 if d_arg is None else d_arg)


@pytest.mark.parametrize("slocal,llocal,exp", [
    (None, None, (200, 600)),
    (500, None, (500, 600)),
    (None, 5000, (200, 5000)),
    (0, 0, (0, 0)),
    (1000, 10000, (1000, 10000)),
])
def test_init_slocal_llocal_arguments(slocal, llocal, exp):
    opt, _ = opt_and_log(smalllocal=200, largelocal=600)
    pd = PeakDetect(opt=opt, slocal=slocal, llocal=llocal)
    assert (pd.sregion, pd.lregion) == exp


def test_init_nolambda_message():
    opt, log = opt_and_log(nolambda=True)
    pd = PeakDetect(opt=opt)
    assert pd.nolambda is True
    assert log.info == ["#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!"]


def test_init_requires_opt():
    with pytest.raises(AttributeError):
        PeakDetect()


def test_init_from_callpeak_options(macs3_argparser, tmp_path):
    """An options namespace from the real callpeak parser and validator,
    completed with what callpeak_cmd sets before peak detection."""
    from MACS3.Utilities.OptValidator import opt_validate_callpeak
    vlog = logging.getLogger("MACS3.Utilities.OptValidator")
    level = vlog.level
    try:
        opt = macs3_argparser.parse_args(
            ["callpeak", "-t", "t.bed", "-c", "c.bed", "-g", "1e6", "-p",
             "0.001", "--outdir", str(tmp_path), "-n", "run"])
        opt = opt_validate_callpeak(opt)
    finally:
        vlog.setLevel(level)
    opt.PE_MODE = False
    opt.tsize = 36
    opt.d = 180
    opt.tocontrol = False
    treat, control = se_dummy()
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    assert (pd.d, pd.maxgap, pd.minlen) == (180, 36, 180)
    assert (pd.sregion, pd.lregion) == (1000, 10000)
    assert pd.gsize == 1e6
    assert pd.log_pvalue == pytest.approx(3.0)
    assert pd.log_qvalue is None
    assert pd.end_shift == 0
    assert pd.nolambda is False


# ------------------------------------
# PeakDetect.call_peaks: arguments given to CallerFromAlignments
# ------------------------------------

@pytest.mark.parametrize("pe", [False, True], ids=["SE", "PE"])
@pytest.mark.parametrize("with_control", [True, False], ids=["ctrl", "noctrl"])
@pytest.mark.parametrize("tocontrol", [False, True], ids=["tolarge", "tosmall"])
@pytest.mark.parametrize("nolambda", [False, True], ids=["lambda", "nolambda"])
def test_caller_arguments(recorder, pe, with_control, tocontrol, nolambda):
    opt, _ = opt_and_log(PE_MODE=pe, tocontrol=tocontrol, nolambda=nolambda,
                         d=150, smalllocal=200, largelocal=600, gsize=10000,
                         shift=7, store_bdg=True, do_SPMR=True)
    treat, control = pe_dummy() if pe else se_dummy()
    control = control if with_control else None
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    pd.call_peaks()
    (caller,) = recorder
    assert caller.treat is treat
    assert caller.control is control
    assert_caller_args(caller.kwargs, expected_caller_args(opt, treat, control))


@pytest.mark.parametrize("ratio", [0.5, 2.0, 1.0])
@pytest.mark.parametrize("tocontrol", [False, True])
def test_caller_arguments_custom_ratio(recorder, ratio, tocontrol):
    opt, _ = opt_and_log(ratio=ratio, tocontrol=tocontrol, d=150,
                         smalllocal=200, largelocal=600, gsize=10000)
    treat, control = se_dummy()
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    pd.call_peaks()
    assert_caller_args(recorder[0].kwargs,
                       expected_caller_args(opt, treat, control))
    assert pd.ratio_treat2control == pytest.approx(
        ratio if ratio != 1.0 else 100 / 40)


@pytest.mark.parametrize("slocal,llocal,ds", [
    (200, 600, [150, 200, 600]),
    (0, 600, [150, 600]),
    (200, 0, [150, 200]),
    (0, 0, [150]),
    (300, 300, [150, 300]),        # llocal not larger than slocal
    (150, 150, [150, 150]),
])
def test_caller_arguments_local_windows(recorder, slocal, llocal, ds):
    opt, _ = opt_and_log(d=150, smalllocal=slocal, largelocal=llocal)
    treat, control = se_dummy()
    PeakDetect(opt=opt, treat=treat, control=control).call_peaks()
    kw = recorder[0].kwargs
    assert kw["ctrl_d_s"] == ds
    assert_caller_args(kw, expected_caller_args(opt, treat, control))


def test_ratio_pe_mode(recorder):
    """PE: fragments count once in the treatment and both ends count in the
    control: ratio = length / (2 * control.total * avg template length)."""
    opt, _ = opt_and_log(PE_MODE=True, d=150, smalllocal=200, largelocal=600)
    treat, control = pe_dummy()
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    pd.call_peaks()
    assert pd.ratio_treat2control == pytest.approx(25000 / (2 * 80 * 250.0))
    assert recorder[0].kwargs["d"] == pytest.approx(250.0)


def test_ratio_pe_mode_control_sum_is_integer(recorder):
    """The PE control sum 2 * 3 * 250.25 = 1501.5 is kept as 1501."""
    opt, _ = opt_and_log(PE_MODE=True, d=150, smalllocal=200, largelocal=600,
                         gsize=10000, tocontrol=True)
    treat = DummyTrack(total=100, length=25001, average_template_length=250.25)
    control = DummyTrack(total=3, length=750, average_template_length=250.0)
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    pd.call_peaks()
    assert pd.ratio_treat2control == 25001 / 1501
    assert recorder[0].kwargs["lambda_bg"] == pytest.approx(1501 / 10000,
                                                            rel=1e-6)


@pytest.mark.parametrize("slocal,llocal,msg", [
    (100, 600, "100 can't be smaller than 150!"),
    (200, 120, "120 can't be smaller than 150!"),
    (400, 300, "300 can't be smaller than 400!"),
])
def test_local_window_assertions(recorder, slocal, llocal, msg):
    opt, _ = opt_and_log(d=150, smalllocal=slocal, largelocal=llocal)
    treat, control = se_dummy()
    pd = PeakDetect(opt=opt, treat=treat, control=control)
    with pytest.raises(AssertionError, match=msg):
        pd.call_peaks()
    assert recorder == []


def test_without_control_llocal_zero(recorder):
    """Without control the lambda window is llocal, so 0 cannot be used."""
    opt, _ = opt_and_log(largelocal=0)
    with pytest.raises(ZeroDivisionError):
        PeakDetect(opt=opt, treat=se_dummy()[0], control=None).call_peaks()


# ------------------------------------
# PeakDetect.call_peaks: which calling method, cutoffs and messages
# ------------------------------------

@pytest.mark.parametrize("with_control", [True, False], ids=["ctrl", "noctrl"])
@pytest.mark.parametrize("log_p,log_q,broad,summits,analysis,expected,msgs", [
    (5.0, None, False, False, False,
     ("call_peaks", ["p"], [5.0], dict(min_length=150, max_gap=200,
                                       call_summits=False,
                                       cutoff_analysis=False)),
     ["#3 Call peaks with given -log10pvalue cutoff: 5.00000 ..."]),
    (5.0, 2.0, False, True, True,
     ("call_peaks", ["p"], [5.0], dict(min_length=150, max_gap=200,
                                       call_summits=True,
                                       cutoff_analysis=True)),
     ["#3 Going to call summits inside each peak ...",
      "#3 Call peaks with given -log10pvalue cutoff: 5.00000 ..."]),
    (None, 2.0, False, False, False,
     ("call_peaks", ["q"], [2.0], dict(min_length=150, max_gap=200,
                                       call_summits=False,
                                       cutoff_analysis=False)),
     []),
    (5.0, None, True, False, True,
     ("call_broadpeaks", ["p"], dict(lvl1_cutoff_s=[5.0], lvl2_cutoff_s=[1.5],
                                     min_length=150, lvl1_max_gap=200,
                                     lvl2_max_gap=800, cutoff_analysis=True)),
     ["#3 Call broad peaks with given level1 -log10pvalue cutoff and "
      "level2: 5.00000, 1.50000..."]),
    (None, 2.0, True, False, False,
     ("call_broadpeaks", ["q"], dict(lvl1_cutoff_s=[2.0], lvl2_cutoff_s=[1.5],
                                     min_length=150, lvl1_max_gap=200,
                                     lvl2_max_gap=800, cutoff_analysis=False)),
     ["#3 Call broad peaks with given level1 -log10qvalue cutoff and "
      "level2: 2.000000, 1.500000..."]),
])
def test_calling_method(recorder, with_control, log_p, log_q, broad, summits,
                        analysis, expected, msgs):
    opt, log = opt_and_log(log_pvalue=log_p, log_qvalue=log_q, broad=broad,
                           call_summits=summits, cutoff_analysis=analysis,
                           log_broadcutoff=1.5, d=150, tsize=200,
                           maxgap=None, minlen=None)
    treat, control = se_dummy()
    pd = PeakDetect(opt=opt, treat=treat,
                    control=control if with_control else None)
    result = pd.call_peaks()
    (caller,) = recorder
    assert caller.calls == [expected, ("destroy",)]
    assert result == ("broad peaks" if broad else "narrow peaks")
    assert pd.peaks == result
    assert log.info == msgs


@pytest.mark.parametrize("with_control", [True, False])
@pytest.mark.parametrize("trackline", [True, False])
def test_trackline_enabled_before_calling(recorder, with_control, trackline):
    opt, _ = opt_and_log(trackline=trackline)
    treat, control = se_dummy()
    PeakDetect(opt=opt, treat=treat,
               control=control if with_control else None).call_peaks()
    calls = [c[0] for c in recorder[0].calls]
    assert calls == (["enable_trackline"] if trackline else []) + \
        ["call_peaks", "destroy"]


@pytest.mark.parametrize("with_control", [True, False])
def test_no_cutoff_given(recorder, with_control):
    """Without a p- or q-value cutoff nothing is called; the caller is still
    destroyed and the missing result raises."""
    opt, _ = opt_and_log(log_pvalue=None, log_qvalue=None)
    treat, control = se_dummy()
    pd = PeakDetect(opt=opt, treat=treat,
                    control=control if with_control else None)
    with pytest.raises(UnboundLocalError):
        pd.call_peaks()
    assert recorder[0].calls == [("destroy",)]


# ------------------------------------
# PeakDetect.call_peaks end to end on small tracks
# ------------------------------------

T_PLUS = [200, 1000, 1005, 1010, 1020, 1030, 1040, 1050, 1060, 1070, 2600,
          3990]
T_MINUS = [1150, 1160, 1170, 1180, 1190, 1200, 1210, 3000, 4500]
C_PLUS = [300, 900, 2500, 4100, 4700]
C_MINUS = [1250, 3300, 4800]


def fw(plus, minus):
    t = FWTrack(fw=50)
    for p in plus:
        t.add_loc(b"chr1", p, 0)
    for p in minus:
        t.add_loc(b"chr1", p, 1)
    t.finalize()
    return t


def _cov(ivs, n):
    a = np.zeros(n + 1, dtype=np.int64)
    for s, e in ivs:
        a[min(max(s, 0), n)] += 1
        a[min(max(e, 0), n)] -= 1
    return np.cumsum(a)[:n]


def se_reference(d, windows, scales, lambda_bg, treat_scale, with_control,
                 shift=0):
    """Per-base treatment pileup and lambda over [0, end of treatment).
    A + tag at p covers [p + shift, p + shift + d), a - tag [p - shift - d,
    p - shift); control window w covers [p - w//2, p + w - w//2) for a +
    tag and [p - (w - w//2), p + w//2) for a - tag."""
    t_iv = [(p + shift, p + shift + d) for p in T_PLUS] + \
        [(p - shift - d, p - shift) for p in T_MINUS]
    n = max(e for _, e in t_iv)
    t = np.maximum(_cov(t_iv, n).astype(F32) * F32(treat_scale), F32(0))
    cp, cm = (C_PLUS, C_MINUS) if with_control else (T_PLUS, T_MINUS)
    c = np.full(n, F32(lambda_bg), dtype=F32)
    for w, s in zip(windows, scales):
        w_iv = [(p - w // 2, p + w - w // 2) for p in cp] + \
            [(p - (w - w // 2), p + w // 2) for p in cm]
        c = np.maximum(c, _cov(w_iv, n).astype(F32) * F32(s))
    return t, c


def rle_bdg(values):
    lines = []
    start = 0
    for i in range(1, len(values) + 1):
        if i == len(values) or values[i] != values[start]:
            lines.append("chr1\t%d\t%d\t%.5f" % (start, i, values[start]))
            start = i
    return lines


def per_base_regions(t, c, cutoff, min_length, max_gap):
    """Regions where -log10 P(X > t) under Poisson(lambda) exceeds cutoff."""
    pairs = {}
    for ti, ci in set(zip(t.astype(int).tolist(), c.tolist())):
        pairs[(ti, ci)] = -np.log10(poisson.sf(ti, ci))
    score = np.array([pairs[(ti, ci)] for ti, ci in
                      zip(t.astype(int).tolist(), c.tolist())])
    above = np.flatnonzero(score > cutoff)
    regions = []
    for x in above:
        if regions and x - regions[-1][1] <= max_gap:
            regions[-1][1] = x + 1
        else:
            regions.append([x, x + 1])
    return [(s, e) for s, e in regions if e - s >= min_length]


@pytest.fixture
def caller_tmp(tmp_path, monkeypatch):
    import tempfile
    d = tmp_path / "pileup_tmp"
    d.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(d))
    return d


def e2e_opt(tmp_path, **kw):
    args = dict(d=100, smalllocal=400, largelocal=1000, gsize=50000,
                tsize=50, maxgap=50, minlen=100, log_pvalue=5.0,
                store_bdg=True, name="run", bdg_treat=str(tmp_path / "t.bdg"),
                bdg_control=str(tmp_path / "c.bdg"),
                cutoff_analysis_file=str(tmp_path / "cut.txt"),
                call_summits=False)
    args.update(kw)
    return opt_and_log(**args)


def read_bdg(path):
    with open(path) as fh:
        return fh.read().splitlines()


@pytest.mark.parametrize("tocontrol", [False, True])
def test_end_to_end_with_control(tmp_path, caller_tmp, tocontrol):
    cutoff = 2.0 if tocontrol else 5.0
    opt, _ = e2e_opt(tmp_path, tocontrol=tocontrol, log_pvalue=cutoff)
    treat, control = fw(T_PLUS, T_MINUS), fw(C_PLUS, C_MINUS)
    peaks = PeakDetect(opt=opt, treat=treat, control=control).call_peaks()
    assert isinstance(peaks, PeakIO)
    ratio = 21 / 8                 # 21 treatment tags, 8 control tags
    if tocontrol:
        lambda_bg, tscale, base = 8 * 100 / 50000, 1 / ratio, 1.0
    else:
        lambda_bg, tscale, base = 21 * 100 / 50000, 1.0, ratio
    t, c = se_reference(100, [100, 400, 1000],
                        [base, base * 100 / 400, base * 100 / 1000],
                        F32(lambda_bg), F32(tscale), True)
    assert read_bdg(tmp_path / "t.bdg") == rle_bdg(t)
    assert read_bdg(tmp_path / "c.bdg") == rle_bdg(c)
    regions = per_base_regions(t, c, cutoff, 100, 50)
    rows = [(p["start"], p["end"]) for p in peaks.get_data_from_chrom(b"chr1")]
    assert rows == regions
    assert regions
    for p in peaks.get_data_from_chrom(b"chr1"):
        tmax = t[p["start"]:p["end"]].max()
        assert p["pileup"] == tmax
        assert t[p["summit"]] == tmax


def test_end_to_end_without_control(tmp_path, caller_tmp):
    opt, _ = e2e_opt(tmp_path)
    peaks = PeakDetect(opt=opt, treat=fw(T_PLUS, T_MINUS),
                       control=None).call_peaks()
    t, c = se_reference(100, [1000], [100 / 1000], F32(100 * 21 / 50000),
                        1.0, False)
    assert read_bdg(tmp_path / "t.bdg") == rle_bdg(t)
    assert read_bdg(tmp_path / "c.bdg") == rle_bdg(c)
    rows = [(p["start"], p["end"]) for p in peaks.get_data_from_chrom(b"chr1")]
    assert rows == per_base_regions(t, c, 5.0, 100, 50)


@pytest.mark.parametrize("shift", [-20, 20])
def test_end_to_end_shift(tmp_path, caller_tmp, shift):
    """opt.shift moves the 5' ends before extension (treatment only)."""
    opt, _ = e2e_opt(tmp_path, shift=shift)
    PeakDetect(opt=opt, treat=fw(T_PLUS, T_MINUS),
               control=fw(C_PLUS, C_MINUS)).call_peaks()
    ratio = 21 / 8
    t, c = se_reference(100, [100, 400, 1000],
                        [ratio, ratio * 100 / 400, ratio * 100 / 1000],
                        F32(21 * 100 / 50000), 1.0, True, shift=shift)
    assert read_bdg(tmp_path / "t.bdg") == rle_bdg(t)
    assert read_bdg(tmp_path / "c.bdg") == rle_bdg(c)


def test_end_to_end_spmr_and_trackline(tmp_path, caller_tmp):
    """--SPMR divides both bedGraphs by the treatment depth in millions
    (21 tags); --trackline adds a track line first."""
    opt, _ = e2e_opt(tmp_path, do_SPMR=True, trackline=True)
    PeakDetect(opt=opt, treat=fw(T_PLUS, T_MINUS),
               control=fw(C_PLUS, C_MINUS)).call_peaks()
    ratio = 21 / 8
    t, c = se_reference(100, [100, 400, 1000],
                        [ratio, ratio * 100 / 400, ratio * 100 / 1000],
                        F32(21 * 100 / 50000), 1.0, True)
    denom = F32(21 / 1e6)
    tl = read_bdg(tmp_path / "t.bdg")
    cl = read_bdg(tmp_path / "c.bdg")
    assert tl[0].startswith('track type=bedGraph name="treatment pileup"')
    assert cl[0].startswith('track type=bedGraph name="control lambda"')
    assert tl[1:] == rle_bdg(t / denom)
    assert cl[1:] == rle_bdg(c / denom)


def test_end_to_end_nolambda(tmp_path, caller_tmp):
    opt, _ = e2e_opt(tmp_path, nolambda=True)
    PeakDetect(opt=opt, treat=fw(T_PLUS, T_MINUS),
               control=fw(C_PLUS, C_MINUS)).call_peaks()
    t, c = se_reference(100, [], [], F32(21 * 100 / 50000), 1.0, True)
    assert read_bdg(tmp_path / "c.bdg") == rle_bdg(c)


def test_end_to_end_matches_direct_caller(tmp_path, caller_tmp):
    """PeakDetect gives the same peaks as CallerFromAlignments called with
    the arguments derived from the definitions."""
    opt, _ = e2e_opt(tmp_path, store_bdg=False)
    treat, control = fw(T_PLUS, T_MINUS), fw(C_PLUS, C_MINUS)
    got = PeakDetect(opt=opt, treat=treat, control=control).call_peaks()
    kw = expected_caller_args(opt, treat, control)
    kw["d"] = int(kw["d"])
    direct = CallerFromAlignments(fw(T_PLUS, T_MINUS), fw(C_PLUS, C_MINUS),
                                  **kw)
    exp = direct.call_peaks(["p"], [5.0], min_length=100, max_gap=50)

    def rows(pk):
        return [(p["start"], p["end"], p["summit"], p["pileup"], p["pscore"],
                 p["qscore"], p["fc"]) for p in pk.get_data_from_chrom(b"chr1")]
    assert rows(got) == rows(exp)
    assert rows(got)


def test_end_to_end_broad(tmp_path, caller_tmp):
    opt, _ = e2e_opt(tmp_path, broad=True, log_broadcutoff=2.0,
                     store_bdg=False)
    treat, control = fw(T_PLUS, T_MINUS), fw(C_PLUS, C_MINUS)
    got = PeakDetect(opt=opt, treat=treat, control=control).call_peaks()
    assert isinstance(got, BroadPeakIO)
    kw = expected_caller_args(opt, treat, control)
    kw["d"] = int(kw["d"])
    exp = CallerFromAlignments(fw(T_PLUS, T_MINUS), fw(C_PLUS, C_MINUS),
                               **kw).call_broadpeaks(
        ["p"], [5.0], [2.0], min_length=100, lvl1_max_gap=50,
        lvl2_max_gap=200)
    keys = ("start", "end", "thickStart", "thickEnd", "blockNum",
            "blockSizes", "blockStarts", "pileup", "pscore", "qscore", "fc")
    assert [[p[k] for k in keys] for p in got.peaks[b"chr1"]] == \
        [[p[k] for k in keys] for p in exp.peaks[b"chr1"]]
    assert got.peaks[b"chr1"]


def test_end_to_end_qvalue(tmp_path, caller_tmp):
    opt, _ = e2e_opt(tmp_path, log_pvalue=None, log_qvalue=2.0)
    peaks = PeakDetect(opt=opt, treat=fw(T_PLUS, T_MINUS),
                       control=fw(C_PLUS, C_MINUS)).call_peaks()
    rows = peaks.get_data_from_chrom(b"chr1")
    assert rows
    assert all(p["qscore"] >= 2.0 for p in rows)


def test_end_to_end_pe(tmp_path, caller_tmp):
    """PE with control: the treatment pileup is fragment coverage."""
    frags = [(100, 300), (950, 1150), (960, 1160), (980, 1200), (1000, 1210),
             (1010, 1190), (2500, 2700), (4000, 4180)]
    cfrags = [(300, 500), (1200, 1400), (2600, 2800), (4300, 4500),
              (4700, 4900)]
    treat, control = PETrackI(), PETrackI()
    for l, r in frags:
        treat.add_loc(b"chr1", l, r)
    for l, r in cfrags:
        control.add_loc(b"chr1", l, r)
    treat.finalize()
    control.finalize()
    opt, _ = e2e_opt(tmp_path, PE_MODE=True, log_pvalue=3.0)
    peaks = PeakDetect(opt=opt, treat=treat, control=control).call_peaks()
    n = 4180
    t = _cov(frags, n).astype(F32)
    assert read_bdg(tmp_path / "t.bdg") == rle_bdg(t)
    # lambda: both ends x of each control fragment cover [x - w//2,
    # x + w//2) for w = d, slocal, llocal; ratio = treatment length /
    # int(2 * control fragments * mean treatment fragment length)
    length = sum(r - l for l, r in frags)
    avg = float(F32(length / len(frags)))
    ratio = length / int(2 * len(cfrags) * avg)   # C long: truncated
    c = np.full(n, F32(length / 50000), dtype=F32)
    for w, s in zip([100, 400, 1000], [ratio, 100 / 400 * ratio,
                                       100 / 1000 * ratio]):
        ends = [x for f in cfrags for x in f]
        c = np.maximum(c, _cov([(x - w // 2, x + w // 2) for x in ends],
                               n).astype(F32) * F32(s))
    assert read_bdg(tmp_path / "c.bdg") == rle_bdg(c)
    regions = per_base_regions(t, c, 3.0, 100, 50)
    assert regions
    assert [(p["start"], p["end"]) for p in
            peaks.get_data_from_chrom(b"chr1")] == regions
