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


# ===========================================================================
# Comprehensive tests for PeakModel and smooth (added below the original
# tests, which are kept unchanged).
#
# Synthetic treatment: binding sites 5000 bp apart.  At a site centred on
# c, k+ reads on the + strand have their 5' end at c - D//2 and k- reads on
# the - strand have their 5' end at c + D - D//2, so the + and - strand
# tag peaks are D bp apart and the fragment length to recover is D.
#
# With bandwidth bw the model uses peaksize P = 2 * bw: each strand's tags
# are extended by bw to both sides, regions above min_tags = round(total *
# lmfold * P / gsize / 2) (strictly) and below max_tags = round(total *
# umfold * P / gsize / 2) (strictly), at least 200 bp long, are strand
# peaks; a + peak and a - peak at most P apart, the + one first, with tag
# counts within a factor of 2, form a pair centred at their midpoint.
# ===========================================================================

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PeakModel import smooth


class ModelLog:
    def __init__(self):
        self.info = []
        self.debug = []
        self.warn = []


def model_opt(gsize=1e6, lmfold=5, umfold=50, bw=100, d_min=20):
    log = ModelLog()
    opt = SimpleNamespace(gsize=gsize, lmfold=lmfold, umfold=umfold, bw=bw,
                          d_min=d_min, info=log.info.append,
                          debug=log.debug.append, warn=log.warn.append)
    return opt, log


def site_track(sites, chroms=(b"chr1",), spacing=5000):
    """sites: list of (D, k_plus, k_minus[, minus_first]); sites are spread
    over ``chroms`` in turn.  Returns (track, list of (chrom, centre))."""
    track = FWTrack(fw=36)
    centres = []
    for i, site in enumerate(sites):
        d, kp, km = site[:3]
        minus_first = len(site) > 3 and site[3]
        chrom = chroms[i % len(chroms)]
        c = 2000 + spacing * (i // len(chroms))
        p_pos, m_pos = c - d // 2, c + d - d // 2
        if minus_first:
            p_pos, m_pos = m_pos, p_pos
        for _ in range(kp):
            track.add_loc(chrom, p_pos, 0)
        for _ in range(km):
            track.add_loc(chrom, m_pos, 1)
        centres.append((chrom, c))
    track.finalize()
    return track, centres


def paired_count(log):
    msgs = [m for m in log.info if m.startswith("#2 Total number of paired peaks: ")]
    assert len(msgs) == 1
    return int(msgs[0].rsplit(" ", 1)[1])


def tag_bounds(total, opt):
    peaksize = 2 * opt.bw
    return (int(round(float(total) * opt.lmfold * peaksize / opt.gsize / 2)),
            int(round(float(total) * opt.umfold * peaksize / opt.gsize / 2)))


def ref_smooth(x, window_len=11, window="hanning"):
    """Weighted moving average centred on each sample.  Beyond the left end
    the signal is mirrored about its first sample (x[-i] = x[i]); beyond
    the right end about the half sample after the last one
    (x[n - 1 + i] = x[n - i]), as in the SciPy cookbook recipe."""
    n = len(x)
    h = window_len // 2
    w = {"flat": np.ones(window_len), "hanning": np.hanning(window_len),
         "hamming": np.hamming(window_len), "bartlett": np.bartlett(window_len),
         "blackman": np.blackman(window_len)}[window]
    w = w / w.sum()

    def at(i):
        if i < 0:
            return x[-i]
        if i >= n:
            return x[2 * n - 1 - i]
        return x[i]
    return np.array([sum(w[j + h] * at(t + j) for j in range(-h, h + 1))
                     for t in range(n)])


def ref_model(plus_line, minus_line, peaksize, d_min):
    """Cross-correlation of the normalised strand profiles for lags
    -P+1 .. P (minus shifted right by the lag), smoothed with an 11-sample
    moving average; d is the lag of the highest local maximum above d_min."""
    m = (minus_line - minus_line.mean()) / (minus_line.std() * len(minus_line))
    p = (plus_line - plus_line.mean()) / (plus_line.std() * len(plus_line))
    n = len(m)
    lags = np.arange(-peaksize + 1, peaksize + 1)
    corr = np.array([np.dot(m[k:], p[:n - k]) if k >= 0 else
                     np.dot(m[:n + k], p[-k:]) for k in lags])
    y = ref_smooth(corr, 11, "flat")
    peaks = [i for i in range(1, len(y) - 1)
             if y[i] > y[i - 1] and y[i] > y[i + 1] and lags[i] > d_min]
    best = max(peaks, key=lambda i: y[i])
    return lags, y, int(lags[best]), sorted(int(lags[i]) for i in peaks)


# ------------------------------------
# NotEnoughPairsException
# ------------------------------------

def test_not_enough_pairs_exception_value_and_str():
    e = NotEnoughPairsException("No enough pairs to build model")
    assert isinstance(e, Exception)
    assert e.value == "No enough pairs to build model"
    assert str(e) == "'No enough pairs to build model'"


@pytest.mark.parametrize("value", [3, None, ("a", 1)])
def test_not_enough_pairs_exception_repr_of_value(value):
    assert str(NotEnoughPairsException(value)) == repr(value)


# ------------------------------------
# PeakModel.__init__ and __str__
# ------------------------------------

def test_init_defaults_before_build():
    opt, log = model_opt()
    model = PeakModel(opt, DummyTreatment(total=10))
    assert (model.d, model.scan_window, model.min_tags) == (0, 0, 0)
    assert model.alternative_d is None
    assert log.info == log.debug == log.warn == []
    assert str(model) == ("\nSummary of Peak Model:\n  Baseline: 0\n"
                          "  Upperline: 0\n  Fragment size: 0\n"
                          "  Scan window size: 0\n")


def test_init_requires_options():
    with pytest.raises(AttributeError):
        PeakModel(SimpleNamespace(gsize=1), DummyTreatment())


# ------------------------------------
# PeakModel.build: tag bounds, pairing and the not-enough-pairs path
# ------------------------------------

@pytest.mark.parametrize("total,lmfold,umfold,bw,gsize", [
    (1000, 10, 30, 100, 1000000),
    (1000, 5, 50, 300, 1000000),
    (12345, 2, 80, 150, 2.7e9),
    (250, 10, 30, 100, 10000),
    (3, 1, 1, 50, 100),
])
def test_build_tag_bounds(total, lmfold, umfold, bw, gsize):
    opt, log = model_opt(gsize=gsize, lmfold=lmfold, umfold=umfold, bw=bw)
    model = PeakModel(opt, DummyTreatment(total=total))
    with pytest.raises(NotEnoughPairsException):
        model.build()
    lo, hi = tag_bounds(total, opt)
    assert model.min_tags == lo
    assert log.debug[0] == "#2 min_tags: %d; max_tags:%d; " % (lo, hi)
    assert "Baseline: %d\n  Upperline: %d\n" % (lo, hi) in str(model)


@pytest.mark.parametrize("n_sites", [0, 1, 50, 99])
def test_build_not_enough_pairs(n_sites):
    opt, log = model_opt(lmfold=5, umfold=10000)
    track, _ = site_track([(100, 5, 5)] * n_sites) if n_sites else (
        DummyTreatment(total=0), None)
    model = PeakModel(opt, track)
    with pytest.raises(NotEnoughPairsException) as exc:
        model.build()
    assert exc.value.value == "No enough pairs to build model"
    assert log.info == ["#2 looking for paired plus/minus strand peaks...",
                        "#2 Total number of paired peaks: %d" % n_sites]
    assert log.warn == [
        "#2 MACS3 needs at least 100 paired peaks at + and - strand to build "
        "the model, but can only find %d! Please make your MFOLD range "
        "broader and try again. If MACS3 still can't build the model, we "
        "suggest to use --nomodel and --extsize 147 or other fixed number "
        "instead." % n_sites,
        "#2 Process for pairing-model is terminated!"]


@pytest.mark.parametrize("n_sites", [100, 101, 150])
def test_build_enough_pairs(n_sites):
    opt, log = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * n_sites)
    lo, hi = tag_bounds(track.total, opt)
    assert lo < 5 < hi
    model = PeakModel(opt, track)
    model.build()
    assert paired_count(log) == n_sites
    assert log.warn == []
    assert log.info[-1] == "#2 Model building with cross-correlation: Done"


def test_build_counts_pairs_on_several_chromosomes():
    opt, log = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 120, chroms=(b"chr1", b"chr2", b"chrX"))
    PeakModel(opt, track).build()
    assert paired_count(log) == 120


def test_build_chromosome_with_one_strand_is_discarded():
    opt, log = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 110)
    for i in range(30):                    # chr2: + strand reads only
        for _ in range(5):
            track.add_loc(b"chr2", 2000 + 5000 * i, 0)
    track.finalize()
    PeakModel(opt, track).build()
    assert paired_count(log) == 110
    assert "Chrom b'chr2' is discarded!" in log.debug


@pytest.mark.parametrize("extra,desc", [
    ((100, 30, 30), "both strands above max_tags"),
    ((100, 1, 1), "both strands not above min_tags"),
    ((100, 5, 2), "+/- tag ratio 2.5"),
    ((100, 2, 5), "+/- tag ratio 0.4"),
    ((100, 5, 5, True), "- strand peak before the + strand peak"),
    ((450, 5, 5), "strand peaks more than 2*bw apart"),
])
def test_build_unpaired_sites(extra, desc):
    """110 good sites plus 30 sites that must not pair."""
    opt, log = model_opt(lmfold=5, umfold=50)
    sites = [(100, 5, 5)] * 110 + [extra] * 30
    track, _ = site_track(sites)
    lo, hi = tag_bounds(track.total, opt)
    assert lo < 2 and 5 < hi < 30
    PeakModel(opt, track).build()
    assert paired_count(log) == 110, desc


def test_build_ratio_inside_bounds_pairs():
    """Tag counts 5 and 3 (ratio 1.67) do pair."""
    opt, log = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 80 + [(100, 5, 3)] * 30)
    PeakModel(opt, track).build()
    assert paired_count(log) == 110


def test_build_mfold_bounds_reject_all():
    """With umfold so low that max_tags <= 5 no strand peak qualifies."""
    opt, log = model_opt(lmfold=1, umfold=4)
    track, _ = site_track([(100, 5, 5)] * 150)
    assert tag_bounds(track.total, opt)[1] <= 5
    with pytest.raises(NotEnoughPairsException):
        PeakModel(opt, track).build()
    assert paired_count(log) == 0


# ------------------------------------
# PeakModel.build: strand profiles, cross-correlation and d
# ------------------------------------

@pytest.mark.parametrize("D,bw", [(100, 100), (150, 100), (80, 100), (51, 100),
                                  (200, 150)])
def test_build_strand_profiles(D, bw):
    """Each tag within P + 5 of a pair centre adds 1 to 10 positions of
    its strand profile, starting P + (tag - centre) (window 2P + 11)."""
    opt, log = model_opt(lmfold=5, umfold=200, bw=bw)
    n = 120
    track, _ = site_track([(D, 5, 5)] * n)
    model = PeakModel(opt, track)
    model.build()
    P = 2 * bw
    W = 1 + 2 * P + 10
    plus = np.zeros(W, dtype="i4")
    minus = np.zeros(W, dtype="i4")
    plus[P - D // 2:P - D // 2 + 10] = 5 * n
    minus[P + D - D // 2:P + D - D // 2 + 10] = 5 * n
    np.testing.assert_array_equal(model.plus_line, plus)
    np.testing.assert_array_equal(model.minus_line, minus)


@pytest.mark.parametrize("D,bw", [(100, 100), (150, 100), (80, 100), (51, 100),
                                  (200, 150)])
def test_build_correlation_profile(D, bw):
    opt, _ = model_opt(lmfold=5, umfold=200, bw=bw)
    track, _ = site_track([(D, 5, 5)] * 120)
    model = PeakModel(opt, track)
    model.build()
    lags, y, d, alt = ref_model(model.plus_line.astype(float),
                                model.minus_line.astype(float), 2 * bw, 20)
    assert len(model.ycorr) == len(model.xcorr) == 4 * bw
    np.testing.assert_allclose(model.ycorr, y, rtol=1e-9, atol=1e-12)
    assert d == D
    assert model.scan_window == 2 * max(model.d, 10)
    assert model.d in model.alternative_d
    assert model.alternative_d == sorted(model.alternative_d)
    assert all(x > 20 for x in model.alternative_d)


@pytest.mark.parametrize("D,bw", [(100, 100), (150, 100), (80, 100),
                                  (200, 150)])
def test_build_d_within_one_of_shift(D, bw):
    opt, log = model_opt(lmfold=5, umfold=200, bw=bw)
    track, _ = site_track([(D, 5, 5)] * 120)
    model = PeakModel(opt, track)
    model.build()
    assert abs(model.d - D) <= 1
    assert "Fragment size: %d\n" % model.d in str(model)


def test_build_d_min_excludes_all_maxima():
    """When d_min is above the only positive-lag maximum no d is found."""
    opt, _ = model_opt(lmfold=5, umfold=200, bw=100, d_min=150)
    track, _ = site_track([(100, 5, 5)] * 120)
    with pytest.raises(AssertionError,
                       match=r"No proper d can be found! Tweak --mfold\?"):
        PeakModel(opt, track).build()


def test_build_two_fragment_sizes_give_alternative_d():
    """Sites with shifts 60 and 160 pooled in one profile give several
    local maxima of the correlation (60, 110 from the cross terms, 160);
    each is reported in alternative_d and d is the highest one."""
    opt, _ = model_opt(lmfold=5, umfold=200, bw=100)
    track, _ = site_track([(160, 5, 5)] * 90 + [(60, 5, 5)] * 60)
    model = PeakModel(opt, track)
    model.build()
    lags, y, d, alt = ref_model(model.plus_line.astype(float),
                                model.minus_line.astype(float), 200, 20)
    assert len(alt) > 1
    assert len(model.alternative_d) == len(alt)
    assert all(abs(a - b) <= 1 for a, b in zip(model.alternative_d, alt))
    assert abs(model.d - d) <= 1


def test_build_with_peakmodel_logging():
    opt, log = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 100)
    PeakModel(opt, track).build()
    assert log.info == ["#2 looking for paired plus/minus strand peaks...",
                        "#2 Total number of paired peaks: 100",
                        "#2 Model building with cross-correlation: Done"]


def test_max_pairnum_is_accepted():
    opt, log = model_opt(lmfold=5, umfold=200)
    track, _ = site_track([(100, 5, 5)] * 120)
    a = PeakModel(opt, track, max_pairnum=10)
    a.build()
    b = PeakModel(opt, track)
    b.build()
    assert a.d == b.d
    np.testing.assert_array_equal(a.ycorr, b.ycorr)


# ------------------------------------
# smooth
# ------------------------------------

SIGNAL = np.array([0.0, 1.0, 3.0, 2.0, 5.0, 4.0, 4.0, 7.0, 1.0, 0.0, 2.0,
                   6.0, 3.0, 3.0, 8.0, 2.0, 1.0, 0.5, 0.0, 9.0])


@pytest.mark.parametrize("window", ["flat", "hanning", "hamming", "bartlett",
                                    "blackman"])
@pytest.mark.parametrize("window_len", [3, 5, 11])
def test_smooth_against_reference(window, window_len):
    y = smooth(SIGNAL, window_len=window_len, window=window)
    assert y.shape == SIGNAL.shape
    np.testing.assert_allclose(y, ref_smooth(SIGNAL, window_len, window),
                               rtol=1e-12, atol=1e-12)


def test_smooth_default_is_hanning_11():
    np.testing.assert_allclose(smooth(SIGNAL), ref_smooth(SIGNAL, 11, "hanning"),
                               rtol=1e-12, atol=1e-12)


def test_smooth_flat_constant_signal():
    x = np.full(30, 2.5)
    np.testing.assert_allclose(smooth(x, 7, "flat"), x)


def test_smooth_flat_interior_is_moving_average():
    x = np.arange(25, dtype=float) ** 2
    y = smooth(x, 5, "flat")
    for t in range(2, 23):
        assert y[t] == pytest.approx(x[t - 2:t + 3].mean())


def test_smooth_window_equal_to_length():
    x = SIGNAL[:11]
    np.testing.assert_allclose(smooth(x, 11, "flat"), ref_smooth(x, 11, "flat"),
                               rtol=1e-12)


def test_smooth_integer_input():
    x = np.array([1, 5, 2, 8, 3, 9, 4], dtype="i4")
    np.testing.assert_allclose(smooth(x, 3, "flat"), ref_smooth(x, 3, "flat"))


@pytest.mark.parametrize("window_len", [0, 1, 2, -3])
def test_smooth_short_window_returns_input(window_len):
    x = SIGNAL.copy()
    assert smooth(x, window_len=window_len) is x


def test_smooth_rejects_2d():
    with pytest.raises(ValueError, match="smooth only accepts 1 dimension arrays."):
        smooth(np.zeros((3, 20)))


@pytest.mark.parametrize("n,window_len", [(10, 11), (0, 3), (2, 3)])
def test_smooth_rejects_short_input(n, window_len):
    with pytest.raises(ValueError,
                       match="Input vector needs to be bigger than window size."):
        smooth(np.zeros(n), window_len=window_len)


def test_smooth_short_input_checked_before_short_window():
    with pytest.raises(ValueError, match="Input vector needs to be bigger"):
        smooth(np.zeros(1), window_len=2)


@pytest.mark.parametrize("window", ["gaussian", "Flat", ""])
def test_smooth_rejects_unknown_window(window):
    with pytest.raises(ValueError, match="Window is on of 'flat', 'hanning', "
                       "'hamming', 'bartlett', 'blackman'"):
        smooth(SIGNAL, 5, window)


def test_smooth_requires_ndarray():
    with pytest.raises(AttributeError):
        smooth(list(SIGNAL))
