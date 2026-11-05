"""Public API and CLI regressions for degenerate HMM models (#814/#757)."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from conftest import MACS3_SCRIPT, macs3_env, no_core_dumps
from MACS3.Signal import HMMR_HMM as hmm
from test_HMMR_HMM import hand_gaussian, hand_poisson, sample_hmm
from test_hmmratac_cmd import atac  # shared synthetic ATAC input fixture


def damage(model, case):
    """Simulate invalid stored parameters left by a final fitting step."""
    if case == "indefinite":
        model.covars_[1] = np.diag([1e12, 1e12, 1e12, -3e-5])
    elif case == "singular":
        model.covars_[1] = 0
    elif case == "asymmetric":
        model.covars_[1, 0, 1] += 2
    elif case == "nan-covariance":
        model.covars_[1, 0, 0] = np.nan
    elif case == "infinite-covariance":
        model.covars_[1, 0, 0] = np.inf
    elif case == "covariance-shape":
        model._covars_ = np.tile(np.eye(3), (3, 1, 1))
    elif case == "zero-diagonal":
        # The public diag getter returns a copy. Mutate the stored
        # representation to simulate fitting bypassing the setter.
        model._covars_[1, 0] = 0
    elif case == "negative-diagonal":
        model._covars_[1, 0] = -1
    elif case == "nan-mean":
        model.means_[1, 0] = np.nan
    elif case == "mean-shape":
        model.means_ = np.ones((3, 3))
    elif case == "zero-transition":
        model.transmat_[1] = 0
    elif case == "negative-transition":
        model.transmat_[1] = [-.1, .5, .6]  # sums to one
    elif case == "nan-transition":
        model.transmat_[1, 0] = np.nan
    elif case == "negative-start":
        model.startprob_ = np.array([-.1, .5, .6])
    elif case == "infinite-start":
        model.startprob_[0] = np.inf
    elif case == "negative-rate":
        model.lambdas_[1, 0] = -1
    elif case == "nan-rate":
        model.lambdas_[1, 0] = np.nan
    else:
        raise AssertionError(case)
    return model


CASES = [
    ("full", "indefinite", "state 1 covariance.*positive-definite"),
    ("full", "singular", "state 1 covariance.*positive-definite"),
    ("full", "asymmetric", "state 1 covariance.*symmetric"),
    ("full", "nan-covariance", "covariance.*finite"),
    ("full", "infinite-covariance", "covariance.*finite"),
    ("full", "covariance-shape", "covariance.*shape"),
    ("diag", "zero-diagonal", "covariance.*positive"),
    ("diag", "negative-diagonal", "covariance.*positive"),
    ("full", "nan-mean", "means_.*finite"),
    ("full", "mean-shape", "means_.*shape"),
    ("full", "zero-transition", "transmat_.*sum to 1"),
    ("full", "negative-transition", "transmat_.*nonnegative"),
    ("full", "nan-transition", "transmat_.*finite"),
    ("full", "negative-start", "startprob_.*nonnegative"),
    ("full", "infinite-start", "startprob_.*finite"),
    ("poisson", "negative-rate", "lambdas_.*nonnegative"),
    ("poisson", "nan-rate", "lambdas_.*finite"),
]


def make_model(kind):
    return hand_poisson() if kind == "poisson" else hand_gaussian(kind)


@pytest.mark.parametrize("operation", ["train", "predict", "save"])
@pytest.mark.parametrize("kind,case,message", CASES)
def test_reject_invalid_parameters(tmp_path, monkeypatch, operation, kind,
                                   case, message):
    model = damage(make_model(kind), case)
    model_type = "poisson" if kind == "poisson" else "gaussian"
    target = tmp_path / "model.json"
    # Rejection must also preserve an existing file, before opening it.
    target.write_text("existing model\n")
    if operation == "train":
        monkeypatch.setattr(type(model), "fit", lambda *a, **kw: model)
        call = lambda: hmm.hmm_training([[1.] * 4] * 3, [3], hmm_type=model_type)
    elif operation == "predict":
        call = lambda: hmm.hmm_predict([[1.] * 4] * 3, [3], model)
    else:
        call = lambda: hmm.hmm_model_save(str(target), model, 10, 2, 1, 0, model_type)
    with pytest.raises(ValueError, match=message) as error:
        call()
    assert "Invalid HMM model:" in str(error.value)
    assert "-l/--lower" in str(error.value)
    assert "-u/--upper" in str(error.value)
    assert target.read_text() == "existing model\n"


@pytest.mark.parametrize("kind,case,message", CASES)
def test_reject_invalid_saved_model(tmp_path, kind, case, message):
    model = make_model(kind)
    model_type = "poisson" if kind == "poisson" else "gaussian"
    path = tmp_path / "model.json"
    hmm.hmm_model_save(str(path), model, 10, 2, 1, 0, model_type)
    document = json.loads(path.read_text())
    damage(model, case)
    document.update(startprob=model.startprob_.tolist(),
                    transmat=model.transmat_.tolist())
    if kind == "poisson":
        document["lambdas"] = model.lambdas_.tolist()
    else:
        document["means"] = model.means_.tolist()
        document["covars"] = model._covars_.tolist()
    path.write_text(json.dumps(document))
    # hmmlearn's covariance setter may reject before our post-load check.
    with pytest.raises(ValueError, match="Invalid HMM model"):
        hmm.hmm_model_init(str(path))


@pytest.mark.parametrize("covariance_type", ["full", "diag", "tied", "spherical"])
def test_valid_single_sequence_covariance_types(covariance_type):
    observations, lengths, _ = sample_hmm(814, n_seq=1, length=240)
    model = hmm.hmm_training(observations, lengths, covar=covariance_type)
    expected = model.predict_proba(observations, lengths)
    actual = hmm.hmm_predict(observations, lengths, model)
    np.testing.assert_array_equal(actual, expected)
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual.sum(axis=1), 1)


def test_valid_small_variance_is_not_rejected():
    model = hand_gaussian()
    model.covars_[1] = np.diag([1e-12, 1., 1., 1.])
    observations = model.means_.tolist()
    np.testing.assert_array_equal(
        hmm.hmm_predict(observations, [3], model),
        model.predict_proba(observations, [3]))


def test_integer_covariance_round_trip_preserves_json(tmp_path):
    model = hand_gaussian()
    model.covars_ = np.tile(np.eye(4, dtype=int), (3, 1, 1))
    original = str(tmp_path / "original.json")
    restored = str(tmp_path / "restored.json")
    hmm.hmm_model_save(original, model, 10, 2, 1, 0, "gaussian")
    loaded = hmm.hmm_model_init(original)[0]
    hmm.hmm_model_save(restored, loaded, 10, 2, 1, 0, "gaussian")
    assert Path(restored).read_text() == Path(original).read_text()


@pytest.mark.parametrize("kind", ["full", "poisson"])
def test_prediction_infers_features_for_manually_configured_models(kind):
    expected_model = make_model(kind)
    model = make_model(kind)
    del model.n_features
    observations = [[1.] * 4] * 3
    np.testing.assert_array_equal(
        hmm.hmm_predict(observations, [3], model),
        expected_model.predict_proba(observations, [3]))


def test_valid_parameters_do_not_require_convergence(monkeypatch):
    model = hand_gaussian()
    assert not model.monitor_.converged
    monkeypatch.setattr(hmm.GaussianHMM, "fit", lambda *a, **kw: model)
    assert hmm.hmm_training([[1.] * 4] * 3, [3]) is model


def test_zero_poisson_rates_and_zero_probabilities_are_valid():
    model = hand_poisson()
    model.lambdas_[0, 0] = 0
    model.startprob_ = np.array([1., 0., 0.])
    model.transmat_ = np.eye(3)
    observations = [[0., 1., 1., 1.]] * 3
    np.testing.assert_array_equal(
        hmm.hmm_predict(observations, [3], model),
        model.predict_proba(observations, [3]))


def test_single_state_poisson_accepts_flat_rates():
    model = hmm.PoissonHMM(n_components=1)
    model.startprob_ = np.ones(1)
    model.transmat_ = np.ones((1, 1))
    model.lambdas_ = np.array([0., 1., 2., 3.])
    np.testing.assert_array_equal(
        hmm.hmm_predict([[0., 1., 2., 3.]] * 3, [3], model), np.ones((3, 1)))


@pytest.mark.parametrize("operation", ["train", "predict"])
@pytest.mark.parametrize("observations,lengths", [
    ([], []), ([[1.] * 4] * 3, [2]), ([[np.nan] * 4] * 3, [3]),
    ([[1.] * 4, [1.] * 3], [2]),
])
def test_empty_and_malformed_observations(operation, observations, lengths):
    with pytest.raises(ValueError):
        if operation == "train":
            hmm.hmm_training(observations, lengths)
        else:
            hmm.hmm_predict(observations, lengths, hand_gaussian())


def run_probe(tmp_path, case, *args):
    return subprocess.run(
        [sys.executable, str(Path(__file__).with_name("hmm_model_probe.py")),
         hmm.__file__, case, *map(str, args)],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True, text=True,
        timeout=30, preexec_fn=no_core_dumps)


def test_constant_training_rejected_in_compiled_build(tmp_path):
    result = run_probe(tmp_path, "constant-training")
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("extra", [[], ["--modelonly"]])
def test_cli_rejects_degenerate_fit_before_saving(tmp_path, atac, extra):
    output = tmp_path / "out"
    result = run_probe(tmp_path, "cli-training", MACS3_SCRIPT, atac.bedpe,
                       output, *extra)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "state 1 covariance must be positive-definite" in result.stderr
    assert "-l/--lower" in result.stderr
    assert "HMM converged:" not in result.stderr
    assert "Write HMM parameters" not in result.stderr
    assert "Traceback" not in result.stderr
    assert not (output / "NA_model.json").exists()
    assert not (output / "NA_accessible_regions.narrowPeak").exists()


@pytest.mark.parametrize("kind,case", [
    ("full", "indefinite"), ("full", "asymmetric"),
    ("diag", "zero-diagonal"), ("poisson", "nan-rate"),
])
def test_cli_rejects_invalid_loaded_model(tmp_path, atac, run_macs3, kind, case):
    model = make_model(kind)
    path = tmp_path / "invalid.json"
    hmm.hmm_model_save(str(path), model, 10, 2, 1, 0,
                       "poisson" if kind == "poisson" else "gaussian")
    document = json.loads(path.read_text())
    damage(model, case)
    document["lambdas" if kind == "poisson" else "covars"] = (
        model.lambdas_ if kind == "poisson" else model._covars_).tolist()
    path.write_text(json.dumps(document))
    output = tmp_path / "out"
    result = run_macs3(["hmmratac", "-i", atac.bedpe, "-f", "BEDPE",
                        "--model", path, "--outdir", output,
                        "-l", "4", "-u", "20", "-c", "2"], timeout=30)
    assert result.returncode == 1, result.stderr
    assert "Invalid HMM model:" in result.stderr
    assert "Retrain if using a saved model" in result.stderr
    assert "Traceback" not in result.stderr
    assert not (output / "NA_accessible_regions.narrowPeak").exists()
