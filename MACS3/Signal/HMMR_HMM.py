# cython: language_level=3
# Time-stamp: <2024-10-04 11:45:30 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
# ------------------------------------
# python modules
# ------------------------------------
import cython
import numpy as np
from cython.cimports import numpy as cnp
from hmmlearn.hmm import GaussianHMM, PoissonHMM
import json

# ------------------------------------
# MACS3 modules
# ------------------------------------
from MACS3.Signal.Prob import pnorm2

# ------------------------------------
# Misc functions
# ------------------------------------

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html


@cython.cfunc
@cython.inline
def get_weighted_density(x: cython.int, m: cython.float,
                         v: cython.float, w: cython.float) -> cython.float:
    """Description:

    parameters:
      1. x: the observed value
      2. m: the mean of gaussian
      3. v: the variance of the gaussian
      4. w: the weight
    return value:
    """
    return w * pnorm2(float(x), m, v)

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# public functions
# ------------------------------------


def _invalid_hmm_model(reason):
    return ValueError(
        f"Invalid HMM model: {reason}. Review the training regions and "
        "adjust -l/--lower or -u/--upper cutoffs, or supply representative "
        "regions with -t/--training. Retrain if using a saved model.")


def _hmm_array(value, shape, name):
    try:
        values = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise _invalid_hmm_model(f"{name} must contain numeric values") from exc
    if values.dtype.kind not in "biuf":
        raise _invalid_hmm_model(f"{name} must contain numeric values")
    if values.shape != shape:
        raise _invalid_hmm_model(f"{name} must have shape {shape}")
    if not np.isfinite(values).all():
        raise _invalid_hmm_model(f"{name} must contain only finite values")
    return values


def _validate_hmm_model(hmm_model):
    """Reject unusable fitted parameters without modifying the model.

    A convergence flag does not guarantee usable emission parameters.
    In particular, hmmlearn's final M-step can leave invalid covariances
    that are only discovered on the next prediction.
    """
    n = getattr(hmm_model, "n_components", None)
    f = getattr(hmm_model, "n_features", None)
    if not hasattr(hmm_model, "n_features"):
        # hmmlearn also accepts manually configured models and infers
        # their feature count from the emissions on first prediction.
        emissions = (getattr(hmm_model, "means_", None)
                     if isinstance(hmm_model, GaussianHMM)
                     else getattr(hmm_model, "lambdas_", None))
        try:
            if isinstance(hmm_model, PoissonHMM):
                emissions = np.atleast_2d(emissions)
            f = np.asarray(emissions).shape[1]
        except (TypeError, ValueError, IndexError):
            pass
    if not isinstance(n, (int, np.integer)) or n <= 0:
        raise _invalid_hmm_model("n_components must be a positive integer")
    if not isinstance(f, (int, np.integer)) or f <= 0:
        raise _invalid_hmm_model("n_features must be a positive integer")
    for name, shape in (("startprob_", (n,)), ("transmat_", (n, n))):
        probabilities = _hmm_array(getattr(hmm_model, name, None), shape, name)
        if (np.any(probabilities < 0)
                or not np.allclose(probabilities.sum(axis=-1), 1)):
            raise _invalid_hmm_model(
                f"{name} must be nonnegative and sum to 1 per distribution")

    if isinstance(hmm_model, PoissonHMM):
        rates = _hmm_array(np.atleast_2d(getattr(hmm_model, "lambdas_", None)),
                           (n, f), "lambdas_")
        if np.any(rates < 0):
            raise _invalid_hmm_model("lambdas_ must be nonnegative")
    elif isinstance(hmm_model, GaussianHMM):
        _hmm_array(getattr(hmm_model, "means_", None), (n, f), "means_")
        kind = hmm_model.covariance_type
        # Inspect hmmlearn's stored representation: the public covars_
        # property expands diagonal/tied/spherical arrays and can conceal
        # malformed shapes (or expand spherical training arrays repeatedly).
        covars = getattr(hmm_model, "_covars_", None)
        if kind == "spherical":
            shape = np.shape(covars)
            if shape not in ((n,), (n, 1), (n, f)):
                raise _invalid_hmm_model("spherical covariance has an invalid shape")
        elif kind in ("full", "diag", "tied"):
            shape = {"full": (n, f, f), "diag": (n, f), "tied": (f, f)}[kind]
        else:
            raise _invalid_hmm_model(f"unknown covariance type {kind}")
        covars = _hmm_array(covars, shape, "covariance")
        if kind in ("diag", "spherical"):
            if np.any(covars <= 0):
                raise _invalid_hmm_model("covariance variances must be positive")
        else:
            matrices = covars if kind == "full" else covars[np.newaxis, :, :]
            for state, covariance in enumerate(matrices):
                label = f"state {state} covariance" if kind == "full" else "shared covariance"
                if not np.allclose(covariance, covariance.T):
                    raise _invalid_hmm_model(f"{label} must be symmetric")
                try:
                    np.linalg.cholesky(covariance.astype(float, copy=False))
                except np.linalg.LinAlgError as exc:
                    raise _invalid_hmm_model(f"{label} must be positive-definite") from exc


@cython.ccall
def hmm_training(training_data: list, training_data_lengths: list,
                 n_states: cython.int = 3, random_seed: cython.int = 12345,
                 covar: str = 'full', hmm_type: str = 'gaussian'):
    # training data should be in array like format: X = np.array([[.5, .3, .1, .1], [.6, .4, 0, 0]])
    # if we do not want init_prob to be updated through learning, set params = 'tmc' and init_prob = initial_state otherwise it will be overwritten
    # according to base documentation, if init_prob not stated, it is set to be equally likely for any state (1/ # of components)
    # if we have other known parameters, we should set these (ie: means_weights, covariance_type etc.)
    rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(random_seed)))
    if hmm_type == "gaussian":
        hmm_model = GaussianHMM(n_components=n_states, covariance_type=covar, random_state=rs, verbose=False)
    if hmm_type == "poisson":
        hmm_model = PoissonHMM(n_components=n_states, random_state=rs, verbose=False)
    hmm_model = hmm_model.fit(training_data, training_data_lengths)
    assert hmm_model.n_features == 4
    _validate_hmm_model(hmm_model)
    return hmm_model


@cython.ccall
@cython.returns(cnp.ndarray)
def hmm_predict(signals: list, lens: list, hmm_model):
    _validate_hmm_model(hmm_model)
    predictions = hmm_model.predict_proba(signals, lens)
    return predictions


@cython.ccall
def hmm_model_save(model_file: str, hmm_model, hmm_binsize: cython.int,
                   i_open_region: cython.int, i_nucleosomal_region: cython.int,
                   i_background_region: cython.int, hmm_type: str):
    if hmm_type == "gaussian":
        if hmm_model.covariance_type not in ("diag", "full"):
            raise Exception(f"Unknown covariance type {hmm_model.covariance_type}")
        _validate_hmm_model(hmm_model)
        if hmm_model.covariance_type == "diag":
            covars = hmm_model.covars_.diagonal(axis1=1, axis2=2)
        elif hmm_model.covariance_type == "full":
            covars = hmm_model.covars_
        with open(model_file, "w") as f:
            json.dump({"startprob": hmm_model.startprob_.tolist(),
                       "transmat": hmm_model.transmat_.tolist(),
                       "means": hmm_model.means_.tolist(),
                       "covars": covars.tolist(),
                       "covariance_type": hmm_model.covariance_type,
                       "n_features": int(hmm_model.n_features),
                       "i_open_region": int(i_open_region),
                       "i_background_region": int(i_background_region),
                       "i_nucleosomal_region": int(i_nucleosomal_region),
                       "hmm_binsize": int(hmm_binsize),
                       "hmm_type": hmm_type}, f)
    if hmm_type == "poisson":
        _validate_hmm_model(hmm_model)
        with open(model_file, "w") as f:
            json.dump({"startprob": hmm_model.startprob_.tolist(),
                       "transmat": hmm_model.transmat_.tolist(),
                       "lambdas": hmm_model.lambdas_.tolist(),
                       "n_features": int(hmm_model.n_features),
                       "i_open_region": int(i_open_region),
                       "i_background_region": int(i_background_region),
                       "i_nucleosomal_region": int(i_nucleosomal_region),
                       "hmm_binsize": int(hmm_binsize),
                       "hmm_type": hmm_type}, f)


@cython.ccall
def hmm_model_init(model_file: str) -> list:
    with open(model_file) as f:
        m = json.load(f)
        hmm_type = m.get("hmm_type", "gaussian")  # if no model_type written, assume older file and use gaussian
        if hmm_type == "gaussian":
            hmm_model = GaussianHMM(n_components=3, covariance_type=m["covariance_type"])
            hmm_model.means_ = np.array(m["means"])
            try:
                covars = np.array(m["covars"])
                hmm_model.covars_ = covars
            except (ValueError, TypeError, IndexError) as exc:
                raise _invalid_hmm_model(f"invalid covariance: {exc}") from exc
            hmm_model.covariance_type = m["covariance_type"]
        if hmm_type == "poisson":
            hmm_model = PoissonHMM(n_components=3)
            hmm_model.lambdas_ = np.array(m["lambdas"])
        hmm_model.startprob_ = np.array(m["startprob"])
        hmm_model.transmat_ = np.array(m["transmat"])
        hmm_model.n_features = m["n_features"]
        _validate_hmm_model(hmm_model)
        return [hmm_model, m["i_open_region"], m["i_background_region"],
                m["i_nucleosomal_region"], m["hmm_binsize"], hmm_type]
