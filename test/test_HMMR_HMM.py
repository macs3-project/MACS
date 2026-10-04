import unittest
import pytest
from MACS3.Signal.HMMR_HMM import hmm_training, hmm_predict, GaussianHMM
import numpy as np
import numpy.testing as npt
import json
from hmmlearn.hmm import PoissonHMM
from scipy.stats import (multivariate_normal,
                         poisson)
from MACS3.Signal.HMMR_HMM import (hmm_model_save,
                                   hmm_model_init)

# ------------------------------------
# Main function
# ------------------------------------
''' This unittest is to check the ouputs of the hmm_training() and hmm_predict() functions
'''


# @pytest.mark.skip(reason="need to refine later")
class Test_HMM_train(unittest.TestCase):
    def setUp(self):
        self.training_data = np.loadtxt("test/large_training_data.txt",
                                        delimiter="\t", dtype="float",
                                        usecols=(2, 3, 4, 5)).tolist()
        self.training_data_lengths = np.loadtxt('test/large_training_lengths.txt', dtype="int").tolist()
        self.expected_converged = True
        self.not_expected_covars = None
        self.not_expected_means = None
        self.not_expected_transmat = None

        self.startprob = [0.01807016, 0.90153727, 0.08039257]
        self.means = [[2.05560411e-01, 1.52959594e+00, 1.73568556e+00, 1.00019720e-04],
                      [1.84467806e-01, 1.46784946e+00, 1.67895745e+00, 1.00016654e-04],
                      [2.06402305e+00, 8.60140461e+00, 7.22907032e+00, 1.00847661e-04]]
        self.covars = [[[1.19859257e-01, 5.33746506e-02, 3.99871507e-02, 1.49805047e-07],
                        [5.33746506e-02, 1.88774896e+00, 7.38204761e-01, 1.70902908e-07],
                        [3.99871507e-02, 7.38204761e-01, 2.34175176e+00, 1.75654357e-07],
                        [1.49805047e-07, 1.70902908e-07, 1.75654357e-07, 1.45312288e-07]],
                       [[1.06135330e-01, 4.16846792e-02, 3.24447289e-02, 1.30393434e-07],
                        [4.16846792e-02, 1.75537103e+00, 6.70848135e-01, 1.49425940e-07],
                        [3.24447289e-02, 6.70848135e-01, 2.22285392e+00, 1.52914017e-07],
                        [1.30393434e-07, 1.49425940e-07, 1.52914017e-07, 1.27205162e-07]],
                       [[5.94746590e+00, 5.24388615e+00, -5.33166471e-01, -1.47228883e-06],
                        [5.24388615e+00, 2.63945986e+01, 3.54212739e+00, -6.03892201e-06],
                        [-5.33166471e-01, 3.54212739e+00, 1.50231166e+01, 1.43141422e-05],
                        [-1.47228883e-06, -6.03892201e-06, 1.43141422e-05, 1.04240673e-07]]]
        self.transmat = [[1.91958645e-03, 9.68166646e-01, 2.99137676e-02],
                         [8.52453717e-01, 1.46924953e-01, 6.21329356e-04],
                         [2.15432113e-02, 6.80080650e-05, 9.78388781e-01]]
        self.n_features = 4

        # for prediction
        self.prediction_data = np.loadtxt("test/small_prediction_data.txt",
                                          delimiter="\t",
                                          dtype="float",
                                          usecols=(2,3,4,5)).tolist()
        self.prediction_data_lengths = np.loadtxt('test/small_prediction_lengths.txt',
                                                  dtype="int").tolist()
        self.predictions = np.loadtxt('test/small_prediction_results.txt',
                                      delimiter="\t",
                                      dtype="float").tolist()

    # @pytest.mark.skip(reason="it may fail with different sklearn+hmmlearn")
    def test_training(self):
        # test hmm_training:
        model = hmm_training(training_data=self.training_data, training_data_lengths=self.training_data_lengths, n_states=3, random_seed=12345, covar='full')
        print(model.startprob_)
        print(model.means_)
        print(model.covars_)
        print(model.transmat_)
        print(model.n_features)
        self.assertEqual(model.monitor_.converged, self.expected_converged)
        self.assertNotEqual(model.covars_.tolist(), self.not_expected_covars)
        self.assertNotEqual(model.means_.tolist(), self.not_expected_means)
        self.assertNotEqual(model.transmat_.tolist(), self.not_expected_transmat)
        npt.assert_allclose(model.startprob_.tolist(), self.startprob, rtol=1e-5)
        npt.assert_allclose(model.means_, self.means, rtol=1e-5)
        npt.assert_allclose(model.covars_, self.covars, rtol=1e-5)
        npt.assert_allclose(model.transmat_, self.transmat, rtol=1e-5)
        npt.assert_allclose(model.n_features, self.n_features, rtol=1e-5)

    @pytest.mark.skip(reason="it may fail with different sklearn+hmmlearn")
    def test_predict(self):
        # test hmm_predict
        hmm_model = GaussianHMM(n_components=3, covariance_type='full')
        hmm_model.startprob_ = np.array(self.startprob)
        hmm_model.transmat_ = np.array(self.transmat)
        hmm_model.means_ = np.array(self.means)
        hmm_model.covars_ = np.array(self.covars)
        hmm_model.covariance_type = 'full'
        hmm_model.n_features = self.n_features
        predictions = hmm_predict(self.prediction_data,
                                  self.prediction_data_lengths,
                                  hmm_model)

        ## This is to write the prediction results into a file for 'correct' answer
        #with open("test/small_prediction_results.txt","w") as f:
        #    for x,y,z in predictions:
        #        f.write(str(x)+"\t"+str(y)+"\t"+str(z)+"\n")

        #npt.assert_allclose(predictions, self.predictions, rtol=1e-5)


# ------------------------------------
# Helpers for the tests below
# ------------------------------------

# emission means (gaussian) or lambdas (poisson) of three well
# separated states for the short, mono, di and tri signals
TRUE_MEANS = np.array([[0.5, 1.0, 0.5, 0.2],
                       [6.0, 9.0, 2.0, 0.5],
                       [12.0, 4.0, 7.0, 3.0]])
TRUE_TRANS = np.array([[0.90, 0.08, 0.02],
                       [0.10, 0.80, 0.10],
                       [0.05, 0.15, 0.80]])
HAND_STARTPROB = np.array([0.6, 0.3, 0.1])


def sample_hmm(seed, n_seq=5, length=80, kind="gaussian"):
    """Sample sequences from a 3-state HMM with 4 emission features.

    Returns (observations as a list of lists, lengths, states).
    """
    rs = np.random.RandomState(seed)
    X = []
    states = []
    for _ in range(n_seq):
        s = rs.randint(3)
        for t in range(length):
            if t:
                s = rs.choice(3, p=TRUE_TRANS[s])
            states.append(s)
            if kind == "poisson":
                X.append(rs.poisson(TRUE_MEANS[s]).tolist())
            else:
                X.append((TRUE_MEANS[s] + rs.normal(0, 0.5, 4)).tolist())
    return X, [length] * n_seq, states


def hand_covars():
    """Three positive definite 4x4 covariance matrices."""
    return np.array([(0.5 + 0.2 * k) * np.eye(4) + 0.1 * np.ones((4, 4))
                     for k in range(3)])


def hand_gaussian(covariance_type="full"):
    """A GaussianHMM with parameters set by hand."""
    m = GaussianHMM(n_components=3, covariance_type=covariance_type)
    m.startprob_ = HAND_STARTPROB.copy()
    m.transmat_ = TRUE_TRANS.copy()
    m.means_ = TRUE_MEANS.copy()
    if covariance_type == "full":
        m.covars_ = hand_covars()
    else:
        m.covars_ = np.array([[0.6, 0.7, 0.8, 0.9],
                              [1.0, 1.1, 1.2, 1.3],
                              [0.5, 0.5, 2.0, 2.0]])
    m.n_features = 4
    return m


def hand_poisson():
    """A PoissonHMM with parameters set by hand."""
    m = PoissonHMM(n_components=3)
    m.startprob_ = HAND_STARTPROB.copy()
    m.transmat_ = TRUE_TRANS.copy()
    m.lambdas_ = TRUE_MEANS + 0.1
    m.n_features = 4
    return m


def gaussian_observations(seed=0, n=12):
    rs = np.random.RandomState(seed)
    rows = [TRUE_MEANS[i % 3] + rs.normal(0, 1.0, 4) for i in range(n)]
    return [r.tolist() for r in rows]


def poisson_observations(seed=0, n=12):
    rs = np.random.RandomState(seed)
    return [rs.poisson(TRUE_MEANS[i % 3] + 0.1).tolist() for i in range(n)]


def ref_posteriors(B, startprob, transmat):
    """Scaled forward-backward: posterior state probabilities of one
    sequence given its emission likelihoods ``B`` (T x K)."""
    T, K = B.shape
    alpha = np.zeros((T, K))
    c = np.zeros(T)
    a = startprob * B[0]
    c[0] = a.sum()
    alpha[0] = a / c[0]
    for t in range(1, T):
        a = (alpha[t - 1] @ transmat) * B[t]
        c[t] = a.sum()
        alpha[t] = a / c[t]
    beta = np.ones((T, K))
    for t in range(T - 2, -1, -1):
        beta[t] = transmat @ (B[t + 1] * beta[t + 1]) / c[t + 1]
    post = alpha * beta
    return post / post.sum(axis=1, keepdims=True)


def gaussian_likelihoods(X, means, covars):
    return np.array([[multivariate_normal.pdf(x, means[k], covars[k])
                      for k in range(len(means))] for x in X])


def poisson_likelihoods(X, lambdas):
    return np.array([[np.prod(poisson.pmf(x, lambdas[k]))
                      for k in range(len(lambdas))] for x in X])


def order_by_sum(a):
    """Rows of ``a`` ordered by their sum (states have no fixed order)."""
    a = np.asarray(a)
    return a[np.argsort(a.sum(axis=1))]


# ------------------------------------
# hmm_training
# ------------------------------------

def test_hmm_training_gaussian_shapes():
    X, lens, _ = sample_hmm(1)
    model = hmm_training(X, lens, random_seed=7)
    assert isinstance(model, GaussianHMM)
    assert model.covariance_type == "full"
    assert model.n_components == 3
    assert model.n_features == 4
    assert model.startprob_.shape == (3,)
    assert model.transmat_.shape == (3, 3)
    assert model.means_.shape == (3, 4)
    assert model.covars_.shape == (3, 4, 4)


@pytest.mark.parametrize("hmm_type", ["gaussian", "poisson"])
def test_hmm_training_probabilities_normalized(hmm_type):
    X, lens, _ = sample_hmm(2, kind=hmm_type)
    model = hmm_training(X, lens, random_seed=7, hmm_type=hmm_type)
    assert model.startprob_.sum() == pytest.approx(1.0, abs=1e-10)
    npt.assert_allclose(model.transmat_.sum(axis=1), np.ones(3), atol=1e-10)
    assert np.all(model.transmat_ >= 0)
    assert np.all(model.startprob_ >= 0)


def test_hmm_training_gaussian_recovers_emission_means():
    # hmmlearn stops after its default 10 Baum-Welch iterations; with this
    # seed the fit converges earlier, at the true parameters (other seeds
    # can stop in a local optimum that merges two states)
    X, lens, _ = sample_hmm(3, n_seq=6, length=100)
    model = hmm_training(X, lens, random_seed=1)
    assert model.monitor_.iter < 10
    npt.assert_allclose(order_by_sum(model.means_), order_by_sum(TRUE_MEANS),
                        atol=0.5)


def test_hmm_training_gaussian_full_covariances_are_symmetric_pd():
    X, lens, _ = sample_hmm(4)
    model = hmm_training(X, lens, random_seed=11)
    for c in model.covars_:
        npt.assert_allclose(c, c.T, atol=1e-12)
        assert np.all(np.linalg.eigvalsh(c) > 0)


def test_hmm_training_poisson_shapes_and_lambdas():
    X, lens, _ = sample_hmm(5, n_seq=6, length=100, kind="poisson")
    model = hmm_training(X, lens, random_seed=11, hmm_type="poisson")
    assert isinstance(model, PoissonHMM)
    assert model.n_features == 4
    assert model.lambdas_.shape == (3, 4)
    npt.assert_allclose(order_by_sum(model.lambdas_),
                        order_by_sum(TRUE_MEANS), atol=1.0)


@pytest.mark.parametrize("hmm_type", ["gaussian", "poisson"])
def test_hmm_training_same_seed_same_model(hmm_type):
    X, lens, _ = sample_hmm(6, kind=hmm_type)
    a = hmm_training(X, lens, random_seed=12345, hmm_type=hmm_type)
    b = hmm_training(X, lens, random_seed=12345, hmm_type=hmm_type)
    npt.assert_array_equal(a.startprob_, b.startprob_)
    npt.assert_array_equal(a.transmat_, b.transmat_)
    if hmm_type == "gaussian":
        npt.assert_array_equal(a.means_, b.means_)
        npt.assert_array_equal(a.covars_, b.covars_)
    else:
        npt.assert_array_equal(a.lambdas_, b.lambdas_)


def test_hmm_training_default_seed_is_12345():
    X, lens, _ = sample_hmm(6)
    a = hmm_training(X, lens)
    b = hmm_training(X, lens, random_seed=12345)
    npt.assert_array_equal(a.means_, b.means_)
    npt.assert_array_equal(a.transmat_, b.transmat_)


def test_hmm_training_diag_covariance():
    X, lens, _ = sample_hmm(7)
    model = hmm_training(X, lens, random_seed=3, covar="diag")
    assert model.covariance_type == "diag"
    # hmmlearn's covars_ expands diagonal covariances to full matrices
    assert model.covars_.shape == (3, 4, 4)
    for c in model.covars_:
        npt.assert_array_equal(c, np.diag(np.diag(c)))


def test_hmm_training_n_states():
    X, lens, _ = sample_hmm(8)
    model = hmm_training(X, lens, n_states=2, random_seed=3)
    assert model.means_.shape == (2, 4)
    assert model.transmat_.shape == (2, 2)


def test_hmm_training_single_sequence():
    X, lens, _ = sample_hmm(9, n_seq=1, length=200)
    model = hmm_training(X, lens, random_seed=3)
    assert model.means_.shape == (3, 4)


def test_hmm_training_requires_four_features():
    X, lens, _ = sample_hmm(10)
    X3 = [row[:3] for row in X]
    with pytest.raises(AssertionError):
        hmm_training(X3, lens, random_seed=3)


def test_hmm_training_unknown_type_raises():
    # no model is built for an unknown type: the function fails when it
    # uses the unassigned model variable
    X, lens, _ = sample_hmm(10)
    with pytest.raises(UnboundLocalError):
        hmm_training(X, lens, random_seed=3, hmm_type="gamma")


def test_hmm_training_requires_lists():
    X, lens, _ = sample_hmm(10)
    with pytest.raises(TypeError):
        hmm_training(np.array(X), lens, random_seed=3)


def test_hmm_training_lengths_must_sum_to_observations():
    X, lens, _ = sample_hmm(10)
    with pytest.raises(ValueError, match="doesn't sum to"):
        hmm_training(X, lens[:-1], random_seed=3)


# ------------------------------------
# hmm_predict
# ------------------------------------

@pytest.mark.parametrize("model_fn, obs_fn", [
    (hand_gaussian, gaussian_observations),
    (hand_poisson, poisson_observations),
])
def test_hmm_predict_shape_and_rows_sum_to_one(model_fn, obs_fn):
    X = obs_fn(n=15)
    pred = hmm_predict(X, [7, 8], model_fn())
    assert isinstance(pred, np.ndarray)
    assert pred.shape == (15, 3)
    npt.assert_allclose(pred.sum(axis=1), np.ones(15), atol=1e-10)
    assert np.all(pred >= 0) and np.all(pred <= 1)


def test_hmm_predict_gaussian_against_forward_backward():
    X = gaussian_observations(seed=1, n=12)
    pred = hmm_predict(X, [12], hand_gaussian())
    B = gaussian_likelihoods(X, TRUE_MEANS, hand_covars())
    npt.assert_allclose(pred, ref_posteriors(B, HAND_STARTPROB, TRUE_TRANS),
                        rtol=1e-6, atol=1e-10)


def test_hmm_predict_gaussian_diag_against_forward_backward():
    X = gaussian_observations(seed=2, n=10)
    model = hand_gaussian("diag")
    pred = hmm_predict(X, [10], model)
    covars = [np.diag(d) for d in [[0.6, 0.7, 0.8, 0.9],
                                   [1.0, 1.1, 1.2, 1.3],
                                   [0.5, 0.5, 2.0, 2.0]]]
    B = gaussian_likelihoods(X, TRUE_MEANS, covars)
    npt.assert_allclose(pred, ref_posteriors(B, HAND_STARTPROB, TRUE_TRANS),
                        rtol=1e-6, atol=1e-10)


def test_hmm_predict_poisson_against_forward_backward():
    X = poisson_observations(seed=3, n=12)
    pred = hmm_predict(X, [12], hand_poisson())
    B = poisson_likelihoods(X, TRUE_MEANS + 0.1)
    npt.assert_allclose(pred, ref_posteriors(B, HAND_STARTPROB, TRUE_TRANS),
                        rtol=1e-6, atol=1e-10)


def test_hmm_predict_single_observation():
    # one observation: posterior is proportional to startprob * likelihood
    X = gaussian_observations(seed=4, n=1)
    pred = hmm_predict(X, [1], hand_gaussian())
    B = gaussian_likelihoods(X, TRUE_MEANS, hand_covars())[0]
    expected = HAND_STARTPROB * B / (HAND_STARTPROB * B).sum()
    npt.assert_allclose(pred[0], expected, rtol=1e-6, atol=1e-12)


def test_hmm_predict_sequences_are_decoded_independently():
    X = gaussian_observations(seed=5, n=20)
    model = hand_gaussian()
    joint = hmm_predict(X, [8, 12], model)
    first = hmm_predict(X[:8], [8], model)
    second = hmm_predict(X[8:], [12], model)
    npt.assert_allclose(joint, np.vstack([first, second]), rtol=1e-12,
                        atol=1e-15)


def test_hmm_predict_multiple_sequences_against_forward_backward():
    X = poisson_observations(seed=6, n=14)
    pred = hmm_predict(X, [5, 9], hand_poisson())
    B = poisson_likelihoods(X, TRUE_MEANS + 0.1)
    expected = np.vstack([ref_posteriors(B[:5], HAND_STARTPROB, TRUE_TRANS),
                          ref_posteriors(B[5:], HAND_STARTPROB, TRUE_TRANS)])
    npt.assert_allclose(pred, expected, rtol=1e-6, atol=1e-10)


def test_hmm_predict_lengths_must_sum_to_observations():
    X = gaussian_observations(n=10)
    with pytest.raises(ValueError, match="doesn't sum to"):
        hmm_predict(X, [4, 4], hand_gaussian())


def test_hmm_predict_empty_input_raises():
    with pytest.raises(ValueError):
        hmm_predict([], [], hand_gaussian())


def test_hmm_predict_requires_lists():
    with pytest.raises(TypeError):
        hmm_predict(np.array(gaussian_observations(n=4)), [4],
                    hand_gaussian())


# ------------------------------------
# hmm_model_save
# ------------------------------------

def test_hmm_model_save_gaussian_full_json(tmp_path):
    model = hand_gaussian()
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 10, 2, 1, 0, "gaussian")
    expected = json.dumps({"startprob": HAND_STARTPROB.tolist(),
                           "transmat": TRUE_TRANS.tolist(),
                           "means": TRUE_MEANS.tolist(),
                           "covars": hand_covars().tolist(),
                           "covariance_type": "full",
                           "n_features": 4,
                           "i_open_region": 2,
                           "i_background_region": 0,
                           "i_nucleosomal_region": 1,
                           "hmm_binsize": 10,
                           "hmm_type": "gaussian"})
    with open(f) as fh:
        assert fh.read() == expected


def test_hmm_model_save_gaussian_diag_saves_diagonal(tmp_path):
    model = hand_gaussian("diag")
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 25, 0, 2, 1, "gaussian")
    with open(f) as fh:
        m = json.load(fh)
    assert m["covariance_type"] == "diag"
    assert m["covars"] == [[0.6, 0.7, 0.8, 0.9], [1.0, 1.1, 1.2, 1.3],
                           [0.5, 0.5, 2.0, 2.0]]
    assert (m["i_open_region"], m["i_nucleosomal_region"],
            m["i_background_region"]) == (0, 2, 1)
    assert m["hmm_binsize"] == 25


def test_hmm_model_save_poisson_json(tmp_path):
    model = hand_poisson()
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 10, 2, 1, 0, "poisson")
    expected = json.dumps({"startprob": HAND_STARTPROB.tolist(),
                           "transmat": TRUE_TRANS.tolist(),
                           "lambdas": (TRUE_MEANS + 0.1).tolist(),
                           "n_features": 4,
                           "i_open_region": 2,
                           "i_background_region": 0,
                           "i_nucleosomal_region": 1,
                           "hmm_binsize": 10,
                           "hmm_type": "poisson"})
    with open(f) as fh:
        assert fh.read() == expected


@pytest.mark.parametrize("hmm_type", ["gaussian", "poisson"])
def test_hmm_model_save_trained_model_structure(tmp_path, hmm_type):
    X, lens, _ = sample_hmm(12, kind=hmm_type)
    model = hmm_training(X, lens, random_seed=5, hmm_type=hmm_type)
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 10, np.int64(2), np.int64(1), np.int64(0),
                   hmm_type)
    with open(f) as fh:
        m = json.load(fh)
    common = ["startprob", "transmat", "n_features", "i_open_region",
              "i_background_region", "i_nucleosomal_region", "hmm_binsize",
              "hmm_type"]
    if hmm_type == "gaussian":
        assert list(m.keys()) == ["startprob", "transmat", "means", "covars",
                                  "covariance_type"] + common[2:]
        assert np.array(m["means"]).shape == (3, 4)
        assert np.array(m["covars"]).shape == (3, 4, 4)
        assert m["covariance_type"] == "full"
    else:
        assert list(m.keys()) == ["startprob", "transmat", "lambdas"] + \
            common[2:]
        assert np.array(m["lambdas"]).shape == (3, 4)
    assert np.array(m["startprob"]).shape == (3,)
    assert np.array(m["transmat"]).shape == (3, 3)
    for k in ("n_features", "i_open_region", "i_background_region",
              "i_nucleosomal_region", "hmm_binsize"):
        assert type(m[k]) is int
    assert m["n_features"] == 4
    assert m["hmm_type"] == hmm_type


@pytest.mark.parametrize("covariance_type", ["tied", "spherical"])
def test_hmm_model_save_unsupported_covariance_type(tmp_path,
                                                    covariance_type):
    model = GaussianHMM(n_components=3, covariance_type=covariance_type)
    with pytest.raises(Exception,
                       match="Unknown covariance type %s" % covariance_type):
        hmm_model_save(str(tmp_path / "m.json"), model, 10, 2, 1, 0,
                       "gaussian")


def test_hmm_model_save_unknown_hmm_type_writes_nothing(tmp_path):
    """Pins the current output.

    hmm_model_save handles only 'gaussian' and 'poisson'; any other
    hmm_type writes no file and raises nothing. The command line limits
    --hmm-type to those two choices, and nothing documents what an
    unknown type should do.
    """
    f = tmp_path / "m.json"
    assert hmm_model_save(str(f), hand_gaussian(), 10, 2, 1, 0,
                          "gamma") is None
    assert not f.exists()


def test_hmm_model_save_requires_str_path(tmp_path):
    with pytest.raises(TypeError):
        hmm_model_save(tmp_path / "m.json", hand_gaussian(), 10, 2, 1, 0,
                       "gaussian")


# ------------------------------------
# hmm_model_init
# ------------------------------------

@pytest.mark.parametrize("kind", ["full", "diag", "poisson"])
def test_hmm_model_init_round_trip(tmp_path, kind):
    if kind == "poisson":
        model, hmm_type, X = hand_poisson(), "poisson", poisson_observations()
    else:
        model, hmm_type, X = (hand_gaussian(kind), "gaussian",
                              gaussian_observations())
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 20, 1, 2, 0, hmm_type)
    ret = hmm_model_init(f)
    assert isinstance(ret, list)
    assert len(ret) == 6
    loaded, i_open, i_bg, i_nuc, binsize, loaded_type = ret
    assert (i_open, i_bg, i_nuc, binsize, loaded_type) == (1, 0, 2, 20,
                                                           hmm_type)
    assert loaded.n_components == 3
    assert loaded.n_features == 4
    npt.assert_array_equal(loaded.startprob_, model.startprob_)
    npt.assert_array_equal(loaded.transmat_, model.transmat_)
    if hmm_type == "gaussian":
        assert isinstance(loaded, GaussianHMM)
        assert loaded.covariance_type == kind
        npt.assert_array_equal(loaded.means_, model.means_)
        npt.assert_array_equal(loaded.covars_, model.covars_)
    else:
        assert isinstance(loaded, PoissonHMM)
        npt.assert_array_equal(loaded.lambdas_, model.lambdas_)
    npt.assert_array_equal(hmm_predict(X, [len(X)], loaded),
                           hmm_predict(X, [len(X)], model))


@pytest.mark.parametrize("hmm_type", ["gaussian", "poisson"])
def test_hmm_model_init_trained_model_round_trip(tmp_path, hmm_type):
    X, lens, _ = sample_hmm(13, kind=hmm_type)
    model = hmm_training(X, lens, random_seed=5, hmm_type=hmm_type)
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 10, 2, 1, 0, hmm_type)
    loaded = hmm_model_init(f)[0]
    npt.assert_array_equal(hmm_predict(X, lens, loaded),
                           hmm_predict(X, lens, model))
    # a second save of the loaded model reproduces the file
    f2 = str(tmp_path / "model2.json")
    hmm_model_save(f2, loaded, 10, 2, 1, 0, hmm_type)
    with open(f) as a, open(f2) as b:
        assert a.read() == b.read()


def test_hmm_model_init_without_hmm_type_is_gaussian(tmp_path):
    # model files written before --hmm-type existed have no 'hmm_type'
    m = {"startprob": HAND_STARTPROB.tolist(),
         "transmat": TRUE_TRANS.tolist(),
         "means": TRUE_MEANS.tolist(),
         "covars": hand_covars().tolist(),
         "covariance_type": "full",
         "n_features": 4,
         "i_open_region": 2,
         "i_background_region": 0,
         "i_nucleosomal_region": 1,
         "hmm_binsize": 10}
    f = tmp_path / "old_model.json"
    f.write_text(json.dumps(m))
    loaded, i_open, i_bg, i_nuc, binsize, hmm_type = hmm_model_init(str(f))
    assert hmm_type == "gaussian"
    assert isinstance(loaded, GaussianHMM)
    assert (i_open, i_bg, i_nuc, binsize) == (2, 0, 1, 10)
    X = gaussian_observations()
    npt.assert_array_equal(hmm_predict(X, [len(X)], loaded),
                           hmm_predict(X, [len(X)], hand_gaussian()))


def test_hmm_model_init_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        hmm_model_init(str(tmp_path / "nothing.json"))


@pytest.mark.parametrize("key", ["means", "covars", "startprob",
                                 "i_open_region", "hmm_binsize"])
def test_hmm_model_init_missing_key(tmp_path, key):
    f = str(tmp_path / "model.json")
    hmm_model_save(f, hand_gaussian(), 10, 2, 1, 0, "gaussian")
    with open(f) as fh:
        m = json.load(fh)
    del m[key]
    with open(f, "w") as fh:
        json.dump(m, fh)
    with pytest.raises(KeyError):
        hmm_model_init(f)


def test_hmm_model_init_not_json(tmp_path):
    f = tmp_path / "model.json"
    f.write_text("not json\n")
    with pytest.raises(json.JSONDecodeError):
        hmm_model_init(str(f))


def test_hmm_model_init_always_builds_three_states(tmp_path):
    # The model file labels exactly three states (i_open_region,
    # i_nucleosomal_region, i_background_region) and hmmratac only trains
    # 3-state models, so hmm_model_init builds 3 components whatever the
    # file holds. A 2-state model written by hmm_model_save loads, but
    # hmmlearn rejects it when it is used.
    X, lens, _ = sample_hmm(14)
    model = hmm_training(X, lens, n_states=2, random_seed=5)
    f = str(tmp_path / "model.json")
    hmm_model_save(f, model, 10, 1, 0, 0, "gaussian")
    loaded = hmm_model_init(f)[0]
    assert loaded.n_components == 3
    assert loaded.startprob_.shape == (2,)
    with pytest.raises(ValueError,
                       match="startprob_ must have length n_components"):
        hmm_predict(X, lens, loaded)
