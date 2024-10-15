import unittest
# from MACS3.Signal.HMMR_HMM import *
import numpy as np
import numpy.testing as npt
from hmmlearn.hmm import PoissonHMM
# from sklearn import cluster
# import json
# class hmmlearn.hmm.PoissonHMM(n_components=1, startprob_prior=1.0, transmat_prior=1.0, lambdas_prior=0.0, lambdas_weight=0.0, 
#                       algorithm='viterbi', random_state=None, n_iter=10, tol=0.01, verbose=False, params='stl', init_params='stl', implementation='log')

# class hmmlearn.hmm.GaussianHMM(n_components=1, covariance_type='diag', min_covar=0.001, startprob_prior=1.0, transmat_prior=1.0, 
#                       means_prior=0, means_weight=0, covars_prior=0.01, covars_weight=1, algorithm='viterbi', random_state=None, n_iter=10, tol=0.01, verbose=False, 
#                       params='stmc', init_params='stmc', implementation='log')


def hmm_training(training_data, training_data_lengths, n_states=3, random_seed=12345):
    rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(random_seed)))
    hmm_model = PoissonHMM(n_components=n_states, random_state=rs, verbose=False)
    hmm_model = hmm_model.fit(training_data, training_data_lengths)
    assert hmm_model.n_features == 4
    return hmm_model


def hmm_predict(signals, lens, hmm_model):
    predictions = hmm_model.predict_proba(signals, lens)
    return predictions


class Test_HMM_train_poisson(unittest.TestCase):
    def setUp(self):
        self.training_data = np.loadtxt("test/large_training_data.txt",
                                        delimiter="\t",
                                        dtype="float",
                                        usecols=(2, 3, 4, 5)).astype(int).tolist()
        self.training_data_lengths = np.loadtxt('test/large_training_lengths.txt',
                                                dtype="int").tolist()
        self.expected_converged = True
        self.not_expected_transmat = None
        self.n_features = 4
        self.startprob = [0.83912489, 0.14996896, 0.01090615]
        self.transmat = [[9.87606722e-01, 1.23932782e-02, 1.75299652e-11],
                         [1.76603580e-02, 9.64232293e-01, 1.81073490e-02],
                         [4.87992301e-14, 2.70319349e-02, 9.72968065e-01]]
        self.lambdas = [[0.03809295, 0.62378578, 0.68739807, 0.],
                        [0.23243362, 3.4420467, 4.256037, 0.],
                        [2.58132377, 11.45924282, 8.13706237, 0.]]
        # for prediction
        self.prediction_data = np.loadtxt("test/small_prediction_data.txt",
                                          delimiter="\t",
                                          dtype="float",
                                          usecols=(2, 3, 4, 5)).astype(int).tolist()
        self.prediction_data_lengths = np.loadtxt('test/small_prediction_lengths.txt',
                                                  dtype="int").tolist()
        self.predictions = np.loadtxt('test/small_prediction_results_poisson.txt',
                                      delimiter="\t", dtype="float").tolist()

    def test_training(self):
        # test hmm_training:
        model = hmm_training(training_data=self.training_data,
                             training_data_lengths=self.training_data_lengths,
                             n_states=3,
                             random_seed=12345)
        # print(model.startprob_)
        # print(model.transmat_)
        # print(model.lambdas_)
        # print(model.n_features)
        self.assertEqual(model.monitor_.converged, self.expected_converged)
        self.assertNotEqual(model.transmat_.tolist(), self.not_expected_transmat)
        npt.assert_allclose(model.startprob_.tolist(), self.startprob)
        npt.assert_allclose(model.transmat_, self.transmat)
        npt.assert_allclose(model.lambdas_, self.lambdas)
        npt.assert_allclose(model.n_features, self.n_features)

    def test_predict(self):
        # test hmm_predict
        hmm_model = PoissonHMM(n_components=3)
        hmm_model.startprob_ = np.array(self.startprob)
        hmm_model.transmat_ = np.array(self.transmat)
        hmm_model.lambdas_ = np.array(self.lambdas)
        hmm_model.n_features = self.n_features
        predictions = hmm_predict(self.prediction_data,
                                  self.prediction_data_lengths,
                                  hmm_model)

        # This is to write the prediction results into a file for 'correct' answer
        # with open("test/small_prediction_results_poisson.txt","w") as f:
        #    for x,y,z in predictions:
        #        f.write(str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
        npt.assert_allclose(predictions, self.predictions)
