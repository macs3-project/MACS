import unittest
import pytest
from MACS3.Signal.HMMR_HMM import *
import numpy as np
import numpy.testing as npt

# ------------------------------------
# Main function
# ------------------------------------
''' This unittest is to check the ouputs of the hmm_training() and hmm_predict() functions
'''

# @pytest.mark.skip(reason="need to refine later")
class Test_HMM_train(unittest.TestCase):
    def setUp( self ):
        self.training_data = np.loadtxt("test/large_training_data.txt", delimiter="\t", dtype="float", usecols=(2,3,4,5)).tolist()
        self.training_data_lengths = np.loadtxt('test/large_training_lengths.txt', dtype="int").tolist()
        self.expected_converged = True
        self.not_expected_covars = None
        self.not_expected_means = None
        self.not_expected_transmat = None
        self.startprob = [0.09411589, 0.82689766, 0.07898644]
        self.means = [[2.02697935e-01, 1.52785266e+00, 1.73790142e+00, 1.00019411e-04],
                      [1.87823916e-01, 1.48213364e+00, 1.69577044e+00, 1.00017125e-04],
                      [2.07360047e+00, 8.63029738e+00, 7.24406955e+00, 1.00852188e-04]]
        self.covars = [[[ 1.18061824e-01,  5.32522674e-02,  4.04981722e-02,  1.43240236e-07],
                        [ 5.32522674e-02,  1.88909221e+00,  7.44040883e-01,  1.64463390e-07],
                        [ 4.04981722e-02,  7.44040883e-01,  2.35914194e+00,  1.69079937e-07],
                        [ 1.43240236e-07,  1.64463390e-07,  1.69079937e-07,  1.38857074e-07]],

                        [[ 1.08338994e-01,  4.38027284e-02,  3.40898529e-02,  1.34873591e-07],
                        [ 4.38027284e-02,  1.78899081e+00,  6.92059837e-01,  1.54578989e-07],
                        [ 3.40898529e-02,  6.92059837e-01,  2.26836145e+00,  1.58248579e-07],
                        [ 1.34873591e-07,  1.54578989e-07,  1.58248579e-07,  1.31639696e-07]],

                        [[ 5.96438746e+00,  5.22590773e+00, -5.59954962e-01, -1.48829290e-06],
                        [ 5.22590773e+00,  2.63829229e+01,  3.49433872e+00, -6.09680431e-06],
                        [-5.59954962e-01,  3.49433872e+00,  1.50531402e+01,  1.43841972e-05],
                        [-1.48829290e-06, -6.09680431e-06,  1.43841972e-05,  1.04838987e-07]]]
        self.transmat = [[3.55718812e-03, 9.71544738e-01, 2.48980738e-02],
                         [9.22578828e-01, 7.32630014e-02, 4.15817043e-03],
                         [2.11090463e-02, 6.34703169e-04, 9.78256251e-01]]
        self.n_features = 4

        # for prediction
        self.prediction_data = np.loadtxt("test/small_prediction_data.txt", delimiter="\t", dtype="float", usecols=(2,3,4,5)).tolist()
        self.prediction_data_lengths = np.loadtxt('test/small_prediction_lengths.txt', dtype="int").tolist()
        self.predictions = np.loadtxt('test/small_prediction_results.txt', delimiter="\t", dtype="float").tolist()

    def test_training( self ):
        # test hmm_training:
        model = hmm_training(training_data = self.training_data, training_data_lengths = self.training_data_lengths, n_states = 3, random_seed = 12345, covar = 'full')
        print(model.startprob_)
        print(model.means_)
        print(model.covars_)
        print(model.transmat_)
        print(model.n_features)
        self.assertEqual( model.monitor_.converged, self.expected_converged )
        self.assertNotEqual( model.covars_.tolist(), self.not_expected_covars )
        self.assertNotEqual( model.means_.tolist(), self.not_expected_means )
        self.assertNotEqual( model.transmat_.tolist(), self.not_expected_transmat )
        npt.assert_allclose( model.startprob_.tolist(), self.startprob )
        npt.assert_allclose(model.means_, self.means)
        npt.assert_allclose(model.covars_, self.covars)
        npt.assert_allclose(model.transmat_, self.transmat)
        npt.assert_allclose(model.n_features, self.n_features)

    def test_predict( self ):
        # test hmm_predict
        hmm_model = hmm.GaussianHMM( n_components=3, covariance_type='full' )
        hmm_model.startprob_ = np.array(self.startprob)
        hmm_model.transmat_ = np.array(self.transmat)
        hmm_model.means_ = np.array(self.means)
        hmm_model.covars_ = np.array(self.covars)
        hmm_model.covariance_type = 'full'
        hmm_model.n_features = self.n_features
        predictions = hmm_predict( self.prediction_data, self.prediction_data_lengths, hmm_model )

        ## This is to write the prediction results into a file for 'correct' answer 
        #with open("test/small_prediction_results.txt","w") as f:
        #    for x,y,z in predictions:
        #        f.write( str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
        
        npt.assert_allclose( predictions, self.predictions )

