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

        self.startprob = [0.01807016, 0.90153727, 0.08039257]
        self.means = [[2.05560411e-01, 1.52959594e+00, 1.73568556e+00, 1.00019720e-04],
                      [1.84467806e-01, 1.46784946e+00, 1.67895745e+00, 1.00016654e-04],
                      [2.06402305e+00, 8.60140461e+00, 7.22907032e+00, 1.00847661e-04]]
        self.covars = [[[ 1.19859257e-01, 5.33746506e-02, 3.99871507e-02, 1.49805047e-07],
                        [ 5.33746506e-02, 1.88774896e+00, 7.38204761e-01, 1.70902908e-07],
                        [ 3.99871507e-02, 7.38204761e-01, 2.34175176e+00, 1.75654357e-07],
                        [ 1.49805047e-07, 1.70902908e-07, 1.75654357e-07, 1.45312288e-07]],
                       [[ 1.06135330e-01, 4.16846792e-02, 3.24447289e-02, 1.30393434e-07],
                        [ 4.16846792e-02, 1.75537103e+00, 6.70848135e-01, 1.49425940e-07],
                        [ 3.24447289e-02, 6.70848135e-01, 2.22285392e+00, 1.52914017e-07],
                        [ 1.30393434e-07, 1.49425940e-07, 1.52914017e-07, 1.27205162e-07]],
                       [[ 5.94746590e+00, 5.24388615e+00, -5.33166471e-01, -1.47228883e-06],
                        [ 5.24388615e+00, 2.63945986e+01, 3.54212739e+00, -6.03892201e-06],
                        [-5.33166471e-01, 3.54212739e+00, 1.50231166e+01, 1.43141422e-05],
                        [-1.47228883e-06, -6.03892201e-06, 1.43141422e-05, 1.04240673e-07]]]
        self.transmat =[[1.91958645e-03, 9.68166646e-01, 2.99137676e-02],
                        [8.52453717e-01, 1.46924953e-01, 6.21329356e-04],
                        [2.15432113e-02, 6.80080650e-05, 9.78388781e-01]]
        self.n_features = 4

        # for prediction
        self.prediction_data = np.loadtxt("test/small_prediction_data.txt", delimiter="\t", dtype="float", usecols=(2,3,4,5)).tolist()
        self.prediction_data_lengths = np.loadtxt('test/small_prediction_lengths.txt', dtype="int").tolist()
        self.predictions = np.loadtxt('test/small_prediction_results.txt', delimiter="\t", dtype="float").tolist()

    @pytest.mark.skip( reason="it may fail with different sklearn+hmmlearn" )
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

    @pytest.mark.skip( reason="it may fail with different sklearn+hmmlearn" )        
    def test_predict( self ):
        # test hmm_predict
        hmm_model = GaussianHMM( n_components=3, covariance_type='full' )
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

