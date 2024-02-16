import unittest
import pytest
from MACS3.Signal.HMMR_HMM import hmm_training, hmm_predict
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
        self.train_data = np.loadtxt("test/large_training_data.txt", delimiter="\t", dtype="float", usecols=(2,3,4,5)).tolist()
        self.training_data_lengths = np.loadtxt('test/large_training_lengths.txt', dtype="int").tolist()
        self.expected_converged = True
        self.not_expected_covars = None
        self.not_expected_means = None
        self.not_expected_transmat = None
        self.means = [[2.05560411e-01, 1.52959594e+00, 1.73568556e+00, 1.00019720e-04],
                      [1.84467806e-01, 1.46784946e+00, 1.67895745e+00, 1.00016654e-04],
                      [2.06402305e+00, 8.60140461e+00, 7.22907032e+00, 1.00847661e-04]]
        self.covars = [[[ 1.19859257e-01,  5.33746506e-02,  3.99871507e-02,  1.49805047e-07],
                        [ 5.33746506e-02,  1.88774896e+00,  7.38204761e-01,  1.70902908e-07],
                        [ 3.99871507e-02,  7.38204761e-01,  2.34175176e+00,  1.75654357e-07],
                        [ 1.49805047e-07,  1.70902908e-07,  1.75654357e-07,  1.45312288e-07]],

                        [[ 1.06135330e-01,  4.16846792e-02,  3.24447289e-02,  1.30393434e-07],
                        [ 4.16846792e-02,  1.75537103e+00,  6.70848135e-01,  1.49425940e-07],
                        [ 3.24447289e-02,  6.70848135e-01,  2.22285392e+00,  1.52914017e-07],
                        [ 1.30393434e-07,  1.49425940e-07,  1.52914017e-07,  1.27205162e-07]],

                        [[ 5.94746590e+00,  5.24388615e+00, -5.33166471e-01, -1.47228883e-06],
                        [ 5.24388615e+00,  2.63945986e+01,  3.54212739e+00, -6.03892201e-06],
                        [-5.33166471e-01,  3.54212739e+00,  1.50231166e+01,  1.43141422e-05],
                        [-1.47228883e-06, -6.03892201e-06,  1.43141422e-05,  1.04240673e-07]]]
        self.transmat = [[1.91958645e-03, 9.68166646e-01, 2.99137676e-02],
                         [8.52453717e-01, 1.46924953e-01, 6.21329356e-04],
                         [2.15432113e-02, 6.80080650e-05, 9.78388781e-01]]
        # self.predictions = np.loadtxt("test/standard_results_hmmratac/predictions.txt", delimiter="\t", dtype="float").tolist()
        # self.cr_data = np.loadtxt("test/standard_results_hmmratac/cr_data.txt", delimiter="\t", dtype="float").tolist()
        # self.cr_data_lengths = np.loadtxt("test/standard_results_hmmratac/cr_data_lengths.txt", dtype="float").astype(np.int64).tolist()



    def test_predict( self ):
        model = hmm_training(training_data = self.train_data, training_data_lengths = self.training_data_lengths, n_states = 3, random_seed = 12345)

        # predictions = hmm_predict(self.cr_data, self.cr_data_lengths, model)
        # predictions = hmm_predict(self.train_data, self.training_data_lengths, model)
        # np.savetxt('test/standard_results_hmmratac/predictions.txt', predictions, delimiter="\t", fmt="%f")

        # test training:
        self.assertEqual( model.monitor_.converged, self.expected_converged )
        self.assertNotEqual( model.covars_.tolist(), self.not_expected_covars )
        self.assertNotEqual( model.means_.tolist(), self.not_expected_means )
        self.assertNotEqual( model.transmat_.tolist(), self.not_expected_transmat )
        npt.assert_allclose(model.means_, self.means)
        npt.assert_allclose(model.covars_, self.covars)
        npt.assert_allclose(model.transmat_, self.transmat)
        
        # test predicition:
        # npt.assert_allclose(predictions, self.predictions) #adjust rtol, atol for tolerance