import unittest
import pytest
import numpy as np
from hmmlearn import hmm
from MACS3.Signal.HMMR_HMM import hmm_training, hmm_predict
import numpy as np
import pandas as pd
# ------------------------------------
# Main function
# ------------------------------------
''' This unittest is to check the ouputs of the hmm_training() and hmm_predict() functions
'''

# @pytest.mark.skip(reason="need to refine later")
class Test_HMM_train(unittest.TestCase):
    def setUp( self ):
        train_data = pd.read_csv('test/large_training/large_training_data.txt', sep='\t', names=['a', 'b', 'c', 'd', 'e', 'f'])
        self.train_data = train_data[['b', 'c', 'd', 'e', 'f']].to_numpy().tolist()
        self.training_data_lengths = np.loadtxt('test/large_training/large_training_lengths.txt', dtype=int).tolist()
        self.expected_converged = True
        self.not_expected_covars = None
        self.not_expected_means = None
        self.not_expected_transmat = None

    def test_predict( self ):
        model = hmm_training(training_data = self.train_data, training_data_lengths = self.training_data_lengths, n_states = 3, random_seed = 12345)

        self.assertEqual( model.monitor_.converged, self.expected_converged )
        self.assertNotEqual( model.covars_.tolist(), self.not_expected_covars )
        self.assertNotEqual( model.means_.tolist(), self.not_expected_means )
        self.assertNotEqual( model.transmat_.tolist(), self.not_expected_transmat )