import unittest
import pytest
import numpy as np
from hmmlearn import hmm
from MACS3.Signal.HMMR_HMM import hmm_training, hmm_predict
# ------------------------------------
# Main function
# ------------------------------------
''' This unittest is to check the ouputs of the hmm_training() and hmm_predict() functions
'''
# def hmm_training( training_data, k):
#     hmm_model = hmm.GaussianHMM( n_components = k )
#     hmm_model.fit( training_data )
#     return hmm_model

# def hmm_predict( signals, hmm_model ):
#     predictions = hmm_model.predict( signals )
#     return predictions

class Test_HMM_train(unittest.TestCase):
    def setUp( self ):
        self.train_data = [[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]] 
        self.test_data = [[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]]
        self.model = hmm_training(training_data = self.train_data, training_data_lengths = [6], n_states = 2)
        self.pred = hmm_predict(signals = self.test_data, lens = [6], hmm_model = self.model)

    #@pytest.mark.skip(reason="state label is random, so we need to fix this test")
    def test_predict( self ):
        mod = hmm_training(self.train_data, [6], 2)
        preds = hmm_predict(self.test_data, [6], mod)
        # print(preds)
        # print(self.pred)
        self.assertEqual( self.pred.tolist(), preds.tolist() )
        #self.assertAlmostEqual()

