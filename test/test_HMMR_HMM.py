import unittest
import numpy as np
from hmmlearn import hmm
# ------------------------------------
# Main function
# ------------------------------------
''' This unittest is to check the ouputs of the hmm_training() and hmm_predict() functions
'''
def hmm_training( training_data, k):
    hmm_model = hmm.GaussianHMM( n_components = k )
    hmm_model.fit( training_data )
    return hmm_model

def hmm_predict( signals, hmm_model ):
    predictions = hmm_model.predict( signals )
    return predictions

class Test_HMM_train(unittest.TestCase):
    def setUp( self ):
        self.train_data = [[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]] 
        self.test_data = [[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]]
        self.correct_predictions = [ 0, 0, 0, 1, 1, 1 ]
        
    def test_predict( self ):
        mod = hmm_training(self.train_data, 2)
        preds = hmm_predict(self.test_data, mod)
        print(preds)
        self.assertEqual( self.correct_predictions[0], preds[0] )
        self.assertEqual( self.correct_predictions[1], preds[1] )
        self.assertEqual( self.correct_predictions[2], preds[2] )
        self.assertEqual( self.correct_predictions[3], preds[3] )
        self.assertEqual( self.correct_predictions[4], preds[4] )
        self.assertEqual( self.correct_predictions[5], preds[5] )