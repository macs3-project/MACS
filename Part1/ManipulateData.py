"""Module Description
@author:  Aseel Awdeh
@contact: araed104@uottawa.ca
"""
# python modules
from __future__ import division
import numpy as np
import pandas as pd
import random
import logging
import os
import sys
import scipy
from scipy import stats
from scipy.optimize import nnls, leastsq
from scipy.linalg import solve
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold
import matplotlib


###############################################################################################################
# ManipulateData.py
#
# Takes as input the controlmatrix as well as the chipseq read count text file.
# Reads them (all or samples from the files)
# Normalizes the data
# Applies nnls to the data
# Returns coefficents and mse
###############################################################################################################
class ManipulateData:
    """
        Input: path to a ChIPseq readcount file.
        Dealing with filtered bam files, where counts per 200bp window with 50bp increments along genome found.

        Step 1: Retrive matrix of controls (nwindows, ncontrols) and chipseq dataset.
        Step 2: Normalize each control to chipseq, by finding total count in control and total count in chipseq and finding ratio
        Step 3: Find weights associated with each control by finding variance of each control. -- W=1/V
        Step 4: Use non negative weighed least squares to find value of coefficients for each control to model background signal for each chipseq.

        chipseq and control files are converted to read counts per 200bp window with 50 bp increments across genome.
        control is a (window, control) matrix

        Attributes:
        ----------
            control_csv_path (str) -- path to csv files containing the counts per window (rows) for each control (column)
            nrows (int) -- number of instances to consider to determine weights per control
            sample (boolean) -- if true, then sample

    """

    def __init__(self, control_csv, nrows = None, control_names = None, sample = None):
        print("#0 Reading control csv file " + str(control_csv))
        self.sample = sample
        self.nrows = nrows
        self.control_names = control_names

        if sample and not self.control_names:
            #https://nikolaygrozev.wordpress.com/2015/06/16/fast-and-simple-sampling-in-pandas-when-loading-data-from-files/
            #Take every N-th row
            #print("Sample from control read count matrix.")
            n = 3
            # Count the lines or use an upper bound
            num_lines = self.nrows
            # The row indices to skip - make sure 0 is not included to keep the header!
            skip_idx = [x for x in range(1, num_lines) if x % n != 0]
            # Read the data
            self.df_control = pd.read_csv(control_csv, sep=',', skiprows=skip_idx)

        elif self.control_names and sample:
            #print("Sample and controls selected.")
            print(self.control_names)
            n = 3
            num_lines = self.nrows#61765409 #sum(1 for l in open(f))
            skip_idx = [x for x in range(1, num_lines) if x % n != 0]
            self.df_control = pd.read_csv(control_csv, sep=',', skiprows=skip_idx, usecols=self.control_names)#, nrows=10000)

        elif self.control_names and self.nrows:
            self.df_control = pd.read_csv(control_csv, sep=',' ,usecols=self.control_names, nrows=self.nrows) #dtype=np.float32, call if control matrix already saved

        elif self.control_names and not self.nrows:
            self.df_control = pd.read_csv(control_csv, sep=',' ,usecols=self.control_names) #dtype=np.float32, call if control matrix already saved

        elif not self.control_names and self.nrows:
            self.df_control = pd.read_csv(control_csv, sep=',' ,nrows=self.nrows) #dtype=np.float32, call if control matrix already saved

        else:
            self.df_control = pd.read_csv(control_csv, sep=',')


    ###############################################################################################################
    # getChIP
    #
    ###############################################################################################################
    def getChIP(self, path):
        """
            Get read count per window for chipseq readcount dataset

            Matrix = (nsamples, nfeatures) = (windows, chip)

            Parameters:
            -----------
            chipPath -- Path for chipseq dataset

            Returns:
            ---------
            Returns dataframe of counts per window for chip dataset
        """

        print("#1 Reading ChIP file")
        if self.sample:
            n = 3
            # Count the lines or use an upper bound
            num_lines = self.nrows#61765409 #sum(1 for l in open(f))
            # The row indices to skip - make sure 0 is not included to keep the header!
            skip_idx = [x for x in range(1, num_lines) if x % n != 0]
            # Read the data
            df_chipseq = pd.read_csv(path, sep='\t', names=['Chr', 'start', 'end', 'count'], skiprows=skip_idx)#, nrows=10000)
        elif self.nrows:
            df_chipseq = pd.read_csv(path, sep='\t', names=['Chr', 'start', 'end', 'count'], nrows=self.nrows)
        else:
            df_chipseq = pd.read_csv(path, sep='\t', names=['Chr', 'start', 'end', 'count'])

        df = pd.DataFrame()
        #NAME.readcount.txt -- get name
        nameChip = path.split("/")
        nameChip = nameChip[len(nameChip) - 1].split(".")[0]
        df[nameChip] = df_chipseq['count']

        return df

    def getControlMatrix(self):
        return self.df_control

    def getControlNames(self):
        return self.df_control.columns.values

    ###############################################################################################################
    # readperbillionNormalizeControl
    #
    ###############################################################################################################
    def readperbillionNormalizeControl(self):
        """
            This is to normalize the read counts for sequencing depth and gene length.
            - Find sum of each control (each column) -Sc
            - Divide Sc/10^9 -- this is called the per million scaling factor -- PMSF
            - Divide each read count in column (control) by PMSF -- RMSF

        """
        print("#0 Read ber billion for Control RPB")
        sum_column = self.df_control.sum(axis=0)
        sum_column = sum_column.as_matrix()
        PMSF = sum_column/1000000000.0
        RPMSF = self.df_control/PMSF
        RPMSF = RPMSF.fillna(0) #if sum of column == 0, then dividing by zero will lead to nan.
        self.df_control = RPMSF

    ###############################################################################################################
    # readperbillionNormalizeChip
    #
    ###############################################################################################################
    def readperbillionNormalizeChip(self, df):
        """
            This is to normalize the read counts for sequencing depth and gene length.
            - Find sum of each control (each column) -Sc
            - Divide Sc/10^9 -- this is called the per million scaling factor -- PMSF
            - Divide each read count in column (control) by PMSF -- RMSF

        """
        print("#1 Read ber billion for ChIP RPB")
        sum_column = df.sum(axis=0)
        sum_column = sum_column.as_matrix()
        PMSF = sum_column/1000000000.0
        RPMSF = df/PMSF
        RPMSF = RPMSF.fillna(0) #if sum of column == 0, then dividing by zero will lead to nan.
        return RPMSF

    ###############################################################################################################
    # Non negative linear regression
    #
    # Input
    ###############################################################################################################
    def nnls_regression(self, y):
        """
            Non negative least square regression:

            y = (instances, 1) -- b target
            X = (instances, features) -- A

            The ordinary least squares minimization function is:
                argmin (Aix - bi)^2
        """
        print("#3 Non negative linear regression")

        y = y.transpose().values.ravel()
        X = self.df_control.as_matrix()

        #Cross validation
        training_error_fold = []
        testing_error_fold = []
        training_r2_fold = []
        kf = KFold(n_splits = 10, shuffle=True)
        for train_index , test_index in kf.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            #Train
            model, rnorm = nnls(X_train, y_train)
            #training error
            y_train_predicted = X_train.dot(model)
            mse_train = mean_squared_error(y_train, y_train_predicted)
            training_error_fold.append(mse_train)

            #Predict on testing test
            y_predicted = X_test.dot(model)
            #Calculate the training error -- compared predicted to actual for test (y_test)
            mse_test = mean_squared_error(y_test, y_predicted)
            testing_error_fold.append(mse_test)

        avg_mse_train = np.array(training_error_fold).mean()
        avg_mse_test = np.array(testing_error_fold).mean()

        #Fit to whole data
        coeff, rnorm = nnls(X, y)
        return coeff, avg_mse_test, avg_mse_train
