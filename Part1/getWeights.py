"""Module Description
@author:  Aseel Awdeh
@contact: araed104@uottawa.ca
"""
# python modules
import os
import sys
import pandas as pd
import numpy as np
from optparse import OptionParser
from scipy.optimize import nnls

from ManipulateData import ManipulateData

###############################################################################################################
# getWeights.py
#
# Takes as input the read count text files for each of the treatment and control samples.
# And produces a weight for each control sample.
###############################################################################################################
#Create a directory if not exists
coeff_dir = 'Coefficients'
if not os.path.exists(coeff_dir):
   os.makedirs(coeff_dir)

mse_dir = 'MSE'
if not os.path.exists(mse_dir):
   os.makedirs(mse_dir)

###############################################################################################################
# main
###############################################################################################################
def main():
    print('Generate weights for each control.')
    usage = 'usage: %prog [options] <control_readcount> <treatment_readcount> <nrows> <sample>'
    parser = OptionParser(usage)
    parser.add_option('-c', '--controlMatrix', help='Matrix of control read counts per window across the genome. Format: columns are the controls, rows are the windows.')
    parser.add_option('-t', '--treatment_readcount', help='Read count file for a treatment sample.')
    parser.add_option('--nrows', default=0, help='Number of rows (windows) in read count text files.')
    parser.add_option('-s', '--sample', default=False, help='Sample rows by taking every third row.')

    (options,args) = parser.parse_args()

    if len(args) != 2 and options.nrows != 0:
    	parser.error('Must provide directory to control read counts, name of output file and number of rows in text files.')
    else:
        controlMatrixcsv = args[0]
        readcountsChipFile = args[1]
        nrows = int(options.nrows)
        sampleBool = options.sample

        chipname = readcountsChipFile.split('/')
        chipname = chipname[len(chipname) - 1]
        chipname = chipname.split('.')[0]
        print('Treatment sample: ' + chipname)

        m = ManipulateData(controlMatrixcsv, sample=sampleBool, nrows=nrows)
        m.readperbillionNormalizeControl()
        runChIP(chipname, readcountsChipFile, m)

###############################################################################################################
# runChIP
#
# Input
#
# Output
###############################################################################################################
def runChIP(file, readcountsChipFile, m):

    if not os.path.isdir(file):
        chipname = file.split('.')[0]
        outputCoefDir = os.path.join(coeff_dir, chipname + '.coefficients.csv')
        outputMSEDir = os.path.join(mse_dir, chipname + '.mse.txt')

        df_chip = m.getChIP(readcountsChipFile)
        print("Normalize ChIP.")
        df_chip = m.readperbillionNormalizeChip(df_chip)

        #Regression
        print("Run regression.")
        coefficients, avg_mse_test, avg_mse_train = m.nnls_regression(df_chip)
        controlNames = m.getControlNames()

        print("Writing coefficents.")
        w = open(outputCoefDir, 'w')
        for i in range(len(coefficients)):
           w.write(controlNames[i] + ', ' + str(coefficients[i]) + '\n')
        w.close()

        print("Writing MSE.")
        w = open(outputMSEDir, 'w')
        w.write(chipname + "," + str(avg_mse_test))
        w.close()

###############################################################################################################
# main
###############################################################################################################
if __name__ == "__main__":
    main()
