"""Module Description
@author:  Aseel Awdeh
@contact: araed104@uottawa.ca
"""
# python modules
import os
import numpy as np
import pandas as pd
import sys
from optparse import OptionParser

###############################################################################################################
# getData.py
#
# Takes a set of control files in the form of read counts (counts per window across genome)
# Combines these files to produce a matrix of controls -- columns are the controls and the rows are the windows.
#
###############################################################################################################

#Create a directory if not exists
sample_dir = 'Data'
if not os.path.exists(sample_dir):
   os.makedirs(sample_dir)

################################################################################
# main
################################################################################
def main():
    print('Combine control read counts.')
    usage = 'usage: %prog [options] <directory_readcount_files> <outputFile_name>'
    parser = OptionParser(usage)
    parser.add_option('-c', '--control_dir', dest='Control', default=None, help='Directory of control read count files.')
    parser.add_option('-o', '--outputFile', dest='Name', default=None, help='Filename for matrix of control read counts per window across the genome.')
    (options,args) = parser.parse_args()
    print(options)
    print(args)

    if len(args) != 2:
    	parser.error('Must provide directory to control read counts and name of output file.')
    else:
        controlDir = args[0]
        outputFile = args[1]
        getControls(controlDir, outputFile)

###############################################################################################################
# getControls -- Takes a path
#
# Input: directory -- pathway to folder of read count files
#
# Output: Output file name
###############################################################################################################
def getControls(directory, outfile_name, nrows=None):
    """
        Get read count per window from directories of controls.
        Need to loop through directories and fill out matrix with read counts for controls

        Matrix = (nsamples, nfeatures) = (windows, controls)

        Keyword arguments:
        directories -- directories of readcounts.txt
        outfile_name -- name of output file
        nrows -- how many rows to consider

        Goal: Save all controls in one csv file for easier access in outfile_name.
    """

    #Data frame to add columns to array
    df_samples = pd.DataFrame()

    #Get list of files on directory -- eg. control files
    dir_Files = os.listdir(directory)
    print("Combine " + str(len(dir_Files)) +" readcount files.")
    #-----------------Loop through directory -------------------------------
    ########################################################################
    # Loop through all the files
    for file in dir_Files:
        #If not a directory -- needs to be readcount file
        if not os.path.isdir(file):
            print(file)
            f_path = os.path.join(directory, file)

            #Read each file
            print("File " + f_path)
            df = pd.read_csv(f_path, sep='\t',names=['Chr', 'start', 'end', 'count'])#, nrows=nrows)#,  dtype={"Chr": object, "count": int})

            #For regression, add controls as cols not rows, to dataframe
            nameBAM = file.split(".")[0]
            df_samples[nameBAM] = df['count']

    df_samples.to_csv(os.path.join(sample_dir, outfile_name), index=False) #write into file format: counts of control

    return df_samples

def getChIP(chipPath, rows):
    """
        Get read count per window for chipseq readcount dataset

        Matrix = (nsamples, nfeatures) = (windows, chip)

        Keyword arguments:
        chipPath -- Path for chipseq dataset
        Rows -- how many rows to consider

        Goal: Returns dataframe of counts per window for chip dataset
    """

    df_chipseq = pd.read_csv(chipPath, sep='\t',names=['Chr', 'start', 'end', 'count'], nrows=rows)
    df = pd.DataFrame()
    #NAME.readcount.txt -- get name
    nameChip = chipPath.split("/")
    nameChip = nameChip[len(nameChip) - 1].split(".")[0]
    df[nameChip] = df_chipseq['count']

    return df

###############################################################################################################
#main
###############################################################################################################
if __name__ == "__main__":
    main()
