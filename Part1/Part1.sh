#!/bin/bash

#########################################
# Make sure Part1.sh, getWeights.py, 
# ManipulateData.py and getData.py are 
# in the same directory.
# Your environment should also be set with
# Python, BEDtools and SAMtools. 
#########################################

#########################################
# Inputs:
# (a) a path to a directory holding all 
# the controls -- required.
#
# (b) full path to the treatment file 
# -- required.
#
# (c) full path to file containing 
# chromosome sizes corresponding to the 
# species genome. 
# File format eg:
# chr1	248956422
# If no chromosome size is provided, by 
# default chromosome sizes corresponding 
# to the human genome will be used. 
# -- optional.
##########################################

#########################################
# The BAM files could either be 
# 1. unsorted or unindexed bam files.
# 2. sorted and indexed bam files.
##########################################

controlDir=$1 #Argument 1
chipFile=$2 #Argument 2
arg3=$3 #Argument 3 for chrom sizes

if [ -z "$1" ]
then
	echo "Control directory not supplied"
	exit 1
fi

if [ -z "$2" ]
then
    echo "ChIP treatment file not supplied"
	exit 1
fi

if [[ ! -d $controlDir ]]
then 
	echo "The first argument should be a directory of the control samples."
	exit 1
fi

if [[ ! -f $chipFile ]]
then
	echo "The second argument should be a file for the treatment sample."
	exit 1
fi

if [[ ! -f "$3" ]]
then
	arg3=""
else
	chrom=$3
fi

#########################################
# Create directories 
#########################################
mainDir=WACS
controlDirFiltered=WACS/Control
controlReadCount=WACS/Control/ReadCount
chipDir=WACS/ChIP

if [[ ! -d  "$mainDir" ]]
then
	mkdir WACS
fi

if [[ ! -d  "$controlDirFiltered" ]]
then
	mkdir WACS/Control
	mkdir WACS/Control/ReadCount
fi


if [[ ! -d  "$chipDir" ]]
then
	mkdir WACS/ChIP
fi


#########################################
# Get chromosome sizes
#########################################

if [ -z "$arg3" ]
then
	if [ ! -f $mainDir/hg38.chrom.sizes ]
	then 
		wget https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -P $mainDir
		chrom=WACS/hg38.chrom.sizes 
	fi
fi

if [ ! -f $mainDir/intervalFile.bed ]
then
	echo "Create interval file with 200 bp windows."
	bedtools makewindows -g $chrom -w 200 > $mainDir/intervalFile.bed
fi
#########################################
# Process treatment sample.
#########################################
echo "Process treatment sample."
FILENAME=`echo $chipFile | rev | cut -d '/' -f 1 | rev`
BASE=`echo $FILENAME | cut -d \. -f 1`
chipReadCount=$chipDir/$BASE.readcount.txt

if [[ $chipFile == *.bam ]]
then
	if [ -f "$chipFile.bai" ]
	then
		echo "File indexed."
	else
		echo "Sort $chipFile"
		samtools sort $chipFile > $chipFile
		echo "Index $FILENAME"
		samtools index $chipFile
	fi

	if [ -s "$chipDir/$BASE.bam" ]
	then
		echo "Filtered file exists."
	else
		echo "Remove unwanted chromosomes from $FILENAME"
		samtools idxstats $chipFile | cut -f 1 | grep -v Un | grep -v random | grep -v M | grep -v E | xargs samtools view -b $chipFile > $chipDir/$BASE.bam
	fi 
	
	if [ -s "$chipDir/$BASE.readcount.txt" ]
	then
		echo "$BASE.readcount.txt file exists and not empty."
	else
		echo "Produce read counts per 200bp windows along genome. Saved in format name.readcount.txt."
		bedtools coverage -a $mainDir/intervalFile.bed -b $chipDir/$BASE.bam -counts > $chipDir/$BASE.readcount.txt
	fi
else
	echo "Incorrect file format. Make sure the file is in BAM format."
	exit 1
fi

#########################################
# Process the control samples. 
#########################################
echo Process the controls in $controlDir
for file in $controlDir*
do
	if [[ $file == *.bam ]]
	then
		FILENAME=`echo $file | rev | cut -d '/' -f 1 | rev`
		BASE=`echo $FILENAME | cut -d \. -f 1`	

		if [ -f "$controlDir$BASE.bam.bai" ]
		then
			echo "File indexed."
		else 
			echo Sort $controlDir$FILENAME
			samtools sort $controlDir$BASE.bam > $controlDir$BASE.bam
			echo Index $controlDir$FILENAME
			samtools index $controlDir$BASE.bam
			fi	

		if [ -s "$controlDirFiltered/$FILENAME" ]
		then
			echo "Filtered file exists."
		else
			echo Remove unwanted chromosomes from $controlDir$FILENAME
			samtools idxstats $controlDir$FILENAME | cut -f 1 | grep -v Un | grep -v random | grep -v M | grep -v E | xargs samtools view -b $controlDir$FILENAME > $controlDirFiltered/$FILENAME 
		fi
			
		if [ -s "$controlReadCount/$BASE.readcount.txt" ]
		then
			echo "$BASE.readcount.txt file exists."
		else
			echo Produce read counts per 200bp windows along genome. Saved in format $BASE.readcount.txt
			bedtools coverage -a $mainDir/intervalFile.bed -b $controlDirFiltered/$FILENAME -counts > $controlReadCount/$BASE.readcount.txt
		fi	
	fi
done

#########################################
# Create control matrix.
########################################## 
echo "Move python files into WACS directory"
if [ -f "getWeights.py" ]
then 
	mv getWeights.py $mainDir
else
	echo "Make sure Part1.sh, getWeights.py, ManipulateData.py and getData.py are in the same directory."
fi
if [ -f "ManipulateData.py" ]
then 
	mv ManipulateData.py $mainDir
else
	echo "Make sure Part1.sh, getWeights.py, ManipulateData.py and getData.py are in the same directory."
fi
if [ -f "getData.py" ]
then 
	mv getData.py $mainDir
else
	echo "Make sure Part1.sh, getWeights.py, ManipulateData.py and getData.py are in the same directory."
fi

cd $mainDir
controlReadCount=Control/ReadCount
echo "Make control matrix and save as outputControl.csv"
python getData.py $controlReadCount outputControl.csv

#########################################
# Get the weights.
########################################## 
echo "Get weights."
echo $chipReadCount
chipReadCount=ChIP/*.readcount.txt
nrows=$(< $chipReadCount wc -l)
python getWeights.py Data/outputControl.csv $chipReadCount -s True --nrows $nrows

