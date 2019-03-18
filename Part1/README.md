# WACS for weighted peak calling approach. 

## Introduction 
Chromatin immunoprecipitation followed by high throughput sequencing (ChIP-seq) allows biologists to identify protein/DNA binding
and histone modifications across the genome. A key step in the analysis of ChIP-seq data is peak calling, which is used to identify 
regions of enrichment. The incorporation of controls is used to account for the bias in the data. However, a recurrent issue is 
the existence of different types of bias in different ChIP-seq experiments. Peak calling can generate different results for the 
same ChIP-seq experiment depending on the controls used. The present study proposed  We introduce WACS (Weighted Analysis of ChIP-seq) 
which is an extension of the well-known peak caller, MACS2. It allows the use of “smart” controls to model the non-signal effect for 
a specific ChIP-seq experiment. WACS first generates the weights per control to model the background distribution per ChIP-seq experiment 
(weights are generated in Part 1). This is then followed by peak calling. 

## Install

- Bioinformatics:
	- SAMtools (http://www.htslib.org/download/)
	- bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
- Python 2.7:
	- pip install numpy 
	- pip install pandas
	- pip install scipy
	
# Make sure Part1.sh, getWeights.py, ManipulateData.py and getData.py are in the same directory.
# Your environment should also be set with Python, BEDtools and SAMtools. 
	
## Usage
#1. To generate the weights per control sample:

./Part.sh controlDir treatmentFile chromSize

Inputs:
(a) a path to a directory holding all the control BAM files -- required.

(b) full path to the treatment BAM file -- required.

(c) full path to file containing chromosome sizes corresponding to the species genome.  -- optional.
File format eg:
chr1	248956422
If no chromosome size is provided, by default chromosome sizes corresponding to the human genome will be used. 

# The BAM files could either be (1) unsorted or unindexed bam files, or (2) sorted and indexed bam files.

#2. Call peaks
For each chip, get control names and their corresponding weights and pass them to macs2 callpeak_wacs.

FILENAME=`basename name.bam`
BASE=`echo $FILENAME | cut -d. -f1`
coefFile=Coefficients/$BASE.coefficients.csv 
ControlBamDir=BAM/Control

if [ -e $coefFile ]
then
	controlNames="`cut -d, -f1 $coefFile`"
	controlWeights="`cut -d, -f2 $coefFile`"

	for i in $controlNames; do controlNames_full+="$ControlBamDir$i.bam "; done

	macs2 callpeak_wacs -t name.bam -c $controlNames_full -w $controlWeights
fi 