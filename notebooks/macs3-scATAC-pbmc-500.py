#!/usr/bin/env python
# coding: utf-8

# # Use MACS3 API to call peaks for each cluster in scATAC-seq data

# ## Load the fragment file and build pileup track

# We will import the fragment file parser from MACS3.

# In[1]:


# %pip install macs3 panda

from MACS3.IO.Parser import FragParser


# Assume the fragment file is at `../test/atac_pbmc_500_v1_fragments.tsv.gz`. Note that MACS3 can directly load gzipped files. We will load the file and build a `PairEndTrack.PETrackII` object, which contains alignment locations, barcodes and counts.

# In[2]:


frag_file = FragParser("../test/atac_pbmc_500_v1_fragments.tsv.gz", buffer_size=100000)
petrack = frag_file.build_petrack(max_count=2)
petrack.finalize()


# Now we can check some basic statistics of the loaded data.

# In[3]:


print(f"average template length is {petrack.average_template_length}")
print(f"total number of bases is {petrack.total}")
print(petrack.get_chr_names())


# ## Call peaks for the entire dataset

# Next, we will call peaks for the entire dataset. The first step is to build a pileup track from the alignments. But before that, we can filter out invalid cells/barcodes. For example, we have a `atac_pbmc_500_v1_singlecell.csv` file downloaded together with the fragment file which is from the standard 10x pipeline. We want to only keep the cell barcodes with at least 500 usable fragments and at least 25% of total fragments are marked as usable.

# In[4]:


import pandas as pd

barcodes = []

# Read the CSV
df = pd.read_csv("../test/atac_pbmc_500_v1_singlecell.csv")

# We keep the cell barcode that 1) is a cell barcode 2) more than 500 usable fragments 3) at least 25% of total fragments are usable
df_pass = df[
    (df["is__cell_barcode"] == 1) & 
    (df["passed_filters"] >= 500) &
    (df["passed_filters"]/df["total"] >= 0.25)]

# Extract barcode values as a list then convert to a set
# Note: we need to encode the string since PETrackII.subset needs bytestrings.
barcodes = set([x.encode() for x in df_pass["barcode"].tolist()])

print("Cell barcodes to be kept:", len(barcodes))

petrack = petrack.subset(barcodes)


# Now we can build the pileup signal track of all valid cell barcodes and their fragments.

# In[5]:


pileup_track = petrack.pileup_bdg()


# We can get the sum, total length, maximum, minimum, mean, and variance from the pileup track by using `summary` function.

# In[6]:


(pileup_sum, pileup_length, pileup_max, pileup_min, pileup_mean, pileup_var) = pileup_track.summary()
print(f"{pileup_sum=}\n{pileup_length=}\n{pileup_max=}\n{pileup_min=}\n{pileup_mean=}\n{pileup_var=}")


# Next, we will call peaks for the entire dataset. Since we are calling peaks for ATAC-seq data, we will use the single whole-genome average pileup value as the background. From above calculation, the average is in the variable `pileup_mean`.

# In[7]:


global_bg_track = pileup_track.set_single_value(pileup_mean) # this will return a new track with all values set as `pileup_mean`


# We construct a score track for comparing observed pileup and background. The score track will contain the observed pileup, the background, and scores for each position in the genome.

# In[8]:


score_track = pileup_track.make_ScoreTrackII_for_macs(global_bg_track, depth1=100, depth2=100) # Note, depth1 and 2 are the same so there is no need for scaling the values


# We will use q-score (-log10 q-value) as scores. MACS3 supports the following methods:
# 
# - p: -log10 pvalue;
# - q: -log10 qvalue;
# - l: log10 likelihood ratio (minus for depletion)
# - s: symmetric log10 likelihood ratio (for comparing two ChIPs)
# - f: log10 fold enrichment
# - F: linear fold enrichment
# - d: subtraction
# - M: maximum
# - m: fragment pileup per million reads
# 
# To use any of the methods, provide the argument to `change_score_method` with `ord`, as shown in the following example.

# In[9]:
score_track.change_score_method(ord('p'))


# In[10]:


a = score_track.get_data_by_chr(b'chr1')[3]


# In[11]:

import numpy as np

print(np.max(a))
print(np.min(a))
