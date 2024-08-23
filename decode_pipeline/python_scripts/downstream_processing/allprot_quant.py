#!/bin/usr python

import os
import pandas as pd
from Bio import SeqIO
import numpy as np
from itertools import groupby
import re
from operator import itemgetter
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from glob import glob
from collections import Counter
import matplotlib as mpl
from matplotlib.lines import Line2D
import random
import copy
import sys

"""
This scripts should be run to obtain normalized protein abundance of all main proteins across patients. The data is used in downstream analyses (e.g. PSEA).

Input:
    1. dataset info (directory, samples, sample_map, output directory)
    2. Validation_search_evidence_dict.p (output of AA_subs_validation2.py)
    3. proteinGroups.txt (MQ output files from validation search).

The MQ results used in this script are from the MQ validation search
    1. read in proteinGroups.txt from each tmt set for quant data
    2. Quantify precursor level abundance of each protein in each label-free sample or TMT set.
    3. Quantify reporter-level abundance of protein in each sample by distributing the precursor intensity by the relative intensity of the reporter ions and normalizing for sample loading (TMT data only)
    4. Alternative normalization by median normalization followed by normalization to reference channel

Returns a dictionary with protein abundances and a table of proteins x samples abundances for precursor and reporter-ion level.
"""

ds = 'LSCC'
proj_dir = '/scratch/tsour.s/Sat_LSCC/'
aa_subs_dir = proj_dir + 'AA_subs_pipeline/'

if ds!='Healthy':
    sample_map = sample_map_list[datasets.index(ds)]
    ref_channel_map = pd.read_excel(outdir+'CPTAC_TMT_reference_channels.xlsx')
    ref_channel = ref_channel_map.loc[ref_channel_map['Dataset']==ds, 'MQ sample'].values[0]
samples = samples_list[datasets.index(ds)]

MQ_protein_dir = aa_subs_dir +'MQ_validation_output/' # directory where validation search results (proteinGroups.txt) are stored
val_evidence_dict = pickle.load(open(aa_subs_dir+'Validation_search_evidence_dict.p', 'rb'))



# distribute precursor intensity by ratios of reporter ions in tmt set. these values used for sample-level quantification of peptide
def distribute_prec_int(protein_row, pg):

    # here we do not use normalized values, because we want the differences in sample loading to be accounted for in distribution of precursor intensity
    row_reporter_int = protein_row[[x for x in pg.columns if 'Reporter intensity corrected' in x]].values
    row_reporter_int_sum = np.sum([x for x in row_reporter_int if ~np.isnan(x)])
    
    prec_int = protein_row['Intensity']
    
    distributed_int = [prec_int*(x/row_reporter_int_sum) if ~np.isnan(x) else np.nan for x in row_reporter_int]
    
    return(distributed_int)

# normalize for sample loading
def median_normalize(sample_raw):
    """
    Input: list of raw intensities from a single TMT channel
    Output: median-normalized list of intensities for TMT channel
    """
    sample_median = np.nanmedian(sample_raw)
    sample_norm = [x/sample_median for x in sample_raw]
    return(sample_norm)
    
# normalize to reference sample (TMT data only)
def reference_normalize(med_norm_data, reference_int):
    ref_norm_data = [x/reference_int for x in med_norm_data]
    return(ref_norm_data)


allprot_abundance_dict = {}
for s in samples: 
    print(s)
    pg = pd.read_csv(MQ_protein_dir+s+'/txt/proteinGroups.txt', '\t', low_memory=False)
    
    prot_ids = pg['Protein IDs'].values
    abund_vals = pg.loc[:, [x for x in pg.columns if 'Reporter intensity corrected' in x]]
    
        
    s_rows = []
    if ds!= 'Healthy':
        s_reporter_int = pg.loc[:,[x for x in pg.columns if 'Reporter intensity corrected' in x]]
        s_reporter_median_norm = s_reporter_int.apply(median_normalize, axis=0)
        s_reporter_ref_norm = s_reporter_median_norm
        for i,row in s_reporter_ref_norm.iterrows():
            ref_int = row['Reporter intensity corrected '+str(ref_channel)+' 1']
            norm_vals = reference_normalize(row.values, ref_int)
            s_reporter_ref_norm.loc[i,:] = norm_vals
            
        for i,row in pg.iterrows():
            precursor = row['Intensity']
            distributed_int = distribute_prec_int(row, pg)
            s_rows.append([precursor]+distributed_int)
        s_df = pd.DataFrame(s_rows, columns = ['Precursor intensity']+['PrecDist '+str(x) for x in list(range(1,len(abund_vals.columns)+1))], index=prot_ids)
        for col in s_df.columns[1:]:
            n = col.split(' ')[-1]
            s_df['Median normalized '+n] = s_reporter_median_norm['Reporter intensity corrected '+n + ' 1']
            s_df['Reference normalized '+n] = s_reporter_ref_norm['Reporter intensity corrected '+n + ' 1']
    
    else:
        for i,row in pg.iterrows():
            precursor = row['Intensity']
            s_rows.append([precursor])
        s_df = pd.DataFrame(s_rows, columns = ['Precursor intensity'], index=prot_ids)
    
        for col in s_df.columns[1:]:
            median_norm_col = median_normalize(s_df[col].values)
            s_df['Median normalized '+col] = median_norm_col

    allprot_abundance_dict[s] = s_df


pickle.dump(allprot_abundance_dict, open(aa_subs_dir+'Allprot_normalized_abundance_dict.p', 'wb'))


# create a spreadsheet from data in Allprot_normalized_abundance_dict.p that is proteins x patients containing reporter-ion level abundances
# for TMT data 
allprots = []
for tmt in tmt_sets:
    tmt_df = prot_abund_dict[tmt]
    tmt_prots = [y for z in [x.split(';') for x in tmt_df.index] for y in z]
    allprots = allprots + tmt_prots
all_prots = list(set(allprots))

tmt_samples = {tmt:list(sample_map.loc[sample_map['sample_ID']==tmt,'sample_name'].values) for tmt in tmt_sets}
sample_mq = {sample:sample_map.loc[sample_map['sample_name']==sample,'MQ'].values[0] for sample in samples}
mainprot_patient_quant_df = pd.DataFrame(index=all_prots, columns=samples)

for p,prot in enumerate(all_prots):
    if p%1000==0:
        print(p)
    for tmt in tmt_sets:
        tmt_df = prot_abund_dict[tmt]
        prot_df = tmt_df.loc[[i for i in tmt_df.index if prot in i.split(';')],:]
        
        tmt_patients = tmt_samples[tmt]
        for patient in tmt_patients:
            patient_mq = sample_mq[patient]
            
            patient_int = np.nansum(prot_df['PrecDist '+str(patient_mq)])
            mainprot_patient_quant_df.loc[prot,patient] = patient_int
            
mainprot_patient_quant_df.to_excel(aa_subs_dir+'Main_protein_patient_reporter_quant_df.xlsx')
