#!/bin/usr python
import pandas as pd
import numpy as np
import pickle
from glob import glob
import re
import scipy as sp
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os
import copy
from collections import Counter
import random
import sys
import multiprocessing as mp

"""
This script will produce a dictionary and spreadsheets containing the precursor and reporter-ion level abundances of all main peptides in the dataset, quantified the same way as SAAP and BPs.

The most efficient way to run this was with multiprocessing on an HPC/cluster environment.

Input:
    1. dataset info (directory, samples, sample_map, output directory)
    2. Validation_search_evidence_dict.p (output of AA_subs_validation2.py)
Returns:
    1. main_peptide_quant_dict.p
    2. main_peptide_precursor_quant_df.xlsx
    3. main_peptide_reporter_quant_df.xlsx
"""

proj_dir = '/scratch/tsour.s/Sat_LSCC/'
aa_subs_dir = proj_dir + 'AA_subs_pipeline/'

val_ev_dict = pickle.load(open(aa_subs_dir+'Validation_search_evidence_dict.p', 'rb'))

# if TMT data:
sample_map = pd.read_excel(aa_subs_dir+'sample_map.xlsx', index_col=0)
samples = list(set(sample_map['sample_name']))
tmt_sets = list(set(sample_map['sample_ID']))

# if label free, just provide list of samples
#samples = open(aa_subs_dir+'proteomic_tissue_list.txt', 'r').read().split('\n')


# function to get precursor abundance of a peptide
def precursor_intensity_quant(pep_idx, s_ev):
    pep_ev_df = s_ev.iloc[pep_idx,:]
    prec_int = pep_ev_df['Intensity'].sum()
    return(prec_int)

# function to get reporter-ion level abundance of a peptide (TMT only)
def distribute_prec_int(pep_precursor, pep_idx, s_ev):
    pep_ev = s_ev.iloc[pep_idx,:]
    pep_reporter_int = pep_ev.loc[:,[x for x in pep_ev.columns if ('Reporter intensity corrected' in x) and ('_Normalized' not in x)]]
    pep_reporter_int_sums = pep_reporter_int.sum(axis=0)
    pep_total_reporter_int = pep_reporter_int_sums.sum()
    pep_distributed = pep_precursor * pep_reporter_int_sums / pep_total_reporter_int
    return(pep_distributed)

# get list of unique peptide sequences
set_peps = []
for s in samples:
    s_ev = val_ev_dict[s]
    peps = list(set(s_ev['Sequence'].values))
    set_peps = set_peps + peps
set_peps = list(set(set_peps))

# function to get the precursor and reporter ion level abundances of a peptide sequence, returned as a dictionary
def get_peptide_dict(pep):
    pep_dict = {pep:{}}
    for s in samples:
        print(s)
        s_ev = val_ev_dict[s]
        s_vals = s_ev.values # sequence is first column
        pep_idx = [i for i in range(len(s_vals)) if s_vals[i][0]==pep]
        pep_precursor = precursor_intensity_quant(pep_idx, s_ev)
        if 'Precursor' not in pep_dict.keys():
            pep_dict['Precursor'] = {s:pep_precursor}
        else:
            pep_dict['Precursor'][s] = pep_precursor
            
        # comment out reporter level quant for label-free data
        pep_distributed = distribute_prec_int(pep_precursor, pep_idx, s_ev)
        s_patients = sample_map.loc[sample_map['sample_ID']==s, 'sample_name'].values
        for pat in s_patients:
            mq_reporter = sample_map.loc[sample_map['sample_name']==pat, 'MQ'].values[0]
            patient_reporter_val = pep_distributed['Reporter intensity corrected '+str(mq_reporter)]
            if 'Patient_dict' not in pep_dict.keys():
               pep_dict['Patient_dict'] = {pat:patient_reporter_val}
            else:
                pep_dict['Patient_dict'][pat] = patient_reporter_val
    pickle.dump(pep_dict, open(aa_subs_dir+pep+'_quant_dict.p','wb'))
    return(pep_dict)

# split up job across all peptide sequences
p = mp.Pool()
pep_dict_list = p.map(get_peptide_dict, set_peps)
p.close()
p.join()

# combine the results of multiprocessing into a single dictionary
main_pep_quant_dict = {i:pep_dict_list[i] for i in range(len(pep_dict_list))}
pickle.dump(main_pep_quant_dict, open(aa_subs_dir+'main_peptide_quant_dict.p', 'wb'))

# create spreadsheets for ease of use in later applications

# precursor level quant spreadsheet
peps = [v['Peptide'] for k,v in main_pep_quant_dict.items()]
samples = list(mainpep_quant[0]['Precursor'].keys())
prec_mainpep_quant_df = pd.DataFrame(index=peps, columns=samples)
for k,v in main_pep_quant_dict.items():
    pep = v['Peptide']
    for t,val in v['Precursor'].items():
        prec_mainpep_quant_df.loc[pep,t]  = val
prec_mainpep_quant_df.to_excel(aa_subs_dir+'main_peptide_precursor_quant_df.xlsx')
     
# TMT only - reporter level quant spreadsheet
patients = sample_map['sample_name'].to_list()
mainpep_quant_df = pd.DataFrame(index=peps, columns=patients)
for k,v in main_pep_quant_dict.items():
    pep = v['Peptide']
    for p,val in v['Patient_dict'].items():
        mainpep_quant_df.loc[pep,p]  = val
mainpep_quant_df.to_excel(aa_subs_dir+'main_peptide_reporter_quant_df.xlsx')


