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
from copy import deepcopy


""" This script is for use with TMT-labeled datasets. It follows the AAS_detection and AAS_validation scripts. It takes in highly confidently identified peptides with AAS, and quantifies them. Output metrics include normalized abundances, and ratios of peptide with AAS to its base peptide, at both precursor and reporter-ion levels, stored in a dict that structured by unique MTP-BP pairs, as opposed to TMT sets.

    This script also generates a dictionary of quant data organized by AAS type
    
    MTP = mistranslated peptide (SAAP or substituted amino acid peptide)
    BP = base peptide (RNA-tempalted peptide)
"""
    
# setting directories and reading in data
code_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/AA_subst_pipeline/'
proj_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Papers/Cao_PDAC/'
aa_subs_dir = proj_dir+'AA_subs_pipeline/'
sample_map = pd.read_excel(aa_subs_dir+'sample_map.xlsx')
samples = ['S'+str(i) for i in list(set(sample_map['TMT plex']))]
MQ_TMT_dict = {'126':1,'127N':2, '127C':3, '128N':4, '128C':5, '129N':6, '129C':7, '130N':8, '130C':9, '131':10, '131C':11}

mtp_dict = pickle.load(open(aa_subs_dir+'Ion_validated_MTP_dict.p', 'rb'))
val_evidence_dict = pickle.load(open(aa_subs_dir+'Validation_search_evidence_dict.p', 'rb'))

print("initializing")
# get unique pairs of MTP-BP seqs across all TMT sets
unq_pairs = []
unq_pair_dicts = []
for s, s_dict in mtp_dict.items():
    print(s)
    for i,v in s_dict['aa subs'].items():
        sub = v
        mtp = s_dict['mistranslated sequence'][i]
        bp = s_dict['DP Base Sequence'][i]
        ev_idx = s_dict['idx_val_evidence'][i]
        pair = [mtp, bp]
        if (pair not in unq_pairs):
            unq_pairs.append(pair)
            unq_pair_dicts.append({'MTP':mtp, 'BP':bp, 'AAS':sub, 'TMT_list':[s], 'ev_idx_list':ev_idx})
            unq_pair_dicts[len(unq_pair_dicts)-1]['TMT_evidence'] = {s:[ev_idx]}
        else:
            curr_dict_idx = [i for i,x in enumerate(unq_pair_dicts) if (x['MTP']==mtp) and (x['BP']==bp)][0]
            unq_pair_dicts[curr_dict_idx]['TMT_list'].append(s)
            [unq_pair_dicts[curr_dict_idx]['ev_idx_list'].append(x) for x in ev_idx]
            if s in unq_pair_dicts[curr_dict_idx]['TMT_evidence'].keys():
                [unq_pair_dicts[curr_dict_idx]['TMT_evidence'][s].append(x) for x in ev_idx]
            else:
                unq_pair_dicts[curr_dict_idx]['TMT_evidence'][s] = ev_idx
print(len(unq_pairs))

# initialize dictionary to store quant info
MTP_quant_dict = {}
for i,pair in enumerate(unq_pair_dicts):
    curr_dict = {'MTP_seq':pair['MTP'], 'BP_seq':pair['BP'], 'sub_index':[i for i,x in enumerate(pair['MTP']) if pair['MTP'][i]!=pair['BP'][i]], 'aa_sub':pair['AAS'],'tmt_sets':list(set(pair['TMT_list'])), 'tmt_evidence':pair['TMT_evidence'],
                 'MTP_PrecInt':{x:np.nan for x in samples}, 'BP_PrecInt':{x:np.nan for x in samples},'Prec_ratio':{x:np.nan for x in samples},'Patient_dict':{}}
    curr_dict['Patient_dict'] = {x:{'MTP_ReportInt':np.nan,'BP_ReportInt':np.nan, 'Reporter_ratio':np.nan} for x in list(set(sample_map['sample_name']))}
    MTP_quant_dict[i] = curr_dict
#pickle.dump(MTP_quant_dict, open(aa_subs_dir+'MTP_quant_dict.p', 'wb'))

""" Functions for quantification"""
# normalize reporter ion intensities
def reporter_ion_normalize(sample_raw):
    """
    Input: list of raw intensities from a single TMT channel
    Output: median-normalized list of intensities for TMT channel
    """
    sample_median = np.median(sample_raw)
    sample_norm = [x/sample_median for x in sample_raw]
    return(sample_norm)

# Quantify precursor intensities and ratio of MTP/BP precursor intensities
# precursor intensities encompass all samples in tmt set
def precursor_intensity_quant(k, tmt_set):
    """
    Input: k=index of MTP in MTP_quant_dict, tmt_set
    Output: precursor intensity of mtp, its bp, their ratio
            Precursor intensities are log-transformed values of those reported in evidence.txt result file from validation MQ search.
            Precursor intensities represent the sum of all precursors mapped to the peptide sequence (across multiple fractions, charge states)
    """
    bp = MTP_quant_dict[k]['BP_seq']
    mtp = MTP_quant_dict[k]['MTP_seq']
    ev_df = val_evidence_dict[tmt_set]
    bp_ev_df = ev_df.loc[ev_df['Sequence']==bp,:]
    bp_prec_int = bp_ev_df['Intensity'].sum()
    
    mtp_ev_df = ev_df.loc[ev_df['Sequence']==mtp,:]
    mtp_prec_int = mtp_ev_df['Intensity'].sum()
  
    prec_ratio = np.log10(mtp_prec_int/bp_prec_int)
    return([mtp_prec_int, bp_prec_int, prec_ratio])


# distribute precursor intensity by ratios of reporter ions in tmt set. these values used for sample-level quantification of peptide
def distribute_prec_int(k, tmt_set):
    """
    Input: k=index of MTP in MTP_quant_dict, tmt set
    Outpt: list of reporter intensity values for TMT channels that reflect precursor intensity distributed by fraction of each reporter ion to total reporter signal
    """
    ev_df = val_evidence_dict[tmt_set]
    bp = MTP_quant_dict[k]['BP_seq']
    bp_ev = ev_df.loc[ev_df['Sequence']==bp,:]
    bp_prec_int = MTP_quant_dict[k]['BP_PrecInt'][tmt_set]
    mtp = MTP_quant_dict[k]['MTP_seq']
    mtp_ev = ev_df.loc[ev_df['Sequence']==mtp,:]
    mtp_prec_int = MTP_quant_dict[k]['MTP_PrecInt'][tmt_set]

    # here we do not use normalized values, because we want the differences in sample loading to be accounted for in distribution of precursor intensity
    bp_reporter_int = bp_ev.loc[:,[x for x in bp_ev.columns if ('Reporter intensity corrected' in x) and ('_Normalized' not in x)]]
    bp_reporter_int_sums = bp_reporter_int.sum(axis=0)
    bp_total_reporter_int = bp_reporter_int_sums.sum()
    bp_distributed = bp_prec_int * bp_reporter_int_sums / bp_total_reporter_int

    mtp_reporter_int = mtp_ev.loc[:,[x for x in mtp_ev.columns if ('Reporter intensity corrected' in x) and ('_Normalized' not in x)]]
    mtp_reporter_int_sums = mtp_reporter_int.sum(axis=0)
    mtp_total_reporter_int = mtp_reporter_int_sums.sum()
    mtp_distributed = mtp_prec_int * mtp_reporter_int_sums / mtp_total_reporter_int


    return(bp_distributed, mtp_distributed)

"""Quantification of precursor and reporter ions"""
print('Quantifying peptides with AAS')
for s in samples:
    print(s)
    tmt_set=s
    ev_df = val_evidence_dict[s]
    for col in [x for x in ev_df if 'Reporter intensity corrected' in x]:
        val_evidence_dict[s][col+'_Normalized'] = reporter_ion_normalize(ev_df[col].values)

    for k,v in MTP_quant_dict.items():
        if s in v['tmt_sets']:
            mtp = v['MTP_seq']
            bp = v['BP_seq']
            bp_ev = ev_df.loc[ev_df['Sequence']==bp,:]
            mtp_ev = ev_df.loc[ev_df['Sequence']==mtp,:]

            mtp_prec_int, bp_prec_int, prec_ratio = precursor_intensity_quant(k, tmt_set)
            v['BP_PrecInt'][tmt_set] = bp_prec_int
            v['MTP_PrecInt'][tmt_set] = mtp_prec_int
            v['Prec_ratio'][tmt_set] = prec_ratio

            bp_distributed, mtp_distributed = distribute_prec_int(k, tmt_set)

            # Now we will assign the distributed reporter ion values to the respective patients
            s_patients = sample_map.loc[sample_map['sample_ID']==s, 'sample_name'].values

            p_dict = v['Patient_dict']
            for pat in s_patients:
                mq_reporter = sample_map.loc[sample_map['sample_name']==pat, 'MQ'].values[0]

                patient_bp_reporter = bp_distributed['Reporter intensity corrected '+str(mq_reporter)]
                patient_mtp_reporter = mtp_distributed['Reporter intensity corrected ' +str(mq_reporter)]
                patient_ratio = np.log10(patient_mtp_reporter/patient_bp_reporter)

                patient_bp_reporter_norm = bp_ev['Reporter intensity corrected '+str(mq_reporter)+'_Normalized'].sum()
                patient_mtp_reporter_norm = mtp_ev['Reporter intensity corrected '+str(mq_reporter)+'_Normalized'].sum()

                p_dict[pat]['BP_ReportInt'] = patient_bp_reporter
                p_dict[pat]['MTP_ReportInt'] = patient_mtp_reporter
                p_dict[pat]['BP_ReportInt_Norm'] = patient_bp_reporter_norm
                p_dict[pat]['MTP_ReportInt_Norm'] = patient_mtp_reporter_norm
                p_dict[pat]['Reporter_ratio'] = patient_ratio

print('saving results to file')
pickle.dump(MTP_quant_dict, open(aa_subs_dir+'MTP_quant_dict.p', 'wb'))
