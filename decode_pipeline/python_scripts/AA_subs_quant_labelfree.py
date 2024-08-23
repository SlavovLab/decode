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


""" This script is for use with label-free datasets. It follows the AAS_detection and AAS_validation scripts. It takes in highly confidently identified peptides with AAS, and quantifies them. Output metrics include normalized abundances, and ratios of peptide with AAS to its base peptide, at precursor ion levels, stored in a dict that structured by unique MTP-BP pairs.

    This script also generates a dictionary of quant data organized by AAS type
"""
    
# setting directories and reading in data
code_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/AA_subst_pipeline/'
proj_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Papers/Wang_HealthyTissues/'
aa_subs_dir = proj_dir+'AA_subs_pipeline/'

mtp_dict = pickle.load(open(aa_subs_dir+'Ion_validated_MTP_dict.p', 'rb'))
samples = list(mtp_dict.keys())
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
            unq_pair_dicts.append({'MTP':mtp, 'BP':bp, 'AAS':sub, 'tissue_list':[s], 'ev_idx_list':ev_idx})
            unq_pair_dicts[len(unq_pair_dicts)-1]['tissue_evidence'] = {s:[ev_idx]}
        else:
            curr_dict_idx = [i for i,x in enumerate(unq_pair_dicts) if (x['MTP']==mtp) and (x['BP']==bp)][0]
            unq_pair_dicts[curr_dict_idx]['tissue_list'].append(s)
            [unq_pair_dicts[curr_dict_idx]['ev_idx_list'].append(x) for x in ev_idx]
            if s in unq_pair_dicts[curr_dict_idx]['tissue_evidence'].keys():
                [unq_pair_dicts[curr_dict_idx]['tissue_evidence'][s].append(x) for x in ev_idx]
            else:
                unq_pair_dicts[curr_dict_idx]['tissue_evidence'][s] = ev_idx
print(len(unq_pairs))

# initialize dictionary to store quant info
MTP_quant_dict = {}
for i,pair in enumerate(unq_pair_dicts):
    curr_dict = {'MTP_seq':pair['MTP'], 'BP_seq':pair['BP'], 'sub_index':[i for i,x in enumerate(pair['MTP']) if pair['MTP'][i]!=pair['BP'][i]], 'aa_sub':pair['AAS'],'tissues':list(set(pair['tissue_list'])), 'tissue_evidence':pair['tissue_evidence'],
                 'MTP_PrecInt':{x:np.nan for x in samples}, 'BP_PrecInt':{x:np.nan for x in samples},'Prec_ratio':{x:np.nan for x in samples}, 'Norm_MTP_PrecInt':{x:np.nan for x in samples}, 'Norm_BP_PrecInt':{x:np.nan for x in samples}}#,'Patient_dict':{}}
  #  curr_dict['Patient_dict'] = {x:{'MTP_ReportInt':np.nan,'BP_ReportInt':np.nan, 'Reporter_ratio':np.nan} for x in list(set(sample_map['sample_name']))}
    MTP_quant_dict[i] = curr_dict
#pickle.dump(MTP_quant_dict, open(aa_subs_dir+'MTP_quant_dict.p', 'wb'))

""" Functions for quantification"""
# normalize reporter ion intensities
def median_normalize(sample_raw):
    """
    Input: list of raw precursor intensities for tissue
    Output: median-normalized list of precursor intensities for tissue
    """
    sample_median = np.median(sample_raw)
    sample_norm = [x/sample_median for x in sample_raw]
    return(sample_norm)

# Quantify precursor intensities and ratio of MTP/BP precursor intensities
# precursor intensities encompass all samples in tmt set
def precursor_intensity_quant(k, tissue):
    """
    Input: k=index of MTP in MTP_quant_dict, tissue
    Output: precursor intensity of mtp, its bp, their ratio
            Precursor intensities are log-transformed values of those reported in evidence.txt result file from validation MQ search.
            Precursor intensities represent the sum of all precursors mapped to the peptide sequence (across multiple fractions, charge states)
    """
    bp = MTP_quant_dict[k]['BP_seq']
    mtp = MTP_quant_dict[k]['MTP_seq']
    ev_df = val_evidence_dict[tissue]
    bp_ev_df = ev_df.loc[ev_df['Sequence']==bp,:]
    bp_prec_int = np.sum([x for x in bp_ev_df['Intensity'].values if ~np.isnan(x)])
    norm_bp_prec_int = np.sum([x for x in bp_ev_df['Intensity'].values if ~np.isnan(x)])
    
    mtp_ev_df = ev_df.loc[ev_df['Sequence']==mtp,:]
    mtp_prec_int = np.sum([x for x in mtp_ev_df['Intensity'].values if ~np.isnan(x)])
    norm_mtp_prec_int = np.sum([x for x in mtp_ev_df['Intensity'].values if ~np.isnan(x)])
  
    prec_ratio = np.log2(mtp_prec_int/bp_prec_int)
    return([mtp_prec_int, bp_prec_int, norm_mtp_prec_int, norm_bp_prec_int, prec_ratio])


"""Quantification of precursor and reporter ions"""
print('Quantifying peptides with AAS')
for s in samples:
    print(s)
    tissue=s
    ev_df = val_evidence_dict[s]
    #for col in [x for x in ev_df if 'Reporter intensity corrected' in x]:
    for col in [x for x in ev_df if x=='Intensity']:
        ev_df[col+'_Normalized'] = median_normalize(ev_df[col].values)
    
    for k,v in MTP_quant_dict.items():
        if s in v['tissues']:
            mtp = v['MTP_seq']
            bp = v['BP_seq']
            bp_ev = ev_df.loc[ev_df['Sequence']==bp,:]
            mtp_ev = ev_df.loc[ev_df['Sequence']==mtp,:]

            mtp_prec_int, bp_prec_int, norm_mtp_prec_int, norm_bp_prec_int, prec_ratio = precursor_intensity_quant(k, tissue)
            v['BP_PrecInt'][tissue] = bp_prec_int
            v['MTP_PrecInt'][tissue] = mtp_prec_int
            v['Norm_BP_PrecInt'][tissue] = norm_bp_prec_int
            v['Norm_MTP_PrecInt'][tissue] = norm_mtp_prec_int
            v['Prec_ratio'][tissue] = prec_ratio


print('saving results to file')
pickle.dump(MTP_quant_dict, open(aa_subs_dir+'MTP_quant_dict.p', 'wb'))
