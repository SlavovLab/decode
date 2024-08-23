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


"""
    This code reads in dictionaries generated in AA_subs_validation1.py.
    Also requires dataset metadata, and evidence.txt result files from validation search.
    Output is MTP dict filtered for MTPs quantified in validation search and for peptides with b/y ion evidence covering site of AAS.
"""

# setting directories and reading in data
code_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/AA_subst_pipeline/'
proj_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Papers/Clark_ccRCC/'
aa_subs_dir = proj_dir+'AA_subs_pipeline/'# aa_subs_dir is the output folder
MQ_dir = aa_subs_dir+'MQ_validation_output//'
sample_map = pd.read_excel(aa_subs_dir+'sample_map.xlsx')
samples = ['S'+str(i) for i in list(set(sample_map['TMT plex']))]

mtp_dict = pickle.load(open(aa_subs_dir+'qMTP_dict.p', 'rb'))
samples = list(mtp_dict.keys())

""" Identify MTPs that are found in the validation MQ search"""

print('Building validation evidence dict')
# initialize and save dictionary of validation evidence data

""" run the chunk below once, then read in data file if running script again"""
val_evidence_dict = {} # create a dict with all DP evidence.txt data from all tmt sets, so do not need to read in multiple files later
for s in samples:
    print(s)
    evidence = pd.read_csv(MQ_dir+s+'/txt/evidence.txt', '\t', engine='python')
    #filter evidence files for PEP < 0.01 and PIF > 0.8
    evidence = evidence.loc[(evidence['PEP']<=0.01) & (evidence['PIF']>=0.8),:]
    val_evidence_dict[s] = evidence
pickle.dump(val_evidence_dict, open(aa_subs_dir+'Validation_search_evidence_dict.p', 'wb'))
#val_evidence_dict = pickle.load(open(aa_subs_dir+'Validation_search_evidence_dict.p', 'rb'))

# function for determining if sequence is identified in validation search
def seq_in_val_search(idx, tmt_set):
    """
    Input: index of mtp_dict[tmt_set], tmt_set (sample)
    Output: if peptide is found, output = [index in mtp_dict, index in evidence file, index in mtp list at idx]
    """
    val_ev_df = val_evidence_dict[tmt_set]
    ev_seqs = list(val_ev_df['Sequence'].values)
    
    mtp_list = mtp_dict[tmt_set]['mistranslated sequence'][idx]
    # some mtp entries have >1 putative AAS. If > identified in validation search, return as separate results
    for i, mtp in enumerate(mtp_list):
        if mtp in ev_seqs:
            return([idx, [i for i,x in enumerate(ev_seqs) if x==mtp], i])
        else:
            return None
        
# apply validation search function to each mtp entry.
print('Identifying mtps that are found in regular search')
val_hit_lists = {s:[] for s in samples}
for s, s_dict in mtp_dict.items():
   # print(s)
    for idx in s_dict['Raw file'].keys():
        result = seq_in_val_search(idx, s)
        if result:
            val_hit_lists[s].append(result)

# Loop through lists of results, create new dict of validated MTPs with link to evidence file index
val_mtp_dict = {}
for s,val_list in val_hit_lists.items():
    print(s)
    val_mtp_dict[s] = {k:{} for k in mtp_dict[s].keys()}
    val_mtp_dict[s]['idx_val_evidence'] = {}

    for i, val in enumerate(val_list):
        mtp_idx = val[0]
        seq_idx = val[2]
        ev_idx = val[1]

        for k in mtp_dict[s].keys():
            if (isinstance(mtp_dict[s][k][mtp_idx], list)) and len(mtp_dict[s][k][mtp_idx])>0: # this is to make sure that we are extracting the correct AAS data and q-values for the mtp found out of list of mtps
                val_mtp_dict[s][k][i] = mtp_dict[s][k][mtp_idx][seq_idx]
            else:
                val_mtp_dict[s][k][i] = mtp_dict[s][k][mtp_idx]
        val_mtp_dict[s]['idx_val_evidence'][i] = ev_idx
    print(len(val_mtp_dict[s]['idx_val_evidence']))
pickle.dump(val_mtp_dict, open(aa_subs_dir+'Validated_MTP_dict.p', 'wb'))
val_mtp_dict = pickle.load(open(aa_subs_dir+'Validated_MTP_dict.p', 'rb'))
        
"""Determine which validated sequences have b/y ion evidence covering site of mistranslation"""
print('Determining b/y ion evidence for peptides with AAS')

# function to determine the number of fragment ions supporting site of AAS
def n_frags_over_MTP(frag_match, mtp, sub_idx, tmt_set):
    """
    Input: fragment matches for peptide from msms.txt (MQ output file), peptide sequence, index of AAS on sequence, tmt_set
    Output: number of fragment ions covering site of AAS
    """
    count = 0
    for f, frag in enumerate(frag_match):
        mtp_frag=0
        if ('NH3' not in frag) and ('H2O' not in frag) and ('(' not in frag) and ('a' not in frag):
            if 'b' in frag:
                frag_start = 0
                frag_end = int(frag[1:])
                frag_seq = mtp[frag_start:frag_end]
                if frag_end>sub_idx:
                    mtp_frag = 1
            elif 'y' in frag:
                frag_start = -int(frag[1:])
                frag_seq = mtp[frag_start:]
                if len(mtp)+frag_start <= sub_idx:
                    mtp_frag=1
            count+=mtp_frag
        
    return(count)

# apply function and annotate val_mtp_dict
"""requires msms.txt files (MQ output)"""
for s in samples:
    print(s)
    ev = val_evidence_dict[s]
    #this takes a while to run. Could speed up by creating dict of msms.txt dataframes and read in.
    msms = pd.read_csv(MQ_dir+'/'+s+'/txt/msms.txt', '\t', low_memory=False)
    
    val_mtp_dict[s]['fragment_evidence'] = {}
    for k,v in val_mtp_dict[s]['aa subs'].items():
        seq = val_mtp_dict[s]['mistranslated sequence'][k]
        bp = val_mtp_dict[s]['DP Base Sequence'][k]
        sub_idx = [i for i,x in enumerate(bp) if seq[i]!=x][0]
        
        ev_idx = val_mtp_dict[s]['idx_val_evidence'][k]
        val_mtp_dict[s]['fragment_evidence'][k] = 0
        for idx in ev_idx:
            row = ev.iloc[idx,:]
            raw_file = row['Raw file']
            scan = row['MS/MS scan number']
            scan_row = msms.loc[(msms['Raw file']==raw_file) & (msms['Scan number']==scan),:]
            if len(scan_row)>0:
                frag_match = scan_row['Matches'].values[0].split(';')
                count_frags = n_frags_over_MTP(frag_match, seq, sub_idx, s)
                if count_frags>val_mtp_dict[s]['fragment_evidence'][k]:
                    val_mtp_dict[s]['fragment_evidence'][k] = count_frags

# filter val_mtp_dict for those with b/y ion evidence
val_ion_mtp_dict = {s:{} for s in samples}
for s in samples:
    ion_idx = [i for i,x in val_mtp_dict[s]['fragment_evidence'].items() if x>1]
    for k,v in val_mtp_dict[s].items():
        val_ion_mtp_dict[s][k] = {i:x for i,x in v.items() if i in ion_idx}
pickle.dump(val_ion_mtp_dict, open(aa_subs_dir+'Ion_validated_MTP_dict.p', 'wb'))
  

