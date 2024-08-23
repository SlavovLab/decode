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
import sys


"""
    This code reads in dictionaries generated in AA_subs_detection_HPC_MP_TMT.py or AA_subs_detection_HPC_MP_labelfree.py.
    Also requires dataset metadata, and evidence.txt result files from DP search
    Output is MTP dict filtered for Bayesian analysis q-values < determined threshold; fasta databases with MTPs appended.
"""

# setting directories and reading in data
proj_dir = '/scratch/tsour.s/Dou_UCEC/'
aa_subs_dir = proj_dir+'AA_subs_pipeline/'# aa_subs_dir is the output folder
MQ_dir = proj_dir+'MQ_output/'
database_dir = proj_dir+'databases/'

mtp_dict = pickle.load(open(aa_subs_dir+'MTP_dict.p', 'rb'))
samples = list(mtp_dict.keys())
fasta_str_list = [database_dir+s+'_noMTP.fasta' for s in samples] # a list of strings to each sample's custom protein database directory
uniprot_fasta_str = database_dir+'Human_uniprot.fasta'

"""Calculate q-values using Bayesian posterior probability and filter MTPs by empirically determined confidence threshold"""
print('calculating q-values with Bayesian analysis for all MTPs')

## function for calculating Bayesian posterior probability from positional probability (prior) and mass error (likelihood)
def posterior_aasub_prob(v, idx): # v is the sample dictionary, e.g. MTP_dict['S1']
    """
    Input: index of mtp in mtp_dict[S]
    Output: posterior probability of observation
    """
    pp = v['aa subs positional probability'][idx]
    merr = v['aa subs mass error (ppm)'][idx]
    PEP = v['DP PEP'][idx]
    Pi = []
    for i in range(len(pp)):
        # multiply the probability of getting the observed mass error by the prior positional probability
        Pi.append(sp.stats.norm.pdf(merr[i], mean, std)*pp[i])
    #normalize each posterior prob. by the sum of posterior probabilities observed for the peptide and multiply by (1-PEP)
    Pi = [(x/sum(Pi)*(1-PEP)) for x in Pi]
    return(Pi)


### generate mass error (likelihood) distribution
""" Requries evidence.txt files from DP search"""
mass_err_ppm = []
DP_evidence_dict = {} # create a dict with all DP evidence.txt data from all tmt sets, so do not need to read in multiple files later#
#DP_evidence_dict = pickle.load(open(aa_subs_dir+'DP_search_evidence_dict.p', 'rb'))
for s in samples:
    print(s)
    evidence = pd.read_csv(MQ_dir+s+'/txt/evidence.txt', '\t', engine='c')
    DP_evidence_dict[s] = evidence
    mass_err_ppm.append(list(evidence['Mass error [ppm]'].values))
    pickle.dump(DP_evidence_dict, open(aa_subs_dir+'DP_search_evidence_dict.p', 'wb'))
mass_err_ppm = [x for y in mass_err_ppm for x in y]
pickle.dump(DP_evidence_dict, open(aa_subs_dir+'DP_search_evidence_dict.p', 'wb'))

# mass error distribution stats
me = [x for x in mass_err_ppm if x >-1000]
mean = np.mean(me)
std = np.std(me)
pdf = sp.stats.norm.pdf(2, mean, std)

# apply posterior probability functional to all MTPs
post_prob = []
for k, v in mtp_dict.items():
    v['Posterior subs probability'] = {i: posterior_aasub_prob(v, i) for i in v['aa subs'].keys()}
    post = v['Posterior subs probability']
    for x in post.values():
        post_prob = post_prob + [i for i in x]

## calculate q values
pval = [1-p for p in post_prob]
ranked_pval = np.sort(pval)
cumsum = np.cumsum(ranked_pval)
ranked_qval = [x/(i+1) for i,x in enumerate(cumsum)]
rank_ind = np.argsort(pval)
qval = [ranked_qval[list(rank_ind).index(i)] for i in range(len(rank_ind))]

# add q-values to mtp_dict
count_q=0
for k,v in mtp_dict.items():
    post = v['Posterior subs probability']
    q_dict = {}
    for i,x in post.items():
        q_dict[i] = list(qval[count_q:count_q+len(x)])
        count_q+=len(x)
    v['q-value'] = q_dict
    
pickle.dump(mtp_dict, open(aa_subs_dir+'qMTP_dict.p', 'wb'))

# Precision-Recall moving threshold analysis to compute optimal confidence threshold
print('Computing PR-curve and optimical confidence threshold')
thresh = max(qval)
n_thresh = len([x for x in qval if x<=thresh])
TP_thresh = np.floor((1-thresh)*n_thresh)

metric_rows = []
for qt in np.logspace(-20,-1, num=100, base=10):
    TP = len([x for x in qval if x<=qt])*(1-qt) # N true pos
    FP = len([x for x in qval if x<=qt])*qt     # N false pos
    FN = TP_thresh - TP                         # N false neg
    TN = len(qval) - TP - FN - FP               # N true neg
    P = TP/(TP+FP)                              # Precision
    R = TP/(TP+FN)                              # Recall
    F_score = (2*P*R)/(P+R)                     # F-score (predictive performance)
    metric_rows.append([qt, TP, FP, FN, TN, P, R, F_score])
metric_df = pd.DataFrame(metric_rows, columns=['q_threshold', 'TP', 'FP', 'FN', 'TN', 'Precision', 'Recall', 'F_score'])
metric_df.to_excel(aa_subs_dir+'tonsil_q-value_Precision_Recall_data.xlsx')

# optimal threshold is at highest F_score
max_F_idx = list(metric_df['F_score'].values).index(max(metric_df['F_score'].values))
q_thresh = metric_df.loc[max_F_idx, 'q_threshold']

# filter mtp_dict for peptides with q-values below determined threshold
filtered_mtp_dict = {}
for s in samples:
    filtered_mtp_dict[s] = {}
    dp_dict = mtp_dict[s]
    hc_idx = [i for i,x in dp_dict['q-value'].items() if any(y<=q_thresh for y in x)]
    for k,v in dp_dict.items():
        filtered_mtp_dict[s][k] = {i:x for i,x in v.items() if i in hc_idx}

pickle.dump(filtered_mtp_dict, open(aa_subs_dir+'qMTP_dict.p', 'wb'))

# filter filtered_mtp_dict for any peptides that have a match in the custom or uniprot fasta file
def in_human_fasta(seq, fasta_seqs):
    if any(re.search(seq, x) for x in fasta_seqs):
        return(True)
    else:
        return(False)
uniprot_fasta = open(uniprot_fasta_str, 'r').read()
uniprot_fasta_entries = uniprot_fasta.split('>')
uniprot_fasta_seqs = [''.join(x.split('\n')[1:-1]) for x in uniprot_fasta_entries]


keys2trim = ['origin aa', 'origin aa index', 'aa subs', 'aa subs positional probability', 'aa subs mass error (ppm)', 'destination aa', 'mistranslated sequence', '1-frame genome substring', '2-frame genome substring', '3-frame genome substring', '4-frame genome substring','5-frame genome substring', '6-frame genome substring', 'Posterior subs probability', 'q-value']

seq_filtered_mtp_dict = {}
for i,s in enumerate(samples):
    fasta_str = fasta_str_list[i]
    fasta = open(fasta_str, 'r').read()
    fasta_entries = fasta.split('>')
    fasta_seqs = [''.join(x.split('\n')[1:-1]) for x in fasta_entries]
    candidate_saap = [x for y in list(filtered_mtp_dict[s]['mistranslated sequence'].values()) for x in y]
    cand_saap = [x for x in candidate_saap if not in_human_fasta(x, fasta_seqs) and not in_human_fasta(x, uniprot_fasta_seqs)]
    
    s_filtered_dict = {k:{} for k in filtered_mtp_dict[s].keys()}
    for idx,seqs in filtered_mtp_dict[s]['mistranslated sequence'].items():
        keep_seq_idx = [i for i,x in enumerate(seqs) if x in cand_saap]
        if len(keep_seq_idx)>0:
            for k,v in filtered_mtp_dict[s].items():
                v = v[idx]
                if k in keys2trim:
                    v = [v[i] for i in keep_seq_idx]
                s_filtered_dict[k][idx] = v
    seq_filtered_mtp_dict[s] = s_filtered_dict
            
pickle.dump(seq_filtered_mtp_dict, open(aa_subs_dir+'qMTP_dict_'+celltype+'.p', 'wb'))

mtp_dict = filtered_mtp_dict

"""Add high-confidence MTPs to TMT set-specific protein fasta databases"""
print('Adding MTPs to fasta databases')
noMTP_fasta = glob(database_dir+'*noMTP.fasta')
for s in samples:
    print(s)
    ff = [x for x in noMTP_fasta if s+'_' in x][0]
    fasta_entries = open(ff, 'r').read().split('>')
    
    sdict = mtp_dict[s]
    base_scans = sdict['DP Base Scan Number']
    dp_scans = sdict['DP Mod Scan Number']
    seqs = sdict['mistranslated sequence']
    
    unq_seqs = list(set([x for y in list(seqs.values()) for x in y]))
    unq_bp_scans = [';'.join([str(base_scans[i]) for i in base_scans.keys() if seq in seqs[i]]) for seq in unq_seqs]
    unq_dp_scans = [';'.join([str(dp_scans[i]) for i in dp_scans.keys() if seq in seqs[i]]) for seq in unq_seqs]
    
    fasta = database_dir + s+'_MTP.fasta'
    
    with open(fasta, 'w') as f:
        for entry in fasta_entries:
            f.write('>'+entry)
        for i,seq in enumerate(unq_seqs):
            seq_idx = ';'.join([str(j) for j,s in seqs.items() if seq in s])
            bp_scan = unq_bp_scans[i]
            dp_scan = unq_dp_scans[i]
            f.write('>MTP|'+seq_idx+'_base'+bp_scan+'_DP'+dp_scan+'\n'+seq+'\n')
    f.close()
""" RE-SEARCH THE RAW PROTEOMICS DATA WITH REGULAR SEARCH IN MAXQUANT. USE FASTA WITH MTPs APPENDED TO END AS DATABASE. """
# transfer appended fasta files to HPC
# use gen_mqpar.py to create sbatch files and xml files for MQ searches on cluster
# search raw proteomics data with fastas that have MTPs
