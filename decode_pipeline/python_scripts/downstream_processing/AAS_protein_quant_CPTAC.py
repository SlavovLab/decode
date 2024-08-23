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


"""
This script takes in results from earlier sections of pipelines to output 2 dictionaries of AAS at the protein level.
Input:
    1. Ion_validation_MTP_quant_dict.p (output of AA_subs_validation2.py)
    2. MTP_quant_dict.p (output of AA_subs_quant.py)
    3. DP_search_evidence_dict.p (output of AA_subs_validation1.py)
    4. sample-specific fasta files (generated in AA_subs_validation1.py)
    5. sample_map.xlsx or list of samples
    6. blast_dict.p (output of compile_blast_results.py)
Output:
    1. MTProtein_isoform_dict.p.
        Data dictionary with protein names (headers in fasta files) as keys, and substitution data as values.
    2. MTProtein_ref_dict.p
        Data dictionary with RefSeq protein names as keys and values of all entries in MTProtein_isoform_dict.p that mapped to the RefSeq protein via blast.
"""

# setting directories and reading in data
code_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/AA_subst_pipeline/'
proj_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Papers/Sat_LSCC/'
aa_subs_dir = proj_dir+'AA_subs_pipeline/'
fasta_dir = proj_dir+'databases/' # need sample-specific fasta files for mapping protein sequences to blast results

sample_map = pd.read_excel(aa_subs_dir+'sample_map.xlsx')
samples = ['S'+str(i) for i in list(set(sample_map['TMT plex']))]
MQ_TMT_dict = {'126':1,'127N':2, '127C':3, '128N':4, '128C':5, '129N':6, '129C':7, '130N':8, '130C':9, '131':10, '131C':11}


mtp_dict = pickle.load(open(aa_subs_dir+'Ion_validated_MTP_dict.p', 'rb'))
mtp_quant_dict = pickle.load(open(aa_subs_dir+'MTP_quant_dict.p', 'rb'))
dp_ev_dict = pickle.load(open(aa_subs_dir+'DP_search_evidence_dict.p', 'rb'))
blast_dict = pickle.load(open(aa_subs_dir+'blast_dict.p', 'rb'))

print('generating dictionary of proteins with AAS')
""" for each uniquely identified razor protein of peptides with AAS, populate dictionary with all MTP-BP data mapping to protein """

# initializing dictionary
mtprot_dict = {}
patients = list(set(sample_map['sample_name'].values))
protein_dict = {'Protein_sequence':'', 'MTP-BP_list':[], 'MTP-BP_pairs':[]}
mtp_bp_pair_dict = {'MTP':'', 'BP':'', 'aa_sub':'', 'aa_sub_idx':'', 'in_TMT_sets':[]}


""" Functions used below to populate mtprot_dict """
# function to get razor protein associated with BP for each MTP
def get_rzr_MTprot(s, mtp_idx):
    """
    Input: tmt_set, index of mtp in mtp_dict
    Output: razor protein of the bp
    """
    rf = mtp_dict[s]['Raw file'][mtp_idx]
    bp_scan = mtp_dict[s]['DP Base Scan Number'][mtp_idx]
    bp_seq = mtp_dict[s]['DP Base Sequence'][mtp_idx]
    try:
        bp_rzr_prot = dp_ev_dict[s].loc[(dp_ev_dict[s]['Sequence']==bp_seq) & (dp_ev_dict[s]['Raw file']==rf), 'Leading razor protein'].values[0]
    except: # if for some reason the base peptide is not in the evidence file, use the first protein in allPeptides 'DP Proteins' column
        bp_rzr_prot = mtp_dict[s]['DP Proteins'][mtp_idx]
    return(bp_rzr_prot)

# function to add quant info 
def annotate_quant(p_dict):
    """
    Input: entry in mtprot_dict
    Output: No output. Annotates mtprot_dict entry with MTP quantification data
    """
    pairs = p_dict['MTP-BP_pairs']
    for pair in pairs:
        mtp = pair['MTP']
        bp = pair['BP']
        quant_dict = [x for i,x in mtp_quant_dict.items() if (x['MTP_seq']==mtp) and (x['BP_seq']==bp)][0]
        p_dict['MTP_PrecInt'] = quant_dict['MTP_PrecInt']
        p_dict['BP_PrecInt'] = quant_dict['BP_PrecInt']
      #  p_dict['Norm_MTP_PrecInt'] = quant_dict['Norm_MTP_PrecInt']
      #  p_dict['Norm_BP_PrecInt'] = quant_dict['Norm_BP_PrecInt']
        p_dict['Prec_ratio'] = quant_dict['Prec_ratio']
        p_dict['Patient_dict'] = quant_dict['Patient_dict']

# function to get sequence of protein - need sequence to link to blast results
def prot_details(p, s_fasta):
    """
    Input: protein name, fasta file directory
    Output: protein sequence
    """
    seqidx = [i for i, line in enumerate(s_fasta) if (p in line) or (all(x in line for x in p.split('|'))) or (all(x in line for x in p.split('_')))]
    # first term in if statement is for reference type headers, second is for STRG headers, third is for customProDB headers
    if len(seqidx)>0:
        seqidx = seqidx[0]
        if re.match('[A-Z]{7,}', s_fasta[seqidx+1]):
            seqseq = s_fasta[seqidx+1]
        else:
            seqseq = s_fasta[seqidx+2]
    else:
        seqseq = ''
    return(seqseq)


# iterate through MTPs, populate mtprot_dict with MTP data
for s in samples:
    print(s)
    for mtp_idx in list(mtp_dict[s]['Raw file'].keys()):
        rzr_prot = get_rzr_MTprot(s, mtp_idx)
        rzr_prot = rzr_prot.strip('|')
        
        mtp = mtp_dict[s]['mistranslated sequence'][mtp_idx]
        bp = mtp_dict[s]['DP Base Sequence'][mtp_idx]
        aas = mtp_dict[s]['aa subs'][mtp_idx]
        aas_idx = [i for i,x in enumerate(bp) if mtp[i]!=x][0]
        
        if rzr_prot not in mtprot_dict.keys():
            mtprot_dict[rzr_prot] = copy.deepcopy(protein_dict)
            mtprot_dict[rzr_prot]['MTP-BP_list'].append([mtp,bp])
            mtprot_dict[rzr_prot]['MTP-BP_pairs'].append(copy.deepcopy(mtp_bp_pair_dict))
            mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['MTP'] = mtp
            mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['BP'] = bp
            mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['aa_sub'] = aas
            mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['aa_sub_idx'] = aas_idx

        else:
            if [mtp,bp] not in mtprot_dict[rzr_prot]['MTP-BP_list']:
                mtprot_dict[rzr_prot]['MTP-BP_list'].append([mtp,bp])
                mtprot_dict[rzr_prot]['MTP-BP_pairs'].append(copy.deepcopy(mtp_bp_pair_dict))
                mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['MTP'] = mtp
                mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['BP'] = bp
                mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['aa_sub'] = aas
                mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['aa_sub_idx'] = aas_idx

        mtp_quant_idx = [i for i,x in mtp_quant_dict.items() if (x['MTP_seq'] ==mtp) and x['BP_seq'] ==bp][0]
        tmt_sets = mtp_quant_dict[mtp_quant_idx]['tmt_sets']
        mtprot_dict[rzr_prot]['MTP-BP_pairs'][-1]['in_TMT_sets'] = tmt_sets
# add quant info
for p, p_dict in mtprot_dict.items():
    annotate_quant(p_dict)

# add protein sequence for mapping to blast
for s in samples:
    print(s)
    s_fasta = open(fasta_dir+s+'_MTP.fasta', 'r').readlines()
    s_fasta = [x for x in s_fasta if x!= '\n']
    for p,pdict in mtprot_dict.items():
        if pdict['Protein_sequence'] == '':
            if (s in pdict['MTP-BP_pairs'][0]['in_TMT_sets']) and (pdict['Protein_sequence']==''):
                seq = prot_details(p, s_fasta)
                pdict['Protein_sequence'] = seq.split('\n')[0]

ickle.dump(mtprot_dict, open(aa_subs_dir+'MTProtein_isoform_dict.p', 'wb'))


"""
For downstream analysis, we want to aggregate data from proteoforms to their reference proteins.
To do this, we blasted all protein sequences against RefSeq all proteins, and annotated each sequence with the top hit.
The results are compiled by compile_blast_dict.py which outputs blast_dict.p and used here to link all proteoforms of a reference protein together.
"""

print('Aggregating data to reference proteins')

# annotate each p_dict with blast result
def annotate_blast_results(prot):
    """
    Input: protein name
    Output: None. Annotates dictionary with blast results
    """
    pdict = mtprot_dict[prot]
    prot_seq = pdict['Protein_sequence'].split('\n')[0]
    if 'Blast' not in pdict.keys():
        pdict['Blast'] = {}
    for s in samples:
        sblast = blast_dict[s]
        seq_idx = [i for i,x in enumerate(sblast['query']) if x==prot_seq]
        if len(seq_idx)>0:
           # print('hit')
            prot_name = sblast['title'][seq_idx[0]]
            frac_pos = sblast['positive_fractions'][seq_idx[0]]
            pdict['Blast']['Ref_id'] = prot_name
            pdict['Blast']['fraction_pos_id'] = frac_pos
            break
            
count = 0
for p in mtprot_dict.keys():
    annotate_blast_results(p)
    if count%1000==0:
        print(count)
    count+=1
pickle.dump(mtprot_dict, open(aa_subs_dir+'MTProtein_isoform_dict.p', 'wb'))

 # create new ref_prot_dict with unique blast result and list of isoform names
unq_blast = list(set([list(pdict['Blast'].values())[0] for pdict in mtprot_dict.values() if len(pdict['Blast'].values())>0]))
#print(len(unq_blast))
blastprot_dict = {blast:{} for blast in unq_blast}
for p, pdict in mtprot_dict.items():
    blast = list(pdict['Blast'].values())
    if len(blast)>0:
        blast = blast[0]
        blastprot_dict[blast][p] = pdict
pickle.dump(blastprot_dict, open(aa_subs_dir+'MTProtein_ref_dict.p', 'wb'))
