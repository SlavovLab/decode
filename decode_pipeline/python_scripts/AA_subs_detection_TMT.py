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
import multiprocessing as mp

"""
    This code reads in results from DP search with MaxQuant and dataset metadata (see sample_map variable) from TMT-labeled data
    and outputs dictionaries of data for putative peptides with amino acid substitutions
"""

# setting directories and reading in data
code_dir = '/home/tsour.s/scripts_templates/AA_subs_pipeline/'
proj_dir = '/scratch/tsour.s/Vas_COAD/'
aa_subs_dir = proj_dir+'AA_subs_pipeline/'# aa_subs_dir is the output folder

""" sample map is the dataset metadata
    columns of sample map are [TMT plex, TMT channel, ParticipantID, Group, MQ, sample_name, sample_ID]
    defined as: TMT plex = tmt experiment #, TMT channel = (e.g.) 130N, ParticipantID = patient study id, Group = Tumor or Normal, MQ = MaxQuant intensity column number (see #MQ_TMT_dict below), sample_name = ParticipantID_Group, sample_ID = S+TMT Plex"""
sample_map = pd.read_excel(aa_subs_dir+'sample_map.xlsx')
samples = ['S'+str(i) for i in list(set(sample_map['TMT plex']))]

#MQ_TMT_dict = {'126':1,'127N':2, '127C':3, '128N':4, '128C':5, '129N':6, '129C':7, '130N':8, '130C':9, '131':10}
mod_df = pd.read_excel(code_dir+'Modifications_table/Modifications_table_091520.xlsx') #all unimod modifications plus TMT label and amidation
path_to_fasta = aa_subs_dir+'GRCh38_latest_genomic.fna'

""" these files can be read in if previously generated a 6-frame translation for species of interest"""
#path_to_w1_translation ='/scratch/tsour.s/human_DNA_fastas/all_frame1_W_aa_ambig.p'
#path_to_sa1_index ='/scratch/tsour.s/human_DNA_fastas/all_frame1_sa_ambig.p'
#path_to_w3_translation ='/scratch/tsour.s/human_DNA_fastas/all_frame3_W_aa_ambig.p'
#path_to_sa3_index ='/scratch/tsour.s/human_DNA_fastas/all_frame3_sa_ambig.p'
#path_to_w6_translation ='/scratch/tsour.s/human_DNA_fastas/all_frame6_W_aa_ambig.p'
#path_to_sa6_index ='/scratch/tsour.s/human_DNA_fastas/all_frame6_sa_ambig.p'


# amino acid molecular weight dictionary
MW_dict = {"G": 57.02147,
            "A" : 71.03712,
            "S" : 87.03203,
            "P" : 97.05277,
            "V" : 99.06842,
            "T" : 101.04768,
            "I" : 113.08407,
            "L" : 113.08407,
            "N" : 114.04293,
            "D" : 115.02695,
            "Q" : 128.05858,
            #"K" : 128.09497,
            "K" : 357.257902, # mass of lysine is adjusted to reflect TMT label
            "E" : 129.0426,
            "M" : 131.04049,
            "H" : 137.05891,
            "F" : 147.06842,
            "R" : 156.10112,
            "C" : 160.030654,
            "Y" : 163.0633,
            "W" : 186.07932,
            }
#list of amino acids
AAs = 'ACDEFGHIKLMNPQRSTVWY'

## dictionary of AAS types and the theoretical mass shift of the substitution
subs_dict = { i+' to '+j : MW_dict[j] - MW_dict[i] for i in MW_dict for j in MW_dict if i!=j}

## update subs_dict to account for the fact that I and L have the same mass
del subs_dict['L to I']
del subs_dict['I to L']
subst_dict={}
for k,v in subs_dict.items(): # unifies I and L
    if k[-1]!='I' and k[-1]!='L':
        subst_dict[k] = v
    elif k[-1]=='I':
        subst_dict[k+'/L']=v
subs_dict = {}
for a in AAs:
    subs_dict[a] = {k:v for k,v in subst_dict.items() if k[0]==a}
# sort AAS by mass shift
sorted_subs, sorted_sub_masses = zip(*sorted(subst_dict.items(), key= lambda x: x[1]))
sites = list(AAs)+['nterm','cterm']

""" Functions adapted from Mordret et al., Mol.Cell, 2019, to identify amino acid substitutions in MQ DP search results """
def prot_nterm(sequence):
    """
    Does the peptide originate at the protein's N-term
    """
    for start in SA_search(sequence, W_aa_ambiguous, sa_ambiguous):
        if W_aa_ambiguous[start-1] == '*':
            return True
        if W_aa_ambiguous[start-2] == '*':
            return True
    return False

def prot_cterm(sequence):
    """
    Does the peptide end at the protein's C-term
    """
    l=len(sequence)
    for start in SA_search(sequence, W_aa_ambiguous, sa_ambiguous):
        end = start+l
        if W_aa_ambiguous[end] == '*':
            return True
    return False

def pep_cterm(modified_sequence):
    """
    Returns the probability that C term AA was modified.
    """
    if modified_sequence[-1] == ')':
        prob = float(modified_sequence[:-1].split('(')[-1])
        if prob >= 0.95:
            return True
        else:
            return False
    else:
        return False

def pep_nterm(modified_sequence):
    """
    Returns the probability that N term AA was modified.
    """
    if modified_sequence[1] == '(':
        prob =  float(modified_sequence[2:].split(')')[0])
        if prob >= 0.95:
            return True
        else:
            return False
    else:
        return False

def is_gene(record):
    if len(record.seq)%3 != 0:
        return False
    if not record.seq[:3] in {'ATG'}:
        return False
    if record.seq[-3:].translate()!='*':
        return False
    return True

def codonify(seq):
    """
    input: a nucleotide sequence (not necessarily a string)
    output: a list of codons
    """
    seq = str(seq)
    l = len(seq)
    return [seq[i:i+3] for i in range(0,l,3)]

def suffix_array(text, _step=16):
    """Analyze all common strings in the text.
    
    Short substrings of the length _step a are first pre-sorted. Then are the
    results repeatedly merged so that the guaranteed number of compared
    characters bytes is doubled in every iteration until all substrings are
    sorted exactly.
    
    Arguments:
        text:  The text to be analyzed.
        _step: Is only for optimization and testing. It is the optimal length
               of substrings used for initial pre-sorting. The bigger value is
               faster if there is enough memory. Memory requirements are
               approximately (estimate for 32 bit Python 3.3):
                   len(text) * (29 + (_size + 20 if _size > 2 else 0)) + 1MB
    
    Return value:      (tuple)
      (sa, rsa, lcp)
        sa:  Suffix array                  for i in range(1, size):
               assert text[sa[i-1]:] < text[sa[i]:]
        rsa: Reverse suffix array          for i in range(size):
               assert rsa[sa[i]] == i
        lcp: Longest common prefix         for i in range(1, size):
               assert text[sa[i-1]:sa[i-1]+lcp[i]] == text[sa[i]:sa[i]+lcp[i]]
               if sa[i-1] + lcp[i] < len(text):
                   assert text[sa[i-1] + lcp[i]] < text[sa[i] + lcp[i]]
    >>> suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])
    
    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """
    tx = text
    size = len(tx)
    step = min(max(_step, 1), len(tx))
    sa = list(range(len(tx)))
    sa.sort(key=lambda i: tx[i:i + step])
    grpstart = size * [False] + [True]  # a boolean map for iteration speedup.
    # It helps to skip yet resolved values. The last value True is a sentinel.
    rsa = size * [None]
    stgrp, igrp = '', 0
    for i, pos in enumerate(sa):
        st = tx[pos:pos + step]
        if st != stgrp:
            grpstart[igrp] = (igrp < i - 1)
            stgrp = st
            igrp = i
        rsa[pos] = igrp
        sa[i] = pos
    grpstart[igrp] = (igrp < size - 1 or size == 0)
    while grpstart.index(True) < size:
        # assert step <= size
        nextgr = grpstart.index(True)
        while nextgr < size:
            igrp = nextgr
            nextgr = grpstart.index(True, igrp + 1)
            glist = []
            for ig in range(igrp, nextgr):
                pos = sa[ig]
                if rsa[pos] != igrp:
                    break
                newgr = rsa[pos + step] if pos + step < size else -1
                glist.append((newgr, pos))
            glist.sort()
            for ig, g in groupby(glist, key=itemgetter(0)):
                g = [x[1] for x in g]
                sa[igrp:igrp + len(g)] = g
                grpstart[igrp] = (len(g) > 1)
                for pos in g:
                    rsa[pos] = igrp
                igrp += len(g)
        step *= 2
    del grpstart
    del rsa
    return sa

def SA_search(P, W, sa):
    """ search for peptide in genome
    Input: P = peptide sequence, W = translated genome, sa = sorted substrings of proteome
    Output: list of indices of string match
    """
    lp = len(P)
    n = len(sa)
    l = 0; r = n
    while l < r:
        mid = int((l+r) / 2)
        a = sa[mid]
        if P > W[a : a + lp]: #either will be in latter half of substrings or former half
            l = mid + 1
        else:
            r = mid
    s = l; r = n
    while l < r:
        mid = int((l+r) / 2)
        a = sa[mid]
        if P < W[a : a + lp]:
            r = mid
        else:
            l = mid + 1
    return [sa[i] for i in range(s, r)] #s is either 0 or sa[mid+1], r is end of substring

def build_ref_translation(path_to_fasta, output_str, f=[1,2,3,4,5,6]):
    """ build reference fasta for removal of homologous sequences """
    """ f = frame """
    record_list = []
    translated_record_list = []
    record_dict = {}
    boundaries_aa = [0]
    
    for i,record in enumerate(SeqIO.parse(open(path_to_fasta,'r+'),'fasta')):
        if f in [1,2,3]:
            record.seq = record.seq.upper()[f:]
        else:
            record.seq = record.seq.upper()[::-1][f:]
        translation = str(record.seq.translate())
        bits = record.description.split(' ')
        record_list.append(record)
        translated_record_list.append(translation)
        boundaries_aa.append(boundaries_aa[-1]+len(translation))
    boundaries_aa = np.array(boundaries_aa[1:]) # an array annotating the genes' cumulative length
    W_aa = ''.join(translated_record_list)
    sa = suffix_array(W_aa)
    W_aa_ambiguous = W_aa.replace('I','L')
    sa_ambiguous = suffix_array(W_aa_ambiguous)

    pickle.dump(W_aa_ambiguous, open(output_str+'W'+str(f)+'_aa_ambig.p', 'wb'))
    pickle.dump(sa_ambiguous, open(output_str+'s'+str(f)+'a_ambig.p', 'wb'))

    return(W_aa_ambiguous, sa_ambiguous)

def find_homologous_peptide(P, W_aa_ambiguous, sa_ambiguous):
    """
    Gets a peptide and returns whether it has homolegous translation in the genome.
    """
    if len(SA_search('K' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
        return(True)
    elif len(SA_search('R' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
        return(True)
    elif len(SA_search('*' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
        return(True)
    else:
        return(False)


def refine_localization_probabilities(modified_seq):
    """
    Input: modified sequence (a string of AA with p of each to contain modification: APKV(.7)ML(.3)L means that V was modified with p = .7 and L with p = .3)
    Output: 2 lists: all candidate residues and their positional probabilities.
    """
    modified_sites = [modified_seq[m.start()-1] for m in re.finditer('\(',modified_seq) ]
    weights = [float(i) for i in re.findall('\(([^\)]+)\)',modified_seq)]
    
    return([modified_sites, weights])


def get_aa_subs(dp_dict, idx, tol=5, pp=0):
    """
    Find DP sequences with mass shift that is within tolerance of theoretical mass shift of AAS. Positional probability value must be reported by MQ for base residue consideration
    Input: idx = iterator value, tol = mass shift error tolerance, pp = positional probability threshold
    Output: AAS information
    """
    cand_residues = dp_dict['DP candidate residues'][idx]
    pos_probs = dp_dict['DP positional probabilities'][idx]
    DP_deltam = dp_dict['DP Mass Difference'][idx]
    DP_mz = dp_dict['m/z'][idx]
    mtol = DP_mz*(tol/1e6)
    
    cand_aa_origin = []
    aa_origin_index = []
    cand_aa_subs = []
    aa_subs_prob = []
    aa_subs_deltam = []
    for i,res in enumerate(cand_residues):
        res_dict = subs_dict[res]
        for s in res_dict:
            delta_m = res_dict[s]
            if ((DP_deltam > 0) and (delta_m > 0)) or ((DP_deltam < 0) and (delta_m < 0)):
                if (DP_deltam > delta_m - mtol) & (DP_deltam < delta_m + mtol) & (pp < pos_probs[i]):
                    cand_aa_origin.append(res)
                    aa_origin_index.append(i)
                    cand_aa_subs.append(s)
                    aa_subs_prob.append(pos_probs[i])
                    aa_subs_deltam.append(np.abs(DP_deltam-delta_m))


    return([cand_aa_origin, aa_origin_index, cand_aa_subs, aa_subs_prob, aa_subs_deltam])


def get_mistranslated_seq(dp_dict, idx):
    """
    Input: idx = iterator value as loop over data dictionary
    Output: list of potential sequences that could explain mass shift
    """
    DP = dp_dict['DP Probabilities'][idx]
    base_seq = dp_dict['DP Base Sequence'][idx]
    cand_residues = dp_dict['DP candidate residues'][idx]
    pos_probs = dp_dict['DP positional probabilities'][idx]
    cand_res_seq_idx = [x.start()-1 for x in re.finditer('\(', DP)]
    sub_origin_idx = dp_dict['origin aa index'][idx]
    sub_dest = dp_dict['destination aa'][idx]
    mps = []
    for i, res in enumerate(cand_residues):
        cand_res_idx = cand_res_seq_idx[i]
        if i in sub_origin_idx:
            sub_idx = sub_origin_idx.index(i)
            dest_aa = sub_dest[sub_idx]
            mod_seq = DP[:cand_res_idx]+dest_aa+DP[cand_res_idx+1:]
            mod_seq = re.sub('\(([^\)]+)\)', '', mod_seq)
            mps.append(mod_seq)
    return(mps)


def find_PTMs(dp_dict, idx, tol=10):
    """
    finds and annotates DPs that can be explained by PTM
    """
    cand_residues = dp_dict['DP candidate residues'][idx]
    pos_probs = dp_dict['DP positional probabilities'][idx]
    DP_deltam = dp_dict['DP Mass Difference'][idx]
    DP_mz = dp_dict['m/z'][idx]
    prot_nterm = dp_dict['Protein N-term'][idx]
    prot_cterm = dp_dict['Protein C-term'][idx]
    pep_nterm = dp_dict['Peptide N-term'][idx]
    pep_cterm = dp_dict['Peptide C-term'][idx]
    
    mtol = DP_mz*(tol/1e6)
    PTM = []
    PTM_aa = []
    PTM_prob = []
    PTM_deltam = []
    for i, res in enumerate(cand_residues):
        #7/21 - added N-term and C-term sites to res_df and rerun
        res_df = mod_df.loc[(mod_df.site == res) | (mod_df.site == 'N-term') | (mod_df.site == 'C-term'), :]
        
        for j, row in res_df.iterrows():
            mod = row['modification']
            pos = row['position']
            delta_m = row['mass shift']
            
            if (DP_deltam > delta_m - mtol) & (DP_deltam < delta_m + mtol):
                if pos == 'Anywhere':
                    term_filter = True
                elif (pos == 'Protein N-term') & (prot_nterm==True):
                    term_filter = True
                elif (pos == 'Any N-term') & (pep_nterm==True):
                    term_filter = True
                elif (pos == 'Protein C-term') & (prot_cterm==True):
                    term_filter = True
                elif (pos == 'Any C-term') & (pep_cterm==True):
                    term_filter = True
                else:
                    term_filter=False
                
                if term_filter == True:
                    PTM.append(mod)
                    PTM_aa.append(res)
                    PTM_prob.append(pos_probs[i])
                    PTM_deltam.append(DP_deltam - delta_m)
    return([PTM, PTM_aa, PTM_prob, PTM_deltam])


""" build reference fasta for removal of homologous sequences """
W1_aa_ambiguous, sa1_ambiguous = build_ref_translation(path_to_fasta, aa_subs_dir, 1)
W2_aa_ambiguous, sa2_ambiguous = build_ref_translation(path_to_fasta, aa_subs_dir, 2)
W3_aa_ambiguous, sa3_ambiguous = build_ref_translation(path_to_fasta, aa_subs_dir, 3)
W4_aa_ambiguous, sa4_ambiguous = build_ref_translation(path_to_fasta, aa_subs_dir, 4)
W5_aa_ambiguous, sa5_ambiguous = build_ref_translation(path_to_fasta, aa_subs_dir, 5)
W6_aa_ambiguous, sa6_ambiguous = build_ref_translation(path_to_fasta, aa_subs_dir, 6)

""" if previously built 6-frame translation, comment out the above and read in files here"""
#W_aa_ambiguous = W1_aa_ambiguous = pickle.load(open(path_to_w1_translation , 'rb'))
#sa_ambiguous = sa1_ambiguous = pickle.load(open(path_to_sa1_index, 'rb'))
#W3_aa_ambiguous = pickle.load(open(path_to_w3_translation , 'rb'))
#sa3_ambiguous = pickle.load(open(path_to_sa3_index, 'rb'))
#W6_aa_ambiguous = pickle.load(open(path_to_w6_translation , 'rb'))
#sa6_ambiguous = pickle.load(open(path_to_sa6_index, 'rb'))


fout = open('aas_detect.out', 'w')
def get_data_dict(s):
    """ Function to identify and annotate peptides with AAS and other canonical PTMs """
    print(s)
    sfout = open(s+'_aas_detect.out', 'w')
    sfout.write(s+'\n')

        # read in MQ DP search results, allPeptides.txt file
    sample_df = pd.read_csv(proj_dir+'MQ_output/'+s+'_DP/combined/txt/allPeptides.txt', '\t', low_memory=False)
    allPep_count = len(sample_df)
    
        # filter allPeptides for peptides with DP that has PEP value<0.01
    # keep only relevant columns
    dp_df = sample_df.loc[~np.isnan(sample_df['DP Mass Difference']) & (sample_df['DP PEP']<=0.01), :]
    cols2keep = ['Raw file', 'Charge', 'm/z', 'Mass', 'Mass precision [ppm]' 'Retention time',
                'Sequence', 'Proteins','Intensity', 'DP Mass Difference', 'DP Time Difference', 'DP PEP',
              'DP Base Sequence', 'DP Probabilities', 'DP Positional Probability', 'DP Base Raw File', 'DP Base Scan Number',
              'DP Mod Scan Number', 'DP Proteins','DP Ratio mod/base']
    dp_df = dp_df.loc[:,cols2keep]
    dp_df = dp_df.reset_index(drop=True)
    dp_count = len(dp_df)
  
        # convert dataframe to dictionary
    dp_dict = dp_df.to_dict()
    dp_dict['count allPeptides'] = allPep_count
    dp_dict['count DP'] = dp_count
   
        # for each DP, identify the candidate base residues and their probabilities of being modified
    candidate_res_dict = {i: refine_localization_probabilities(v) for i,v in dp_dict['DP Probabilities'].items()}
    dp_dict['DP candidate residues'] = {i: candidate_res_dict[i][0] for i in candidate_res_dict.keys()}
    dp_dict['DP positional probabilities'] = {i: candidate_res_dict[i][1] for i in candidate_res_dict.keys()}
    dp_dict['count candidate residues per peptide'] = {i: len(v) for i,v in dp_dict['DP candidate residues'].items()}
   
    dp_dict['Protein N-term'] = {i: prot_nterm(v) for i,v in dp_dict['DP Base Sequence'].items()}
    dp_dict['Protein C-term'] = {i: prot_cterm(v) for i,v in dp_dict['DP Base Sequence'].items()}
    dp_dict['Peptide N-term'] = {i: pep_nterm(v) for i,v in dp_dict['DP Probabilities'].items()}
    dp_dict['Peptide C-term'] = {i: pep_cterm(v) for i,v in dp_dict['DP Probabilities'].items()}
   
        # Find and annotate DPs that can be explained by AAS
    print('Finding AA subs')
    sfout.write('Finding AA subs\n')
    aa_subs_dicts = {i: get_aa_subs(dp_dict, i) for i in dp_dict['DP candidate residues'].keys()}
    dp_dict['origin aa'] = {i: aa_subs_dicts[i][0] for i in aa_subs_dicts.keys()}
    dp_dict['origin aa index'] = {i: aa_subs_dicts[i][1] for i in aa_subs_dicts.keys()}
    dp_dict['aa subs'] = {i: aa_subs_dicts[i][2] for i in aa_subs_dicts.keys()}
    dp_dict['aa subs positional probability'] = {i: aa_subs_dicts[i][3] for i in aa_subs_dicts.keys()}
    dp_dict['aa subs mass error (ppm)'] = {i: aa_subs_dicts[i][4] for i in aa_subs_dicts.keys()}
    dp_dict['destination aa'] = {i: [x[-1] for x in v] if len(v)>0 else [] for i,v in dp_dict['aa subs'].items()}
    dp_dict['mistranslated sequence'] = {i: get_mistranslated_seq(dp_dict, i) for i,v in dp_dict['aa subs'].items()}
    mp_count = sum([1 for x in list(dp_dict['origin aa'].values()) if len(x)>0])
    dp_dict['count mistranslated peptides'] = mp_count
    dp_dict['count aa subs per peptide'] = {i: len(v) for i,v in dp_dict['origin aa'].items()}

        # Find and annotate DPs that can be explained by canonical PTM
    print('Finding PTMs')
    sfout.write('Finding PTMs\n')
    PTM_dict = {i: find_PTMs(dp_dict, i) for i in dp_dict['DP candidate residues'].keys()}
    dp_dict['PTM'] = {i: PTM_dict[i][0] for i in PTM_dict.keys()}
    dp_dict['PTM site'] = {i: PTM_dict[i][1] for i in PTM_dict.keys()}
    dp_dict['PTM positional probability'] = {i: PTM_dict[i][2] for i in PTM_dict.keys()}
    dp_dict['PTM mass error [observed-expected] (ppm)'] = {i: PTM_dict[i][3] for i in PTM_dict.keys()}
   
        # Annotate whether or not peptide with AAS has homologous sequence that could have arisen from elsewhere in genome
    dp_dict['1-frame genome substring'] = {i: [find_homologous_peptide(x, W1_aa_ambiguous, sa1_ambiguous) for x in v] for i,v in dp_dict['mistranslated sequence'].items()}
    dp_dict['2-frame genome substring'] = {i: [find_homologous_peptide(x, W2_aa_ambiguous, sa2_ambiguous) for x in v] for i,v in dp_dict['mistranslated sequence'].items()}
    dp_dict['3-frame genome substring'] = {i: [find_homologous_peptide(x, W3_aa_ambiguous, sa3_ambiguous) for x in v] for i,v in dp_dict['mistranslated sequence'].items()}
    dp_dict['4-frame genome substring'] = {i: [find_homologous_peptide(x, W4_aa_ambiguous, sa4_ambiguous) for x in v] for i,v in dp_dict['mistranslated sequence'].items()}
    dp_dict['5-frame genome substring'] = {i: [find_homologous_peptide(x, W5_aa_ambiguous, sa5_ambiguous) for x in v] for i,v in dp_dict['mistranslated sequence'].items()}
    dp_dict['6-frame genome substring'] = {i: [find_homologous_peptide(x, W6_aa_ambiguous, sa6_ambiguous) for x in v] for i,v in dp_dict['mistranslated sequence'].items()}

    sfout.close()    
    return(dp_dict)

""" split AAS and PTM search for each sample/TMT with multiprocessing"""
p = mp.Pool()
data_dict_list = p.map(get_data_dict, samples)

p.close()
p.join()

data_dict = {s:data_dict_list[i] for i,s in enumerate(samples)}
pickle.dump(data_dict, open(aa_subs_dir+'DP_dict.p', 'wb'))

""" "split data_dict into dict for DPs with canonical PTMs (for potential further analysis outside scope of this study) and dict with filtered MTPs """
fout.write('filtering data dict\n')
mtp_dict = {}
ptm_dict = {}

#mtp_dict = pickle.load(open(aa_subs_dir+'MTP_dict.p', 'rb'))
#ptm_dict = pickle.load(open(aa_subs_dir+'PTM_dict.p', 'rb'))
                                                                                                             
for s in samples:
    if s not in mtp_dict.keys():
        dp_dict = data_dict[s]
        fout.write(s)
        curr_mtp_dict = {}
        curr_ptm_dict = {}
        ptm_idx = [i for i,x in dp_dict['DP PEP'].items() if len(dp_dict['PTM'][i])>0]
    
        # MTP dict are peptides with potential AAS that cannot be explained by PTM and do not have homologous sequence in genome
        mtp_idx = [i for i, x in dp_dict['DP PEP'].items() if (len(dp_dict['PTM'][i])==0) and (len(dp_dict['aa subs'][i])>0) and (dp_dict['1-frame genome substring'][i][0]==False) and (dp_dict['2-frame genome substring'][i][0]==False) and (dp_dict['3-frame genome substring'][i][0]==False) and (dp_dict['4-frame genome substring'][i][0]==False) and (dp_dict['5-frame genome substring'][i][0]==False) and (dp_dict['6-frame genome substring'][i][0]==False)]
        for k,v in dp_dict.items():
            if k!= 'count mistranslated peptides' and k!='count allPeptides' and k!='count DP':
                curr_mtp_dict[k] = {i:v[i] for i in mtp_idx}
                curr_ptm_dict[k] = {i:v[i] for i in ptm_idx}
        mtp_dict[s] = curr_mtp_dict
        ptm_dict[s] = curr_ptm_dict
        print('N mtps = '+str(len(mtp_dict[s]['mistranslated sequence'].keys())))
    
fout.write('saving files\n')
fout.close()
pickle.dump(mtp_dict, open(aa_subs_dir+'MTP_dict.p', 'wb'))
pickle.dump(ptm_dict, open(aa_subs_dir+'PTM_dict.p', 'wb'))
