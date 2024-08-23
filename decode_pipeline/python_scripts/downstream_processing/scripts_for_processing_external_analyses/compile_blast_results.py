#!/usr/bin python

import sys
from Bio.Blast import NCBIXML
import pickle
import pandas as pd
import re
from glob import glob
import os
import copy

""" This script takes in the output of a blastp search of substituted peptides against human RefSeq protein database, parses it, and compiles a dictionary of the top hit for each SAAP. It also creates a mapping between each main protein (ensembl), the blast reference protein, uniprot and gene names using uniprot data.

Input:
    1. directory containing blast results
    2. Allprot_normalized_abundance_dict.p (output of allprot_quant.py - a dictionary containing the abundances of all the proteins found in the validation search)
    3. ensebml_map; uniprot_map (files containing ID mapping data from uniprot)
    4. directory containing sample-specific protein fasta (output of AA_subs_validation1.py)
    
Output:
    1. TMT set/label-free sample specific compilations of top blast result for each protein sequence
    2. Dataset compilation of output 1.
    3. map of blast reference protein to protein name (fasta header) for each TMT set
    4. map of blast reference protein to protein name (fasta header), ensembl ID, uniprot ID and gene ID for dataset
    
Some parts of this code (indicated below) are more efficiently run for each TMT set or label-free sample separately, on a HPC or cluster environment
"""


proj_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Papers/Sat_LSCC/'
blast_dir = proj_dir+'blast/'
aa_subs_dir = proj_dir + 'AA_subs_pipeline/'
samples = ['S'+str(i) for i in list(range(1,26))]
genome_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/genomes/'


ensembl_map = pd.read_table(genome_dir + 'HUMAN_9606_idmapping_selected.tab') # human protein ID map downloaded from Uniprot. supplied in pipeline files.
ensembl_map_vals = ensembl_map.values
uniprot_map = pickle.load(open(genome_dir + 'human_uniprot_10Oct2022_dict.p', 'rb')) # a dictionary of human uniprot database entries. supplied in pipeline files.

allprot_abundance_dict = pickle.load(open(aa_subs_dir+'Allprot_normalized_abundance_dict.p', 'rb'))

### compile top blast result for each protein sequence into a dictionary for each sample
# to speed up, run separately for each sample, on cluster
for s in samples:
    set_dict = {'title':[],'length':[], 'id_fractions':[], 'positive_fractions':[],'gap_fractions':[], 'query':[]}
    blast_files = glob(proj_dir+'blast/'+s+'_*.xml')
    for blast in blast_files:
        print(blast)
        blast_handle = open(blast)
        blast_records = NCBIXML.parse(blast_handle)
        try:
            for record in blast_records:
                for alignment in record.alignments:
                    title = alignment.title
                    length = alignment.length
                    hsp = alignment.hsps[0]
                    ids = hsp.identities
                    id_fraction = ids/length
                    pos_id = hsp.positives
                    pos_fraction = pos_id/length
                    gaps = hsp.gaps
                    gap_fraction = gaps/length
                    query=hsp.query

                    set_dict['title'].append(title)
                    set_dict['length'].append(length)
                    set_dict['id_fractions'].append(id_fraction)
                    set_dict['positive_fractions'].append(pos_fraction)
                    set_dict['gap_fractions'].append(gap_fraction)
                    set_dict['query'].append(query)
            blast_dict[s] = set_dict
        except:
            #print(blast)
            pass
    pickle.dump(set_dict, open(blast_dir+s+'_blast_dict.p', 'wb'))

### combine all compiled results from each sample/TMT set into a single dictionary
blast_dict = {}
for s in samples:
    set_dict = pickle.load(open(blast_dir+s+'_blast_dict.p', 'rb'))
    blast_dict[s] = set_dict
pickle.dump(blast_dict, open(aa_subs_dir+'blast_dict.p', 'wb'))


### map blast reference protein to fasta protein names via their sequence
# to speed up, run separately for each sample, on cluster

# function to get protein sequence
def prot_details(prot_name, s_fasta):
    """
    input: protein name as in fasta header; TMT/label-free sample-specific fasta
    output: protein sequence
    """
    seqidx = [i for i, line in enumerate(s_fasta) if (prot_name in line) or (all(x in line for x in prot_name.split('|'))) or (all(x in line for x in prot_name.split('_')))]
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


# function to get the blast reference protein id from prot name
def blast_prot(blast_dict, prot, s, s_fasta):
    '''
    input: blast_dict.p, protein name as in fasta header; TMT or label-free sample name; TMT/label-free sample-specific fasta
    output: blast reference protein
    '''
    prot_seq = prot_details(prot, s_fasta).split('\n')[0]
    prot_name = ''
    sblast = blast_dict[s]
    seq_idx = [i for i,x in enumerate(sblast['query']) if x==prot_seq]
    if len(seq_idx)>0:
        prot_name = sblast['title'][seq_idx[0]]
    return(prot_name)
    

for s in samples:
    s_fasta = open(proj_dir+'databases/'+s+'_MTP.fasta', 'r').readlines()
    s_fasta = [x for x in s_fasta if x!= '\n']
    prot_df = allprot_abundance_dict[s]
    blast_names = []
    prot_names = []
    for i, row in prot_df.iterrows():
        prots = row.name
        prots = prots.split(';')
        for prot in prots:
            prot_name = blast_prot(blast_dict, prot, s, s_fasta):
            hit=0
            if len(prot_name)>0:
                blast_names.append(prot_name)
                prot_names.append(prot)
                hit=1
                break
        if hit==0:
            blast_names.append('Unknown_'+prots[0])
            prot_names.append(prots[0])
            
blast_name_map = pd.DataFrame(zip(prot_names, blast_names), columns=['Protein name', 'Blast protein'])
blast_name_map.to_excel(blast_dir+s+'_blast_map.xlsx')


### compile all blast_maps into one per dataset and map blast protein to Ensembl, Uniprot and gene name
blast_maps = glob(blast_dir+'*_blast_map.xlsx')

for blast_map in blast_maps:
    s = blast_map.split('/')[-1].split('_')[0]
    bmap = pd.read_excel(blast_map, index_col=0)
    all_bmaps.append(bmap)

all_bmaps_df = pd.concat(all_bmaps)
all_bmaps_df = all_bmaps_df.drop_duplicates()

all_ds_maps = {}
gene_names = []
uniprot_ids = []
for i,row in all_bmaps_df.iterrows():
   # if i%1000==0:
    #    print(i)
    prot_genes = []
    prot_up = []
    protein = row['Protein name']
    if 'ENSP' in protein:
        ens = protein.split('_')[0]
        ens_rows = [k for k,i in enumerate(ensembl_map_vals) if (isinstance(i[-2], str)) and (ens in i[-2])]
        ens_df = ensembl_map.loc[ens_rows]
        uniprot_acc = ens_df['UniProtKB-AC'].values
        for up in uniprot_acc:
            prot_up.append(uniprot_acc)
            if up in uniprot_map.keys():
                [prot_genes.append(x) for x in uniprot_map[up]['Gene_names']]
    gene_names.append(prot_genes)
    uniprot_ids.append(prot_up)
all_bmaps_df['Gene names'] = gene_names
all_bmaps_df['Uniprot IDs'] = gene_names
all_bmaps_df.to_excel(aa_subs_dir_dir+'blast_map_w_gene.xlsx')

