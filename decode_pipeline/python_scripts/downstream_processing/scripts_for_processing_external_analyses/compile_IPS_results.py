import pandas as pd
import numpy as np
import pickle
from glob import glob
import re
import scipy as sp
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

""" This script will take in the results of mapping protein sequences containing AAS to known biological/structural domains using InterProScan (https://www.ebi.ac.uk/interpro/download/InterProScan/)
Input:
    1. MTProtein_isoform_dict.p (output of AAS_protein_quant.py)
    2. MTProt_InterProScan.tsv (output by InterProScan)
Returns:
    1. IPS_results.p : dictionary with compiled InterProScan results
    2. MTProtein_isoform_dict.p annotated with InterProScan results
    3. InterProScan_domains_toProteins_map.xlsx : excel sheet mapping proteins to domains
"""

proj_dir = '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Papers/Sat_LSCC/'
aa_subs_dir = proj_dir + 'AA_subs_pipeline/'
mtprot_dict = pickle.load(open(aa_subs_dir+'MTProtein_isoform_dict.p', 'rb'))
ips_df = pd.read_csv(aa_subs_dir+'MTProt_InterProScan.tsv', '\t',header=None, names=list(header))

# compile interpro results into dict
header=['Accession', 'Sequence MD5 digest', 'Sequence length', 'Analysis',
       'Signature accessioon', 'Signature description', 'Start location',
       'Stop location', 'Score', 'Status', 'Date',
       'InterPro annotations - accession',
       'InterPro annotations - description', 'GO annotations'] #determined from InterPro documentation

ips_dict = {}
for i, row in ips_df.iterrows():
    db = row['Analysis']
    if db not in ips_dict.keys():
        ips_dict[db] = {}
    desc = row['Signature description']
    if desc not in ips_dict[db].keys():
        ips_dict[db][desc] = {'Proteins':{},'Domain_MTP_count':0}
    ips_dict[db][desc]['Proteins'][row['Accession']] = row.to_dict()
pickle.dump(ips_dict, open(aa_subs_dir+'IPS_results.p', 'wb'))

# link ips results to allprot_dict
prot_ips_rows = []
count_fail = 0
for analysis, dom_dict in ips_dict.items():
    print(analysis)
    for dom in dom_dict.keys():
        prot_dicts = dom_dict[dom]['Proteins']
        #print(dom)
        prot_ips_rows.append([analysis, dom, list(prot_dicts.keys())])
        for p, prot_dict in prot_dicts.items():
            try:
                pdict = mtprot_dict[p]
                pdict['InterProScan'] = {analysis:{'Signature description':prot_dict['Signature description'],
                                                      'Signature accession':prot_dict['Signature accessioon']}}
            except:
                count_fail+=1
                
prot_ips_df = pd.DataFrame(prot_ips_rows, columns=['Analyis', 'Domain', 'Proteins'])
prot_ips_df.to_excel(aa_subs_dir+'InterProScan_domains_toProteins_map.xlsx')

pickle.dump(mtprot_dict, open(aa_subs_dir+'MTProtein_isoform_dict.p', 'wb'))
