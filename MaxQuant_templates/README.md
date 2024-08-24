This directory contains sample MaxQuant parameter templates for searching raw MS-proteomics data with MaxQuant dependent peptide search and standard database search. 

gen_mqpar.py is a script that can be used to generate new MaxQuant parameter files with user-specified raw files, fasta database, and other parameters using the provided files as templates. 

## Running the dependent peptide search in MaxQuant
1. Download raw proteomics data from relevant databases e.g. with wget
2. Edit MQ parameter file using gen_mqpar.py with file paths of custom protein database and raw proteomics data files. Ensure dependent peptide parameter is set to True.
3. Run MQ with edited parameter file as input. For a good tutorial on running MQ in Linux see [MaxQaunt in Linux](https://atchen.me/research/2019/03/21/mq-linux.html)

### Notable MaxQuant parameter settings for dependent peptide search
- MQ v1.6.17.0
- dependentPeptides = True
- dependentPeptideFDR = 0.01
- Met oxidation and N-term acetylation variable modifications
- If searching TMT data, set lcmsRunType = Standard and add TMT labels as fixed modifications
- Cys carbamidomethyl fixed modification         
