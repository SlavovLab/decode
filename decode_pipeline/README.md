# Pipeline for identifying, validating and quantifying amino acid substitutions from alternate RNA decoding

## Getting started 

The python scripts in [python_scripts](python_scripts/) are used for identifying, validating and quantifying amino acid substitutions in results from a dependent peptide search in MaxQuant. 

The scripts are dependent on one another and intended to be run in a specific order. After identifying candidate peptides with substitutions, the user must run a validation search with MaxQuant or another search engine. 

### Prerequisites
1. Dependent peptide search in MaxQuant, ideally with sample-specific protein databases
2. Table with sample metadata (specific format for TMT-labeled CPTAC datasets)

### CPTAC sample metadata 
Columns: 
1. TMT plex (TMT experiment #, e.g. [1-22] for LSCC
2. TMT channel (e.g. 130N)
3. ParticipantID (patient study ID, e.g. C3L-02655)
4. Group (Tumor or Normal)
5. MQ (MaxQuant intensity column # corresponding to TMT channel defined in MQ_TMT_dict)
> MQ_TMT_dict = {'126':1,'127N':2, '127C':3, '128N':4, '128C':5, '129N':6, '129C':7, '130N':8, '130C':9, '131':10, ‘131C’:11}
6. Sample_name (ParticipantID_Group, e.g. C3L-02655_Tumor)
7. sample_ID (S+TMT plex, e.g. S1)
Example can be found in dataset directories in [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs)

# Usage 

## 1. Filter modified peptides from DP search to find candidate substituted peptides
### Run **AA_subs_detection_TMT.py** or **AA_subs_detection_labelfree.py** for TMT-labeled or label-free data, respectively. 

These scripts analyze DP search results to find peptides with known modifications (e.g. PTMs) and peptides with potential AAS.

### Dependencies:
1. allPeptides.txt (MaxQuant output)
2. evidence.txt (MaxQuant output)
3. sample_map.xlsx (dataset metadata - described above)
4. Modifications_table_091520.xlsx (Downloaded all known modifications from unimod.org, appended with addition and loss of TMT label. Can be found in [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs)
5. Genome fasta file

### Outputs:
1. DP_dict.p : Dictionary of allPeptides.txt annotated with modifications and/or AA substitution data (e.g. modification/substitution type, mass shift, new sequence, modification error, genome substring data)
2. PTM_dict.p : Subset of DP_dict corresponding to known PTMs (biological or experimental)
3. MTP_dict.p : Subset of DP_dict corresponding to candidate peptides with AAS
4. DP_search_evidence_dict.p, a dictionary with TMT sets as keys and evidence.txt dataframes as values.
5.  6-frame translations of genome used to assess if peptide can arise from noncanonical translation
These can be used as input if previously generated

Output files for data analyzed in [LINK TO PREPRINT] can be found in [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs)

## 2. Filter candidate peptides by 1% FDR and append identified candidate peptides to custom fasta files

### Run **AA_subs_validation1.py** 

This script computes q-values for candidate SAAP, appends candidate SAAP to custom fasta files for validation search.

### Dependencies:
1. sample_map.xlsx (dataset metadata)
2. MTP_dict.p (from AA_subs_detection.py)
3. Protein fasta databases (used for DP search)
4. Evidence.txt (MaxQuant DP output)

### Outputs:
1. DP_search_evidence_dict.p : A dictionary with TMT sets/samples as keys and corresponding DP search evidence.txt dataframes as values
2. qMTP_dict.p : MTP_dict.p annotated with posterior substitution probabilities and q-values and filtered for candidate SAAP that pass confidence threshold
3. custom_fasta_MTP.fasta : Custom fasta file with candidate SAAP appended

## 3. Run validation database search in MaxQuant
See MaxQuant_templates/README.md

## 4. Filter candidate substituted peptides to those that are identified in validation search 

### Run **AA_subs_validation2.py** 

This script generates a version of MTP_dict.p filtered for those that are quantified in the validation search and that have 2+ fragment ions providing evidence at site of substitution

### Dependencies:
1. sample_map.xlsx (dataset metadata)
2. qMTP_dict.p (from AA_subs_validation1.py)
4. Evidence.txt (MaxQuant validation search output)
5. Msms.txt (MaxQuant validation search output)

### Outputs:
1. Validated_MTP_dict.p : qMTP_dict.p filtered for those peptides quantified in validation search
2. Ion_validation_MTP_dict.p : Validated_MTP_dict.p filtered for peptides with 2+ fragment ions providing evidence for substitution
3. Validation_search_evidence_dict.p : A dictionary with TMT sets/samples as keys and corresponding validation search evidence.txt dataframes as values

## 5. Quantify validated substituted peptides

### Run **AA_subs_quant_TMT.py** or **AA_subs_quant_labelfree.py** 

These scripts quantifies validated substituted peptides and their corresponding base peptides and computes ratios between the 2. AA_subs_quant_TMT will quantify peptides at the precursor (MS1) and reporter ion (MS2) level, while the label-free version quantifies peptides only at the precursor level. 

### Dependencies:
1. sample_map.xlsx (dataset metadata)
2. Ion_validated_MTP_dict.p (from AA_subs_validation2.py)
4. Validation_search_evidence_dict.p (from AA_subs_validation2.py)

### Outputs:
1. MTP_quant_dict.p : This dictionary is structured differently from the previously generated dictionaries. It contains indices as keys, and each value is a dictionary for a unique SAAP-BP pair identified in the dataset. Each peptide dictionary contains quantification data for the peptide, including SAAP abundance, BP abundance, and RAAS, at the precursor (experiment)-level and reporter ion (patient sample)-level.


## Further analysis
The data files produced with this pipeline can be used for futher analysis of the substituted peptides. The results from running this pipeline as described in [Tsour *et al*](https://doi.org/10.1101/2024.08.26.609665) can be found in [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs) and used to generate the analyses and figures described in [decode_analysis](decode_analysis) and [decode_figures](decode_figures)

Additional scripts for quantifying the non-modified peptides, proteins, and modified proteins from the MaxQuant search data are provided in ~/decode_pipeline/python_script/downstream_processing/. The output from these scripts used in generating publication figures are also provided in [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs)

Other external analyses, such as mapping substitutions to structural domains with InterProScan and blasting all sample-specific protein sequences were also run. Scripts used to compile the output from those analyses are provided in ~/decode_pipeline/python_script/downstream_processing/scripts_for_processing_external_datasets/. The output from these scripts used in generating publication figures are also provided in [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs)

