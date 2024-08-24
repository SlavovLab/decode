# decode
Computational pipeline for quantifying amino acid substitutions from alternate RNA decoding

The code in this repository can be used to identify, validate, and quantify amino acid substitutions in LC-MS proteomics data that arise from alternate decoding of RNA. 

The pipeline has multiple steps and multiple sections of the pipeline require input from external analysis, such as database searches of proteomics data. AS such, this is not a standalone software package that can be run automatically. Detailed instructions for running the pipeline are provided in the various README.md files in this repository.

## Overview of pipeline
1. Custom protein sequence database generation by *in-silico* translation of RNA-seq data
2. Dependent peptide search with MaxQuant (external software analysis)
3. Identifying candidate peptides with amino acid substitutions
4. Validation database search (external software analysis)
5. Quantifying validated peptides with amino acid substitutions
6. Downstream data analysis

## Usage
### Step 1: Custom protein databases
Use RNA-seq data matched to LC-MS proteomics data to create sample-specific protein databases. 
The code for this step is in custom_protein_database_pipeline/ and the README.md in that directory contains detailed instructions for running the code.
If no matched RNA-seq data is available, this step can be skipped, but caution should be taken in interpreting quantified amino acid substitutions as there is lower confidence that they are not encoded in the genome. 
