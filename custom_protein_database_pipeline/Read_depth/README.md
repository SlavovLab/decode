# This directory contains code to get the read depth at each base from aligned RNA-seq data. 
# Getting started 

### Prerequisites 
RNA-seq data in .bam or .fastq format

If using CPTAC data, download RNA-seq data from GDC (NIH Genomic Data Commons) using GDC transfer tool client, gdc_manifest and gdc_user_token (restricted access required). Use gdc_sample_sheet (downloaded from GDC cart) to map file name to sample name. Sequencing files will be saved as FileID/FileName.bam

### Software package dependencies
1. samtools v1.10
2. fastp v0.23.4
3. hisat2 v2.2.1

# Usage
## 1. Run PGDB_N_reads.sh in terminal or Linux after installing dependencies. 

This script was written to be run in a Slurm environment on a single file downloaded from GDC. The preamble may need editing. Update directory to location of dataset RNA-seq data files. 

Input is FileID, FileName, SampleName where FileID specifies the directory of the specific RNA-seq data file in the dataset directory, FileName is the name of the .bam file (i.e. FileName.bam) and SampleName is the desired name of the output files (i.e. SampleName.fasta)

Example command line call to run script:
> PGDB_N_reads.sh FileID FileName SampleName

### Consideration for large datasets with many RNA-seq files 
To run multiple in-silico translations at once in a HPC environment, use gdc_sample_sheet to create id_list.txt (list of File IDs), file_list.txt (list of File Names), sample_list.txt (list of desired output sample names). Corresponding lines in each file should reference the same sample. 

To execute for all files in a dataset:
1. Create arrays on the command line from id_list.txt, file_list.txt, sample_list.txt. E.g.
> idArray=(); while IFS= read -r line; do idArray+=("$line"); done<id_list.txt
2. Loop over arrays to submit all jobs for dataset simultaneously. E.g.
> for i in {0..185}; do  sbatch --job-name=${sampleArray[$i]} --output=PGDB_${sampleArray[$i]}.out PGDB_main.sh ${idArray[$i]} ${fileArray[$i]} ${sampleArray[$i]}; done

## 2. Extract N reads pertaining to locations with identified SAAP with extract_SAAP_reads.py



