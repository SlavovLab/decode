#!/bin/bash


# Ensure R is installed
if ! command -v R &> /dev/null; then
    echo "R needs to be isntalled  and in your PATH."
    exit 1
fi


## Note, that all scripts also load functions defined in
## raas_utils.R

## Initialize BP/SAAP and RAAS TABLES
## NOTE, that this script is called redundantly from each script below.
## It loads all required libraries, and defines all input files and output
## data paths. This script also generates the
## "Supplemental_Data_7.SAAP_coordinates.tsv" as well as tables for
## unique protein sites used in the random forest model by Andrew Leduc.
R --vanilla <  raas_init.R &> $DECODE/log/raas_init.txt  

### CALCULATE RAAS PROFILES and GENERATE PUBLICATION FIGURES
R --vanilla < codons.R &> $DECODE/log/codons.txt  # Fig. 2j, EFig. 3k, 4
R --vanilla < aminoacids.R &> $DECODE/log/aminoacids.txt  # Fig. 3b,c,e, EFig. 5a,c,f
R --vanilla < function.R &> $DECODE/log/function.txt # Fig. 4c, EFig. 7
R --vanilla < motifs.R &> $DECODE/log/motifs.txt  # Fig. 4c,d, Fig. 5c, EFig. 8a,c
R --vanilla < kraq.R &> $DECODE/log/kraq.txt  # EFig. 8b
R --vanilla < proteins.R &> $DECODE/log/proteins.txt # Fig. 5a,b, EFig. 10a
R --vanilla < structure.R &> $DECODE/log/structure.txt  # EFig. 10b,c
R --vanilla < rna.R &> $DECODE/log/rna.txt # SFig. 2

## TODO: add 1d protein plots and chimeraX generating code
