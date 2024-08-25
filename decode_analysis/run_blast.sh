#!/bin/bash

## ANALYSIS OF AMINO ACID SUBSTITUTION SITES
##  as provided by Shiri Tsour.

## SOFTWARE DEPENDENCIES

## additionally to standard bash tools (cat, grep, cut, awk, gunzip,
## rsync), this script additionally requires NCBI's blast, and a
## recent R version. 

if [[ -z "${DECODE}" ]]; then
    echo "Please define a path for input and output of this script:" 
    echo "export DECODE=<YOURDATAPATH>"
fi

## NUMBER OF CPUs used for blast
nthreads=7

## OUTPUT DATA PATHS
## SET THIS PATH TO WHERE YOU WANT TO SAVE INPUT and OUTPUT

## generate output paths
mkdir -p $DECODE/log
mkdir -p $DECODE/figures
mkdir -p $DECODE/processedData
mkdir -p $DECODE/originalData

### DATA COLLECTION

## 1) REQUIRED MANUAL DOWNLOADS

## download additionalData.zip which provides human genome data files
## generated in by the genomeBrowser project AND unzip in $DECODE/

## cp -a ~/data/mistrans/additionalData.zip $DECODE/
## cd $DECODE; unzip additionalData.zip

echo "testing existence of required input data"

required=${DECODE}/additionalData/features_GRCh38.110.tsv.gz
if [ ! -f "$required" ]; then
    echo "MANUAL STEPS REQIRED: Please download the additionalData.zip"
    echo "and unzip as instructed."
    exit 1
fi

echo "testing availabilty of NCBI blast"

# Ensure blast is installed
if ! command -v blastp &> /dev/null; then
    echo "blastp needs to be isntalled  and in your PATH."
    exit 1
fi
# Ensure R is installed
if ! command -v R &> /dev/null; then
    echo "R needs to be isntalled  and in your PATH."
    exit 1
fi


## 2) COLLECT BP/SAAP SEQUENCES
## produces $DECODE/processedData/unique_bp.fas and
##          $DECODE/processedData/unique_saap.tsv
echo "extracting BP and SAAP sequences"

## write BP into a fasta file for use with blast
gunzip -c $DECODE/additionalData/All_SAAP_TMTlevel_quant_df.txt.gz \
    | cut -f 5 | sort | uniq | grep -v "^BP$" | \
    grep -v -e '^$' > $DECODE/processedData/tmp.txt
gunzip -c $DECODE/additionalData/All_SAAP_patient_level_quant_df.txt \
    | cut -f 5 | sort | uniq | grep -v "^BP$" | \
    grep -v -e '^$' > $DECODE/processedData/tmp2.txt
cat $DECODE/processedData/tmp.txt $DECODE/processedData/tmp2.txt | \
    sort |uniq | awk '{print ">" $0 ORS $0}' - \
		     > $DECODE/processedData/unique_bp.fas
rm -f $DECODE/processedData/tmp.txt $DECODE/processedData/tmp2.txt


## BP/SAAP as simple table, basis for search in proteins
gunzip -c $DECODE/additionalData/All_SAAP_TMTlevel_quant_df.txt.gz \
    | cut -f 4,5 | sort | uniq | grep -v -P "^SAAP\tBP" | \
    grep -v -e '^$' > $DECODE/processedData/tmp.tsv
gunzip -c $DECODE/additionalData/All_SAAP_patient_level_quant_df.txt \
    | cut -f 4,5 | sort | uniq | grep -v -P "^SAAP\tBP" | \
    grep -v -e '^$' > $DECODE/processedData/tmp2.tsv
cat $DECODE/processedData/tmp.tsv $DECODE/processedData/tmp2.tsv | \
    sort | uniq > $DECODE/processedData/unique_saap.tsv
rm -f $DECODE/processedData/tmp.txt $DECODE/processedData/tmp2.txt


## 3) GET Ensembl PROTEIN SEQUENCES
echo "downloading Ensembl protein fasta"
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz $DECODE/originalData/

## 4) ADD PATIENT-SPECIFIC SINGLE AMINO REPLACEMENTS
## produces $DECODE/processedData/all_proteins.fa
echo "adding mutations to fasta"
R --vanilla < get_mutated_proteins.R  &> $DECODE/log/mutated_proteins.txt 


## 5) BLAST BP IN Ensembl+MUTATION PROTEINS
## produces $DECODE/processedData/unique_bp_blast.tsv
echo "RUNNING BLAST SEARCH OF BP"

makeblastdb -in $DECODE/processedData/all_proteins.fa -parse_seqids \
		      -title "ensembl hg38 proteins" -dbtype prot  \
    &> $DECODE/log/blast_database.txt
## NOTE: blast hits are pre-filtered for query length=alignment length
##       and >75% identity
format="6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"
blastp  -num_threads ${nthreads} -task blastp-short \
	   -query  $DECODE/processedData/unique_bp.fas \
	   -db $DECODE/processedData/all_proteins.fa  \
	   -outfmt "$format" 2> $DECODE/log/unique_bp_blast.txt \
    | awk '{if($5==$6 && $3>75) print}' \
    | grep -v "^#" > $DECODE/processedData/unique_bp_blast.tsv 


## 6) SELECT BEST-MATCHING PROTEIN HIT FOR EACH BP
## produces $DECODE/processedData/bp_mapped.tsv
echo "getting best blast hit for each BP"

R --vanilla < get_protein_match.R &> $DECODE/log/protein_match.txt

## 7) COLLECT ALL DATA FOR UNIQUE BP/SAAP PAIRS
## produces $DECODE/processedData/saap_mapped.tsv

## NOTE: this script requires large additional data produced
## by genomeBrowser project and partially
## requiring high performance computing. The script can be run
## after doing these calculations and providing the correct path
## to the resulting data directory; see
## https://gitlab.com/raim/genomeBrowser/-/blob/master/data/mammary/setup.sh

## Here, we provide the produced file saap_mapped.tsv in the
## downloaded additionalData folder instead, to allow running
## the R analysis scripts that produce the final published figures, below.
echo "collecting BP protein, transcript and genome coordinates" \
     "and various structural data"

R --vanilla < map_peptides.R &> $DECODE/log/map_peptides.txt



