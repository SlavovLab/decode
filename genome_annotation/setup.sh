#!/bin/bash

## DATA SETUP and CALCULATIONS for ANALYSIS of AMINO ACID SUBSTITUTIONS

## NOTE that these scripts were extracted from their original place in a
## larger multi-project script base by Rainer Machne ((C) 2025), accessible
## at https://gitlab.com/raim/genomeBrowser/-/tree/master/data/mammary

## This is a log file more than a script, manual steps are
## required in all cases within `if false;` conditions:
## some data downloads, data transfer to and calculations at a HPC platform.

## SETUP DATA DIRECTORIES FOR HUMAN GENOME DATA

## PATH WHERE DATA IS STORED - NOTE THAT THIS PATH IS
## REQUIRED TO BE DEFINED IN R scripts
export MAMDATA=/home/{$USER}/data/mistrans_data/ 

## PATH TO SCRIPTS (e.g. this one; git@github.com:SlavovLab/decode.git
## subdirectory genome_analysis)
SRC=<SOURCECODE_PATH>

## required additional tools
bioawk=~/scripts/bioawk

## main directory structure
ORIGDATA=$MAMDATA/originalData
PROCDATA=$MAMDATA/processedData

## generate directory structure
mkdir $MAMDATA
mkdir $ORIGDATA
mkdir $PROCDATA
mkdir $MAMDATA/chromosomes
mkdir $MAMDATA/log


## download data to originalData
## includes trivial preprocessing: format conversions, data extractions, etc.
## NOTE: do this manually, since many files
## are behind wget/rsync firewalls
if false; then
    $SRC/download.sh
fi

## ANALYZE GFF3 FILE STRUCTURE

## count of annotated features
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | cut -f 3  | grep -v "\\#"  | sort| uniq -c \
    | sort -nr > $PROCDATA/features_count.txt

## PROBLEM: ca. 10k more proteins in protein fasta than mRNA in gff3.
## These are annotated to scaffold sequences instead of properly assigned
## chromosome positions.

## analyze gff3 vs. protein fasta
## get all proteins in protein fasta
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.pep.all.fa.gz \
    | grep "^>" | cut -f 1 -d ' ' | sed 's/^>//;s/\..*//' \
    | sort |uniq > $PROCDATA/proteinIDs_from_fasta.txt

## get all proteins in gff3
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep ENSP  | sed 's/.*protein_id=//' \
    | sort | uniq  > $PROCDATA/proteinIDs_from_gff3.txt 

## manual check: most proteins in fasta that do not appear in gff3
## are annotated to a scaffold
## diff
diff  $PROCDATA/proteinIDs_from_fasta.txt \
      $PROCDATA/proteinIDs_from_gff3.txt 



## CHROMOSOME TABLE
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tchromosome\t" | cut -f 1,5 | grep "^[0-9]" \
    | sort -n > $MAMDATA/chromosomes/sequenceIndex.tmp
## add X,Y and MT
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tchromosome\t" | cut -f 1,5 \
    | grep "^[X|Y]" >> $MAMDATA/chromosomes/sequenceIndex.tmp
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tchromosome\t" | cut -f 1,5 \
    | grep "^MT" >> $MAMDATA/chromosomes/sequenceIndex.tmp

echo -e "#ID\tname\tlength" > $MAMDATA/chromosomes/sequenceIndex.csv
nl -w2 $MAMDATA/chromosomes/sequenceIndex.tmp \
   >> $MAMDATA/chromosomes/sequenceIndex.csv

\rm $MAMDATA/chromosomes/sequenceIndex.tmp


## FEATURE TABLE: 21k genes + 15k pseudogenes
## NOTE: rRNA etc. are TRANSCRIPTS!
##       type of gene is in description column as biotype
## get genes: everything with \tID=gene
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tID=gene" \
	   > $ORIGDATA/Homo_sapiens.GRCh38.110_genes.gff3
gzip $ORIGDATA/Homo_sapiens.GRCh38.110_genes.gff3

## get transcripts: every mRNA and everything with transcript\t
## and everything with gene_segment
#gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
#    | grep -P "\tmRNA\t" > $ORIGDATA/Homo_sapiens.GRCh38.110_transcripts.gff3
#gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
#    | grep -P "transcript\t" \
#ď	   >> $ORIGDATA/Homo_sapiens.GRCh38.110_transcripts.gff3

## get transcripts: everything with \tID=transcript
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tID=transcript" \
	   > $ORIGDATA/Homo_sapiens.GRCh38.110_transcripts.gff3
gzip $ORIGDATA/Homo_sapiens.GRCh38.110_transcripts.gff3

## get proteins: every CDS
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tCDS\t" \
	   > $ORIGDATA/Homo_sapiens.GRCh38.110_proteins.gff3
gzip $ORIGDATA/Homo_sapiens.GRCh38.110_proteins.gff3

## get transcript: exons and 5/3'UTR
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\texon\t" \
	   > $ORIGDATA/Homo_sapiens.GRCh38.110_exons.gff3
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "prime_UTR\t" \
	   >> $ORIGDATA/Homo_sapiens.GRCh38.110_exons.gff3
gzip $ORIGDATA/Homo_sapiens.GRCh38.110_exons.gff3

## get 5' and 3'UTR
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "prime_UTR\t" \
	   > $ORIGDATA/Homo_sapiens.GRCh38.110_utr.gff3
gzip $ORIGDATA/Homo_sapiens.GRCh38.110_utr.gff3

### GENERATE SOME USEFUL ANNOTATIONS

## protein and transcript lengths
$bioawk -c fastx '{ print $name, length($seq) }' \
       < $ORIGDATA/Homo_sapiens.GRCh38.pep.all.fa.gz \
       > $ORIGDATA/Homo_sapiens.GRCh38.pep.all_lengths.tsv
$bioawk -c fastx '{ print $name, length($seq) }' \
       < $ORIGDATA/Homo_sapiens.GRCh38.cdna.all.fa.gz \
       > $ORIGDATA/Homo_sapiens.GRCh38.cdna.all_lengths.tsv

## a simple protein<->transcript map
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tCDS\t" | sed 's/.*transcript://;s/;protein_id=/\t/' \
    | sort | uniq >  $ORIGDATA/protein_transcript_map.tsv


## ID MAPPINGS

## grep uniprot/ensembl/name protein mappings
gunzip -c $ORIGDATA/HUMAN_9606_idmapping.dat.gz \
    | grep Ensembl_PRO | cut -f 1,3 - \
			     > $ORIGDATA/uniprot_ensembl.dat
gunzip -c $ORIGDATA/HUMAN_9606_idmapping.dat.gz \
    | grep Gene_Name | cut -f 1,3 - > $ORIGDATA/uniprot_name.dat

## get human gene name synonyms
## only human
gunzip -c $ORIGDATA/gene_info.gz | grep -P "^9606\t" \
    | cut -f 3,5  > $ORIGDATA/gene_synonyms.tsv



## ANNOTATE tRNA
## TODO: get tRNA set and integrate into feature table
##http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
R --vanilla < $SRC/annotation_tRNA.R


### GENERATE FEATURE and other ANNOTATION TABLES 

## genome feature (transcripts, proteins) annotation file
R --vanilla < $SRC/annotation.R >  $MAMDATA/log/annotation.log 


### GENOME ANNOTATION FILES FOR CODON MAPPING

## CODING SEQUENCES: generate transcript fasta with only the coding
## regions: reads Homo_sapiens.GRCh38.cdna.all.fa.gz and
## Homo_sapiens.GRCh38.110_utr.gff3.gz, cuts out only the coding
## sequences, and generates processedData/coding.fa and a QC figure in
## processedData/annotation/
R --vanilla < $SRC/map_utr.R >  $MAMDATA/log/utr.log

## MAP CDS per protein
## PROTEIN COORDINATES and SORTED (!) CDS LENGTHS
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.110.gff3.gz \
    | grep -P "\tCDS\t" | cut -f 1,4,5,7,9 | sed 's/\tID.*protein_id=/\t/g'  \
    | sort -k5,5 -k2,2n | awk ' {print $5, $1, $2, $3, $4 }' \
	> $PROCDATA/protein_cds_coordinates.tsv
cat $PROCDATA/protein_cds_coordinates.tsv \
    | awk ' {print $1, $4-$3 +1, $5 } ' \
	  > $PROCDATA/protein_cds_lengths.tsv

## TODO: finish this to write out for mapping protein sites to the genome
## via 2nd codon in transcript -> exon number, and exon coordinates
R --vanilla < $SRC/map_cds.R



## length of coding regions
$bioawk -c fastx '{ print $name, length($seq) }' \
       < $PROCDATA/coding.fa \
       | sed 's/,//g' > $PROCDATA/coding_length.tsv
## length of proteins
$bioawk -c fastx '{ print $name, length($seq) }' \
       < $ORIGDATA/Homo_sapiens.GRCh38.pep.all.fa.gz \
       | sed 's/,//g' > $PROCDATA/protein_length.tsv

## CODON FREQUENCY
R --vanilla < $SRC/codon_frequency.R 



### required for PROTEIN SECONDARY STRUCTURE PRED.

## get a cleaned protein fasta file: remove all proteins <6 AA
gunzip -c $ORIGDATA/Homo_sapiens.GRCh38.pep.all.fa.gz \
    | $bioawk -c fastx '{if(length($seq) > 5){print ">"$name; print $seq}}' \
	      - > $PROCDATA/Homo_sapiens.GRCh38.pep.large.fa

## get a list of all IDs in this fasta file
grep ">" $PROCDATA/Homo_sapiens.GRCh38.pep.large.fa \
    | sed 's/> *//' > $PROCDATA/Homo_sapiens.GRCh38.pep.large_ids.txt



## RUN CALCULATIONS: hmmer/PFAM, S4pred, iupred3
## NOTE: follow the manual steps in the following scripts
if false; then
    $SRC/setup_hpc.sh
    $SRC/run_calculations.sh
fi
