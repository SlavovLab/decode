#!/bin/bash

### RUN VARIOUS PREDICTION ALGORITHM BASED ON HUMAN GENOME
## PROTEIN, RNA and DNA SEQUENCES

## NOTE: PATHs are set in setup.sh

## required tools and paths
bioawk=/home/raim/scripts/bioawk

### S4PRED : PREDICT SECONDARY STRUCTURE FOR ALL ENSEMBL PROTEINS 

## NOTE: run at NEU Discovery HPC platform, and data copied to main
## data/mistrans/processedData from there
if false; then
    ## set up directories at the `discovery` HPC of NEU
    DSRC=<SOURCECODE_PATH_ON_HPC>
    ##$SRC/setup_hpc.sh
    mkdir ~/s4pred
    cd ~/s4pred/
    ## run s4pred secondary structure prediction on all human proteins
    $DSRC/s4pred_batch.sh
    ## check results: <102691 of 121449 after first run
    ls ~/data/mistrans/processedData/s4pred/ -1 |wc -l
    grep ">" ~/data/mistrans/processedData/Homo_sapiens.GRCh38.pep.large.fa|wc -l
    ## repeat CANCELLED jobs - NOTE: run twice with ncpu=40 to catch last 4000
    $DSRC/s4pred_repeat.sh
    ## collect missing in chunks of 1000
    ncpu=20
    seq_per_call=1000
    sbatch --job-name=s4pred_missing --output=s4pred_missing.out -p short \
	   --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB \
	   $DSRC/s4pred_missing.sh $seq_per_call $ncpu
    sbatch --job-name=s4pred_missing2 --output=s4pred_missing2.out -p short \
	   --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB \
	   $DSRC/s4pred_missing.sh $seq_per_call $ncpu
    sbatch --job-name=s4pred_missing3 --output=s4pred_missing3.out -p short \
	   --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB \
	   $DSRC/s4pred_missing.sh $seq_per_call $ncpu

    ## clean up
    rm *.fas chunk.fasta missing_ids.txt -f
fi

## copy s4pred data to exon
rsync -avz r.machne@login.discovery.neu.edu:data/mistrans/processedData/s4pred \
      $PROCDATA/
## COLLECT S4pred predictions, and write output to fasta
$SRC/s4pred_collect.sh
## check completeness
grep ">" $PROCDATA/Homo_sapiens.GRCh38.pep.large_s4pred.fas |wc -l
grep ">" $PROCDATA/Homo_sapiens.GRCh38.pep.large.fa|wc -l
## compress!
gzip $PROCDATA/Homo_sapiens.GRCh38.pep.large_s4pred.fas


### IUPRED3

## iupred3 requires to be run on each sequence separately
## * loop through ID file, check if it's already done,
##   then extract single fasta and run iupred3
## * FAST: 24 h in simple loop on a desktop computer. NOTE:
##   can be run twice as fast when the two loops in script
##   are run manually in two parallel terminals
$SRC/iupred3_run.sh
## NOTE: columns are
# POS   RES     IUPRED2 ANCHOR2



### SCAN FOR PFAM DOMAINS

hmmpress $ORIGDATA/pfam/Pfam-A.hmm
hmmscan --noali --notextw \
	-o $PROCDATA/Homo_sapiens.GRCh38.pep.large_annotations.txt \
	--domtblout $PROCDATA/Homo_sapiens.GRCh38.pep.large_annotations.tbl \
	-E 0.01 --cpu 4 \
	$ORIGDATA/pfam/Pfam-A.hmm \
	$PROCDATA/Homo_sapiens.GRCh38.pep.large.fa \
    &> $PROCDATA/Homo_sapiens.GRCh38.pep.large_annotations.log

## replace spaces with ; and omit descript to make parsing easier
sed 's/ \+/;/g'  \
    $PROCDATA/Homo_sapiens.GRCh38.pep.large_annotations.tbl \
    | cut -d';' -f 1-22 \
	  > $PROCDATA/Homo_sapiens.GRCh38.pep.large_annotations.csv



