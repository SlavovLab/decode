
# GENOME DOWNLOAD, PROCESSING and ANNOTATION FOR MISTRANSLATION MAPPING and ANALYSIS

by Rainer Machné, 2023-2025.

The scripts in this folder download human genome data from ENSEMBL and
various other sources (supplemental material of publications), and
provide various mapping files, such as exon structures of protein
transcripts, for mapping amino acid substitution sites to the
positions of their codons in ENSEMBL transcripts and chromosomes. The
output of these scripts is used as an input for analysis of amino acid
substitutions.

The scripts are written in `bash` and `R` and depend on various
additionally installed tools. The bash scripts were never actually run
as one script, they are used as a documentation of downloading and
data processing steps run manually in a bash console. 

**All files generated here and required for downstream analysis are
also provided in a separte folder. See README file in
`../decode_analysis`.**

Note, that these scripts were initially hosted in
https://gitlab.com/raim/genomeBrowser/-/tree/master/data/mammary and
copied/adapted from there on 2025-07-28. 

1. `setup.sh` is the master script that sets up a directory structure,
   and calls all other scripts. It further contains a series of 
   processing steps that extracts information from downloaded genome files.
   The steps must be reproduced in the order of that script.
2. `download.sh` contains `wget` calls to download original
   data. *NOTE** that several download commands (`wget`) will not work
   since some download URLs do not support scripted downloads. These
   files must be looked up and downloaded manually.
3. R scripts: these take the pre-processed data files (`setup.sh`) as input
   to do some more involved mapping and coordinate calculations.
4. Bash scripts: additional bash scripts `setup_hpc.sh` and
   `run_calculations.sh` also require manual steps and guide the
   calculation of protein structure predictions, partially on
   Northeastern's Discovery HPC platform.
