#!/bin/bash


cd $MAMDATA/originalData/

## required utils (besides bash standards)
pdftotext=/usr/bin/pdftotext
myBWTBG=~/programs/ucsc_utils/bigWigToBedGraph

## ENSEMBL GENOME ANNOTATION FILE
#rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_embl/homo_sapiens .
## NOTE: downloaded on 2023-11-27
echo "downloading genome gff3 file"
rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz  .


echo "downloading chromosome fasta"
# # NOTE: downloaded on 20240807
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz .

## TODO: compare with chromosome lengths in data/mammary/chromosomes/sequenceIndex.csv
gunzip -c Homo_sapiens.GRCh38.dna.toplevel.fa.gz |grep ">"

## download ensembl human proteins
echo "downloading protein fasta file"
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz .

echo "downloading protein coding transcript (cDNA) fasta file"
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz .

echo "downloading noncoding transcript fasta file"
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz .

echo "downloading tRNAs from tRNAscan-SE"
wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
tar zxvf hg38-tRNAs.tar.gz hg38-tRNAs.bed


echo "downloading various data sets"

## liftOver files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -P $MAMDATA/originalData/

## tRNA DATA

## @Gingold2014: tRNAs differentially upregulated in growing and differentiating
## cells; NOTE: auto-download forbidden; stored manually
wget https://www.cell.com/cms/10.1016/j.cell.2014.08.011/attachment/67879827-afd7-4c7c-924e-f20ec5c4bf69/mmc2.xls -P $MAMDATA/originalData/ -O gingold14_mmc2.xls
## @Torres2019: tRNA classification
## NOTE: auto-download forbidden; stored manually
wget https://www.pnas.org/doi/suppl/10.1073/pnas.1821120116/suppl_file/pnas.1821120116.sd03.xlsx -P $MAMDATA/originalData/ -O torres19_s3.xlsx

### describePROT data
## NOTE: These files were extracted from the json file of describePROT by
## Andrew Leduc, 20240517.

## ID MAPPING

## TODO: compare uniprot mapping
## UNIPROT ID MAPPING
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz

## NCBI GENE INFO
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz

## ENSEMBL BIOMART: uniprot/swissprot, GO, slimGO for ensembl peptide ID
wget -O ensembl_gene_uniprot_$(date +"%Y%m%d").tsv 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_peptide_id" />
		<Attribute name = "uniprot_gn_id" />
		<Attribute name = "uniprot_isoform" />
		<Attribute name = "uniprotswissprot" />
	</Dataset>
</Query>'
gzip ensembl_gene_uniprot_$(date +"%Y%m%d").tsv

now=$(date +"%Y%m%d");wget -O ensembl_gene_goslim_$(date +"%Y%m%d").tsv 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_peptide_id" />
		<Attribute name = "goslim_goa_accession" />
		<Attribute name = "goslim_goa_description" />
	</Dataset>
</Query>'
gzip ensembl_gene_goslim_$(date +"%Y%m%d").tsv
now=$(date +"%Y%m%d");wget -O ensembl_gene_go_$(date +"%Y%m%d").tsv 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_peptide_id" />
        	<Attribute name = "go_id" />		    
		<Attribute name = "go_linkage_type" />
	</Dataset>
</Query>'
gzip ensembl_gene_go_$(date +"%Y%m%d").tsv
## TODO: also save a obo file to get term defintions etc.

## ensembl/refseq protein mapping
now=$(date +"%Y%m%d");wget -O ensembl_refseq_$(date +"%Y%m%d").tsv 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_peptide_id" />
		<Attribute name = "refseq_protein" />
	</Dataset>
</Query>'
gzip ensembl_refseq_$(date +"%Y%m%d").tsv

## ensembl/refseq protein mapping
now=$(date +"%Y%m%d");wget -O ensembl_refseq_transcripts_$(date +"%Y%m%d").tsv 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "refseq_mrna" />
	</Dataset>
</Query>'
gzip ensembl_refseq_transcripts_$(date +"%Y%m%d").tsv

## Guttierez-Monreal et al. 2016: circadian expression in MCF10A
## in geo soft format
wget https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21250/soft/GPL21250_family.soft.gz

## Janes et al. 2011: stochastic profiling of MCF10A organoid
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120030/soft/GSE120030_family.soft.gz

## human GO annotation
wget http://geneontology.org/gene-associations/goa_human.gaf.gz
## goslim generic: recommended for humans
wget https://current.geneontology.org/ontology/subsets/goslim_generic.obo



## PHASTCONS TRACK
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw
## extract regions to bed
## TODO: use this to extract data from all exons
##~/programs/ucsc/bigWigToBedGraph $MAMDATA/originalData/hg38.phastCons7way.bw -chrom=chr1 -start=1200 -end=2000  out.bed

## EBI pfam annotation of human proteins
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/proteomes/9606.tsv.gz

## PFAM 
mkdir ./pfam
cd ./pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
cd -

## Hu.MAP 2.0 protein complexes
wget http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt

## CORUM
## NOTE: manually downloaded since wget doesn't work
wget -O humanComplexes_corum_v4.1.txt.zip https://mips.helmholtz-muenchen.de/corum/download/releases/current/humanComplexes.txt.zip -
#wget https://mips.helmholtz-muenchen.de/corum/download/releases/current/uniprot_corum_mapping.txt
unzip humanComplexes_corum_v4.1.txt.zip

## UNIPROT ALPHA FOLD
mkdir ./afold
cd ./afold
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/v4/UP000005640_9606_HUMAN_v4.tar
cd -

## PROTEIN HALF-LIVES - @Mathieson2018
## all proteins for different cell types
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM5_ESM.xlsx
## proteasome
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM6_ESM.xlsx
## nuclear pore
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM7_ESM.xlsx

## STRUCTURE PREDICTIONS at DescribePROT
wget http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/download_database/9606_database.json # raw data 1.5 GB
wget http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/download_database_value/9606_value.csv # protein level summaries

## in vitro S20 proteasome targets
wget https://www.embopress.org/doi/suppl/10.1038/s44320-024-00015-y/suppl_file/44320_2024_15_moesm1_esm.xlsx

## @Watson2023 proteome and phospoproteome of temperature and osmotic
## shock or adaptation - NOTE: mouse cells - TODO: map mouse orthologs
## and/or start mouse genome browser project.
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06626-z/MediaObjects/41586_2023_6626_MOESM4_ESM.xlsx # proteome
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06626-z/MediaObjects/41586_2023_6626_MOESM5_ESM.xlsx # phosphoproteome

### PDB FILE PROTEASOME
wget https://files.rcsb.org/download/6KWY.pdb
##wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/6kwy.xml.gz

## SIFTS : PDB gene ID mapping
##wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_ensembl.csv
##wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_pfam.csv
wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_ensembl.csv.gz
wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_pfam.csv.gz

## @Yang2022: thermal stability prediction
## https://structure-next.med.lu.se/ProTstab2/
wget https://structure-next.med.lu.se/ProTstab2/download/download_human_predict_result -P $MAMDATA/originalData/ -O ProTstab2_human.csv

## @Savitski2014: Tracking cancer drugs in living cells by thermal
## profiling of the proteome
## NOTE: forbidden, requires manual download
wget https://www.science.org/doi/suppl/10.1126/science.1255784/suppl_file/table_s11_thermal_profiling_jurkat_intact_cells.xlsx \
     -O $MAMDATA/originalData/savitski14_tableS11.xlsx
wget https://www.science.org/doi/suppl/10.1126/science.1255784/suppl_file/table_s3_thermal_profiling_atp_cell_extract.xlsx \
     -O $MAMDATA/originalData/savitski14_tableS3.xlsx

### RNA MOD DATA

## @Zhang2023: pseudouridylation
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41589-023-01304-7/MediaObjects/41589_2023_1304_MOESM3_ESM.xlsx \
     -O $MAMDATA/originalData/zhang23_stables_1-6.xlsx

## @Dai2023: pseudouridylation
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01505-w/MediaObjects/41587_2022_1505_MOESM3_ESM.xlsx \
     -O $MAMDATA/originalData/dai23_stables_1-23.xlsx

## m6A-Atlas modification database - dl 20240811
wget http://rnamd.org/m6a/download/high/hg38/hg38_CL_Tech.txt.gz \
     -P $MAMDATA/originalData/ 
wget http://rnamd.org/m6a/download/high/hg38/human_all_logFC.txt.gz \
     -P $MAMDATA/originalData/ 

## PSI-SITE collection
wget http://180.208.58.19/%ce%a8-WHISTLE/Download/Human_Experiment.txt \
     -O $MAMDATA/originalData/piano_psi_human.txt

### RIBOSOME PROFILING DATA

## Ribo-seq data resources
## https://riboseq.org/
## * download bigwig format of riboseq data from:
##   https://gwips.ucc.ie/downloads/index.html
## * Taiwanese database: mapped to refseq transcripts:
## http://cosbi4.ee.ncku.edu.tw/HRPDviewer/user/download_p1

## UCSC Ribo-seq track from GWIPS-viz

## For each study used to generate this track, raw fastq files were
## downloaded from a repository (e.g., NCBI GEO datasets). Cutadapt was
## used to trim the relevant adapter sequence from the reads, after which
## reads below 25 nt in length were discarded. The trimmed reads were
## aligned to ribosomal RNA using Bowtie and aligning reads were
## discarded. The remaining reads were then aligned to the hg38 (GRCh38)
## genome assembly using Bowtie. An offset of 15 nt (to infer the
## position of the A-site) was added to the most 5' nucleotide coordinate
## of each uniquely-mapped read.

wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/gwipsvizRiboseq/gwipsvizRiboseq.bw \
     -P $MAMDATA/originalData/ 
$myBWTBG $MAMDATA/originalData/gwipsvizRiboseq.bw $MAMDATA/processedData/gwipsvizRiboseq.bed
## save disk space
gzip  $MAMDATA/processedData/gwipsvizRiboseq.bed
\rm $MAMDATA/originalData/gwipsvizRiboseq.bw

### NUCLEOSOME DATA
## https://compbio-zhanglab.org/NUCOME/


## AMINO ACID MASSES for MS
wget https://proteomicsresource.washington.edu/docs/amino_acid_masses.xlsx  \
     -P $MAMDATA/originalData/ 


## Tissue Proteomes

## NOTE: access via url forbidden, MANUAL DOWNLOAD
##wget https://www.embopress.org/doi/suppl/10.15252/msb.20188503/suppl_file/msb188503-sup-0003-tableev1.zip \
##     -O $MAMDATA/originalData/wang19_table_EV1.zip 
unzip $MAMDATA/originalData/wang19_table_EV1.zip -d $MAMDATA/originalData/
mv $MAMDATA/originalData/Table_EV1.xlsx $MAMDATA/originalData/wang19_table_EV1.xlsx

## Tissue Transcriptomes

## NOTE: access via url forbidden, MANUAL DOWNLOAD
##wget https://www.embopress.org/doi/suppl/10.15252/msb.20188513/suppl_file/msb188513-sup-0004-tableev2.zip -O $MAMDATA/originalData/eraslan19_table_EV2.zip 
## mv ~/Downloads/msb188513-sup-0004-tableev2.zip ~/data/mammary/originalData/eraslan19_table_EV2.zip
unzip $MAMDATA/originalData/eraslan19_table_EV2.zip -d $MAMDATA/originalData/
mv $MAMDATA/originalData/Table_EV2/Table_EV2.tsv $MAMDATA/originalData/eraslan19_table_EV2.tsv
mv $MAMDATA/originalData/Table_EV2/README.txt $MAMDATA/originalData/eraslan19_table_EV2_readme.txt
rm $MAMDATA/originalData/Table_EV2/ -rf


## MSigDB Hallmark Gene sets
## NOTE: no access via wget, downloaded manually from
##
## wget https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt \
##      -P $MAMDATA/originalData/ 
