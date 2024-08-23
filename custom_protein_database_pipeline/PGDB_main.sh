#!/bin/bash

#SBATCH --job-name=${3}
#SBATCH --output=PGDB_${3}.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80Gb
#SBATCH --time=24:00:00
#SBATCH --partition=short

module load freebayes/1.3.4
module load gffcompare/0.11.5
module load stringtie/2.1.1
module load samtools/1.10

'''
This script takes in a raw RNA-sequening data file in .bam format that was downloaded from GDC and in-silico translates it into a sample-specific protein database using a two-pronged appproach: SNV/InDel calling and de-novo transcriptome assembly.

Input:
1 = File ID (GDC, string)
2 = File name (GDC, string)
3 = Sample name (GDC-CaseID_Tumor or GDC-CaseID_Normal, string)

Output:
Custom fasta file of the in-silico translated coding region transcripts to be used for proteomics data searches

Usage:
- run on HPC cluster (Discovery) with Slurm. Checkpoints before each step to determine if step was previously completed and target file exists. If not, step is run. Can run for a single file or submit as an array for all files in a dataset.
- example command to submit single job to slurm: sbatch PGDB_Jul2022.sh 00904565-9953-412f-9d55-4ff719afbc5b a7f02909-cfb7-4d53-9773-9d26ef365609.rna_seq.transcriptome.gdc_realn.bam C3L-02968_Normal
- example command to submit multiple jobs using sample sheet info (see pipeline instructions part 1.2): for i in {0..185}; do  sbatch --job-name=${sampleArray[$i]} --output=PGDB_${sampleArray[$i]}.out PGDB_Jul2022.sh ${idArray[$i]} ${fileArray[$i]} ${sampleArray[$i]}; done

'''


echo $1 #FileID. File ID in GDC.
echo $2 #FileName. File name in GDC.
echo $3 #SampleName. Desired file name. All new files will have sampleName as head of file name

# change to directory where RNAseq data is stored
cd /scratch/tsour.s/Clark_ccRCC/RNAseq_data/$1

''' start here with 1 .bam file '''
#sort gdc bam file by read name
if [[ ! -f $3.collate.bam ]]
then 
	srun samtools collate -o $3.collate.bam $2
	echo 'sorted bam'
fi 
 
#write paired end reads to 2 separate files
if [[ ! -f $3.paired_1.fq ]]
then 
	srun samtools fastq -1 $3.paired_1.fq -2 $3.paired_2.fq -n $3.collate.bam
	echo 'fasta conversion complete'
fi

''' start here with 2 paired end .fq files '''
# trim adapters, read quality filtering, make QC outputs [adapted from Spritz]
if [[ $3.paired.trim_1.fq ]]
then 
	srun /home/tsour.s/fastp -q 28 -i $3.paired_1.fq -I $3.paired_2.fq -o $3.paired.trim_1.fq -O $3.paired.trim_2.fq -h $3.fastp_QC_report.html -j $3.fastp_QC_report.json -w 16 -R $3.fastp_report --detect_adapter_for_pe
	echo "adapter trimming complete"
fi

#align using hisat2
if [[ ! -f $3.hisat2_summary.txt ]]
then 
	srun /home/tsour.s/hisat2-2.2.1/hisat2 -x /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic -1 $3.paired.trim_1.fq -2 $3.paired.trim_2.fq --dta -S $3.hisat2.sam --summary-file $3.hisat2_summary.txt -p 28
	echo 'hisat2 alignment complete'
fi
  
#convert to bam
if [[ ! -f $3.hisat2.header.bam ]]
then 
	srun samtools view -bS $3.hisat2.sam > $3.hisat2.bam
	srun samtools view -bS -h $3.hisat2.sam > $3.hisat2.header.bam #with header for customProDB
	echo 'bam conversion complete'
fi
   
#sort
if [[ ! -f $3.sorted.header.bam ]]
then 
	srun samtools sort $3.hisat2.bam -o $3.sorted.bam
	srun samtools sort $3.hisat2.header.bam -o $3.sorted.header.bam
	echo 'sorted bam'
fi
  
#create index
if [[ ! -f $3.sorted.header.bai ]]
then 
	srun samtools index -b $3.sorted.bam $3.sorted.bai
	srun samtools index -b $3.sorted.header.bam $3.sorted.header.bai
	echo 'index created'
fi

#assemble with string tie
if [[ ! -f $3.stringtie.gtf ]]
then 
	srun stringtie $3.sorted.bam -o $3.stringtie.gtf -G /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.100.gff3 -f 0.15 -c 2 -p 28
	echo 'stringtie assembly complete'
fi

#compare and annotate assembly with reference gff
if [[ ! -f $3.gffcmp ]]
then 
	srun gffcompare -r /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.100.gff3 -o $3.gffcmp $3.stringtie.gtf
	echo 'gffcompare complete'
fi
   
#filter for CDS (exons) only, and filter for all else and write these 2 transcript fastas to file
if [[ ! -f $3.gffcmp.annotated.transcript.gtf ]]
then 
	srun /home/tsour.s/gffread/gffread $3.gffcmp.annotated.gtf --no-pseudo --force-exons -M -T -o $3.gffcmp.annotated.transcript.CDS.gtf
	srun /home/tsour.s/gffread/gffread -w $3.CDS.DNA.fasta -g /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa $3.gffcmp.annotated.transcript.CDS.gtf 
	srun /home/tsour.s/gffread/gffread $3.gffcmp.annotated.gtf --no-pseudo -M -T -o $3.gffcmp.annotated.transcript.gtf
	srun /home/tsour.s/gffread/gffread -w $3.all.DNA.fasta -g /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa $3.gffcmp.annotated.transcript.gtf
	echo 'genomic fastas written'
fi

#convert off compare output gtf to bed format, high confident and all but unknown/opposite strand/repeat
if [[ ! -f $3.gffcmp.annotated.all.bed  ]]
then 
	srun python /home/tsour.s/scripts_templates/gffcompare_to_bed.py $3.gffcmp.annotated.gtf $3.gffcmp.annotated.HC.bed -C "=,c,k,m,n,j,e"
	srun python /home/tsour.s/scripts_templates/gffcompare_to_bed.py $3.gffcmp.annotated.gtf $3.gffcmp.annotated.all.bed 
	echo 'bed conversion complete'
fi
   
#translate bed
if [[ ! -f $3.translation.fa ]]
then 
	srun python /home/tsour.s/scripts_templates/translate_bed.py $3.gffcmp.annotated.HC.bed --twobit=/home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.2bit --cds --min_length 10 --reference GRCh38 --bed $3.CDS.translation.bed --fasta $3.CDS.translation.fa -v
	srun python /home/tsour.s/scripts_templates/translate_bed.py $3.gffcmp.annotated.all.bed --twobit=/home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.2bit --min_length 10 --reference GRCh38 --bed $3.translation.bed --fasta $3.translation.fa -v
	echo 'translation complete'
fi

#freeBayes variant calling
if [[ ! -f $3.vcf ]]
then 
	srun freebayes -b $3.sorted.bam -f /home/tsour.s/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa  -v $3.vcf
	echo 'variant calling complete'
fi
  
#customProDB R script
if [[ ! -f $3.indel.fasta ]]
then 
	srun Rscript /home/tsour.s/scripts_templates/customProDB_ensembl_discovery.R $3.sorted.header.bam $3.vcf /home/tsour.s/genomes/customProDB_files/ensembl/ $3.
	echo 'customProDB complete'
   
fi

#merge fastas
if [[ ! -f $3.CDS.final.fasta ]]
then 
	srun python /home/tsour.s/scripts_templates/fasta_merge_files_and_filter_unique_sequences.py $3.CDS.final.fasta sequence '^>([^ ]+).*$' $3.indel.fasta $3.rpkm.fasta $3.SNV.fasta $3.CDS.translation.fa
	echo 'fasta merge complete'
fi
