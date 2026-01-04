#!/bin/bash

#SBATCH --job-name=$3
#SBATCH --output=PGDB_$3.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=50Gb
#SBATCH --time=24:00:00
#SBATCH --partition=short

module load samtools/1.10

'''
This script takes in a raw RNA-sequening data file in .bam format that was downloaded from GDC and aligns it, and saves the read depth at each base location.

Input:
1 = File ID (GDC, string)
2 = File name (GDC, string)
3 = Sample name (GDC-CaseID_Tumor or GDC-CaseID_Normal, string)

Output:
Table with Chromosome, Position, and N reads columns.

Usage:
- run on HPC cluster (Discovery) with Slurm. Checkpoints before each step to determine if step was previously completed and target file exists. If not, step is run. Can run for a single file or submit as an array for all files in a dataset.
- example command to submit single job to slurm: sbatch PGDB_N_Reads.sh 00904565-9953-412f-9d55-4ff719afbc5b a7f02909-cfb7-4d53-9773-9d26ef365609.rna_seq.transcriptome.gdc_realn.bam C3L-02968_Normal
- example command to submit multiple jobs using sample sheet info (see pipeline instructions part 1.2): for i in {0..185}; do  sbatch --job-name=${sampleArray[$i]} --output=PGDB_${sampleArray[$i]}.out PGDB_N_reads.sh ${idArray[$i]} ${fileArray[$i]} ${sampleArray[$i]}; done

'''


echo $1 #FileID. File ID in GDC.
echo $2 #FileName. File name in GDC.
echo $3 #SampleName. Desired file name. All new files will have sampleName as head of file name

# change to directory where RNAseq data is stored
cd /scratch/tsour.s/Gillette_LUAD/RNAseq_data/$1 # 

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
    
# trim adapters, read quality filtering, make QC outputs [adapted from Spritz]
if [[ ! -f $3.paired.trim_1.fq ]]
then 
	srun /home/tsour.s/fastp -q 28 -i $3.paired_1.fq -I $3.paired_2.fq -o $3.paired.trim_1.fq -O $3.paired.trim_2.fq -h $3.fastp_QC_report.html -j $3.fastp_QC_report.json -w 16 -R $3.fastp_report --detect_adapter_for_pe
	echo "adapter trimming complete"
fi

#align using hisat2
if [[ ! -f $3.hisat2_summary.txt ]]
then 
	srun /home/tsour.s/hisat2-2.2.1/hisat2 -x grch38_110_hisat2_index -1 $3.paired.trim_1.fq -2 $3.paired.trim_2.fq --dta -S $3.hisat2.sam --summary-file $3.hisat2_summary.txt -p 28
	echo 'hisat2 alignment complete'
fi
  
#convert to bam
if [[ ! -f $3.hisat2.bam ]]
then 
	srun samtools view -bS $3.hisat2.sam > $3.hisat2.bam
	echo 'bam conversion complete'
fi

#sort
if [[ ! -f $3.sorted.bam ]]
then 
	srun samtools sort $3.hisat2.bam -o $3.sorted.bam
	echo 'sorted bam'
fi
  
#create index
if [[ ! -f $3.sorted.bai ]]
then 
	srun samtools index -b $3.sorted.bam $3.sorted.bai
	echo 'index created'
fi

#get read depth at each base pair 
if [[ ! -f $3.per_base_coverage.txt ]]
then 
	run samtools depth -a -H $3.sorted.bam -o $3.per_base_coverage.txt
	echo 'depth file created'
fi

