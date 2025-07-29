#!/bin/bash


module load python/3.6.6

ncpu=20
batch=/home/${USER}/work/mistrans/
fasta=/home/${USER}/data/mammary/processedData/Homo_sapiens.GRCh38.pep.large.fa

## 10min for 100 sequences in tests with $ncpu=20
## 600 per hour -> do 600 per run to split into ca. 200 jobs

total=`grep ">" $fasta |wc -l`
seq_per_call=2000
start=0
end=`echo $total $seq_per_call | awk '{print int($1/$2)}' -`

for (( i=$start; i<=$end; i++ )); do

    echo running loop $i
    
    ## run script
    sbatch -a $i --job-name=s4pred_${i} --output=s4pred_${i}.out -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB $batch/scripts/s4pred_run.sh $i $seq_per_call $ncpu
done

