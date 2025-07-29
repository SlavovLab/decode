
## Record of setting up data and tools for high performance computating
## at Northeastern's discovery cluster

mkdir ~/programs
mkdir data
mkdir data/mistrans
mkdir work

## copy required data: originalData and processedData (run from source
## computer; see setup.sh for folders)
## rsync -avz ${MYDATA}/mistrans_data/processedData r.machne@login.discovery.neu.edu:data/mistrans/
## rsync -avz ${MYDATA}/mistrans_data/originalData r.machne@login.discovery.neu.edu:data/mistrans/


## clone THIS code base
cd ~/work
git clone git@github.com:SlavovLab/decode.git


## DOWNLOAD AND COMPILE REQUIRED TOOLS

## bioawk
cd ~/programs
git clone https://github.com/lh3/bioawk 
cd bioawk
make

## S4pred: protein secondary structure
cd ~/programs
git clone https://github.com/psipred/s4pred
cd s4pred
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights.tar.gz
tar -xvzf weights.tar.gz

## HMMER: PFAM predictions
cd ~/programs
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxvf hmmer.tar.gz
cd hmmer-3.4

## compile hmmer
srun -N 1 -n 28 --constraint=ib --pty /bin/bash
module load gcc/9.2.0
./configure --prefix /home/r.machne/
make
make check
make install

## get Pfam
srun -N 1 -n 28 --constraint=ib --pty /bin/bash
mkdir ~/data/mistrans/originalData/pfam
cd ~/data/mistransoriginalData/pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/userman.txt
echo "Pfam-A.hmm.gz downloaded on 2024-04-11 as the current release (36.0) from EBI" > README.md
gunzip Pfam-A.hmm.gz
hmmpress ~/data/mistrans/originalData/pfam/Pfam-A.hmm
