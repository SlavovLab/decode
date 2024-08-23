#!/bin/bash
#SBATCH --job-name=mouse_lung_DP
#SBATCH --output=mouse_lung_DP.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100Gb
#SBATCH --time=120:0:0
#SBATCH --partition=long
source /home/tsour.s/.bashrc
srun mono $MQ_1_6_17_0 /home/tsour.s/MQ_templates/mouse_lung_DP.xml
