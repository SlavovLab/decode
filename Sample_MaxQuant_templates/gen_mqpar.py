#!/usr/bin/env python3
#coding: utf-8

import argparse
import os
import sys
import re
import time

"""
This script takes in a MaxQuant parameter file template (xml format) and customizes it for a new MaxQuant search.
Input: template_paramters.xml, directory with raw files, desired output file names, MaxQuant version, fasta database and number of threads.
Output: new parameters.xml, directory for storing MaxQuant output files, slurm script for running MaxQuant with the new parameters.xml. Slurm parameters such as cpu, partition, etc. should be updated in the code below.
"""

parser = argparse.ArgumentParser(description='MaxQuant fun')
parser.add_argument('input_xml', type=argparse.FileType('r'))
parser.add_argument('raw_file_folders', type=str, nargs='+')
parser.add_argument('-o', '--outfile', type=str, required=True)
parser.add_argument('-O', '--outdir', type=str, required=True)
parser.add_argument('-f', '--fasta', type=str, required=True)
parser.add_argument('-mq', '--mq-version', type=str, default='1_6_17_0', choices=['1_5_5_1', '1_6_2_3', '1_6_2_10', '1_6_3_2', '1_6_3_3', '1_6_3_4', '1_6_7_0', '1_6_17_0'], help='MaxQuant version')
parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to specify in MaxQuant configuration file')

args = parser.parse_args()

mqpar = open(args.input_xml.name, 'r')
mqpar_text = mqpar.read()
mqpar.close()


# replace fasta path
fasta_path = args.fasta
fasta_path = ('<fastaFilePath>' + fasta_path + '</fastaFilePath>')
mqpar_text = re.sub(r'\<fastaFilePath\>(.|\n|\r)*\<\/fastaFilePath\>', fasta_path, mqpar_text)


file_counter = 0
file_path_repl_text = '<filePaths>\n'

for folder in args.raw_file_folders:
  files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
  # only take raw files
  #files = [f for f in files if f[-5:] == '.wiff']
  files = [f for f in files if f[-4:] == '.raw']

  for file in files:
    file_path_repl_text += ('<string>' + os.path.join(folder, file) + '</string>\n')
    file_counter += 1

file_path_repl_text += '</filePaths>'


mqpar_text = re.sub(r'\<filePaths\>(.|\n|\r)*\<\/filePaths\>', file_path_repl_text, mqpar_text)

# change experiments, fractions, ptms, paramGroupIndices to reflect # of raw files added
experiments_text = '<experiments>\n'
fractions_text = '<fractions>\n'
ptms_text = '<ptms>\n'
group_inds_text = '<paramGroupIndices>\n'
reference_channels_text = '<referenceChannel>\n'

for i in range(0, file_counter):
  experiments_text += '<string>1</string>\n'
  fractions_text += '<short>'+str(i+1)+'</short>\n'
  ptms_text += '<boolean>False</boolean>\n'
  group_inds_text += '<int>0</int>\n'
  reference_channels_text += '<string></string>\n'

experiments_text += '</experiments>\n'
fractions_text += '</fractions>\n'
ptms_text += '</ptms>\n'
group_inds_text += '</paramGroupIndices>\n'
reference_channels_text += '</referenceChannel>'

mqpar_text = re.sub(r'\<experiments\>(.|\n|\r)*\<\/referenceChannel\>', \
  experiments_text+fractions_text+ptms_text+group_inds_text+reference_channels_text, mqpar_text)


# name the output folder after the named xml output
output_folder = os.path.basename(args.outfile)
# remove the .xml, if it exists
output_folder = re.sub(r'\.xml', '', output_folder)
# remove the beginning "mqpar_", if it exists
output_folder = re.sub(r'mqpar_', '', output_folder)
# append the scratch folder
output_folder = (args.outdir + output_folder)

# create the folder
if not os.path.exists(output_folder):
  os.makedirs(output_folder)

# add xml around
mqpar_text = re.sub(r'\<fixedCombinedFolder\>(.|\n|\r)*\<\/fixedCombinedFolder\>',\
  ('<fixedCombinedFolder>' + output_folder + '</fixedCombinedFolder>'), mqpar_text)

# replace andromeda index -- folder inside of the output folder
andromeda_index_path = os.path.join(output_folder, 'andromeda_index')
andromeda_index_path = ('<fixedSearchFolder>' + andromeda_index_path + '</fixedSearchFolder>')
mqpar_text = re.sub(r'\<fixedSearchFolder\>(.|\n|\r)*\<\/fixedSearchFolder\>', andromeda_index_path, mqpar_text)

# replace temp folder -- folder inside of the output folder

temp_path = os.path.join(output_folder, 'temp')
temp_path = ('<tempFolder>' + temp_path + '</tempFolder>')
mqpar_text = re.sub(r'\<tempFolder\>(.|\n|\r)*\<\/tempFolder\>', temp_path, mqpar_text)

# replace number of threads
threads_tag = ('<numThreads>' + str(args.threads) + '</numThreads>')
mqpar_text = re.sub(r'\<numThreads\>(.|\n|\r)*\<\/numThreads\>', threads_tag, mqpar_text)


### create the slurm script
"NOTE: this should be updated to reflect desired memory, cpu, partition, etc."
slurm_script = ('#!/bin/bash\n'
'#SBATCH --job-name={JOBNAME}\n'
'#SBATCH --output={JOBNAME}.out\n'
'#SBATCH --ntasks=1\n'
'#SBATCH --cpus-per-task=24\n'
'#SBATCH --mem=100Gb\n'
'#SBATCH --time=72:0:0\n'
'#SBATCH --partition=long\n'
'source /home/tsour.s/.bashrc\n'
'srun mono ${MQ_VERSION} {MQPAR}\n'
)

slurm_script = re.sub(r'{MQ_VERSION}', ('MQ_' + args.mq_version), slurm_script)
slurm_script = re.sub(r'{MQPAR}', os.path.abspath(args.outfile), slurm_script)
slurm_script = re.sub(r'{JOBNAME}', os.path.basename(output_folder), slurm_script)

# write slurm script - same format as the output folder
slurm_script_path = os.path.join(os.path.expanduser('~'), 'scripts_templates', os.path.basename(output_folder)) + '.sh'
slurm_script_file = open(slurm_script_path, 'w')
slurm_script_file.write(slurm_script)
slurm_script_file.close()

# write XML file
out_file = open(args.outfile, 'w')
out_file.write(mqpar_text)
out_file.close()
print('success!')

