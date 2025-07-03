#!/bin/bash -l

# usage 
# sbatch scriptname.sh python-qc-dir file_prefix

# !!!!!!! These modules have to be load once ever per user on the head or login node !!!!!!!!
# module load perl
# module load gcc
# module load make
# cpanm install XML::Parser

# Define job name
#SBATCH -J PXD_qc
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 120:00:00
# Define array length - not used here
# Define cores
#SBATCH -c 10

module load gcc
module load perl
module load python

# List all modules
module list
#
echo =========================================================
echo SLURM job: submitted date = $(date)
date_start=$(date +%s)
hostname
echo Current directory: $(pwd)
echo "Print the following environmental variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"

#############################################


#### These are post-processing scripts that I have written (Andy), here I am not doing an array job, so I can process a full run in one go
#FILE_PREFIX=${PWD##*/}
pathtopythonqc=$1
FILE_PREFIX=$2

python3 ${pathtopythonqc}/Convert_pepXML_toCSV.py interact-ipro.pep.xml ${FILE_PREFIX}
python3 ${pathtopythonqc}/CalculateFDR_and_threshold_v2.py -f ${FILE_PREFIX}_interact-ipro.pep.tsv  -d rev_ -t 0.01
python3 ${pathtopythonqc}/calculate_psm_stats.py ${FILE_PREFIX}_interact-ipro.pep_thresholded.tsv  ${FILE_PREFIX}

