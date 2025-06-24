#!/bin/bash -l

# usage 
# sbatch run-combine-and-threshold-BATCH.sh pxdlist.txt outsuffix

# !!!!!!! These modules have to be load once ever per user on the head or login node !!!!!!!!
# module load perl
# module load gcc
# module load make
# cpanm install XML::Parser

# Define job name
#SBATCH -J PXD_threshold_combine
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 120:00:00
# Define array length - not used here
# Define cores
#SBATCH -c 5

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
pxdlist=$1
outsuffix=$2

#### combine and run things

python3 1-python-qc-scripts/Combine-PXD-thresholded-tsv.py -l ${pxdlist} -o ${outsuffix}
python3 1-python-qc-scripts/Collapse_by_max_prob_v2.py ${outsuffix}_PepAtlas-COMBINED-thresholded.tsv
python3 1-python-qc-scripts/CalculateFDR_and_threshold_untested_v2.py -f ${outsuffix}_PepAtlas-COMBINED-thresholded_peptide_final_collapsed.tsv -d rev_ -t 0.01 -p


