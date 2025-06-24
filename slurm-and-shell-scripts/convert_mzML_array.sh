#!/bin/bash


# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 120:00:00
### Cores sent by owning script

module load mono
module load thermorawfileparser


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

FILE=`sed -n ${SLURM_ARRAY_TASK_ID}p samples.txt`

mono $THERMORAWFILEPARSER_BIN/ThermoRawFileParser.exe -i ${FILE}
mono $THERMORAWFILEPARSER_BIN/ThermoRawFileParser.exe -i ${FILE} -m 1	#Generate metadata stats for all files

