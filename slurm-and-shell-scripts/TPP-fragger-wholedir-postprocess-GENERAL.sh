#!/bin/bash -l

# usage sbatch scriptname.sh

# !!!!!!! These modules have to be load once ever per user on the head or login node !!!!!!!!
# module load perl
# module load gcc
# module load make
# cpanm install XML::Parser

# Define job name
#SBATCH -J PXD_fragger
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 120:00:00
# Define array length - not used here
# Define cores
#SBATCH -c 80

module load mono
module load thermorawfileparser
module load gcc
module load perl
module load tpp
module load adoptopenjdk
module load msfragger
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

#Convert whole folder of raw to mzML
# mono $THERMORAWFILEPARSER_BIN/ThermoRawFileParser.exe -d ./
#Make metadata files
# mono $THERMORAWFILEPARSER_BIN/ThermoRawFileParser.exe -d ./ -m 1	

#There are two ways to run MS-Fragger with or without array jobs, see notes elsewhere for pros and cons - here I am used non-array job to run on all mzML files
paramsfile=$1
java -jar $MSFRAGGER_BIN/MSFragger-3.5.jar ${paramsfile} *.mzML

# run TPP prophet commands using all comet output .pep.xml files
#### !!!! check correct DECOY prefix !!!!! If using FragPipe GUI to make database, then default is rev_
xinteract -OPAEd -drev_ -PPM -p0 -Ninteract-prob.pep.xml *pepXML -THREADS=$SLURM_CPUS_ON_NODE
InterProphetParser THREADS=$SLURM_CPUS_ON_NODE interact-prob.pep.xml interact-ipro.pep.xml 
#tpp2mzid interact-ipro.ptm.pep.xml #not currently using mzid output

#### These are post-processing scripts that I have written (Andy), here I am not doing an array job, so I can process a full run in one go
FILE_PREFIX=${PWD##*/}
pathtopythonqc=$2
python3 ${pathtopythonqc}/Convert_pepXML_toCSV.py interact-ipro.pep.xml ${FILE_PREFIX}
python3 ${pathtopythonqc}/CalculateFDR_and_threshold_v2.py ${FILE_PREFIX}_interact-ipro.pep.tsv rev_ 0.01
python3 ${pathtopythonqc}/calculate_psm_stats.py ${FILE_PREFIX}_interact-ipro.pep_thresholded.tsv $FILE_PREFIX

#Tidy up afterwards for file moving!
#gzip -k ${FILE_PREFIX}_interact-ipro.pep_thresholded.tsv
#gzip ${FILE_PREFIX}_interact-ipro.pep.tsv
#gzip *.pepXML
#gzip *.xml
#gzip *.mzML



