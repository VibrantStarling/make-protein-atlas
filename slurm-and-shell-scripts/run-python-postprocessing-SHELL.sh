#! /bin/bash

# usage 
# bash run-python-postprocessing-SHELL.sh pxd-list.txt path/to/run-python-postprocessing-BATCH.sh python/qc/directory

for i in $(cut -f1 ${1})
do cd ${i}; 
sbatch --array=1-1 ${2} ${3} ${i}
cd ../
done
