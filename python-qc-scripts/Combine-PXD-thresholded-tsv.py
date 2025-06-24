# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:27:01 2025

@author: telma
"""

import argparse
import pandas as pd
import numpy as np 
from glob import glob
from pathlib import Path
import re


def get_files(pxd_list, modpeptide):
    PXDs = pd.read_csv(pxd_list, sep='\t', names = ['pxd_code','type'])
    files = []
    for pxd in PXDs['pxd_code']:
        if modpeptide is not None:
            files = files + glob(pxd+'/'+pxd+'*_modpeptide_final_collapsed.tsv')
        else: 
            files = files + glob(pxd+'/'+pxd+'*_peptide_final_collapsed.tsv')
    print("There are "+ str(len(files)) + " files to merge.")
    return files


def count_pep_occurance(files):
    # a dictionary of {peptide: number_of_occurances}
    peptide_occurance = {}
    for file in files:
        print("Adding counts for " + file)
        file_peptides = []
        for chunk in pd.read_csv(file, sep = '\t', chunksize=2):
            file_peptides += list(chunk['peptide'])
        for peptide in set(file_peptides):
                peptide_occurance[peptide] = peptide_occurance.get(peptide, 0) + 1
    return peptide_occurance

def combine_tsvs(files, output_file, peptide_occurance):
    # combine a list of tsv files into one
    headers = pd.read_csv(files[0], nrows=1, header=None, sep = "\t")
    headers[headers.shape[1]+1] = 'occurances'
    headers.to_csv(output_file, mode="a", index=False, header=None, sep="\t")
    for tsv_file_name in files:
        chunk_container = pd.read_csv(tsv_file_name, chunksize=2, sep="\t")
        for chunk in chunk_container:
            chunk['occurances'] = chunk['peptide'].map(peptide_occurance)
            chunk.to_csv(output_file, mode="a", index=False, header=None, sep="\t")
    #
    print("Combined files saved to " + output_file )
    return

def main():
    '''
    executes the arguments and functions to run and optimise liftoff
    '''
    parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
    Combines thresholded tsvs produced by 'CalculateFDR_and_threshold.py' scripts.
    -------------------------------------------------------------
    
         ''',
         epilog="written by Helen Rebecca Davison") 
    parser.add_argument('-l', '--pxd_list', \
                        help="Tab delimited list of PXD codes and the type of analysis (LF, TMT,iTRAQ)\
                        with one entry on each line",
                        required=True)
    parser.add_argument('-o', '--outname', \
                        help="Prefix to give the output",
                        required=True)
    parser.add_argument('-m','--modpeptide',  action='store', nargs='*',\
                        help="Combine .tsv results files that have been collapsed by mod_peptide.")

    args = parser.parse_args()
    
    
    pxd_list = args.pxd_list
    output_file = str(Path(args.outname+"_PepAtlas-COMBINED-thresholded.tsv"))
    modpeptide = args.modpeptide
    
    # error handling
    if Path(output_file).exists():
        raise Exception(output_file + " exists. Please remove or choose a different prefix name.")

              
    ## get files for processing
    files = get_files(pxd_list, modpeptide)
    
    if len(files) == 0:
        raise Exception("No files found. Check if you need to specify '-m' and make sure your .tsv have been run through Collapse_by_max_prob_v2.py")
    
    peptide_occurance = count_pep_occurance(files)
    combine_tsvs(files, output_file, peptide_occurance)
    

main()