### Based on DECOY string, want to calculate FDR and only write out global FDR below threshold ###

import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import sys
import argparse
from pathlib import Path

#decoy_string = "DECOY"
#QTHRESHOLD = 0.01

def get_is_decoy(row, decoy_string):

    #print("check decoy string",decoy_string)
    all_proteins = [row.loc["protein"]]
    if str(row.loc["alt_protein"]) != "nan": #Gives nan if empty for some reason
        all_proteins.extend(str(row.loc["alt_protein"]).split(";"))

    is_decoy = 1
    for protein in all_proteins:
        #print("p",protein)
        if protein.startswith(decoy_string) == False:   #Any proteins not a decoy, then it is not a decoy
            is_decoy = 0
            #print("\tFalse")
            break

    return pd.Series({'is_decoy': is_decoy})



#calculate FDR - Kerry's method
def calculateFDR(results_file, decoy_string, q_thresh,prob_column_header):
    #extract results to df
    #df=Extract_DF.extract_PTMprophet_IDent_df(results_file,PXD)
    df = pd.read_csv(results_file,sep="\t",index_col=None,low_memory=False)

    df[prob_column_header]=df[prob_column_header].astype(float)
    df=df.sort_values(by=prob_column_header,ascending=False)
    df=df.reset_index(drop=True)
    #set decoy and target

    #Old method was flawed (especially if : in protein names but possibly generally flawed)
    #df['Protein_count'] = df['protein'].str.count(":")+1
    #df['Decoy_count'] = df['protein'].str.count(decoy_string)
    #df['Decoy'] = np.where(df['Protein_count']==df['Decoy_count'],1,0)
    #df=df.drop(columns=['Protein_count', 'Decoy_count'])
    df['Decoy'] = df.apply(get_is_decoy,axis=1, args=[decoy_string])

    #target_count = row_count-decoy_count
    df['decoy_count']=df['Decoy'].cumsum()
    df['row_count']=(df.index)+1
    df['target_count']=df['row_count']-df['decoy_count']
    #FDR = decoy_count/target_count
    df['FDR']= df['decoy_count']/df['target_count']
    df['FDR'] = df['FDR'].astype(float)
    #q_val = min FDR at this position or lower
    df['q_value']=df['FDR']
    df['q_value']=df.iloc[::-1]['FDR'].cummin()
    #return as csv
    output = results_file[:-4] + "_withfdr.tsv"
    df.to_csv(output,index=False,sep="\t")
    print("Output written to ",output)

    df_thresholded = df[df["q_value"] < q_thresh]
    output = results_file[:-4] + "_thresholded.tsv"

    df_thresholded.to_csv(output, index=False, sep="\t")
    print("Output written to ", output)

    #FDR plots
    fig,ax = plt.subplots(figsize=(12,6))
    #ax1=df.plot.scatter(x='FDR',y='row_count')
    #ax2=df.plot.line(x='q_value',y='row_count')
    #plt.ylabel("PSM count")
    #plt.savefig('FDR.jpg',dpi=300)

    ax.plot(df[prob_column_header],df["q_value"])

    #df.plot.line(x='prob', y='FDR',ylim=(0,1),xlim=(1,0))
    plt.xlabel("Peptide Probability")
    plt.ylabel("q-value")
    output = results_file[:-4] + "FDR_score.jpg"
    plt.savefig(output,dpi=300)
    plt.close('all')
    print("Finished - Output written to ", output)

def calculateFDR_Peptide_Level(results_file, decoy_string, q_thresh,prob_column_header):
    #extract results to df
    #df=Extract_DF.extract_PTMprophet_IDent_df(results_file,PXD)
    df = pd.read_csv(results_file,sep="\t",index_col=None)

    df[prob_column_header]=df[prob_column_header].astype(float)
    df=df.sort_values(by=[prob_column_header,'count_of_PSMs'],ascending=False)
    df=df.reset_index(drop=True)
    #set decoy and target
    df['Decoy'] = df.apply(get_is_decoy, axis=1, args=[decoy_string])
    #target_count = row_count-decoy_count
    df['decoy_count']=df['Decoy'].cumsum()
    df['row_count']=(df.index)+1
    df['target_count']=df['row_count']-df['decoy_count']
    #FDR = decoy_count/target_count
    df['FDR']= df['decoy_count']/df['target_count']
    df['FDR'] = df['FDR'].astype(float)
    #q_val = min FDR at this position or lower
    df['q_value']=df['FDR']
    df['q_value']=df.iloc[::-1]['FDR'].cummin()
    #return as csv
    output = results_file[:-4] + "_pep_level_withfdr.tsv"
    df.to_csv(output,index=False,sep="\t")
    print("Output written to ",output)

    df_thresholded = df[df["q_value"] < q_thresh]
    output = results_file[:-4] + "_pep_level_thresholded.tsv"

    df_thresholded.to_csv(output, index=False, sep="\t")
    print("Output written to ", output)

    #FDR plots
    fig,ax = plt.subplots(figsize=(12,6))
    #ax1=df.plot.scatter(x='FDR',y='row_count')
    #ax2=df.plot.line(x='q_value',y='row_count')
    #plt.ylabel("PSM count")
    #plt.savefig('FDR.jpg',dpi=300)

    ax.plot(df[prob_column_header],df["q_value"])

    #df.plot.line(x='prob', y='FDR',ylim=(0,1),xlim=(1,0))
    plt.xlabel("Peptide Probabilty")
    plt.ylabel("q-value")
    output = results_file[:-4] + "_pep_level_FDR_score.jpg"
    plt.savefig(output,dpi=300)
    plt.close('all')
    print("Finished - Output written to ", output)



### TO DO Need to add in keyword for pp_prob or ip_prob - to switch between column used

def main():
    '''
    executes the arguments and functions to run and optimise liftoff
    '''
    parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
    Calculate FDR and thresholds for a tsv output from MSfragger.
    -------------------------------------------------------------
     ''',
     epilog="written by Andy Jones, edited by Helen Rebecca Davison") 
    parser.add_argument('-f', '--results_file', \
                    help="The unthresholded results tsv from peptide prophet \
                        called 'PXD000000_interact-ipro.pep.tsv' or 'PXD000000_interact-ipro.pep_thresholded_peptide_final_collapsed.tsv'",
                    required=True)
    parser.add_argument('-d','--decoy', \
                    help="the decoy string, like 'rev_'",
                    required=True)
    parser.add_argument('-t','--threshold', nargs='?', const=0.01, default=0.01, type=float, \
                    help="The cutoff value for q-value. The default is 0.01.")
    parser.add_argument('-p','--peptide_level',  action='store', nargs='*',\
                    help="Calculate FDR at the peptide level")
    args = parser.parse_args()
    
    results_file = args.results_file
    decoy_string = args.decoy
    threshold = args.threshold
    peptide_level = args.peptide_level
    
    #sys.tracebacklimit = 0

    if not Path(results_file).exists():
            raise Exception("Results file " + results_file +" does not exist")

    print("\033[33m{}\033[0;0m".format("Using results file: "))
    print(str(results_file))
    print("\033[33m{}\033[0;0m".format("Using decoy string: "))
    print(decoy_string)
    print("\033[33m{}\033[0;0m".format("Using threshold: "))
    print(threshold)

    if peptide_level is not None:
        calculateFDR_Peptide_Level(results_file, decoy_string, threshold,"ip_prob")
        print("Calculating peptide-level stats and thresholding at " + str(threshold))
    else:
        calculateFDR(results_file, decoy_string, threshold,"pp_prob")
        print("Calculating PSM stats and thresholding at " + str(threshold))
    

main()