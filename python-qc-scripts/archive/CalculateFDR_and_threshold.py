### Based on DECOY string, want to calculate FDR and only write out global FDR below threshold ###

import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import sys

#decoy_string = "DECOY"
#QTHRESHOLD = 0.01

def get_is_decoy(row):

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
    df['Decoy'] = df.apply(get_is_decoy,axis=1)

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
    df['Protein_count'] = df['protein'].str.count(":")+1
    df['Decoy_count'] = df['protein'].str.count(decoy_string)
    df['Decoy'] = np.where(df['Protein_count']==df['Decoy_count'],1,0)
    df=df.drop(columns=['Protein_count', 'Decoy_count'])
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
    output = results_file[:-4] + "_pep_level__thresholded.tsv"

    df_thresholded.to_csv(output, index=False, sep="\t")
    print("Output written to ", output)

    #FDR plots
    fig,ax = plt.subplots(figsize=(12,6))
    #ax1=df.plot.scatter(x='FDR',y='row_count')
    #ax2=df.plot.line(x='q_value',y='row_count')
    #plt.ylabel("PSM count")
    #plt.savefig('FDR.jpg',dpi=300)

    ax.plot(df["prob"],df["q_value"])

    #df.plot.line(x='prob', y='FDR',ylim=(0,1),xlim=(1,0))
    plt.xlabel("Peptide Probabilty")
    plt.ylabel("q-value")
    output = results_file[:-4] + "_pep_level_FDR_score.jpg"
    plt.savefig(output,dpi=300)
    plt.close('all')
    print("Finished - Output written to ", output)


def testing():
    file_to_process = "D:/temp/PXD028712/PXD028712_pep_prophet.tsv"
    calculateFDR(file_to_process,"DECOY",0.01)

### TODO Need to add in keyword for pp_prob or ip_prob - to switch between column used
if len(sys.argv)!= 4 and len(sys.argv)!=5:
    print("Exit - expected usage\npython ",  sys.argv[0]," PXD00000_interact-ipro.pep.tsv decoy_string qvaluethreshold [TRUE]\n"
          "e.g. python3 ",  sys.argv[0]," D:/temp/PXD028712/PXD028712_pep_prophet.tsv DECOY 0.01\n"
                                        "Purpose of the script is to calculate q-values via DECOY string FDR method and threshold at given threshold\n"
                                        "writing out an output with and without the threshold\n"
                                        "There is a 4th optional argument to specify to process a peptide-level tsv file")
else:

    decoy_string = sys.argv[2]
    print("Using decoy string",decoy_string)
    if len(sys.argv) == 4:
        calculateFDR(sys.argv[1], sys.argv[2], float(sys.argv[3]),"pp_prob")
        print("Calculating PSM stats and thresholding at ",sys.argv[3])
    else:
        calculateFDR_Peptide_Level(sys.argv[1], sys.argv[2], float(sys.argv[3]),"ip_prob")
        print("Calculating peptide-level stats and thresholding at ",sys.argv[3])



