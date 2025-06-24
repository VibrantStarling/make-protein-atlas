### Read local CSV format and take best scoring PEPTIDE ID, and collapse - (keeping count of PSMs with Prob > threshold)
### Might want to trial use of pandas? ###

# usage:
# python Collapse_by_max_prob_v2.py inputlocation/PXD00000_interact-ipro.pep_thresholded.tsv [-m]
# -m tells the script to collapse on modified peptides instead of peptides

import pandas as pd
import numpy as np
import sys

def read_and_collapse(infile, pep_type):

    df = pd.read_csv(infile,sep="\t")
    #pivot_df = df.pivot_table(index=["peptide"],values = ["prob"], aggfunc={np.max,np.count_nonzero})
    df_only_cols = df[[pep_type,"pp_prob"]]
    df_only_cols = df_only_cols.rename(columns={"pp_prob":"count_of_PSMs"})
    group_df = df_only_cols.groupby(by=pep_type).count()
    outfile = infile[:-4] + "_" + pep_type + "_count.tsv"
    #group_df.to_csv(outfile,sep="\t")

    #outfile = infile[:-4] + "_collapsed.tsv"
    #pivot_df.to_csv(outfile,sep="\t",index=False)

    df = df.sort_values("pp_prob", ascending=False)    #check sorted, I think it will have been sorted by the Prophets anyway
    df_no_dups = df.drop_duplicates(subset=pep_type, keep='first', inplace=False)
    df_no_dups = df_no_dups.set_index(pep_type)
    outfile = infile[:-4] + "_" + pep_type + "_no_dups.tsv"
    #df_no_dups.to_csv(outfile,sep="\t")

    final_df = pd.concat([df_no_dups,group_df],axis=1)
    outfile = infile[:-4] + "_" + pep_type.replace('_','') + "_final_collapsed.tsv"
    final_df.to_csv(outfile, sep="\t",index_label=pep_type)
    print("Finished - output written to",outfile)




    #df['Protein_loc'] = df.groupby('peptide')['peptide'].transform('count')


def testing():
    outfolder = "D:/temp/PXD028712/"
    #infile= outfolder + "PXD028712_pep_prophet_mini3_thresholded.tsv"
    infile= outfolder + "PXD028712_pep_prophet_thresholded.tsv"
    read_and_collapse(infile)
    

if len(sys.argv) < 2:
    print("Exit - expected usage\npython ",  sys.argv[0]," inputlocation/PXD00000_interact-ipro.pep_thresholded.tsv [-m] \n"
                                                         "Script collapses based on max probability, and does PSM count with"
                                                         "unique modified peptide string being the unit of grouping. \n"
                                                         "specify -m to collapse on modified peptides instead of peptides")

else:
    if (len(sys.argv) > 2) and (sys.argv[2] == "-m"):
        pep_type = "mod_peptide"
        print("collapsing on modified peptides")
    else:
        pep_type = "peptide"
        print("collapsing on peptides")
    read_and_collapse(sys.argv[1], pep_type)
    

    







