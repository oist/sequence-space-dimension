#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from tqdm import tqdm
import os
import pandas as pd
import sys
import pickle

sys.path.append("/bucket/KondrashovU/seq_space/scripts/seq_space_lib/")
sys.path.append("../scripts/seq_space_lib/")
import sequence_space_lib as seqsp

sys.path.append("/bucket/KondrashovU/seq_space/scripts/simulate_prot_evol_lib/")
sys.path.append("../scripts/simulate_prot_evol_lib/")
import simulate_prot_evol_lib as simevol


if __name__ == "__main__":

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-t', '--tree', help='list of files with trees',
                        required=True)
    PARSER.add_argument('-if', '--input_folder', help='folder with input trees',
                        required=False, default="")
    PARSER.add_argument('-of', '--output_folder', help='folder to write to simulation outputs to',
                        required=False, default="")
    PARSER.add_argument('-a', '--return_all', type=eval, choices=[True, False],  
                        help='set to True to include the list with all dn/ds values for each pair of seqs',
                        required=False, default=False)
    PARSER.add_argument('-s', '--keep_suffix', type=eval, choices=[True, False],  
                        help='set to True to include the file extension in the output file name',
                        required=False, default=False)
    #Optional

    ARGS = vars(PARSER.parse_args())

    tree=ARGS['tree']
    in_folder=ARGS['input_folder']
    out_folder=ARGS["output_folder"]
    return_all=ARGS["return_all"]
    keep_suffix=ARGS["keep_suffix"]

    #initialize empty arrays for outputs
    results=[]

    if '.pickle' in tree:
        #input is one file
        inlist=[tree]
    else: #input is a list of names of tree files
        inlist=open(tree, 'r').read().split('\n')[:-1]

    #iterate over trees 
    for treefile in tqdm(inlist):
        ogname='ali_'+re.findall('tree_(.+).pickle', os.path.basename(treefile))[0]
        with open(in_folder+treefile, "rb") as intree:
            atree = pickle.load(intree)
        result=simevol.dn_ds_from_tree(atree, return_all=return_all)
        results.append([ogname]+list(result))

    if return_all:
        colnames=['OG_id', 'dnds_mean_tree', 'dnds_median_tree', 'dnds_var_tree', 'dnds_list_tree']
    else:
        colnames=['OG_id', 'dnds_mean_tree', 'dnds_median_tree', 'dnds_var_tree']

    df=pd.DataFrame(results, columns=colnames).round(5)

    if keep_suffix:
        outfn=os.path.basename(tree)
    else:
        outfn=os.path.splitext(os.path.basename(tree))[0]
    with open(out_folder+'dnds_from_tree_'+outfn+'.tsv', 'w+') as csvout:
        df.to_csv(csvout, sep='\t', index=False)
