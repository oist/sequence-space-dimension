#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from tqdm import tqdm
import os 
import pandas as pd
import numpy as np
import sys
import pickle
import re

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
    
    ARGS = vars(PARSER.parse_args())
    tree=ARGS['tree']
    in_folder=ARGS['input_folder']
    out_folder=ARGS["output_folder"]

    
    if '.pickle' in tree:
        #input is one file
        inlist=[tree]
    else: #input is a list of names of tree files
        inlist=open(tree, 'r').read().split('\n')[:-1]

    rows=[]
    #iterate over trees 
    for file in tqdm(inlist):
        ogname='ali_'+re.findall('tree_(.+).pickle', os.path.basename(file))[0]

        leaves, intnodes, leavesancs = simevol.get_ancestors_seqs_from_tree(in_folder+file, internal_to_fasta=False)
        alpha_anc, ninvar_anc, _= simevol.get_usage_aa_per_site(intnodes, alphabet='nt', colfreqs_out=False)
        alpha_ancleaf, ninvar_ancleaf, _= simevol.get_usage_aa_per_site(leavesancs, alphabet='nt', colfreqs_out=False)
        rows.append([ogname, np.mean(alpha_anc), np.median(alpha_anc), ninvar_anc, np.mean(alpha_ancleaf),np.median(alpha_ancleaf),  ninvar_ancleaf])
    
leafancdf=pd.DataFrame(rows, columns=['OG_id', 'mean_alphabet_size_anc', 'median_alphabet_size_anc', 'num_invar_sites_anc', 
                                      'mean_alphabet_size_leaves_anc',  'median_alphabet_size_leaves_anc',  'num_invar_sites_leaves_anc'])
leafancdf['num_var_sites_anc']=300-leafancdf['num_invar_sites_anc']
leafancdf['num_var_sites_leaves_anc']=300-leafancdf['num_invar_sites_leaves_anc']
leafancdf.to_csv(out_folder+os.path.splitext(os.path.basename(tree))[0]+'_usage_leaves_anc.csv', index=False)

