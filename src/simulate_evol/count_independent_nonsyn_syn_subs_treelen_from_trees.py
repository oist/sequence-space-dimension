#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from tqdm import tqdm
import os 
import pandas as pd
import numpy as np
import sys

sys.path.append("../scripts/seq_space_lib/")
import sequence_space_lib as seqsp

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
    for treefile in tqdm(inlist):
        rows.append(simevol.count_independent_nonsyn_syn_subs_treelen(in_folder+treefile))
    df=pd.DataFrame(np.array(rows).reshape(len(rows),4), columns=['OG_id', 'nonsyn_subs_count', 'syn_subs_count', 'total_treelength'])
    df.to_csv(out_folder+os.path.splitext(os.path.basename(tree))[0]+'_total_nonsyn_subs_count.csv', index=False)