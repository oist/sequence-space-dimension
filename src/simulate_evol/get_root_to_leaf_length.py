#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from tqdm import tqdm
import pandas as pd
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
    PARSER.add_argument('-t', '--list_file', help='list of lists of files with trees',
                        required=True)
    PARSER.add_argument('-if', '--input_folder', help='folder with input trees',
                        required=False, default="")
    PARSER.add_argument('-of', '--output_folder', help='folder to write the outputs to',
                        required=False, default="")
    
    ARGS = vars(PARSER.parse_args())
    list_file=ARGS['list_file']
    in_folder=ARGS['input_folder']
    out_folder=ARGS["output_folder"]

    lists_by_gr=pd.read_csv(list_file, header=None)
    allconds=list(lists_by_gr[0])
    alllens_list=[]
    listpath='/bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/'
    for cond in range(4):
        #input is a list of names of tree files
        inlist=open(listpath+lists_by_gr[1][cond], 'r').read().split('\n')[:-1]
        alllens_list.append(simevol.get_root_to_leaf_length(inlist, path=in_folder))

    df=pd.DataFrame({'tree':allconds, 'root_to_leaf_lengths':alllens_list})
    df.to_csv(out_folder+'root_to_leaf_lengths_by_tree_type.tsv', index=False, sep='\t')
