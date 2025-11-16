#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from tqdm import tqdm
import os
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
    PARSER.add_argument('-t', '--tree', help='list of files with trees',
                        required=True)
    PARSER.add_argument('-if', '--input_folder', help='folder with input trees',
                        required=False, default="")
    PARSER.add_argument('-of', '--output_folder', help='folder to write to simulation outputs to',
                        required=False, default="")
    #Optional

    #output files parameters
    PARSER.add_argument('-a', '--include_ancestors', type=eval, help='set to true if dnds and average distance to be calculated with and without ancestors, (False - only the leaves)',
                         required=False, default=False)
    PARSER.add_argument('-df', '--dist_to_file', type=eval, help='set to True to write pairwise distances to matrix, False by default',
                         required=False, default=False)
    PARSER.add_argument('-d', '--dim_only', type=eval, choices=[True, False], help='only get the dimensionality into the final summary table, by default: False', 
                        required=False, default=False)
    PARSER.add_argument('-w', '--win_size', type=float,  help='window size too use when calculating dimensionality', 
                        required=False, default=0.3)
    PARSER.add_argument('-r', '--dim_range_from_data', type=eval,  help='set to True to use distance range for dimensionality calculation from the data (not 0-1), False by default', 
                        required=False, default=False)
    PARSER.add_argument('-do', '--dim_out', type=eval, help='set to True to write dimensionality output file, False by default',
                         required=False, default=False)
    ARGS = vars(PARSER.parse_args())

    tree=ARGS['tree']
    in_folder=ARGS['input_folder']
    out_folder=ARGS["output_folder"]
    dim_only=ARGS["dim_only"]
    include_internal=ARGS['include_ancestors']
    dist_to_file_bool=ARGS['dist_to_file']
    win_size=ARGS['win_size']
    dim_out_bool=ARGS['dim_out']
    dim_range_from_data=ARGS['dim_range_from_data']

    #initialize empty arrays for outputs
    results=[]

    if '.pickle' in tree:
        #input is one file
        inlist=[tree]
    else: #input is a list of names of tree files
        inlist=open(tree, 'r').read().split('\n')[:-1]

    #iterate over trees 
    for treefile in tqdm(inlist):
        params_from_name=re.findall('tree_(\S+)_(\S+)_(\S+)_(\d+)_(continuousF\S+run\d+)_f(0\.\d+)_fnum(\d+)_0.pickle',
                                           treefile)[0]
        if dist_to_file_bool:
            dist_to_file=out_folder+'dist_'+os.path.splitext(treefile)[0][5:]
        else:
            dist_to_file=None
        
        if dim_out_bool:
            dim_out=out_folder+'dim_'+os.path.splitext(treefile)[0][5:]
        else:
            dim_out=None

        if params_from_name[0]=='notree' or params_from_name[0]=='no_tree':
            include_internal=False
        result=list(simevol.get_stats_from_simulated_tree(in_folder+treefile, dim_only=dim_only, verbose=True, include_internal=include_internal, 
                                   dist_to_file=dist_to_file, win_size_perc=win_size, out_dist=False, dim_out=dim_out, 
                                   dim_range_from_data=dim_range_from_data))
        
        result.extend(params_from_name)
        results.append(result)

    colnames_param=['mode', 'tree_len_scale', 'branch_length_distr', 'gamma', 'suffix', 'fraction_allowed_subs', 'replicate']

    if dim_only:
        if include_internal:
            colnames=['mean_dist',
                      'dimensionality',
                      'mean_dist_anc',
                      'dimensionality_anc',
                      'mean_dist_leaves_anc',
                      'dimensionality_leaves_anc'] + colnames_param
            out_file_name='dist_dim_w_anc'
        else:
            colnames=['mean_dist',
                      'dimensionality'] + colnames_param
            out_file_name='dist_dim'

    else:
        if include_internal:
            colnames=['dnds_mean',
                      'mean_dist',
                      'dimensionality',
                      'dnds_mean_anc',
                      'mean_dist_anc',
                      'dimensionality_anc',
                      'dnds_mean_leaves_anc',
                      'mean_dist_leaves_anc',
                      'dimensionality_leaves_anc'] + colnames_param
            
            out_file_name='dnds_dist_dim_w_anc'
        else:
            colnames=['dnds_mean',
                      'mean_dist',
                      'dimensionality'] + colnames_param
            out_file_name='dnds_dist_dim'

    df=pd.DataFrame(results, columns=colnames).round(5)
    
    if dim_range_from_data:
        dimsuffix='_dim_range_from_data'
    else:
        dimsuffix=''

    outfn=os.path.splitext(os.path.basename(tree))[0]
    with open(out_folder+out_file_name+'_'+outfn+'_winsize{}.csv'.format(str(win_size)+dimsuffix), 'w+') as csvout:
            df.to_csv(csvout, index=False)
