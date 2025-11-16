#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from tqdm import tqdm
import numpy as np
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
    PARSER.add_argument('-p', '--params', help='simulation parameters and output file suffix as underscore separated string',
                        required=True)
    PARSER.add_argument('-of', '--output_folder', help='folder to write to simulation outputs to',
                        required=False, default="")
    #Optional
    
    PARSER.add_argument('-r', '--replicates', type=int, help='number of replicates to simulate for every value of F, if f_type==discrete or custom, number of simulations in total if f_type==continuous',
                        required=False, default=10)
    PARSER.add_argument('-ft', '--f_type', help='type of F values to use: discrete -(0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5), continuous - uniformly distributed between 0.07 and 0.5, or custom _ delimited  list',
                        required=False, default='discrete')
    PARSER.add_argument('-fu', '--f_upper',  type=float, help='upper limit of F values for continuousF',
                        required=False, default=0.5)
    PARSER.add_argument('-fl', '--f_lower',  type=float, help='lower limit of F values for continuousF',
                        required=False, default=0.07)
    PARSER.add_argument('-kt', '--keep_tree', type=eval, choices=[True, False], help='keep the tree the same for different values of F, by default: False - every simulation uses a new tree', 
                        required=False, default=False)
    
    #Model parameters
    PARSER.add_argument('-s', '--seq_len',type=int, help='sequence length in nucleotides',
                        required=False, default=900)
    PARSER.add_argument('-l', '--n_leaves', type=int, help='number of leaves in the tree / sequences in the protein family',
                        required=False, default=300)
    PARSER.add_argument('-er', '--exp_rate', type=float, help = 'rate of exponential distributions used to simulate branch lengths when they decrease with tree depth', 
                        required=False, default=0.05349)
    PARSER.add_argument('-bs', '--branch_len_scaling', help='way to scale branch lengths as a function of the distance from root, default = False (no scaling)', 
                        required=False, default='None')
    PARSER.add_argument('-bf', '--scale_factor', type=float, help='Factor to scale branch lengths by depth, only when increasing with distance from root',
                         required=False, default=0)

    
    #output files parameters
    PARSER.add_argument('-a', '--include_ancestors', type=eval, help='set to true if dnds and average distance to be calculated with and without ancestors, (False - only the leaves)',
                         required=False, default=False)
    PARSER.add_argument('-df', '--dist_to_file', type=eval, help='set to True to write pairwise distances to matrix, False by default',
                         required=False, default=False)
    PARSER.add_argument('-d', '--dim_only', type=eval, choices=[True, False], help='only get the dimensionality into the final summary table, by default: False', 
                        required=False, default=False)
    PARSER.add_argument('-w', '--win_size', type=float,  help='window size too use when calculating dimensionality', 
                        required=False, default=0.4)
    
    ARGS = vars(PARSER.parse_args())

    params = ARGS["params"]
    out_folder=ARGS["output_folder"]
    replicates = ARGS["replicates"]
    f_type=ARGS["f_type"]
    f_upper=ARGS['f_upper']
    f_lower=ARGS['f_lower']
    keep_tree=ARGS["keep_tree"]
    ntseqlen=ARGS["seq_len"]
    n_leaves=ARGS["n_leaves"]
    dim_only=ARGS["dim_only"]
    exp_rate=ARGS['exp_rate']
    shortleaves=ARGS['branch_len_scaling']
    scale_factor=ARGS['scale_factor']
    include_internal=ARGS['include_ancestors']
    dist_to_file_bool=ARGS['dist_to_file']
    win_size=ARGS['win_size']

    aaseqlen=int(int(ntseqlen)/3)
    symparams=params.split('_') #params in the format tree_1.2_exp_gamma{_batchID}

    #initialize empty arrays for outputs
    realfs=[]
    results=[]
    newtree=False #by default don't simulates a new tree for each replicate

    #get f values
    if f_type=='discrete': #if the f values are discrete
        f_values=[0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5] #standard set of f values
    elif '0' in f_type: #if the f values are cusom set
            f_values=[float(i) for i in f_type.split('_')]
    elif f_type=='continuous': #if the f values are randomly uniformly sampled
        f_values=np.random.uniform(f_lower, f_upper, int(replicates))
        #if the f values are randomly uniformly sampled, only simulate 1 replicate per f value 
        #(replicates then refer to the number of random f values sampled)
        replicates=1
        
    #choose if we want to use the same tree for all the replicates and f values or a new one for each replicate and f value
    if bool(keep_tree): #simulate the tree that will be used in all replicates
        simtree = simevol.simulate_tree_with_depth_scaled_branch_lengths(int(n_leaves), scale_factor=float(scale_factor),
                                                                 length_scale=float(symparams[1]), distribution=symparams[2], 
                                                                 exp_rate=float(exp_rate), shortleaves=shortleaves,
                                                                 starlike=symparams[0], starlike_limit=symparams[2])

    else: #set "simulate a new tree for each replicate"-switch to True
        simtree=None
        newtree=True
    
    for i in range(int(replicates)):
        print("Working on iteration ", i, flush=True)
        for findex, f in enumerate(f_values):
            outtree=out_folder+'tree_'+params+'_f'+"{:.3f}".format(f)+ '_fnum'+ str(findex)+'_'+str(i) + '.pickle'
            outfasta=out_folder+'ali_'+params+'_f'+"{:.3f}".format(f) + '_fnum'+ str(findex)+'_'+str(i)+'.fasta'
            if dist_to_file_bool:
                dist_to_file=out_folder+'dist_'+params+'_f'+"{:.3f}".format(f)+ '_fnum'+ str(findex)+'_'+str(i)
            else:
                dist_to_file=None
            result=simevol.simulate_evolve_tree_get_stats(symparams, int(ntseqlen), aaseqlen, int(n_leaves), f, 
                                                  exp_rate=float(exp_rate), shortleaves=shortleaves, 
                                                  simtree=simtree, newtree=newtree, dim_only=dim_only, 
                                                  outtree=outtree, outfasta=outfasta, 
                                                  scale_factor=float(scale_factor), out_dist=False,  
                                                  include_internal=include_internal, dist_to_file=dist_to_file, 
                                                  win_size_perc=win_size)
            results.append(result)
        realfs.extend(f_values)
        print("Finished", flush=True)


    if dim_only:
        if include_internal:
            colnames={0:'mean_dist', 1:'dimensionality', 2:'mean_dist_anc', 
                      3:'dimensionality_anc', 4:'mean_dist_leaves_anc', 
                      5:'dimensionality_leaves_anc'}
            out_file_name='dist_dim_w_anc'
        else:
            colnames={0:'mean_dist', 1:'dimensionality' }
            out_file_name='dist_dim'
        df=pd.DataFrame(results).rename(columns=colnames)

    else:
        if include_internal:
            colnames={0:'dnds_mean', 1:'mean_dist', 2:'dimensionality', 
                      3:'dnds_mean_anc', 4:'mean_dist_anc', 5:'dimensionality_anc', 
                      6:'dnds_mean_leaves_anc', 7:'mean_dist_leaves_anc', 8:'dimensionality_leaves_anc'}
            out_file_name='dnds_dist_dim_w_anc'
        else:
            colnames={0:'dnds_mean', 1:'mean_dist', 2:'dimensionality' }
            out_file_name='dnds_dist_dim'
        df=pd.DataFrame(results).rename(columns=colnames)

    flen=len(realfs)
    df['fraction_allowed_subs']=realfs
    df['mode']=[symparams[0]]*flen
    df['tree_len_scale']=[symparams[1]]*flen
    df['exp_rate']=[exp_rate]*flen
    df['branch_len_scaling']=[shortleaves]*flen
    df['branch_length_distr']=[symparams[2]]*flen
    df['gamma']=[symparams[3]]*flen
    df['scale_factor']=[scale_factor]*flen
    if  f_type!='continuous':
        df['replicate']= simevol.get_replicate_indices(int(replicates), len(f_values))
    if len(symparams)>4:
        df['suffix']=[symparams[4]]*flen

    meta_params='_'.join([f_type,'keep_tree', str(keep_tree), 'reps'+str(replicates), 
                         shortleaves, 'exp_rate'+str(exp_rate), 'scale'+str(scale_factor)])
    
    with open(out_folder+out_file_name+'_'+meta_params +'_'+params + '.csv', 'w') as csvout:
            df.to_csv(csvout, index=False)
