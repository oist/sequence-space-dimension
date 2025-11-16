
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
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
    PARSER.add_argument('-n', '--nseq', type=int, help='number of sequences per family',
                        required=True)
    PARSER.add_argument('-of', '--out_folder', help='folder to write to simulation outputs to',
                        required=False, default="")
    PARSER.add_argument('-m', '--fmatrpath', help='input folder with matrices',
                        required=False, default='params/fracs_and_matr_allowed_for_nseq_10linspace.pickle')
    PARSER.add_argument('-r', '--repl', type=int, help='number of replicates per matrix and number of sequences',
                        required=False, default=5)
    PARSER.add_argument('-w', '--win_size', type=float,  help='window size too use when calculating dimensionality', 
                        required=False, default=0.3)
    PARSER.add_argument('-nm', '--newmatr', type=bool,  help='generate new matrices', required=False, default=True)
    PARSER.add_argument('-s', '--starlike', type=bool, help='set to true if starlike tree', required=False, default=False)

    ARGS = vars(PARSER.parse_args())

    nseq=ARGS['nseq']
    repl=ARGS['repl']
    out_folder=ARGS['out_folder']
    fmatrpath=ARGS['fmatrpath']
    win_size_perc=ARGS['win_size']
    newmatr=ARGS['newmatr']
    starlike=ARGS['starlike']

    if starlike==True:
        lenscale=0.7
        file_suffix='_starlike'
    else:
        lenscale=1
        file_suffix=''

    rows=[]
    seqlen=300
    with open(fmatrpath, "rb") as f:
        fmatr = pickle.load(f)
    colnames=['fraction_allowed_subs', 'nseq', 'rep', 'dist', 'dist_var', 'dim', 'dim_ci']
    print('res:'+'\t'.join(colnames))
    for f, matrix in fmatr:
        for k in tqdm(range(repl)):
            try:
                if newmatr==True:
                    matrix=simevol.generate_matrix_AA(seqlen, fraction=f)
                else:
                    pass
                
                random_sequence = simevol.generate_nt_seq(900, 1, seqfrom='matrix', matrix=matrix)
                simtree = simevol.simulate_tree_with_depth_scaled_branch_lengths(nseq, starlike=starlike, length_scale=lenscale)
                leaves, nonsyn_subs, syn_subs, subs_attempts, tree, statslist = simevol.evolve_tree(random_sequence,matrix, simtree,
                                                                                                neighbors={}, limit_on='all_subs', gamma=0, 
                                                                                                equal_f_split=True, universal_neighbors=True)
                print('sequences done')
                dist_vect, dist, (dim, dimr2, dim_ci)=simevol.get_dnds_dist_dim_stats_from_seq(leaves,
                                            dist_to_file=out_folder+'simu_tree{}_f{:.3f}_nseq{}_rep{}.dist'.format(file_suffix, f, nseq, k),
                                            win_size_perc=win_size_perc, dim_only=True, out_dist=True, nt=True, 
                                            dim_out=None, dim_stats_out=True)
                rows.append(f)
                rows.append(nseq)
                rows.append(k)
                rows.append(dist)
                rows.append(np.var(dist_vect))
                rows.append(dim)
                rows.append(dim_ci)
            except ValueError:
                rows.append(f)
                rows.append(nseq)
                rows.append(k)
                rows.append(np.nan)
                rows.append(np.nan)
                rows.append(np.nan)
                rows.append(np.nan)
            print('res:{:.4f}\t{:.0f}\t{:.0f}\t{:.5f}\t{:.5f}\t{:.4f}\t{:.4f}'.format(*rows[-7:]))
    df=pd.DataFrame(np.array(rows).reshape(len(fmatr)*repl,7), columns=colnames)
    df.to_csv(out_folder+'simu_tree{}_nseq{}_repl{}.csv'.format(file_suffix, nseq, repl), float_format='%.5f', index=False)
