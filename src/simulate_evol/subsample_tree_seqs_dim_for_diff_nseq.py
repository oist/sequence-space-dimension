#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import re
import os

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
    PARSER.add_argument('-f', '--flist', help='list of pdistm files (matrices) to subsample',
                        required=True)
    PARSER.add_argument('-if', '--input_folder', help='folder with input trees',
                        required=False, default="")
    PARSER.add_argument('-of', '--output_folder', help='folder to write to simulation outputs to',
                        required=False, default="")
    PARSER.add_argument('-w', '--win_size', type=float,  help='window size too use when calculating dimensionality', 
                        required=False, default=0.3)
    PARSER.add_argument('-m', '--matr_to_file', type=bool, help='set to True to write subsampled matrices to separate pdistm files',
                         required=False, default=False)
    PARSER.add_argument('-n', '--nsamples', type=int, help='number of replicate subsamplings of the matrix',
                        required=False, default=5)
    PARSER.add_argument('-r', '--dim_range_from_data', type=eval,  help='set to True to use distance range for dimensionality calculation from the data (not 0-1), False by default', 
                        required=False, default=False)
    PARSER.add_argument('-s', '--nseqs', type=eval, help='number of seqeunces to subsample from the matrix',
                        required=False, default='[100,1000]')
    
    
    ARGS = vars(PARSER.parse_args())
    flist=ARGS['flist']
    in_folder=ARGS['input_folder']
    out_folder=ARGS["output_folder"]
    win_size_perc=ARGS['win_size']
    matr_to_file=ARGS['matr_to_file']
    nsamples=ARGS['nsamples']
    dim_range_from_data=ARGS['dim_range_from_data']
    nseqs=ARGS['nseqs']
    
    inlist=open(flist, 'r').read().split('\n')[:-1]

    rows=[]
    colnames=['fraction_allowed_subs', 'nseq', 'rep', 'sample', 'dist', 'dist_var', 'dim', 'dim_ci']
    print('res:'+'\t'.join(colnames))

    #iterate over matrices 
    for file in tqdm(inlist):
        #get f and replicate number from this matrix
        starlike, f, rep=re.findall('simu_(.*)_f(.+)_nseq10000_rep(.+).dist.pdistm', os.path.basename(file))[0]
        for nseq in nseqs:
            subsampled_matr=seqsp.submatrices(in_folder+file,nseq,nsamples)
            for sample_ind, pdistm in enumerate(subsampled_matr):
                try:
                    #record matrix to file
                    if matr_to_file:
                        out_matr_name=out_folder+'simu_{}_f{}_nseq{}_rep{}_sample{}.dist'.format(starlike, f, nseq, rep, sample_ind)
                        distdf=pd.DataFrame.from_records(pdistm)
                        with open('%s.pdistm' % (out_matr_name), 'w+') as outm:
                            distdf.to_csv(outm, index=False, header=False)
                    #get flat matrix
                    arrlen= (nseq*nseq-nseq)/2
                    pdistm_flat, _=seqsp.matrix_to_array(pdistm)
                    #get dim
                    dim, dimr2, dim_ci=seqsp.max_k((pdistm_flat, arrlen), nparray=False, from_file=False, 
                                    flat=True, plot=False, win_size_perc=win_size_perc, out_file=False, 
                                    fname=None, ci_r2_out_flat=True, reg_range_from_data=dim_range_from_data)

                    rows.append(f)
                    rows.append(nseq)
                    rows.append(rep)
                    rows.append(sample_ind)
                    rows.append(np.mean(pdistm_flat))
                    rows.append(np.var(pdistm_flat))
                    rows.append(dim)
                    rows.append(dim_ci)
                    #subsample this matrix to get 1000 and 100 sequences 
                except ValueError:
                    rows.append(f)
                    rows.append(nseq)
                    rows.append(rep)
                    rows.append(sample_ind)
                    rows.append(np.nan)
                    rows.append(np.nan)
                    rows.append(np.nan)
                    rows.append(np.nan)
                print('res:{}\t{:.0f}\t{}\t{:.0f}\t{:.5f}\t{:.5f}\t{:.4f}\t{:.4f}'.format(*rows[-8:]))
        
    if dim_range_from_data:
        dimsuffix='_dim_range_from_data'
    else:
        dimsuffix=''

    nseqsuffix='_'.join([str(x) for x in nseqs])

    df=pd.DataFrame(np.array(rows).reshape(len(inlist)*len(nseqs)*nsamples,len(colnames)), columns=colnames)
    df.to_csv(out_folder+'simu_{}_nseq{}_subsampled_from10000{}.csv'.format(starlike, nseqsuffix, dimsuffix), float_format='%.5f', index=False)

