#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to get the length of sequences with gaps (alignment), the fraction of invariable sites, 
the average and median fraction of gaps and the alphabet size (number of different characters in every column).
"""
from Bio import AlignIO, SeqIO
from tqdm import tqdm
import argparse
import os
import pandas as pd
import numpy as np
import re

def ali_len_nvar_sites_gap_content(file_list, outdir, indir, flist=True):
    if flist==False:
        fl=[file_list]
    else:
        fl=open(file_list, 'r').read().split("\n")[:-1]

    ogs=[]
    seq_lengths=[]
    frac_gaps_across_col_median=[]
    frac_gaps_across_col_mean=[]
    invar_sites=[]
    alpha_size_median=[]
    alpha_size_mean=[]
    for i in tqdm(fl):
        #get OG ID
        og_id=re.findall("^[A-Za-z0-9]*", os.path.split(i)[-1])[0]
        #open alignment
        ali=AlignIO.read(open(indir+i), 'fasta')
        #get the length of sequences with gaps (alignment)
        length =  ali.get_alignment_length()
        #get the number of sequences in the alignment
        anseq=len(ali)
        #count the number of sites with only one possible character
        n_invar_sites=0
        #get the fraction of gaps at every site
        gap_fraction=[]
        #get the alphabet size for every site (including gaps)
        col_alpha_size=[]
        #iterate over every column in the alignment and get mentionned values for every site
        for col in range(length):
            #get a column in the alignment as a string
            currcol=ali[:, col]
            #count what fraction of this column constitute gaps
            gap_fraction.append(currcol.count('-')/anseq)
            #get the list of the characters in this column
            colvars=set(currcol)
            #get the number of characters in this column
            col_alpha_size.append(len(colvars))
            #if a column has the same character in all sequences count it as an invariant site
            if col_alpha_size[-1]==1:
                #this is an invariable site
                n_invar_sites=n_invar_sites+1
        #record average characteristics for the whole alignment
        ogs.append(og_id)
        seq_lengths.append(length)
        frac_gaps_across_col_median.append(np.median(gap_fraction))
        frac_gaps_across_col_mean.append(np.mean(gap_fraction))
        invar_sites.append(n_invar_sites)
        alpha_size_median.append(np.median(col_alpha_size))
        alpha_size_mean.append(np.mean(col_alpha_size))

    df=pd.DataFrame({'OG_name':ogs, 'Length_of_ali':seq_lengths, 'num_invar_sites':invar_sites, 
                  'median_frac_gaps':frac_gaps_across_col_median,
                  'mean_frac_gaps':frac_gaps_across_col_mean, 'median_alphabet_size':alpha_size_median,
                'mean_alphabet_size':alpha_size_mean})
    
    df['median_frac_gaps']=df['median_frac_gaps'].round(decimals=3)
    df['mean_frac_gaps']=df['mean_frac_gaps'].round(decimals=3)
    df['median_alphabet_size']=df['median_alphabet_size'].round(decimals=3)
    df['mean_alphabet_size']=df['mean_alphabet_size'].round(decimals=3)

    #write dataframe with summary stats to file
    with open(outdir+os.path.split(file_list)[-1].replace('.list', '.')+'ali_length_gap_content_var_sites.tsv', 'w+') as of:
        df.to_csv(of, sep='\t', index=False)
                  

if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file_list', help='list of OGs IDs',
                        required=True)
    #Optional
    PARSER.add_argument('-o', '--outdir', help='folder to write the output files (only sequences having nt)',
                        required=False, default="")
    PARSER.add_argument('-i', '--indir', help='folder where the input files are located',
                        required=False, default="")
    
    PARSER.add_argument('-l', '--list_bool', type=eval, choices=[True, False], 
                        help='boolean,  True - if FILE is a list of file names, False - if File is a file name',
                        required=False, default='True')
    
    ARGS = vars(PARSER.parse_args())
                  

    file_list = ARGS['file_list']              
    outdir = ARGS['outdir']
    list_bool = ARGS['list_bool']
    indir = ARGS['indir']

                  
    ali_len_nvar_sites_gap_content(file_list, outdir, indir, flist=list_bool)
