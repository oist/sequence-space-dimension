#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to get the length of sequences with gaps (ali_len) and without gaps 
(mean_seq_len, median_seq_len), the fraction of invariable sites (num_invar_sites, num_invar_sites_no_gaps), 
the average and median fraction of gaps and the alphabet size (number of 
different characters in every column).

median_usage_per_site mean_usage_per_site (same as median_alphabet_size	mean_alphabet_size?)
frac_of_pdist_over_dnds frac_of_pdist_over_dnds_unique
num_invar_sites	median_frac_gaps	mean_frac_gaps	median_alphabet_size	mean_alphabet_size (alphabet size without gaps)

USE ON UNIQUE ALIGNEMNTS ONLY!
"""

from Bio import AlignIO
from tqdm import tqdm
import argparse
import os
import pandas as pd
import numpy as np

def p_entropy(a):
    """ Entropy. """
    return (-1) * np.sum(a * np.log(a))

def keff_entropy(a):
    """ Effective number of amino acids (entropy). """
    return np.exp(p_entropy(a))

def keff_homoplasy(a):
    """Effective number of amino acids (homoplasy)."""
    return 1.0/np.sum(np.square((a)))

def usage_ali_seq_len_nvar_sites_gap_content(file_list, outdir, indir, suffix, flist=True):
    if flist==False:
        fl=[file_list]
    else:
        fl=open(file_list, 'r').read().split("\n")[:-1]

    ogs=[]
    ali_lengths=[]
    median_seq_len=[]
    mean_seq_len=[]
    nsites_02gaps=[]
    nsites_05gaps=[]
    alpha_size_median=[]
    alpha_size_mean=[]
    alpha_size_median02=[]
    alpha_size_mean02=[]
    alpha_size_median05=[]
    alpha_size_mean05=[]
    invar_sites=[]
    invar_sites_02gaps=[]
    invar_sites_05gaps=[]
    invar_sites_no_gaps=[]
    frac_gaps_across_col_median=[]
    frac_gaps_across_col_mean=[]
    #k_eff and entropy
    keff_entropy_mean=[]
    keff_homoplasy_mean=[]
    keff_entropy_02gaps_mean=[]
    keff_homoplasy_02gaps_mean=[]
    keff_entropy_05gaps_mean=[]
    keff_homoplasy_05gaps_mean=[]

    for i in tqdm(fl):
        #get OG ID
        og_id=os.path.split(i)[-1].removesuffix(suffix)
        #open alignment
        ali=AlignIO.read(open(indir+i), 'fasta')
        #get the length of sequences with gaps (alignment)
        length =  ali.get_alignment_length()
        #get the number of sequences in the alignment
        anseq=len(ali)
        #get mean and median sequence lengths
        indseqlen=[]
        for sind in range(anseq):
            #get a column in the alignment as a string
            seq=ali[sind, :]
            indseqlen.append(len(str(seq.seq).replace('-','')))
        mean_seq_len.append(np.mean(indseqlen))
        median_seq_len.append(np.median(indseqlen))

        #count the number of sites with only one possible character
        n_invar_sites=0
        n_invar_sites_no_gaps=0
        n_invar_sites_02gaps=0
        n_invar_sites_05gaps=0
        #get the fraction of gaps at every site
        gap_fraction=[]
        #get the alphabet size for every site (including gaps)
        col_alpha_size=[]
        col_alpha_size_02gaps=[]
        col_alpha_size_05gaps=[]
        site_keff_entropy=[]
        site_keff_homoplasy=[]
        site_keff_entropy_02gaps=[]
        site_keff_homoplasy_02gaps=[]
        site_keff_entropy_05gaps=[]
        site_keff_homoplasy_05gaps=[]
        #iterate over every column in the alignment and get mentionned values for every site
        for col in range(length):
            #get a column in the alignment as a string
            currcol=ali[:, col]
            #count what fraction of this column constitute gaps
            gap_fraction.append(currcol.count('-')/anseq)

            #get amino acid frequences per site ignoring gap frequencies)
            ungapped_currcol = currcol[currcol != '-']
            colfreqs=np.unique(currcol, return_counts=True)[1]/len(ungapped_currcol)

            #get entropy and effective number of amino acids
            site_keff_entropy.append(keff_entropy(colfreqs))
            site_keff_homoplasy.append(keff_homoplasy(colfreqs))

            #get the list of the characters in this column
            colvars=list(set(currcol))
            if '-' in colvars:
                gapin=True
                colvars.remove('-')
            else:
                gapin=False

            #get the number of characters in this column (excluding gaps)
            col_alpha_size.append(len(colvars))
            if gap_fraction[-1]<=0.2:
                col_alpha_size_02gaps.append(len(colvars))
                site_keff_entropy_02gaps.append(site_keff_entropy[-1])
                site_keff_homoplasy_02gaps.append(site_keff_homoplasy[-1])
                if col_alpha_size[-1]==1:
                    n_invar_sites_02gaps=n_invar_sites_02gaps+1
            if gap_fraction[-1]<=0.5:
                col_alpha_size_05gaps.append(len(colvars))
                site_keff_entropy_05gaps.append(site_keff_entropy[-1])
                site_keff_homoplasy_05gaps.append(site_keff_homoplasy[-1])
                if col_alpha_size[-1]==1:
                    n_invar_sites_05gaps=n_invar_sites_05gaps+1
            #if a column has the same character in all sequences count it as an invariant site
            if col_alpha_size[-1]==1:
                if gapin==True:
                    #this is an invariable site with respect to amino acids
                    n_invar_sites=n_invar_sites+1
                else:
                    #this is an invariable site with respect to amino acids AND gaps
                    n_invar_sites_no_gaps=n_invar_sites_no_gaps+1
                    n_invar_sites=n_invar_sites+1
    
        #record average characteristics for the whole alignment
        ogs.append(og_id)
        ali_lengths.append(length)
        frac_gaps_across_col_median.append(np.median(gap_fraction))
        frac_gaps_across_col_mean.append(np.mean(gap_fraction))
        invar_sites.append(n_invar_sites)
        invar_sites_no_gaps.append(n_invar_sites_no_gaps)
        invar_sites_02gaps.append(n_invar_sites_02gaps)
        invar_sites_05gaps.append(n_invar_sites_05gaps)
        nsites_02gaps.append(len(col_alpha_size_02gaps))
        nsites_05gaps.append(len(col_alpha_size_05gaps))
        alpha_size_median.append(np.median(col_alpha_size))
        alpha_size_mean.append(np.mean(col_alpha_size))
        alpha_size_mean02.append(np.mean(col_alpha_size_02gaps))
        alpha_size_median02.append(np.median(col_alpha_size_02gaps))
        alpha_size_mean05.append(np.mean(col_alpha_size_05gaps))
        alpha_size_median05.append(np.median(col_alpha_size_05gaps))
        keff_entropy_mean.append(np.mean(site_keff_entropy))
        keff_homoplasy_mean.append(np.mean(site_keff_homoplasy))
        keff_entropy_02gaps_mean.append(np.mean(site_keff_entropy_02gaps))
        keff_homoplasy_02gaps_mean.append(np.mean(site_keff_homoplasy_02gaps))
        keff_entropy_05gaps_mean.append(np.mean(site_keff_entropy_05gaps))
        keff_homoplasy_05gaps_mean.append(np.mean(site_keff_homoplasy_05gaps))

    df=pd.DataFrame({'OG_id':ogs, 'ali_len':ali_lengths, 'mean_seq_len':mean_seq_len, 'median_seq_len':median_seq_len,
                      'num_invar_sites':invar_sites, 'num_invar_sites_no_gaps':invar_sites_no_gaps,
                      'num_invar_sites_less02gaps':invar_sites_02gaps, 'num_invar_sites_less05gaps':invar_sites_05gaps,
                      'num_sites_less02gaps':nsites_02gaps, 'num_sites_less05gaps':nsites_05gaps,
                  'median_frac_gaps_per_site':frac_gaps_across_col_median,
                  'mean_frac_gaps_per_site':frac_gaps_across_col_mean, 'median_alphabet_size':alpha_size_median,
                'mean_alphabet_size':alpha_size_mean, 'median_alphabet_size_less02gaps':alpha_size_median02,
                'mean_alphabet_size_less02gaps':alpha_size_mean02, 'median_alphabet_size_less05gaps':alpha_size_median05,
                'mean_alphabet_size_less05gaps':alpha_size_mean05, 
                'mean_keff_entropy':keff_entropy_mean, 'mean_keff_homoplasy':keff_homoplasy_mean,
                'mean_keff_entropy_less02gaps':keff_entropy_02gaps_mean, 'mean_keff_homoplasy_less02gaps':keff_homoplasy_02gaps_mean,
                'mean_keff_entropy_less05gaps':keff_entropy_05gaps_mean, 'mean_keff_homoplasy_less05gaps':keff_homoplasy_05gaps_mean})
    
    df['num_var_sites']=df['ali_len']-df['num_invar_sites']
    df['num_var_sites_less02gaps']=df['num_sites_less02gaps']-df['num_invar_sites_less02gaps']
    df['num_var_sites_less05gaps']=df['num_sites_less05gaps']-df['num_invar_sites_less05gaps']
    df=df.round(3)

    #write dataframe with summary stats to file
    with open(outdir+os.path.split(file_list)[-1].replace('.list', '_')+'ali_seq_length_alpha_size_var_invar_sites_keff.tsv', 'w+') as of:
        df.to_csv(of, sep='\t', index=False)
                  

if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file_list', help='list of alignments files (unique seqs only!)',
                        required=True)
    PARSER.add_argument('-s', '--suffix', help='file name after OG_id',required=True)
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
    suffix = ARGS['suffix']

                  
    usage_ali_seq_len_nvar_sites_gap_content(file_list, outdir, indir, suffix, flist=list_bool)
