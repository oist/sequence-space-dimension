#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to calculate pairwise distances between sequences from multiple alignments

Inputs:
- FILENAME - protein/ncleotide alignmnet file in fasta format
- IN_FOLDER
- OUT_FOLDER
- REM_DUPLICATES - boolean type, True - remove duplicates sequences, False - keep all sequences
- AA_SITES_ONLY - boolean variable, True - consider only sites with non-gap characters in both sequences for distance calculation, 
False - consider all sites in which there is a non-gap character in at least one of the sequences. Default: False
- SUFFIX - (optional) suffix to add to the output file name (pairwise distances matrix)

Usage:
python3 hamming_distance_multiple_ali_universal.py -a natural_FPs_jun2020_aligned.fasta -if gfp/  -of gfp/
"""

import argparse
from Bio import SeqIO
import pandas as pd
import re
import sys 
import numpy as np
import itertools
from scipy import stats
from sklearn.linear_model import LinearRegression
from scipy.spatial.distance import squareform

def pairwise_distance_from_multiple_ali(filename, remdupl=True, in_folder='', out_folder='', aa_sites_only=False, suffix=''):
    """ 
    This is a function calculate hamming distacnes between all pairs of sequences from the multiple alignment.
    Input file should be an AA alignment file or a multifasta file of aligned AA sequences in fasta format
    Output is a matrix of pairwise distances normalized on the alignment lenght of each pair
    Sequences are NOT realigned. This script is not sensitive to the columns consisting only of gaps.
    """
    
    fn=re.findall('(.+)\.', filename)[0]
    
    with open(in_folder+filename, 'r') as infile:
        aligned_seq=list(SeqIO.parse(infile, 'fasta'))
        ids=[rec.id for rec in aligned_seq]
        sequences=[str(rec.seq).upper() for rec in aligned_seq]
        if remdupl==True: #filtering out duplicates (identical sequences in different species)
            
            filtdf=pd.DataFrame({'id':ids, 'seq':sequences})
            filtdf.drop_duplicates(subset=['seq'], keep='first', inplace=True)
            d=dict(zip(filtdf['id'], filtdf['seq']))
        elif remdupl==False:
            
            d=dict(zip(ids, sequences))
        
    # getting the list of all possible combinations of sequences in this alignments    
    comb=list(itertools.combinations(list(d.keys()), 2))
    distances=[]  
    
    # getting the list (l) of non-gapped sites for each sequence in the alignment
    for name, seq in d.items():
        l=[]
        for ind in range(len(seq)):
            if seq[ind]!='-':
                l.append(ind)
        d[name]=[seq, l]
        
    for i in comb:
        
        if aa_sites_only==False:
            # getting the list of common positions in the alignment for a pair of sequences (all positions where at least one of two sequences has a non-gap character)
            inter=set(d[i[0]][1]+d[i[1]][1])
            
        elif aa_sites_only==True:
            #getting the list of positions where both sequences have a non-gap character
            inter= set(d[i[0]][1]).intersection(d[i[1]][1])
            
        # calculating hamming distance between 2 sequences
        disnorm=np.nan
        distance = 0
        if len(inter)!=0:
            for k in inter:
            # Add 1 to the distance if these two characters are not equal
                if d[i[0]][0][k] != d[i[1]][0][k]:
                    distance += 1
            disnorm=distance/len(inter) # normalize on a alignment length in common
        else: #if there are no common sites, the distance is 1
            disnorm=1
        distances.append(disnorm)
    
    #transform the list of distances into a matrix
    distdf=pd.DataFrame.from_records(np.triu(squareform(distances))).round(5)
    #record matrix to file
    if suffix!='':
        suffix='_'+suffix
    with open('%s%s%s.pdistm' % (out_folder, fn, suffix), 'w+') as outm:
        distdf.to_csv(outm, index=False, header=False)
            
        
if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-a', '--alignment', help='multiple sequence alignment in fasta format',
                        required=True)
    #Required
    PARSER.add_argument('-if', '--input_folder', help='folder where the multiple sequence alignment in fasta format are',
                        required=False, default="")
    PARSER.add_argument('-of', '--output_folder', help='folder where the multiple sequence alignment in fasta format are',
                        required=False, default="")
    #Optional
    
    PARSER.add_argument('-aa', '--aa_sites_only', type=eval, choices=[True, False], help='boolean type, True - only sites with non-gap characters in both seqs, False - all sites with a non-gap character in at least one of the sequences',
                        required=False, default='False')
    PARSER.add_argument('-s', '--suffix', help='suffix to add to the resulting pairwise distance matrix',
                        required=False, default='')
    PARSER.add_argument('-d', '--rem_duplicates', type=eval, choices=[True, False], help='boolean type, True - remove duplicates sequences, False - keep all sequences',
                        required=False, default='True')
    
    ARGS = vars(PARSER.parse_args())

    alignment = ARGS["alignment"]
    aa_sites_only = ARGS['aa_sites_only']
    suffix = ARGS['suffix']
    in_folder = ARGS["input_folder"]
    out_folder=ARGS["output_folder"]
    remdupl = ARGS["rem_duplicates"]

    pairwise_distance_from_multiple_ali(alignment, remdupl=remdupl, in_folder=in_folder, out_folder=out_folder, aa_sites_only= aa_sites_only, suffix=suffix)
