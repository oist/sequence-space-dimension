#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to calculate mean and median sequence length 

Usage:
python3 hamming_distance_pairwise_ali_universal.py -a natural_FPs_jun2020_aligned.fasta -if gfp/ -wd gfp/

"""

import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import re
from tqdm import tqdm 


def get_unique_seq_len(file_list, in_folder='', out_file='', filtlist=None, filtlen=None):
    res=[]
    fl=open(file_list, 'r').read().split("\n")[:-1]
    for filename in fl:
        print(filename)
        fn=re.findall('(.+)\.', filename)[0]
        with open(in_folder+filename, 'r') as infile:
            aligned_seq=list(SeqIO.parse(infile, 'fasta'))
            ids=[rec.id for rec in aligned_seq]
            sequences=[str(rec.seq).replace('-','').upper() for rec in aligned_seq]
            #filtdf=pd.DataFrame({'id':ids, 'seq':sequences})
            #filtdf.drop_duplicates(subset=['seq'], keep='first', inplace=True)
            #sequences_u=list(filtdf['seq'])
        lens=[]
        for i in sequences:
            lens.append(len(i))
        if filtlist:
            if fn in filtlist.split(','):
                lens=[i for i in lens if i <= int(filtlen)]
        res.append([fn, np.mean(lens), np.median(lens)])

    df=pd.DataFrame(res, columns=['OG_id', 'mean_seq_len_orig', 'median_seq_len_orig'])
    df.to_csv(out_file, sep=',', index=False)

if __name__=='__main__':
    if len(sys.argv)>4:
        get_unique_seq_len(sys.argv[1], in_folder=sys.argv[2], out_file=sys.argv[3], filtlist=sys.argv[4], filtlen=sys.argv[5])
    else:
        get_unique_seq_len(sys.argv[1], in_folder=sys.argv[2], out_file=sys.argv[3])
