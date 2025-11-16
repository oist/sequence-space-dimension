#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to list to the output file only the files with more than ARGV[2] sequences in them 
Inputs:
- ARGV[1] - a table (txt file) with file names in first column and sequecnes counts in second,  separated by colon ':'
- ARGV[2] - minimum number of sequences a file should contain in order to be kept in the output file 
Outputs:
- tsv (tabulated table called ARGV[1][:-4]_seqs_or_more.txt) with the file names in the first column and the number of sequences in the file in the
second column (table contains only the files with ARGV[2] or more sequences)
Usage:
python3 filter_fasta_by_num_of_seq.py all_gammaproteo_eggNOGs_nt_aligned_num_of_seq.txt 200
"""

import sys
import pandas as pd

lengs = pd.read_csv(sys.argv[1], header=None, sep=":", dtype ={0:str, 1:int})
filt_lengs=lengs[lengs[1]>=int(sys.argv[2])]

filt_lengs=filt_lengs.rename(columns={0:'OG_name', 1:'num_of_seq'})
with open(sys.argv[1][:-4]+'_'+str(sys.argv[2])+'_seqs_or_more.tsv', 'w+') as of:
	filt_lengs.to_csv(of, index=False, sep='\t')
