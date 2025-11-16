#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import sys
from tqdm import tqdm
import os
import re
import pandas as pd
import sys

"""
Function to check that the translations of the nt sequences of gammaproteo eggNOGs
correspond exactly to the aa sequences in gammaproteo eggNOGs. 
First argument is an aa file list, second - nucleotiede file list, 
third - folder name to write the output list of ids and the table with the number 
of non-matching sequences per eggnog (optional)
"""

if __name__ == "__main__":
    #open lists with the names of files to compare
    #(files should be in the same order)
    if sys.argv[1][-5:]=='.list':
        aa=pd.read_csv(sys.argv[1], header=None)
        nt=pd.read_csv(sys.argv[2], header=None)
    else: #inputs are individual files to compare
        aa=[[sys.argv[1]]]
        nt=[[sys.argv[2]]]
        nonlist=True

    idsonly=[]
    
    if len(sys.argv)==4:
        folder=sys.argv[3]
    else:
        folder=''

    oglist=[]
    countlist=[]
    idlist=[]
    matching=[]
    match_index=[]
    #open each file and read sequences
    for ntfile, aafile in tqdm(zip(list(nt[0]), list(aa[0]))):
        aas=list(SeqIO.parse(aafile, 'fasta'))
        nts=list(SeqIO.parse(ntfile, 'fasta'))
        length=len(aas)
        match_index.append(length)
        #get the id/name of the eggnog from the file name 
        og_id=re.findall("^[A-Za-z0-9]*", os.path.split(ntfile)[-1])[0]
        counter=0
        #iterate by every pair of corresponding sequences from two files
        for index, i, k in enumerate(zip(nts, aas)):
            #check that sequences' ids and sequences are identical
            if i.id==k.id and i.seq==k.seq:
                if nonlist:
                    matching.append(i.id)
                    match_index.append(index)
                continue
            else: #record sequence records that differ in id or in sequence
                idlist.append(og_id+'_'+i.id)
                idsonly.append(i.id)
                counter+=1
        oglist.append(og_id)
        countlist.append(counter)
    if nonlist:
        with open(folder+og_id+'_nt_aa_mismatch_ids.list', 'w+') as idso:
            idso.write('\n'.join(idsonly)+'\n')
        with open(folder+og_id+'_nt_aa_match_ids.list', 'w+') as idso:
            idso.write('\n'.join(matching)+'\n')
        with open(folder+og_id+'_nt_aa_match_index.list', 'w+') as idso:
            idso.write('\n'.join(match_index)+'\n')
    
    else:
        #write the list of all sequence ids woth their og_id if their nt sequence translation differs from eggnog aa sequence
        with open(folder+'eggnogs_nt_aa_mismatch_ids.list', 'w+') as idso:
            idso.write('\n'.join(idlist)+'\n')

        #table with the number of non-matching sequences per eggnog
        countdf=pd.DataFrame({'OG_name':oglist, 'nt_aa_mismatch_count':countlist})
        with open(folder+'eggnogs_nt_aa_mismatch_counts.tsv', 'w') as vo:
            countdf.to_csv(vo, sep='\t', index=False)