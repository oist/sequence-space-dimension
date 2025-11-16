#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to retreive from aa alignments only the sequences that are present in nt alignments
and remove all the gaps.
"""
from Bio import SeqIO
from tqdm import tqdm
import argparse
import pandas as pd

def get_nt_corresponding_to_aa(file_list, outdir, aa_folder_suf= '1236_raw_algs/1236/*.raw_alg.faa', nt_folder_suf='eggNOGs_nt_aligned/*.fasta'):
    fl=open(file_list, 'r').read().split("\n")[:-1]
    for i in tqdm(fl):
        #read both nt and aa files
        aa=list(SeqIO.parse(aa_folder_suf.replace('*', i), 'fasta'))
        nt=list(SeqIO.parse(nt_folder_suf.replace('*', i), 'fasta'))
        
        #get seqs ids from nt files
        ntids=[]
        for ntrec in nt:
            ntids.append(ntrec.id)
        
        #retreive corresponding seqs from aa files
        aa_cor_nt=[]
        for aarec in aa:
            if aarec.id in ntids:
                aarec.seq=aarec.seq.ungap("-")
                aa_cor_nt.append(aarec)
        
        #write aa seqs corresponding to nt seqs to output file
        with open(outdir+i+'.aa_corr_to_nt.faa', 'w+') as of:
            SeqIO.write(aa_cor_nt, of, "fasta")
                  
                  
if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file_list', help='list of OGs IDs',
                        required=True)
    PARSER.add_argument('-o', '--outdir', help='folder to write the output files (only sequences having nt)',
                        required=False, default="")
    #Optional
    
    PARSER.add_argument('-a', '--aa_folder_suf', help='folder and the file format/ending of the aa file with *  in the place of OG ID',
                        required=False, default='1236_raw_algs/1236/*.raw_alg.faa')
    PARSER.add_argument('-n', '--nt_folder_suf', help='folder and the file format/ending of the nt file with *  in the place of OG ID',
                        required=False, default='eggNOGs_nt_aligned/*.fasta')
    
    ARGS = vars(PARSER.parse_args())
                  

    file_list = ARGS['file_list']              
    outdir = ARGS['outdir']
    aa_folder_suf = ARGS['aa_folder_suf']
    nt_folder_suf = ARGS['nt_folder_suf']              
                  
    get_nt_corresponding_to_aa(file_list, outdir, aa_folder_suf=aa_folder_suf , nt_folder_suf=nt_folder_suf)
