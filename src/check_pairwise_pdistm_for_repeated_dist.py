#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to reset incorrect matrices from pairwise alignment after muscle malfunction

Usage:
python3 hamming_distance_pairwise_ali_universal.py -a natural_FPs_jun2020_aligned.fasta -if gfp/ -wd gfp/

"""

import argparse
from Bio import SeqIO
import itertools
import os
import pandas as pd
import numpy as np
import re
from tqdm import tqdm 
import sys
from scipy.spatial.distance import squareform
import subprocess
from io import StringIO

def upper_tri_index(x, y, n):
    if x > y:
        x, y = y, x
    return (x * (2*n - x - 1)) // 2 + (y - x - 1)

def equal_tail_length(lst):
    if not lst:
        return 0
    last = lst[-1]
    k = 1
    for v in reversed(lst[:-1]):
        if v == last:
            k += 1
        else:
            break
    return k


def equal_tail_length_np(arr):
    if arr.size == 0:
        return 0
    last = arr[-1]
    # Compare reversed array to last value
    matches = np.flip(arr == last)
    # Find first False → break point
    k = np.argmax(~matches)
    return arr.size if k == 0 else k

def check_matrix_for_repeated_dist(filename, in_folder='', out_folder='', suffix='', 
                                   remdupl=True, outlen=False, write_partial=None):
    """ 
    This is a function to check if the matrix doesn't contain erroneusly repeated 
    the same pairwise distance multiple times because of the alignment software malfunciton

    Args:
        write_partial (None or int): if int write the distances to a file every INT pairs, recommended int=10000 or 5000  
    """
    
    fn=re.findall('(.+)\.', filename)[0]
    
    if os.path.isfile('%s%s%s.pdistm' % (out_folder, fn, suffix)): 
        #check if the ids of the last aligned seqs are the last ids pair of the unique sequences combinations
        #open original fasta with sequences
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

        comb=list(itertools.combinations(list(d.keys()), 2))
        num_of_seq=len(d)
        #open last aligned fasta 
        aligned=list(SeqIO.parse(out_folder+fn+'.ali', 'fasta'))
        if str(aligned[0].id) in comb[-1] and str(aligned[1].id) in comb[-1]:
            print('\n'+fn+' is fine')
            return
        else:
            os.remove('%s%s%s.pdistm' % (out_folder, fn, suffix))
            index_last_aligned_0=list(d.keys()).index(str(aligned[0].id))
            index_last_aligned_1=list(d.keys()).index(str(aligned[0].id))
            index_last_pair=upper_tri_index(index_last_aligned_0, index_last_aligned_1, num_of_seq)
            calculated = np.loadtxt('%s%s%s.pdistm.partial' % (out_folder, fn, suffix), delimiter=",").ravel()

    elif os.path.isfile('%s%s%s.pdistm.partial' % (out_folder, fn, suffix)):
        #check if there is a stretch of identical entries
        #load the values from the progress file
        calculated = np.loadtxt('%s%s%s.pdistm.partial' % (out_folder, fn, suffix), delimiter=",").ravel()
        #check for repeating last entries
        num_repeats=equal_tail_length_np(calculated)
        if num_repeats==1:
            print('\n'+fn+' is fine')
            return 
        else:
            print(num_repeats, len(calculated))
            index_last_pair=len(calculated)-num_repeats

    else:
        print('\n'+fn+' not started')
        return

    print('\nFixing '+fn+', last calculated pair:'+str(index_last_pair))

    #reset all partial files (*.pdistm.partial and *.plenm) to the last aligned pair 
    nrows=index_last_pair//write_partial
    os.rename('%s%s%s.pdistm.partial' % (out_folder, fn, suffix), '%s%s%s.pdistm.partial' % (out_folder+'partial_broken_repeated/', fn, suffix))
    calc_write=calculated[:nrows*write_partial]
    with open('%s%s%s.pdistm.partial' % (out_folder, fn, suffix), 'a') as outm:
        arr = calc_write.reshape(nrows, write_partial)
        np.savetxt(outm, arr, delimiter=",", fmt="%.5f")
        #for i in range(nrows):
        #    np.savetxt(outm, [calculated[write_partial*i:write_partial*(i+1)]], delimiter=",", fmt='%.5f')

    if outlen and os.path.isfile('%s%s%s.plenm' % (out_folder, fn, suffix)):
        lengths=np.loadtxt('%s%s%s.plenm' % (out_folder, fn, suffix), delimiter=",").ravel()
        lengths_write=lengths[:nrows*write_partial]
        os.rename('%s%s%s.plenm' % (out_folder, fn, suffix), '%s%s%s.plenm' % (out_folder+'partial_broken_repeated/', fn, suffix))
        with open('%s%s%s.plenm' % (out_folder, fn, suffix), 'a') as outl:
            arrl = lengths_write.reshape(nrows, write_partial)
            np.savetxt(outl, arrl, delimiter=",", fmt="%.5f")
            #for i in range(nrows):
            #    np.savetxt(outm, [lengths[write_partial*i:write_partial*(i+1)]], delimiter=",", fmt='%.0f')

    
if __name__=='__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--fasta', help='sequences in fasta format',
                        required=True)
    PARSER.add_argument('-if', '--input_folder', help='folder where the multiple sequence alignment in fasta format is',
                        required=False, default="")
    PARSER.add_argument('-of', '--output_folder', help='folder to write pairwise distance matrices outputs',
                        required=False, default="")
    #Optional

    PARSER.add_argument('-s', '--suffix', help='suffix to add to the resulting pairwise distance matrix',
                        required=False, default='')
    PARSER.add_argument('-d', '--rem_duplicates', type=eval, choices=[True, False], help='boolean type, True - remove duplicates sequences, False - keep all sequences',
                        required=False, default='True')
    #Optional
    PARSER.add_argument('-l', '--out_len',  type=eval, choices=[True, False], help='set to True to record pairwise length matrices', 
                        required=False, default=False)
    PARSER.add_argument('-w', '--write_partial', type=int, help='number of pairwise distances to write to file to save partial progress, set to None if not needed', required=False, default=None)

   
    ARGS = vars(PARSER.parse_args())

    fasta = ARGS["fasta"]
    suffix = ARGS['suffix']
    in_folder = ARGS["input_folder"]
    out_folder=ARGS["output_folder"]
    remdupl = ARGS["rem_duplicates"]
    out_len = ARGS['out_len']
    write_partial=ARGS['write_partial']


    check_matrix_for_repeated_dist(fasta, remdupl=remdupl, in_folder=in_folder, 
                                   out_folder=out_folder, suffix=suffix, outlen=out_len,
                                    write_partial=write_partial)


