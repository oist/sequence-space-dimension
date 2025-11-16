#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to calculate pairwise distances between sequences in multiple alignments

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

def align_mafft_wo_file(in_fasta, threads=1, extra_args=None):
    """
    Run MAFFT

    Args:
        in_fasta (str): path to input fasta
        extra_args (list): extra command-line args (list of strings)

    Returns:
        alignment (Bio.Align.MultipleSeqAlignment)
        ali_text (str)  - alignment as FASTA/Clustal text
        log_text (str)  - stderr text
    """
    extra_args = list(extra_args) if extra_args else []

    cmd = ["mafft"]
    cmd += ["--thread", str(threads)] + extra_args + [in_fasta]

    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"mafft failed ({r.returncode}).\n{r.stderr}")

    return(r.stdout)

def pairwise_distance_pairwise_ali(filename, in_folder='', out_folder='', suffix='', 
                                   remdupl=True, outlen=False, alialg='muscle', threads=1, 
                                   write_partial=None, startfrom_pair=0, remgaps=False, filtlong=None):
    """ 
    This is a function to pairwise align (using MUSCLE) and calculate hamming distacnes between all pairs of sequences in the alignment.
    Input file should be an AA alignment file or a multifasta file of non-aligned AA sequences in fasta format (1 seq = 1 line)
    Output is an array and a matrix of pairwise distances normalized on the alignment lenght of each pair
    All pairs of sequences are realigned.

    Args:
        write_partial (None or int): if int write the distances to a file every INT pairs, recommended int=10000 or 5000  
    """
    
    fn=re.findall('(.+)\.', filename)[0]
    
    if os.path.isfile('%s%s%s.pdistm' % (out_folder, fn, suffix)):
        if outlen and os.path.isfile('%s%s%s.plenm' % (out_folder, fn, suffix)):
            #check that both length and distances were successfully saved
            dist_len = len(np.loadtxt('%s%s%s.pdistm' % (out_folder, fn, suffix), delimiter=",").ravel().tolist())
            len_len = len(np.loadtxt('%s%s%s.plenm' % (out_folder, fn, suffix), delimiter=",").ravel().tolist())
            if len_len==dist_len:
                print('%s%s%s.pdistm' % (out_folder, fn, suffix)+' is already calculated')
                return
            else:
                print('\nplenm and pdistm files do not correspond, resetting pdistm file for '+fn)
                os.remove('%s%s%s.pdistm' % (out_folder, fn, suffix))
                print('\n'+fn+'\n')
        else:
            print('%s%s%s.pdistm' % (out_folder, fn, suffix)+' is already calculated')
            return
    else:
        print('\n'+fn+'\n')

    with open(in_folder+filename, 'r') as infile:
        aligned_seq=list(SeqIO.parse(infile, 'fasta'))
        ids=[rec.id for rec in aligned_seq]
        if remgaps: #if removing gaps (realigning sequences from MSA)
            sequences=[str(rec.seq).replace('-','').upper() for rec in aligned_seq]
        else:
            sequences=[str(rec.seq).upper() for rec in aligned_seq]

        filtdf=pd.DataFrame({'id':ids, 'seq':sequences})

        if filtlong: #filter sequnces longer than filtlong
            filtdf=filtdf[filtdf['seq'].str.len()<=int(filtlong)]
            print('%s: Filtered out %d sequences longer than %s\n' % (fn, len(ids)-len(filtdf), str(filtlong)))

        if remdupl==True: #filtering out duplicates (identical sequences in different species)
            filtdf.drop_duplicates(subset=['seq'], keep='first', inplace=True)

        d=dict(zip(filtdf['id'], filtdf['seq']))
        
    num_of_seq=len(d)
    print('%d unique seqeunces\n' % num_of_seq)

    distances=[]  
    lengths=[]

    if write_partial is not None: #if writing the progress file is required
        write_partial=int(write_partial) 
        if startfrom_pair=='from_data': #if the pair to start from should be set based on what's been calculated
            if os.path.isfile('%s%s%s.pdistm.partial' % (out_folder, fn, suffix)):
                #load the values from the progress file
                calculated = np.loadtxt('%s%s%s.pdistm.partial' % (out_folder, fn, suffix), delimiter=",").ravel().tolist()
                count_pairs=len(calculated)
                #make sure there are no broken rows
                if count_pairs%write_partial==0: #full row
                    startfrom_pair=count_pairs
                    print('partial progress detected - start from pair %d\n' % startfrom_pair)
                else:
                    raise ValueError("the number of pairs in pdist is {} - non-divisible by row length".format(count_pairs))
                
                if outlen:
                    if os.path.isfile('%s%s%s.plenm' % (out_folder, fn, suffix)):
                        #load the values from the progress file
                        calculated_len = np.loadtxt('%s%s%s.plenm' % (out_folder, fn, suffix), delimiter=",").ravel().tolist()
                        count_pairs_len=len(calculated_len)
                        if count_pairs!=count_pairs_len: 
                            raise ValueError("the number of pairs in pdistm is different from plenm")
                    else:
                        raise ValueError("the {}.plenm doesn't exist".format(out_folder+fn+suffix))
            else:
                startfrom_pair=0

        elif int(startfrom_pair)==0:
            #reset partial files
            if os.path.isfile('%s%s%s.pdistm.partial' % (out_folder, fn, suffix)):
                print('resetting partial files for '+fn)
                os.remove('%s%s%s.pdistm.partial' % (out_folder, fn, suffix))
            if outlen:
                if os.path.isfile('%s%s%s.plenm' % (out_folder, fn, suffix)):
                    os.remove('%s%s%s.plenm' % (out_folder, fn, suffix))


    comb=list(itertools.combinations(list(d.keys()), 2))[int(startfrom_pair):]

    for ind, i in tqdm(enumerate(comb)):

        with open(out_folder+fn+'.faa', 'w+') as out:
            out.write('\n'.join(['>'+i[0],d[i[0]],'>'+i[1],d[i[1]]])+'\n')
          
        if os.path.isfile(out_folder+fn+'.ali'):
            os.remove(out_folder+fn+'.ali')

        if  alialg =='mafft':
            #uncomment and use this to run mafft the usual way by letting it write the alignment into file and reading it with python
            #os.system('mafft --thread %d %s.faa > %s.ali 2> %s.log' % (int(threads), out_folder+fn, out_folder+fn, out_folder+fn)) 
            #the default way takes the alignment from the mafft subprocess directly and avoids file writing-reading steps
            ali_fa = align_mafft_wo_file(out_folder+fn+'.faa', threads=threads)
            aligned=list(SeqIO.parse(StringIO(ali_fa), 'fasta'))
            os.remove(out_folder+fn+'.ali') #remove the alignment file after reading it so that it doesn't get read again if muscle malfunctions

        elif alialg == 'muscle':
            os.system('muscle -align %s.faa -output %s.ali -quiet 2> %s.log' % (out_folder+fn, out_folder+fn, out_folder+fn)) #run muscle
            aligned=list(SeqIO.parse(out_folder+fn+'.ali', 'fasta'))
            os.remove(out_folder+fn+'.ali') #remove the alignment file after reading it so that it doesn't get read again if muscle malfunctions
        
        distance = 0
        
        length=len(str(aligned[0].seq))
        lengths.append(length)
        
        for k in range(length):
        # Add 1 to the distance if these two characters are not equal
            if str(aligned[0].seq)[k] != str(aligned[1].seq)[k]:
                distance += 1
            disnorm=distance/length # normalize on a alignment length in common
        distances.append(disnorm)

        if write_partial:
            if ind%write_partial==0 and ind!=0:
                with open('%s%s%s.pdistm.partial' % (out_folder, fn, suffix), 'a') as outm:
                    np.savetxt(outm, [distances[-write_partial:]], delimiter=",", fmt='%.5f')
                if outlen:
                    with open('%s%s%s.plenm' % (out_folder, fn, suffix), 'a') as outm:
                        np.savetxt(outm, [lengths[-write_partial:]], delimiter=",", fmt='%.0f')

    #if some distances were calculated previously - read them from the partial progress file 
    # and add to the list so that the final pdistm always has all pairs  
    if int(startfrom_pair)!=0:
        distances_all=calculated+distances
    else:
        distances_all=distances

    #transform the list of distances into a matrix
    distdf=pd.DataFrame.from_records(np.triu(squareform(distances_all))).round(5)

    #record matrix to file
    if suffix!='':
        suffix='_'+suffix

    #check that the length is correct
    dflen=len(distdf)
    if num_of_seq==dflen:
        with open('%s%s%s.pdistm' % (out_folder, fn, suffix), 'w+') as outm:
            distdf.to_csv(outm, index=False, header=False)
        if outlen:
            if int(startfrom_pair)!=0:
                lengths_all=calculated_len+lengths
            else:
                lengths_all=lengths
            lendf=pd.DataFrame.from_records(np.triu(squareform(lengths_all)))
            with open('%s%s%s.plenm' % (out_folder, fn, suffix), 'w+') as outm:
                lendf.to_csv(outm, index=False, header=False)
    else:
        raise ValueError('the length of the matrix %d differs from the number of sequences %d' % (dflen, num_of_seq))
    
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
    PARSER.add_argument('-a', '--alialg', help='alignment algorithm, mafft or muscle', required=False, default='muscle')
    PARSER.add_argument('-t', '--threads', help='number of threads to use for mafft', required=False, default=1)
    PARSER.add_argument('-w', '--write_partial', help='number of pairwise distances to write to file to save partial progress, set to None if not needed', required=False, default=None)
    PARSER.add_argument('-p', '--startfrom_pair', help='pair number to start from', required=False, default='from_data')
    PARSER.add_argument('-g', '--remgaps', type=eval, choices=[True, False], help='set to True to remove gaps from the imput seqeunces', 
                        required=False, default=False)
    PARSER.add_argument('-fl', '--filtlong', help='maximum allowed sequence length, set to None if not needed', required=False, default=None)

    ARGS = vars(PARSER.parse_args())

    fasta = ARGS["fasta"]
    suffix = ARGS['suffix']
    in_folder = ARGS["input_folder"]
    out_folder=ARGS["output_folder"]
    remdupl = ARGS["rem_duplicates"]
    out_len = ARGS['out_len']
    alialg=ARGS['alialg']
    threads=ARGS['threads']
    write_partial=ARGS['write_partial']
    startfrom_pair=ARGS['startfrom_pair']
    remgaps=ARGS['remgaps']
    filtlong=ARGS['filtlong']
    

    pairwise_distance_pairwise_ali(fasta, remdupl=remdupl, in_folder=in_folder, 
                                   out_folder=out_folder, suffix=suffix, outlen=out_len,
                                   alialg=alialg, threads=threads, write_partial=write_partial,
                                   startfrom_pair=startfrom_pair, remgaps=remgaps, filtlong=filtlong)


