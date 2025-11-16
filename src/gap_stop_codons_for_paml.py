#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to replace with gaps all instances of stop codons in frame and frameshift signs ("!") in nucleotide alignments
Inputs:
- FILELIST - list of the input files
- INFOLDER - directory with the input alignment files
- OUTFOLDER - directory for output alignmnets
- PARALEL - set to False if you want to run all files consequentially on one node

Usage:
python3 gap_stop_codons_for_paml.py -f dnds_NT_ali_list.txt -if alignments/ -of dnds/NT_alignments_no_stop_codons/ 

"""

from Bio import AlignIO
import argparse 
import re
from tqdm import tqdm

def replacer(s, newstring, index):
    """
    Function to replace part of the string
    Inserted string replaces the substring of the same length starting at the given index
    """
    
    # raise an error if index is outside of the string
    if index not in range(len(s)):
        raise ValueError("index outside given string")
        
    # insert the new string between "slices" of the original
    return s[:index] + newstring + s[index + len(newstring):]

def gap_stop_codons(filename, indir, outdir):
    """
    A function to gap all stop codons ('TGA', 'TAA', 'TAG') in the nucleotide sequences alignments (for PAML).
    Stop codon signs (*) in the corresponding amino acid alignmnets remain.
    """
    with open(indir+filename, 'r') as infile:
        s=AlignIO.read(infile, 'fasta')
    l=s.get_alignment_length()

    assert l%3==0, 'sequence length is not a multiple of 3'
    
    stop_codons=['TGA', 'TAA', 'TAG']

    #iterating by each sequence in the alignment
    newali=[]
    for record in s:
        seq=str(record.seq).upper()
        
        #replace all "!" (franeshift signs) with gaps
        seq=seq.replace('!', '-')
        
        #iterating by codons
        for n in range(int(l/3)): #n is a codon number 
            codon=seq[3*n:3*n+3] #current codon
            
            #checing if the codon is in the list of stop codons
            if codon in stop_codons:
                if codon== 'TGA':
                    seq = replacer(seq, '---', 3*n)
                elif codon== 'TAA':
                    seq = replacer(seq, '---', 3*n)
                elif codon== 'TAG':
                    seq = replacer(seq, '---', 3*n)
            else:
                continue
        newali.append('>'+str(record.id)+'\n'+seq+'\n')
    
    fformat = re.findall('(\.[^\.]+)$',filename)[0]
    with open(outdir+filename.replace(fformat, '_nostop'+fformat), 'w+') as outfile:
        outfile.write(''.join(newali))
        
if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file', help='filename (if paralel=True) or list of the names of the fasta files with the nucleotide coding sequence alignments (if paralel=False)',
                        required=True)
    PARSER.add_argument('-if', '--infolder', help='folder where the multiple sequence alignment in fasta format are',
                        required=False, default="")
    PARSER.add_argument('-of', '--outfolder', help='folder where the output alignmnets with gapped stop codons should be placed',
                        required=False, default="")
    #Optional
    PARSER.add_argument('-p', '--paralel', help='True - file is a file name, False - file is a list of file names',
                        required=False, default=True)
   
    ARGS = vars(PARSER.parse_args())

    file = ARGS["file"]
    infolder = ARGS["infolder"]
    outfolder=ARGS["outfolder"]
    paralel=ARGS['paralel']

    if paralel==True:
        gap_stop_codons(file, infolder, outfolder)
    else:
        with open(file) as flist:
            files = flist.readlines()

        for ffile in tqdm(files):
            gap_stop_codons(ffile[:-1], infolder, outfolder)
