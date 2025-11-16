import random
import pandas as pd
import itertools
import numpy as np
import dnds
from Bio import SeqIO
import sys
import re
import os

#usage:
# hamming_distance_wo_gaps.py FASTA_WITH_PATH FRACTION_OF_DISTS_TO_SAMPLE OUTPUT_PATH

def fast_sample_unique_pairs(n, k):
    seen = set()
    pairs = []
    while len(pairs) < k:
        a = random.randint(0, n - 1)
        b = random.randint(0, n - 1)
        if a == b:
            continue
        pair = tuple(sorted((a, b)))  # enforce order so (a, b) ≡ (b, a)
        if pair not in seen:
            seen.add(pair)
            pairs.append(pair)
    return pairs

def get_distance_array(sequence_list, fraction=1.0, upper=True):
    """
    Function that calculates pairwise distances between sequences in a list.

    IMPORTANT: pairwise distances calculated here are only valid if the sequences 
                don't have any gaps (all of the same length). Otherwise, a standalone 
                function hamming_distance_multiple_ali_universal.py should be used to get
                a pairwise distance matrix. 
                Pairwise distances here are calculated using dnds.hamming() function, 
                and it's output is divided by the length of the first sequnence in the 
                family (since the lengths are the same). This function doesn't remove gaps.
                
    Input:
    - list of sequeces (strings)
    - fraction of possible unique pairs to sample

    Output:
    - np.ndarray of pairwise distances
    """
    nseqs = len(sequence_list)
    seqlenaa = len(sequence_list[0])
    total_pairs = nseqs * (nseqs - 1) // 2

    if upper:
        seqlist=[i.upper() for i in sequence_list]
        sequence_list=seqlist

    if not (0 < fraction <= 1):
        raise ValueError("fraction must be between 0 and 1")

    if fraction == 1.0:
        # use all possible pairs
        sampled_pairs = list(itertools.combinations(range(nseqs), 2))
    else:
        num_sampled_pairs = max(1, int(total_pairs * fraction))
        #we need this function because we need unique pairs, not unique indices
        sampled_pairs = fast_sample_unique_pairs(nseqs, num_sampled_pairs)

    frac_dist_evolved = np.zeros(len(sampled_pairs))
    for ind, (i, j) in enumerate(sampled_pairs):
        frac_dist_evolved[ind] = dnds.hamming(sequence_list[i], sequence_list[j]) / seqlenaa

    return frac_dist_evolved


fn=re.findall('(.+)\.', os.path.basename(sys.argv[1]))[0]
    
with open(sys.argv[1], 'r') as infile:
    aligned_seq=list(SeqIO.parse(infile, 'fasta'))
    sequences=[str(rec.seq) for rec in aligned_seq]

seqlen=len(sequences[0])
finalseq=[]
for i in sequences:
    if len(i)!=seqlen:
        continue
    else:
        finalseq.append(i)

dist_arr=get_distance_array(finalseq, float(sys.argv[2]))

np.save(sys.argv[3] + fn+'_sampled'+sys.argv[2]+'.npy', dist_arr)
