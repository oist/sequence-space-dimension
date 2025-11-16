#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to remove second (and all subsequent) duplicate sequences from
the pdistm matrices.

Inputs:
- PDISTM

Outputs:
- PDISM_UNIQ

Usage: filter_identical_seqs_from_pdistm.py -m -if ./ -of ./
"""

import argparse
import re
import os
import pandas as pd
import numpy as np

def filter_matr_for_duplicates(matr, infolder='', outfolder=''):
  #reading the matrix from cav file
  pdist=pd.read_csv(infolder+'/'+matr, header=None).to_numpy()
  #get the indices of the columns of the elements in the upper triangular
  #part of the matrix
  _, tcolsi = np.triu_indices(np.shape(pdist)[0], k = 1)
  #get the elements in the upper triangular part of the matrix (all that have values)
  tval=pdist[np.triu_indices(np.shape(pdist)[0], k = 1)]
  #get the indices of the elements in the upper triangular part of the matrix 
  #whose values is zero (they indicate the duplicate sequences we want ro get rid of)
  col_to_remove = tcolsi[tval==0]
  #remove the rows of the all repeating sequences (keep the first sequence
  #and remove all subsequencet)
  pdist=np.delete(pdist, col_to_remove, axis=0)
  #remove the columns of the all repeating sequences (keep the first sequence
  #and remove all subsequencet)
  pdist=np.delete(pdist, col_to_remove, axis=1)
  
  df=pd.DataFrame.from_records(np.array(pdist))
  fn=os.path.splitext(matr)[0]
  with open(outfolder+'/'+fn+'_unique.pdistm', 'w+') as outm:
    df.to_csv(outm, index=False, header=False)

if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-m', '--matrix', help='pairwise distance/dissimilarity matrix name to filter out duplicates from',
                        required=True)
    #Optional
    PARSER.add_argument('-if', '--input_folder', help='folder where the matrix to filter is',
                        required=False, default=".")
    PARSER.add_argument('-of', '--output_folder', help='folder where to write filtered matrix',
                        required=False, default=".")
    
    ARGS = vars(PARSER.parse_args())

    matr = ARGS["matrix"]
    infolder = ARGS["input_folder"]
    outfolder = ARGS["output_folder"]   

    filter_matr_for_duplicates(matr, infolder=infolder, outfolder=outfolder)
