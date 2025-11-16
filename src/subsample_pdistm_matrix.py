#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import random
import numpy as np
import pandas as pd
import argparse

def submatrices(matrix_f,n,x): 
    """
    Inputs: 
    matrix – initial matrix as a csv file
    n – size of the matrices we get
    x – number of subsamples
    """
    matrix= pd.read_csv(matrix_f,header = None).to_numpy()
    subsamples=[]
    for j in range(0,x):
        #Makes a random list of numbers within the interval (0, len(matrix)) so the numbers in list do not repeat.
        #We take k = ((len(matrix)) - (len of matrix we are getting)) of those numbers
        to_remove = random.sample(list(range(0,len(matrix))), k = len(matrix)-n)
        #Cuts of the particular columns and rows from the initial matrix Number of rows are from the list avobe.
        subsamples.append(np.delete(np.delete(matrix,to_remove,0),to_remove,1))
    return(subsamples)
  
def multiple_subsamples(matrix, sizes=[5000, 6000, 7000, 8000], reps=5):
    for size in sizes:
        sub_list=submatrices(matrix,size,reps)
        for ind, subm in enumerate(sub_list):
            with open(matrix[:-7]+'_'+str(size)+'_'+str(ind)+'.pdistm', 'w+') as outs:
                np.savetxt(outs, subm, delimiter=",", fmt='%.7f')
              
if __name__ == '__main__':	
    
    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-m', '--matrix', help='pairwise distance/dissimilarity matrix name to subsample',
                        required=True)
    #Optional
    PARSER.add_argument('-s', '--sizes', type=eval, help='a list of sizes of the subsampled matrices',
                        required=False, default="[100]")
    PARSER.add_argument('-r', '--reps',  type=eval, help='number of repetitions of subsampling for each size',
                        required=False, default="5")
    
    ARGS = vars(PARSER.parse_args())

    matrix = ARGS["matrix"]
    sizes = ARGS["sizes"]
    reps = ARGS["reps"]   
    
    multiple_subsamples(matrix, sizes=sizes, reps=reps)
