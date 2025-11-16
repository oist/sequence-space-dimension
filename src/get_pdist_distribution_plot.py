#importing modules
import argparse
import pandas as pd
import pickle
import re
import numpy as np
import ast
import random
import math 
from tqdm import tqdm
import matplotlib.pyplot as plt

#defining functions

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier

def matrix_to_array(matrix, nparray=False, from_file = False):
    """
    Function to convert the upper triangular matrix into a 1-dimensional array with matrix values (from upper half of the input matrix) 
    Possible input file format examples:
    when nparray==False:
    Square matrix in csv format (either upper triangular or symmetric but NOT lower triangular):
    0,1,2
    0,0,3
    0,0,0
    
    when nparray==True:
    Square matrix in np.array format (either upper triangular or symmetric but NOT lower triangular):
    array([[0, 1, 2],
           [1, 0, 3],
           [2, 3, 0]]))
    """
    #first, we open an input matrix in n-dimensional numpy array format, for example:
    #    array([[0, 1, 2],
    #           [0, 0, 3],
    #           [0, 0, 0]]))
    
    if from_file==True:
        if nparray==True: #if input file is a np.array in a txt file (nparray variable should be set to True) 
            pdist=np.array(ast.literal_eval(open(matrix).read()))
        else: #if input file is a csv
            pdist=pd.read_csv(matrix, header=None).to_numpy()
            return(list(pdist[np.triu_indices(np.shape(pdist)[0], k = 1)]), len(pdist))

    #return only the elements above diagonal with np.triu_indices function as a 1d array given the matrix dimension (number of rows/columns)    
    else:
        return(list(matrix[np.triu_indices(np.shape(matrix)[0], k = 1)]), len(matrix))
      
        
def sample_plot_pdist(file, out_folder='', many=False, fraqplot=1):
    """
    Function to draw a histogram of pairwise distances in
    - a file (matrix) -> use MANY=False flag, MATRIX argument then refers to the matrix file name 
    - a set of files (matrices) -> use MANY=True flag, MATRIX argument then refers to the txt file with the list of matrix file names
    """
    
    fn=re.findall('(.+)\.', file.split('/')[-1])[0]
    dists=[]
    if many == True: #input is a list of files
        with open(file) as inp:
            files = inp.read().splitlines()
        for matr in tqdm(files):
            arrm=tuple(matrix_to_array(matr, from_file=True))[0]
            if fraqplot!=1:
                fraqkeep=int(round_half_up(fraqplot*len(arrm)))  #get the number of pairwise distances to keep
                arr_to_plot=random.sample(arrm, fraqkeep)
            else: #plot all the values in the array
                arr_to_plot=list(arrm)
            dists.extend(arr_to_plot)
                    
    elif many == False: #input is one matrix
        arrm=tuple(matrix_to_array(file, from_file=True))[0]
        if fraqplot!=1:
            fraqkeep=int(round_half_up(fraqplot*len(arrm))) #get the number of pairwise distances to keep
            dists=random.sample(arrm, fraqkeep)
        else: #plot all the values in the array
            dists=list(arrm)  
    #plotting
    print('plotting now')
    counts_bins=plt.hist(dists, bins=50)
    with open(out_folder+fn+'_'+str(int(fraqplot*100))+"_pairs_hist.pickle", 'wb') as histbin:
        pickle.dump(counts_bins, histbin)
    plt.xlim(-0.04,1)
    plt.ylabel('pairs of sequences')
    plt.xlabel('pairwise distance')
    plt.suptitle(fn+', '+str(int(fraqplot*100))+'% of all pairs')
    plt.savefig(out_folder+fn+'_'+str(int(fraqplot*100))+"_pairs.png", dpi=200)
    
    
if __name__ == '__main__':	

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file', help='input file, pdistm matrix or list of matrices',
                        required=True)
    #Required
    PARSER.add_argument('-of', '--output_folder', help='folder to save the plots',
                        required=False, default="")
    #Optional
    
    PARSER.add_argument('-m', '--many', type=eval, choices=[True, False], help='boolean type, True - input file is a list of matrices, False - input file is a matrix',
                        required=False, default='False')
    PARSER.add_argument('-fr', '--fraqplot', help='fraction of pairwise distances to plot, default = 1 (all)',
                        required=False, default=1)
    
    ARGS = vars(PARSER.parse_args())

    file = ARGS["file"]
    out_folder=ARGS["output_folder"]
    many = ARGS["many"]
    fraqplot= ARGS["fraqplot"]
    
    sample_plot_pdist(file, out_folder=out_folder, many=many, fraqplot=float(fraqplot))
    
    
