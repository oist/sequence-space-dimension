#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
This is the script to calculate dimension (k coefficient) from pairwise Levenstein (Hamming) distances matrix (*.pdistm file)
Script calculates:
    (a) max_k - the slope of the steepest region of the curve of the size set by WIN_SIZE (10 points - by default) 
    with 95% confidence interval CI_max_k
    CI = ts*standard error on the slope , there ts - values from t-student inverse dictribution corresponding to 95% confidence level and 
    N-1 degrees of freedom (with N - number of rows in the matrix)
    (b) breakpoint - the first point included in region used to calculate a max_k
    
This script uses:
- regular spacing (in log scale) between points on the pairs of points-distance plot for regression
- normalization by the number of distinct points (sequences) in the matrix
  k = log((pairs of points at dist d)/(# points)^2)/log(d)
    
Supports two modes of calculation: 
    - MODE="max" : dimension is the highest slope of the linear regression slopes done on all windows of size WIN_SIZE_PERC*NUMBINS
    - MODE="both": function outputs the "max" dimension value and "0-50%" dimension 
      (the slope of the linear regression done on the first WIN_SIZE_PERC*NUMBINS (0-WIN_SIZE_PERC% divergence))

Input: 
- DIST_FILE - upper triangular/symmetric matrix with pairwise distances between sequences
Example matrix:
0,1,2    0,1,2
0,0,3 or 1,0,3
0,0,0    2,3,0
Output:
*.coefnr file: csv table with header containing k coefficient, its R^2 value (goodness of linear regression approximation) and confidence interval for dimension
coefficient k calculated for a range of pairwise distance between sequences pairs with shighest slope (max_k), as well as the breakpoint (in % of sequence distance)
(default "max" mode)

"""

#importing modules
import argparse
import pandas as pd
import re
import os
import numpy as np
import itertools
from scipy import stats
import ast
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import math
import sys

sys.path.append("/bucket/KondrashovU/seq_space/scripts/seq_space_lib/")
sys.path.append("../scripts/seq_space_lib/")
import sequence_space_lib as seqsp

from sequence_space_lib import t_student_one_sided

if __name__=='__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file', help='upper triangular matrix with pairwise distances',
                        required=True)

    #Optional
    PARSER.add_argument('-i', '--in_folder', help='folder to find the input file', required=False, default="")
    PARSER.add_argument('-of', '--out_folder', help='folder to wtite the output file', required=False, default="")
    PARSER.add_argument('-op', '--out_plots', help='folder to wtite the output plots', required=False, default="")
    
    PARSER.add_argument('-l', '--log', type=eval, choices=[True, False], help='scaling to use for regression - True=log (equal), False=non_log (unequal)', required=False, default="False")
    PARSER.add_argument('-w', '--win_size_perc', type=eval, help="window size to get the max slope as a percentage of NUMBINS, with default numbins=20, +1 == +5% sequence divergence", 
                        required=False, default=0.4)
    PARSER.add_argument('-n', '--numbins', type=eval, help="number of bins to use when binning the pairwise distances", required=False, default=20)
    
    PARSER.add_argument('-p', '--plots',  type=eval, choices=[True, False], help='set to True to plot the regression', required=False, default="True")
    PARSER.add_argument('-s', '--save_plot',  type=eval, choices=[True, False], help='set to True to save the regression plot', required=False, default="True")
    PARSER.add_argument('-a', '--nparray', type=eval, choices=[True, False],  help='set to True if input matrix is in np.array format instead of csv table format', required=False, default='False')
    PARSER.add_argument('-m', '--mode', choices=['max', 'min_max', 'min_max_50'],  help="how to calculate dimension, 'max' - on the range of max slope, 'both' - on the range of max slope and on 0-50 perc distance range", 
                        required=False, default='max')
    PARSER.add_argument('-fl', '--flat', type=eval, choices=[True, False],  help="indicates if the input is a flattened 1d np.array/matrix (usually from create_normal_matrix function - a simulated matrix)", 
                        required=False, default='False')
    PARSER.add_argument('-r', '--dim_range_from_data', type=eval,  help='set to True to use distance range for dimensionality calculation from the data (not 0-1), False by default', 
                        required=False, default=False)

    ARGS = vars(PARSER.parse_args())
    
    file = ARGS["file"]
    in_folder = ARGS["in_folder"]
    out_folder = ARGS["out_folder"]
    out_plots = ARGS["out_plots"]
    log=ARGS["log"]
    win_size_perc = ARGS["win_size_perc"]
    numbins=ARGS['numbins']
    bool_plots = ARGS["plots"]
    save_plot=ARGS['save_plot']
    nparray = ARGS["nparray"]
    mode=ARGS["mode"]
    simu_flat=ARGS["flat"]
    dim_range_from_data=ARGS['dim_range_from_data']

    #get the values from inverse T-Student distribution for confidence level 0.05 (95%) 
    #(precalculated for some regression window sizes)
    #one sided T Student test is performed because our curve is cumulative, therefore slope will be > 0
    if win_size_perc*numbins==10:
        t_stud=1.859548
    elif win_size_perc*numbins==20:
        t_stud=1.734064
    elif win_size_perc*numbins==8:
        t_stud=1.943180
    elif win_size_perc*numbins==16:
        t_stud=1.761310
    elif win_size_perc*numbins==6:
        t_stud=2.131847
    elif win_size_perc*numbins==12:
        t_stud=1.812461
    else:
        t_stud=t_student_one_sided(win_size_perc*numbins)
        
    
    seqsp.max_k(file, log=log, win_size_perc=win_size_perc, numbins=numbins, t_stud=t_stud, mode=mode,plot=bool_plots, in_folder=in_folder, 
          out_folder=out_folder, out_plots=out_plots, nparray=nparray, from_file=True,save_picture=save_plot, flat=simu_flat, label = '', 
          reg_range_from_data=dim_range_from_data) 

