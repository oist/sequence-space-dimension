#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This is a main sequence space dimensionality project library file
#
# Function sections:
# - generic functions
# - MDS plotting 
# - scatterplots with density shown by color
# - dimensionality 
# - effective topological dimension
# - simulate and subsample matrices
# - reading yn00/PAML output table with dS, dN and dN/dS estimations
#
# How to use this library:
# import sys
# add the path to the sequence_space_lib.py and  __init__.py files to your path 
# sys.path.append("scripts/seq_space_lib/")
# import like normal library
# import sequence_space_lib as seqsp 

import pandas as pd
import numpy as np
import os
import re
import ast
import random
import math
import matplotlib.pyplot as plt
import itertools
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO 
#sklearn
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.linear_model import LinearRegression
#scipy
from scipy import stats
from scipy.spatial.distance import squareform, pdist
import scipy.spatial.distance as scipy_dist
from scipy.special import comb

#generic functions
def round_half_up(n, decimals=0):
    """
    Function whic rounds a float number half up 
    Source: get_pdist_distribution_plots_test.ipynb
    """
    multiplier=10**decimals
    return(math.floor(n*multiplier + 0.5) / multiplier)

def matrix_to_array(matrix, nparray=False, in_folder='', from_file = False):
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
           
    Source: get_pdist_distribution_plots_test.ipynb
    """
    #first, we open an input matrix in n-dimensional numpy array format, for example:
    #    array([[0, 1, 2],
    #           [0, 0, 3],
    #           [0, 0, 0]]))
    
    if from_file==True:
        if nparray==True: #if input file is a np.array in a txt file (nparray variable should be set to True) 
            pdist=np.array(ast.literal_eval(open(in_folder+matrix).read()))
        else: #if input file is a csv
            pdist=pd.read_csv(in_folder+matrix, header=None).to_numpy()
            return(list(pdist[np.triu_indices(np.shape(pdist)[0], k = 1)]), len(pdist))

    #return only the elements above diagonal with np.triu_indices function as a 1d array given the matrix dimension (number of rows/columns)    
    else:
        return(list(matrix[np.triu_indices(np.shape(matrix)[0], k = 1)]), len(matrix))


def t_student_one_sided(sample_size, p=0.05):
    # One-sided inverse Students t-distribution
    # p - probability, df - degrees of freedom
    tinv = lambda p, df: abs(stats.t.ppf(p, df))
    return(tinv(p, sample_size-2))

# MDS plotting 

def MDS_graph(matrix, pcx='PC1', pcy='PC2', n_components=2,max_range=10, 
              stress=False, axs=[], n=None, save=False, path='', k='', metric=True,
              og_name='', title=None, dot_colors=[], alpha=0.5, size=7, noplot=False):
    """
    Function to perform multidimensional scaling (MDS) on pairwise distance matrix
    Source: seq_space_visualization+phylip.ipynb
    """
    
    if type(matrix)==str: #if the input file is a triangular distance matrix in a csv file
        data = pd.read_csv(matrix, header=None)
        data = data + data.T
        og=re.findall('\d+', matrix)[0]
    elif type(matrix)==np.ndarray:
        data=matrix
        og=og_name
    else:
        raise  TypeError('Unkown type of input file')
        
    length=str(len(data))
    
    if stress==False:
        #visualize in 2D without stress plot
        if n_components==2:
            mds = MDS(metric=metric, n_components=n_components, dissimilarity='precomputed', random_state=0)
            transform=mds.fit_transform(data)
            if noplot:
                return(transform)
            if dot_colors!=[]:
                plt.scatter(transform[:,0],  transform[:,1], alpha=alpha, c=dot_colors, s=size)
            else:
                plt.scatter(transform[:,0],  transform[:,1], alpha=alpha, s=size)
            plt.ylabel(pcx)
            plt.xlabel(pcy)
            if title==None:
                plt.title('NCBI %s, pairwise ali, n=%s, k=%s' % (og,length, str(k)))
            else:
                plt.title(title)
            if save==True:
                plt.savefig(path+'NCBI_%s_pairwise_ali_projection.png' % og,  dpi=300)
                
        #visualize in 3D       
        if n_components==3:
            mds = MDS(metric=metric, n_components=n_components, dissimilarity='precomputed', random_state=0)
            transform=mds.fit_transform(data)
            if noplot:
                return(transform)
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(transform[:,0],  transform[:,1], transform[:,2], s=size)
    
    #visualize in 2D with stress plot        
    elif stress==True:
        #calculate stress graph 
        stress = []
        #create list to write point coordinates in different number of dimensions
        transform=[]
        # Max value for n_components
        for dim in range(1, max_range):
            # Set up the MDS object
            mds = MDS(n_components=dim, dissimilarity='precomputed', random_state=0)
            transform.append(mds.fit_transform(data))
            # Retrieve the stress value
            stress.append(mds.stress_)
        # Plot 
        if axs==[]:
            fig, ax=plt.subplots(1,2, figsize=(12,5))
            # Plot MDS scatter
            #plot MDS in 2D (second entry with index 1 in transform list)
            ax[0].scatter(transform[1][:,0],  transform[1][:,1], alpha=alpha, s=size)
            ax[0].set_ylabel(pcx)
            ax[0].set_xlabel(pcy)
            # Plot stress vs. n_components
            ax[1].plot(range(1, max_range), stress)
            ax[1].set_xticks(range(1, max_range, 2))
            ax[1].set_xlabel('n_components')
            ax[1].set_ylabel('stress')

            plt.suptitle('NCBI %s, pairwise ali, n=%s, k=%s' % (og,length, str(k)))

            if save==True:
                plt.savefig(path+'NCBI_%s_pairwise_ali_projection_stress.png' % og,  dpi=300)
        else:
            #plot MDS in 2D (second entry with index 1 in transform list)
            axs[n, 0].scatter(transform[1][:,0],  transform[1][:,1], alpha=alpha, s=size)
            #axs[n, 0].set_ylabel(pcx)
            #axs[n, 0].set_xlabel(pcy)
            axs[n, 0].set_xlim([-0.6, 0.6])
            axs[n, 0].set_ylim([-0.7, 0.65])
            # Plot stress vs. n_components
            axs[n, 1].plot(range(1, max_range), stress)
            axs[n, 1].set_xticks(range(1, max_range, 2))
            axs[n, 1].set_xlabel('n_components')
            axs[n, 1].set_ylabel('stress')
            
            axs[n, 0].set_title('NCBI %s, pairwise ali, n=%s, k=%s' % (og,length, str(k)))
            
    return(data, transform)

#plot things

def sample_plot_pdist(file, in_folder='', out_folder='', many=False, fraqplot=1):
    """
    Function to draw a histogram of pairwise distances in
    - a file (matrix) -> use MANY=False flag, MATRIX argument then refers to the matrix file name 
    - a set of files (matrices) -> use MANY=True flag, MATRIX argument then refers to the txt file with the list of matrix file names
    
    Source: get_pdist_distribution_plots_test.ipynb
    """
    
    fn=re.findall('(.+)\.', file)[0]
    dists=[]
    if many == True: #input is a list of files
        with open(file) as inp:
            files = inp.read().splitlines()
        for matr in files:
            arrm=tuple(matrix_to_array(matr, in_folder=in_folder, from_file=True))[0]
            if fraqplot!=1:
                fraqkeep=int(round_half_up(fraqplot*len(arrm)))  #get the number of pairwise distances to keep
                arr_to_plot=random.sample(arrm, fraqkeep)
            else: #plot all the values in the array
                arr_to_plot=list(arrm)
            dists.extend(arr_to_plot)
                    
    elif many == False: #input is one matrix
        arrm=tuple(matrix_to_array(file, in_folder=in_folder, from_file=True))[0]
        if fraqplot!=1:
            fraqkeep=int(round_half_up(fraqplot*len(arrm))) #get the number of pairwise distances to keep
            dists=random.sample(arrm, fraqkeep)
        else: #plot all the values in the array
            dists=list(arrm)  
    #plotting
    print('plotting now')
    plt.hist(dists, bins=50)
    plt.ylabel('pairs of sequences')
    plt.xlabel('pairwise distance')
    plt.suptitle(fn+', '+str(int(fraqplot*100))+'% of all pairs')
    plt.savefig(out_folder+fn+'_'+str(int(fraqplot*100))+"_pairs.png", dpi=200)


def plot_density(x, y, xlabel, ylabel, xscale='linear', yscale='linear', xlim=None, ylim=None, 
                 big=False, regline=True, size=0.4, density=True, dropna=False, inpdf=False, ax=None, 
                 regcolor='r', palette='viridis', legend=''):
    
    """
    Make scatterplots with point density shown by color
    """
    if not ax:
        fig, ax=plt.subplots()

    if dropna==True: #x is dataframe and xlabel and ylabel are column names
        dropdf=x.dropna(subset=[xlabel, ylabel])
        x=dropdf[xlabel]
        y=dropdf[ylabel]
    elif type(y)==pd.core.series.Series or type(y)==list:
        pass
    elif inpdf or y=='' or y==None:
        x=x[xlabel]
        y=x[ylabel]

    if big== True:
        fig.figure(figsize=(6, 4), dpi= 200)

    if density==True:
        xy = np.vstack([x,y])
        z = stats.gaussian_kde(xy)(xy)
        ax.scatter(x,y, s=size, c=z, cmap=palette)
    else:
        ax.scatter(x,y, s=size, alpha=0.5, cmap=palette)

    if regline==True:
        m, b, r_value, p_value, std_err = stats.linregress(x,y)   
        if b >= 0:
            sign = "+"
        else:
            sign = "-"
        abs_b = abs(b)
        ax.plot(x, m*x + b,label=legend+"y={:.2f}x".format(m)+sign+"{:.2f}\nR^2={:.4f}".format(abs_b,r_value**2), c=regcolor)
        ax.legend()
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if ylim:
        ax.set_ylim(ylim[0], ylim[1])
    if xlim:
        ax.set_xlim(xlim[0], xlim[1])
    
    ax.set_ylabel(ylabel, size=14)
    ax.set_xlabel(xlabel, size=14)
    
    
#dimensionality 

def get_curves_LD_all_kcoefs_ranges(matrix, log=False, norm=False,  numbins=20, nparray=False, 
                                    in_folder='', from_file=False, flat=False, reg_range_from_data=False):
    """Function to plot the number of pairs at a given distance against the distance
    
    Inputs:
    - matrix: 2d matrix of pairwise distances between sequences, if Flat=True,
        it can be an 1d array (flattened matrix) or a tuple (1d array, array size) that 
        also includes number of pairwisw distances
    - log: set to True to plot the curve in log scale (False by default, not used)
    - norm: set to True to normalize the curve (False by default, not used)
    - numbins: number of bins used to calculate the curves (20 by default)
    - nparray: set to True if the input is a 2d matrix that has to be read from a file 
        and was written to this file as a numpy array (False by default)
    - in_folder: if the input is a 2d matrix that has to be read from a file, it's the path 
        to this file (empty string by default)
    - from_file: set to True if the input is a 2d matrix that has to be read from a file 
        (False by default)
    - flat: set to True if the input is a 1d array (False by default)
    - reg_range_from_data: set to True to set the regression range equal to 
    [min(distances), max(distances)], instead of (0, 1) by default. Works
    only with log=False
    
    Outputs:
    - x: range of pairwise distances
    - y: number of pairs at a given distance
    """
    
    #read pairwise distance matrix as an 1d array
    if flat==True:
        if type(matrix)==tuple:
            k, size= matrix
        elif type(matrix)==np.ndarray:
            k= matrix
            size=len(k)
    else:
        k,size=matrix_to_array(matrix, nparray=nparray, in_folder=in_folder, from_file=from_file)
    
    if log==True:
        # using logspace so that the points appear to be equally spaced in the log scale
        binn=np.histogram(k, bins=list(np.append(np.zeros(1), np.logspace(-2.99573227,0,numbins,base=math.e))))
    elif reg_range_from_data:
        custombins=list(np.linspace(min(k), max(k), numbins+1))
        binn=np.histogram(k, bins=custombins)
    else:
        custombins=list(np.linspace(0,1, numbins+1))
        binn=np.histogram(k, bins=custombins)
    #get x axis (range of pairwise distances)
    dista=np.log(binn[1][1:])
    #get y axis (cumulative pairs counts)
    cumulative=np.cumsum(binn[0])
    
    #check if there are any empty bins at small distances
    #if there are any, these bins will be excluded (because np.log(0)=-inf which is a problem during dimension calculation)
    if 0 in cumulative:
        last0_ind=np.where(cumulative==0)[0][-1]+1
        cumulative=cumulative[last0_ind:]
        dista=dista[last0_ind:]
    #normalize y scale or not: DOES NOT alter dimension value in the end, only usefull for visualization purposes
    if norm==True:
        # applying the Normalization: log(E)/N^2
        return(dista, np.log(cumulative/size**2))
    else:
        return(dista, np.log(cumulative))
    
def max_k(matrix, log=False, win_size_perc=0.4, numbins=20, t_stud=None, mode='max', 
          norm=False, plot=False, in_folder='', out_folder='',out_plots='',out_file=True, xlim=None, ylim=None,
          nparray=False,from_file=True,save_picture=False, flat=False, label='', ax=None, regcolor='r', scattercolor='tab:blue',
          fname=None, ci_r2_out_flat=False, range_out_flat=False, reg_range_from_data=False, skipbin=False, 
          norange=False,legendsize=12, ylabel='log(N)', xlabel='log(Pairwise distance)', fontsize=14):
    
    """
    Function to calculate dimension from the matrix of pairwise distances.
    Supports three modes of calculation: 
    - MODE="max": dimension is the highest slope of the linear regression slopes done on all windows of size WIN_SIZE_PERC*NUMBINS
    - MODE="min_max": same as "max" mode + it gives a min_k value - the slope of the curve from 0% sequence distance to (but not including) the 
      first point included in region used to calculate max_k
    - MODE="min_max_50": function outputs the "max" and "min" dimension values and "0-50%" dimension 
      (the slope of the linear regression done on the first WIN_SIZE_PERC*NUMBINS (0-WIN_SIZE_PERC% divergence))
    Source: all_important_functions.ipynb
    """
    #get the file name prefix 
    if flat==False:
        fn=re.findall('(.+)\.', matrix)[0]
    elif fname!=None:
        fn=fname
    else:
        fn='simulated_matrix'
    
    #get x and y for regression, where
    #x - pairwise distance in log scale (0-100% divergence)
    #y - log of number of pairs per bin at a a given pairwise distance
    if skipbin:
        (x, y)=matrix
    else:
        (x, y)=get_curves_LD_all_kcoefs_ranges(matrix, log=log, norm =norm, numbins=numbins, 
                                               nparray=nparray, in_folder=in_folder, from_file = from_file, 
                                               flat=flat, reg_range_from_data=reg_range_from_data)

    slopes = []
    x_range = []
    rvals= []
    sterrs = []
    bs = []
    
    win_size=int(win_size_perc*numbins)
    div_incr=100/numbins #divergence increment
    length=len(x)
    
    #to find maximun slope check all intervals of the curve of length WIN_SIZE
    for i in range(len(x)-win_size+1):
        win_x=x[i:i+win_size]
        win_y=y[i:i+win_size]
        slope, b, r_value, p_value, std_err =  stats.linregress(win_x, win_y)
        #record only positive slopes
        if slope <= 0:
            continue
        slopes.append(slope)
        rvals.append(r_value)
        x_range.append(win_x)
        bs.append(b)
        sterrs.append(std_err)
        
    if len(slopes) == 0:
        raise ValueError(
            """
            No positive slopes or no windows
            with all values above threshold found!
            """
            )

    max_k = max(slopes)
    ind=slopes.index(max_k)
    max_x = x_range[ind]
    max_r = rvals[ind]
    max_b = bs[ind]
    
    #adjustment in case the actual length of the x axis is different from numbins (empty bins removed in get_curves_LD_all_kcoefs_ranges function)
    difflen=numbins-length    

    #get the range of the distances used to get the highest slope
    if log==True:
        x_coor_nolog=list(np.append(np.zeros(1), np.logspace(-2.99573227,0,numbins,base=math.e)))[difflen:]
        log_str='log'
    elif reg_range_from_data and flat and not skipbin:
        x_coor_nolog=list(np.linspace(np.exp(x[0])-(np.exp(x[1])-np.exp(x[0])), np.exp(x[-1]), numbins+1))[difflen:]
        log_str='non_log_range_from_data'
    else:
        x_coor_nolog=list(np.linspace(0,1, numbins+1))[difflen:]
        log_str='non_log'

    range_start=x_coor_nolog[ind]
    range_end=x_coor_nolog[ind+win_size]
    
    if mode=='min_max' or mode=='min_max_50':
        flat_x=x[:ind]
        flat_y=y[:ind]
        if len(flat_y)==0 or len(flat_x) ==0:
            fslope=fr_value=flat_ci=flat_end=flat_start=fstd_err=fb=np.nan
        else:
            fslope, fb, fr_value, fp_value, fstd_err =  stats.linregress(flat_x, flat_y)     
            flat_ci=fstd_err*t_student_one_sided(ind+1)
            
            flat_start=x_coor_nolog[0]
            flat_end=x_coor_nolog[ind-1]
        
    #calculate confidence interval for the slope
    #one sided T Student test is performed because our curve is cumulative, therefore slope will be > 0
    if t_stud!=None and difflen==0:
        #use precalculated t_student value
        max_ci=t_stud*sterrs[ind]
    else:
        t_stud=t_student_one_sided(win_size)
        max_ci=t_stud*sterrs[ind]
    
    if out_file==True: #write a text output file 
        #0-50 curve
        if mode=='min_max_50':
            win50_x=x[0:int(numbins*0.5)-difflen]
            win50_y=y[0:int(numbins*0.5)-difflen]
            k50, b50, r50, _, sterr50 =  stats.linregress(win50_x, win50_y)
            x50 =win50_x
            if win_size_perc==0.5 and difflen==0:
                ci50=t_stud*sterr50 
            else:
                t_stud50=t_student_one_sided(len(x[0:int(numbins*0.5)-difflen]))
                ci50=t_stud50*sterr50
                
            #writing outputs to file    
            text=['OG_id,0-50_k,0-50_R^2,0-50_CI,range_start,range_end,max_k,max_k_R^2,max_CI,min_k,min_k_R^2,min_CI']
            text.append('%s,{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n'.format(k50, r50**2, ci50, range_start, range_end, max_k, max_r**2, max_ci,fslope, fr_value**2, flat_ci) % fn)
            with open(out_folder+'%s_%s_%s_%s.coefnr' % (fn, log_str, str(win_size_perc), str(numbins)), 'w+') as out:
                out.write('\n'.join(text))
                
        elif mode=="max":
            #writing outputs to file
            text=['OG_id,range_start,range_end,max_k,max_k_R^2,max_CI']
            text.append('%s,{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n'.format(range_start, range_end, max_k, max_r**2, max_ci) % fn)
            with open(out_folder+'%s_%s_%s_%s.coefnr' % (fn, log_str, str(win_size_perc), str(numbins)), 'w+') as out:
                out.write('\n'.join(text))
                
        elif mode=="min_max":
            #writing outputs to file
            text=['OG_id,range_start,range_end,max_k,max_k_R^2,max_CI,min_k,min_k_R^2,min_CI']
            text.append('%s,{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n'.format(range_start, range_end, max_k, max_r**2, max_ci, fslope, fr_value**2, flat_ci) % fn)
            with open(out_folder+'%s_%s_%s_%s.coefnr' % (fn, log_str, str(win_size_perc), str(numbins)), 'w+') as out:
                out.write('\n'.join(text))
    
    
    #saving plot    
    if plot ==True and flat==False:
        if ax:
            pass
        else:
            fig, ax=plt.subplots()
        ax.scatter(x, y, label= label, c=scattercolor)
        #plot max slope
        if norange:
            ax.plot(max_x, max_k*max_x + max_b, label=r"$D_{cor}=$"+"{:.2f}$\pm${:.2f}".format(max_k, max_ci)+label, c=regcolor)
        else:
            ax.plot(max_x, max_k*max_x + max_b, label="{:.2f}-{:.2f}, y={:.2f}x+{:.2f}\nR^2={:.3f}, 95% CI={:.3f}".format(range_start, range_end, max_k, max_b,max_r**2, max_ci), 
                c=regcolor
                )
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        #plot min slope (regression on all the points from 0 (or the first non-empty bin) to but not including the first point of max slope)
        if mode=='min_max' or mode=='min_max_50':
            ax.plot(flat_x, fslope*flat_x + fb, label="{:.2f}-{:.2f}, y={:.2f}x+{:.2f}\nR^2={:.3f}, 95% CI={:.3f}".format(flat_start, flat_end, fslope, fb, fr_value**2, flat_ci), c='orange')
        #plot slope from 0 (or the first non-empty bin) to 50% distance
        if mode=='0-50' or mode=='min_max_50':
            ax.plot(x50, k50*x50 + b50, label="y={:.2f}x+{:.2f}\nR^2={:.3f}, 95% CI={:.3f}".format(k50, b50,r50**2, ci50), c='k')
        ax.legend(loc='lower right',fontsize=legendsize) 
        ax.set_ylabel(ylabel, size=fontsize)
        ax.set_xlabel(xlabel, size=fontsize)
        if save_picture==True:        
            plt.savefig(out_plots+'%s_%s_%s_%s.png' % (fn, log_str, str(win_size_perc), str(numbins)))
    elif plot ==True and flat==True:
        if ax:
            pass
        else:
            fig, ax=plt.subplots()
        ax.scatter(x, y, c=scattercolor)
        #plot max slope
        if norange:
            ax.plot(max_x, max_k*max_x + max_b, label=r"$D_{cor}=$"+"{:.2f}$\pm${:.2f}".format(max_k, max_ci)+label, c=regcolor)
        else:
            ax.plot(max_x, max_k*max_x + max_b, label="{:.2f}-{:.2f}, y={:.2f}x+{:.2f}\nR^2={:.3f}, 95% CI={:.3f}".format(range_start, range_end, max_k, max_b,max_r**2, max_ci)+label, c=regcolor)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        #plot min slope (regression on all the points from 0 (or the first non-empty bin) to but not including the first point of max slope)
        if mode=='min_max' or mode=='min_max_50':
            ax.plot(flat_x, fslope*flat_x + fb, label="{:.2f}-{:.2f}, y={:.2f}x+{:.2f}\n".format(flat_start, flat_end, fslope, fb)+label)
        #plot slope from 0 (or the first non-empty bin) to 50% distance
        if mode=='0-50' or mode=='min_max_50':
            ax.plot(x50, k50*x50 + b50, label="y={:.2f}x+{:.2f}\n".format(k50, b50)+label)
        ax.legend(loc='lower right',fontsize=legendsize) 
        ax.set_ylabel(ylabel, size=fontsize)
        ax.set_xlabel(xlabel, size=fontsize)

        if save_picture==True:        
            plt.savefig(out_plots+'%s_%s_%s_%s.png' % (fn, log_str, str(win_size_perc), str(numbins)))
    if flat==True:
        if ci_r2_out_flat:
            if range_out_flat:
                return(max_k, max_r**2, max_ci, range_start, range_end)
            else:
                return(max_k, max_r**2, max_ci)
        else:
            if range_out_flat:
                return(max_k, range_start, range_end)
            else:
                return(max_k)
    
#effective topological dimension

def add_topol_dim_cols(df, base_alpha=True, ci=True, seqlen_col='median_seq_len_orig', 
                       alpha_col='median_alphabet_size_less02gaps', 
                       dimcol='dim_ur', dimCIcol='dim_CI_ur', prefix=''):

    """
    Adds columns for effective topological dimension, 
    calculated from correlation dimension, maximum distance and sequence length
    """

    if seqlen_col=='num_var_sites_less02gaps' or seqlen_col=='num_var_sites':
        maxdist=df[seqlen_col]
    else:
        maxdist=df[seqlen_col]*df['max_ld']
    df[prefix+'n_seq_from_dim_log10']= df[dimcol]*np.log10(maxdist)
    df[prefix+'n_seq_from_dim_log20']= df[prefix+'n_seq_from_dim_log10'] / np.log10(20)
    if base_alpha:
        df[prefix+'n_seq_from_dim_log_alpha']= df[prefix+'n_seq_from_dim_log10'] / np.log10(df[alpha_col])
    if ci:
        # Confidence intervals for dim
        df[prefix+'log10_nseq_ci'] = df[dimCIcol] * np.log10(maxdist)
        df[prefix+'log20_nseq_ci'] = df[prefix+'log10_nseq_ci'] / np.log10(20)  

#truncated binomial distribution
def truncated_binomial_log_product_base20(L, p, n=20, base=20):
    """
    Compute log-base-20 of the product ∏ u^{n_u} where n_u = L * f_u,
    and f_u is the truncated binomial PMF with Bin(n, p), u ≥ 1.
    """
    q0 = (1 - p) ** n
    norm = 1 - q0  # normalization to truncate at u >= 1

    log_base_20 = 0.0
    for u in range(1, n + 1):
        prob_u = comb(n, u) * (p ** u) * ((1 - p) ** (n - u)) / norm
        log_base_20 += L * prob_u * (np.log(u) / np.log(base))  # change of base

    return log_base_20

    
def add_expected_topol_dim_cols(df, dimfrom='usage', base=20, use_col='median_alphabet_size_less02gaps', binomial=False, seqlen_col='median_seq_len_orig'):
    """
    Function that calculates expected effective topological dimension based on amino 
    acid usage or dN/dS.
    """
    if type(base)==int:
        base_name= str(base)
    else:
        base_name='_alpha'

    if dimfrom=='usage':
        df['n_seq_expected_from_{}_log{}'.format(dimfrom, base_name)] = df[seqlen_col]*np.log(df[use_col])/np.log(base)
        if binomial:
            df['n_seq_bin_expected_from_{}_log{}'.format(dimfrom, base_name)] = [truncated_binomial_log_product_base20(L, p, 20) for L, p in zip(df[seqlen_col], np.array(df[use_col])/20)]
    elif dimfrom=='dn/ds': 
        df['n_seq_expected_from_{}_log{}'.format(dimfrom, base_name)] = df[seqlen_col]*np.log(df[use_col]*20)/np.log(base)
        if binomial:
            df['n_seq_bin_expected_from_{}_log{}'.format(dimfrom, base_name)] = [truncated_binomial_log_product_base20(L, p, 20) for L, p in zip(df[seqlen_col], np.array(df[use_col]))]

def fraction_of_explored_space(df, prefix=''):
    df[prefix+'n_seq_diff_usage_log20']=(df['n_seq_bin_expected_from_usage_log20']-df[prefix+'n_seq_from_dim_log20'])
    df[prefix+'n_seq_diff_dnds_log20']= (df['n_seq_bin_expected_from_dn/ds_log20']-df[prefix+'n_seq_from_dim_log20'])
    df[prefix+'n_seq_diff_usage_log10']=df[prefix +'n_seq_diff_usage_log20']* math.log10(20)
    df[prefix+'n_seq_diff_dnds_log10']= df[prefix +'n_seq_diff_dnds_log20']* math.log10(20)
    return df

#simulate and subsample matrices    
    
def create_normal_matrix(dims, size=4000):
    '''
    Creates a matrix of the different distances of points normally distributed on 'dims' dimensions and of size 'size'.
    Source: all_important_functions.ipynb
    '''
    eucl_coords_by_dim=[]    
    for dim in range(1, dims+1):
        eucl_coords_by_dim.append(np.random.normal(0,1,size))
    eucl_coords_by_pt=np.array(eucl_coords_by_dim).T
    dist_arr=scipy_dist.pdist(eucl_coords_by_pt)
    dist_norm=dist_arr/max(dist_arr)
    
    return(dist_norm, size)

def submatrices(matrix_f,n,x): 
    """
    Inputs: 
    matrix - initial matrix as a csv file
    n - size of the matrices we get
    x - number of subsamples
    Source: all_important_functions.ipynb
    """
    if type(matrix_f)==str: #matrix is a csv file
        matrix= pd.read_csv(matrix_f,header = None).to_numpy()
    else: #matrix is a numpy array (NxN shaped)
        matrix=matrix_f
        
    # Make an exception if the size of the subsampled matrix equals
    # to the size of the matrix to subsample from. Then return the original one
    if matrix.shape[0]==n:
        subsamples=[matrix]*x
        return(subsamples)
    subsamples=[]
    for j in range(0,x):
        #Makes a random list of numbers within the interval (0, len(matrix)) so the numbers in list do not repeat.
        #We take k = ((len(matrix)) - (len of matrix we are getting)) of those numbers
        to_remove = random.sample(list(range(0,len(matrix))), k = len(matrix)-n)
        #Cuts of the particular columns and rows from the initial matrix Number of rows are from the list avobe.
        subsamples.append(np.delete(np.delete(matrix,to_remove,0),to_remove,1))
    return(subsamples)
  
def multiple_subsamples(matrix, sizes=[5000, 6000, 7000, 8000], reps=5):
    """
    Creates subsample-matrices of different sizes from a big matrix
    Source: all_important_functions.ipynb
    """
    for size in sizes:
        sub_list=submatrices(matrix,size,reps)
        for ind, subm in enumerate(sub_list):
            with open(matrix[:-7]+'_'+str(size)+'_'+str(ind)+'.pdistm', 'w+') as outs:
                np.savetxt(outs, subm, delimiter=",", fmt='%.7f')
                
def subsample_matrix_given_points(matrix, indices_to_keep):
    """
    Function to select from a matrix only rows and columns with given indices
    Indices of the rows and columns to keep are specified by indices_to_keep argument
    Source: cluster_seqs_by_ds_dev.ipynb
    """
    to_remove=np.delete(np.array(range(len(matrix))), indices_to_keep)
    return(np.delete(np.delete(matrix,to_remove,0),to_remove,1))


#reading yn00/PAML output table with dS, dN and dN/dS estimations
def open_yn00_rst(file):

    """
    Function to open PAML yn00 rst file and to be able to read correctly 
    the merged values from different columns in table (like 1.888112.7890)
    Source: cluster_seqs_by_ds_dev.ipynb
    """    
     
    def insert_space(match):
        return(' '.join(re.findall('\d+\.\d{4}', match.group(0))))
    
    with open(file, 'r') as opened_file:
        newtable=re.sub('\d+\.\d+\.\d+\S*', insert_space, opened_file.read())
    if newtable[:8]=="\tt\tkappa" or newtable[:7]=="t\tkappa" or newtable[:7]=="t kappa": #enable opening rst table with the header and with removed unneeded columns
        t=pd.read_csv(StringIO(newtable), sep='\s+', header=0)
    else:
        t=pd.read_csv(StringIO(newtable), sep='\s+', header=None)
        #for rst files produced with dnds_pairwise_ndata_n_allinfile_universal.py in "sites" mode (with estimated number of Syn and Nonsyn sites)
        if len(t.columns)==13: 
            t=t.drop([0, 1, 2, 3, 8, 11], axis=1)
        elif newtable[:3]==' YN': #for rst files produced with yn00 ndata = 1 option 
            t=t.drop([0, 5, 8], axis=1)
        else: #for rst files produced with yn00 ndata = n (pairs) option
            t=t.drop([0, 1, 6, 9], axis=1)
        t.columns=['t',  'kappa', 'dNdS', 'dN', 'dNSE', 'dS', 'dSSE']
    t.apply(pd.to_numeric)

    return(t)


def get_ds_for_clustering(dndstable, output='dS', input='yn00'):
    
    """Function which opens ynout table with open_yn00_rst() fuction, replaces Nan's and -0's 
    with 99's and 0's respectively and outputs either np.ndarray of dS (if output=='dS') or
    the whole rst dataframe (if output=='rst').
    Input: if input=='yn00': name of the rst file from yn00 program of PAML package (ynout file)
           if input=='ng86': the name of the table of Nei&Gojobori 1986 dnds estimates (3 columns: dnds, dn, ds, space-delimited)
    Output: either np.ndarray with dS or whole rst dataframe depending on the output parameter
    
    Source: cluster_seqs_by_ds_dev.ipynb
    """

    if input=='yn00':
        rawt=open_yn00_rst(dndstable)
        fixedt=rawt.replace(to_replace={np.nan:float(99), -0:0})
    elif input=='ng86':
        fixedt=pd.read_csv(dndstable, sep=' ', names=['dNdS', 'dN', 'dS'])

    if output=='dS':
        return(fixedt.dS)
    elif output=='rst':
        return(fixedt)
    

def multiple_files_generator(dslist, fastalist):
    """
    Function to run sequence clustering on multiple files consequently
    """
    #open lists of corresponding (by order) rst/dS files and fastas with sequences
    dsfiles=open(dslist).readlines()
    fastafiles=open(fastalist).readlines()
    
    #output pairs of ds and fasta
    for dsf, fastaf in tqdm(zip(dsfiles, fastafiles)):
        yield(dsf.strip(), fastaf.strip())
		