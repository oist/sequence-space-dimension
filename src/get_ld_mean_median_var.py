#!/usr/bin/env python3

import pandas as pd
import os
import re
from tqdm import tqdm
import numpy as np
import sys

if len(sys.argv)>1:
    folder=sys.argv[1]
    out_suffix=sys.argv[2]
else:
    folder='/nfs/scistore08/kondrgrp/lisakova/cogs/cut_cogs/pdistm_matrices/'
    out_suffix='cut_cogs'

def matrix_to_array(matrix):
    array=[]
    pdist=pd.read_csv(matrix, header=None).to_numpy()
    for i, row in enumerate(pdist):
        array+=list(row[i+1:])
    return(array)


if __name__=='__main__':
    cld=os.listdir(folder)
    cmeanld=[]
    cmedianld=[]
    cvarld=[]
    fn=[]
    for ci in tqdm(cld):
        if ci[-7:]=='.pdistm':
            fn.append(re.findall('(.+)\.', ci)[0])
            cldcsv=matrix_to_array(folder+ci)
            cmeanld.append(np.mean(cldcsv))
            cmedianld.append(np.median(cldcsv))
            cvarld.append(np.var(cldcsv))
    p=pd.DataFrame({'OG_id':fn,'mean_ld':cmeanld, 'median_ld':cmedianld, 'var_ld':cvarld})
    p.to_csv('allstats_mean_median_var_ld_%s.csv' % out_suffix, index=False)
