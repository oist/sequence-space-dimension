#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np
from tqdm import tqdm

with open(sys.argv[1], 'r') as filelist:
	l=filelist.read().split("\n")[:-1]

means=[]
medians=[]
vars=[]
for i in tqdm(l):
	pdist=pd.read_csv(sys.argv[2]+i, header=None).to_numpy()
	dist=list(pdist[np.triu_indices(np.shape(pdist)[0], k = 1)])
	means.append(np.mean(dist))
	medians.append(np.median(dist))
	vars.append(np.var(dist))

df=pd.DataFrame({'OG_name':l, 'LD_mean':means, 'LD_median':medians, 'LD_var':vars})

df.to_csv(sys.argv[1]+'_LD.tsv', sep='\t', index=False)
