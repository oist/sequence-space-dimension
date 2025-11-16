#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse 
import os
import pandas as pd
import numpy as np
import re
from tqdm import tqdm
from io import StringIO 
import math
import sys
import itertools 
from Bio import SeqIO
sys.path.append("/bucket/KondrashovU/seq_space/scripts/seq_space_lib/")
import sequence_space_lib as seqsp

def filter_dnds(ynout, dndsalg='yn00', dnclose=True, seq_file='', dnds_more1_filt=False):


    if dndsalg=='yn00':
        t=seqsp.open_yn00_rst(ynout)
        len0=len(t.index)
        
        if seq_file!='':
            # open fasta with sequences
            with open(seq_file, 'r') as infile:
                allseq = list(SeqIO.parse(infile, 'fasta'))   
                # get sequences ids as list 

            fastaids = [rec.id for rec in allseq]
            combn = list(itertools.combinations(fastaids, 2))

            t=t.join(pd.DataFrame.from_records(combn, columns =['seq1', 'seq2']))


        #omit rows with nan or NA
        t.dropna(axis=0, how='any', subset=['dNdS', 'dN', 'dNSE','dS', 'dSSE'], inplace=True)
        t=t[t['dN']<90] #remove dn=99 (similar to nan)
        #get % of pairs filtered out because they contained NA 
        lenwona=len(t.index)
        percna=(len0-lenwona)/len0

    elif dndsalg=='ng86':
        t=pd.read_csv(ynout, sep=' ', names=['dNdS', 'dN', 'dS'])
        len0=lenwona=len(t.index)
        percna=None

        if seq_file!='':
            # open fasta with sequences
            with open(seq_file, 'r') as infile:
                allseq = list(SeqIO.parse(infile, 'fasta'))   
            fastaids = [rec.id for rec in allseq]
            combn = list(itertools.combinations(fastaids, 2))
            t=t.join(pd.DataFrame.from_records(combn, columns =['seq1', 'seq2']))

    #filter out extreme values of dS
    t=t[t['dS']<0.8]
    #get % of pairs filtered out because the dS was > 0.8 
    lendshigh=len(t.index)
    percdshigh=(lenwona-lendshigh)/len0
    
    #get % of pairs filtered out because the dS was < 0.1
    t=t[t['dS']>0.1]
    lendslow=len(t.index)
    percdslow=(lendshigh-lendslow)/len0
    
    if dndsalg=='yn00':
        if dnclose==True:
            #removing rows if dSSE >= dS 
            t=t[t['dS']>t['dSSE']]
            #removing rows if dNSE >= dN 
            t=t[t['dN']>t['dNSE']]
            percsehigh=(lendslow-len(t.index))/len0
    elif dndsalg=='ng86':
        percsehigh=None
        #check if there are any dn/ds equal to -1 with dN != 0
        if (t[(t[['dNdS']] < 0).all(1)].dN==0).all()!=True:
            print('negative dN/dS, nonzero dN detected in %s' % ynout)

        #filter out rows with negative dN
        t=t[t['dN']>=0]
        #set all dnds equal to -1 to 0
        t['dNdS'] = t['dNdS'].replace(to_replace={-1:0})
        t= t[t['dNdS'] >=0]
    else:
        percsehigh=None
    
    if dnds_more1_filt==True:
        t= t[t['dNdS']<=1]

    return(t, len0, len(t.index), percna, percdshigh, percdslow, percsehigh)

if __name__=='__main__':
    
    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--file_list', help='list with names of yn00 rst files (both ndata=1 and ndata=n formats are supported)',
                        required=True)
    
    #Optional
    PARSER.add_argument('-a', '--dndsalg', help='PAML algorithm used to get the dnds table: yn00 or ng86', required=False, default='yn00')
    PARSER.add_argument('-l', '--levels', type=eval, help='number of levels or groupings on which to calculate summary stats, 1 (lowest) lewel by default, 2 max', 
                        required=False, default='1')
    PARSER.add_argument('-id', '--input_dir', help='folder with input .ynout files', required=False, default='')
    PARSER.add_argument('-od', '--output_dir', help='path to write the output table', required=False, default='')
    
    PARSER.add_argument('-dn', '--dnclose', type=eval, choices=[True, False],
                        help='set to False if the removal of rows with dN values close to dN SE is not desired, True by default', 
                        required=False, default='True')
    PARSER.add_argument('-s', '--suffix', 
                        help='part of the file name to exclude from the final table e.g if for file AAA_XXX_YY.ynout you only want to keep the AAA, set it to XXX_YY', 
                        required=False, default='')
    PARSER.add_argument('-sf', '--seq_file', 
                        help='file with list of the fasta file names with sequences used for dn/ds calculation. Note that the order and number of sequences in these fastas should be identical to the ones usdd for PAML', 
                        required=False, default='')
    PARSER.add_argument('-yd', '--ynout_dir', 
                        help='folder to record filtered ynout tables (if ynout_bool==True)', 
                        required=False, default='')
    PARSER.add_argument('-yo', '--ynout_bool', 
                        help='file with list of the fasta file names with sequences used for dn/ds calculation. Note that the order and number of sequences in these fastas should be identical to the ones usdd for PAML', 
                        type=eval, choices=[True, False], required=False, default='False')
    PARSER.add_argument('-t', '--sum_tables', 
                        help='boolean, set to true if outputting summary tables is desired, set to False if you only need filtered ynouts', 
                        type=eval, choices=[True, False], required=False, default='True')
    PARSER.add_argument('-p', '--no_pos_sel', 
                        help='boolean, set to true if want to remove all dnds>1 from final table, default - False', 
                        type=eval, choices=[True, False], required=False, default='False')
    
 
    ARGS = vars(PARSER.parse_args())
    
    file_list = ARGS['file_list']
    dndsalg = ARGS['dndsalg']
    input_dir = ARGS['input_dir']
    output_dir = ARGS['output_dir']
    dnclose = ARGS['dnclose']
    suffix = ARGS['suffix']
    levels= ARGS['levels']
    seq_file=ARGS['seq_file']
    ynout_dir= ARGS['ynout_dir']
    ynout_bool = ARGS['ynout_bool']
    sum_tables = ARGS['sum_tables']
    no_pos_sel = ARGS['no_pos_sel']

    #main
    with open(file_list, 'r') as filel:
        fl=filel.readlines()
    if seq_file!='':
        with open(seq_file, 'r') as sfilel:
            sfl=sfilel.read().splitlines()
    else:
        sfl=['']*len(fl)

    if no_pos_sel:
        no_pos_sel_suffix='dnds_less_than_1_only_'
        no_pos_sel_suffix_yn='_dnds_less_than_1_only'
    else:
        no_pos_sel_suffix=no_pos_sel_suffix_yn=''

    ids=[]
    ogids=[]
    ogdfs=[]
    medians=[]
    means=[]
    var=[]
    lengths1=[]
    lengths2=[]
    dnmeans=[]
    dsmeans=[]
    dnmedians=[]
    dsmedians=[]
    dnvars=[]
    dsvars=[]
    pdshigh=[]
    pdslow=[]
    pna=[]
    psehigh=[]

    if sum_tables==True:
        for ind, i in tqdm(enumerate(fl)):
            og_dfs=[]
            df, len1, len2, percna, percdshigh, percdslow, percsehigh = filter_dnds(input_dir+i[:-1], dndsalg=dndsalg, 
                                                                                    dnclose=dnclose, seq_file=sfl[ind],
                                                                                    dnds_more1_filt=no_pos_sel)
            if df.empty:
                print(i+': no pairs after filtering')
                continue
            if suffix=='_': #means that the id is the part of a file name before the first underscore
                ids.append(re.findall('(^[^_]+)_.*\.ynout', i)[0])
            else:
                ids.append(re.findall('(\S+)%s\.ynout' % suffix, i)[0])
            medians.append(np.median(df.dNdS))
            means.append(np.mean(df.dNdS))
            var.append(np.var(df.dNdS))
            lengths1.append(len1)
            lengths2.append(len2)
            dnmeans.append(np.mean(df.dN))
            dsmeans.append(np.mean(df.dS))
            dnmedians.append(np.median(df.dN))
            dsmedians.append(np.median(df.dS))
            dnvars.append(np.var(df.dN))
            dsvars.append(np.var(df.dS))
            pdshigh.append(percdshigh)
            pdslow.append(percdslow)
            if dndsalg!='ng86':
                pna.append(percna)
                psehigh.append(percsehigh)

            if levels==2:
                og, qid = re.findall('([^_]+)_([^/]+)%s\.ynout' % suffix, i)[0]
                ogids.append(og)
                df['OG_id']=[og]*len2
                df['query_id']=[qid]*len2
                ogdfs.append(df)

            elif ynout_bool==True : #record ynout table for each level 1 group after filtering
                with open(ynout_dir+i[:-1].replace('.ynout', no_pos_sel_suffix_yn+'_parsed_filtered.ynout'), 'w+') as ogdfrecord:
                    df.to_csv(ogdfrecord,  index=False, sep='\t')
            
        
        if dndsalg!='ng86':
            mvdf=pd.DataFrame({'OG_id':ids, 'dnds_median':medians, 'dnds_mean':means, 'dnds_var':var, 
                        'orig_length':lengths1,'filt_length':lengths2, 
                        'dn_mean':dnmeans, 'ds_mean':dsmeans, 'dn_median':dnmedians, 'ds_median':dsmedians, 'dn_var':dnvars, 'ds_var':dsvars,
                        'frac_NA':pna, 'frac_dS>0.8':pdshigh, 'frac_dS<0.1':pdslow, 'frac_SE_high':psehigh})
        else:
            mvdf=pd.DataFrame({'OG_id':ids, 'dnds_median':medians, 'dnds_mean':means, 'dnds_var':var, 
                        'orig_length':lengths1,'filt_length':lengths2, 
                        'dn_mean':dnmeans, 'ds_mean':dsmeans, 'dn_median':dnmedians, 'ds_median':dsmedians, 'dn_var':dnvars, 'ds_var':dsvars,
                        'frac_dS>0.8':pdshigh, 'frac_dS<0.1':pdslow})
            
        file_list_no_path=re.findall('([^/]+$)',file_list)[0]

        with open(output_dir+'dnds_dn_ds_means_medians_var_fracfilt_'+no_pos_sel_suffix+file_list_no_path+'.csv', 'w+') as handle:
            mvdf.to_csv(handle, index=False, float_format='%.4f')

    else:
        for ind, i in tqdm(enumerate(fl)):
            og_dfs=[]
            df, len1, len2, percna, percdshigh, percdslow, percsehigh = filter_dnds(input_dir+i[:-1], dndsalg=dndsalg, 
                                                                                    dnclose=dnclose, seq_file=sfl[ind],
                                                                                    dnds_more1_filt=no_pos_sel)
            if df.empty:
                print(i+': no pairs after filtering')
                continue

            og, qid = re.findall('([^_]+)_([^/]+)%s\.ynout' % suffix, i)[0]
            ogids.append(og)
            df['OG_id']=[og]*len2
            df['query_id']=[qid]*len2
            ogdfs.append(df)

    if levels==2:

        if sum_tables==True:
            ogsummary=[]
            alldfs=pd.concat(ogdfs,  ignore_index=True)
            for oggr in list(set(ogids)):
                #select only the rows from a given OG
                targetdf=alldfs[alldfs['OG_id']==oggr]
                
                #record ynout table for each level 1 group after filtering
                if ynout_bool==True:
                    with open(ynout_dir+oggr+no_pos_sel_suffix_yn+'_parsed_filtered.ynout', 'w+') as ogdfrecord:
                        targetdf.to_csv(ogdfrecord,  index=False, sep='\t')

                #get summary statistics of selected columns
                ssdf=targetdf[['dNdS', 'dN', 'dS']].describe().loc[['mean', 'std', '50%']]
                #rename "50%" into median
                ssdf=ssdf.rename(index={'50%':'median'})
                #flatten the table into 1 row
                v = ssdf.unstack().to_frame().T
                v.columns = v.columns.map('_'.join)
                #add number of pairs into the table
                v['length']=len(targetdf.index)
                v['OG_id']=oggr
                ogsummary.append(v)
            
            ogfinal=pd.concat(ogsummary,  ignore_index=True)
            ogfinal.to_csv(output_dir+'dnds_dn_ds_means_medians_var_by_OG_'+no_pos_sel_suffix+file_list_no_path+'.csv', index=False, float_format='%.4f')

        else:
            alldfs=pd.concat(ogdfs,  ignore_index=True)
            for oggr in list(set(ogids)):
                #select only the rows from a given OG
                targetdf=alldfs[alldfs['OG_id']==oggr]
                
                #record ynout table for each level 1 group after filtering
                if ynout_bool==True:
                    with open(ynout_dir+oggr+no_pos_sel_suffix_yn+'_parsed_filtered.ynout', 'w+') as ogdfrecord:
                        targetdf.to_csv(ogdfrecord,  index=False, sep='\t')

