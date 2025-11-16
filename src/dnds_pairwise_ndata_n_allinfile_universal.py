#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    This is a function to get the dN, dS and dN/dS values for all pairs of sequences in the alignments.
    Progrm runs alignment on all sequences in orthogroup at once. 
    First, sequences are recorded in the temporary phylip_paml.phy file in phylip format. 
    Second, dn/ds are calculated using yn00 program of PAML. Only rst file with dN +-SE, dS +-SE and dN/dS 
    values (calculated using Yang & Nielsen (2000) method) is recorded, all other outputs of yn00 are deleted.
    
    "split" mode - writes all combinations of pairs to several separate files (each - 50 000 pairs max), 
    runs yn00 consequtively on each of these files
    
    "sites" mode - adds the number of synonymous and non-synonymous sites in the alignment of the pair of sequences to the standard rst output 
    Output table will look like:  "seq. seq.     S       N        t   kappa   omega (dN/dS)     dN +- SE    dS +- SE"
    
    Arguments and inputs:
    - FASTA - nucleotide alignment name
    - WORKING_DIR - temporary folders (with yn00 input and original outputs) for dnds calculation for each 
    alignment will be created here; this is also a directory where should be the output directories ynout_dir/ 
    (for dn and ds values) and logs/ (for yn00 log files)
    - DNDSALG - {'yn00', 'ng86'} - PAML algorithm to use when estimating dn and ds 
    - FASTADIR
    - KEEP_ALL {False, True} - if True, keeps all PAML outputs
    - SPLIT - split all sequence pairs in batches of SPLIT_SIZE so that program will not crash on big datasets
    - SPLIT_SIZE - number to pairs of sequences on which probgram calculates dn/ds in one run, default = 1000
    - SITES - {False, True} - False (default) - no counts (MLE) of syn and non-syn sites in the output, True - S and N counts estimation in output
    
    Outputs:
    *.ynout - resulting output file is a table with columns:  YN:   t   kappa   omega (dN/dS)     dN +- SE    dS +- SE
    
    Required folder structure:
    - FASTADIR/ (location of the folder with original fasta alignments)
    - WORKING_DIR/
      \-logs/ (yn00 stdout)
      \-ynout_dir/ (renamed rst files for each alignment with dn/ds values)
      \-ctl file (created by the programm itself, it's a yn00 control file, see PAML manual)
    
    ynout_dir is not reqired if KEEP_ALL is set to True (program keeps all PAML outputs)
    
    Usage: python3 ../scripts/dnds_pairwise_ndata_n_allinfile_universal.py -f 1RM7J_short_nostop.fasta 
    -wd ~/Documents/kondrashov_lab/seq_space/dnds/eggntest_sites/ 
    -fd ~/Documents/kondrashov_lab/seq_space/dnds/eggntest_sites/ -k True -s False -c True
    """

from Bio import AlignIO, Align, SeqIO
import argparse
import shutil
import itertools
import re
import os

def format_phylip_run_yn00(fasta, working_dir, dndsalg='yn00', keep_all=False, fastadir='', split=False, split_size=1000, sites=False):
    
    """The function to get the dN, dS and dN/dS values for all pairs of sequences in the alignments."""
    
    #entering the working directory 
    os.chdir(working_dir) 
    
    name=re.findall('(.+)\.', fasta)[0]
    
    #creating temporary file-specific directory to run PAML (it does not allow setting output file names
    #separate directory is needed to run PAML on different files in parallel
    if os.path.isdir(name):
        pass
    else:
        os.mkdir(name)
        
    #reading fasta
    s=list(SeqIO.parse(fastadir+fasta, format='fasta'))
    
    comb=list(itertools.combinations(range(len(s)), 2))
    
    #if split==True, the program writes a certain number of sequence pairs to file and runs yn00 on this sequences until the dn/ds on all pairs are calulated
    if split==True:
        split_size=int(split_size)
        #mini-function to split all pairs of sequences in bathces of given size X
        final_list= lambda test_list, x: [test_list[i:i+x] for i in range(0, len(test_list), x)]
        comb=final_list(comb, split_size)

        #entering temporary file-specific directory
        os.chdir(name)

        ali_len=len(s[0].seq)

        #creating list where all pairwise dnds values will be added 

        for ind, k in enumerate(comb):
            combali=[]
            for i in k:
            #rewriting alignment in PAML-legible phylip format for each pair of sequences
            #writing alignment in PAML-legible phylip format to temporary file (universal name so that yn00.ctl file is common for all alignments-files)
                combali.append(' 2 %d\n' % ali_len +'\n'.join(['nameseq1', str(s[i[0]].seq), 'nameseq2', str(s[i[1]].seq)])+'\n\n')

            with open('phylip_paml.phy', 'w+') as out:
                out.write(''.join(combali))

            if ind==0:
                ctltext=['seqfile = phylip_paml.phy\n', 'outfile = pamlout\n', 'verbose = 1\n', 'icode = 0\n', 'weighting = 0\n', 'commonf3x4 = 0\n']
                ctltext.append('ndata = %s\n' % str(len(k)))

                with open('yn00.ctl', 'w+') as ynctl:
                    ynctl.write(''.join(ctltext))

            #running yn00 function from PAML package
            os.system('yn00 yn00.ctl > ../logs/log_%s_%d' % (name, ind))
            #record the results
            if dndsalg=='yn00' or dndsalg=='both':
                if sites==False:
                    #moving output file with dn, ds, and dn/ds values from temporary directory to output directory and renaming it from general name to alignment-specific
                    #resulting output file is a table with columns:  YN:   t   kappa   omega     dN +- SE    dS +- SE
                    shutil.move('rst', '../ynout_dir/'+name+'_'+str(ind)+'.ynout')
                    
                else:
                    #parsing the pamlout file to get the rst table but with the S and N values (also no YN and pair number in the beginning of the line)
                    os.system("grep 'seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE' -A2  pamlout | grep -v -e 'seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE\|--\|^$' > %s" % '../ynout_dir/'+name+'_'+str(ind)+'.ynout')
            #if we need the dn and ds estimated using Nei & Gojobori 1986 algorithm
            if dndsalg=='ng86' or dndsalg=='both':
                os.system("grep 'Nei & Gojobori 1986. dN/dS (dN, dS)' -A5  pamlout | grep  'nameseq2'  |sed 's/nameseq2//g'| sed -e 's/(//g' -e 's/)//g' -e 's/^[ \t]*//' > %s" % '../ynout_dir/'+name+'_'+str(ind)+'.ng86.ynout')

        #exiting the temporary file-specific directory and deleting it and its content
        os.chdir('../')
        shutil.rmtree(name)
    
    #if split==False program writes all combinations of sequences to 1 file and runs yn00 on all combinations at once
    elif split==False:
        #entering temporary file-specific directory
        os.chdir(name)

        ali_len=len(s[0].seq)

        #creating list where all pairwise dnds values will be added 
        combali=[]
        for i in comb:
            #rewriting alignment in PAML-legible phylip format for each pair of sequences
            #writing alignment in PAML-legible phylip format to temporary file (universal name so that yn00.ctl file is common for all alignments-files)
            combali.append(' 2 %d\n' % ali_len +'\n'.join(['nameseq1', str(s[i[0]].seq), 'nameseq2', str(s[i[1]].seq)])+'\n\n')

        with open('phylip_paml.phy', 'w+') as out:
            out.write(''.join(combali))

        ctltext=['seqfile = phylip_paml.phy\n', 'outfile = pamlout\n', 'verbose = 1\n', 'icode = 0\n', 'weighting = 0\n', 'commonf3x4 = 0\n']
        ctltext.append('ndata = %s\n' % str(len(comb)))

        with open('yn00.ctl', 'w+') as ynctl:
            ynctl.write(''.join(ctltext))

        #running yn00 function from PAML package
        os.system('yn00 yn00.ctl > ../logs/log_%s' % name)

        if keep_all==False:
            #if we need the dn and ds estimated using Yang and Nielsen 2000 algorithm
            if dndsalg=='yn00' or dndsalg=='both':
                #exiting the temporary file-specific directory and deleting it and its content
                #moving output file with dn, ds, and dn/ds values from temporary directory to output directory and renaming it from general name to alignment-specific
                #resulting output file is a table with columns:  YN:   t   kappa   omega     dN +- SE    dS +- SE
                if sites==False:
                    #getting the output without syn and non-syn sites counts
                    shutil.move('rst', '../ynout_dir/'+name+'.ynout')
                else:
                    #parsing the pamlout file to get the rst table but with the S and N values (also no YN and pair number in the beginning of the line)
                    os.system("grep 'seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE' -A2  pamlout | grep -v -e 'seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE\|--\|^$' > %s" % '../ynout_dir/'+name+'.ynout')

            #if we need the dn and ds estimated using Nei & Gojobori 1986 algorithm
            if dndsalg=='ng86'  or dndsalg=='both':
                os.system("grep 'Nei & Gojobori 1986. dN/dS (dN, dS)' -A5  pamlout | grep  'nameseq2'  |sed 's/nameseq2//g'| sed -e 's/(//g' -e 's/)//g' -e 's/^[ \t]*//' > %s" % '../ynout_dir/'+name+'.ng86.ynout')
            
            os.chdir('../')
            shutil.rmtree(name)        
        else:    
            if dndsalg=='ng86' or dndsalg=='both':
                
                os.system("grep 'Nei & Gojobori 1986. dN/dS (dN, dS)' -A5  pamlout | grep  'nameseq2'  |sed 's/nameseq2//g'| sed -e 's/(//g' -e 's/)//g' -e 's/^[ \t]*//' > %s" % '../ynout_dir/'+name+'.ng86.ynout')

if __name__=='__main__':
    
    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-f', '--fasta', help='fasta file with alignment',
                        required=True)
    PARSER.add_argument('-wd', '--working_dir', help='working directory with ctl file', required=True)
    
    #Optional
    PARSER.add_argument('-a', '--dndsalg', choices=['yn00', 'ng86', 'both'], 
                        help='PAML algorithm used to estimate dn and ds, yn00 - Yang & Nielsen 2000, ng89 - Nei & Gojobori 1986',
                        required=False, default='yn00')

    PARSER.add_argument('-fd', '--input_fasta_dir', help='folder with input alignments', required=False, default='')
    
    PARSER.add_argument('-k', '--keep_all', type=eval, choices=[True, False], 
                        help='boolean variable, False - delete all outputs except for rst file, True - keep all yn00 outputs',
                        required=False, default='False')
    
    PARSER.add_argument('-ss', '--split_size', 
                        help='number of pairs of sequences (alignment chuncks size) on which to separately run yn00, default - 1000 pairs', 
                        required = False, default=1000)
    
    PARSER.add_argument('-s', '--split', type=eval, choices=[True, False],
                        help='boolean variable indicating whether to run yn00 on all pairs of sequences at once (False) or on chunks of length split_size',
                        required = False, default='False')
    
    PARSER.add_argument('-c', '--sites', type=eval, choices=[True, False],
                        help='boolean variable, False (default) - no counts (MLE) of syn and non-syn sites in the output, True - S and N counts estimation in output',
                        required = False, default='False')
    
    ARGS = vars(PARSER.parse_args())
    
    fasta = ARGS['fasta']
    working_dir = ARGS['working_dir']
    dndsalg = ARGS['dndsalg']
    fastadir = ARGS['input_fasta_dir']
    keep_all = ARGS['keep_all']
    split_size=ARGS['split_size']
    split=ARGS['split']
    sites=ARGS['sites']

    format_phylip_run_yn00(fasta, working_dir, dndsalg=dndsalg, keep_all=keep_all, 
                           fastadir=fastadir, split_size=split_size, split=split, sites=sites)
    

