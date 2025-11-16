#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=hmm
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/hmm_entero_gamma_remaining2181.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=3-20:00:00

#Specify memory in gigabytes
#SBATCH --mem-per-cpu=10G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-4%4

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/hmmer_link_gamma_to_entero/
SINP=/bucket/KondrashovU/seq_space/ncbi_enterobacterales/genomes_db/
OUT=/flash/KondrashovU/ladaisa/entero/hmmer/
LOUT=/flash/KondrashovU/ladaisa/logs/entero/
FILE=all_eggnog_gamma546_remaining$SLURM_ARRAY_TASK_ID

#load hmmer 
module load bioinfo-ugrp-modules  DebianMed/12.0 hmmer/3.3.2+dfsg-1

hmmsearch --cpu 25 --noali --tblout $OUT/$FILE.tsv -o $OUT/$FILE.out -E 1e-10 \
--tformat fasta $INP/$FILE.hmm $SINP/Enterobacterales.prt &> $LOUT/$FILE_hmmsearch.log

