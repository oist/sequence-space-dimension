#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=cont_paml
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/paml/dnds_paml_finalcontinuousF_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=1-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=60G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-27%27

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

INP=/bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/
OUT=/flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
LOGOUT=/flash/KondrashovU/ladaisa/logs/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_list_final_all_fastas_5244.list | tail -1)

ml bioinfo-ugrp-modules  DebianMed/12.0
module load python/3.11.4

for i in $(seq 1 200); do python3 /bucket/KondrashovU/seq_space/scripts/dnds_pairwise_ndata_n_allinfile_universal.py \
-f $(head -$i $INP/lists/$INFILE | tail -1) -fd $INP/fasta/ -wd $OUT/paml/ -a 'both' -k False -s False -c True ; \
 done >> $LOGOUT/simulate_evol/paml/$INFILE.log 2>&1
