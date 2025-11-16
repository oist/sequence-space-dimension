#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=align
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/dnds/dnds_yn00_ng86_nosplit_enterobacterales200seqs_2597_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=60G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-650%650

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

INP=/bucket/KondrashovU/seq_space/ncbi_enterobacterales/
OUT=/flash/KondrashovU/ladaisa/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_entero_2597_nt_ali_wo_stops.list   | tail -1)

ml bioinfo-ugrp-modules  DebianMed/12.0
module load python/3.11.4

for i in $(seq 1 4); do python3 /bucket/KondrashovU/seq_space/scripts/dnds_pairwise_ndata_n_allinfile_universal.py \
-f $(head -$i $INP/lists_split/$INFILE | tail -1) -fd $INP/dnds/NT_alignments_no_stop_codons/ -wd $OUT/entero/dnds/ -a 'both' -k False -s False -c False ; \
 done >> $OUT/logs/entero/dnds/$INFILE.log 2>&1
