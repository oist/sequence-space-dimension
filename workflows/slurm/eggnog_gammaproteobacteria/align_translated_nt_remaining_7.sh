#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=align_rem_big
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/align_translated_nt/align_translated_nt_eggnog_gamma_remaining_7_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=90G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-7%7

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

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/
OUT=/flash/KondrashovU/ladaisa/eggnoggamma/align_translated_nt_main_1585/
OUTL=/flash/KondrashovU/ladaisa/logs/eggnoggamma/align_translated_nt/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/lists/eggNOGs_nt_translated_aligned_aa_remaining_7_to_report_gaps.list  | tail -1)

ml bioinfo-ugrp-modules  DebianMed/12.0
module load muscle/5.1.0-1

#align
muscle -super5 $INP/eggNOGs_nt_aligned_ungapped_translated/$INFILE -output $OUT$INFILE.ali > $OUTL/$INFILE.log 2>&1
