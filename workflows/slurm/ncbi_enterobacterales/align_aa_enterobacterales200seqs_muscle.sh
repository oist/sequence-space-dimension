#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=align
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/align_muscle_200seqs/align_aa_enterobacterales200seqs_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=2-20:00:00

#Specify memory in gigabytes
#SBATCH --mem=20G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-97%97

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
OUT=/flash/KondrashovU/ladaisa/entero/align_muscle_200seqs/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_enterobacterales_2605_fastas_200seqs.list  | tail -1)

ml bioinfo-ugrp-modules  DebianMed/12.0
module load muscle/5.1.0-1

#align
while read -r line; do muscle -align $INP/fastas_200seqs/$line -output $OUT$line.ali; done < $INP/$INFILE >> /flash/KondrashovU/ladaisa/logs/entero/align_muscle_200seqs/$INFILE.log 2>&1



