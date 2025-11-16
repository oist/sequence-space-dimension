#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=pali_gam
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/pairali/pairwise_ali_muscle_continue_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=10G

#SBATCH --array=1-134%25

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

ml python/3.11.4
ml bioinfo-ugrp-modules  DebianMed
module load muscle/5.1.0-1

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/
OUT=/flash/KondrashovU/ladaisa/eggnoggamma/pairali/

INL=random_0.05_og_ids_gamma_for_pairali.list

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/$INL | tail -1)

python3 /bucket/KondrashovU/seq_space/scripts/hamming_distance_pairwise_ali_universal.py -f $INFILE -l True -d True -if $INP/1236_raw_algs/1236/ -of $OUT/pdistm_matrices_unique_pairwise_ali/ -a muscle -w 5000 -p from_data -g True  &>>  $OUT/gamma_pairwise_ali_dist_new_$INFILE.log


