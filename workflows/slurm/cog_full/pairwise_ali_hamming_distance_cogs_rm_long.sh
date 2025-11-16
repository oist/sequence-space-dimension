#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=pali_cogs_rmlong
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/cogs/full/pairwise_ali_muscle_rmlong_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=25G

#SBATCH --array=0-7%8

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

INP=/bucket/KondrashovU/seq_space/cogs/full/
OUT=/flash/KondrashovU/ladaisa/cogs/full/

INL=(COG1196.fa COG2931.fa COG3170.fa COG3210.fa COG3319.fa COG3321.fa COG5276.fa COG5281.fa)

python3 /bucket/KondrashovU/seq_space/scripts/hamming_distance_pairwise_ali_universal.py -f "${INL[$SLURM_ARRAY_TASK_ID]}" -l True -d True -if $INP/sequences/ -of $OUT/muscle/pdistm_matrices_unique_pairwise_ali/ -a muscle -w 10000 -p 0 -fl 10000 &>>  $OUT/muscle/cogs_pairwise_ali_dist_rmlong_$SLURM_ARRAY_TASK_ID.log
