#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=pali_vert
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/ncbivert/pairali/pairwise_ali_muscle_continue_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=10G

#SBATCH --array=1-616%80

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

INP=/bucket/KondrashovU/seq_space/ncbi_vertebrates/
OUT=/flash/KondrashovU/ladaisa/ncbivert/pairali/

INL=random_0.05_og_ids_vert_for_pairali_orig_aa.list

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/$INL | tail -1)

python3 /bucket/KondrashovU/seq_space/scripts/hamming_distance_pairwise_ali_universal.py -f $INFILE -l True -d True -if $INP/making_ali/sequences_in_all_diam_folder/OG_AA_sample/ -of $OUT/pdistm_matrices_unique_pairwise_ali/ -a muscle -w 5000 -p from_data  &>>  $OUT/vert_pairwise_ali_dist_new_$INFILE.log


