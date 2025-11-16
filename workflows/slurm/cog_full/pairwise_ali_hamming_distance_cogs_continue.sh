#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=pali_cogs_full
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/cogs/full/pairwise_ali_muscle_continue_$SLURM_JOB_ID.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=16G

#SBATCH --array=1-296%296

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
ml bioinfo-ugrp-modules  DebianMed/12.0
module load muscle/5.1.0-1

INP=/bucket/KondrashovU/seq_space/cogs/full/
OUT=/flash/KondrashovU/ladaisa/cogs/full/

INLIST=$(head -$SLURM_ARRAY_TASK_ID $INP/lists/list_of_lists_cogs_num_of_seq_200_seqs_or_more_split.list | tail -1)

while read -r INFILE; do 
	python3 /bucket/KondrashovU/seq_space/scripts/hamming_distance_pairwise_ali_universal.py -f $INFILE -l True -d True -if $INP/sequences/ -of $OUT/muscle/pdistm_matrices_unique_pairwise_ali/ -a muscle -w 10000 -p from_data  &>>  $OUT/muscle/cogs_pairwise_ali_dist_$SLURM_ARRAY_TASK_ID.log
done < $INP/lists/${INLIST}

