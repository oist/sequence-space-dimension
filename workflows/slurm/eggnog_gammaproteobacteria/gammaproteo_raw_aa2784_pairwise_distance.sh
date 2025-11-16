#!/bin/bash

#SBATCH --partition=short

#Give your job a name
#SBATCH --job-name=hdist_gamma
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/hdir/pdist_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=2:00:00

#Specify memory in gigabytes
#SBATCH --mem=20G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-2784%100

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=lada.isakova@oist.jp
#SBATCH --mail-type=ALL

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

module load python/3.10.2

INP=/bucket/KondrashovU/seq_space/
OUT=/flash/KondrashovU/ladaisa/eggnoggamma/hdir/

#Run a binary that takes nth line of input.txt as input

python3 $INP/scripts/hamming_distance_multiple_ali_universal.py -a $(head -$SLURM_ARRAY_TASK_ID $INP/eggnog/gammaproteo/lists/all_gammaproteo_eggNOGs_raw_aa_200_seqs_or_more.list | tail -1) -if $INP/eggnog/gammaproteo/1236_raw_algs/1236/ -of $OUT/ -d True
