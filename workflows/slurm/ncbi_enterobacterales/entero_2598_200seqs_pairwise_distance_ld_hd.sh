#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=hdist_ent
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/pdist/pdist_aa_enterobacterales200seqs_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=20:00:00

#Specify memory in gigabytes
#SBATCH --mem=20G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-10%10

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

module load python/3.11.4

INP=/bucket/KondrashovU/seq_space/
OUTH=/flash/KondrashovU/ladaisa/entero/pdistm_matrices_unique_hamming/
OUTL=/flash/KondrashovU/ladaisa/entero/pdistm_matrices_unique_levenshtein/

#Run a binary that takes nth line of input.txt as input

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/ncbi_enterobacterales/list_of_lists_enterobacterales_2598_aa_ali_200seqs.list | tail -1)

for i in $(seq 1 261); do python3 $INP/scripts/hamming_distance_multiple_ali_universal.py -a $(head -$i $INP/ncbi_enterobacterales/$INFILE | tail -1) \
-if $INP/ncbi_enterobacterales/aa_ali_200seqs/ -of $OUTL/ -d True -aa False ; done

for i in $(seq 1 261); do python3 $INP/scripts/hamming_distance_multiple_ali_universal.py -a $(head -$i $INP/ncbi_enterobacterales/$INFILE | tail -1) \
-if $INP/ncbi_enterobacterales/aa_ali_200seqs/ -of $OUTH/ -d True -aa True -s aa_only ; done


