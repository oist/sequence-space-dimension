#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=k_ent_range
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/dim_reg_range_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=20:00:00

#Specify memory in gigabytes
#SBATCH --mem=10G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-97%97

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

INP=/bucket/KondrashovU/seq_space/ncbi_enterobacterales/
OUT=/flash/KondrashovU/ladaisa/

#Run a binary that takes nth line of input.txt as input

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_enterobacterales_2605_aligned_fastas_200seqs.list | tail -1)

for i in $(seq 1 27); do
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $(head -$i $INP/$INFILE | tail -1).pdistm -i $INP/pdistm_matrices_unique_levenshtein/ -of $OUT/entero/kcoefs_mult_ali_max_k_range_from_data/  -w 0.5 -n 20 -r True -p True -s True -m 'max' -l False -op $OUT/entero/kcoefs_mult_ali_max_k_range_from_data_plots/;
done 

