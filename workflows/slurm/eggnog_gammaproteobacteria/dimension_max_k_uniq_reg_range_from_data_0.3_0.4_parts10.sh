#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=k_gamma
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/pdistm_logs/dim_reg_range_0304_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=20:00:00

#Specify memory in gigabytes
#SBATCH --mem=50G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-2%2

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

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/
OUT=/flash/KondrashovU/ladaisa/eggnoggamma/

#Run a binary that takes nth line of input.txt as input

INLIST=$(head -$SLURM_ARRAY_TASK_ID $INP/lists/egnogg_gamma_raw_pdistm_lists.list  | tail -1)

while read -r INFILE; do 
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $INFILE -i $INP/pdistm_matrices_unique/ -of $OUT/kcoefs_mult_ali_max_k_range_from_data/ -op $OUT/kcoefs_mult_ali_max_k_range_from_data_plots/ -w 0.4 -n 20 -r True -p True -s True -m 'max' -l False
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $INFILE -i $INP/pdistm_matrices_unique/ -of $OUT/kcoefs_mult_ali_max_k_range_from_data/ -op $OUT/kcoefs_mult_ali_max_k_range_from_data_plots/ -w 0.3 -n 20 -r True -p True -s True -m 'max' -l False
done < $INP/lists/$INLIST
