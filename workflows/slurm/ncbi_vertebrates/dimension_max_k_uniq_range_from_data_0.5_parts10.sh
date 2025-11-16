#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=vert_uk_reg_range
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/ncbivert/kcoefs_mult_ali_max_k/dimension_reg_range_uniq_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=15:00:00

#Specify memory in gigabytes
#SBATCH --mem=10G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-10%10

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

ml python/3.10.2

INP=/bucket/KondrashovU/seq_space/ncbi_vertebrates/
OUT=/flash/KondrashovU/ladaisa/ncbivert/

INLIST=$(head -$SLURM_ARRAY_TASK_ID $INP/lists/pdistm_matrices_unique_lists.list | tail -1)

while read -r INFILE; do 
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $INFILE -i $INP/pdistm_matrices_unique/ -of $OUT/kcoefs_mult_ali_max_k_range_from_data_unique/ -w 0.5 -n 20 -r True -p False -s False -m "max" -l False
done < $INP/lists/$INLIST