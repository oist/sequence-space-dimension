#!/bin/bash

#SBATCH --partition=short

#Give your job a name
#SBATCH --job-name=dim_ur_cogs_full
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/cogs/full/pairali_dimension_ur_uniq.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=2:00:00

#Specify memory in gigabytes
#SBATCH --mem=50G

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

INP=/bucket/KondrashovU/seq_space/cogs/full/
OUT=/flash/KondrashovU/ladaisa/cogs/full/


while read -r INFILE; do 
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $INFILE -i $OUT/muscle/pdistm_matrices_unique_pairwise_ali/ -of $OUT/kcoefs_mult_ali_max_k_range_from_data_pairali/ -op $OUT/kcoefs_mult_ali_max_k_range_from_data_pairali_plots/ -w 0.5 -n 20 -r True -p True -s True -m 'max' -l False
done < $INP/pairali_done30102025.list


