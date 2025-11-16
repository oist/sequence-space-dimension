#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=simu_dim_unique
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/dim_unique_03_04_and_reg_range_w_anc.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=50G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Submit a job array of N jobs limiting the number of simultateously running ones to K

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

OUT=/flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
INP=/bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/

for f in $INP/pdistm_unique_renamed/*pdistm; do 
	BASE=$(basename $f)
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $BASE -i $INP/pdistm_unique_renamed/ -of $OUT/kcoefs_range_from_data_unique/ -w 0.5 -n 20 -r True -p False -s False -m 'max' -l False
	python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f $BASE -i $INP/pdistm_unique_renamed/ -of $OUT/kcoefs_unique/ -w 0.3 -n 20 -r False -p False -s False -m 'max' -l False
done 
