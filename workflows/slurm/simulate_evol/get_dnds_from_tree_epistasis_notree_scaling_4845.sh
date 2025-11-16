#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=dnds_tree
#SBATCH --output=/flash/KondrashovU/ladaisa/simulate_evol/epistasis/dnds_from_tree_epistasis_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=50G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#SBATCH --array=1-100%100

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

OUT=/flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
INP=/bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/

FILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_final_trees_4845_notreefinal.list | tail -1)


module load python/3.11.4

#get dnds from tree
python /bucket/KondrashovU/seq_space/simulate_evol/get_dnds_from_simulated_tree.py -t $INP/lists/${FILE} -a False -if $INP/trees/ -of $OUT -s True