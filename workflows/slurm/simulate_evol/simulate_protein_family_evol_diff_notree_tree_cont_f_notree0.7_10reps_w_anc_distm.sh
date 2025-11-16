#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=notree10contF
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/uniform_frac/simulate_final_notree10contF_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=450G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#SBATCH --array=1-48%48

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

OUT=/flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
INP=/bucket/KondrashovU/seq_space/simulate_evol/

ALLPARAMS=$(head -$SLURM_ARRAY_TASK_ID $INP/params/evol_all_tree_notree_w_epistasis_notree0.7_10reps_cont_f.list | tail -1)

IFS=';' read -r params replicates f_type exp_rate branch_len_scaling scale_factor <<< $ALLPARAMS

module load python/3.11.4

python $INP/simulate_protein_family_diff_tree_w_epistasis_evol_stats_and_ali.py -p $params -of $OUT \
-r $replicates -ft $f_type -fu 0.7 -kt False -er $exp_rate -bs $branch_len_scaling -bf $scale_factor \
-d True -a True -df True > $OUT/pyout_$params$branch_len_scaling.log 2>&1
