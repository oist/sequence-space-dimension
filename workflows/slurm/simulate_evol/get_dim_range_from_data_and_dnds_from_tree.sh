#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=dim_range
#SBATCH --output=/flash/KondrashovU/ladaisa/simulate_evol/epistasis/dnds_dim_range_from_tree_epistasis_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=100G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#SBATCH --array=1-5%5

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

OUT=/flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
INP=/bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/

FLIST="None_exp_rate0.05349_scale0.0_notree_0.7_gammaFitAll_any_continuousF.notree.list\nNone_exp_rate0.05349_scale0.0_notree_1_gammaFitAll_any_continuousF.notree.list\nNone_exp_rate0.05349_scale0.0_tree_1_exp_any_continuousF.list\nexp_node_order_rate_exp_rate0.1_scale0.0_tree_2.3_exp_any_continuousF.shortleaves.list\nexp_node_order_rate_rev_exp_rate0.05_scale1.0tree_0.05_exp_any_continuousF.longleaves.list"

FILE=$(printf "$FLIST\n" | head -$SLURM_ARRAY_TASK_ID | tail -1)

module load python/3.11.4

#get dimensionality 
echo $FILE
python /bucket/KondrashovU/seq_space/simulate_evol/get_dim_dist_from_simulated_tree.py -t $INP${FILE} -a False -d True -df False -do True -r True -w 0.5 -if $INP/trees/ -of $OUT
echo "done with dimensionality, proceeding to the dnds from tree"
#get dnds from tree
python /bucket/KondrashovU/seq_space/simulate_evol/get_dnds_from_simulated_tree.py -t $INP${FILE} -a False -if $INP/trees/ -of $OUT
