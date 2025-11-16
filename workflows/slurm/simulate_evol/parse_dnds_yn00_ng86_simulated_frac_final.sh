#!/bin/bash

#SBATCH --partition=short

#Give your job a name
#SBATCH --job-name=cont_paml_parse
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/paml/parse_yn00_ng86_contf_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=2:00:00

#Specify memory in gigabytes
#SBATCH --mem=60G

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

INP=/bucket/KondrashovU/seq_space/simulate_evol/
OUT=/flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
LOGOUT=/flash/KondrashovU/ladaisa/logs/

ml bioinfo-ugrp-modules  DebianMed/12.0
module load python/3.11.4

python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
 -f $INP/data/uniform_f_w_dist/paml_dnds_yn00_simu.list -id $OUT/paml/ynout_dir/ -od $OUT/paml/ \
 --ynout_dir $OUT/paml/ynout_dir_yn00_parsed_filtered/ -a 'yn00' -l 1 -t True -yo True
 
 python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
 -f $INP/data/uniform_f_w_dist/paml_dnds_ng86_simu.list -id $OUT/paml/ynout_dir/ -od $OUT/paml/ \
 --ynout_dir $OUT/paml/ynout_dir_ng86_parsed_filtered/ -s '.ng86' -a 'ng86' -l 1 -t True -yo True 

