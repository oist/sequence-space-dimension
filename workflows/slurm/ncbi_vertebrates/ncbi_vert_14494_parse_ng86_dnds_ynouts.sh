#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=parse_dnds
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/ncbivert/recalculate_dnds_ng86_yn00/parse_ng86_dnds14494.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4:00:00

#Specify memory in gigabytes
#SBATCH --mem=100G

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

module load python/3.11.4

INP=/bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/
OUT=/flash/KondrashovU/ladaisa/

#Run a binary that takes nth line of input.txt as input

python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
-f $INP/ncbi_vert_dnds_ng86_ynout.list -id $INP/ynout_dir_ng86/ -od $OUT/ncbivert/recalculate_dnds_ng86_yn00/ \
-s '_' -a 'ng86' -l 1

