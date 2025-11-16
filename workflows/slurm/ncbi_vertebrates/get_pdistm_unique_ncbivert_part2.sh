#!/bin/bash

#SBATCH --partition=short

#Give your job a name
#SBATCH --job-name=ncbi_unique2
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/ncbivert/get_unique_pdistm_ncbi_2.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=02:00:00

#Specify memory in gigabytes
#SBATCH --mem=10G

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

#Run a binary that takes nth line of input.txt as input

while read -r line; do python3 /bucket/KondrashovU/seq_space/scripts/filter_identical_seqs_from_pdistm.py -m $line -if $INP/pdistm_matrices/ -of $OUT/pdistm_matrices_unique/ ; echo $line >> $OUT/get_unique_pdistm_2.log; done < $INP/lists/pdistm_matrices_all_wo_folder_part2.list
