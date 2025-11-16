#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=from_nt_dnds
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/yn00_from_nt/dnds2771_nosplit_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=60G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l-isakova@oist.jp
#SBATCH --mail-type=ALL

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-278%278

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

module load python/3.11.4

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/
OUT=/flash/KondrashovU/ladaisa/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_nt_alignments_no_stop_codons_from_nt_2771.list | tail -1)

#Run a binary that takes nth line of input.txt as input

for i in $(seq 1 10); do python3 /bucket/KondrashovU/seq_space/scripts/dnds_pairwise_ndata_n_allinfile_universal.py \
-f $(head -$i $INP/$INFILE | tail -1) -fd $OUT/eggnoggamma/NT_alignments_no_stop_codons_from_nt/ -wd $OUT/eggnoggamma/yn00_from_nt/ -a 'yn00' -k False -s False -c True ; \
 done >> $OUT/logs/eggnoggamma/yn00_from_nt/log_$SLURM_ARRAY_TASK_ID.log 2>&1




