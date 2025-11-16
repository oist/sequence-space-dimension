#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=vt_dnds
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/ncbivert/recalculate_dnds_ng86_yn00/dnds_yn00_ng86_14495_nosplit_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=4-00:00:00

#Specify memory in gigabytes
#SBATCH --mem=80G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l-isakova@oist.jp
#SBATCH --mail-type=ALL

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-806%806

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

module load python/3.11.4

INP=/bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/
OUT=/flash/KondrashovU/ladaisa/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_dnds_NT_ali_list_no_stop_codons_no_paralogs.list | tail -1)

#Run a binary that takes nth line of input.txt as input

for i in $(seq 1 18); do python3 /bucket/KondrashovU/seq_space/scripts/dnds_pairwise_ndata_n_allinfile_universal.py \
-f $(head -$i $INP/lists_split/$INFILE | tail -1) -fd $INP/NT_alignments_no_stop_codons/ -wd $OUT/ncbivert/recalculate_dnds_ng86_yn00/ -a 'both' -k False -s False -c False ; \
 done >> $OUT/logs/ncbivert/recalculate_dnds_ng86_yn00/$INFILE.log 2>&1
 
 
