#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=parseent
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/dnds/parse_yn00_ng86_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=20:00:00

#Specify memory in gigabytes
#SBATCH --mem=40G

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

INP=/bucket/KondrashovU/seq_space/ncbi_enterobacterales/
OUT=/flash/KondrashovU/ladaisa/

#Run a binary that takes nth line of input.txt as input

python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
 -f $INP/dnds/entero_ynout_yn00_2597.list -id $INP/dnds/ynout_dir_yn00/ -od $OUT/entero/dnds/ \
 --ynout_dir $OUT/entero/dnds/ynout_dir_yn00_parsed_filtered_seqnames/ -s '.fasta.nt_nostop' -a 'yn00' -l 1 -t True -yo True \
 -sf $INP/dnds/entero_ynout_corr_nt_ali_for_seqnames_2597.list
 
 python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
 -f $INP/dnds/entero_ynout_ng86_2597.list -id $INP/dnds/ynout_dir_ng86/ -od $OUT/entero/dnds/ \
 --ynout_dir $OUT/entero/dnds/ynout_dir_ng86_parsed_filtered_seqnames/ -s '.fasta.nt_nostop.ng86' -a 'ng86' -l 1 -t True -yo True \
 -sf $INP/dnds/entero_ynout_corr_nt_ali_for_seqnames_2597.list


