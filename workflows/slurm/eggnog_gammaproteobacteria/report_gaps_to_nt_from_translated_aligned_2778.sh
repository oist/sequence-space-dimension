#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=repgaps
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/NT_alignments_no_stop_codons_from_nt/report_gaps_to_nt_2778_%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=20:00:00

#Specify memory in gigabytes
#SBATCH --mem=20G

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-96%96

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

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/
OUT=/flash/KondrashovU/ladaisa/eggnoggamma/

INFILE=$(head -$SLURM_ARRAY_TASK_ID $INP/list_of_lists_eggNOGs_nt_translated_aligned_aa_2778.list  | tail -1)

ml bioinfo-ugrp-modules  DebianMed/12.0
ml python/3.11.4

#insert gaps from aa alignemnts into nucleotides to create the nt alignments
for i in $(seq 1 29); do java -jar /apps/unit/KondrashovU/ladaisa/macse_v2.07.jar -prog reportGapsAA2NT \
-align_AA $INP/eggNOGs_nt_translated_aligned_aa/$(head -$i $INP/$INFILE | tail -1).ali -seq $INP/eggNOGs_nt_aligned_ungapped/$(head -$i $INP/$INFILE | tail -1) \
-out_NT $OUT/eggNOGs_nt_translated_aligned_nt/$(head -$i $INP/$INFILE | tail -1).nt.ali; \
python3 /bucket/KondrashovU/seq_space/scripts/gap_stop_codons_for_paml.py -f $(head -$i $INP/$INFILE | tail -1).nt.ali \
-if $OUT/eggNOGs_nt_translated_aligned_nt/ -of $OUT/NT_alignments_no_stop_codons_from_nt/ ; done >> /flash/KondrashovU/ladaisa/logs/eggnoggamma/NT_alignments_no_stop_codons_from_nt/log_$SLURM_ARRAY_TASK_ID.NT.log 2>&1
