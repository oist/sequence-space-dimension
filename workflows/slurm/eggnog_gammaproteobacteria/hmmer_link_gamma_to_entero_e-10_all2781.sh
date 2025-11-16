#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=hmm
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/entero/hmm_entero_gamma.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=3-20:00:00

#Specify memory in gigabytes
#SBATCH --mem=100G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=51

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

INP=/bucket/KondrashovU/seq_space/eggnog/gammaproteo/hmmer_link_gamma_to_entero/
SINP=/bucket/KondrashovU/seq_space/ncbi_enterobacterales/genomes_db/
OUT=/flash/KondrashovU/ladaisa/entero/hmmer/
LOUT=/flash/KondrashovU/ladaisa/logs/entero/

#load hmmer 
module load bioinfo-ugrp-modules  DebianMed/12.0 hmmer/3.3.2+dfsg-1

hmmsearch --cpu 50 --noali --tblout $OUT/all_eggnog_gamma2781.tsv -o $OUT/all_eggnog_gamma2781.out -E 1e-10 \
--tformat fasta $INP/all_eggnog_gamma2781.hmm $SINP/Enterobacterales.prt &> $LOUT/all_eggnog_gamma2781_hmmsearch.log


