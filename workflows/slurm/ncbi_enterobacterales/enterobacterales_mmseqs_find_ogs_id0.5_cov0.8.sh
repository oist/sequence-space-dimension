#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=enterobact_clust
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/mmseqs_entrobacterales/mmseqss_cluster.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=48:00:00

#Specify memory in gigabytes
#SBATCH --mem=100G

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=l.isakova@oist.jp
#SBATCH --mail-type=ALL

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#Load modules, if needed

INP=/bucket/KondrashovU/seq_space/ncbi_enterobacterales/genomes_db/
OUT=/flash/KondrashovU/ladaisa/mmseqs_entrobacterales/

ml bioinfo-ugrp-modules  DebianMed/12.0  mmseqs2/14-7e284+ds-1+b2

mmseqs cluster $INP/Enterobacterales.db $OUT/Enterobacterales_clustered $OUT/tmp --min-seq-id 0.5 --cov-mode 0 -c 0.8 --threads ${SLURM_CPUS_PER_TASK}
