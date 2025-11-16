#!/bin/bash

#SBATCH --partition=compute

#Give your job a name
#SBATCH --job-name=filt_rename_simu_pdistm
#SBATCH --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/filt_rename_simu_pdistm.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=10:00:00

#Specify memory in gigabytes
#SBATCH --mem=50G

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

ml python/3.11.4

INP=/bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/pdistm/
OUT=/flash/KondrashovU/ladaisa/simulate_evol/

#Run a binary that takes nth line of input.txt as input

for ind in $(seq 1 12444); do  line=$(ls $INP/*pdistm|head -${ind}|tail -1); BASE=$(basename $line '.pdistm');\
 python3 /bucket/KondrashovU/seq_space/scripts/filter_identical_seqs_from_pdistm.py -m ${BASE}.pdistm -if $INP \
 -of $OUT/pdistm_unique_renamed/ ; mv $OUT/pdistm_unique_renamed/${BASE}_unique.pdistm $OUT/pdistm_unique_renamed/SimuOG${ind}.pdistm; \
 echo "${BASE},SimuOG${ind}" >> $OUT/simu_og_pdistm_names_mapping.csv; done 

zip -r $OUT/simu_og_unique_renamed_pdistm.zip $OUT/pdistm_unique_renamed/

