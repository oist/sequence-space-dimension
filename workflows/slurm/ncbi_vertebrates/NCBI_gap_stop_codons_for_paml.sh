#!/bin/bash

#Give your job a name
#SBATCH --job-name=dnds_ncbi
#SBATCH --output=dnds/tmp/logs/%A_%a.log

#Specify time limit; max is 10 days, i.e. 240 hours
#SBATCH --time=240:00:00

#Specify memory in gigabytes
#SBATCH --mem=5G

#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-type=ALL

#Submit a job array of N jobs limiting the number of simultateously running ones to K
#SBATCH --array=1-14583%200

#Do not restart the job if it fails
#SBATCH --no-requeue

#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

#For single-CPU jobs, make sure that they use a single thread
export OMP_NUM_THREADS=1

#Load modules, if needed

module load python3

#Run a binary that takes nth line of input.txt as input
python3 gap_stop_codons_for_paml.py -f $(head -$SLURM_ARRAY_TASK_ID dnds_NT_ali_list.txt | tail -1) -if alignments_nt/ -of dnds/NT_alignments_no_stop_codons/
